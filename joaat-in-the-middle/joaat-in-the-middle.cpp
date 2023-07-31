#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <future>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#include <immintrin.h>

#include "parallel.h"

using usize = std::size_t;
using u8 = std::uint8_t;
using u32 = std::uint32_t;

using String = std::string;
using StringView = std::string_view;

template <typename T>
using Vec = std::vector<T>;

template <typename Key>
using HashSet = std::unordered_set<Key>;

template <typename T, typename U>
using Pair = std::pair<T, U>;

static u32 joaat_partial(u32 hash, const void* data, usize length)
{
    for (usize i = 0; i != length; ++i) {
        hash += static_cast<const unsigned char*>(data)[i];
        hash += hash << 10; // hash *= 0x401
        hash ^= hash >> 6;
    }

    return hash;
}

static u32 joaat_partialr(u32 hash, const void* data, usize length)
{
    for (usize i = length; i != 0; --i) {
        // inverse of hash ^= hash >> 6;
        hash ^= (hash >> 6) ^ (hash >> 12) ^ (hash >> 18) ^ (hash >> 24) ^ (hash >> 30);

        // inverse of hash += hash << 10;
        hash *= 0xC00FFC01;

        // inverse of hash += data[i];
        hash -= static_cast<const unsigned char*>(data)[i - 1];
    }

    return hash;
}

static u32 joaat_definalize(u32 hash)
{
    // inverse of hash += hash << 15;
    hash *= 0x3FFF8001;

    // inverse of hash ^= hash >> 11;
    hash ^= (hash >> 11) ^ (hash >> 22);

    // inverse of hash += hash << 3;
    hash *= 0x38E38E39;

    return hash;
}

using FilterWord = u32;

struct Collider {
    Vec<Vec<String>> Parts {};
    Vec<u32> Suffixes {};

    Vec<Vec<u32>> Prefixes {};

    Vec<Pair<const String*, usize>> CurrentParts;

    usize PrefixPos {};
    usize SuffixPos {};

    Vec<FilterWord> Filter {};

    Vec<u32> HashIndices;
    Vec<u32> HashBuckets;
    Vec<u8> SubHashes;

    Collider(u32 seed, Vec<u32> hashes, Vec<Vec<String>> parts);
    ~Collider();

    void PushPrefix(const String* suffixes, usize suffix_count);
    void PopPrefix();

    void PushSuffix(const String* prefixes, usize prefix_count);

    String GetPrefix(usize index);
    String GetSuffix(usize index);

    void Compile(usize prefix_table_size, usize suffix_table_size);

    usize Match(HashSet<String>& found);
    usize Collide(HashSet<String>& found);
};

template <typename T>
static inline bool bit_test(const T* bits, usize index)
{
    constexpr usize Radix = sizeof(T) * CHAR_BIT;

#if defined(_MSC_VER) && !defined(__clang__) && (defined(_M_IX86) || defined(_M_X64))
    // `index % Radix` should be a no-op on x86, as it is masked off by shl/shr/bt
    // Yet MSVC still generates a horrible `movzx, and, movzx` sequence
    // Note also, _bittest{64} is not recommended due to its high latency
    return (bits[index / Radix] & (T(1) << (index /*% Radix*/))) != 0;
#else
    return (bits[index / Radix] & (T(1) << (index % Radix))) != 0;
#endif
}

template <typename T>
static inline void bit_set(T* bits, usize index)
{
    constexpr usize Radix = sizeof(T) * CHAR_BIT;

    bits[index / Radix] |= (T(1) << (index % Radix));
}

Collider::Collider(u32 seed, Vec<u32> hashes, Vec<Vec<String>> parts)
{
    Parts = std::move(parts);

    Prefixes.resize(Parts.size() + 1);
    Prefixes[0] = { seed };

    Suffixes = std::move(hashes);

    PrefixPos = 0;
    SuffixPos = Parts.size();

    CurrentParts.resize(Parts.size());
}

Collider::~Collider()
{
}

static void ExpandHashes(const u32* __restrict input, u32* __restrict output, usize count, const void* __restrict suffix, usize suffix_length)
{
    parallel_partition(count, 0x10000, 0, [=](usize start, usize count) {
        if (count >= 4) {
            do {
                __m128i hashes = _mm_loadu_si128((const __m128i*)&input[start]);

                for (usize i = 0; i != suffix_length; ++i) {
                    hashes = _mm_add_epi32(hashes, _mm_set1_epi32(static_cast<const unsigned char*>(suffix)[i]));
                    hashes = _mm_add_epi32(hashes, _mm_slli_epi32(hashes, 10));
                    hashes = _mm_xor_si128(hashes, _mm_srli_epi32(hashes, 6));
                }

                _mm_storeu_si128((__m128i*)&output[start], hashes);

                start += 4;
                count -= 4;
            } while (count >= 4);
        }

        for (usize i = 0; i != count; ++i) {
            usize j = start + i;

            output[j] = joaat_partial(input[j], suffix, suffix_length);
        }

        return true;
    });
}

void Collider::PushPrefix(const String* suffixes, usize suffix_count)
{
    const Vec<u32>& prefixes = Prefixes[PrefixPos];
    CurrentParts[PrefixPos] = { suffixes, suffix_count };

    ++PrefixPos;
    Vec<u32>& hashes = Prefixes[PrefixPos];

    usize prefix_count = prefixes.size();
    hashes.resize(prefix_count * suffix_count);

    for (usize i = 0; i < suffix_count; ++i) {
        StringView suffix = suffixes[i];
        ExpandHashes(prefixes.data(), &hashes[i * prefix_count], prefix_count, suffix.data(), suffix.size());
    }
}

void Collider::PopPrefix()
{
    --PrefixPos;
}

static void ShrinkHashes(const u32* __restrict input, u32* __restrict output, usize count, const void* __restrict prefix, usize prefix_length)
{
    parallel_partition(count, 0x10000, 0, [=](usize start, usize count) {
        if (count >= 4) {
            do {
                __m128i hashes = _mm_loadu_si128((const __m128i*)&input[start]);

                for (usize i = prefix_length; i != 0; --i) {
                    __m128i value = _mm_set1_epi32(static_cast<const unsigned char*>(prefix)[i - 1]);

                    hashes = _mm_xor_si128(hashes, _mm_srli_epi32(hashes, 6));
                    hashes = _mm_xor_si128(hashes, _mm_srli_epi32(hashes, 12));
                    hashes = _mm_xor_si128(hashes, _mm_srli_epi32(hashes, 24));

                    hashes = _mm_sub_epi32(hashes, _mm_slli_epi32(hashes, 10));
                    hashes = _mm_add_epi32(hashes, _mm_slli_epi32(hashes, 20));

                    hashes = _mm_sub_epi32(hashes, value);
                }

                _mm_storeu_si128((__m128i*)&output[start], hashes);

                start += 4;
                count -= 4;
            } while (count >= 4);
        }

        for (usize i = 0; i != count; ++i) {
            usize j = start + i;

            output[j] = joaat_partialr(input[j], prefix, prefix_length);
        }

        return true;
    });
}

void Collider::PushSuffix(const String* prefixes, usize prefix_count)
{
    usize suffix_count = Suffixes.size();

    Vec<u32> suffixes;
    suffixes.resize(suffix_count * prefix_count);

    for (usize i = 0; i < prefix_count; ++i) {
        StringView prefix = prefixes[i];
        ShrinkHashes(Suffixes.data(), &suffixes[i * suffix_count], suffix_count, prefix.data(), prefix.size());
    }

    Suffixes.swap(suffixes);

    --SuffixPos;
    CurrentParts[SuffixPos] = { prefixes, prefix_count };
}

String Collider::GetPrefix(usize index)
{
    String result;

    usize prefix_count = Prefixes[PrefixPos].size();

    for (usize i = PrefixPos; i != 0; --i) {
        const auto& [suffixes, suffix_count] = CurrentParts[i - 1];
        prefix_count /= suffix_count;

        result.insert(0, suffixes[index / prefix_count]);
        index %= prefix_count;
    }

    return result;
}

String Collider::GetSuffix(usize index)
{
    String result;

    usize suffix_count = SubHashes.size();

    for (usize i = SuffixPos; i != Parts.size(); ++i) {
        const auto& [prefixes, prefix_count] = CurrentParts[i];
        suffix_count /= prefix_count;

        result.append(prefixes[index / suffix_count]);
        index %= suffix_count;
    }

    return result;
}

static void SortHashesWithIndices(u32* hashes, u32* indices, usize count, u32 bit)
{
    if (count < 16) {
        for (usize i = 1; i < count; ++i) {
            u32 hash = hashes[i];
            u32 index = indices[i];

            usize j = i;

            for (; (j != 0) && (hashes[j - 1] > hash); --j) {
                hashes[j] = hashes[j - 1];
                indices[j] = indices[j - 1];
            }

            hashes[j] = hash;
            indices[j] = index;
        }

        return;
    }

    usize pivot = count;

    for (usize i = 0; i < pivot; ++i) {
        u32 hash = hashes[i];

        if ((hash >> bit) & 0x1) {
            u32 index = indices[i];

            do {
                --pivot;

                if (i == pivot)
                    break;

                std::swap(hash, hashes[pivot]);
                std::swap(index, indices[pivot]);
            } while ((hash >> bit) & 0x1);

            hashes[i] = hash;
            indices[i] = index;
        }
    }

    if (bit) {
        bit -= 1;

        const auto sort_lower = [=] { SortHashesWithIndices(hashes, indices, pivot, bit); };
        const auto sort_upper = [=] { SortHashesWithIndices(hashes + pivot, indices + pivot, count - pivot, bit); };

        if ((count > 0x10000) && (bit > 25)) {
            auto future = std::async(std::launch::async, sort_lower);
            sort_upper();
            future.wait();
        } else {
            sort_lower();
            sort_upper();
        }
    }
}

using Stopwatch = std::chrono::high_resolution_clock;

void Collider::Compile(usize prefix_table_size, usize suffix_table_size)
{
    printf("Compiling...\n");

    while (PrefixPos != SuffixPos) {
        const Vec<String>& next_prefix = Parts[PrefixPos];
        const Vec<String>& next_suffix = Parts[SuffixPos - 1];

        usize next_prefix_size = Prefixes[PrefixPos].size() * next_prefix.size();
        usize next_suffix_size = Suffixes.size() * next_suffix.size();

        bool more_prefixes = next_prefix_size < prefix_table_size;
        bool more_suffixes = next_suffix_size < suffix_table_size;

        if (more_prefixes && more_suffixes) {
            more_prefixes = next_prefix_size < next_suffix_size;
            more_suffixes = !more_prefixes;
        }

        if (more_prefixes) {
            printf("Expanding Prefixes %zu\n", PrefixPos);
            PushPrefix(next_prefix.data(), next_prefix.size());
        } else if (more_suffixes) {
            printf("Expanding Suffixes %zu\n", SuffixPos - 1);
            PushSuffix(next_suffix.data(), next_suffix.size());
        } else {
            break;
        }
    }

    usize prefix_count = Prefixes[PrefixPos].size();
    usize suffix_count = Suffixes.size();

    HashIndices.resize(suffix_count);

    for (usize i = 0; i < suffix_count; ++i)
        HashIndices[i] = static_cast<u32>(i);

    auto start = Stopwatch::now();

    printf("Building suffix lookup...\n");
    SortHashesWithIndices(Suffixes.data(), HashIndices.data(), suffix_count, 31);

    printf("Building suffix filter...\n");

    // Create using sorted hashes to improve cache hits
    constexpr usize FilterRadix = sizeof(FilterWord) * CHAR_BIT;
    Filter.resize(1 + (UINT32_MAX / FilterRadix));

    for (usize i = 0; i < suffix_count; ++i)
        bit_set(Filter.data(), Suffixes[i]);

    printf("Building suffix buckets...\n");

    HashBuckets.resize(0x1000000);
    SubHashes.resize(Suffixes.size());

    usize here = 0;

    for (usize i = 0; i < HashBuckets.size(); ++i) {
        for (; here < Suffixes.size(); ++here) {
            u32 hash = Suffixes[here];

            if ((hash >> 8) > i)
                break;

            SubHashes[here] = hash & 0xFF;
        }

        HashBuckets[i] = static_cast<u32>(here);
    }

    auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(Stopwatch::now() - start).count();

    printf("Built lookup in %lli ms\n", delta);

    Suffixes.clear();

    printf("Compiled: %zu/%zu (%zu/%zu)\n", PrefixPos, SuffixPos, prefix_count, suffix_count);
}

static void FindMatches(const u32* hashes, usize count, const FilterWord* suffix_filter, void (*callback)(void* context, usize index), void* context)
{
    parallel_partition(count, 0x10000, 0, [=](usize start, usize count) {
        for (usize i = 0; i != count; ++i) {
            usize j = start + i;

            if (u32 hash = hashes[j]; bit_test(suffix_filter, hash))
                callback(context, j);
        }

        return true;
    });
}

usize Collider::Match(HashSet<String>& found)
{
    const u32* hashes = Prefixes[PrefixPos].data();
    usize hash_count = Prefixes[PrefixPos].size();

    // printf("Matching %zu\n", hash_count);

    usize total = 0;

    std::mutex match_lock;
    std::mutex buffer_lock;
    std::mutex found_lock;
    std::condition_variable cond;

    const usize buffer_size = 0x1000;
    auto match_buffer = std::make_unique<std::string[]>(buffer_size);
    auto found_buffer = std::make_unique<std::string[]>(buffer_size);

    std::atomic<usize> head = 0;
    std::atomic<usize> tail = 0;

    auto flush_matches = [&] {
        std::unique_lock buffer_guard { buffer_lock };

        std::unique_lock match_guard { match_lock };

        match_buffer.swap(found_buffer);
        usize count = head.exchange(0);
        tail = 0;

        match_guard.unlock();
        cond.notify_all();

        HashSet<String> sub_found(std::make_move_iterator(&found_buffer[0]), std::make_move_iterator(&found_buffer[count]));
        buffer_guard.unlock();

        std::lock_guard found_guard { found_lock };
        found.merge(std::move(sub_found));
        total += count;
    };

    const auto add_match = [&](String match) {
        usize index = tail.fetch_add(1);

        while (index >= buffer_size) {
            std::unique_lock guard { match_lock };

            index = tail.fetch_add(1);

            if (index < buffer_size)
                break;

            cond.wait(guard);
        }

        match_buffer[index] = std::move(match);

        if (head.fetch_add(1) == buffer_size - 1) {
            flush_matches();
        }
    };

    auto process_match = [&](usize index) {
        u32 hash = hashes[index];

        usize hash_bucket = hash >> 8;
        auto bucket_start = hash_bucket ? HashBuckets[hash_bucket - 1] : 0;
        auto bucket_end = HashBuckets[hash_bucket];

        const u8* subs = SubHashes.data();

        auto find = std::equal_range(&subs[bucket_start], &subs[bucket_end], hash & 0xFF);

        for (; find.first != find.second; ++find.first) {
            String prefix = GetPrefix(index);
            String suffix = GetSuffix(HashIndices[find.first - subs]);
            add_match(prefix + suffix);
        }
    };

    FindMatches(
        hashes, hash_count, Filter.data(), [](void* context, usize index) {
            (*static_cast<decltype(&process_match)>(context))(index);
        },
        &process_match);

    flush_matches();

    return total;
}

usize Collider::Collide(HashSet<String>& found)
{
    printf("Collide %zu/%zu\n", PrefixPos, SuffixPos);

    if (PrefixPos == SuffixPos) {
        usize count = Match(found);
        printf("Matches %zu/%zu\n", count, found.size());
        return count;
    }

    const Vec<String>& sub_parts = Parts[PrefixPos];

    usize total = 0;

    for (usize i = 0; i < sub_parts.size(); ++i) {
        PushPrefix(&sub_parts[i], 1);
        total += Collide(found);
        PopPrefix();
    }

    return total;
}

Vec<String> LoadFile(std::string name)
{
    Vec<String> results;

    std::ifstream input(name);

    for (String line; std::getline(input, line);)
        results.push_back(line);

    return results;
}

Vec<String> LoadFileOrLiteral(const char* path)
{
    if (path[0] == '$')
        return { path + 1 };

    return LoadFile(path);
}

int main(int argc, char** argv)
{
    if (argc < 4) {
        printf("Usage: %s <list_of_hashes.txt> <output_file.txt> <string_parts.txt or $string_literal>\n", argv[0]);
        return 0;
    }

    printf("Loading hashes\n");

    Vec<u32> hashes;

    for (String str : LoadFileOrLiteral(argv[1])) {
        hashes.push_back(joaat_definalize(std::stoul(str, 0, 16)));
    }

    std::ofstream output(argv[2]);

    Vec<Vec<String>> parts;

    for (int i = 3; i < argc; ++i) {
        parts.push_back(LoadFileOrLiteral(argv[i]));
    }

    auto start = Stopwatch::now();

    Collider collider(0, hashes, parts);
    collider.Compile(usize(1) << 28, usize(1) << 31);

    printf("Searching\n");

    HashSet<String> found;
    usize total = collider.Collide(found);

    auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(Stopwatch::now() - start).count();

    printf("Found %zu/%zu results in %lli ms\n", found.size(), total, delta);

    for (String str : found) {
        output << str << "\n";
    }
}
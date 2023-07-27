#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <future>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "parallel.h"

using usize = std::size_t;
using u8 = std::uint8_t;
using u32 = std::uint32_t;

using String = std::string;
using StringView = std::string_view;

template <typename T>
using Vec = std::vector<T>;

template <typename Key, typename Value>
using HashMap = std::unordered_map<Key, Value>;

template <typename Key>
using HashSet = std::unordered_set<Key>;

template <typename T, typename U>
using Pair = std::pair<T, U>;

static u8 to_lower_8(u8 value)
{
    return (value >= 'A' && value <= 'Z') ? (value + ('a' - 'A')) : value;
}

static u32 joaat_partial(u32 hash, const void* data, usize length)
{
    for (usize i = 0; i != length; ++i) {
        hash += to_lower_8(static_cast<const unsigned char*>(data)[i]);
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
        hash -= to_lower_8(static_cast<const unsigned char*>(data)[i - 1]);
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
    Vec<u32> HashStringMapper {};

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
    parallel_partition(count, 0x100000, 0, [=](usize start, usize count) {
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
    parallel_partition(count, 0x100000, 0, [=](usize start, usize count) {
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

    usize suffix_count = Suffixes.size();

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
    if (count < 32) {
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

        if (pivot > 0x100000) {
            auto future = std::async(std::launch::async, sort_lower);
            sort_upper();
            future.wait();
        } else {
            sort_lower();
            sort_upper();
        }
    }
}

#include <chrono>

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

    printf("Compiled: %zu/%zu (%zu/%zu)\n", PrefixPos, SuffixPos, prefix_count, suffix_count);

    printf("Building suffix lookup...\n");

    HashStringMapper.reserve(suffix_count);

    for (usize i = 0; i < suffix_count; ++i)
        HashStringMapper.push_back(static_cast<u32>(i));

    auto start = Stopwatch::now();

    SortHashesWithIndices(Suffixes.data(), HashStringMapper.data(), suffix_count, 31);

    auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(Stopwatch::now() - start).count();

    printf("Built lookup in %lli\n", delta);

    // exit(0);

    printf("Building suffix filter...\n");

    // Create using sorted hashes to improve cache hits

    constexpr usize FilterRadix = sizeof(FilterWord) * CHAR_BIT;
    Filter.resize(1 + (UINT32_MAX / FilterRadix));

    for (usize i = 0; i < suffix_count; ++i)
        bit_set(Filter.data(), Suffixes[i]);

    printf("Compiled\n");
}

static void FindMatches(const u32* hashes, usize count, const FilterWord* suffix_filter, void (*callback)(void* context, usize index), void* context)
{
    parallel_partition(count, 0x1000000, 0, [=](usize start, usize count) {
        for (usize i = 0; i != count; ++i) {
            usize j = start + i;

            if (u32 hash = hashes[j]; bit_test(suffix_filter, hash))
                callback(context, j);
        }

        return true;
    });
}

#include <shared_mutex>

usize Collider::Match(HashSet<String>& found)
{
    const u32* hashes = Prefixes[PrefixPos].data();
    usize hash_count = Prefixes[PrefixPos].size();

    printf("Matching %zu\n", hash_count);

    std::shared_mutex lock;
    usize total = 0;

    auto process_match = [&](usize index) {
        u32 hash = hashes[index];

        auto find = std::equal_range(Suffixes.begin(), Suffixes.end(), hash);

        for (; find.first != find.second; ++find.first) {
            String prefix = GetPrefix(index);
            String suffix = GetSuffix(HashStringMapper[find.first - Suffixes.begin()]);

            String match = prefix + suffix;

            std::lock_guard guard { lock };
            found.emplace(std::move(match));
            ++total;
        }
    };

    FindMatches(
        hashes, hash_count, Filter.data(), [](void* context, usize index) {
            (*static_cast<decltype(&process_match)>(context))(index);
        },
        &process_match);

    return total;
}

usize Collider::Collide(HashSet<String>& found)
{
    printf("Collide %zu/%zu %zu/%zu\n", PrefixPos, SuffixPos, Prefixes[PrefixPos].size(), Suffixes.size());

    if (PrefixPos == SuffixPos) {
        return Match(found);
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

#include <fstream>

Vec<String> LoadFile(std::string name)
{
    Vec<String> results;

    std::ifstream input(name);

    for (String line; std::getline(input, line);)
        results.push_back(line);

    return results;
}

int main()
{
    String BASE_PATH = R"(W:\Projects\joaat-in-the-middle\joaat-in-the-middle\strings\)";
    Vec<Vec<String>> PARTS;
    Vec<u32> TARGET_HASHES;

    Vec<String> TARGET_STRINGS = LoadFile(BASE_PATH + "cam_targets.txt");

    for (StringView string : TARGET_STRINGS) {
        TARGET_HASHES.push_back(joaat_partial(0, string.data(), string.size()));
    }

    PARTS.push_back({ "cam" });
    PARTS.push_back({ "Cinematic" });

    Vec<String> MIDDLES = LoadFile(BASE_PATH + "cam_middles.txt");

    PARTS.push_back(MIDDLES);
    PARTS.push_back(MIDDLES); // 11/12
    PARTS.push_back(MIDDLES); // 35/57
    PARTS.push_back(MIDDLES); // 57/182
    PARTS.push_back(MIDDLES); // 99/497
    PARTS.push_back(MIDDLES); // 1606/2725
    PARTS.push_back(MIDDLES); // 106236/117888
    PARTS.push_back(MIDDLES); // 7238934/8017458
    // PARTS.push_back(MIDDLES);

    PARTS.push_back({ "Metadata" });

    Collider collider(0, TARGET_HASHES, PARTS);

    collider.Compile(usize(1) << 32, usize(1) << 30);

    HashSet<String> found;
    usize total = collider.Collide(found);

    printf("Found %zu/%zu\n", found.size(), total);

    for (String str : TARGET_STRINGS) {
        if (found.find(str) == found.end())
            printf("Didn't find %s\n", str.c_str());
    }

    for (String str : found) {
        if (std::find(TARGET_HASHES.begin(), TARGET_HASHES.end(), joaat_partial(0, str.data(), str.size())) == TARGET_HASHES.end())
            printf("Invalid String: %s\n", str.c_str());
    }

    printf("Cleanup...\n");
}
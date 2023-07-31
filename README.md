# Joaat Meet-In-The-Middle Attack

An implementation of a meet-in-the-middle attack on the jenkins-one-at-a-time hashing algorithm. For some more information, see [my blog post](https://0x1f9f1.github.io/2023/07/29/hash-mitm.html).

The code is focused more on using fast algorithms, rather than micro-optimisations. This includes:
* Using a meet-in-the-middle attack.
* Sorting the pre-computed hashes using a multi-threaded in-place hybrid MSD-radix/insertion sort, while maintaining a mapping between the sorted hashes and their original position.
* Reducting memory usage by compressing hashes into buckets.
* Using a 2<sup>32</sup>-bit bitset to check for valid hashes (this could have used a bloom filter, but as the number of hashes increases, the target size of a bloom filter becomes larger than 2<sup>32</sup> bits anyway).
* SIMD hashing
* Multi-threaded hashing and matching
* Using the original position of the hash to encode its string value (avoiding the need to store the string representation of each hash).
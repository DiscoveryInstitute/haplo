#pragma once

#include <cstddef>
#include <functional>
#include <random>
#include <string>
#include <utility>
#include <vector>

class BlocksRange {

	template<typename V> friend class Blocks;
    friend class FileCheckpointer;
    friend class FileMemorySaver;
    friend class HaploBlockBoundaries;

  private:

    std::vector<long> m_toggles;

    BlocksRange(std::vector<long>&& toggles) : m_toggles(toggles) {}
    BlocksRange combine (const BlocksRange &range_A, const BlocksRange &range_B,
            std::function<bool (bool, bool)> rule) const;

  public:

    /** Default constructor */
    BlocksRange() {}
    /** Create BlocksRange where the whole range is empty/false. */
    static BlocksRange createEmpty() { return BlocksRange(std::vector<long>{}); }
    /** Create BlocksRange where the whole range is full/true. */
    static BlocksRange createFull() { return BlocksRange(std::vector<long>{0}); }
    /** Create a BlocksRange with a single flip at the specified location. */
    static BlocksRange createOneToggle(long b) { return BlocksRange(std::vector<long>{b}); }
    /** Create a random BlocksRange which randomly starts as true or false,
     *  and then has a number of flips distributed randomly at the given rate. */
    static BlocksRange createRandomToggles(long blocks, double rate, std::mt19937_64& rng);

    /** Is the whole range empty/false? */
    bool empty() const { return m_toggles.empty(); }
    /** Is the whole range full/true? */
    bool full() const { return m_toggles.size() == 1 && m_toggles[0] == 0; }

    /** Number of blocks included in this range.
     *  Must also pass max_range because BlocksRange does not end of range (for efficiency). */
    long getSizeOfRange(long max_range) const;

    /** Return new BlocksRange, where every subrange that was true is now false, and vice versa. */
    BlocksRange inverse() const;
    /** Return new BlocksRange that is union of 'this' and 'other'. */
    BlocksRange unionWith(const BlocksRange& other) const;
    /** Return new BlocksRange that is union of 'this' and the inverse of 'other'. */
    BlocksRange unionWithInverse(const BlocksRange& other) const;
    /** Return new BlocksRange that is intersection of 'this' and 'other'. */
    BlocksRange intersectionWith(const BlocksRange& other) const;
    /** Return new BlocksRange that is intersection of 'this' and the inverse of 'other'. */
    BlocksRange intersectionWithInverse(const BlocksRange& other) const;
    /** Return new BlocksRange that is symmetric-difference (XOR) between 'this' and 'other'. */
    BlocksRange symmetricDifference(const BlocksRange& other) const;

	bool operator==(const BlocksRange&) const;
    std::string toString() const;
    size_t getSizeInMemory() const;
};

#pragma once

#include <cstddef>
#include <functional>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include "BlocksRange.h"

template<typename V>
class Blocks {

    friend class FileCheckpointer;
    friend class FileMemorySaver;
    friend class HaploBlockVariants;

  private:

    std::vector<std::pair<long, V>> m_intervals;

  public:

    /** Default constructor */
    Blocks<V>() {}

    /** Create Blocks where the whole range is one value. */
    Blocks<V>(V value) : m_intervals({std::pair<long,V>(0,value)}) {}

    /** Combine two Blocks according to 'combinePattern'. */
    static Blocks<V> combine(const Blocks<V>& a, const Blocks<V>& b, const BlocksRange& combinePattern);

    /** Filter Blocks according to 'filterPattern'. */
    Blocks<V> filtered(const BlocksRange& filterPattern) const;

    Blocks<V> compressed();

    bool operator==(const Blocks<V>& other) const;

    V getValue(long position) const;

    std::string toString() const;
    size_t getSizeInMemory() const;

};

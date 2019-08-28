#include <algorithm>
#include <string>
#include "Blocks.h"
#include "BlocksRange.h"

using namespace std;

template<typename V>
Blocks<V> Blocks<V>::combine(
    const Blocks<V>& blocks_A, const Blocks<V>& blocks_B, const BlocksRange& pattern) {

    auto it0 = blocks_A.m_intervals.cbegin();   // 0 for off/unselected
    auto it1  = blocks_B.m_intervals.cbegin();  // 1 for on/selected
    auto end0 = blocks_A.m_intervals.cend();
    auto end1  = blocks_B.m_intervals.cend();
    auto tog_it = pattern.m_toggles.cbegin();
    auto tog_end = pattern.m_toggles.cend();
    V current0 = V();
    V current1 = V();

	Blocks<V> blocks_C;
    auto& intervals_C = blocks_C.m_intervals;

    for ( ; tog_it != tog_end; ++tog_it) {

        for ( ; it1 != end1 && it1->first < *tog_it; ++it1) {
            intervals_C.push_back(*it1);
            current1 = it1->second;
        }

        swap(it0, it1);
        swap(end0, end1);
        swap(current0, current1);

        for ( ; it1 != end1 && it1->first <= *tog_it; ++it1) {
            current1 = it1->second;
        }
        if (current1 != current0) {
            intervals_C.emplace_back(*tog_it, current1);
        }

    }
    for ( ; it1 != end1; ++it1) {
        intervals_C.push_back(*it1);
        current1 = it1->second;
    }

    return blocks_C;
}

template<typename V>
Blocks<V> Blocks<V>::filtered(const BlocksRange& filterPattern) const {
    return combine(*this, Blocks<V>(), filterPattern);
}

template<typename V>
Blocks<V> Blocks<V>::compressed() {
    const auto& old_intervals = this->m_intervals;
    Blocks<V> new_blocks;
    auto& new_intervals = new_blocks.m_intervals;
    auto it = old_intervals.begin();
    auto end_it = old_intervals.end();
    new_intervals.emplace_back(*it);
    for ( ; it != end_it; ++it) {
        if (new_intervals.back().second != it->second) {
            new_intervals.emplace_back(*it);
        }
    }
    return new_blocks;
}

template<typename V>
bool Blocks<V>::operator==(const Blocks<V>& other) const {
    return (this->m_intervals.size() == other.m_intervals.size()) &&
            equal(this->m_intervals.cbegin(), this->m_intervals.cend(), other.m_intervals.cbegin());
}



template<typename V>
V Blocks<V>::getValue(long position) const {
    auto itBegin = m_intervals.cbegin();
    auto itEnd = m_intervals.cend();
    typedef pair<long,V> P;
    auto it = upper_bound(
                itBegin, itEnd, make_pair(position, V()),
                [](const P& a, const P& b) { return a.first < b.first; });
    return it==itBegin ? V() : prev(it)->second;
}

template<typename V>
string Blocks<V>::toString() const {
	string s = "";
    const auto& v = m_intervals;
	long n = v.size();
    if (n==0) return "[0,end): " + to_string(V()) + ", ";
	long last = n-1;
    if (v[0].first != 0) {
        s += "[0," + to_string(v[0].first) + "):" + to_string(V()) + ", ";
    }
	for(int i=0; i<last; i++) {
		s += "[" + to_string(v[i].first) + "," + to_string(v[i+1].first) + "):" + to_string(v[i].second) + ", ";
	}
	s += "[" + to_string(v[last].first) + ",end):" + to_string(v[last].second) + ", ";
	return s;
}

template<typename V>
size_t Blocks<V>::getSizeInMemory() const {
    return sizeof(vector<pair<long,V>>) + m_intervals.size() * sizeof(pair<long,V>);
}



// Specify which class(es) we want to pre-compile.
template class Blocks<int>;
template class Blocks<long>;

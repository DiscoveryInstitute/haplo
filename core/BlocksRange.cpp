#include <algorithm>
#include <functional>
#include <iterator>
#include <random>
#include <string>
#include <vector>
#include "BlocksRange.h"

using namespace std;

long BlocksRange::getSizeOfRange(long max_blocks) const {
    long sum = 0L;
    int n = m_toggles.size();
    for (int i=1; i<n; i+=2) {
        sum += m_toggles[i] - m_toggles[i-1];
    }
    if (n % 2 == 1) {
        sum += max_blocks - m_toggles[n-1];
    }
    return sum;
}

BlocksRange BlocksRange::createRandomToggles(
        long total_blocks, double toggle_rate, mt19937_64& rng) {

    // Generate how many random toggles there will be.
    int n_random = poisson_distribution<int>(toggle_rate)(rng);

    // Generate the toggles
    uniform_int_distribution<long> random(0, total_blocks-1);
    long toggles[n_random+1];
    int it = 0;
    if (2 * random(rng) >= total_blocks) toggles[it++] = 0;       // Start true or false.
    for (int i=0; i<n_random; i++) toggles[it++] = random(rng);   // Add random toggles.

    // Sort toggles into order
    sort(toggles, toggles+it);

    // Remove duplicates: if there are n toggles the same value, then n <- n % 2
    int prev = -1;
    int n = it;
    it = 0;
    for (int i=0; i<n; i++) {
        if (toggles[i] == prev) {
            --it; prev = -1;
        } else {
            prev = toggles[it++] = toggles[i];
        }
    }

    // Copy the results into a new vector and return new BlocksRange.
    return BlocksRange(vector<long>(toggles,toggles+it));
}


BlocksRange BlocksRange::inverse() const {
    int n = m_toggles.size();
    vector<long> new_toggles;
    // If toggles start with 0, remove it. If not, add it.
    if (!m_toggles.empty() && m_toggles[0] == 0) {
        new_toggles.reserve(n-1);
        for (int i=1; i<n; i++) new_toggles.push_back(m_toggles[i]);
        return BlocksRange(move(new_toggles));
    } else {
        new_toggles.reserve(n+1);
        new_toggles.push_back(0);
        for (int i=0; i<n; i++) new_toggles.push_back(m_toggles[i]);
        return BlocksRange(move(new_toggles));
    }
}

void toggle(bool& value) { value = !value; }

BlocksRange BlocksRange::combine (
        const BlocksRange &range_A, const BlocksRange &range_B,
        function<bool (bool, bool)> rule) const
{

    auto iterator_A = range_A.m_toggles.cbegin();
    auto iterator_B = range_B.m_toggles.cbegin();
    auto end_A = range_A.m_toggles.cend();
    auto end_B = range_B.m_toggles.cend();

    BlocksRange range_C;
    auto& toggle_positions_C = range_C.m_toggles;

    long block = 0;
    bool toggle_state_A = false;
    bool toggle_state_B = false;
    bool toggle_state_C = false;

    if (iterator_A != end_A && *iterator_A == block) { iterator_A++; toggle(toggle_state_A); }
    if (iterator_B != end_B && *iterator_B == block) { iterator_B++; toggle(toggle_state_B); }
    // Check if we need to add a toggle to C at block 0.
    if (toggle_state_C != rule(toggle_state_A, toggle_state_B)) {
        toggle_positions_C.push_back(block);
        toggle(toggle_state_C);
    }
    while (iterator_A != end_A && iterator_B != end_B) {
        // Advance A and/or B to next block where there is a toggle.
        long next_block_A = *iterator_A;
        long next_block_B = *iterator_B;
        if (next_block_A <= next_block_B) { block = *iterator_A++; toggle(toggle_state_A); }
        if (next_block_B <= next_block_A) { block = *iterator_B++; toggle(toggle_state_B); }
        // Check if we need to add a toggle to C.
        if (toggle_state_C != rule(toggle_state_A, toggle_state_B)) {
            toggle_positions_C.push_back(block);
            toggle(toggle_state_C);
        }
    }
    while (iterator_A != end_A) {
        // Advance A to next block where there is a toggle.
        block = *iterator_A++; toggle(toggle_state_A);
        // Check if we need to add a toggle to C.
        if (toggle_state_C != rule(toggle_state_A, toggle_state_B)) {
            toggle_positions_C.push_back(block);
            toggle(toggle_state_C);
        }
    }
    while (iterator_B != end_B) {
        // Advance B to next block where there is a toggle.
        block = *iterator_B++; toggle(toggle_state_B);
        // Check if we need to add a toggle to C.
        if (toggle_state_C != rule(toggle_state_A, toggle_state_B)) {
            toggle_positions_C.push_back(block);
            toggle(toggle_state_C);
        }
    }

    return range_C;
}

BlocksRange BlocksRange::unionWith(const BlocksRange& other) const {
    return combine(*this, other, [](bool a, bool b){ return a || b; });
}
BlocksRange BlocksRange::unionWithInverse(const BlocksRange& other) const {
    return combine(*this, other, [](bool a, bool b){ return a || !b; });
}
BlocksRange BlocksRange::intersectionWith(const BlocksRange& other) const {
    return combine(*this, other, [](bool a, bool b){ return a && b; });
}
BlocksRange BlocksRange::intersectionWithInverse(const BlocksRange& other) const {
    return combine(*this, other, [](bool a, bool b){ return a && !b; });
}
BlocksRange BlocksRange::symmetricDifference(const BlocksRange& other) const {
    return combine(*this, other, [](bool a, bool b){ return !a != !b; });
}


bool BlocksRange::operator==(const BlocksRange& other) const {
	return equal(this->m_toggles.cbegin(), this->m_toggles.cend(), other.m_toggles.cbegin());
}


string BlocksRange::toString() const {
	if (m_toggles.empty()) return "[empty]";
    string s = "";
	int last = m_toggles.size() - 1;
	for (int i=0; i<last; i+=2) {
		s += "[" + to_string(m_toggles[i]) + "," + to_string(m_toggles[i+1]) + ") ";
	}
	if (last%2==0) {
		s += "[" + to_string(m_toggles[last]) + ",end)";
    }
    return s;
}

size_t BlocksRange::getSizeInMemory() const {
    return sizeof(vector<long>) + m_toggles.size() * sizeof(long);
}


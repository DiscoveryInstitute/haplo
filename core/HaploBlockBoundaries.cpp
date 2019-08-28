#include "HaploBlockBoundaries.h"
#include <algorithm>
#include "BlocksRange.h"

using namespace std;


HaploBlockBoundaries::HaploBlockBoundaries(long chromosome_length, long maximum_blocks, mt19937_64& rng) {

    set<LocusId> boundary_set {0, chromosome_length};
    auto nblocks = [&]() -> long { return boundary_set.size() - 1; };

    const long first=1L;
    const long last=chromosome_length-1L;
    uniform_int_distribution<LocusId> randomLocus(first, last); // exclude first and last
    if (maximum_blocks > chromosome_length) maximum_blocks = chromosome_length;
    if (maximum_blocks < chromosome_length/2) {
        while (nblocks() < maximum_blocks) {
            boundary_set.insert(randomLocus(rng));
        }
    } else {
        for(LocusId loc=first; loc<=last; ++loc) boundary_set.insert(loc);
        while (nblocks() > maximum_blocks) {
            boundary_set.erase(randomLocus(rng));
        }
    }

    m_block_boundaries = vector<LocusId>(boundary_set.begin(), boundary_set.end());

    for (LocusId locus : m_block_boundaries) m_block_boundary_used[locus] = 0L;
    m_block_boundary_used[0] = 1L;
    m_block_boundary_used[chromosome_length] = 1L;

}

BlocksRange HaploBlockBoundaries::createRandomRecombinationPattern(
        double nucleotide_rate, mt19937_64& rng) const {

    // If only one block (2 boundaries), there is no recombination.
    if (m_block_boundaries.size()==2) return BlocksRange();

    // Generate how many random toggles there will be.
    long chromosome_length = m_block_boundaries.back() - 1;
    double total_rate = nucleotide_rate * (chromosome_length - 1);
    int n_random = poisson_distribution<int>(total_rate)(rng);
    vector<long> toggles; toggles.reserve(n_random+1);

    // Start true or false.
    if (bernoulli_distribution(0.5)(rng)) toggles.push_back(0);

    // Generate the toggles
    uniform_int_distribution<BlockId> randomBlock(1, m_block_boundaries.size()-2); // exclude first and last
    for (int i=0; i<n_random; i++) {
        toggles.push_back(randomBlock(rng));   // Add random toggles.
    }
    sort(toggles.begin(), toggles.end());

    // Remove duplicate pairs (they cancel each other)
    size_t i=0, j=0; // start at second element
    for( ; i<toggles.size(); ++i, ++j) {
        toggles[j] = toggles[i];
        if (j > 0 && toggles[j] == toggles[j-1]) j -= 2;
    }
    toggles.resize(j);

    for (auto& toggle : toggles) toggle = m_block_boundaries[toggle];

    // Move the result into a new BlocksRange.
    return toggles.empty() ? BlocksRange() : BlocksRange(move(toggles));
}

void HaploBlockBoundaries::recordBlockBoundaries(const BlocksRange& pattern) {
    for(LocusId locus : pattern.m_toggles) m_block_boundary_used[locus]++;
}

vector<LocusId> HaploBlockBoundaries::getRecordedBlockBoundaries() const {
    vector<LocusId> boundaries;
    for (const auto& entry : m_block_boundary_used) if (entry.second > 0) boundaries.push_back(entry.first);
    return boundaries;
}

vector<LocusId> HaploBlockBoundaries::getRandomBlockBoundaries(long nblocks, mt19937_64& rng) const {
    const LocusId maximum_blocks = m_block_boundaries.size()-1;
    set<long> block_set {0, maximum_blocks};
    const long first = 1L;
    const long last = maximum_blocks - 1;
    uniform_int_distribution<BlockId> randomBlock(first, last); // exclude first and last
    auto out_nblocks = [&]() -> long { return block_set.size()-1; };
    if (nblocks < maximum_blocks/2) {
        while (out_nblocks() < nblocks) {
            block_set.insert(randomBlock(rng));
        }
    } else {
        for(LocusId loc=first; loc<=last; ++loc) block_set.insert(loc);
        while (out_nblocks() > nblocks) {
            block_set.erase(randomBlock(rng));
        }
    }

    vector<LocusId> block_vec; block_vec.reserve(block_set.size());
    for (const auto& blockid : block_set) { block_vec.push_back(m_block_boundaries[blockid]); }
    return block_vec;

}

string HaploBlockBoundaries::toString() const {
	string s = "";
	for (auto& b : m_block_boundaries) {
		s += to_string(b) + " ";
	}
	return s;
}

bool HaploBlockBoundaries::operator==(const HaploBlockBoundaries& o) const {
    return m_block_boundaries == o.m_block_boundaries
        && m_block_boundary_used == o.m_block_boundary_used;
}


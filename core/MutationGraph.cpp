#include <iterator>
#include <map>
#include <vector>
#include "Chromosome.h"
#include "MutationGraph.h"

using namespace std;

MutationGraph::MutationGraph(vector<vector<LocusId>> founder_alleles, const bool use_mutation_loci)
        : m_use_mutation_loci(use_mutation_loci)
{
    int n_founder_alleles = founder_alleles.size();
    m_next_id = 1;
    for (AlleleId i = 0; i < n_founder_alleles; ++i) {
        vector<LocusId>& loci = founder_alleles.at(i);
        long n_mutations = loci.size();
        vector<pair<BlockId, long>> mutations;
        if (m_use_mutation_loci) {
            mutations.reserve(loci.size());
            for (auto locus : loci) mutations.emplace_back(locus, UNKNOWN_GENERATION);
        }
        m_nodes[m_next_id++] = Node { UNKNOWN_ALLELE, n_mutations, move(mutations), 0, 0 };
    }
}


AlleleId MutationGraph::addAlleleMutation(AlleleId parent_id, long gen, const vector<LocusId>& loci) {
    if (parent_id == UNKNOWN_ALLELE) throw logic_error("Cannot mutate UNKNOWN_ALLELE");

    long n_mutations = loci.size();
    vector<pair<BlockId, long>> mutations;
    if (m_use_mutation_loci) {
        mutations.reserve(loci.size());
        for (auto locus : loci) mutations.emplace_back(locus, gen);
    }
    m_nodes.at(parent_id).n_child_nodes++;

    AlleleId id = m_next_id++;
    m_nodes[id] = Node { parent_id, n_mutations, move(mutations), 0, 0 };
    return id;
}


AlleleId MutationGraph::getFounderId(AlleleId id) const {
    while (m_nodes.at(id).parent_id) {
        id = m_nodes.at(id).parent_id;
    }
    return id;
}

bool MutationGraph::isDerivativeAllele(AlleleId parent_id, AlleleId id) const {
    while (true) {
        if (id == parent_id) return true;
        if (id == UNKNOWN_ALLELE) return false;
        id = m_nodes.at(id).parent_id;
    }
}


int MutationGraph::getMutationDistance(AlleleId id_a, AlleleId id_b) const {
    if (id_a == UNKNOWN_ALLELE || id_b == UNKNOWN_ALLELE)
        throw logic_error("Cannot get mutation distance of with UNKNOWN_ALLELE");

    const Node* node_a;
    const Node* node_b;
    int distance_a, distance_b;
    int overlap_a, overlap_b;

    node_a = &m_nodes.at(id_a); distance_a = 0;
    node_b = &m_nodes.at(id_b); distance_b = 0;

    // Calculate distance to top
    while (node_a->parent_id) {
        distance_a += node_a->n_mutations;
        node_a = &m_nodes.at(node_a->parent_id);
    }
    while (node_b->parent_id) {
        distance_b += node_b->n_mutations;
        node_b = &m_nodes.at(node_b->parent_id);
    }
    if (node_a != node_b) {
        distance_a += node_a->n_mutations;
        distance_b += node_b->n_mutations;
        return distance_a + distance_b;
    }

    // Calculate distance of overlap
    node_a = &m_nodes.at(id_a); overlap_a = distance_a;
    node_b = &m_nodes.at(id_b); overlap_b = distance_b;
    while (node_a != node_b) {
        bool up_a = overlap_a >= overlap_b;
        bool up_b = overlap_a <= overlap_b;
        if (up_a) {
            overlap_a -= node_a->n_mutations;
            node_a = &m_nodes.at(node_a->parent_id);
        }
        if (up_b) {
            overlap_b -= node_b->n_mutations;
            node_b = &m_nodes.at(node_b->parent_id);
        }
    }

    return distance_a + distance_b - overlap_a - overlap_b;
}


void MutationGraph::simplify() {
    // Remove all nodes on dead branches.
    vector<AlleleId> dead_branch_ids;
    for (auto it = m_nodes.rbegin(); it != m_nodes.rend(); ++it) {
        auto* node = &it->second;
        if (node->n_child_nodes == 0 && node->frequency == 0) {
            AlleleId parent_id = node->parent_id;
            dead_branch_ids.push_back(it->first);
            if (parent_id) {
                auto* parent_node = &m_nodes.find(parent_id)->second;
                --parent_node->n_child_nodes;
            }
        }
    }
    for (AlleleId id : dead_branch_ids) m_nodes.erase(id);

    // Compress other unforked branches by removing empty nodes.
    vector<AlleleId> unforked_branch_empty_node_ids;
    string s = toString();
    for (auto it = m_nodes.rbegin(); it != m_nodes.rend(); ++it) {
        auto* node = &it->second;
        while(node->parent_id) {
            auto parent_it = m_nodes.find(node->parent_id);
            auto* parent_node = &parent_it->second;
            if (parent_node->n_child_nodes==1 && parent_node->frequency==0) {
                if (parent_node->parent_id == UNKNOWN_ALLELE) break; // leave founder nodes alone
                node->parent_id = parent_node->parent_id;
                node->n_mutations += parent_node->n_mutations;
                if (m_use_mutation_loci) {
                    copy(parent_node->mutations.begin(), parent_node->mutations.end(),
                         back_inserter(node->mutations));
                }
                unforked_branch_empty_node_ids.push_back(parent_it->first);
            } else {
                break;
            }
        }
    }
    for (AlleleId id : unforked_branch_empty_node_ids) m_nodes.erase(id);

    // Assign to a copy of itself (hopefully aligns the memory)
    auto copy_nodes = m_nodes;
    m_nodes = move(copy_nodes);

    // CHECK !!!
    for (auto& entry : m_nodes) {
        auto id = entry.first;
        auto* node = &entry.second;
        if (node->parent_id != UNKNOWN_ALLELE) {
            auto parent_it = m_nodes.find(node->parent_id);
            if (parent_it==m_nodes.end() || id <= node->parent_id || node->parent_id <0) {
                throw logic_error("invalid parent_id");
            }
        }
        if (node->n_child_nodes == 0 && node->frequency == 0) {
            throw logic_error("dead nodes still present");
        }
    }
}


vector<AlleleId> MutationGraph::getAlleleIds() const {
    vector<AlleleId> ret;
    for (auto& entry : m_nodes) ret.push_back(entry.first);
    return ret;
}

string MutationGraph::toString() const {

    string ret;
    ret += to_string(m_nodes.size()) + " ";
    for (auto entry : m_nodes) {
        auto& id = entry.first;
        auto& node = entry.second;
        long n_mutations = node.n_mutations;
        auto& mutations = node.mutations;
        long generation = (mutations.empty()) ? UNKNOWN_GENERATION : mutations.back().second;
        ret += "(" + to_string(id) + " "
                + to_string(getParentId(id)) + " "
                + to_string(getFounderId(id)) + " "
                + to_string(n_mutations) + " "
                + to_string(generation) + " "
                + to_string(node.frequency) + ") ";
    }
    return ret;
}

bool MutationGraph::Node::operator==(const Node& o) const {
    return parent_id     == o.parent_id
        && n_mutations   == o.n_mutations
        && mutations     == o.mutations
        && n_child_nodes == o.n_child_nodes
        && frequency     == o.frequency;
}

bool MutationGraph::operator==(const MutationGraph& o) const {
    return m_use_mutation_loci == o.m_use_mutation_loci
        && m_nodes             == o.m_nodes
        && m_next_id           == o.m_next_id;
}

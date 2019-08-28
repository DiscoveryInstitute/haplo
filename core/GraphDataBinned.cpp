#include "GraphDataBinned.h"

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

using namespace std;

string GraphDataBinned::toString(bool bar_end) const {
    string s;
    size_t n = m_weight_bins.size();
    double delta = (m_max - m_min) / n;
    for(size_t i = 0; i < n; i++) {
        double bin_min = m_min + i * delta;
        double bin_max = i+1==n ? m_max : bin_min + delta;
        double bin_weight = m_weight_bins.at(i);
        double bin_value = m_value_bins.at(i);
        double bin_average = bin_weight ? bin_value / bin_weight : 0.0;
        s += to_string(bin_min); s += " ";
        if (bar_end) { s += to_string(bin_max); s += " "; }
        s += to_string(bin_weight); s += " ";
        s += to_string(bin_average); s += "\n";
    }
    return s;
}

void assertBinsIdentical(const GraphDataBinned& a, const GraphDataBinned& b) {
    if (a.m_weight_bins.size() != b.m_weight_bins.size() ||
        a.m_value_bins.size() != b.m_value_bins.size() ||
        a.m_min != b.m_min ||
        a.m_max != b.m_max ||
        a.m_delta_inv != b.m_delta_inv ||
        a.m_bin_includes_upper_bound != b.m_bin_includes_upper_bound)
        throw invalid_argument("Attempt to combine data with different bins");
}

GraphDataBinned& GraphDataBinned::merge(const GraphDataBinned& other) {
    assertBinsIdentical(*this, other);
    size_t n = m_weight_bins.size();
    for (size_t i = 0; i < n; ++i) {
        this->m_weight_bins[i] += other.m_weight_bins[i];
        this->m_value_bins[i] += other.m_value_bins[i];
    }
    return *this;
}

GraphDataBinned& GraphDataBinned::divide(const GraphDataBinned& divisor) {
    assertBinsIdentical(*this, divisor);
    size_t n = m_weight_bins.size();
    for (size_t i = 0; i < n; ++i) {
        double weight = divisor.m_weight_bins[i];
        double value = divisor.m_value_bins[i];
        if (weight == 0.0) continue;
        if (value == 0.0) this->m_value_bins[i] = numeric_limits<double>::max();
        else this->m_value_bins[i] *= weight / value;
    }
    return *this;

}

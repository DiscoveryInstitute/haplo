#include "GraphDataBinned2D.h"

#include <cstddef>
#include <stdexcept>
#include <string>

using namespace std;

string GraphDataBinned2D::toString(bool bar_end) const {
    string s;
    size_t n = m_weight_bins.size();
    double delta_x = (m_range_x.second - m_range_x.first) / n;
    double delta_y = (m_range_y.second - m_range_y.first) / n;
    for(size_t ix = 0; ix < n; ++ix) {
        double bin_min_x = m_range_x.first + ix * delta_x;
        double bin_max_x = ix+1==n ? m_range_x.second : bin_min_x + delta_x;
        for(size_t iy = 0; iy < n; ++iy) {
            double bin_min_y = m_range_y.first + iy * delta_y;
            double bin_max_y = iy+1==n ? m_range_y.second : bin_min_y + delta_y;

            double bin_weight = m_weight_bins.at(ix).at(iy);
            double bin_value = m_value_bins.at(ix).at(iy);
            double bin_average = bin_weight ? bin_value / bin_weight : 0.0;
            s += to_string(bin_min_x) + " " + to_string(bin_min_y) + " ";
            if (bar_end) s += to_string(bin_max_x) + " " + to_string(bin_max_y) + " ";
            s +=  to_string(bin_weight) + " " + to_string(bin_average) + "\n";
        }
    }
    return s;
}

GraphDataBinned2D& GraphDataBinned2D::merge(GraphDataBinned2D& other) {
    if (this->m_weight_bins.size() != other.m_weight_bins.size() ||
        this->m_value_bins.size() != other.m_value_bins.size() ||
        this->m_range_x != other.m_range_x ||
        this->m_range_y != other.m_range_y ||
        this->m_delta_inv_x != other.m_delta_inv_x ||
        this->m_delta_inv_y != other.m_delta_inv_y ||
        this->m_bin_includes_upper_bound != other.m_bin_includes_upper_bound)
        throw invalid_argument("Attempt to merge data with different bins");

    size_t n = m_weight_bins.size();
    for (size_t ix = 0; ix < n; ++ix) {
        for (size_t iy = 0; iy < n; ++iy) {
            this->m_weight_bins[ix][iy] += other.m_weight_bins[ix][iy];
            this->m_value_bins[ix][iy] += other.m_value_bins[ix][iy];
        }
    }
    return *this;
}


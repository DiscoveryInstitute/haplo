#pragma once

#include <cmath>
#include <string>
#include <vector>

class GraphDataBinned2D {

    public:

        std::vector<std::vector<double>> m_weight_bins;
        std::vector<std::vector<double>> m_value_bins;
        std::pair<double,double> m_range_x;
        std::pair<double,double> m_range_y;
        double m_delta_inv_x;
        double m_delta_inv_y;
        bool m_bin_includes_upper_bound;

        /** Create a data-structure for summarising data in bins for graphs. */
        GraphDataBinned2D(std::pair<double,double> range_x, std::pair<double,double> range_y,
                          std::size_t nbins, bool bin_includes_upper_bound = false)
            : m_weight_bins(nbins, std::vector<double>(nbins)),
              m_value_bins(nbins, std::vector<double>(nbins)),
              m_range_x(range_x),
              m_range_y(range_y),
              m_delta_inv_x(nbins/(range_x.second-range_x.first)),
              m_delta_inv_y(nbins/(range_y.second-range_y.first)),
              m_bin_includes_upper_bound(bin_includes_upper_bound) {}

        // For making copies, but with data zeroed, for use in omp parallel for loops.
        GraphDataBinned2D zeroedCopy() {
            return GraphDataBinned2D(m_range_x, m_range_y, m_value_bins.size(), m_bin_includes_upper_bound);
        }

        /** Sample the value at a particular position. */
        inline void sampleValue(double x, double y, double value, double multiplier = 1.0) {
            const long n = (long)m_weight_bins.size();
            double ifltx = m_delta_inv_x * (x - m_range_x.first);
            double iflty = m_delta_inv_y * (y - m_range_y.first);
            long ix = m_bin_includes_upper_bound ? ceil(ifltx) - 1 : ifltx;
            long iy = m_bin_includes_upper_bound ? ceil(iflty) - 1 : iflty;
            if (ix < 0 || ix >= n) return;
            if (iy < 0 || iy >= n) return;
            m_weight_bins.at(ix).at(iy) += multiplier;
            m_value_bins.at(ix).at(iy)  += value * multiplier;
        }

        /** Sample weight at a particular position. */
        inline void sampleWeight(double x, double y, double multiplier = 1.0) {
            sampleValue(x, y, 0.0, multiplier);
        }

        std::string toString(bool endbar=false) const;

        // For reduction in omp parallel for loops.
        GraphDataBinned2D& merge(GraphDataBinned2D& other);
};

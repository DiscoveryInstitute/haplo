#pragma once

#include <cmath>
#include <string>
#include <vector>

class GraphDataBinned {

    public:

        std::vector<double> m_weight_bins;
        std::vector<double> m_value_bins;
        double m_min;
        double m_max;
        double m_delta_inv;
        bool m_bin_includes_upper_bound;

        /** Create a data-structure for summarising data in bins for graphs. */
        GraphDataBinned(double mn, double mx, std::size_t nbins, bool bin_includes_upper_bound = false)
            : m_weight_bins(nbins), m_value_bins(nbins),
              m_min(mn), m_max(mx), m_delta_inv(nbins / (mx-mn)),
              m_bin_includes_upper_bound(bin_includes_upper_bound) {}

        // For making copies, but with data zeroed, for use in omp parallel for loops.
        GraphDataBinned zeroedCopy() {
            return GraphDataBinned(m_min, m_max, m_value_bins.size(),m_bin_includes_upper_bound);
        }

        /** Sample the value (y-position) at a particular x-position (for line graphs). */
        inline void sampleValue(double x, double y, double multiplier = 1.0) {
            double iflt = m_delta_inv * (x - m_min);
            long i = m_bin_includes_upper_bound ? ceil(iflt) - 1 : iflt;
            if (i < 0 || i >= (long)m_weight_bins.size()) return;
            m_weight_bins.at(i) += multiplier;
            m_value_bins.at(i) += y * multiplier;
        }

        /** Sample weight at a particular x-position (for histograms). */
        inline void sampleWeight(double x, double multiplier = 1.0) {
            sampleValue(x, 0.0, multiplier);
        }

        std::string toString(bool endbar=false) const;

        // For reduction in omp parallel for loops.
        GraphDataBinned& merge(const GraphDataBinned& other);

        GraphDataBinned& divide(const GraphDataBinned& other);
};

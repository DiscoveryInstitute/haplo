#include "FileOutputs.h"

#include <fstream>
#include <string>
#include "FileWriter.h"

using namespace std;


FileOutputs::FileOutputs(const string& pfx) : m_file_prefix(pfx) {
    ofstream ofs(pfx + "test");
    if (!ofs) throw runtime_error("Directory not available to write files: " + pfx);
}


void FileOutputs::assertNoExistingOutputFiles() const {

    auto filenames = {m_file_prefix + "Q.txt",
                      m_file_prefix + "S.txt",
                      m_file_prefix + "Pi.txt",
                      m_file_prefix + "Phi.txt",
                      m_file_prefix + "snp_dist.txt",
                      m_file_prefix + "r2.txt",
                      m_file_prefix + "dprime.txt",
                      m_file_prefix + "CLDP.txt",
                      m_file_prefix + "aCLDP.txt",
                      m_file_prefix + "nKLD.txt",
                      m_file_prefix + "pair_spectrum.txt"};

    for (string filename : filenames) {
        ifstream ifs(filename);
        if (ifs) throw runtime_error("Output file already exists: " + filename);
    }
}


void FileOutputs::writeOutputFiles(const Analysis& analysis) const {

    FileWriter::write(m_file_prefix + "Q.txt",
                      to_string(analysis.total_blocks),
                      "Q - total generated haploblocks");
    FileWriter::write(m_file_prefix + "S.txt",
                      to_string(analysis.total_snp_sites),
                      "S - total number of SNPs (mutations)");
    FileWriter::write(m_file_prefix + "Pi.txt",
                      to_string(analysis.sn_diversity),
                      "Pi - total diversity per SN site");

    FileWriter::write(m_file_prefix + "Phi.txt",
                      analysis.snp_freq_spectrum.toString(true),
                      "Phi - SNP frequency spectrum");
    FileWriter::write(m_file_prefix + "snp_dist.txt",
                      analysis.snp_distribution.toString(true),
                      "SNP distribution along chromosome");
    FileWriter::write(m_file_prefix + "r2.txt",
                      analysis.correlation_function.toString(true),
                      "r2 - correlation vs genome distance");
    FileWriter::write(m_file_prefix + "dprime.txt",
                      analysis.dprime_function.toString(true),
                      "D' - correlation vs genome distance");
    FileWriter::write(m_file_prefix + "CLDP.txt",
                      analysis.complete_linkage_function.toString(true),
                      "CLDP - complete linkage-disequilibrium proportion vs genome distance");
    FileWriter::write(m_file_prefix + "aCLDP.txt",
                      analysis.almost_complete_linkage_function.toString(true),
                      "aCLDP - almost-complete linkage-disequilibrium proportion vs genome distance");
    FileWriter::write(m_file_prefix + "nKLD.txt",
                      analysis.normalized_kullback_leibler_function.toString(true),
                      "nKLD - normalised Kullback-Leibler function vs genome distance");
    FileWriter::write(m_file_prefix + "sigma2.txt",
                      analysis.sigma_squared_function.toString(true),
                      "sigma2 - average(numerator(r2))/average(denominator(r2))");

    FileWriter::write(m_file_prefix + "pair_spectrum.txt",
                      analysis.freq_pair_spectrum.toString(true),
                      "pair_spectrum - SNPs frequency pair spectrum");


}

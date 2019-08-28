#include <cstdio>
#include <cstdlib>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include "core/Analysis.h"
#include "core/AncestralChromosomes.h"
#include "core/AncestralRecombinationGraph.h"
#include "core/Chromosome.h"
#include "core/Parser.h"
#include "core/PopulationStructureHistory.h"
#include "hpc/FileCheckpointer.h"
#include "hpc/FileMemorySaver.h"
#include "hpc/FileOutputs.h"
#include "hpc/FileReader.h"
#include "hpc/FileSNPsWriter.h"
#include "hpc/FileWriter.h"
#include "hpc/LoggerSmart.h"

using namespace std;

constexpr bool BACKWARDS = false;
constexpr bool FORWARDS = true;

int main(int argc, char** argv) {

    //--------------------------------------------------------------------
    // State.
    mt19937_64 rng;
    PopulationStructureHistory population_history;
    MutationRateHistory mutation_history;
    HaploBlockBoundaries boundaries;
    AncestralRecombinationGraph graph;
    HaploBlockVariants variants;
    AncestralChromosomes chromosomes;

    long gen;
    long nGenerations;
    long extantGeneration;
    long foundingGeneration;


    try {

    //--------------------------------------------------------------------
    // Read command line arguments.
    bool restart_at_checkpoint = false;
    string input_filename;
    for (int iarg=1; iarg<argc; ++iarg) {
        string argstr(argv[iarg]);
        if (argstr=="--restart") restart_at_checkpoint = true;
        else {
            if(input_filename=="") input_filename = argstr;
            else throw invalid_argument("Need one input filename");
        }
    }

    //--------------------------------------------------------------------
    // Read input file.
    Parser inputs(FileReader::parseTOML(input_filename));
    // Set up logger.
    const bool verbose = inputs.getBool("verbose_logging");
    LoggerSmart logger(verbose, 10); // log progress every 10s
    logger.logEvent("Inputs read:\n" + inputs.toString());
    // Set up memory unloader.
    const bool use_memory_unloading = inputs.getBool("use_memory_unloading");
    FileMemorySaver memory_saver(inputs.getString("temporary_file_prefix"));
    logger.logEvent("Memory saver created.");
    // Set up memory unloader.
    FileCheckpointer checkpointer("checkpoint.bin", 60); // checkpoint every 60s
    logger.logEvent("Checkpointer created.");
    // Before we start, check that there are no output files already there.
    //   Don't want to overwrite results that were expensive to compute!
    const FileOutputs file_outputs(inputs.getString("output_file_prefix"));
    file_outputs.assertNoExistingOutputFiles();
    logger.logEvent("Checked no output files.");

    const double alpha = inputs.getNonNegativeDouble("fertility_parameter_alpha");
    const double beta  = inputs.getNonNegativeDouble("mating_parameter_beta");
    const double recombination_rate = inputs.getNonNegativeDouble("recombination_rate");
    const bool cull_non_parents = inputs.getBool("cull_non_parents");
    const bool cull_na_parents = inputs.getBool("cull_nonancestral_parents");
    const bool hide_na_blocks = inputs.getBool("hide_nonancestral_blocks");

    if (restart_at_checkpoint) {
        bool direction;
        checkpointer.load(direction, gen, rng, population_history, mutation_history, boundaries, graph, variants, chromosomes);
        nGenerations = population_history.getNumberOfGenerations();
        extantGeneration = 0L;
        foundingGeneration = nGenerations - 1L;
        if (direction==BACKWARDS) goto continue_backwards;
        if (direction==FORWARDS)  goto continue_forwards;
    }

    //--------------------------------------------------------------------
    // Set Random Seed.
    {
        const string random_seed = inputs.getString("random_seed");
        seed_seq seed(random_seed.begin(), random_seed.end());
        rng = mt19937_64(seed);
    }

    //--------------------------------------------------------------------
    // Read or construct PopulationStructureHistory.
    {
        const string population_history_filename = inputs.getString("population_structure_history_file");
        population_history = PopulationStructureHistory(
                (population_history_filename != "")
                ? FileReader::readPopulationStructureHistory(population_history_filename)
                : PopulationStructureHistory::create(
                        inputs.getString("population_growth_type"),
                        inputs.getPositiveLong("population_number_of_generations"),
                        inputs.getPositiveLong("population_initial"),
                        inputs.getPositiveLong("population_final"),
                    inputs.getPositiveLong("population_change_parameter")));

        nGenerations = population_history.getNumberOfGenerations();
        extantGeneration = 0L;
        foundingGeneration = nGenerations - 1L;

        if (logger.verbose()) {
            string s = (population_history_filename != "")
                ? "Population History was read from file:\n"
                : "Population History was constructed:\n";
            for (long gen = nGenerations-1; gen >= 0L; --gen) {
                long nmales = population_history.getNumberOfMales(gen);
                long nfemales = population_history.getNumberOfFemales(gen);
                s += to_string(gen) + " " + to_string(nmales) + " " + to_string(nfemales) + "\n";
            }
            logger.logEvent(s);
        }
    }

    //--------------------------------------------------------------------
    // Read or construct MutationRateHistory.
    {
        const string mutation_history_filename = inputs.getString("mutation_rate_history_file");
        mutation_history = MutationRateHistory(
                (mutation_history_filename != "")
                ? FileReader::readMutationRateHistory(mutation_history_filename)
                : MutationRateHistory::Constant(nGenerations, inputs.getNonNegativeDouble("mutation_rate")));
    }

    //--------------------------------------------------------------------
    // Generate AncestralRecombinationGraph.
    {
        boundaries = HaploBlockBoundaries(
                inputs.getPositiveLong("chromosome_length"),
                inputs.getPositiveLong("maximum_blocks"),
                rng);

        graph = AncestralRecombinationGraph(&population_history, &boundaries,
                                            inputs.getChromosomeType("chromosome_type"))
                .initialiseExtantGeneration(inputs.getPositiveLong("population_sample"));

        logger.logEvent("HBB and ARG initialised.\n");
        FileWriter::write("log_ancestral.txt","# gen size_of_ancestral_haploblocks");
    }

    for (gen = extantGeneration; gen<foundingGeneration; ++gen) {

        {
            graph.backwardsCoalesceOneGeneration(
                    gen, alpha, beta, recombination_rate, rng,
                    cull_non_parents, cull_na_parents, hide_na_blocks);
        }

        logger.logProgress("Generation " + to_string(gen) + " of ARG done.");
        FileWriter::append("log_ancestral.txt", to_string(gen) + " " + to_string(graph.getSizeOfAncestralHaploBlocks(gen)) + "\n");
        if (use_memory_unloading) memory_saver.unloadGenerationToFile(graph, gen);
        if (checkpointer.timeDue()) {
            checkpointer.save(BACKWARDS, gen, rng, population_history, mutation_history, boundaries, graph, variants, chromosomes);
            memory_saver.deleteDiscardedGenerations();
        }
        continue_backwards:
        FileReader::stopIfFindStopFile();
    }

    //--------------------------------------------------------------------
    // Propagate AncestralRecombinationGraph and generate mutations
    //   to create AncestralChromosomes.
    {

        const long max_nodes_per_block = inputs.getPositiveLong("population_sample") * 6;
        // 2x haploid samples -> 4x nodes + leeway = 6x                                ^
        variants = HaploBlockVariants(
                        boundaries,
                        max_nodes_per_block,
                        inputs.getPositiveLong("primordial_block_length"),
                        inputs.getNonNegativeDouble("primordial_probability_dimorphic"),
                        inputs.getNonNegativeDouble("primordial_probability_tetramorphic"),
                        inputs.getNonNegativeDouble("primordial_diversity"),
                        inputs.getBool("use_mutation_loci"),
                        rng);

        chromosomes = AncestralChromosomes(&graph, &variants)
                        .initialiseFoundingGeneration();

        logger.logEvent("\nHBV and AC initialised.\n");
    }

    for (gen = foundingGeneration-1; gen>=extantGeneration; --gen) {

        if (use_memory_unloading) memory_saver.reloadGenerationFromFile(graph, gen);

        {
            const double mutation_rate = mutation_history.getMutationRate(gen);
            chromosomes.forwardsDropOneGeneration(gen, mutation_rate, rng);
        }

        logger.logProgress("Generation " + to_string(gen) + " of AC done.");
        if (use_memory_unloading) memory_saver.discardGeneration(chromosomes, gen+1);
        if (use_memory_unloading) memory_saver.discardGeneration(graph, gen);
        if (variants.sampleChromosomesAndSimplify(chromosomes.getGeneration(gen))) logger.logEvent("HaploBlockVariants simplified.");
        if (checkpointer.timeDue() || gen == extantGeneration) {
            checkpointer.save(FORWARDS, gen, rng, population_history, mutation_history, boundaries, graph, variants, chromosomes);
            memory_saver.deleteDiscardedGenerations();
        }
        continue_forwards:
        FileReader::stopIfFindStopFile();
    }

    logger.logEvent("\nSimulation complete.");

    //--------------------------------------------------------------------
    // Write SNPs
    if (inputs.getBool("output_final_SNP_sites")) {
        logger.logEvent("Writing SNPs.");

        FileSNPsWriter::WriteSNPs(inputs.getString("output_file_prefix")+"SNPs.sites",
                                  variants, chromosomes.getExtantGeneration(), rng);

        logger.logEvent("Wrote SNPs.");
    }

    //--------------------------------------------------------------------
    // Calculate final output statistics.

    auto analysis = Analysis::Calculate(logger, inputs, variants, chromosomes.getExtantGeneration());

    file_outputs.writeOutputFiles(analysis);
    logger.logEvent("Results written to files.");

    const string unused_inputs = inputs.getUnusedKeysAsString();
    if (unused_inputs != "") {
        printf("There were unused inputs. See log_unused.txt.\n");
        FileWriter::write("log_unused.txt",
                          unused_inputs.c_str(),
                          "Keys from the input file that were not used.\n");
    }

    return EXIT_SUCCESS;

    //--------------------------------------------------------------------

    }
    catch(const exception& e) {
        fprintf(stderr, "%s\n", e.what());
        return EXIT_FAILURE;
    }

}

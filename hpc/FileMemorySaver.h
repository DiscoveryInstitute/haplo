#pragma once
#include <string>
#include <vector>
#include "core/AncestralChromosomes.h"
#include "core/AncestralRecombinationGraph.h"

class FileMemorySaver {

    const std::string m_file_prefix; // can be or include a directory name
    std::vector<std::string> m_discard_list; // list of files to delete (at checkpoint)

public:

    FileMemorySaver(const std::string& pfx="");

    void unloadGenerationToFile(AncestralRecombinationGraph& arg, long gen) const;
    void reloadGenerationFromFile(AncestralRecombinationGraph& arg, long gen) const;
    void discardGeneration(AncestralRecombinationGraph& arg, long gen);

    void unloadGenerationToFile(AncestralChromosomes& ac, long gen) const;
    void reloadGenerationFromFile(AncestralChromosomes& ac, long gen) const;
    void discardGeneration(AncestralChromosomes& ac, long gen);

    void deleteDiscardedGenerations();
};

#pragma once

#include <map>
#include <string>
#include "core/MutationRateHistory.h"
#include "core/PopulationStructureHistory.h"

class FileReader {

    public:

        static std::map<std::string,std::string> parseTOML(std::string input_filename);

        static PopulationStructureHistory readPopulationStructureHistory(std::string input_filename);

        static MutationRateHistory readMutationRateHistory(std::string input_filename);

        static void stopIfFindStopFile();
};

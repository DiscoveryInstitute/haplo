#pragma once

#include <string>
#include "core/Analysis.h"

class FileOutputs {

    const std::string m_file_prefix;

public:

    FileOutputs(const std::string& pfx="");
    void assertNoExistingOutputFiles() const;
    void writeOutputFiles(const Analysis& analysis) const;

};

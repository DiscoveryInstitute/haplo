#pragma once

#include <string>

class FileWriter {

public:

    static void write(const std::string& output_filename,
                      const std::string& content,
                      const std::string& comment="");

    static void append(const std::string& output_filename,
                       const std::string& content,
                       const std::string& comment="");
};

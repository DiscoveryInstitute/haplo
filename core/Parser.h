#pragma once

#include <map>
#include <set>
#include <string>
#include "Chromosome.h"

class Parser {

    friend class FileCheckpointer;

private :

    std::map<std::string, std::string> m_map;
    std::set<std::string> m_used;

public:

    Parser() {}// for creating unassigned variables

    static bool parseBool (
            const std::string&  value, const std::string& key = "");
    static double parseDouble (
            const std::string&  value, const std::string& key = "");
    static long parseLong (
            const std::string&  value, const std::string& key = "");
    static ChromosomeType parseChromosomeType (
            const std::string&  value, const std::string& key = "");

    static double parseNonNegativeDouble(
            const std::string&  value, const std::string& key = "");
    static long parsePositiveLong(
            const std::string&  value, const std::string& key = "");

    Parser(std::map<std::string, std::string> map) : m_map(map) {}

    const std::string& getString(const std::string& key);

    bool getBool(const std::string&  key);
    double getDouble(const std::string&  key);
    long getLong(const std::string&  key);
    ChromosomeType getChromosomeType(const std::string&  key);

    double getNonNegativeDouble(const std::string&  key);
    long getPositiveLong(const std::string& key);

    std::string getUnusedKeysAsString() const;

    std::string toString() const;

    bool operator==(const Parser& o) const { return m_map == o.m_map && m_used == o.m_used; }

};

#include "Parser.h"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>

using namespace std;


string msg(const string& type, const string& value, const string& key) {
    return "Not a " + type + ": " + value + " (key: " + key + ")";
}

// Static methods.

bool Parser::parseBool(const string& value, const string& key) {
    string s = value;
    transform(s.begin(), s.end(), s.begin(), ::tolower);
    if (s == "true") return true;
    if (s == "false") return false;
    throw invalid_argument(msg("boolean", value, key));
}

double Parser::parseDouble(const string& value, const string& key) {
    try {
        return stod(value);
    } catch (const exception& e) {
        throw invalid_argument(msg("double", value, key));
    }
}

long Parser::parseLong(const string& value, const string& key) {
    try {
        double dval = stod(value);
        if (fabs(dval) > (double)numeric_limits<long>::max())
            throw invalid_argument("Number too large: " + value + " (key: " + key + ")");
        return (long) dval;
    } catch (const exception& e) {
        throw invalid_argument(msg("long", value, key));
    }
}

ChromosomeType Parser::parseChromosomeType(
        const string& value, const string& key) {
    string s = value;
    transform(s.begin(), s.end(), s.begin(), ::tolower);
    if (s == "autosome") return AUTO;
    if (s == "x chromosome") return X;
    if (s == "y chromosome") return Y;
    if (s == "mitochondrial dna") return MITO;
    throw invalid_argument(msg("ChromosomeType", value, key));
}

double Parser::parseNonNegativeDouble(
        const string& value, const string& key) {
    double dval = parseDouble(value);
    if (0.0 <= dval) return dval;
    throw invalid_argument(msg("non-negative double", value, key));
}

long Parser::parsePositiveLong(
        const string& value, const string& key) {
    long lval = parseLong(value);
    if (0 < lval) return lval;
    throw invalid_argument(msg("positive long", value, key));
}

// Member methods.

const string& Parser::getString(const string& key) {
    try {
        m_used.emplace(key);
        return m_map.at(key);
    } catch (const out_of_range& oor) {
        throw out_of_range("Key not found: " + key);
    }
}

bool Parser::getBool(const string& key) {
    return parseBool(getString(key), key);
}

double Parser::getDouble(const string& key) {
    return parseDouble(getString(key), key);
}

long Parser::getLong(const string& key) {
    return parseLong(getString(key), key);
}

ChromosomeType Parser::getChromosomeType(const string& key) {
    return parseChromosomeType(getString(key), key);
}

double Parser::getNonNegativeDouble(const string& key) {
    return parseNonNegativeDouble(getString(key), key);
}

long Parser::getPositiveLong(const string& key) {
    return parsePositiveLong(getString(key), key);
}

string Parser::getUnusedKeysAsString() const {
    string s;
    for (const auto& entry : m_map) {
        const string& key = entry.first;
        if (m_used.count(key)==0) {
            s += key + "\n";
        }
    }
    return s;
}

string Parser::toString() const {
    string s;
    for (const auto& entry : m_map) {
        s += entry.first + " : " + entry.second + "\n";
    }
    return s;
}


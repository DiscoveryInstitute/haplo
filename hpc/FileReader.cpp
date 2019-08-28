#include "FileReader.h"

#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "core/Parser.h"

using namespace std;

inline long nearestLong(double x) { return (long)(x>0 ? x+0.5 : x-0.5); }

/**
 *  Reads a file using a very simple TOML-like format.
 *  Returns a map/dictionary from string keys to string values.
 *    Key/value pairs are given line by line:
 *      keyA = valueA
 *      keyB = valueB
 *      ...
 *    Except that
 *      whitespace at the start of a line means it continues previous line
 *      # means comment - it and the rest of the line will be ignored
 *      other whitespace is condensed into single spaces
 */
map<string,string> FileReader::parseTOML(string input_filename) {
    ifstream ifs(input_filename, ifstream::in);
    if (ifs.fail()) throw invalid_argument("Cannot open " + input_filename);
    stringstream raw_ss;
    raw_ss << ifs.rdbuf();
    string s = raw_ss.str();
    s = regex_replace(s, regex(R"([\n\r]+)"), "\n");   // replace native line endings
    s = regex_replace(s, regex(R"(#[^\n]*\n)"), "\n"); // remove comments
    s = regex_replace(s, regex(R"(\s*\n)"), "\n");     // remove wspace before newline
    s = regex_replace(s, regex(R"(\n\s+)"), " ");      // continue line if starts with space
    s = regex_replace(s, regex(R"([^\S\n]+)"), " ");   // remove extra space (except newlines)
    s = regex_replace(s, regex(R"([^\S\n]*=[^\S\n]*)"), "="); // remove space around '='
    map<string,string> kvmap;
    istringstream lines_ss(s);
    for (string line; getline(lines_ss,line); ) {
        if (line.length()==0) continue;
        istringstream tokens_ss(line);
        vector<string> tokens;
        for (string token; getline(tokens_ss,token,'='); ) {
            tokens.push_back(token);
        }
        if (tokens.size()==0 || tokens[0].empty())
                throw invalid_argument("No key in this entry: " + line);
        if (tokens.size()==1) tokens.push_back("");
        if (tokens.size()>2)
                throw invalid_argument("Extra values in this entry: " + line);
        kvmap.emplace(tokens[0], tokens[1]);
    }
    return kvmap;
}

/**
 *  Read PopulationStructureHistory from a text file.
 *  The format of each line is: 'gen nmales nfemales' .
 *  Generations are generations before present, so the first line gives
 *    the total number of generations to simulate back.
 *  Generations must be in descending order (forward in time to present).
 *  The first and last generation must be specified.
 *  Any generations not explicitly specified, are interpolated
 *    piecewise-linearly, and non-integer values are rounded down.
 *  Example:
 *     4 4 4
 *     3 1 1
 *     0 4 4
 *  becomes:
 *     4 4 4
 *     3 1 1
 *     2 2 2
 *     1 3 3
 *     0 4 4
 */
PopulationStructureHistory FileReader::readPopulationStructureHistory(string input_filename) {

    vector<pair<long,long>> history;
    long previousGen=-1L; // to keep stupid compilers happy

    ifstream ifs(input_filename, ifstream::in);
    if (ifs.fail()) throw invalid_argument("Cannot open " + input_filename);

    // For each line in the text file.
    for (string line; getline(ifs,line); ) {
        line = regex_replace(line, regex(R"([\s]*#.*)"), ""); // remove comments (and trim space)
        line = regex_replace(line, regex(R"([\s]+)"), " "); // trim whitespace to single space
        if (line.length()==0) continue;

        // Read the point.
        istringstream tokens_ss(line);
        vector<string> tokens;
        for (string token; getline(tokens_ss, token, ' '); ) {
            if (!token.empty()) tokens.push_back(token);
        }
        if (tokens.size() != 3) {
            throw invalid_argument("Need 3 args per line: gen nmales nfemales. "
                                   "Got: " + line);
        }
        long nextGen = (long)Parser::parseNonNegativeDouble(tokens[0]);
        long nextNMales = Parser::parsePositiveLong(tokens[1]);
        long nextNFemales = Parser::parsePositiveLong(tokens[2]);

        // If first point, create the vector and enter point.
        if (history.empty()) {
            const long ngens = nextGen+1;
            history = vector<pair<long,long>>(ngens);
            history[nextGen] = { nextNMales, nextNFemales };
            previousGen = nextGen;
        } else {
        // Otherwise, interpolate from the last point.
            if (nextGen >= previousGen)
                throw invalid_argument(
                        "Need descending generation numbers. Got" + line);
            auto& previousEntry = history[previousGen];
            long previousNMales = previousEntry.first;
            long previousNFemales = previousEntry.second;
            double norm = 1.0 / (previousGen - nextGen);
            history[nextGen] = { nextNMales, nextNFemales };
            for (long gen = nextGen+1; gen < previousGen; ++gen) {
                double t = (gen - nextGen) * norm;
                long nMales = nearestLong(t * previousNMales + (1-t) * nextNMales);
                long nFemales = nearestLong(t * previousNFemales + (1-t) * nextNFemales);
                history[gen] = { nMales, nFemales };
            }
            previousGen = nextGen;
        }
    }

    // Check the history is not empty.
    auto& extantEntry = history[0];
    if (extantEntry.first==0L && extantEntry.second==0L)
        throw invalid_argument("No extant generation!");

    return PopulationStructureHistory(history);
}

/**
 *  Read MutationRateHistory from a text file.
 *  The format of each line is: 'gen rate'.
 *  See readPopulationStructureHistory for details about interpolation.
 */
MutationRateHistory FileReader::readMutationRateHistory(string input_filename) {

    vector<double> history;
    long previousGen=-1L; // to keep stupid compilers happy

    ifstream ifs(input_filename, ifstream::in);
    if (ifs.fail()) throw invalid_argument("Cannot open " + input_filename);

    // For each line in the text file.
    for (string line; getline(ifs,line); ) {
        line = regex_replace(line, regex(R"([\s]*#.*)"), ""); // remove comments (and trim space)
        line = regex_replace(line, regex(R"([\s]+)"), " "); // trim whitespace to single space
        if (line.length()==0) continue;

        // Read the point.
        istringstream tokens_ss(line);
        vector<string> tokens;
        for (string token; getline(tokens_ss, token, ' '); ) {
            if (!token.empty()) tokens.push_back(token);
        }
        if (tokens.size() != 2) {
            throw invalid_argument("Need 2 args per line: gen rate. "
                                   "Got: " + line);
        }
        long nextGen = (long)Parser::parseNonNegativeDouble(tokens[0]);
        double nextRate = Parser::parseNonNegativeDouble(tokens[1]);

        // If first point, create the vector and enter point.
        if (history.empty()) {
            const long ngens = nextGen+1;
            history = vector<double>(ngens);
            history[nextGen] = nextRate;
            previousGen = nextGen;
        } else {
        // Otherwise, interpolate from the last point.
            if (nextGen >= previousGen)
                throw invalid_argument(
                        "Need descending generation numbers. Got" + line);
            double previousRate = history[previousGen];
            double norm = 1.0 / (previousGen - nextGen);
            history[nextGen] = nextRate;
            for (long gen = nextGen+1; gen < previousGen; ++gen) {
                double t = (gen - nextGen) * norm;
                double rate = t * previousRate + (1-t) * nextRate;
                history[gen] = rate;
            }
            previousGen = nextGen;
        }
    }

    return MutationRateHistory(history);
}

void FileReader::stopIfFindStopFile() {
    ifstream ifs("stop", ifstream::in);
    if (!ifs.fail()) {
        printf("Stopping because 'stop' file found\n");
        exit(0);
    }
}


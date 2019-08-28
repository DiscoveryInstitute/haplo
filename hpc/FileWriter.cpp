#include "FileWriter.h"

#include <fstream>
#include <stdexcept>
#include <string>

using namespace std;

void FileWriter::write(
        const string &output_filename,
        const string &content,
        const string &comment) {

    ofstream ofs(output_filename, ofstream::out | ofstream::trunc);
    if (ofs.fail()) throw invalid_argument("Cannot open " + output_filename);

    if (comment != "") ofs << "# " << comment << endl;
    ofs << content << endl;

}

void FileWriter::append(
        const string &output_filename,
        const string &content,
        const string &comment) {

    ofstream ofs(output_filename, ofstream::out | ofstream::app);
    if (ofs.fail()) throw invalid_argument("Cannot open " + output_filename);

    if (comment != "") ofs << "# " << comment << endl;
    ofs << content << endl;
}

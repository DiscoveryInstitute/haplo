#pragma once
#include <fstream>
#include <string>
#include <utility>
#include <vector>

template<typename V>
inline void writePrimitive(std::ofstream& ofs, const V& p) {
    ofs.write((char*)&p, sizeof(p));
}
template<typename V>
inline void readPrimitive(std::ifstream& ifs, V& p) {
    ifs.read((char*)&p, sizeof(p));
}

inline void writeString(std::ofstream& ofs, const std::string& s) {
    writePrimitive(ofs, s.size());
    ofs.write(s.c_str(),s.size());
}
inline void readString(std::ifstream& ifs, std::string& s) {
    size_t sz;
    readPrimitive(ifs,sz);
    s.resize(sz);
    ifs.read(&s[0],s.size());
}

/** VECTOR */
template<typename V>
inline void writePrimitiveVector(std::ofstream& ofs, const std::vector<V>& v) {
    long n_units = v.size();
    writePrimitive(ofs, n_units);
    char* ptr = (char*)&v[0];
    ofs.write(ptr, n_units * sizeof(V));
}
template<typename V>
inline void readPrimitiveVector(std::ifstream& ifs, std::vector<V>& v) {
    long n_units;
    readPrimitive(ifs, n_units);
    v.resize(n_units);
    char* ptr = (char*)&v[0];
    ifs.read(ptr, n_units * sizeof(V));
}

/** MAP */
template<typename U, typename V>
inline void writePrimitiveMap(std::ofstream& ofs, const std::map<U,V>& m) {
    long n_units = m.size();
    writePrimitive(ofs, n_units);
    for (const auto& e : m) {
        writePrimitive(ofs, e);
    }
}
template<typename U, typename V>
inline void readPrimitiveMap(std::ifstream& ifs, std::map<U,V>& m) {
    long n_units;
    readPrimitive(ifs, n_units);
    for (int i=0; i<n_units; ++i) {
       std::pair<U,V> e;
       readPrimitive(ifs, e);
       m.emplace(e.first, e.second);
    }
}



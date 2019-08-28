#pragma once
#include <cstdio>

class Logger {
    const bool m_verbose;
protected:
    inline void log(const std::string& msg) {
        if (m_verbose) {
            printf("%s\n", msg.c_str());
            fflush(stdout);
        }
    }
public:
    Logger(bool vb) : m_verbose(vb) {}
    bool verbose() { return m_verbose; }
    virtual void logEvent(const std::string& msg) {
        log(msg);
    }
    virtual void logProgress(const std::string& msg) {
        log(msg);
    }
};


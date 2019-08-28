#pragma once
#include <ctime>
#include "../core/Logger.h"

class LoggerSmart : public Logger {

    private:

        const long minimum_interval;
        long previous_log_time;

    public:

        LoggerSmart(bool vb, long min_int)
            : Logger(vb),
              minimum_interval(min_int*CLOCKS_PER_SEC),
              previous_log_time(-min_int*CLOCKS_PER_SEC) {}

        virtual void logEvent(const std::string& msg) {
             log(msg);
             previous_log_time = -minimum_interval;
        }

        virtual void logProgress(const std::string& msg) {
            long current_time = clock();
            if (current_time < previous_log_time + minimum_interval) return;
            previous_log_time = current_time;
            log(msg);
        }

};


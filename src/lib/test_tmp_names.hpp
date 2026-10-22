#pragma once

#include <unistd.h>
#include <string>

// Unique temporary file names for the tests (process id + counter), so that
// concurrent runs do not clobber each other and a stale file left by an
// interrupted run is never silently reused.
inline std::string unique_tmp_name(const std::string &base) {
    static int counter = 0;
    return "mutselc_" + std::to_string(::getpid()) + "_" + std::to_string(counter++) + "_" + base;
}

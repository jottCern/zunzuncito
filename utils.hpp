#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <string>
#include <cmath>

#include "fwd.hpp"

// perform glob filename expansion of given pattern.
// The result will be sorted by filename.
// Throws a runtime_error exception in case of an error.
std::vector<std::string> glob(const std::string& pattern);

// print the event content to standard out
void dump(const event & evt, int njetsmax = 4, bool with_mc_info = true);


inline float hypotf(float x, float y){
    return sqrtf(x*x + y*y);
}

#endif

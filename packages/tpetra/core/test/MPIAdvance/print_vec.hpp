#pragma once

#include <iostream>

template <typename V>
void print_vec(const std::string &label, const V &vec) {
    std::stringstream ss;
    ss << label << ":";
    for (const auto &e : vec) {
        ss << " " << e;
    }
    std::cerr << ss.str() << std::endl;
}

template <typename V>
void print_vec(const int rank, const std::string &label, const V &vec) {
    std::stringstream ss;
    ss << "[" << rank << "]" << label;
    print_vec(ss.str(), vec);
}
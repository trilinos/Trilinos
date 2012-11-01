#ifndef STRING_UTILS_H
#define STRING_UTILS_H

inline
std::string intToString(int i)
{
    std::stringstream ss;
    ss << i;
    return ss.str();
}

#endif

#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <sstream>
#include <vector>
namespace stringUtils
{

inline
std::string intToString(int i)
{
    std::stringstream ss;
    ss << i;
    return ss.str();
}

inline
void convertCharArrayToStringVector(int numArgs, const char** charArray, std::vector<std::string> &stringVector)
{
    stringVector.resize(numArgs);
    for(int i=0; i < numArgs; i++)
    {
        stringVector[i] = std::string(charArray[i]);
    }
}

inline
double stringToDouble(const std::string &doubleString)
{
    double value = 0;
    std::istringstream inputStringStream(doubleString);
    inputStringStream >> value;
    return value;
}

}
#endif

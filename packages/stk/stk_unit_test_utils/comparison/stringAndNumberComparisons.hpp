#ifndef STRING_AND_NUMBER_COMPARISONS_HPP
#define STRING_AND_NUMBER_COMPARISONS_HPP

#include <string>
#include <iostream>

namespace unitTestUtils
{

inline bool isNear(double a, double b, double tolerance)
{
    bool isNear = false;
    double diff = a - b;
    if(diff > -tolerance && diff < tolerance)
    {   
        isNear = true;
    }   
    return isNear;
}

bool approximatelyEqualAsNumbers(const std::string &expectedWord, const std::string &actualWord, double tol);

bool areStringsEqualWithToleranceForNumbers(const std::string &expectedString, const std::string &actualString, double tol);

}
#endif

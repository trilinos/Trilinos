#include <stk_unit_test_utils/stringAndNumberComparisons.hpp>

#include <iostream>
#include <sstream>
#include <string>

namespace stk {
namespace unit_test_util {

bool approximatelyEqualAsNumbers(const std::string &expectedWord, const std::string &actualWord, double tol)
{
    std::istringstream expectedStream(expectedWord);
    std::istringstream actualStream(actualWord);
    double expectedDouble, actualDouble;
    expectedStream >> expectedDouble;
    actualStream >> actualDouble;
    bool isEqual = false;
    bool didStringsConvertToNumbers = !expectedStream.fail() && !actualStream.fail();
    if(didStringsConvertToNumbers && isNear(expectedDouble, actualDouble, tol))
    {
        isEqual = true;
    }
    return isEqual;
}

bool areStringsEqualWithToleranceForNumbers(const std::string &expectedString, const std::string &actualString, double tol)
{
    std::istringstream expectedStream(expectedString);
    std::istringstream actualStream(actualString);
    bool isEqual = true;
    while(!expectedStream.fail())
    {
        std::string expectedWord;
        expectedStream >> expectedWord;
        std::string actualWord;
        actualStream >> actualWord;
        if ( expectedWord.substr(0,4) != "----" ) // only compare non dash lines
        {
            if(!expectedStream.fail() && expectedWord != "SKIP")
            {
                if(expectedWord != actualWord && !approximatelyEqualAsNumbers(expectedWord, actualWord, tol))
                {
                    std::cerr << "Expected: '" << expectedWord << "' but got '" << actualWord << "'" << std::endl;
                    isEqual = false;
                }
            }
        }
    }
    return isEqual;
}

namespace simple_fields {

bool isNear(double a, double b, double tolerance)
{
  return stk::unit_test_util::isNear(a, b, tolerance);
}

bool approximatelyEqualAsNumbers(const std::string &expectedWord, const std::string &actualWord, double tol) {
  return stk::unit_test_util::approximatelyEqualAsNumbers(expectedWord, actualWord, tol);
}

bool areStringsEqualWithToleranceForNumbers(const std::string &expectedString, const std::string &actualString, double tol) {
  return stk::unit_test_util::areStringsEqualWithToleranceForNumbers(expectedString, actualString, tol);
}

} // namespace simple_fields

}
}

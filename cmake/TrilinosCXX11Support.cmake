INCLUDE(CheckCXXSourceRuns)

FUNCTION(CHECK_CXX11_SUPPORT VARNAME)

  # support for >> in addition to > > when closing double templates
  SET(SOURCE_CXX11_CONSECUTIVE_RIGHT_ANGLE_BRACKETS
  "
#include <vector>
int main() {
  // check >> closing brackets
  std::vector<std::vector<float>> vecvecfloat(1);
  vecvecfloat[0].resize(1);
  vecvecfloat[0][0] = 0.0f;
  return 0;
}
  "
  )
  CHECK_CXX_SOURCE_RUNS("${SOURCE_CXX11_CONSECUTIVE_RIGHT_ANGLE_BRACKETS}" CXX11_CONSECUTIVE_RIGHT_ANGLE_BRACKETS)

  # support for auto and typedecl()
  SET(SOURCE_CXX11_AUTOTYPEDVARIABLES
  "
#include <vector>
int main() {
  std::vector<int> vec(10);
  auto b        = vec.begin();
  decltype(b) e = vec.end();
  std::fill(b,e,1);
  return 0;
}
  "
  )
  CHECK_CXX_SOURCE_RUNS("${SOURCE_CXX11_AUTOTYPEDVARIABLES}" CXX11_AUTOTYPEDVARIABLES)

  # support for lambda expressions
  SET(SOURCE_CXX11_LAMBDAS
  "
#include <vector>
#include <algorithm>
int main() {
  // two examples, taken from the wikipedia article on C++0x :)
  std::vector<int> some_list;
  int total = 0;
  int value = 5;
  std::for_each(some_list.begin(), some_list.end(), [&total](int x) {
      total += x;
      });
  std::for_each(some_list.begin(), some_list.end(), [&, value](int x) {
      total += x * value;
      });
  // 
  return 0;
}
  "
  )
  CHECK_CXX_SOURCE_RUNS("${SOURCE_CXX11_LAMBDAS}" CXX11_LAMBDAS)

  IF (NOT CXX11_CONSECUTIVE_RIGHT_ANGLE_BRACKETS OR NOT CXX11_AUTOTYPEDVARIABLES OR NOT CXX11_LAMBDAS)
    SET(${VARNAME} FALSE PARENT_SCOPE)
  ELSE()
    SET(${VARNAME} TRUE PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()

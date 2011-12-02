# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER

INCLUDE(CheckCXXSourceRuns)

FUNCTION(TRIBITS_CHECK_CXX11_SUPPORT VARNAME)

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

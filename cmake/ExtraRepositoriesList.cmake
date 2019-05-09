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


#
# List of extra repositories that contain extra Trilinos packages.
#
# See documentation for the function
# TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES() in the TriBITS Guide and
# Reference:
#
#  https://tribits.org/doc/TribitsDevelopersGuide.html#tribits-project-define-extra-repositories
#

IF (NOT "$ENV{Trilinos_REPOS_URL_BASE}" STREQUAL "")
  SET(Trilinos_REPOS_URL_BASE  $ENV{Trilinos_REPOS_URL_BASE}
    CACHE STRING "Set from env var Trilinos_REPOS_URL_BASE" FORCE)
ELSE()
  SET(Trilinos_REPOS_URL_BASE_DEFAULT git@github.com:trilinos/)
  SET(Trilinos_REPOS_URL_BASE  ${Trilinos_REPOS_URL_BASE_DEFAULT}
    CACHE STRING "Base URL to Trilinos repos <url-base><repo-name>")
ENDIF()
MARK_AS_ADVANCED(Trilinos_REPOS_URL_BASE)

TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES(
  MOOCHO_repo  packages/moocho  GIT  ${Trilinos_REPOS_URL_BASE}moocho  NOPACKAGES  Nightly
  CTrilinos_repo  packages/CTrilinos  GIT  ${Trilinos_REPOS_URL_BASE}CTrilinos  NOPACKAGES  Nightly
  ForTrilinos_repo  packages/ForTrilinos  GIT  ${Trilinos_REPOS_URL_BASE}ForTrilinos  NOPACKAGES  EX
  Optika_repo  packages/optika  GIT  ${Trilinos_REPOS_URL_BASE}optika  NOPACKAGES  Nightly
  Mesquite_repo  packages/mesquite  GIT  ${Trilinos_REPOS_URL_BASE}mesquite  NOPACKAGES  Nightly
  preCopyrightTrilinos  ""  GIT  software.sandia.gov:/space/git/preCopyrightTrilinos  ""  Continuous
  TerminalApplication  ""  GIT  software.sandia.gov:/space/git/TerminalApplication  ""   EX 
  )

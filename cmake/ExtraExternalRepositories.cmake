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
# List of extra external repositories that contain extra Trilinos packages.
#
# Extra repositories must be listed here if they are going to be used in
# Continuous (CI) or Nightly testing driven through the
# TribitsCTestDriverCore.cmake script.
#
# The six fields for each extra external repos are:
#
# NAME: Gives the name of the extra repo (and also the directory name if
# not specified below)
#
# DIR: The relative direcotry location that will be created under the base
# Trilinos base source directory.  This is the name that will be passed to the
# VC tool when cloning/checkingout the repository.
#
# REPOTYPE: The typeof the VC repo (GIT, SVN, or CVS)
#
# REPOURL: This is the URL the extra repository will be cloned (or chekced
# out) from.
#
# PACKSTAT: If the repo has packages, will be "".  If it is NO_PACKAGES, the
# repository does not provide any packges.  A value of NO_PACKAGES is for
# cases where the repo provides code but not add-on Trilinos packages
# (e.g. such as is the case with Dakota used by TriKota).
#
# CATEGORY: This is category of tests where extra repo will be pulled in for.
# Valid categories include:
#
#   Continuous: Continuous integration testing run throughout the development day
#
#   Nightly: Nightly testing run at the end of every development day
#
#   EX: Experimental, not run implicitly for any type of testing, including
#     for 'Experimental' builds.
#
# The reason that you would want to list 'EX' repositories, is that it leaves
# the option to explicitly add the repositories in specialized nightly
# testing.
#
# NOTE: The extra repositories must be listed in assending order of package
# dependencies, just like Trilinos packages in the file
# TrilinosPackages.cmake.  This also means that there can not be any cicular
# dependencies between groups of packages within an extra repository.
#
# NOTE: The packages in a downstream extra repo can depend on packages in an
# upstream repo.
#

SET( Trilinos_EXTRAREPOS_DIR_REPOTYPE_REPOURL_PACKSTAT_CATEGORY
  preCopyrightTrilinos  ""  GIT  software.sandia.gov:/space/git/preCopyrightTrilinos  ""  Continuous
  TerminalApplication  ""  GIT  software.sandia.gov:/space/git/TerminalApplication  ""   EX 
  )

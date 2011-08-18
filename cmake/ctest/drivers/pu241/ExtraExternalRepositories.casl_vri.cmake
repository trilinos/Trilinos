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
# List of extra external repositories for CASL.
#
# See documentation in Trilinos preCopyrightTrilinos/ExtraExternalRepositories.cmake
#

SET( Trilinos_EXTRAREPOS_DIR_REPOTYPE_REPOURL_PACKSTAT_CATEGORY
  #Dakota  packages/TriKota/Dakota SVN
  #   https://software.sandia.gov/svn/public/dakota/public/trunk  NOPACKAGES  Nightly  
  StarCCMExt      ""                                  GIT  casl-dev.ornl.gov:/git-root/StarCCMClient  ""  Continuous
  DeCARTExt       ""                                  GIT  casl-dev.ornl.gov:/git-root/casl_decart    ""  Continuous
  Panzer          ""                                  GIT  software.sandia.gov:/space/git/Panzer      ""  Continuous
  NeutronicsExt   "NeutronicsExt"                     GIT  casl-dev.ornl.gov:/git-root/denovoExt      ""  Continuous
  Denovo          "NeutronicsExt/denovo"              GIT  casl-dev.ornl.gov:/git-root/denovo         NOPACKAGES  Continuous
  Nemesis         "NeutronicsExt/denovo/src/nemesis"  GIT  casl-dev.ornl.gov:/git-root/nemesis        NOPACKAGES  Continuous
  CASLBOA         ""                                  GIT  casl-dev.ornl.gov:/git-root/casl_boa       ""  Continuous
  CASLRAVE        ""                                  GIT  casl-dev.ornl.gov:/git-root/casl_rave      ""  Continuous
  RELAP5Ext       ""                                  GIT  casl-dev.ornl.gov:/git-root/casl_relap5    ""  Continuous
  LIMEExt         ""                                  GIT  software.sandia.gov:/space/git/LIMEExt     ""  Continuous
  PSSDriversExt   ""                                  GIT  casl-dev.ornl.gov:/git-root/casl_vripss    ""  Continuous
  )

#!/usr/bin/env python

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
# Read in the commandline arguments
#

usageHelp = r"""dump-package-deps-table.py [OPTIONS]

Tool that dumps an XML file that can be read by CTest/CDash to specify the
Trilinos package dependenices in a way that is independent of Trilinos.  In
CTest/CDash terminology, a Trilinos package is a "Subproject".

By default, if you just run:

   $ SOME_DIR/dump-cdash-deps-xml-file.py

then the XML file will get written into the main Trilinos source directory
where it can be checked in on the next checkin.

You can also change what XML input file is used and what XML file is written.
This is maily to facilitate unit testing of this script.

Have fun looking through all of the Trilinos dependencies!

"""


from TrilinosDependencies import defaultTrilinosDepsXmlInFile, \
  defaultTrilinosDepsHtmlOutFile, defaultCDashDepsXmlFile

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--input-xml-deps-file", dest="inputXmlDepsFile", type="string",
  default=defaultTrilinosDepsXmlInFile,
  help="Input XML file giving the Trilinos dependencies "+\
    "(default = "+defaultTrilinosDepsXmlInFile+")." )

clp.add_option(
  "--output-cdash-deps-xml-file", dest="outputCDashDepsXmlFile", type="string",
  default=defaultCDashDepsXmlFile,
  help="Output XML file giving the Trilinos dependices CDash language"+\
    "(default = "+defaultCDashDepsXmlFile+")." )

(options, args) = clp.parse_args()


#
# Execute the commands
#


from TrilinosDependencies import getTrilinosDependenciesFromXmlFile

trilinosDependencies = getTrilinosDependenciesFromXmlFile(
  options.inputXmlDepsFile)

trilinosDependencies.writeCDashXmlDepsFile(
  options.outputCDashDepsXmlFile)

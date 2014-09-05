// @HEADER
// ***********************************************************************
// 
//         Trilinos: An Object-Oriented Solver Framework 
//                Copyright (2014) Sandia Corporation
// 
// Copyright (2014) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000, there is a non-exclusive license for use of this
// work by or on behalf of the U.S. Government.  Export of this program
// may require a license from the United States Government.
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// NOTICE:  The United States Government is granted for itself and others
// acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
// license in this data to reproduce, prepare derivative works, and
// perform publicly and display publicly.  Beginning five (5) years from
// July 25, 2001, the United States Government is granted for itself and
// others acting on its behalf a paid-up, nonexclusive, irrevocable
// worldwide license in this data to reproduce, prepare derivative works,
// distribute copies to the public, perform publicly and display
// publicly, and to permit others to do so.
//
// NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
// OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
// ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
// RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
// INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
// THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
//
// ************************************************************************
// @HEADER

// Teuchos includes
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_ParameterListExceptions.hpp"
#include "Teuchos_VerboseObject.hpp"
#include <fstream>

// Flatten the namespaces of certain classes
using std::string;
using Teuchos::ParameterList;
using Teuchos::CommandLineProcessor;
using Teuchos::RCP;

int main(int argc, char* argv[])
{

  std::ofstream out;
  out.open("stratimikos.xml", std::ofstream::out);

  Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
  //
  // We will print to standard out, and that output will be valid XML
  // describing the validated ParameterList.  For the purposes of
  // generating nicely-formatted HTML documentation for this
  // ParameterList, we also need to include an XSL header line.  This
  // bool will control whether we include this header line, which can
  // be controlled at the command line.
  //
  bool xsl_header_flag = true;

  //
  // Set up the command line processor.  All versions of this
  // executable should support the add-xsl-header/suppress-xsl-header
  // command line options.  If you want a single executable to support
  // multiple ParameterLists, you could put additional options here to
  // control which ParameterList to output.
  //
  CommandLineProcessor clp(false);  //don't throw exceptions
  clp.recogniseAllOptions(true);
  clp.setOption("add-xsl-header",
		"suppress-xsl-header",
		&xsl_header_flag, 
		"XSL header flag");

  //
  // Parse the command line and quit if not successful
  //
  CommandLineProcessor::EParseCommandLineReturn parse_return =
    clp.parse(argc, argv);
  if(parse_return != CommandLineProcessor::PARSE_SUCCESSFUL)
    return parse_return;

  //
  // Print the XSL header line if requested
  //
  if (xsl_header_flag )
    out << "<?xml-stylesheet type=\"text/xsl\" "
        << "href=\"common/parameterList/parameterList.xsl\"?>\n";

  //
  // Obtain the validated ParameterList and write it to the fancy
  // output stream.  If you wanted to support multiple ParameterLists,
  // this is where the logic would go to choose between them.  Note
  // that Domi has a function that returns the valid ParameterList,
  // but that a more common use case will be to construct a class
  // (that supports the construct-then-init paradigm) with a default
  // constructor and then call its getValidParameters() method.
  //

  Teuchos::writeParameterListToXmlOStream(*linearSolverBuilder.getValidParameters(), out);

}

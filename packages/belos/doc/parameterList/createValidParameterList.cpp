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
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_ParameterListExceptions.hpp"
#include "Teuchos_VerboseObject.hpp"
#include <fstream>

// Add #include statements here to pull in the code that will be
// needed to generate one or more validated ParameterLists
#include "Tpetra_MultiVector.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosSolverFactory.hpp"

// Added for self-registration of the required solvers for the new DII setup
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"

typedef Tpetra::MultiVector<> multivector_type;
typedef Tpetra::Operator<> operator_type;
typedef multivector_type::scalar_type scalar_type;
typedef Belos::SolverManager<scalar_type, multivector_type, operator_type> solver_type;

// Flatten the namespaces of certain classes
using std::string;
using Teuchos::ParameterList;
using Teuchos::CommandLineProcessor;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::parameterList;

void writeXSLHeader( std::ofstream& out )
{
      out << "<?xml-stylesheet type=\"text/xsl\" "
          << "href=\"common/parameterList/parameterList.xsl\"?>\n";
}

// This test is a special case with no libs to provide the new DII registration.
// SolverFactoryParent is now going to run a specialized form of
// registerSolverFactoryForLib based on the template parameters. For this test
// that will pick up Tpetra and attempt to call the following method:
// Belos::Details::Tpetra::registerSolverFactory()
// For this particular test that won't be available so I have manually recreated
// it here and then register just the managers actually used in this test.
// TODO: Decide how we want to handle this  in the new DII setup.
namespace Belos {
namespace Details {
namespace Tpetra {
  // Note this is going to be called premain due to SolverFactoryParent having
  // a link to this method via its constructor and a specialized form of
  // registerSolverFactoryForLib().
  void registerSolverFactory() {
    Belos::Impl::registerSolverSubclassForTypes<
        Belos::PseudoBlockGmresSolMgr<scalar_type, multivector_type, operator_type>,
        scalar_type, multivector_type, operator_type> ("BLOCK GMRES");
    Belos::Impl::registerSolverSubclassForTypes<
        Belos::PseudoBlockGmresSolMgr<scalar_type, multivector_type, operator_type>,
        scalar_type, multivector_type, operator_type> ("PSEUDOBLOCK GMRES");
    Belos::Impl::registerSolverSubclassForTypes<
        Belos::PseudoBlockCGSolMgr<scalar_type, multivector_type, operator_type>,
        scalar_type, multivector_type, operator_type> ("BLOCK CG");
    Belos::Impl::registerSolverSubclassForTypes<
        Belos::PseudoBlockCGSolMgr<scalar_type, multivector_type, operator_type>,
        scalar_type, multivector_type, operator_type> ("PSEUDOBLOCK CG");
  }
} // namespace Tpetra
} // namespace Details
} // namespace Belos

int main(int argc, char* argv[])
{

  bool success = false;
  bool verbose = false;
  try
  {
    //
    // Open an output file stream for writing our XML file
    //
    std::ofstream out;
   
    //
    // We will print to the 'out' ofstream, and that output will be
    // valid XML describing the validated ParameterList.  For the
    // purposes of generating nicely-formatted HTML documentation for
    // this ParameterList, we also need to include an XSL header line.
    // This bool will control whether we include this header line,
    // which can be controlled at the command line.
    //
    bool xsl_header_flag = true;

    //
    // Set up the command line processor.  All versions of this
    // executable should support the add-xsl-header /
    // suppress-xsl-header command line options.  If you want a single
    // executable to support multiple ParameterLists, you could put
    // additional options here to control which ParameterList to
    // output.
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
    // Here is where code should go that is required to generate the
    // validated XML file.  If your class uses a construct-then-init
    // idiom, then this is where the default constructor would be
    // called.  If this executable supports ParameterLists for more
    // than one class, then this is where the logic to pick the
    // appropriate class would go.
    //

	Belos::SolverFactory<scalar_type, multivector_type, operator_type> sFactory;
	RCP<solver_type> solver;
	for ( int i = 0; i < 4; ++i )
	{
		RCP<ParameterList> nullParams = parameterList();
		std::cout << "writing parameter list " << i+1 << " of 9." << std::endl;
		if ( i == 0 ) 
		{
			solver = sFactory.create("Block GMRES", nullParams);
			out.open("belos_BlockGmres.xml", std::ofstream::out);
			if ( xsl_header_flag )
				writeXSLHeader( out );
			RCP<const ParameterList> gmresParams = solver -> getValidParameters();
			Teuchos::writeParameterListToXmlOStream( *gmresParams, out);
			out.close();
		}
		else if ( i == 1 ) 	
		{
			solver = sFactory.create("Pseudo Block GMRES", nullParams);
			out.open("belos_PseudoBlockGmres.xml", std::ofstream::out);
			if ( xsl_header_flag )
				writeXSLHeader( out );
			RCP<const ParameterList> pseudoGmresParams = solver -> getValidParameters();
			Teuchos::writeParameterListToXmlOStream( *pseudoGmresParams, out);
			out.close();
		}
		else if ( i == 2 ) 
		{
			solver = sFactory.create("Block CG", nullParams);
			out.open("belos_BlockCG.xml", std::ofstream::out);
			if ( xsl_header_flag )
				writeXSLHeader( out );
			RCP<const ParameterList> blockCgParams = solver -> getValidParameters();
			Teuchos::writeParameterListToXmlOStream( *blockCgParams, out);
			out.close();
		}
		else if ( i == 3 ) 
		{
			solver = sFactory.create("Pseudo Block CG", nullParams);
			out.open("belos_PseudoBlockCG.xml", std::ofstream::out);
			if ( xsl_header_flag )
				writeXSLHeader( out );
			RCP<const ParameterList> pseudoBlockCgParams = solver -> getValidParameters();
			Teuchos::writeParameterListToXmlOStream( *pseudoBlockCgParams, out);
			out.close();
		}
	}

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );

}

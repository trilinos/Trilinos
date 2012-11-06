/*
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Version.hpp"

#include "AlgorithmA.hpp"


//
// Here is a simple driver function that I call over and over to show
// different features of FancyOStream
//

void doAlgorithmStuff( Teuchos::ParameterList *algoParams = 0 )
{

  // Here I just create the algorithm object that derives from VerboseObject.
  // By default, this object will print to *Verbose::getDefaultOStream()
  AlgorithmA algoA;
  if (algoParams)
    algoA.setParameterList(Teuchos::rcp(algoParams,false));
  // Note that here I could change the stream just this object prints to
  // by calling algoA.setOStream(...).
  
  // Now I call the algorithm which will print to its default output stream
  algoA.doAlgorithm();
  
  *algoA.getOStream() << std::endl;

  TEUCHOS_ASSERT(algoA.getParameterList().getRawPtr() == algoParams);

}

//
// Test that static initialization of VerboseObjectBase and VerboseObject works!
//

class TestVerboseObjectBaseInitialization {
public:
  TestVerboseObjectBaseInitialization()
    {
      // Get the verbosity level for AlgorithmA
      Teuchos::EVerbosityLevel verbLevel = Teuchos::VerboseObject<AlgorithmA>::getDefaultVerbLevel();
      TEUCHOS_TEST_FOR_EXCEPT_PRINT(verbLevel!=Teuchos::VERB_DEFAULT,&std::cerr);
      // Print to the default default OStream to make sure that the initialization
      // trick worked!
      *Teuchos::VerboseObjectBase::getDefaultOStream()
        << "\n***\n*** Printing to default OStream before main() even starts!\n***\n\n"
        << std::flush;
    }
};

static TestVerboseObjectBaseInitialization testVerboseObjectBaseInitialization;

//
// Main driver program
//

int main(int argc, char* argv[])
{

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;
  using Teuchos::OSTab;
  using Teuchos::dyn_cast;
  using Teuchos::CommandLineProcessor;

  bool success = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  const int numProcs = Teuchos::GlobalMPISession::getNProc();

  try {

    // Get some commandline options
		CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);
		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    // Here I am just grabbing the default output stream
    RCP<FancyOStream>
      out = VerboseObjectBase::getDefaultOStream();
    // Note that the VerboseObject manages FancyOStream objects and not just
    // std::ostream objects.  This is important to the design and very
    // resonable I think.

    *out << std::endl << Teuchos::Teuchos_Version() << std::endl << std::endl;

    //
    // Now I call doAlgorithmStuff() a bunch of times with different setups to
    // show the different kinds of line prefix options
    //
  
    *out << "\n***\n*** Testing VerboseObject base class use\n***\n";
  
    *out << "\n*** Algorithm output with default formatting\n\n";
    doAlgorithmStuff();
  
    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Algorithm output with no front matter\n\n";
    out->setShowAllFrontMatter(false);
    doAlgorithmStuff();
  
    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Algorithm output with processor ranks\n\n";
    out->setShowAllFrontMatter(false).setShowProcRank(true);
    doAlgorithmStuff();
  
    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Algorithm output with line prefix names\n\n";
    out->setShowAllFrontMatter(false).setShowLinePrefix(true);
    doAlgorithmStuff();
  
    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Algorithm output with tab counts\n\n";
    out->setShowAllFrontMatter(false).setShowTabCount(true);
    doAlgorithmStuff();
  
    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Algorithm output with line prefix names and tab counts\n\n";
    out->setShowAllFrontMatter(false).setShowLinePrefix(true).setShowTabCount(true);
    doAlgorithmStuff();
  
    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Algorithm output with processor ranks and line prefix names\n\n";
    out->setShowAllFrontMatter(false).setShowProcRank(true).setShowLinePrefix(true);
    doAlgorithmStuff();
  
    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Algorithm output with processor ranks and tab counts\n\n";
    out->setShowAllFrontMatter(false).setShowProcRank(true).setShowTabCount(true);
    doAlgorithmStuff();
  
    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Algorithm output with processor ranks, line prefix names, and tab counts\n\n";
    out->setShowAllFrontMatter(false).setShowProcRank(true).setShowLinePrefix(true).setShowTabCount(true);
    doAlgorithmStuff();
  
    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Algorithm output with processor ranks, line prefix names, and tab counts but no output for AlgorithmA\n\n";
    Teuchos::VerboseObject<AlgorithmA>::setDefaultVerbLevel(Teuchos::VERB_NONE);
    out->setShowAllFrontMatter(false).setShowProcRank(true).setShowLinePrefix(true).setShowTabCount(true);
    doAlgorithmStuff();
    Teuchos::VerboseObject<AlgorithmA>::setDefaultVerbLevel(Teuchos::VERB_DEFAULT);

    *out << "\n*** Running the algorithm by setting parameters in the parameter list ...\n";

    Teuchos::ParameterList algoParams("AlgorithmA");

    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Set AlgorithmA verbosity level to extreme through a parameter list\n\n";
    algoParams.sublist("VerboseObject").set("Verbosity Level","extreme");
    algoParams.set("Algo Type","Harry");
    algoParams.set("Algo Tol",0.3);
    doAlgorithmStuff(&algoParams);

    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Set AlgorithmA verbosity level to medium and the output file \"AlgorithmA.out\" through a parameter list\n\n";
    algoParams.sublist("VerboseObject").set("Verbosity Level","medium");
    algoParams.sublist("VerboseObject").set("Output File","AlgorithmA.out");
    algoParams.set("Algo Type","John");
    algoParams.set("Algo Tol",10);
    doAlgorithmStuff(&algoParams);

    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Set AlgorithmA verbosity level to low and the output back to default through a parameter list\n\n";
    algoParams.sublist("VerboseObject").set("Verbosity Level","low");
    algoParams.sublist("VerboseObject").set("Output File","none");
    algoParams.set("Algo Tol","20");
    doAlgorithmStuff(&algoParams);

    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n***\n*** Do some more simple tests to make sure things work correctly\n***\n\n";

    //
    // Now I do some other simple tests just to see that FancyOStream is working
    // correctly
    //

    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1).setShowTabCount(true);
    out->setProcRankAndSize(mpiSession.getRank(),mpiSession.getNProc());
    
    *out << "\n***\n*** Testing basic FancyOStream and OSTab classes\n***\n\n";
    
    *out << "\nThis is very good output\nand I like it a lot!\n";
    *out << "";
    *out << "\n";
    *out << "This should";
    *out << " all be";
    *out << " printed on";
    *out << " the same";
    *out << " line two lines below the above output!\n";
    RCP<FancyOStream>
      out2 = rcp(new FancyOStream(rcp(new std::ostringstream),"  "));
    {
      OSTab tab1(out);
      *out << "This should be indented one tab!\n";
      {
        OSTab tab2(out);
        *out << "This should be indented two tabs!\n";
        *out2 << "This should be indented zero tabs from out2!\n";
        {
          OSTab tab3(out2);
          *out << "This should be indented two tabs!\n";
          *out2 << "This should be indented one tab from out2!\n";
        }
      }
      *out << "This should be indented one tab!\n";
    }
    *out << "This should be indented zero tabs!\n";
    
    *out << std::endl; // This required overflow() to be overridden!

    *out << "\n***\n*** Now outputting the latent output that was sent to out2\n***\n\n"
         << dyn_cast<std::ostringstream>(*out2->getOStream()).str();

    if(success)
      *out << "\nEnd Result: TEST PASSED" << std::endl;
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);
    
  return ( success ? 0 : 1 );
  
}

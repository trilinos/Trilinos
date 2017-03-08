// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//
//  This test tests some of the Anasazi status test classes
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziSVQBOrthoManager.hpp"

#include "AnasaziStatusTestCombo.hpp"
#include "AnasaziStatusTestMaxIters.hpp"
#include "AnasaziStatusTestWithOrdering.hpp"
#include "AnasaziStatusTestResNorm.hpp"
#include "AnasaziStatusTestOutput.hpp"

#include "AnasaziLOBPCG.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "MyMultiVec.hpp"
#include "MyOperator.hpp"

using namespace Teuchos;
using namespace Anasazi;

class get_out : public std::logic_error {
  public: get_out(const string &whatarg) : std::logic_error(whatarg) {}
};

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  bool verbose = false;
  bool success = false;
  try {
    bool debug = false;

    // number of global elements
    const int blockSize = 5;
    const int dim = 99;

    CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("debug","nodebug",&debug,"Print debugging info.");
    if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
      MPI_Finalize();
#endif
      return -1;
    }

    //
    // Issue several useful typedefs;
    typedef double                    ST;
    typedef MultiVec<ST>              MV;
    typedef Operator<ST>              OP;
    typedef ScalarTraits<ST>         SCT;
    typedef SCT::magnitudeType        MT;

    //
    // Create an output manager
    RCP<OutputManager<ST> > printer
      = rcp( new BasicOutputManager<ST>() );
    int verbosity = Errors;
    if (verbose || debug) {
      verbosity += Warnings;
    }
    if (debug) {
      verbosity += Debug;
    }
    printer->setVerbosity( verbosity );
    //
    // Create a sort manager
    RCP< SortManager<MT> > sorter =
      rcp( new BasicSort<MT>("LM") );
    //
    // Create an orthogonalization manager
    RCP< MatOrthoManager<ST,MV,OP> > ortho =
      rcp( new SVQBOrthoManager<ST,MV,OP>() );

    printer->stream(Warnings) << Anasazi_Version() << std::endl << std::endl;

    //
    // Create an identity matrix
    std::vector<ST> diag(dim);
    for (int i=0; i<dim; i++) diag[i] = 1.0;
    RCP<MyOperator<ST> > I = rcp( new MyOperator<ST>(diag) );
    //
    // Create the solution eigenvectors
    std::vector<SCT::magnitudeType> v(blockSize);
    RCP< MyMultiVec<ST> > ivec = rcp( new MyMultiVec<ST>(dim,blockSize) );
    for (int i=0; i<blockSize; i++) (*ivec)(i,i) = 1.0;
    //
    // Create the solution eigenvalues
    //
    std::vector<SCT::magnitudeType> T(blockSize);
    for (int i=0; i<blockSize; i++) T[i] = 1.0;
    //
    // Create the residual vectors
    RCP< MyMultiVec<ST> > R = rcp( new MyMultiVec<ST>(dim,blockSize) );
    //
    // Create an eigenvalue problem
    RCP< Eigenproblem<ST,MV,OP> > problem =
      rcp( new BasicEigenproblem<ST,MV,OP>(I,ivec) );
    problem->setHermitian(true);
    problem->setNEV(blockSize);
    problem->setProblem();

    //
    // Create a StatusTestCombo (this will wrap around all tests)
    StatusTestCombo<ST,MV,OP> stcombo;
    //
    // Create a StatusTestOutput (this will also wrap around all tests)
    StatusTestOutput<ST,MV,OP> stoutput(printer,null,1,Passed+Failed+Undefined);
    //
    // Create a StatusTestMaxIters
    StatusTestMaxIters<ST,MV,OP> stmaxiter(1);
    //
    // Create a StatusTestResNorm
    StatusTestResNorm<ST,MV,OP> stresnorm(SCT::zero());

    //
    // Create a parameter list
    ParameterList pls;
    pls.set("Block Size",blockSize);
    //
    // Create an eigensolver for testing
    LOBPCG<ST,MV,OP> lobpcg(problem,sorter,printer,rcp(&stoutput,false),ortho,pls);
    //
    // Initialize the solver with the exact eigensolution
    {
      LOBPCGState<ST,MV> state;
      state.X = ivec;
      state.R = R;
      state.T = rcp( &T, false );
      lobpcg.initialize(state);
    }


    ////////////////////////////////////////////////////////////////////////////////
    //
    // Perform tests
    //
    try {

      //
      // test StatusTestOutput
      {
        //
        // stoutput currently has null child pointer
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.getChild() != null, get_out, "StatusTestOutput::getChild() should have returned Teuchos::null.");
        //
        // calling checkStatus() with null child pointer should result in a StatusTestError exception
        bool threw_expected_exception;
        try {
          stoutput.checkStatus(&lobpcg);
          threw_expected_exception = false;
        }
        catch (const StatusTestError &ste) {
          threw_expected_exception = true;
        }
        TEUCHOS_TEST_FOR_EXCEPTION( threw_expected_exception == false, get_out, "StatusTestOutput::checkStatus() should have thrown exception.");
      }

      //
      // test StatusTestResNorm
      {
        stoutput.setChild(rcp(&stresnorm,false));
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.getStatus()          != Undefined, get_out, "StatusTestOutput::setChild() should reset status to Undefined.");
        //
        // solver has residual norms == 0 < SCT::prec()
        printer->print(Warnings,"*** StatusTestResNorm: 0 < -prec: Failed.\n");
        stresnorm.setTolerance(-SCT::prec());  // 0 < -prec() == false
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus() != Undefined, get_out, "StatusTestResNorm::setTolerance() should reset status to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg) != Failed, get_out, "StatusTestResNorm::checkStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.getStatus()          != Failed, get_out, "StatusTestOutput::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus()          != Failed, get_out, "StatusTestResNorm::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.howMany() != 0, get_out, "StatusTestResNorm::howMany() should have returned 0.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.whichVecs().size() != 0, get_out, "StatusTestResNorm::whichVecs() should have been empty.");
        stoutput.clearStatus();
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.getStatus()  != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus() != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");

        printer->print(Warnings,"*** StatusTestResNorm: 0 < prec: Passed.\n");
        stresnorm.setTolerance(SCT::prec());  // 0 < prec() == true
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus() != Undefined, get_out, "StatusTestResNorm::setTolerance() should reset status to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg) != Passed, get_out, "StatusTestResNorm::checkStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.getStatus()          != Passed, get_out, "StatusTestOutput::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus()          != Passed, get_out, "StatusTestResNorm::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.howMany() != blockSize, get_out, "StatusTestResNorm::howMany() should have returned blockSize.");
        TEUCHOS_TEST_FOR_EXCEPTION( (int)stresnorm.whichVecs().size() != blockSize, get_out, "StatusTestResNorm::whichVecs() should have had length blockSize.");
        std::vector<int> whch(stresnorm.whichVecs());
        for (int i=0; i<(int)whch.size(); i++) {
          TEUCHOS_TEST_FOR_EXCEPTION( whch[i] != i, get_out, "StatusTestResNorm::howMany() should have contained {0,blockSize-1}.");
        }
        stoutput.clearStatus();
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.getStatus()  != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus() != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
      }

      //
      // test StatusTestMaxIters
      {
        stoutput.setChild(rcp(&stmaxiter,false));
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.getStatus()          != Undefined, get_out, "StatusTestOutput::setChild() should reset status to Undefined.");
        //
        // solver has numIters() == 0
        printer->print(Warnings,"*** StatusTestMaxIters: 0 >= 1: Failed.\n");
        stmaxiter.setMaxIters(1); // 0 >= 1 == false
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus() != Undefined, get_out, "StatusTestMaxIters::setMaxIters() should reset status to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg) != Failed, get_out, "StatusTestMaxIters::checkStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.getStatus()          != Failed, get_out, "StatusTestResNorm::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()          != Failed, get_out, "StatusTestResNorm::getStatus() unexpected return.");
        stoutput.clearStatus();
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.getStatus()  != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus() != Undefined, get_out, "StatusTestMaxIters::clearStatus() should reset status to Undefined.");

        printer->print(Warnings,"*** StatusTestMaxIters: 0 >= 0: Passed.\n");
        stmaxiter.setMaxIters(0); // 0 >= 0 == true
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus() != Undefined, get_out, "StatusTestMaxIters::setMaxIters() should reset status to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg) != Passed, get_out, "StatusTestMaxIters::checkStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.getStatus()          != Passed, get_out, "StatusTestResNorm::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()          != Passed, get_out, "StatusTestResNorm::getStatus() unexpected return.");
        stoutput.clearStatus();
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.getStatus()  != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus() != Undefined, get_out, "StatusTestMaxIters::clearStatus() should reset status to Undefined.");

        printer->print(Warnings,"*** StatusTestMaxIters: 0 < 0: Failed.\n");
        stmaxiter.setNegate(true); // 0 < 0 == false
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus() != Undefined, get_out, "StatusTestMaxIters::setMaxIters() should reset status to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg) != Failed, get_out, "StatusTestMaxIters::checkStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.getStatus()          != Failed, get_out, "StatusTestResNorm::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()          != Failed, get_out, "StatusTestResNorm::getStatus() unexpected return.");
        stoutput.clearStatus();
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.getStatus()  != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus() != Undefined, get_out, "StatusTestMaxIters::clearStatus() should reset status to Undefined.");

        printer->print(Warnings,"*** StatusTestMaxIters: 0 < 1: Passed.\n");
        stmaxiter.setMaxIters(1); // 0 < 1 == true
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus() != Undefined, get_out, "StatusTestMaxIters::setMaxIters() should reset status to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg) != Passed, get_out, "StatusTestMaxIters::checkStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.getStatus()          != Passed, get_out, "StatusTestResNorm::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()          != Passed, get_out, "StatusTestResNorm::getStatus() unexpected return.");
        stoutput.clearStatus();
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.getStatus()  != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus() != Undefined, get_out, "StatusTestMaxIters::clearStatus() should reset status to Undefined.");
      }

      //
      // test StatusTestCombo(AND)
      //
      // also test clearStatus()
      {
        stoutput.setChild(rcp(&stcombo,false));
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.getStatus()          != Undefined, get_out, "StatusTestOutput::setChild() should reset status to Undefined.");
        stcombo.setTests( tuple<RCP<StatusTest<ST,MV,OP> > >(rcp(&stresnorm,false),rcp(&stmaxiter,false)) );
        TEUCHOS_TEST_FOR_EXCEPTION( stcombo.getTests().size() != 2, get_out, "StatusTestCombo::getTests() should have two tests.");
        TEUCHOS_TEST_FOR_EXCEPTION( stcombo.getStatus()    != Undefined, get_out, "StatusTestCombo::setTests() should reset status to Undefined.");
        stcombo.setComboType( stcombo.AND );
        TEUCHOS_TEST_FOR_EXCEPTION( stcombo.getComboType() != stcombo.AND, get_out, "StatusTestCombo::getComboType() should be AND.");
        TEUCHOS_TEST_FOR_EXCEPTION( stcombo.getStatus()    != Undefined, get_out, "StatusTestCombo::setComboType() should reset status to Undefined.");

        stresnorm.setTolerance(SCT::prec());  // 0 < prec() == true (first test)
        stmaxiter.setNegate(false);
        stmaxiter.setMaxIters(0);             // 0 >= 0     == true (second test)
        // test that T & T => T
        printer->print(Warnings,"*** StatusTestCombo(AND): T & T: Passed.\n");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg)   != Passed, get_out, "StatusTestOutput::checkStatus(): unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus()           != Passed, get_out, "StatusTestOutput::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus()           != Passed, get_out, "StatusTestCombo::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus()           != Passed, get_out, "StatusTestResNorm::clearStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()           != Passed, get_out, "StatusTestMaxIters::clearStatus() unexpected return.");
        stoutput.clearStatus();
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus() != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus() != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus() != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus() != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");

        stresnorm.setTolerance(-SCT::prec()); // 0 < -prec() == false (first test)
        // test that F & T => F
        printer->print(Warnings,"*** StatusTestCombo(AND): F & T: Failed.\n");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg)   != Failed, get_out, "StatusTestOutput::checkStatus(): unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus()           != Failed, get_out, "StatusTestOutput::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus()           != Failed, get_out, "StatusTestCombo::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus()           != Failed, get_out, "StatusTestResNorm::clearStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()           != Passed, get_out, "StatusTestMaxIters::clearStatus() unexpected return.");
        stoutput.clearStatus();
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus() != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus() != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus() != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus() != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");

        stresnorm.setTolerance(SCT::prec());  // 0 < prec() == true (first test)
        stmaxiter.setNegate(true);  // 0 < 0 == false   (second test)
        // test that T & F => F
        printer->print(Warnings,"*** StatusTestCombo(AND): T & F: Failed.\n");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg)   != Failed , get_out, "StatusTestOutput::checkStatus(): unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus()           != Failed , get_out, "StatusTestOutput::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus()           != Failed , get_out, "StatusTestCombo::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus()           != Passed , get_out, "StatusTestResNorm::clearStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()           != Failed, get_out, "StatusTestMaxIters::clearStatus() unexpected return.");
        stoutput.clearStatus();
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus() != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus() != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus() != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus() != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");

        stresnorm.setTolerance(-SCT::prec()); // 0 < -prec() == false (first test)
        // test that F & F => F
        printer->print(Warnings,"*** StatusTestCombo(AND): F & F: Failed.\n");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg)   != Failed , get_out, "StatusTestOutput::checkStatus(): unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus()           != Failed , get_out, "StatusTestOutput::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus()           != Failed , get_out, "StatusTestCombo::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus()           != Failed , get_out, "StatusTestResNorm::clearStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()           != Failed , get_out, "StatusTestMaxIters::clearStatus() unexpected return.");
        stoutput.clearStatus();
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus() != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus() != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus() != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus() != Undefined, get_out, "StatusTestOutput::clearStatus() should reset all to Undefined.");
      }

      //
      // test StatusTestCombo(OR)
      {
        stcombo.setComboType( stcombo.OR );
        TEUCHOS_TEST_FOR_EXCEPTION( stcombo.getComboType() != stcombo.OR, get_out, "StatusTestCombo::getComboType() should be OR.");
        TEUCHOS_TEST_FOR_EXCEPTION( stcombo.getStatus()    != Undefined, get_out, "StatusTestCombo::setComboType() should reset status to Undefined.");

        stresnorm.setTolerance(SCT::prec());  // 0 < prec() == true (first test)
        stmaxiter.setNegate(false);
        stmaxiter.setMaxIters(0);             // 0 >= 0     == true (second test)
        // test that T | T => T
        printer->print(Warnings,"*** StatusTestCombo(OR): T | T: Passed.\n");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg)   != Passed, get_out, "StatusTestOutput::checkStatus(): unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus()           != Passed, get_out, "StatusTestOutput::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus()           != Passed, get_out, "StatusTestCombo::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus()           != Passed, get_out, "StatusTestResNorm::clearStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()           != Passed, get_out, "StatusTestMaxIters::clearStatus() unexpected return.");

        stresnorm.setTolerance(-SCT::prec()); // 0 < -prec() == false (first test)
        // test that F | T => T
        printer->print(Warnings,"*** StatusTestCombo(OR): F | T: Passed.\n");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg)   != Passed, get_out, "StatusTestOutput::checkStatus(): unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus()           != Passed, get_out, "StatusTestOutput::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus()           != Passed, get_out, "StatusTestCombo::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus()           != Failed, get_out, "StatusTestResNorm::clearStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()           != Passed, get_out, "StatusTestMaxIters::clearStatus() unexpected return.");

        stresnorm.setTolerance(SCT::prec());  // 0 < prec() == true (first test)
        stmaxiter.setNegate(true);  // 0 < 0 == false   (second test)
        // test that T | F => T
        printer->print(Warnings,"*** StatusTestCombo(OR): T | F: Passed.\n");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg)   != Passed, get_out, "StatusTestOutput::checkStatus(): unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus()           != Passed, get_out, "StatusTestOutput::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus()           != Passed, get_out, "StatusTestCombo::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus()           != Passed, get_out, "StatusTestResNorm::clearStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()           != Failed, get_out, "StatusTestMaxIters::clearStatus() unexpected return.");

        stresnorm.setTolerance(-SCT::prec()); // 0 < -prec() == false (first test)
        // test that F | F => F
        printer->print(Warnings,"*** StatusTestCombo(OR): F | F: Failed.\n");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg)   != Failed, get_out, "StatusTestOutput::checkStatus(): unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus()           != Failed, get_out, "StatusTestOutput::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus()           != Failed, get_out, "StatusTestCombo::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus()           != Failed, get_out, "StatusTestResNorm::clearStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()           != Failed, get_out, "StatusTestMaxIters::clearStatus() unexpected return.");
      }

      //
      // test StatusTestCombo(SEQAND)
      {
        stcombo.setComboType( stcombo.SEQAND );
        TEUCHOS_TEST_FOR_EXCEPTION( stcombo.getComboType() != stcombo.SEQAND, get_out, "StatusTestCombo::getComboType() should be SEQAND.");
        TEUCHOS_TEST_FOR_EXCEPTION( stcombo.getStatus()    != Undefined, get_out, "StatusTestCombo::setComboType() should reset status to Undefined.");

        stresnorm.setTolerance(SCT::prec());  // 0 < prec() == true (first test)
        stmaxiter.setNegate(false);
        stmaxiter.setMaxIters(0);             // 0 >= 0     == true (second test)
        // test that T && T => T
        printer->print(Warnings,"*** StatusTestCombo(SEQAND): T && T: Passed.\n");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg)   != Passed, get_out, "StatusTestOutput::checkStatus(): unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus()           != Passed, get_out, "StatusTestOutput::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus()           != Passed, get_out, "StatusTestCombo::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus()           != Passed, get_out, "StatusTestResNorm::clearStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()           != Passed, get_out, "StatusTestMaxIters::clearStatus() unexpected return.");

        stresnorm.setTolerance(-SCT::prec()); // 0 < -prec() == false (first test)
        // test that F && U => F
        printer->print(Warnings,"*** StatusTestCombo(SEQAND): F && U: Failed.\n");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg)   != Failed, get_out, "StatusTestOutput::checkStatus(): unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus()           != Failed, get_out, "StatusTestOutput::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus()           != Failed, get_out, "StatusTestCombo::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus()           != Failed, get_out, "StatusTestResNorm::clearStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()           != Undefined, get_out, "StatusTestMaxIters::clearStatus() unexpected return.");

        stresnorm.setTolerance(SCT::prec());  // 0 < prec() == true (first test)
        stmaxiter.setNegate(true);  // 0 < 0 == false   (second test)
        // test that T && F => F
        printer->print(Warnings,"*** StatusTestCombo(SEQAND): T && F: Failed.\n");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg)   != Failed, get_out, "StatusTestOutput::checkStatus(): unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus()           != Failed, get_out, "StatusTestOutput::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus()           != Failed, get_out, "StatusTestCombo::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus()           != Passed, get_out, "StatusTestResNorm::clearStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()           != Failed, get_out, "StatusTestMaxIters::clearStatus() unexpected return.");

        stresnorm.setTolerance(-SCT::prec()); // 0 < -prec() == false (first test)
        // test that F && U => F
        printer->print(Warnings,"*** StatusTestCombo(SEQAND): F && U: Failed.\n");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg)   != Failed, get_out, "StatusTestOutput::checkStatus(): unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus()           != Failed, get_out, "StatusTestOutput::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus()           != Failed, get_out, "StatusTestCombo::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus()           != Failed, get_out, "StatusTestResNorm::clearStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()           != Undefined, get_out, "StatusTestMaxIters::clearStatus() unexpected return.");
      }

      //
      // test StatusTestCombo(SEQOR)
      {
        stcombo.setComboType( stcombo.SEQOR );
        TEUCHOS_TEST_FOR_EXCEPTION( stcombo.getComboType() != stcombo.SEQOR, get_out, "StatusTestCombo::getComboType() should be SEQOR.");
        TEUCHOS_TEST_FOR_EXCEPTION( stcombo.getStatus()    != Undefined, get_out, "StatusTestCombo::setComboType() should reset status to Undefined.");

        stresnorm.setTolerance(SCT::prec());  // 0 < prec() == true (first test)
        stmaxiter.setNegate(false);
        stmaxiter.setMaxIters(0);             // 0 >= 0     == true (second test)
        // test that T || U => T
        printer->print(Warnings,"*** StatusTestCombo(SEQOR): T || U: Passed.\n");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg)   != Passed, get_out, "StatusTestOutput::checkStatus(): unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus()           != Passed, get_out, "StatusTestOutput::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus()           != Passed, get_out, "StatusTestCombo::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus()           != Passed, get_out, "StatusTestResNorm::clearStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()           != Undefined, get_out, "StatusTestMaxIters::clearStatus() unexpected return.");

        stresnorm.setTolerance(-SCT::prec()); // 0 < -prec() == false (first test)
        // test that F || T => T
        printer->print(Warnings,"*** StatusTestCombo(SEQOR): F || T: Passed.\n");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg)   != Passed, get_out, "StatusTestOutput::checkStatus(): unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus()           != Passed, get_out, "StatusTestOutput::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus()           != Passed, get_out, "StatusTestCombo::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus()           != Failed, get_out, "StatusTestResNorm::clearStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()           != Passed, get_out, "StatusTestMaxIters::clearStatus() unexpected return.");

        stresnorm.setTolerance(SCT::prec());  // 0 < prec() == true (first test)
        stmaxiter.setNegate(true);  // 0 < 0 == false   (second test)
        // test that T || U => T
        printer->print(Warnings,"*** StatusTestCombo(SEQOR): T || U: Passed.\n");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg)   != Passed, get_out, "StatusTestOutput::checkStatus(): unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus()           != Passed, get_out, "StatusTestOutput::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus()           != Passed, get_out, "StatusTestCombo::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus()           != Passed, get_out, "StatusTestResNorm::clearStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()           != Undefined, get_out, "StatusTestMaxIters::clearStatus() unexpected return.");

        stresnorm.setTolerance(-SCT::prec()); // 0 < -prec() == false (first test)
        // test that F || F => F
        printer->print(Warnings,"*** StatusTestCombo(SEQOR): F || F: Failed.\n");
        TEUCHOS_TEST_FOR_EXCEPTION( stoutput.checkStatus(&lobpcg)   != Failed, get_out, "StatusTestOutput::checkStatus(): unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(  stoutput.getStatus()           != Failed, get_out, "StatusTestOutput::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION(   stcombo.getStatus()           != Failed, get_out, "StatusTestCombo::getStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stresnorm.getStatus()           != Failed, get_out, "StatusTestResNorm::clearStatus() unexpected return.");
        TEUCHOS_TEST_FOR_EXCEPTION( stmaxiter.getStatus()           != Failed, get_out, "StatusTestMaxIters::clearStatus() unexpected return.");
      }
      success = true;
    } // end of try
    catch (const get_out &go) {
      printer->stream(Warnings) << go.what() << std::endl;
    }

    printer->print(Warnings,"\n");

    if (success)
      printer->print(Warnings,"End Result: TEST PASSED\n");
    else
      printer->print(Warnings,"End Result: TEST FAILED\n");
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}

/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
//@HEADER
*/

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_Comm.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_LinearProblem.h"
#include "Epetra_CrsMatrix.h"

#include "AztecOO.h"

#include "ParameterHelper.hpp"

#include "build_problem.hpp"
#include "build_solver.hpp"

void process_command_line(int argc, char*argv[], std::string& xml_file);

int main(int argc, char*argv[])
{
  Teuchos::GlobalMPISession mpisess(&argc,&argv,&std::cout);

  Teuchos::Time timer("total");
  timer.start();

#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  //Just get one parameter from the command-line: the name of an xml file
  //to get parameters from.

  std::string xml_file("xml_params.xml");
  process_command_line(argc, argv, xml_file);

  //Read the contents of the xml file into a ParameterList. That parameter list
  //should specify a matrix-file and optionally which Belos solver to use, and
  //which Tifpack preconditioner to use, etc. If there are sublists of parameters
  //for Belos and Tifpack, those will be passed to the respective destinations
  //from within the build_problem and build_solver functions.

  std::cout << "Every proc reading parameters from xml_file: "
            << xml_file << std::endl;
  Teuchos::ParameterList test_params =
      Teuchos::ParameterXMLFileReader(xml_file).getParameters();

  //The build_problem function is declared in build_problem.hpp.
  //Note that build_problem calls build_precond and sets a preconditioner on the
  //linear-problem, if a preconditioner is specified.

  Teuchos::RCP<Epetra_LinearProblem> problem = build_problem(test_params, Comm);

  //The build_solver function is declared in build_solver.hpp:

  Teuchos::RCP<AztecOO> solver = build_solver(test_params, problem);

  int max_iters = solver->GetAllAztecOptions()[AZ_max_iter];
  double tol    = solver->GetAllAztecParams()[AZ_tol];

  Teuchos::Time prec_time("precond");
  prec_time.start();
  double cond_est = 0;
  int err_code = -1;
  if (solver->GetAllAztecOptions()[AZ_precond] != AZ_user_precond) {
    err_code = solver->ConstructPreconditioner(cond_est);
  }
  prec_time.stop();
  if (Comm.MyPID() == 0 && err_code == 0) {
    std::cout << "Time to compute preconditioner: " << prec_time.totalElapsedTime() << "s"<<std::endl;
  }

  int ret = solver->Iterate(max_iters, tol);

  if (err_code == 0) {
    solver->DestroyPreconditioner();
  }

  int actual_iters = (int)solver->GetAztecStatus()[AZ_its];

  if (Comm.MyPID() == 0) {
    std::cout << "Converged in " << actual_iters << " iterations." << std::endl;
  }

  if (problem->GetLHS()->NumVectors() > 1) {
    throw std::runtime_error("ERROR: MultiVector->NumVectors()>1.");
  }

  Epetra_MultiVector* x = problem->GetLHS();
  Epetra_MultiVector* b = problem->GetRHS();
  Epetra_Operator* A = problem->GetOperator();
  Teuchos::RCP<Epetra_MultiVector> r = Teuchos::rcp(new Epetra_MultiVector(*b));
  A->Apply(*x, *r);
  r->Update(1.0, *b, -1.0);
  double norm = 0;
  r->Norm2(&norm);

  if (Comm.MyPID() == 0) {
    std::cout << "2-Norm of residual vec: " << norm << std::endl;
  }

  //Clean up by destroying heap-allocated objects that are not
  //held by Teuchos::RCP:
  delete A;
  delete x;
  delete b;
  Epetra_Operator* prec = solver->GetPrecOperator();
  delete prec;

  //If the xml file specified a number of iterations to expect, then we will
  //use that as a test pass/fail criteria.

  if (test_params.isParameter("expectNumIters")) {
    int expected_iters = 0;
    helper::GetParameter(test_params, "expectNumIters", expected_iters);
    double eps = 2.e-6;
    if (ret == 0 && actual_iters <= expected_iters && norm < eps) {
      if (Comm.MyPID() == 0) {
        std::cout << "End Result: TEST PASSED" << std::endl;
      }
    }
    else {
      if (Comm.MyPID() == 0) {
        std::cout << "Actual iters("<<actual_iters
           <<") > expected number of iterations ("
              <<expected_iters<<"), or resid-norm(" << norm << ") >= "<<eps <<std::endl;
      }
    }
  }

  timer.stop();
  if (Comm.MyPID() == 0) {
    std::cout << "proc 0 total program time: " << timer.totalElapsedTime() << std::endl;
  }

  return 0;
}

void process_command_line(int argc, char*argv[], std::string& xml_file)
{
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("xml_file", &xml_file, "XML Parameters file");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    throw std::runtime_error("Error parsing command-line.");
  }
}


/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
//  
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
// 
// ***********************************************************************
// 
// @HEADER

*/

#include <sys/types.h>
#include <unistd.h>

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Thyra includes
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"

// include basic Epetra information
#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
   #include "mpi.h"
#else
   #include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Export.h"

// EpetraExt 
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"

// Thyra-Epetra adapter includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_SIMPLEPreconditionerFactory.hpp"
#include "Teko_LSCPreconditionerFactory.hpp"
#include "Teko_StridedEpetraOperator.hpp"
#include "Teko_EpetraBlockPreconditioner.hpp"

// Aztec includes
#include "AztecOO.h"
#include "AztecOO_Operator.h"

// #include <EcUtils++/directory.h>

#include <iostream>
#include <fstream>
#include <cmath>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::null;
using Teuchos::ParameterList;

int main(int argc,char * argv[])
{
   // calls MPI_Init and MPI_Finalize
   Teuchos::GlobalMPISession mpiSession(&argc,&argv);

   // Handles some I/O to the output screen
   RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

   // build global (or serial communicator
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif


   std::string solveName = "Amesos";
   if(argc>1)
      solveName = argv[1];

   std::cout << "Using \"" << solveName << "\" for approximate solve" << std::endl;
 
   // get process information
   int numProc = Comm.NumProc();
   int myPID   = Comm.MyPID();

   std::cout << "MPI_PID = " << myPID << ", UNIX_PID = " << getpid() << std::endl;

   // output garbage
   *out << "Approaching Barrier: proc = " << numProc << ", pid = " << myPID << std::endl;
   Comm.Barrier();

   RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile("solverparams.xml");

   Epetra_Map map(15444,0,Comm);
   Epetra_CrsMatrix * ptrA = 0;
   Epetra_Vector * ptrf = 0;
   Epetra_Vector * ptrx = 0;

   std::cout << "Reading matrix market file" << std::endl;
   EpetraExt::MatrixMarketFileToCrsMatrix("../data/nsjac_test.mm",map,map,map,ptrA);
   EpetraExt::MatrixMarketFileToVector("../data/nsrhs_test.mm",map,ptrf);
   EpetraExt::MatrixMarketFileToVector("../data/nslhs_test.mm",map,ptrx);

   RCP<Epetra_CrsMatrix> A = rcp(ptrA);
   RCP<Epetra_Vector> b = rcp(ptrf);
   RCP<Epetra_Vector> x = rcp(ptrx);

   std::cout << "Building strided operator" << std::endl;
   std::vector<int> vec(2); vec[0] = 2; vec[1] = 1;
   Teuchos::RCP<Teko::Epetra::StridedEpetraOperator> sA 
         = Teuchos::rcp(new Teko::Epetra::StridedEpetraOperator(vec,A));

   double alpha = 0.9;
   RCP<Teko::InverseFactory> inverse = Teko::invFactoryFromParamList(*paramList,solveName);
#if 1
   RCP<Teko::BlockPreconditionerFactory> precFact 
         = rcp(new Teko::NS::SIMPLEPreconditionerFactory(inverse,alpha));
#else
   RCP<Teko::NS::LSCStrategy> precStrat = rcp(new Teko::NS::InvLSCStrategy(inverse));
   RCP<Teko::BlockPreconditionerFactory> precFact = rcp(new Teko::NS::LSCPreconditionerFactory(precStrat));
#endif

   std::cout << "Preconditioner factory built" << std::endl;
   Teko::Epetra::EpetraBlockPreconditioner prec(precFact);
   prec.buildPreconditioner(sA);

   std::cout << "Preconditioner built" << std::endl;

   Epetra_LinearProblem problem(&*A,&*x,&*b);

   // build solver
   std::cout << "Setting solver parameters" << std::endl;
   AztecOO solver(problem);
   solver.SetAztecOption(AZ_solver,AZ_gmres);
   solver.SetAztecOption(AZ_precond,AZ_none);
   solver.SetAztecOption(AZ_kspace,50);
   solver.SetAztecOption(AZ_output,10);
   solver.SetPrecOperator(&prec);

   std::cout << "Solving" << std::endl;
   solver.Iterate(1000,1e-5);

   return 0;
}

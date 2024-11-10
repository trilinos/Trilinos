// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <string>

#include "MockModelEval_A.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

int main(int argc, char *argv[]) {

  int status=0; // 0 = pass, failures are incremented
  bool success=true;

  // Initialize MPI 
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  int Proc=mpiSession.getRank();
#ifdef HAVE_MPI
  MPI_Comm appComm = MPI_COMM_WORLD;
#else
  int appComm=0;
#endif

  using Teuchos::RCP;
  using Teuchos::rcp;

  // Command-line argument for input file
  //char* defaultfile="input_1.xml";
  
  try {

    RCP<EpetraExt::ModelEvaluator> Model = rcp(new MockModelEval_A(appComm));

    // Set input arguments to evalModel call
    EpetraExt::ModelEvaluator::InArgs inArgs = Model->createInArgs();

    RCP<Epetra_Vector> x = rcp(new Epetra_Vector(*(Model->get_x_init())));
    inArgs.set_x(x);

    int num_p = inArgs.Np();     // Number of *vectors* of parameters
    RCP<Epetra_Vector> p1;
    if (num_p > 0) {
      p1 = rcp(new Epetra_Vector(*(Model->get_p_init(0))));
      inArgs.set_p(0,p1);
    }
    int numParams = p1->MyLength(); // Number of parameters in p1 vector

    // Set output arguments to evalModel call
    EpetraExt::ModelEvaluator::OutArgs outArgs = Model->createOutArgs();

    RCP<Epetra_Vector> f = rcp(new Epetra_Vector(x->Map()));
    outArgs.set_f(f);

    int num_g = outArgs.Ng(); // Number of *vectors* of responses
    RCP<Epetra_Vector> g1;
    if (num_g > 0) {
      g1 = rcp(new Epetra_Vector(*(Model->get_g_map(0))));
      outArgs.set_g(0,g1);
    }

    RCP<Epetra_Operator> W_op = Model->create_W();
    outArgs.set_W(W_op);

    RCP<Epetra_MultiVector> dfdp = rcp(new Epetra_MultiVector(
                             *(Model->get_x_map()), numParams));
    outArgs.set_DfDp(0, dfdp);

    RCP<Epetra_MultiVector> dgdp = rcp(new Epetra_MultiVector(g1->Map(), numParams));
    outArgs.set_DgDp(0, 0, dgdp);

    RCP<Epetra_MultiVector> dgdx = rcp(new Epetra_MultiVector(x->Map(), g1->MyLength()));
    outArgs.set_DgDx(0, dgdx);

    // Now, evaluate the model!
    Model->evalModel(inArgs, outArgs);

    // Print out everything
   if (Proc == 0)
    std::cout << "Finished Model Evaluation: Printing everything {Exact in brackets}" 
         << "\n-----------------------------------------------------------------"
         << std::setprecision(9) << std::endl;
    x->Print(std::cout << "\nSolution vector! {3,3,3,3}\n");
    if (num_p>0) p1->Print(std::cout << "\nParameters! {1,1}\n");
    f->Print(std::cout << "\nResidual! {8,5,0,-7}\n");
    if (num_g>0) g1->Print(std::cout << "\nResponses! {2}\n");
    RCP<Epetra_CrsMatrix> W = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_op, true);
    W->Print(std::cout << "\nJacobian! {6 on diags}\n");
    dfdp->Print(std::cout << "\nDfDp sensitivity MultiVector! {-1,0,0,0}{0,-4,-6,-8}\n");
    dgdp->Print(std::cout << "\nDgDp response sensitivity MultiVector!{2,2}\n");
    dgdx->Print(std::cout << "\nDgDx^T response gradient MultiVector! {-2,-2,-2,-2}\n");

    if (Proc == 0)
     std::cout <<
      "\n-----------------------------------------------------------------\n";
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  if (!success) status+=1000;

  if (Proc==0) {
    if (status==0) 
      std::cout << "TEST PASSED" << std::endl;
    else 
      std::cout << "TEST Failed" << std::endl;
  }

  return status;
}

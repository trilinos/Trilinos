//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include <iostream>
#include <sstream>

// FEApp is defined in Trilinos/packages/sacado/example/FEApp
#include "twoD_diffusion_ME.hpp"

// Epetra communicator
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// AztecOO solver
#include "AztecOO.h"

// Stokhos Stochastic Galerkin
#include "Stokhos_Epetra.hpp"

// Timing utilities
#include "Teuchos_TimeMonitor.hpp"

int main(int argc, char *argv[]) {
  int n = 32;
  int num_KL = 2;
  int p = 5;
  bool full_expansion = false;

// Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  int MyPID;

  try {

    {
    TEUCHOS_FUNC_TIME_MONITOR("Total PCE Calculation Time");

    // Create a communicator for Epetra objects
    Teuchos::RCP<Epetra_Comm> Comm;
#ifdef HAVE_MPI
    Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    Comm = Teuchos::rcp(new Epetra_SerialComm);
#endif

    MyPID = Comm->MyPID();

    // Create Stochastic Galerkin basis and expansion
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(num_KL); 
    for (int i=0; i<num_KL; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(p,true));
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
    int sz = basis->size();
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk;
    if (full_expansion)
      Cijk = basis->computeTripleProductTensor(sz);
    else
      Cijk = basis->computeTripleProductTensor(num_KL+1);
    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion = 
      Teuchos::rcp(new Stokhos::AlgebraicOrthogPolyExpansion<int,double>(basis,
									 Cijk));
    std::cout << "Stochastic Galerkin expansion size = " << sz << std::endl;

    // Create application
    Teuchos::RCP<twoD_diffusion_ME> model = 
      Teuchos::rcp(new twoD_diffusion_ME(Comm, n, num_KL));
    
    // Set up stochastic parameters
    Epetra_LocalMap p_sg_map(num_KL, 0, *Comm);
    Teuchos::Array<Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> >sg_p_init(1);
    sg_p_init[0]= Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(basis, p_sg_map));
    for (int i=0; i<num_KL; i++) {
      sg_p_init[0]->term(i,0)[i] = 0.0;
      sg_p_init[0]->term(i,1)[i] = 1.0;
    }

    // Setup stochastic initial guess
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_x_init = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(basis, 
						       *(model->get_x_map())));
    sg_x_init->init(0.0);
    
    // Setup stochastic Galerkin algorithmic parameters
    Teuchos::RCP<Teuchos::ParameterList> sgParams = 
      Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::ParameterList& sgOpParams = 
      sgParams->sublist("SG Operator");
    Teuchos::ParameterList& sgPrecParams = 
      sgParams->sublist("SG Preconditioner");
    if (!full_expansion) {
      sgParams->set("Parameter Expansion Type", "Linear");
      sgParams->set("Jacobian Expansion Type", "Linear");
    }
    sgOpParams.set("Operator Method", "Matrix Free");
    sgPrecParams.set("Preconditioner Method", "Mean-based");
    sgPrecParams.set("Mean Preconditioner Type", "ML");
    Teuchos::ParameterList& precParams = 
      sgPrecParams.sublist("Mean Preconditioner Parameters");
    precParams.set("default values", "SA");
    precParams.set("ML output", 10);
    precParams.set("max levels",5);
    precParams.set("increasing or decreasing","increasing");
    precParams.set("aggregation: type", "Uncoupled");
    precParams.set("smoother: type","ML symmetric Gauss-Seidel");
    precParams.set("smoother: sweeps",2);
    precParams.set("smoother: pre or post", "both");
    precParams.set("coarse: max size", 200);
#ifdef HAVE_ML_AMESOS
    precParams.set("coarse: type","Amesos-KLU");
#else
    precParams.set("coarse: type","Jacobi");
#endif

    // Create stochastic Galerkin model evaluator
    Teuchos::RCP<Stokhos::SGModelEvaluator> sg_model =
      Teuchos::rcp(new Stokhos::SGModelEvaluator(model, basis, Teuchos::null,
						 expansion, Cijk, sgParams,
						 Comm, sg_x_init, sg_p_init));

    // Create vectors and operators
    Teuchos::RCP<const Epetra_Vector> sg_p = sg_model->get_p_init(1);
    Teuchos::RCP<Epetra_Vector> sg_x = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_x_map())));
    *sg_x = *(sg_model->get_x_init());
    Teuchos::RCP<Epetra_Vector> sg_f = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_f_map())));
    Teuchos::RCP<Epetra_Vector> sg_dx = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_x_map())));
    Teuchos::RCP<Epetra_Operator> sg_J = sg_model->create_W();
    Teuchos::RCP<Epetra_Operator> sg_M = sg_model->create_WPrec()->PrecOp;

    
    // Setup InArgs and OutArgs
    EpetraExt::ModelEvaluator::InArgs sg_inArgs = sg_model->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs sg_outArgs = sg_model->createOutArgs();
    sg_inArgs.set_p(1, sg_p);
    sg_inArgs.set_x(sg_x);
    sg_outArgs.set_f(sg_f);
    sg_outArgs.set_W(sg_J);
    sg_outArgs.set_WPrec(sg_M);

    // Evaluate model
    sg_model->evalModel(sg_inArgs, sg_outArgs);

    // Print initial residual norm
    double norm_f;
    sg_f->Norm2(&norm_f);
    std::cout << "\nInitial residual norm = " << norm_f << std::endl;

    // Setup AztecOO solver
    AztecOO aztec;
    aztec.SetAztecOption(AZ_solver, AZ_gmres);
    aztec.SetAztecOption(AZ_precond, AZ_none);
    aztec.SetAztecOption(AZ_kspace, 20);
    aztec.SetAztecOption(AZ_conv, AZ_r0);
    aztec.SetAztecOption(AZ_output, 1);
    aztec.SetUserOperator(sg_J.get());
    aztec.SetPrecOperator(sg_M.get());
    aztec.SetLHS(sg_dx.get());
    aztec.SetRHS(sg_f.get());

    // Solve linear system
    aztec.Iterate(1000, 1e-12);

    // Update x
    sg_x->Update(-1.0, *sg_dx, 1.0);

    // Compute new residual & response function
    EpetraExt::ModelEvaluator::OutArgs sg_outArgs2 = sg_model->createOutArgs();
    sg_f->PutScalar(0.0);
    sg_outArgs2.set_f(sg_f);
    sg_model->evalModel(sg_inArgs, sg_outArgs2);

    // Print initial residual norm
    sg_f->Norm2(&norm_f);
    std::cout << "\nFinal residual norm = " << norm_f << std::endl;

    if (norm_f < 1.0e-10)
      std::cout << "Test Passed!" << std::endl;

    }

    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();

  }
  
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
  catch (string& s) {
    std::cout << s << std::endl;
  }
  catch (char *s) {
    std::cout << s << std::endl;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" <<std:: endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

}

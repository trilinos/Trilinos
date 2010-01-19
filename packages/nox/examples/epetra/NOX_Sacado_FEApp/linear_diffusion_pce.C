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
#include "FEApp_ModelEvaluator.hpp"

// Epetra communicator
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// AztecOO solver
#include "AztecOO.h"

// Stokhos Stochastic Galerkin
#include "Stokhos.hpp"

// Timing utilities
#include "Teuchos_TimeMonitor.hpp"

int main(int argc, char *argv[]) {
  int nelem = 100;
  double h = 1.0/nelem;
  int num_KL = 3;
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
    
    // Create mesh
    std::vector<double> x(nelem+1);
    for (int i=0; i<=nelem; i++)
      x[i] = h*i;

    // Set up application parameters
    Teuchos::RCP<Teuchos::ParameterList> appParams = 
      Teuchos::rcp(new Teuchos::ParameterList);

    // Problem
    Teuchos::ParameterList& problemParams = 
      appParams->sublist("Problem");
    problemParams.set("Name", "Heat Nonlinear Source");

    // Boundary conditions
    problemParams.set("Left BC", 0.0);
    problemParams.set("Right BC", 0.0);

    // Source function
    Teuchos::ParameterList& sourceParams = 
      problemParams.sublist("Source Function");
    sourceParams.set("Name", "Constant");
    sourceParams.set("Constant Value", 1.0);

    // Material
    Teuchos::ParameterList& matParams = 
      problemParams.sublist("Material Function");
    matParams.set("Name", "KL Exponential Random Field");
    matParams.set("Mean", 1.0);
    matParams.set("Standard Deviation", 0.5);
    matParams.set("Number of KL Terms", num_KL);
    Teuchos::Array<double> a(1), b(1), L(1);
    a[0] = 0.0; b[0] = 1.0; L[0] = 1.0;
    matParams.set("Domain Lower Bounds", a);
    matParams.set("Domain Upper Bounds", b);
    matParams.set("Correlation Lengths", L);

    // Response functions
    Teuchos::ParameterList& responseParams =
      problemParams.sublist("Response Functions");
    responseParams.set("Number", 1);
    responseParams.set("Response 0", "Solution Average");

    // Free parameters (determinisic, e.g., for sensitivities)
    Teuchos::RefCountPtr< Teuchos::Array<std::string> > free_param_names =
	Teuchos::rcp(new Teuchos::Array<std::string>);
    free_param_names->push_back("Constant Source Function Value");
    
    // Create Stochastic Galerkin basis and expansion
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(num_KL); 
    for (int i=0; i<num_KL; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(p));
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
    appParams->set("Enable Stochastic Galerkin", true);
    appParams->set("Stochastic Galerkin expansion", expansion);
    appParams->set("SG Method", "AD");
    Teuchos::RCP<FEApp::Application> app = 
      Teuchos::rcp(new FEApp::Application(x, Comm, appParams, false));
    
    // Set up stochastic parameters
    Teuchos::Array< Stokhos::VectorOrthogPoly<Epetra_Vector> > sg_params(1);
    Epetra_LocalMap p_sg_map(num_KL, 0, *Comm);
    sg_params[0].reset(basis, Stokhos::EpetraVectorCloner(p_sg_map));
    for (int i=0; i<num_KL; i++) {
      sg_params[0].term(i,0)[i] = 0.0;
      sg_params[0].term(i,1)[i] = 1.0;
    }
    Teuchos::RefCountPtr< Teuchos::Array<std::string> > sg_param_names =
      Teuchos::rcp(new Teuchos::Array<std::string>);
    for (int i=0; i<num_KL; i++) {
      std::stringstream ss;
      ss << "KL Exponential Function Random Variable " << i;
      sg_param_names->push_back(ss.str());
    }
    
    // Create application model evaluator
    Teuchos::RCP<EpetraExt::ModelEvaluator> model = 
      Teuchos::rcp(new FEApp::ModelEvaluator(app, free_param_names,
					     sg_param_names));
    
    // Setup stochastic Galerkin algorithmic parameters
    Teuchos::RCP<Teuchos::ParameterList> sgParams = 
      Teuchos::rcp(&(appParams->sublist("SG Parameters")),false);
    if (!full_expansion) {
      sgParams->set("Parameter Expansion Type", "Linear");
      sgParams->set("Jacobian Expansion Type", "Linear");
    }
    sgParams->set("Jacobian Method", "Matrix Free");
    sgParams->set("Mean Preconditioner Type", "ML");
    Teuchos::ParameterList& precParams = 
      sgParams->sublist("Preconditioner Parameters");
    precParams.set("default values", "DD");

    // Create stochastic Galerkin model evaluator
    Teuchos::Array<int> sg_p_index(1);
    Teuchos::Array<int> sg_g_index(1);
    sg_p_index[0] = 1;
    sg_g_index[0] = 0;
    Teuchos::RCP<Stokhos::SGModelEvaluator> sg_model =
      Teuchos::rcp(new Stokhos::SGModelEvaluator(model, basis, Cijk, sg_p_index,
						 sg_g_index, sgParams,
						 Comm, sg_params));

    // Create vectors and operators
    Teuchos::RCP<const Epetra_Vector> sg_p = sg_model->get_p_init(0);
    Teuchos::RCP<Epetra_Vector> sg_x = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_x_map())));
    *sg_x = *(sg_model->get_x_init());
    Teuchos::RCP<Epetra_Vector> sg_f = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_f_map())));
    Teuchos::RCP<Epetra_Vector> sg_dx = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_x_map())));
    Teuchos::RCP<Epetra_Operator> sg_J = sg_model->create_W();
    Teuchos::RCP<Epetra_Operator> sg_M = sg_model->create_M();

    // Setup InArgs and OutArgs
    EpetraExt::ModelEvaluator::InArgs sg_inArgs = sg_model->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs sg_outArgs = sg_model->createOutArgs();
    sg_inArgs.set_p(0, sg_p);
    sg_inArgs.set_x(sg_x);
    sg_outArgs.set_f(sg_f);
    sg_outArgs.set_W(sg_J);
    sg_outArgs.set_M(sg_M);

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
    aztec.SetAztecOption(AZ_output, 10);
    aztec.SetUserOperator(sg_J.get());
    aztec.SetPrecOperator(sg_M.get());
    aztec.SetLHS(sg_dx.get());
    aztec.SetRHS(sg_f.get());

    // Solve linear system
    aztec.Iterate(100, 1e-12);

    // Update x
    sg_x->Update(-1.0, *sg_dx, 1.0);

    // Compute new residual & response function
    Teuchos::RCP<Epetra_Vector> sg_g = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_g_map(0))));
    EpetraExt::ModelEvaluator::OutArgs sg_outArgs2 = sg_model->createOutArgs();
    sg_outArgs2.set_f(sg_f);
    sg_outArgs2.set_g(0, sg_g);
    sg_model->evalModel(sg_inArgs, sg_outArgs2);

    // Print initial residual norm
    sg_f->Norm2(&norm_f);
    std::cout << "\nFinal residual norm = " << norm_f << std::endl;

    // Print mean and standard deviation
    EpetraExt::BlockVector sg_g_block(View, *(model->get_g_map(0)), *sg_g);
    Stokhos::VectorOrthogPoly<Epetra_Vector> sg_g_vec_poly(
      basis, Stokhos::EpetraVectorCloner(sg_g_block));
    Stokhos::OrthogPolyApprox<int,double> sg_g_poly(basis);
    for (int i=0; i<sz; i++)
      sg_g_poly[i] = sg_g_vec_poly[i][0];
    std::cout << "\nResponse Expansion = " << std::endl;
    std::cout.precision(12);
    sg_g_poly.print(std::cout);
    std::cout << "\nResponse Mean =      " << sg_g_poly.mean() << std::endl;
    std::cout << "Response Std. Dev. = " << sg_g_poly.standard_deviation() 
		<< std::endl;

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

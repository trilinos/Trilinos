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

// NOX
#include "NOX.H"
#include "NOX_Epetra.H"
#include "NOX_Epetra_LinearSystem_SGGS.hpp"
#include "NOX_Epetra_LinearSystem_SGJacobi.hpp"

// Stokhos Stochastic Galerkin
#include "Stokhos.hpp"

// Timing utilities
#include "Teuchos_TimeMonitor.hpp"

//The probability distribution of the random variables.
double uniform_weight(const double& x){
 return 1;
}

//Given the computed sfem solution, computes the mean solution.
void computeMeanSoln(const Epetra_Vector& u, Epetra_Vector& Eofu, Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis){


  const Epetra_Comm& Comm = u.Map().Comm();

  //Roll up solution vector into a multivector containing deterministic components.
  int N_x = u.MyLength()/basis->size();
  int N_xi = basis->size();

  Eofu.PutScalar(0.0);

  Epetra_Map Map(N_x, 0, Comm);
  // Form x and y into block vectors.
  Epetra_MultiVector uBlock(Map,N_xi);

  //Get the triple product tensor.
  Teuchos::RCP< const Stokhos::Sparse3Tensor<int, double> > Cijk =
    basis->computeTripleProductTensor(basis->size());
  double val;
  int i, j;
  int n = Cijk->num_values(0);

  for( int l = 0; l<n; l++){
    Cijk->value(0,l,i,j,val);
    if(i==0 && j == 0) break;
  }
  std::cout << "val = " << val << "\n";
  for(int i = 0; i< N_x; i++){
    Eofu[i] = val*u[i];
  }

}

//Given the computed SFEM solution and the mean solution, computes the variance via <u^2> - <u>^2
void computeVarianceSoln(const Epetra_Vector& u, const Epetra_Vector& Eofu, Epetra_Vector& varOfu, Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis){

  const Epetra_Comm& Comm = u.Map().Comm();

  int N_x = u.MyLength()/basis->size();
  int N_xi = basis->size();


  Epetra_Map Map(N_x, 0, Comm);
  varOfu.Multiply(-1.0, Eofu, Eofu, 0.0);
  Epetra_MultiVector uBlock(Map,N_xi);

  int MyLength = uBlock.MyLength();
  for( int c=0; c<N_xi ; c++){
    for( int i=0; i<MyLength; i++){
      uBlock[c][i] = (u)[c*N_x + i];
    }
  }

  uBlock.Multiply(1.0, uBlock, uBlock, 0.0);
  Teuchos::Array< double > norms = basis->norm_squared();
  for(int c = 0; c<N_xi; c++){
    varOfu.Update(norms[c],*uBlock(c),1.0);
  }

}

int main(int argc, char *argv[]) {   
Teuchos::Time TotalTimer("Total Timer",false);
TotalTimer.start();
//TotalTimer.stop();
  int n, num_KL, p;
  double sigma, mean, weightCut;
  std::string solve_method, precMethod, randField;
  if(argc < 10){
    n = 32; //Number of mesh points
    p = 5; //Polynomial degree
    num_KL = 3;  //Terms in KL expansion
    sigma = .1;
    mean = .2;
    weightCut = 1;   // Support for distribution is +-weightCut
    solve_method = "SG_GMRES";
    precMethod = "Mean-based";
    //precMethod = "Gauss-Seidel";
    //precMethod = "Approx-Gauss-Seidel";
    randField = "UNIFORM";
    //randField = "LOG-NORMAL";
  }else{
    n = atoi(argv[1]);
    p = atoi(argv[2]);
    num_KL = atoi(argv[3]);
    sigma = atof(argv[4]);
    mean = atof(argv[5]);
    weightCut = atof(argv[6]);
    solve_method = argv[7];
    precMethod = argv[8];
    randField = argv[9];
  }
std::cout<< "sigma = " << sigma << " mean = " << mean << "\n";
 
//  int n = 32;
//  int num_KL = 2; 
//int p = 5;
  bool full_expansion;
  if (randField == "UNIFORM")
    full_expansion = false;
  else if (randField == "LOG-NORMAL")
    full_expansion = true;
  bool matrix_free = true;
  bool scaleOP = true; 
//  bool write_linear_system = false;

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
    
    // Create application
    Teuchos::RCP<twoD_diffusion_ME> model =
      Teuchos::rcp(new twoD_diffusion_ME(Comm, n, num_KL, sigma, 
                                          mean, full_expansion));

    // Set up NOX parameters
//    Teuchos::RCP<Teuchos::ParameterList> noxParams =
  //    Teuchos::rcp(sg_model->sublist("NOX"),false);
    Teuchos::RCP<Teuchos::ParameterList> noxParams = 
      Teuchos::rcp(new Teuchos::ParameterList);

    // Set the nonlinear solver method
    noxParams->set("Nonlinear Solver", "Line Search Based");

    // Set the printing parameters in the "Printing" sublist
    Teuchos::ParameterList& printParams = noxParams->sublist("Printing");
    printParams.set("MyPID", MyPID); 
    printParams.set("Output Precision", 3);
    printParams.set("Output Processor", 0);
    printParams.set("Output Information", 
                    NOX::Utils::OuterIteration + 
                    NOX::Utils::OuterIterationStatusTest + 
                    NOX::Utils::InnerIteration +
		    NOX::Utils::Parameters + 
		    NOX::Utils::Details + 
		    NOX::Utils::LinearSolverDetails +
                    NOX::Utils::Warning + 
                    NOX::Utils::Error);

    // Create printing utilities
    NOX::Utils utils(printParams);

    // Sublist for line search 
    Teuchos::ParameterList& searchParams = noxParams->sublist("Line Search");
    searchParams.set("Method", "Full Step");

    // Sublist for direction
    Teuchos::ParameterList& dirParams = noxParams->sublist("Direction");
    dirParams.set("Method", "Newton");
    Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");

    // Sublist for linear solver for the Newton method
    Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
    // Sublist for convergence tests
    Teuchos::ParameterList& statusParams = noxParams->sublist("Status Tests");
    statusParams.set("Test Type", "Combo");
    statusParams.set("Number of Tests", 2);
    statusParams.set("Combo Type", "OR");
    Teuchos::ParameterList& normF = statusParams.sublist("Test 0");
    normF.set("Test Type", "NormF");
    normF.set("Tolerance", 1e-12);
    normF.set("Scale Type", "Scaled");
    Teuchos::ParameterList& maxIters = statusParams.sublist("Test 1");
    maxIters.set("Test Type", "MaxIters");
    maxIters.set("Maximum Iterations", 1);

    Teuchos::ParameterList det_lsParams;
    det_lsParams.set("Aztec Solver", "GMRES");  
    det_lsParams.set("Max Iterations", 5000);
    det_lsParams.set("Size of Krylov Subspace", 100);
    det_lsParams.set("Tolerance", 3e-13); 
   // det_lsParams.set("Max Iterations", 1);
  //  det_lsParams.set("Size of Krylov Subspace", 1);
    //det_lsParams.set("Tolerance", 1e-4); 
    det_lsParams.set("Output Frequency", 0);
    det_lsParams.set("Preconditioner", "ML");    
    det_lsParams.set("Zero Initial Guess", true);    
    Teuchos::ParameterList& det_ML = 
      det_lsParams.sublist("ML");
    ML_Epetra::SetDefaults("SA", det_ML);
    det_ML.set("ML output", 0);
    det_ML.set("max levels",5);
    det_ML.set("increasing or decreasing","increasing");
    det_ML.set("aggregation: type", "Uncoupled");
    det_ML.set("smoother: type","ML symmetric Gauss-Seidel");
    det_ML.set("smoother: sweeps",2);
    det_ML.set("smoother: pre or post", "both");
    det_ML.set("coarse: max size", 200);
#ifdef HAVE_ML_AMESOS
    det_ML.set("coarse: type","Amesos-KLU");
#else
    det_ML.set("coarse: type","Jacobi");
#endif
   
//   det_ML.set("ML output", 10);
//    det_lsParams.set("Write Linear System", write_linear_system);

    // Create NOX interface
    Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> det_nox_interface = 
       Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(model));

     // Create NOX linear system object
    Teuchos::RCP<const Epetra_Vector> det_u = model->get_x_init();
    Teuchos::RCP<Epetra_Operator> det_A = model->create_W();
    Teuchos::RCP<NOX::Epetra::Interface::Required> det_iReq = det_nox_interface;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> det_iJac = det_nox_interface;
    Teuchos::ParameterList det_printParams;
    det_printParams.set("MyPID", MyPID); 
    det_printParams.set("Output Precision", 3);
    det_printParams.set("Output Processor", 0);
    det_printParams.set("Output Information", NOX::Utils::Error);
    Teuchos::RCP<NOX::Epetra::LinearSystem> det_linsys =
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(det_printParams, 
							det_lsParams,
                                                        det_iReq, 
							det_iJac, 
							det_A,
                                                        *det_u));

    // Create Stochastic Galerkin basis and expansion
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(num_KL); 
    for (int i=0; i<num_KL; i++)
      if (randField == "UNIFORM")
        bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(p,true));
      else if (randField == "LOG-NORMAL")      
        bases[i] = Teuchos::rcp(new Stokhos::HermiteBasis<int,double>(p,true));

    //  bases[i] = Teuchos::rcp(new Stokhos::DiscretizedStieltjesBasis<int,double>("beta",p,&uniform_weight,-weightCut,weightCut,true));
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

    if (!full_expansion) {
      sgParams->set("Parameter Expansion Type", "Linear");
      sgParams->set("Jacobian Expansion Type", "Linear");
    }
    sgParams->set("Jacobian Method", "Matrix Free");
    if (precMethod == "Mean-based")  {
      sgParams->set("Preconditioner Method", "Mean-based");
      sgParams->set("Mean Preconditioner Type", "ML");
      Teuchos::ParameterList& precParams =
	sgParams->sublist("Preconditioner Parameters");
      precParams = det_ML;
    }
    else if(precMethod == "Gauss-Seidel") {
      sgParams->set("Preconditioner Method", "Gauss-Seidel");
      Teuchos::ParameterList& GS_params = sgParams->sublist("Gauss-Seidel");
      GS_params.sublist("Deterministic Krylov Solver") = det_lsParams;
      GS_params.set("Deterministic Solver", det_linsys);
      GS_params.set("Max Iterations", 2);
      GS_params.set("Tolerance", 3e-13);
      GS_params.set("Save MatVec Table", false);
    }
    else if (precMethod == "Approx-Gauss-Seidel")  {
      sgParams->set("Preconditioner Method", "Approximate Gauss-Seidel");
      sgParams->set("Symmetric Gauss-Seidel", false);
      sgParams->set("Mean Preconditioner Type", "ML");
      Teuchos::ParameterList& precParams =
	sgParams->sublist("Preconditioner Parameters");
      precParams = det_ML;
    }
    else if (precMethod == "Approx-Jacobi")  {
      sgParams->set("Preconditioner Method", "Approximate Jacobi");
      sgParams->set("Symmetric Gauss-Seidel", false);
      sgParams->set("Mean Preconditioner Type", "ML");
      Teuchos::ParameterList& precParams =
        sgParams->sublist("Preconditioner Parameters");
      precParams = det_ML;
    }
    else
      sgParams->set("Preconditioner Method", "Jacobi");

   // Create stochastic Galerkin model evaluator
    Teuchos::RCP<Stokhos::SGModelEvaluator> sg_model =
      Teuchos::rcp(new Stokhos::SGModelEvaluator(model, basis, Teuchos::null,
                                                 expansion, Cijk, sgParams,
                                                 Comm, sg_x_init, sg_p_init,scaleOP));

     // Create NOX interface
    Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> nox_interface =
       Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(sg_model));

    // Parameter list for SG Jacobi iterations 
    lsParams.set("Max Iterations", 5000);
    lsParams.set("Tolerance", 1e-12);

    // Create NOX stochastic linear system object
    Teuchos::RCP<const Epetra_Vector> u = sg_model->get_x_init();
    Teuchos::RCP<Epetra_Operator> A = sg_model->create_W();
    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = nox_interface;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = nox_interface;

    // Build NOX group
    Teuchos::RCP<NOX::Epetra::Group> grp;
 
    if (solve_method=="SG_GMRES") {
      lsParams.set("Aztec Solver", "GMRES");
      lsParams.set("Max Iterations", 5000);
      lsParams.set("Size of Krylov Subspace", 100);
      lsParams.set("Tolerance", 1e-12); 
      lsParams.set("Output Frequency", 1);
      lsParams.set("Preconditioner", "ML");
      lsParams.set("Zero Initial Guess", true);    
      Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linsys;
    if (matrix_free) {
      Teuchos::RCP<Epetra_Operator> M = sg_model->create_WPrec()->PrecOp;
      Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = nox_interface;
      lsParams.set("Preconditioner", "User Defined");
      linsys = 
	Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
							  iJac, A, iPrec, M,
							  *u));
    }
    else {
      linsys = 
	Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
							  iReq, iJac, A, 
							  *u));
    }
    grp = Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, *u, linsys));
    std::cout << "Solver is SG GMRES " << std::endl;
   }
  else if (solve_method=="SG_GS") {
    det_lsParams.set("Tolerance", 3e-13); 
    lsParams.sublist("Deterministic Krylov Solver") = det_lsParams;
    Teuchos::RCP<NOX::Epetra::LinearSystemSGGS> linsys =
      Teuchos::rcp(new NOX::Epetra::LinearSystemSGGS(printParams, lsParams,
                                                         det_linsys, Cijk,
                                                         iReq, iJac, A, *u, *det_u));
    grp = Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, *u, linsys));
    std::cout << "Solver is SG GS " << std::endl;
  }
  else {
    det_lsParams.set("Tolerance", 3e-13); 
    lsParams.sublist("Deterministic Krylov Solver") = det_lsParams;
    Teuchos::RCP<NOX::Epetra::LinearSystemSGJacobi> linsys =
      Teuchos::rcp(new NOX::Epetra::LinearSystemSGJacobi(printParams, lsParams,
                                                         det_linsys, Cijk,
                                                          iReq, iJac, A, *u, *det_u));
    grp = Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, *u, linsys));
    std::cout << "Solver is SG JACOBI " << std::endl;
 }

//    // Build NOX group
//    Teuchos::RCP<NOX::Epetra::Group> grp = 
//      Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, *u, linsys));

    // Create the Solver convergence test
    Teuchos::RCP<NOX::StatusTest::Generic> statusTests =
      NOX::StatusTest::buildStatusTests(statusParams, utils);

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver = 
      NOX::Solver::buildSolver(grp, statusTests, noxParams);

    Teuchos::Time SolutionTimer("Total Timer",false);
    SolutionTimer.start();
    // Solve the system
    NOX::StatusTest::StatusType status = solver->solve();
    SolutionTimer.stop();

    // Get final solution
    const NOX::Epetra::Group& finalGroup = 
      dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
    const Epetra_Vector& finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

//    std::cout << "finalSolution" << finalSolution <<std::endl;

//////////////////////////////////////////////////////////////////////
//Post process and output the results.
////////////////////////////////////////////////////////////////////
Epetra_Vector Eofu(*(model->get_x_map()),true);
Epetra_Vector varOfu(*(model->get_x_map()),true);
Epetra_Vector specificSol(*(model->get_x_map()),true);
computeMeanSoln(finalSolution, Eofu, basis);
computeVarianceSoln(finalSolution, Eofu, varOfu, basis);
TotalTimer.stop();

std::ofstream mean;
mean.open("mean.txt");
Eofu.Print(mean);
mean.close();

std::ofstream var;
var.open("var.txt");
varOfu.Print(var);
var.close();

std::ofstream time;
time.open("gal_time.txt");
//Teuchos::ParameterList ML_Output = det_ML.GetOutputList();
//double precon_time = ML_Output.get("time: total apply",-2000.0);
time << TotalTimer.totalElapsedTime(false) << "\n";
time << SolutionTimer.totalElapsedTime(false) << "\n";
//time << precon_time << "\n";
//time << system.ApplyTime()<< "\n";
time.close();

std::ofstream dof;
dof.open("gal_dof.txt");
dof << basis->size() << "\n";
dof.close();

/*int iters = aztec_solver.NumIters();
std::ofstream iters_file;
iters_file.open("gal_iters.txt");
iters_file << iters << "\n";
iters_file.close();
*/
      
/*    // Evaluate SG responses at SG parameters
    EpetraExt::ModelEvaluator::InArgs sg_inArgs = sg_model->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs sg_outArgs = 
      sg_model->createOutArgs();
    Teuchos::RCP<const Epetra_Vector> sg_p = sg_model->get_p_init(1);
    Teuchos::RCP<Epetra_Vector> sg_g = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_g_map(1))));
    sg_inArgs.set_p(1, sg_p);
    sg_inArgs.set_x(Teuchos::rcp(&finalSolution,false));
    sg_outArgs.set_g(1, sg_g);
    sg_model->evalModel(sg_inArgs, sg_outArgs);

    // Print mean and standard deviation
    Stokhos::EpetraVectorOrthogPoly sg_g_poly(basis, View, 
					      *(model->get_g_map(0)), *sg_g);
    Epetra_Vector mean(*(model->get_g_map(0)));
    Epetra_Vector std_dev(*(model->get_g_map(0)));
    sg_g_poly.computeMean(mean);
    sg_g_poly.computeStandardDeviation(std_dev);
    std::cout << "\nResponse Expansion = " << std::endl;
    std::cout.precision(12);
    sg_g_poly.print(std::cout);
    std::cout << "\nResponse Mean =      " << std::endl << mean << std::endl;
    std::cout << "Response Std. Dev. = " << std::endl << std_dev << std::endl;
  */    

    if (status == NOX::StatusTest::Converged) 
      utils.out() << "Test Passed!" << std::endl;

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

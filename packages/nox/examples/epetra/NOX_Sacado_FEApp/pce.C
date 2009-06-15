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

#include "NOX.H"
#include "NOX_Epetra.H"

// FEApp is defined in Trilinos/packages/sacado/example/FEApp
#include "FEApp_ModelEvaluator.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Stokhos.hpp"
#include "Stokhos_SGModelEvaluator.hpp"
#include "Stokhos_SGQuadModelEvaluator.hpp"
#include "Sacado_PCE_OrthogPoly.hpp"
#include "Ifpack.h"
#include "Teuchos_TimeMonitor.hpp"

// The preconditioner we will use for PCE
class IfpackPreconditionerFactory : public Stokhos::PreconditionerFactory {
public:
  IfpackPreconditionerFactory(const Teuchos::RCP<Teuchos::ParameterList>& p) :
    precParams(p) {}
  virtual ~IfpackPreconditionerFactory() {}
  virtual Teuchos::RCP<Epetra_Operator> 
  compute(const Teuchos::RCP<Epetra_Operator>& op) {
    Teuchos::RCP<Epetra_RowMatrix> mat = 
      Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(op, true);
    Ifpack Factory;
    std::string prec = precParams->get("Ifpack Preconditioner", "ILU");
    int overlap = precParams->get("Overlap", 0);
    ifpackPrec = Teuchos::rcp(Factory.Create(prec, mat.get(), overlap));
    ifpackPrec->SetParameters(*precParams);
    int err = ifpackPrec->Initialize();   
    err = ifpackPrec->Compute();
    return ifpackPrec;
  }
protected:
  Teuchos::RCP<Teuchos::ParameterList> precParams;
  Teuchos::RCP<Ifpack_Preconditioner> ifpackPrec;
};


int main(int argc, char *argv[]) {
  unsigned int nelem = 100;
  double h = 1.0/nelem;
  double alpha = 1.0;
  double leftBC = 0.0;
  double rightBC = 0.1;
  //double rightBC = 0.0;
  unsigned int numalpha = 2;

  bool do_pce = true;
  bool do_dakota = false;

  int MyPID;

  try {

    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
#endif

    // Create a communicator for Epetra objects
    Teuchos::RCP<Epetra_Comm> Comm;
#ifdef HAVE_MPI
    Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    Comm = Teuchos::rcp(new Epetra_SerialComm);
#endif

    // Get filenames for Dakota runs
    if (argc > 1 && argc != 3) {
      std::cout << "Usage:  pce.exe [input_file] [output file]" << std::endl;
      exit(-1);
    }

    std::string input_filename, output_filename;
    if (argc == 3) {
      input_filename = std::string(argv[1]);
      output_filename = std::string(argv[2]);
      do_pce = false;
      do_dakota = true;
    }

    MyPID = Comm->MyPID();
    
    // Create mesh
    vector<double> x(nelem+1);
    for (unsigned int i=0; i<=nelem; i++)
      x[i] = h*i;

    // Set up application parameters
    Teuchos::RCP<Teuchos::ParameterList> appParams = 
      Teuchos::rcp(new Teuchos::ParameterList);

    // Problem
    Teuchos::ParameterList& problemParams = 
      appParams->sublist("Problem");
    problemParams.set("Name", "Heat Nonlinear Source");

    // Boundary conditions
    problemParams.set("Left BC", leftBC);
    problemParams.set("Right BC", rightBC);

    // Source function
    Teuchos::ParameterList& sourceParams = 
      problemParams.sublist("Source Function");
    sourceParams.set("Name", "Multi-Variate Exponential");
    sourceParams.set("Nonlinear Factor Dimensions", numalpha);
    for (unsigned int i=0; i<numalpha; i++) {
      std::stringstream ss;
      ss << "Nonlinear Factor " << i;
      sourceParams.set(ss.str(), alpha/numalpha);
    }

    // Material
    Teuchos::ParameterList& matParams = 
      problemParams.sublist("Material Function");
    matParams.set("Name", "Constant");
    matParams.set("Constant Value", 1.0);

    // Response functions
    Teuchos::ParameterList& responseParams =
      problemParams.sublist("Response Functions");
    responseParams.set("Number", 1);
    responseParams.set("Response 0", "Solution Average");
    
    // Read in parameter values from input file
    if (do_dakota) {
      std::ifstream input_file(input_filename.c_str());
      int nvar;
      std::string name;
      input_file >> nvar >> name;
      std::vector<double> vals(nvar);
      for (int i=0; i<nvar; i++) {
        input_file >> vals[i] >> name;
        std::stringstream ss;
        ss << "Nonlinear Factor " << i;
        sourceParams.set(ss.str(), vals[i]);
      }
      input_file.close();
    }

    Teuchos::RefCountPtr< Teuchos::Array<std::string> > free_param_names =
	Teuchos::rcp(new Teuchos::Array<std::string>);
    free_param_names->push_back("Constant Function Value");

    // Create application
    Teuchos::RCP<FEApp::Application> app = 
      Teuchos::rcp(new FEApp::Application(x, Comm, appParams, false));

    // Create initial guess
    Teuchos::RCP<const Epetra_Vector> u = app->getInitialSolution();
    NOX::Epetra::Vector nox_u(*u);

    // Create model evaluator
    Teuchos::RCP<EpetraExt::ModelEvaluator> model = 
      Teuchos::rcp(new FEApp::ModelEvaluator(app, free_param_names));

    // Create NOX interface
    Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> interface =
      Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(model));

    // Set up NOX parameters
    Teuchos::RCP<Teuchos::ParameterList> noxParams =
      Teuchos::rcp(&(appParams->sublist("NOX")),false);

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
    lsParams.set("Aztec Solver", "GMRES");  
    lsParams.set("Max Iterations", 20);
    lsParams.set("Size of Krylov Subspace", 20);
    lsParams.set("Tolerance", 1e-4); 
    lsParams.set("Output Frequency", 50);
    lsParams.set("Preconditioner", "Ifpack");
    lsParams.set("RCM Reordering", "Enabled");

    // Create the Jacobian matrix
    Teuchos::RCP<Epetra_Operator> A = model->create_W(); 

    // Create the linear system
    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linsys = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, 
                                                        lsParams,
                                                        iReq, iJac, A, nox_u));
    
    // Create the Group
    Teuchos::RCP<NOX::Epetra::Group> grp =
      Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, nox_u, linsys)); 

    // Create the Solver convergence test
    Teuchos::RCP<NOX::StatusTest::NormF> wrms = 
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-10, 
					   NOX::StatusTest::NormF::Unscaled));
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = 
      Teuchos::rcp(new NOX::StatusTest::MaxIters(10));
    Teuchos::RCP<NOX::StatusTest::Combo> combo = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(wrms);
    combo->addStatusTest(maxiters);

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver = 
      NOX::Solver::buildSolver(grp, combo, noxParams);

    // Solve
    NOX::StatusTest::StatusType status = solver->solve();
    if (status == NOX::StatusTest::Converged) 
      utils.out() << "Nonlinear solver converged!" << endl;
    else {
      utils.out() << "Nonlinear solver failed to converge!" << endl;
    }

    // Get the Epetra_Vector with the final solution from the solver
    const NOX::Epetra::Group& finalGroup = 
      dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
    Epetra_Vector finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

    // Compute our objective function
    double norm_u = 0.0;
    finalSolution.MeanValue(&norm_u);

    // Print objective function to file for Dakota
    if (do_dakota) {
      std::ofstream output_file(output_filename.c_str());
      output_file.precision(12);
      output_file.setf(ios::scientific);
      output_file << norm_u << " " << norm_u << std::endl;
      output_file.close();
    }

    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();

    if (do_pce) {

      TEUCHOS_FUNC_TIME_MONITOR("PCE Time");

      unsigned int d = numalpha;
      unsigned int p = 5;
    
      // Create SG basis and expansion
      typedef Stokhos::LegendreBasis<int,double> basis_type;
      std::vector< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d); 
      for (unsigned int i=0; i<d; i++)
        bases[i] = Teuchos::rcp(new basis_type(p));
      Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
        Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
      Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
        Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
      unsigned int sz = basis->size();
      Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion = 
	Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(basis, 
								      quad));
      Sacado::PCE::OrthogPoly<double>::initExpansion(expansion);

      bool sg_local = true;
      if (sg_local) {
	// Create new app for Stochastic Galerkin solve
	appParams->set("Enable Stochastic Galerkin", true);
	appParams->set("Stochastic Galerkin basis", basis);
	appParams->set("Stochastic Galerkin quadrature", quad);
	appParams->set("SG Method", "AD");
	//appParams->set("SG Method", "Gauss Quadrature");
	app = Teuchos::rcp(new FEApp::Application(x, Comm, appParams, false));
      }

      // Set up stochastic parameters
      Teuchos::Array< Teuchos::Array< Teuchos::RCP<Epetra_Vector> > > sg_p(1);
      sg_p[0].resize(sz);
      Epetra_LocalMap p_sg_map(d, 0, *Comm);
      for (unsigned int i=0; i<sz; i++)
	sg_p[0][i] = Teuchos::rcp(new Epetra_Vector(p_sg_map));
      for (unsigned int i=0; i<d; i++)
	(*(sg_p[0][i+1]))[i] = 1.0;
      Teuchos::RefCountPtr< Teuchos::Array<std::string> > sg_param_names =
	Teuchos::rcp(new Teuchos::Array<std::string>);
      for (unsigned int i=0; i<d; i++) {
	std::stringstream ss;
	ss << "Exponential Source Function Nonlinear Factor " << i;
	sg_param_names->push_back(ss.str());
      }
      std::vector<int> sg_p_index(1);
      
      if (sg_local) {
	sg_p_index[0] = 1;
	model = Teuchos::rcp(new FEApp::ModelEvaluator(app, free_param_names,
						       sg_param_names));

      }
      else {
	sg_p_index[0] = 0;
	std::vector<int> sg_p_index_quad(1);
	sg_p_index_quad[0] = 0;
	std::vector<int> sg_g_index_quad(1);
	sg_g_index_quad[0] = 0;
	Teuchos::RCP<EpetraExt::ModelEvaluator> underlying_model = 
	  Teuchos::rcp(new FEApp::ModelEvaluator(app, sg_param_names));
	model =
	  Teuchos::rcp(new Stokhos::SGQuadModelEvaluator(underlying_model, 
							 quad, sg_p_index_quad,
							 sg_g_index_quad));
      }

      Teuchos::RCP<Teuchos::ParameterList> sgParams = 
        Teuchos::rcp(&(appParams->sublist("SG Parameters")),false);
      sgParams->set("Jacobian Method", "Matrix Free");
      //sgParams->set("Jacobian Method", "Fully Assembled");
      Teuchos::RCP<Teuchos::ParameterList> precParams = 
        Teuchos::rcp(&(sgParams->sublist("SG Preconditioner")),false);
      precParams->set("Ifpack Preconditioner", "ILU");
      precParams->set("Overlap", 0);
      Teuchos::RCP<Stokhos::PreconditionerFactory> sg_prec = 
	Teuchos::rcp(new IfpackPreconditionerFactory(precParams));
      sgParams->set("Preconditioner Factory", sg_prec);
      std::vector<int> sg_g_index(1);
      sg_g_index[0] = 0;
      Teuchos::RCP<Stokhos::SGModelEvaluator> sg_model =
	Teuchos::rcp(new Stokhos::SGModelEvaluator(model, basis, sg_p_index,
						   sg_g_index,
						   sg_p, sgParams));
      interface = 
        Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(sg_model));
      iReq = interface;
      iJac = interface;
      Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = interface; 
      Teuchos::RCP<const EpetraExt::BlockVector> sg_init = 
        Teuchos::rcp_dynamic_cast<const EpetraExt::BlockVector>(sg_model->get_x_init());
      EpetraExt::BlockVector sg_u(*sg_init);
      sg_u.LoadBlockValues(finalSolution, 0);
      NOX::Epetra::Vector sg_nox_u(sg_u);
      A = sg_model->create_W(); 
      if (sgParams->get<std::string>("Jacobian Method") == "Matrix Free") {
        Teuchos::RCP<Epetra_Operator> M = sg_model->create_prec();
        lsParams.set("Preconditioner", "User Defined");
        linsys = 
          Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, 
                                                            lsParams,
                                                            iJac, A, 
                                                            iPrec, M,
                                                            sg_nox_u));
      }
      else {
        //lsParams.set("Write Linear System", true);
        linsys = 
          Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, 
                                                            lsParams,
                                                            iReq, iJac, A, 
                                                            sg_nox_u));
      }
      
      grp = Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, sg_nox_u, 
						linsys));

      // Solve SG nonlinear system
      solver = NOX::Solver::buildSolver(grp, combo, noxParams);
      NOX::StatusTest::StatusType status2 = solver->solve();
      if (status2 == NOX::StatusTest::Converged) 
        utils.out() << "SG Nonlinear solver converged!" << endl;
      else
        utils.out() << "SG Nonlinear solver failed to converge!" << endl;
      
      // Get the Epetra_Vector with the final solution from the solver
      const NOX::Epetra::Group& sg_group = 
        dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
      const Epetra_Vector& sg_solution = 
        (dynamic_cast<const NOX::Epetra::Vector&>(sg_group.getX())).getEpetraVector();

      // Compute objective function -- average of u over domain x
      EpetraExt::ModelEvaluator::InArgs inArgs = sg_model->createInArgs();
      Teuchos::RCP<const Epetra_Vector> x_sg = Teuchos::rcp(&sg_solution,false);
      inArgs.set_x(x_sg);
      EpetraExt::ModelEvaluator::OutArgs outArgs = sg_model->createOutArgs();
      Teuchos::RCP<Epetra_Vector> g_sg = 
	Teuchos::rcp(new Epetra_Vector(*(sg_model->get_g_map(0))));
      outArgs.set_g(0, g_sg);
      sg_model->evalModel(inArgs, outArgs);
      g_sg->Print(std::cout);

      // Print mean and standard deviation
      utils.out().precision(12);
      double mean = (*g_sg)[0];
      double std_dev = 0.0;
      const std::vector<double> nrm2 = basis->norm_squared();
      for (int i=1; i<basis->size(); i++)
        std_dev += (*g_sg)[i]*(*g_sg)[i]*nrm2[i];
      std_dev = std::sqrt(std_dev);

      utils.out() << "Mean =      " << mean << std::endl;
      utils.out() << "Std. Dev. = " << std_dev << std::endl;

      if (status2 == NOX::StatusTest::Converged) 
	utils.out() << "Test Passed!" << endl;
    }

    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();

#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif

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

}

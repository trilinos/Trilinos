//@HEADER
// ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#include "Teuchos_MPISession.hpp"
#include "TSFVector.hpp"
#include "TSFLinearCombination.hpp"
#include "TSFVectorType.hpp"
#include "TSFVectorSpace.hpp"
#include "TSFEpetraVectorType.hpp"
#include "TSFNonlinearOperator.hpp"
#include "TSFBICGSTABSolver.hpp"
#include "TSFAztecSolver.hpp"
#include "TSFGhostView.hpp"
#include "TSFGhostImporter.hpp"
#include "Teuchos_Time.hpp"
#include <cmath>
#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_TSF_Group.H"
#include "NOX_Parameter_Teuchos2NOX.H"

using namespace Teuchos;
using namespace TSFExtended;
using namespace TSFExtendedOps;

/*
 * Solves the first-order nonlinear ODE \f$ x'(t)=-\sin(x) \f$using FD on a uniform mesh
 * with initial conditions \f$x(0)=A\f$. In parallel, the FD stencil has overlap so
 * we must use ghost points obtained using the TSF GhostView object.
 *
 * The exact solution is \f$x(t) = 2.0*\tan^{-1} \left[ e^{-t} \tan(A/2) \right]\f$.
 */


namespace TSFExtended
{
  class ODE : public NonlinearOperatorBase<double>
  {
  private:
    double T_;
    double x0_;
    int N_;
    double h_;
    VectorType<double> type_;
    RefCountPtr<GhostImporter<double> > importer_;
    mutable RefCountPtr<GhostView<double> > ghostView_;
    Array<int> localRows_;
    
  public:
    /** */
    ODE(const double& T, const double& x0, int n, const VectorType<double>& type)
      : NonlinearOperatorBase<double>(),
        T_(T),
        x0_(x0),
        N_(n*MPISession::getNProc()), 
        h_(0.0), 
        type_(type), 
        importer_(), 
        ghostView_(),
        localRows_()
    {
      int nLocal = n;
      h_ = T_/((double) N_-1);
      int nProc = MPISession::getNProc();
      int rank = MPISession::getRank();

      
      /* create a vector space in which n elements live on the local processor */
      int low = n*rank;
      int high = n*(rank+1);

      localRows_.resize(n);
      for (int i=0; i<n; i++)
        {
          localRows_[i] = low + i;
        }

      VectorSpace<double> space = type.createSpace(N_, n, &(localRows_[0]));
      setDomainAndRange(space, space);

      /* create an importer that will give us access to the ghost elements */

      if (nProc > 1)
        {
          int nGhosts = 2;
          if (rank==0 || rank==nProc-1) nGhosts = 1;
          
          Array<int> ghosts(nGhosts);
          if (rank==0)
            {
              ghosts[0] = n;
            }
          else if (rank==nProc-1)
            {
              ghosts[0] = localRows_[0]-1;
            }
          else
            {
              ghosts[0] = localRows_[0]-1;
              ghosts[1] = localRows_[n-1]+1;
            }
          
          importer_ = type.createGhostImporter(space, nGhosts, &(ghosts[0]));
        }
      else
        {
          importer_ = type.createGhostImporter(space, 0, 0);
        }
    }

    
    /** */
    LinearOperator<double> computeJacobianAndFunction(Vector<double>& f) const 
    {
      /* create a new operator */
      LinearOperator<double> J = type_.createMatrix(domain(), range());
      /* get a "view" of a loadable matrix underneath the operator */
      RefCountPtr<LoadableMatrix<double> > matview = J.matrix();

      /* get the current evaluation point */
      const Vector<double>& u = currentEvalPt();

      /* */
      f = range()->createMember();

      /* we'll need to look at ghost elements, so create a ghost view */
      importer_->importView(u, ghostView_);

      /* fill in the Jacobian and function value */
      int nProc = MPISession::getNProc();
      int rank = MPISession::getRank();
      for (int i=0; i<localRows_.size(); i++)
        {
          Array<int> colIndices;
          Array<double> colVals;
          int globalIndex = localRows_[i];
          const double& ui = ghostView_->getElement(globalIndex);

          if (rank==0 && i==0)
            {
              colIndices = tuple(globalIndex);
              colVals = tuple(1.0);
              f.setElement(globalIndex, ui-x0_);
            }
          else if (globalIndex == N_-1)
            {
              colIndices = tuple(globalIndex, globalIndex-1, globalIndex-2);
              colVals = tuple(cos(ui) + 1.5/h_, -2.0/h_, +0.5/h_);
              const double& u1 = ghostView_->getElement(globalIndex-1);
              const double& u2 = ghostView_->getElement(globalIndex-2);
              f.setElement(globalIndex, sin(ui) + (1.5*ui - 2.0*u1 + 0.5*u2)/h_);
            }
          else 
            {
              colIndices = tuple(globalIndex-1, globalIndex, globalIndex+1);
              colVals = tuple(-0.5/h_, cos(ui), 0.5/h_);
              const double& uPlus = ghostView_->getElement(globalIndex+1);
              const double& uMinus = ghostView_->getElement(globalIndex-1);
              f.setElement(globalIndex, sin(ui) + (uPlus-uMinus)/2.0/h_);              
            }
          matview->setRowValues(localRows_[i], colIndices.size(), 
                            &(colIndices[0]), &(colVals[0]));
        }

      matview->freezeValues();

      
      return J.ptr();
    }

    /** */
    Vector<double> getInitialGuess() const 
    {
      Vector<double> rtn = domain()->createMember();

      for (int i=0; i<localRows_.size(); i++)
        {
          double t = localRows_[i] / ((double) N_-1);

          rtn.setElement(localRows_[i], x0_*exp(-t));
        }
      return rtn;
    }

    /** */
    Vector<double> exactSoln() const 
    {
      Vector<double> rtn = domain()->createMember();

      for (int i=0; i<localRows_.size(); i++)
        {
          double t = localRows_[i] / ((double) N_-1);

          rtn.setElement(localRows_[i], 2.0*atan(exp(-t)*tan(x0_/2.0)));
        }
      return rtn;
    }

    /* */
    GET_RCP(NonlinearOperatorBase<double>);


  };
}



int main(int argc, void *argv[]) 
{
  try
    {
      int verbosity = 2;

      MPISession::init(&argc, &argv);

      VectorType<double> type = new EpetraVectorType();

      int n = 10;
      double T = 1.0;
      double A = 1.0;
      ODE* ode = new ODE(T, A, n, type);

      NonlinearOperator<double> F = ode;
      
      Vector<double> x0 = ode->getInitialGuess();

      cerr << "initial guess: " << x0 << endl;

      std::map<int,int> azOptions;
      std::map<int,double> azParams;

      azOptions[AZ_solver] = AZ_gmres;
      azOptions[AZ_precond] = AZ_dom_decomp;
      azOptions[AZ_subdomain_solve] = AZ_ilu;
      azOptions[AZ_graph_fill] = 1;
      azParams[AZ_max_iter] = 100;
      azParams[AZ_tol] = 1.0e-12;

      LinearSolver<double> linSolver = new AztecSolver(azOptions,azParams);

      cerr << "solver = " << linSolver << endl;
      
      NOX::TSF::Group grp(x0, F, linSolver);

      // Set up the status tests
      NOX::StatusTest::NormF statusTestA(grp, 1.0e-10);
      NOX::StatusTest::MaxIters statusTestB(20);
      NOX::StatusTest::Combo statusTestsCombo(NOX::StatusTest::Combo::OR, statusTestA, statusTestB);

      ParameterList solverParameters;
      solverParameters.set("Nonlinear Solver", "Line Search Based");
      
      ParameterList& lineSearchParameters 
        = solverParameters.sublist("Line Search");

      // Set the line search method
      lineSearchParameters.set("Method", "More'-Thuente");

      // Create the linear solver parameters sublist
      ParameterList& linearSolverParameters 
        = solverParameters.sublist("Linear Solver");

      // Set the line search method
      linearSolverParameters.set("Tolerance","1.0e-12");

      // Create the printing parameter sublist
      ParameterList& printingParameters = solverParameters.sublist("Printing");

      // Set the output precision
      printingParameters.set("Output Precision",8);

      NOX::Parameter::Teuchos2NOX converter;
      NOX::Parameter::List noxParameters = converter.toNOX(solverParameters);

      // Create the solver
      NOX::Solver::Manager solver(grp, statusTestsCombo, noxParameters);

      // Solve the nonlinear system
      NOX::StatusTest::StatusType status = solver.solve();

      // Print the answer
      cout << "\n" << "-- Parameter List From Solver --" << "\n";
      solver.getParameterList().print(cout);

      // Get the answer
      grp = solver.getSolutionGroup();

      // Print the answer
      cout << "\n" << "-- Final Solution From Solver --" << "\n";
      grp.print();

      cerr << F.currentEvalPt() << endl;

      cerr << ode->exactSoln() << endl;

      double error = (F.currentEvalPt() - ode->exactSoln()).norm2();
      cerr << "error norm = " << error / sqrt((float) n * MPISession::getNProc()) << endl;
      
    }
  catch(std::exception& e)
    {
      cerr << "Caught exception: " << e.what() << endl;
    }
  MPISession::finalize();

}

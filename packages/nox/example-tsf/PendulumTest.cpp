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
 * Solves the equation of motion of a pendulum as a boundary value problem. 
 *
 * The equation is u'' + 4.0*pi^2 sin(u) = 0.
 * Initial and final conditions are u(0)=0, u(T/4)=A.
 */


namespace TSFExtended
{
  class Pendulum : public NonlinearOperatorBase<double>
  {
  private:
    double T_;
    
  public:
    /** */
    Pendulum(const double& T, int n, const VectorType<double>& type)
      : NonlinearOperatorBase<double>(),
        T_(T), type_(type)
    {
      int nLocal = n;
      int nProc = MPISession::getNProc();
      int rank = MPISession::getRank();
      int dimension = nt*nProc;
      
      /* create a vector space in which n elements live on the local processor */
      int low = n*rank;
      int high = n*(rank+1);
      std::vector<int> localRows(n);
      for (int i=0; i<n; i++)
        {
          localRows[i] = low + i;
        }

      VectorSpace<double> space = type.createSpace(dimension, nt, localRows);
      setDomainAndRange(space, space);

      /* create an importer that will give us access to the ghost elements */

      int nGhosts = 2;
      if (rank==0 || rank==nProc-1) nGhosts = 1;

      std::vector<int> ghosts(nGhosts);
      if (rank==0)
        {
          ghosts[0] = n;
        }
      else if (rank==nProc-1)
        {
          ghosts[0] = localRows[0]-1;
        }
      else
        {
          ghosts[0] = localRows[0]-1;
          ghosts[1] = localRows[n-1]+1;
        }

      importer_ = type.createGhostImporter(space, nGhosts, &(ghosts[0]));
    }

    /** */
    void setT(double T) {T_ = T;}

    
    /** */
    LinearOperator<double> computeJacobianAndFunction(Vector<double>& f) const 
    {
      /* create a new operator */
      LinearOperator<double> J = type_.createMatrix(domain(), range());
      /* get a "view" of a loadable matrix underneath the operator */
      RefCountPtr<LoadableMatrix<double> > matview = J.matrix();


      /* get the current evaluation point */
      const Vector<double>& u = currentEvalPt();

      /* we'll need to look at ghost elements, so create a ghost view */
      importer_->importView(u, ghostView_);

      /* fill in the Jacobian and function value */
      int nProc = MPISession::getNProc();
      int rank = MPISession::getRank();
      for (int i=0; i<nLocalRows; i++)
        {
          Array<int> colIndices;
          Array<double> colVals;
          int globalIndex = localRows[i];
          const double& ui = ghostView_->getElement(globalIndex);

          if ((rank==0 && i==0) || (rank==nProc-1 && i==nLocalRows-1))
            {
              colIndices = tuple(globalIndex);
              colVals = tuple(1.0);
              if (rank==0) f[globalIndex] = ui-u0_;
              else f[globalIndex] = ui-u1_;
            }
          else
            {
              colIndices = tuple(globalIndex-1, globalIndex, globalIndex+1);
              colVals = tuple(-1.0/h_/h_, 2.0/h_/h_ + sin(ui), -1.0/h_/h_);
              const double& uPlus = ghostView_->getElement(globalIndex+1);
              const double& uMinus = ghostView_->getElement(globalIndex-1);
              f[globalIndex] = sin(ui) + (2.0*ui - uMinus - uPlus)/h_/h_;              
            }
          mat->setRowValues(localRows[i], colIndices.size(), 
                            &(colIndices[0]), &(colVals[0]));
        }

      matview->freezeValues();

      
      return J.ptr();
    }

    /** */
    Vector<double> getInitialGuess() const 
    {
      Vector<double> rtn = domain()->createMember();
      rtn.setElement(0, m_);
      return rtn;
    }

    /* */
    GET_RCP(NonlinearOperatorBase<double>);

  private:
    /* eccentricity of orbit */
    double e_;
    /* mean anomaly */
    double m_;
    /* vector type */
    VectorType<double> type_;
  };
}


int main(int argc, void *argv[]) 
{
  try
    {
      int verbosity = 2;

      MPISession::init(&argc, &argv);

      VectorType<double> type = new EpetraVectorType();

      /* a number of occasional interest */
      const double pi = 4.0*atan(1.0);

      /* eccentricity */
      double e = 0.1;
      Kepler* kepler = new Kepler(e, pi/4.0, type);
      // Kepler* kepler = new Kepler(e, pi/4.0, type);
      NonlinearOperator<double> F  = kepler;
      
      
      Vector<double> x0 = F.domain().createMember();
      x0.setElement(0, 0.1);

      ParameterList linSolverParams;

      linSolverParams.set(LinearSolverBase<double>::verbosityParam(), 2);
      linSolverParams.set(IterativeSolver<double>::maxitersParam(), 100);
      linSolverParams.set(IterativeSolver<double>::tolParam(), 1.0e-10);

      LinearSolver<double> linSolver 
        = new BICGSTABSolver<double>(linSolverParams);

      cerr << "solver = " << linSolver << endl;
      
      NOX::TSF::Group grp(x0, F, linSolver);

      // Set up the status tests
      NOX::StatusTest::NormF statusTestA(grp, 1.0e-10);
      NOX::StatusTest::MaxIters statusTestB(50);
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
      linearSolverParameters.set("Tolerance","1.0e-14");

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
    }
  catch(std::exception& e)
    {
      cerr << "Caught exception: " << e.what() << endl;
    }
  MPISession::finalize();

}

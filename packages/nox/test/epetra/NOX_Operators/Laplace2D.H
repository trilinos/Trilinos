//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
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

// NOX Test Lacplace 2D Problem
// ----------------------------
// Simple nonlinear PDE problem.
// This file provides the problem 
//
// -\Delta u + \lambda e^u = 0  in \Omega = (0,1) \times (0,1)
//                       u = 0  on \partial \Omega
//
// for use as a driver for various nox-epetra tests

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "NOX.H"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "NOX_Epetra_LinearSystem_AztecOO.H"

namespace Laplace2D {

// this is required to know the number of lower, upper, left and right
// node for each node of the Cartesian grid (composed by nx \timex ny 
// elements)

static void
get_myNeighbours( const int i, const int nx, const int ny,
			     int & left, int & right, 
			     int & lower, int & upper);


// This function creates a CrsMatrix, whose elements corresponds
// to the discretization of a Laplacian over a Cartesian grid,
// with nx grid point along the x-axis and and ny grid points 
// along the y-axis. For the sake of simplicity, I suppose that
// all the nodes in the matrix are internal nodes (Dirichlet
// boundary nodes are supposed to have been already condensated)

Epetra_CrsMatrix * 
CreateLaplacian( const int nx, const int ny, const Epetra_Comm * Comm);

// ==========================================================================
// This class contians the main definition of the nonlinear problem at
// hand. A method is provided to compute F(x) for a given x, and another
// method to update the entries of the Jacobian matrix, for a given x.
// As the Jacobian matrix J can be written as
//    J = L + diag(lambda*exp(x[i])),
// where L corresponds to the discretization of a Laplacian, and diag
// is a diagonal matrix with lambda*exp(x[i]). Basically, to update
// the jacobian we simply update the diagonal entries. Similarly, to compute
// F(x), we reset J to be equal to L, then we multiply it by the
// (distributed) vector x, then we add the diagonal contribution
// ==========================================================================

} // namespace Laplace2D


class PDEProblem 
{

public:

  // constructor. Requires the number of nodes along the x-axis
  // and y-axis, the value of lambda, and the Epetra_Communicator
  // (to define a Map, which is a linear map in this case)
  PDEProblem(const int nx, const int ny, const double lambda,
	     const Epetra_Comm * Comm);

  // destructor
  ~PDEProblem();

  // compute F(x)
  void ComputeF(const Epetra_Vector & x, Epetra_Vector & f);

  // update the Jacobian matrix for a given x
  void UpdateJacobian(const Epetra_Vector & x);

  // returns a pointer to the internally stored matrix
  Epetra_CrsMatrix * GetMatrix()
  {
    return Matrix_;
  }
  
private:
  
  int nx_, ny_;
  double hx_, hy_;
  Epetra_CrsMatrix * Matrix_;
  double lambda_;

}; /* class PDEProblem */
 
// ==========================================================================
// This is the main NOX class for this example. Here we define
// the interface between the nonlinear problem at hand, and NOX.
// The constructor accepts a PDEProblem object. Using a pointer
// to this object, we can update the Jacobian and compute F(x),
// using the definition of our problem. This interface is bit
// crude: For instance, no PrecMatrix nor Preconditioner is specified.
// ==========================================================================

class SimpleProblemInterface : public NOX::Epetra::Interface::Required      ,
                               public NOX::Epetra::Interface::Jacobian      ,
                               public NOX::Epetra::Interface::Preconditioner
{

public:
 
  //! Constructor
  SimpleProblemInterface( PDEProblem * Problem ) :
    Problem_(Problem) {};

  //! Destructor
  ~SimpleProblemInterface() 
  {
  };

  bool computeF(const Epetra_Vector & x, Epetra_Vector & f,
                NOX::Epetra::Interface::Required::FillType F )
  {
    Problem_->ComputeF(x,f);
    return true;
  };
  
  bool computeJacobian(const Epetra_Vector & x, Epetra_Operator & Jac)
  {

    Problem_->UpdateJacobian(x);
    return true;
  }

  bool computePreconditioner(const Epetra_Vector & x, Epetra_Operator & Op, Teuchos::ParameterList *)
  {

    Problem_->UpdateJacobian(x);
    return true;
  }

  bool computePrecMatrix(const Epetra_Vector & x, Epetra_RowMatrix & M) 
  {
    cout << "*ERR* SimpleProblem::preconditionVector()\n";
    cout << "*ERR* don't use explicit preconditioning" << std::endl;
    throw 1;
  }  
  
private:
  
  PDEProblem * Problem_;
  
}; /* class SimpleProblemInterface */


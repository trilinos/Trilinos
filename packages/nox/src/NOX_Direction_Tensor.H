#ifdef WITH_PRERELEASE
//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef NOX_DIRECTION_TENSOR_H
#define NOX_DIRECTION_TENSOR_H

#include "NOX_Direction_Generic.H" // base class
#include "NOX_Parameter_List.H"	   // class data element 

namespace NOX {
namespace Direction {

//! %Tensor direction computation
/*! 
  Computes the Tensor direction by finding the root or smallest magnitude 
  minimizer of the local model
  \f[ M_T(x_c+d) = F_c + Jd + a_c(s^Td)^2. \f]
  
  To use this direction, specify that the "Method" is "%Tensor" in the
  "%Direction" sublist of the parameters that are passed to the solver
  (see NOX::Direction::Manager for more information on choosing the
  search direction).  Be aware that the tensor step is not guaranteed
  to be a descent direction on the function.  Thus, not all "basic"
  line search approaches would be appropriate.  Instead, use the
  NOX::Solver::TensorBased solver, which approriately handles these
  considerations.

  This class also computes the curvilinear step, which is a parametric
  step that spans the directions of the tensor step and the %Newton
  step.  At lambda=1, the curvilinear step equals the full tensor
  step, and as lambda nears 0, the curvilinear step approaches the
  %Newton direction.  This step provides a monotonic decrease in the
  norm of the local tensor model as lambda varies from 0 to 1.

  \author Brett Bader (CU-Boulder) bader&nbsp;\@&nbsp;cs.colorado.edu

*/

class Tensor : public Generic {

public:

  //! Constructor
  Tensor(Parameter::List& params);

  //! Destructor
  virtual ~Tensor();

  // derived
  virtual bool reset(Parameter::List& params);

  // derived
  virtual bool compute(Abstract::Vector& dir, Abstract::Group& grp, 
			  const Solver::Generic& solver);

  //! Compute the curvilinear step dT(lambda) and return in dir.
  bool computeCurvilinearStep(Abstract::Vector& dir, Abstract::Group& grp, 
			      const Solver::Generic& solver, double lambda);


protected:

private:

  // Various functions that are taken from C to at least get a tensor
  // method implemented.  At some point these functions (as well as
  // the non-abstract vector/matrix data structures) should be
  // replaced with a class or abstraction layer.  I realize this is
  // not elegant, but it is quick and dirty and will work okay for the
  // moment.  The descriptions here are mostly placeholders.

  //!  Allocate storage for a matrix
  void** allocate_matrix(int rows, int cols, double**& a);

  //!  Allocate storage for a vector of ints
  void* allocate_vector(int n, int*& x);

  //!  Allocate storage for a vector of doubles
  void* allocate_vector(int n, double*& x);

  //!  Delete a matrix previosly allocated
  void delete_matrix(double** A);

  //!  Print a matrix
  void print_matrix(int rows, int cols, double** A);

  //!  Print a vector of doubles
  void print_vector(int n, double* x);

  //!  Print a vector of ints
  void print_vector(int n, int* x);

  //!  Calculate x'*y
  double inner_product(int n, double* x, double* y);

  //!  Calculate 2-norm of x
  double norm(int n, double* x);

  //!  Apply Givens rotations to vector v
  void givapp(double* c, double* s, double* v, int k);

  //!  Backsolve an upper triangular system Ux=b.
  double* backsolve(double** U, double* b, int* p, int n, int dim=0);
  
private:

  //! Types of preconditioning 
  enum preconditioningType {None, Left, Right};

  //! Types of reorthogonalization
  enum reorthogonalizeType {Never, AsNeeded, Always};

  //! Flag for type of preconditioning (left/right)
  preconditioningType precondition;

  //! Flag for reorthogonalizing basis after first pass of Gram-Schmidt
  reorthogonalizeType reorth;

  //! "Direction" sublist with parameters for the direction vector
  Parameter::List* paramsPtr;

  //! Vector containing the predicted RHS
  Abstract::Vector* predf;

  //! Temporary storage vectors
  Abstract::Vector* vecw;

  //! Tensor step
  Abstract::Vector* dTensor;

  //! Curvilinear step at the value lambdaBar
  Abstract::Vector* dTLambda;

  //! Newton step
  Abstract::Vector* dNewton;  // maybe get rid of this once dn gets into group

  //! Array of orthonormal basis vectors
  Abstract::Vector** basisVptr;

  //! Flag for "continuing" with reset() inside compute() where n is known
  bool isFreshlyReset;

  //! Flag for noting whether dt has been calculated
  bool isTensorCalculated;

  //! Flag for noting whether dt(lambdaBar) and params have been calculated
  bool isCLParamsCalculated;

  //! Problem size
  int probSize;

  //! Maximum number of Krylov iterations (kmax <= n)
  int kmax;

  //! Maximum number of starting vectors in initial basis
  int maxp;

  //! Frequency for printing iteration history
  int outputFreq;

  //! Maximum dimension of data structures (maxDim = kmax + p)
  int maxDim;

  //! Threshold tolerance for tensor solve
  double tol;

  //! Band Hessenberg matrix from Arnoldi process
  double** hess;

  //! Matrix of coefficients for Givens rotations
  double** givensC;

  //! Matrix of coefficients for Givens rotations
  double** givensS;

  //! Working vector for the new column of Hessenberg matrix
  double* newHessCol;

  //! Projected right hand side
  double* vecg;

  //! Projected right hand side
  double* vecq;

  //! Temporary working vector
  double* vecz;

  //! Vector for calculating the residual error
  double* terrvec; 

  //! Number of Arnoldi iterations for last solve
  int iterations;

  //! Counter for the number of total Arnoldi iterations since reset
  int arnoldiIters;

  //! Threshold for a real root in quadratic equation
  double lambdaBar;

  //! Parameter for the curvilinear linesearch
  double stJinvF;

  //! Parameter for the curvilinear linesearch
  double stJinvA;

  //! Norm of previous step direction, s = xp - xc
  double normS;

  //! Approximate value of s'*dt
  double beta;

  //! Smallest magnitude root/minimizer of quadratic equation
  double y1;

};
} // namespace Direction
} // namespace NOX
#endif
#endif

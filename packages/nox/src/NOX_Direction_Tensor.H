// $Id$ 
// $Source$ 

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

#ifndef NOX_DIRECTION_TENSOR_H
#define NOX_DIRECTION_TENSOR_H

#include "NOX_Common.H"

#ifdef WITH_PRERELEASE

#include "NOX_Direction_Generic.H" // base class
#include "NOX_Direction_Utils_InexactNewton.H" // class data element
#include "Teuchos_ParameterList.hpp"	   // class data element 
#include "NOX_Abstract_Group.H"    // returnType
#include "Teuchos_RCP.hpp"

// Forward declaration
namespace NOX {
class Utils;
}

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
  line search approaches would be appropriate.  Fortunately, the
  NOX::Solver::TensorBased solver approriately handles these
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
  Tensor(const Teuchos::RCP<NOX::GlobalData>& gd,
	 Teuchos::ParameterList& params);

  //! Destructor
  virtual ~Tensor();

  // derived
  virtual bool reset(const Teuchos::RCP<NOX::GlobalData>& gd,
		     Teuchos::ParameterList& params);

  // derived
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, 
		       const NOX::Solver::Generic& solver);

  // derived
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, 
		       const NOX::Solver::LineSearchBased& solver);

  //! Compute the curvilinear step dT(lambda) and return in dir.
  bool computeCurvilinearStep(NOX::Abstract::Vector& dir,
			      const NOX::Abstract::Group& soln, 
			      const NOX::Solver::Generic& solver,
			      double lambda) const;

  //! Compute the curvilinear step dT(lambda) and return in dir.
  bool computeCurvilinearStep2(NOX::Abstract::Vector& dir,
			       const NOX::Abstract::Group& soln, 
			       const NOX::Solver::Generic& solver,
			       double lambda) const;
  
  //! Return the approximate Newton step generated alongside tensor step 
  const NOX::Abstract::Vector& getNewton() const;

protected:

private:

  //! Types of steps 
  enum stepType {NewtonStep, TensorStep3, TensorStep2, TensorStepFP};

  //! Types of preconditioning 
  enum preconditioningType {None, Left, Right};

  //! Types of reorthogonalization
  enum reorthogonalizeType {Never, AsNeeded, Always};

  //! Print an error message
  void throwError(const std::string& functionName, const std::string& errorMsg) const;

  // Compute approximate Newton and tensor steps by solving the local models
  bool solveModels(NOX::Abstract::Vector& dir, NOX::Abstract::Group& soln, 
		   const NOX::Solver::Generic& solver, double& error,
		   double*& yn);

  // Compute the tensor step
  bool computeTensorStep(const NOX::Abstract::Group& soln, 
			 const NOX::Solver::Generic& solver,
			 double lambda) const;

  
private:

  // Various functions that are taken from C to at least get a tensor
  // method implemented.  At some point these functions (as well as
  // the non-abstract vector/matrix data structures) should be
  // replaced with a class or abstraction layer.  I realize this is
  // not elegant, but it is quick and dirty and will work okay for the
  // moment.  The descriptions here are mostly placeholders.

  //!  Allocate storage for a matrix
  void** allocate_matrix(int rows, int cols, double**& a) const;

  //!  Allocate storage for a vector of ints
  void* allocate_vector(int n, int*& x) const;

  //!  Allocate storage for a vector of doubles
  void* allocate_vector(int n, double*& x) const;

  //!  Delete a matrix previosly allocated
  void delete_matrix(double** A) const;

  //!  Print a matrix
  void print_matrix(int rows, int cols, double** A) const;

  //!  Print a vector of doubles
  void print_vector(int n, double* x) const;

  //!  Print a vector of ints
  void print_vector(int n, int* x) const;

  //!  Calculate x'*y
  double inner_product(int n, double* x, double* y) const;

  //!  Calculate 2-norm of x
  double norm(int n, double* x) const;

  //!  Compute Givens sin and cos values
  void computeGivens(double a, double b,
		     double& c, double& s) const;

  //!  Apply Givens rotations to vector [a; b]
  void applyGivens(double c, double s,
		   double& a, double& b) const;

  //!  Compute Householder vector
  void computeHouseholder(double** a, int j, 
			  int i1, int i2,
			  double* z) const;

  //!  Apply Householder rotations to a column vector in a matrix
  void applyHouseholder(const double* z, double** a,
			int j, int i1, int i2) const;

  //!  Apply Householder rotations to a vector
  void applyHouseholder(const double* z, double* a,
			int i1, int i2) const;

  //!  Find the smallest magnitude minimizer of quartic equation
  double minQuartic(double q1, double q2, double q3,
		    double q4, double q5,
		    bool chooseGlobal, double& qval) const;
  
  //!  Calculate beta by minimizing quartic equation
  double calculateBeta(double* vecg, double* vecq,
		       double h1, int kk0, int kk,
		       int pp, double normS, double& qval);

  //!  Calculate beta by solving the quadratic equation
  double calculateBeta(double qa, double qb, double qc,
		       double& qval, double& lambdaBar,
		       double normS = 1) const;

  //!  Calculate norm of local model residual
  double getNormModelResidual(const NOX::Abstract::Vector& dir,
			      const NOX::Abstract::Group& soln,
			      const NOX::Solver::Generic& solver,
			      bool isTensorModel) const;

  //!  Calculate directional derivative
  double getDirectionalDerivative(const NOX::Abstract::Vector& dir,
				  const NOX::Abstract::Group& soln) const;

  //!  Compute residuals and print out information
  void printDirectionInfo(const std::string& dirName,
			  const NOX::Abstract::Vector& dir,
			  const NOX::Abstract::Group& soln,
			  const NOX::Solver::Generic& solver,
			  bool isTensorModel=false) const;
  
  //!  Apply the preconditioner to a vector
  NOX::Abstract::Group::ReturnType 
  applyPreconditioner(bool useTranspose,
		      const NOX::Abstract::Group& soln,
		      Teuchos::ParameterList& params,
		      const NOX::Abstract::Vector& input,
		      NOX::Abstract::Vector& result,
		      const std::string& errLocation) const;
    
  
  //!  Backsolve an upper triangular system Ux=b.
  double* backsolve(double** U, double* b, int* p, int n, int dim=0,
		    bool isTranspose=false) const;
  

private:

  //! Global data pointer.  Keep this so the parameter list remains valid.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  //! Printing Utilities
  Teuchos::RCP<NOX::Utils> utils;

  //! Flag for the base direction to compute after the first iteration 
  stepType requestedBaseStep;

  //! Flag for the direction to be computed this iteration 
  stepType requestedStep;
  
  //! Flag for type of preconditioning (left/right)
  preconditioningType precondition;

  //! Flag for reorthogonalizing basis after first pass of Gram-Schmidt
  reorthogonalizeType reorth;

  //! "Direction" sublist with parameters for the direction vector
  /*!
    \note This is pointer rather than a reference to allow for the
    reset function.
  */
  Teuchos::ParameterList* paramsPtr;

  //! "Linear Solver" sublist inside "Direction"->"Tensor" sublist
  Teuchos::ParameterList* localParamsPtr;
  
  bool isAugmentSubspace;

  //! Determined based on "Use Shortcut Method".
  bool useShortcutMethod;

  //! Determined based on "Rescue Bad %Newton Solve".
  bool doRescue;

  //! Vector containing the tensor term sc
  Teuchos::RCP<NOX::Abstract::Vector> scPtr;

  //! Vector containing the tensor term ac
  Teuchos::RCP<NOX::Abstract::Vector> acPtr;

  //! Vector containing the product Jc*sc/norm(sc)
  Teuchos::RCP<NOX::Abstract::Vector> jcxshatPtr;

  //! Vector containing M'\sc
  Teuchos::RCP<NOX::Abstract::Vector> mtinvscPtr;

  //! Temporary storage vectors
  Teuchos::RCP<NOX::Abstract::Vector> vecw;

  //! Tensor step
  Teuchos::RCP<NOX::Abstract::Vector> dTensor;

  //! Initial guess of the step
  Teuchos::RCP<NOX::Abstract::Vector> dInitial;

  //! Curvilinear step at the value lambdaBar
  Teuchos::RCP<NOX::Abstract::Vector> dTLambda;

  //! Newton step
  Teuchos::RCP<NOX::Abstract::Vector> dNewton;  // maybe get rid of this once dn gets into group

  //! Array of orthonormal basis vectors
  std::vector< Teuchos::RCP<NOX::Abstract::Vector> > basisVptr;

  //! Flag for "continuing" with reset() inside compute() where n is known
  bool isFreshlyReset;

  //! Flag for noting if tensor step has been calculated (for curvilinear step)
  mutable bool isTensorCalculated;

  //! Flag for noting whether dt(lambdaBar) and params have been calculated
  mutable bool isCLParamsCalculated;

  //! Flag for noting whether to augment the subspace with direction sc
  bool isSubspaceAugmented;

  //! Flag for noting whether the transpose of the preconditioner is available
  bool isMinvTransAvailable;

  //! Flag for noting whether initial guess is nonzero
  bool isDir0;

  //! Value of inner product of initial guess and sc
  double dir0xsc;

  //! Problem size
  int probSize;

  //! Maximum number of times to restart local solver
  int maxRestarts;

  //! Maximum number of Krylov iterations (kmax <= n)
  int kmax;

  //! Maximum number of starting vectors in initial basis
  int pmax;

  //! Index holding either 0, k, or k+1
  int k1;
  int k;
  int p;

  //! Frequency for printing iteration history
  int outputFreq;

  //! Maximum dimension of data structures (maxDim = kmax + pmax)
  int maxDim;

  //! Threshold tolerance for tensor solve
  double tol;

  //! Absolute threshold tolerance 
  double errTol;

  //! Band Hessenberg matrix from Arnoldi process
  double** hess;

  //! Matrix of coefficients for Givens rotations
  double** givensC;

  //! Matrix of coefficients for Givens rotations
  double** givensS;

  //! Matrix of coefficients for Householder reflections
  double** hhZ;

  //! Temporary vector for Hessenberg matrix H*Q1
  double* hq1save;

  //! Temporary vector for right-hand side g
  double* gsave;

  //! Temporary vector for right-hand side q
  double* qsave;
  
  //! Vector for orthogonal Hessenberg matrix Q1, kth column
  double* qk;
  
  //! Vector for orthogonal Hessenberg matrix Q1, subdiagonal
  double* q1subdiag;
  
  //! Vector for orthogonal Hessenberg matrix Q1, norm-squared of column
  double* qNorm2;
  
  //! Vector for ghost column of H
  double* hk;
  
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

  //! Permutation index vector
  int* pindex; 

  //! Number of Arnoldi iterations for last solve
  int iterations;

  //! Counter for the number of total Arnoldi iterations since reset
  int arnoldiIters;

  //! Threshold for a real root in quadratic equation
  double lambdaBar;

  //! Parameter for the curvilinear linesearch
  double ata;

  //! Parameter for the curvilinear linesearch
  double atb;

  //! Parameter for the curvilinear linesearch
  double atc;

  //! Parameter for the curvilinear linesearch
  double btb;

  //! Parameter for the curvilinear linesearch
  double btc;

  //! Parameter for the curvilinear linesearch
  double ctc;

  //! Parameter for the curvilinear linesearch
  mutable double stJinvF;

  //! Parameter for the curvilinear linesearch
  mutable double stJinvA;

  //! Norm of previous step direction, s = xp - xc
  double normS;

  //! Approximate value of s'*dt
  mutable double beta;

  //! Smallest magnitude root/minimizer of quadratic equation
  mutable double y1;

  //! Counter for number of Jacobian-vector products
  mutable int multsJv;
  
  //! Counter for number of Preconditioner-vector products
  mutable int multsMv;

  //! Inexact Newton solve utilities
  NOX::Direction::Utils::InexactNewton inexactNewtonUtils;
  
};
} // namespace Direction
} // namespace NOX


#endif

#endif  // WITH_PRERELEASE

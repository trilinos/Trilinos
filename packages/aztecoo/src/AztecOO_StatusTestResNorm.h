
/* Copyright (2001) Sandia Corportation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 *
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef AZTECOO_STATUSTESTRESNORM_H
#define AZTECOO_STATUSTESTRESNORM_H

#include "AztecOO_StatusTest.h"
#include <vector>
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_Operator;

//! AztecOO_StatusTestResNorm: An implementation of AztecOO_StatusTest using a family of residual norms.

/*! AztecOO_StatusTestResNorm is an implementation of AztecOO_StatusTest that allows a user to construct
   one of a family of residual tests for use as a status/convergence test for AztecOO.  The form of the
   test is
   \f[
   \frac{\|r\|}{\sigma} \le \tau
   \f]
   where 
   <ul>
   <li> \f$r\f$ is the residual vector, implicitly or explicitly computed (determined by enum ResType),
   <li> \f$\|r\|\f$ is the residual norm determined by the enum NormType  (1-norm, 2-norm or inf-norm), 
   <li> \f$\sigma\f$ is the scale factor that can be passed in as a precomputed double precision number, 
   or can be selected from by the enum ScaleType (norm of RHS, norm of initial residual).
   <li> \f$\tau\f$ is the tolerance that is passed in as a double precision number to the constructor.  
   The value of \f$\tau\f$ can be reset using the ResetTolerance() method.
   </ul>

*/

class AztecOO_StatusTestResNorm: public AztecOO_StatusTest {

 public:

  //@{ \name Enums.
  /*! 
    \brief Select how the residual vector is produced.
  */
  enum ResType {Implicit, /*!< Use the residual vector produced by the iterative solver. */
		Explicit  /*!< Explicitly compute the residual vector r = b - A*x using the linear problem.
			    NOTE:  If if the ResType is set to Explicit, if the SolutionUpdated argument
			  to CheckStatus() is false, we will not be able to compute an explicit result and we
			  will proceed with the residual implicitly defined. */
  };
  /*! 
    \brief Select the norm type, where \f$ x\f$ is the norm argument and \f$ w \f$ is an optionally defined 
    weighting vector.
  */
  enum NormType {OneNorm, /*!< Use the one-norm \f$\sum_{i=1}^{n}(|x_i w_i|)\f$. */
		 TwoNorm, /*!< Use the two-norm \f$\sqrt(\sum_{i=1}^{n}((x_i w_i)^2)\f$. */
		 InfNorm  /*!< Use the inf-norm \f$(\max_{i=1}^{n}\{|x_i w_i|\})\f$. */
  };
  /*! 
    \brief Select the scale type.
  */
  enum ScaleType {NormOfRHS,     /*!< Use the norm of the right-hand-side. */
		  NormOfInitRes, /*!< Use the initial residual vector (exlicitly computed). */
		  None,          /*!< Use unscaled residual. */
		  UserProvided   /*!< User provides an explicit value that the norm of the residual will be divided by. */
  };
  //@}

  //@{ \name Constructors/destructors.
  //! Constructor
  /*! The constructor takes a single argument specifying the tolerance (\f$\tau\f$).  
    If none of the form definition methods are called, we use \f$\|r\|_2/\|r^{(0)}\|_2 \le \tau\f$ as the
    stopping criterion, where \f$\|r\|_2\f$ uses the least costly form of the 2-norm of 
    residual available from the iterative method and \f$\|r^{(0)}\|_2\f$ is the corresponding norm of the 
    initial residual.  The least costly form of the 2-norm depends on the chosen iterative method.  Most Krylov
    methods produce the preconditioned residual vector in a form that would be exact in infinite precision arithmetic.
    This vector may be different from the true residual either because left scaling or preconditioning was used, or
    because round-off error has introduced significant error, or both.

    \param Operator (In) The original linear operator that was passed in to the AztecOO solver object.
    \param LHS (In) The original left hand side vector that was passed in to the AztecOO solver object.  NOTE: AztecOO accepts
    multivector objects, but AztecOO_StatusTestResNorm does not handle residual tests for multiple vectors.  Most AztecOO users
    tend to have Epetra_Vectors, even though they are passing them in as Epetra_MultiVectors, so this should not be an issue.  If
    you are truly using Epetra_MultiVector objects, remember that for a multivector object mv, the ith vector can be extracted as
    mv(i).
    \param RHS (In) The original right hand side vector that was passed in to the AztecOO solver object.  See note for LHS.
    \param Tolerance (In) A double value that is used to test the convergence criterion.  Can be reset using ResetTolerance().
  */
  AztecOO_StatusTestResNorm(const Epetra_Operator & Operator, 
			    const Epetra_Vector & LHS, const Epetra_Vector & RHS,double Tolerance);

  //! Destructor
  virtual ~AztecOO_StatusTestResNorm();
  //@}

  //@{ \name Form and parameter definition methods.

  //! Define form of the residual, its norm and optional weighting vector.
  /*! This method defines the form of \f$\|r\|\f$.  We specify:
    <ul>
    <li> Whether the residual vector should be explicitly computed, or taken from the iterative method.
    <li> The norm to be used on the residual (this may be different than the norm used in DefineScaleForm()).
    <li> A weighting vector that will be multiplied, element by element, with the residual vector prior to 
    computing the norm.  This argument defaults to a zero vector and no weighting will be done.
    </ul>
  */
  int DefineResForm( ResType TypeOfResidual, NormType TypeOfNorm, Epetra_Vector * Weights = 0);
  
  //! Define form of the scaling, its norm, its optional weighting vector, or, alternatively, define an explicit value.
  /*! This method defines the form of how the residual is scaled (if at all).  It operates in two modes:
    <ol>
    <li> User-provided scaling value:
    <ul> 
    <li> Set argument TypeOfScaling to UserProvided.
    <li> Set ScaleValue to a non-zero value that the residual norm will be divided by.
    <li> TypeOfNorm and Weights arguments will be ignored.
    <li> Sample use:  Define ScaleValue = \f$\|A\|_{\infty}\f$ where \f$ A \f$ is the matrix of the linear problem.
    </ul>
    
    <li> Use a supported Scaling Form:
    <ul>
    <li> Define TypeOfScaling to be the norm of the right hand side, the initial residual vector, or to none.
    <li> Define norm to be used on the scaling vector (this may be different than the norm used in DefineResForm()).
    <li> Define a weighting vector that will be multiplied, element by element, with the residual vector prior to 
    computing the norm.  This argument defaults to a zero vector and no weighting will be done.
    </ul>
    </ol>
  */
  int DefineScaleForm( ScaleType TypeOfScaling, NormType TypeOfNorm, Epetra_Vector * Weights = 0, 
		       double ScaleValue = 1.0);
  //! Reset the value of the tolerance
  /*! We allow the tolerance to be reset for cases where, in the process of testing the residual, we find that
    the initial tolerance was too tight or too lax.
  */
  int ResetTolerance(double Tolerance) {tolerance_ = Tolerance; return(0);};
  //@}

  //@{ \name Methods that implement the AztecOO_StatusTest interface.
  //! Indicates if residual vector is required by this convergence test.
  /*! The value returned by this method will depend on several factors.  Once an AztecOO_StatusTestResNorm object
    is constructed and the DefineResForm and DefineScaleForm methods are optionally called, this method can tested.  
    For most Krylov solvers, there is no extra cost to providing the residual vector.  However, GMRES and Transpose-free
    QMR will need to explicitly compute this vector if ResidualVectorRequired() returns true, so this is an extra cost for
    these two iterative methods.
  */
  virtual bool ResidualVectorRequired() const;

  //! Check convergence status: Unconverged, Converged, Failed.
  /*! This method checks to see if the convergence criteria are met.  Depending on how the residual test
    is constructed this method will return the appropriate status type.

    \param CurrentIter (In) Current iteration of iterative method.

    \param CurrentResVector (In) The current residuals of the iterative process.  

    \param CurrentResNormEst (In) Estimate of the two-norm of the residual.  The value will be
    set to -1.0 if no estimate is available.

    \param SolutionUpdated (In) If this argument is true, then the solution vector that is part 
    of the Epetra_LinearProblem
    object being solved is consistent with the residual. 

    \return AztecOO_StatusType: Unconverged, Converged or Failed.
  */
  virtual AztecOO_StatusType CheckStatus(int CurrentIter, Epetra_MultiVector * CurrentResVector, 
				 double CurrentResNormEst,
				 bool SolutionUpdated);
  virtual AztecOO_StatusType GetStatus() const;

  virtual ostream& Print(ostream& stream, int indent = 0) const;
  //@}

 protected:
  //@{ \name Internal Methods.
  
  double ComputeNorm(const Epetra_Vector & vec, NormType typeofnorm);
  //@}
 private:

  //@{ \name Private data members.
  
  //! User's linear operator
  const Epetra_Operator & operator_;
  
  //! User's initial guess and solution vector
  const Epetra_Vector & lhs_;
   
  //! User's right hand side vector
  const Epetra_Vector & rhs_;
   
  //! Tolerance used to determine convergence
  double tolerance_;
 
  //! Type of residual to use (explicit or implicit)
  ResType restype_;
  
  //! Type of norm to use on residual (one, two or infinity)
  NormType resnormtype_;
  
  //! Type of scaling to use (Norm of RHS, Norm of Initial Residual, None or User provided)
  ScaleType scaletype_;
  
  //! Type of norm to use on scaling vector (one, two or infinity)
  NormType scalenormtype_;
  
  //! Vector of weights
  Epetra_Vector * weights_;
  
  //! Scale value.
  double scalevalue_;
  
  //! Residual norm value.
  double resvalue_;

  //! Test value = resvalue_ / scalevalue_
  double testvalue_;
  
  //! Status
  AztecOO_StatusType status_;
  
  //! Is residual vector required?
  bool resvecrequired_;

  //! Is this the first time CheckStatus is called?
  bool firstcallCheckStatus_;

  //! Is this the first time DefineResForm is called?
  bool firstcallDefineResForm_;

  //! Is this the first time DefineScaleForm is called?
  bool firstcallDefineScaleForm_;

  Epetra_Vector * localresvector_;

  //! Is this the current residual vector explicitly computed?
  bool curresvecexplicit_;


  //@}

};

#endif /* AZTECOO_STATUSTESTRESNORM_H */

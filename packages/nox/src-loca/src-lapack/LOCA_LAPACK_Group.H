// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

#ifndef LOCA_LAPACK_GROUP_H
#define LOCA_LAPACK_GROUP_H

#include "LOCA_Abstract_Group.H"	       // base class
#include "LOCA_Abstract_TransposeSolveGroup.H" // base class
#include "NOX_LAPACK_Group.H"                  // base class
#include "LOCA_Parameter_Vector.H"             // class data element
#include "LOCA_LAPACK_Interface.H"             // class data element

#include "Teuchos_ConfigDefs.hpp"              // for complex<double>

namespace LOCA {

  //! %LOCA BLAS/LAPACK support
  namespace LAPACK {

    //! Extension of the NOX::LAPACK::Group to %LOCA.  
    /*!
      This class is derived both from the NOX::LAPACK::Group and
      LOCA::Abstract::Group classes and therefore inherits the implementation
      of the NOX::Abstract::Group interface provided by NOX::LAPACK::Group.

      This class provides implementations of %LOCA AbstractGroup virtual 
      methods specific to the %LAPACK group.  It stores a parameter vector
      for setting/retrieving parameter values 
      (LOCA::Continuation::AbstractGroup), provides a facility for computing
      eigenvalues (LOCA::Continuation::AbstractGroup) using the %LAPACK
      routines DGEEV and DGGEV, augements the Jacobian matrix for homotopy
      (LOCA::Homotopy::AbstractGroup), and stores and manipulates a mass
      matrix (LOCA::Bifurcation::HopfBord::AbstractGroup).  Since it is
      derived from the LOCA::Abstract::Group (which is in-turn derived
      from all FiniteDifference groups), this group implicitly uses the
      finite-difference implementations of parameter and second derivatives
      provided by the FiniteDifference groups.  This group
      can therefore be used as an underlying group for all of LOCA's 
      continuation and bifurcation algorithms.

      The computeF() and computeJacobian() methods of the NOX::LAPACK::Group
      parent class are overloaded here.  They both set the entire contents
      of the parameter vector in the problem interface before calling the
      NOX::LAPACK::Group computeF() and computeJacobian().

      This group has several constructors supplying different information.
      All require a LOCA::LAPACK::Interface object to link to the 
      application code.  Set hasMassMat to true if the system has a mass matrix
      (only relevant for eigenvalue and Hopf calculations).  Finally,
      separate used and allocated row/column dimensions can be specified.  
      This functionality exists primarily to link with Fortran codes which
      preallocate all arrays to a fixed size but only use a portion of that
      array.
    */
    
    class Group : 
      public NOX::LAPACK::Group, 
      public LOCA::Abstract::Group,
      public virtual LOCA::Abstract::TransposeSolveGroup {

    public:

      //! Constructor 
      Group(const Teuchos::RCP<LOCA::GlobalData>& global_data,
	    LOCA::LAPACK::Interface& i);

      //! Copy constructor
      Group(const LOCA::LAPACK::Group& source, 
	    NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      ~Group();

      //! Assignment operator
      LOCA::LAPACK::Group& operator=(const LOCA::LAPACK::Group& source);

      /*! 
       * @name Overloaded NOX::LAPACK::Group  methods.
       */
      //@{

      //! Assignment operator
      NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);
      
      //! Assignment operator
      NOX::LAPACK::Group& operator=(const NOX::LAPACK::Group& source);

      //! Cloning function
      Teuchos::RCP<NOX::Abstract::Group> 
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Overloaded computeF()
      /*!
	Calls LOCA::LAPACK::Interface::setParams before evalulating F.
      */
      NOX::Abstract::Group::ReturnType computeF();

      //! Overloaded computeJacobian()
      /*!
	Calls LOCA::LAPACK::Interface::setParams before evalulating J.
      */
      NOX::Abstract::Group::ReturnType computeJacobian();

      //@}

      /*! 
       * @name Implementation of LOCA::Abstract::TransposeSolveGroup methods.
       */
      //@{

      //! Solve Jacobian-tranpose system 
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianTransposeInverse(Teuchos::ParameterList& params, 
				    const NOX::Abstract::Vector& input, 
				    NOX::Abstract::Vector& result) const;

      //! Solve Jacobian-tranpose system with multiple right-hand sides
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianTransposeInverseMultiVector(
				    Teuchos::ParameterList& params, 
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result) const;

      //@}
      
      /*! 
       * @name Implementation of LOCA::MultiContinuation::AbstractGroup virtual methods.
       */
      //@{

      //! Copy
      virtual void copy(const NOX::Abstract::Group& source);

      //! Set the parameter vector
      void setParams(const LOCA::ParameterVector& p);

      //! Set parameter indexed by paramID
      virtual void setParam(int paramID, double val);

      //! Set parameter indexed by paramID
      virtual void setParam(std::string paramID, double val);

      //! Return a const reference to the parameter vector owned by the group. 
      const LOCA::ParameterVector& getParams() const;

      //! Return copy of parameter indexed by paramID
      virtual double getParam(int paramID) const;

      //! Return copy of parameter indexed by paramID
      virtual double getParam(std::string paramID) const;

      //! Projects solution to a few scalars for multiparameter continuation
      /*! 
       * This method is called every time a solution is saved by the 
       * multiparameter continuation code MF for later visualization
       * and should project the solution vector down to a few scalars.
       * The array \c px will be preallocated to the proper length
       * given by projectToDrawDimension().
       *
       * The implementation here is to call the corresponding method
       * in the interface.
       */
      virtual void projectToDraw(const NOX::Abstract::Vector& x,
				 double *px) const;

      //! Returns the dimension of the project to draw array
      /*!
       * The implementation here is to call the corresponding method 
       * in the interface.
       */
      virtual int projectToDrawDimension() const;

      //! Compute a scaled dot product
      /*! 
       * The implementation here is a.dot(b) / a.length()
       */
      virtual double
      computeScaledDotProduct(const NOX::Abstract::Vector& a,
			      const NOX::Abstract::Vector& b) const;

      //! Print out the solution vector and continuation parameter
      void printSolution(const double conParam) const;

      //! Print out a vector and a parameter
      void printSolution(const NOX::Abstract::Vector& x_, 
			 const double conParam) const;

      //! Scales a vector
      /*! 
       * The implementation here is x.scale(1.0/sqrt(x.length))
       */
      virtual void
      scaleVector(NOX::Abstract::Vector& x) const;

      //@}

      /*! 
       * @name Implementation of LOCA::TimeDependent::AbstractGroup virtual methods.
       */
      //@{

      //! Compute the shifted matrix
      virtual NOX::Abstract::Group::ReturnType
      computeShiftedMatrix(double alpha, double beta);

      //! Multiply the shifted matrix by a vector.
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrix(const NOX::Abstract::Vector& input,
			 NOX::Abstract::Vector& result) const;

      //! Multiply the shifted matrix by a multi-vector.
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrixMultiVector(
				const NOX::Abstract::MultiVector& input,
				NOX::Abstract::MultiVector& result) const;
      
      /*! 
       * \brief Apply the inverse of the shifted matrix by a multi-vector, as 
       * needed by the shift-and-invert and generalized Cayley transformations.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrixInverseMultiVector(
			        Teuchos::ParameterList& params, 
				const NOX::Abstract::MultiVector& input,
				NOX::Abstract::MultiVector& result) const;

      //@}

      /*! 
       * @name Implementation of LOCA::Hopf::MooreSpence::AbstractGroup virtual methods.
       */
      //@{

      //! Is  \f$J+i\omega B\f$ valid
      virtual bool isComplex() const;

      //! Compute \f$J+i\omega B\f$
      /*!
       * The argument \b frequency stores \f$\omega\f$.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeComplex(double frequency);

      //! Compute \f$(J+i\omega B)(y+iz)\f$
      virtual NOX::Abstract::Group::ReturnType
      applyComplex(const NOX::Abstract::Vector& input_real,
		   const NOX::Abstract::Vector& input_imag,
		   NOX::Abstract::Vector& result_real,
		   NOX::Abstract::Vector& result_imag) const;

      //! Compute \f$(J+i\omega B)(y+iz)\f$
      virtual NOX::Abstract::Group::ReturnType
      applyComplexMultiVector(const NOX::Abstract::MultiVector& input_real,
			      const NOX::Abstract::MultiVector& input_imag,
			      NOX::Abstract::MultiVector& result_real,
			      NOX::Abstract::MultiVector& result_imag) const;

      //! Solve \f$(J+i\omega B)(y+iz) = a+ib\f$
      virtual NOX::Abstract::Group::ReturnType
      applyComplexInverseMultiVector(
			    Teuchos::ParameterList& params,
			    const NOX::Abstract::MultiVector& input_real,
			    const NOX::Abstract::MultiVector& input_imag,
			    NOX::Abstract::MultiVector& result_real,
			    NOX::Abstract::MultiVector& result_imag) const;

      //@}

      /*! 
       * @name Implementation of LOCA::Hopf::MinimallyAugmented::AbstractGroup virtual methods.
       */
      //@{

      /*!
       * Computes conjugate-tranpose matrix vector product 
       * \f$ (J+i\omega B)^H (x + iy) \f$.
       */
      virtual NOX::Abstract::Group::ReturnType 
      applyComplexTranspose(const NOX::Abstract::Vector& input_real,
			    const NOX::Abstract::Vector& input_imag,
			    NOX::Abstract::Vector& result_real,
			    NOX::Abstract::Vector& result_imag) const;

      /*!
       * Computes conjugate-tranpose matrix vector product 
       * \f$ (J+i\omega B)^H (x + iy) \f$.
       */
      virtual NOX::Abstract::Group::ReturnType 
      applyComplexTransposeMultiVector(
			   const NOX::Abstract::MultiVector& input_real,
			   const NOX::Abstract::MultiVector& input_imag,
			   NOX::Abstract::MultiVector& result_real,
			   NOX::Abstract::MultiVector& result_imag) const;

      //! Solve \f$(J+i\omega B)^H (x + iy) = a+ib\f$
      virtual NOX::Abstract::Group::ReturnType
      applyComplexTransposeInverseMultiVector(
			    Teuchos::ParameterList& params,
			    const NOX::Abstract::MultiVector& input_real,
			    const NOX::Abstract::MultiVector& input_imag,
			    NOX::Abstract::MultiVector& result_real,
			    NOX::Abstract::MultiVector& result_imag) const;

      //@}

      /*! 
       * @name Implementation of LOCA::Homotopy::AbstractGroup virtual methods.
       */
      //@{
      /*! 
       * \brief Replace Jacobian \f$J\f$ by \f$aJ+bI\f$ where \f$I\f$ is 
       * the identity matrix.
       */
      virtual NOX::Abstract::Group::ReturnType 
      augmentJacobianForHomotopy(double a, double b);

      //@}

      //! Return reference to Jacobian matrix
      NOX::LAPACK::Matrix<double>& getJacobianMatrix() {
	return jacSolver.getMatrix();
      }

      //! Return reference to Jacobian matrix
      const NOX::LAPACK::Matrix<double>& getJacobianMatrix() const {
	return jacSolver.getMatrix();
      }

      //! Return reference to shifted matrix
      NOX::LAPACK::Matrix<double>& getShiftedMatrix() {
	return shiftedSolver.getMatrix();
      }

      //! Return reference to shifted matrix
      const NOX::LAPACK::Matrix<double>& getShiftedMatrix() const {
	return shiftedSolver.getMatrix();
      }

#ifdef HAVE_TEUCHOS_COMPLEX

      //! Return reference to complex matrix
      NOX::LAPACK::Matrix< std::complex<double> >& getComplexMatrix() {
	return complexSolver.getMatrix();
      }

      //! Return reference to complex matrix
      const NOX::LAPACK::Matrix< std::complex<double> >& 
      getComplexMatrix() const {
	return complexSolver.getMatrix();
      }

#endif

    protected:

      //! resets isValid flags
      void resetIsValid();

    protected:

      //! LOCA Global data object
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Referece to current problem
      LOCA::LAPACK::Interface& locaProblemInterface;

      //! vector of parameters
      ParameterVector params;

      //! Shifted matrix (alpha*J+beta*M)
      mutable NOX::LAPACK::LinearSolver<double> shiftedSolver;

      //! Frequency for Hopf calculations
      double freq;

      //! Flag indicating whether complex matrix is valid
      bool isValidComplex;

#ifdef HAVE_TEUCHOS_COMPLEX
      //! Complex matrix
      mutable NOX::LAPACK::LinearSolver< std::complex<double> > complexSolver;
#endif
      
    };

  } // namespace LAPACK
} // namespace LOCA


#endif

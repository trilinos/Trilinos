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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_LAPACK_GROUP_H
#define LOCA_LAPACK_GROUP_H

#include "LOCA_Abstract_Group.H"	// base class
#include "NOX_LAPACK_Group.H"           // base class
#include "LOCA_Parameter_Vector.H"      // class data element
#include "LOCA_LAPACK_Interface.H"      // class data element

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
      (only relevant for eigenvalue and Hopf calculations).  Pass a 
      parameter list to change the algorithm used for singular Jacobian
      matrix solves (see LOCA::Bifurcation::TPBord::SingularSolveGroup,
      only relevant for turning point and pitchfork calculations).  Finally,
      separate used and allocated row/column dimensions can be specified.  
      This functionality exists primarily to link with Fortran codes which
      preallocate all arrays to a fixed size but only use a portion of that
      array.
    */
    
    class Group : public NOX::LAPACK::Group, public LOCA::Abstract::Group {

    public:

      //! Constructor 
      /*!
       * Set hasMassMat to true if the system has a mass matrix
       * (Only relevant for eigenvalue and Hopf calculations).
       */
      Group(LOCA::LAPACK::Interface& i, bool hasMassMat = false);

      //! Constructor with parameter list 
      /*!
       * Pass a parameter for singular matrix solves. (see 
       * LOCA::Bifurcation::TPBord::SingularSolveGroup)
       */
      Group(NOX::Parameter::List& params,
	    LOCA::LAPACK::Interface& i, bool hasMassMat = false);

      //! Constructor with used and allocated dimensions of matrices
      Group(LOCA::LAPACK::Interface& i, int m, int n, bool hasMassMat = false);

      //! Constructor with parameter list and with used, allocated dims
      Group(NOX::Parameter::List& params,
	    LOCA::LAPACK::Interface& i, int m, int n, bool hasMassMat = false);

      //! Copy constructor
      Group(const Group& source, 
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
      NOX::Abstract::Group* clone(NOX::CopyType type = NOX::DeepCopy) const;

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
       * @name Implementation of LOCA::Abstract::Group virtual methods.
       */
      //@{

      //! Assignment operator
      LOCA::Abstract::Group& operator=(const LOCA::Abstract::Group& source);

      //@}

      /*! 
       * @name Implementation of LOCA::Continuation::AbstractGroup virtual methods.
       */
      //@{

      //! Set the parameter vector
      void setParams(const LOCA::ParameterVector& p);

      //! Set parameter indexed by paramID
      virtual void setParam(int paramID, double val);

      //! Set parameter indexed by paramID
      virtual void setParam(string paramID, double val);

      //! Return a const reference to the parameter vector owned by the group. 
      const LOCA::ParameterVector& getParams() const;

      //! Return copy of parameter indexed by paramID
      virtual double getParam(int paramID) const;

      //! Return copy of parameter indexed by paramID
      virtual double getParam(string paramID) const;

      //! Compute a scaled dot product
      /*! 
       * The implementation here is a.dot(b) / a.length()
       */
      virtual double
      computeScaledDotProduct(const NOX::Abstract::Vector& a,
			      const NOX::Abstract::Vector& b) const;

      //! Scales a vector
      /*! 
       * The implementation here is x.scale(1.0/sqrt(x.length))
       */
      virtual void
      scaleVector(NOX::Abstract::Vector& x) const;

      //! applyJacobianInverse for multiple right-hand sides
      virtual NOX::Abstract::Group::ReturnType
      applyJacobianInverseMulti(NOX::Parameter::List& params,
			   const NOX::Abstract::Vector* const* inputs,
			   NOX::Abstract::Vector** outputs, int nVecs) const;

      //! Compute eigenvalues/eigenvectors
      virtual NOX::Abstract::Group::ReturnType
      computeEigenvalues(NOX::Parameter::List& params);

      //! Print out the solution vector and continuation parameter
      void printSolution(const double conParam) const;

      //! Print out a vector and a parameter
      void printSolution(const NOX::Abstract::Vector& x_, 
			 const double conParam) const;

      //@}

      
      /*! 
       * @name Implementation of LOCA::MultiContinuation::AbstractGroup virtual methods.
       */
      //@{

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

      //@}

      /*! 
       * @name Implementation of LOCA::Bifurcation::TPBord::AbstractGroup virtual methods.
       */
      //@{

      /*!
       * \brief Solve a bordered system of equations
       */
      /*!
       * Compute the solution to the bordered system of equations
       * \f[
       *     \begin{bmatrix}
       *        J & a \\
       *        b^T & 0
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *        \text{vResult} \\
       *        \text{sResult}
       *     \end{bmatrix} = 
       *     \begin{bmatrix}
       *        \text{vInput} \\
       *        \text{sInput}
       *     \end{bmatrix}
       * \f]
       */
      virtual NOX::Abstract::Group::ReturnType
      applyBorderedJacobianInverse(bool trans,
				   NOX::Parameter::List& params,
				   const NOX::Abstract::Vector& a,
				   const NOX::Abstract::Vector& b,
				   const NOX::Abstract::Vector& vInput,
				   double sInput,
				   NOX::Abstract::Vector& vResult,
				   double& sResult) const;

      //@}

      /*! 
       * @name Implementation of LOCA::TimeDependent::AbstractGroup virtual methods.
       */
      //@{

      //! Compute the mass matrix of the system
      virtual NOX::Abstract::Group::ReturnType
      computeMassMatrix();

      //! Multiply the mass matrix by a vector
      virtual NOX::Abstract::Group::ReturnType
      applyMassMatrix(const NOX::Abstract::Vector& input,
		      NOX::Abstract::Vector& result) const;

      //! Returns true if mass matrix is valid
      virtual bool isMassMatrix() const;

      //@}

      /*! 
       * @name Implementation of LOCA::Bifurcation::HopfBord::AbstractGroup virtual methods.
       */
      //@{

      //! Solve \f$(J+i\omega B)(y+iz) = a+ib\f$
      virtual NOX::Abstract::Group::ReturnType
      applyComplexInverse(NOX::Parameter::List& params,
			  const NOX::Abstract::Vector& input_real,
			  const NOX::Abstract::Vector& input_imag,
			  double frequency,
			  NOX::Abstract::Vector& result_real,
			  NOX::Abstract::Vector& result_imag) const;

      //! applyComplexInverse for multiple right-hand-sides
      virtual NOX::Abstract::Group::ReturnType
      applyComplexInverseMulti(NOX::Parameter::List& params,
			       const NOX::Abstract::Vector* const* inputs_real,
			       const NOX::Abstract::Vector* const* inputs_imag,
			       double frequency,
			       NOX::Abstract::Vector** results_real,
			       NOX::Abstract::Vector** results_imag,
			       int nVecs) const;

      //@}

      /*! 
       * @name Implementation of LOCA::Homotopy::AbstractGroup virtual methods.
       */
      //@{

      /*! 
       * \brief Replace Jacobian \f$J\f$ by \f$J+pI\f$ where \f$I\f$ is the 
       * identity matrix and \f$p\f$ is a scalar (\b conParamValue).
       */
      virtual NOX::Abstract::Group::ReturnType 
      augmentJacobianForHomotopy(double conParamValue);

      //@}

      //! Print out a vector and a parameter
      void printSolution(const NOX::LAPACK::Vector& x_, 
			 const double conParam) const;

      //! %LAPACK interface to applyMassMatrix
      virtual NOX::Abstract::Group::ReturnType
      applyMassMatrix(const NOX::LAPACK::Vector& input,
		      NOX::LAPACK::Vector& result) const;
      
      //! Returns true if system has a mass matrix
      bool hasMass() const;

    protected:

      //! resets isValid flags
      void resetIsValid();

    protected:
      //! Referece to current problem
      LOCA::LAPACK::Interface& locaProblemInterface;

      //! vector of parameters
      ParameterVector params;

      //! Mass matrix
      NOX::LAPACK::Matrix massMatrix;

      //! Flag indicating whether we have a mass matrix
      bool hasMassMatrix;

      //! Flag indicating whether mass matrix is valid
      bool isValidMass;
      
    };

  } // namespace LAPACK
} // namespace LOCA


#endif

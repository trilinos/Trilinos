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

#ifndef LOCA_ABSTRACT_GROUP_H
#define LOCA_ABSTRACT_GROUP_H

#include "LOCA_Continuation_AnasaziGroup.H"                   // base class
#include "LOCA_Bifurcation_HopfBord_FiniteDifferenceGroup.H"  // base class
#include "LOCA_Bifurcation_TPBord_SingularSolveGroup.H"       // base class
#include "LOCA_Homotopy_AbstractGroup.H"                      // base class
#include "LOCA_MultiContinuation_FiniteDifferenceGroup.H"     // base class

namespace LOCA { 

  //! %LOCA abstract interface namespace
  namespace Abstract { 

    /*!
     * \brief Compatiblity class for AbstractGroup hierarchy.
     */
    /*!
     * This class is derived from all %LOCA AbstractGroup abstract base
     * classes as well as all FiniteDifference groups and any other groups
     * that provided default implementations for AbstractGroup pure virtual
     * methods.  This class provides definitions for all needed assignment
     * operators and provides definitions for some pure virtual methods by
     * printing error messages.  This class exists primarily for compatiblity
     * to an older class hierarchy and will most likely be removed in the
     * future.
     */
    class Group : 
      public LOCA::Continuation::AnasaziGroup,
      public LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup,
      public LOCA::Bifurcation::TPBord::SingularSolveGroup,
      public virtual LOCA::Homotopy::AbstractGroup,
      public LOCA::MultiContinuation::FiniteDifferenceGroup
    {

    public:

      //! Default constructor.
      Group(const DerivUtils& deriv = LOCA::DerivUtils());
  
      /*! 
       * \brief Constructor with parameter list 
       * (for LOCA::Bifurcation::TPBord::SingularSolveGroup).
       */
      Group(NOX::Parameter::List& params,
	    const DerivUtils& deriv = LOCA::DerivUtils());
 
      //! Copy constructor
      Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~Group();

      //! Assignment operator
      virtual Group& operator=(const Group& source);

      /*! 
       * @name Implementation of LOCA::Continuation::AbstractGroup virtual methods.
       */
      //@{

      //! Assignment operator
      virtual LOCA::Continuation::AbstractGroup& 
      operator=(const LOCA::Continuation::AbstractGroup& source);

      //@}

      /*! 
       * @name Implementation of LOCA::Continuation::AnasaziGroup virtual methods.
       */
      //@{

      //! Assignment operator
      virtual LOCA::Continuation::AnasaziGroup& 
      operator=(const LOCA::Continuation::AnasaziGroup& source);

      //@}

      
      /*! 
       * @name Implementation of LOCA::Continuation::FiniteDifferenceGroup virtual methods.
       */
      //@{

      //! Assignment operator
      virtual LOCA::Continuation::FiniteDifferenceGroup& 
      operator=(const LOCA::Continuation::FiniteDifferenceGroup& source);

      //@}

      /*! 
       * @name Implementation of LOCA::Bifurcation::TPBord::AbstractGroup virtual methods.
       */
      //@{

      //! Assignment operator
      virtual LOCA::Bifurcation::TPBord::AbstractGroup& 
      operator=(const LOCA::Bifurcation::TPBord::AbstractGroup& source);

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
       * The default implementation here is to throw an error.
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
       * @name Implementation of LOCA::Bifurcation::TPBord::FiniteDiffernceGroup virtual methods.
       */
      //@{

      //! Assignment operator
      virtual LOCA::Bifurcation::TPBord::FiniteDifferenceGroup& 
      operator=(const LOCA::Bifurcation::TPBord::FiniteDifferenceGroup& source);
      //@}

      /*! 
       * @name Implementation of LOCA::Bifurcation::TPBord::SingularSolveGroup virtual methods.
       */
      //@{

      //! Assignment operator
      virtual LOCA::Bifurcation::TPBord::SingularSolveGroup& 
      operator=(const LOCA::Bifurcation::TPBord::SingularSolveGroup& source);

      //@}

      /*! 
       * @name Implementation of LOCA::TimeDependent::AbstractGroup virtual methods.
       */
      //@{

      //! Assignment operator
      virtual LOCA::TimeDependent::AbstractGroup& 
      operator=(const LOCA::TimeDependent::AbstractGroup& source);

      //! Compute the mass matrix of the system
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeMassMatrix();

      //! Multiply the mass matrix by a vector
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyMassMatrix(const NOX::Abstract::Vector& input,
		      NOX::Abstract::Vector& result) const;

      //! Is mass matrix valid
      /*!
       * Implementation here returns false.
       */
      virtual bool isMassMatrix() const;

      //@}

      /*! 
       * @name Implementation of LOCA::Bifurcation::HopfBord::AbstractGroup virtual methods.
       */
      //@{
      
      //! Assignment operator
      virtual LOCA::Bifurcation::HopfBord::AbstractGroup& 
      operator=(const LOCA::Bifurcation::HopfBord::AbstractGroup& source);

      //! Solve \f$(J+i\omega B)(y+iz) = a+ib\f$
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyComplexInverse(NOX::Parameter::List& params,
			  const NOX::Abstract::Vector& input_real,
			  const NOX::Abstract::Vector& input_imag,
			  double frequency,
			  NOX::Abstract::Vector& result_real,
			  NOX::Abstract::Vector& result_imag) const;

      //@}

      /*! 
       * @name Implementation of LOCA::Bifurcation::HopfBord::FiniteDiffernceGroup virtual methods.
       */
      //@{

      //! Assignment operator
      virtual LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup& 
      operator=(const LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup& source);

      //@}

      /*! 
       * @name Implementation of LOCA::Homotopy::AbstractGroup virtual methods.
       */
      //@{

      //! Assignment operator
      virtual LOCA::Homotopy::AbstractGroup& 
      operator=(const LOCA::Homotopy::AbstractGroup& source);

      /*! 
       * \brief Replace Jacobian \f$J\f$ by \f$J+pI\f$ where \f$I\f$ is the 
       * identity matrix and \f$p\f$ is a scalar (\b conParamValue).
       */
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType 
      augmentJacobianForHomotopy(double conParamValue);

      //@}

      /*! 
       * @name Implementation of LOCA::MultiContinuation::AbstractGroup virtual methods.
       */
      //@{

      //! Assignment operator
      virtual LOCA::MultiContinuation::AbstractGroup& 
      operator=(const LOCA::MultiContinuation::AbstractGroup& source);

      //! Set parameters indexed by (integer) paramIDs
      virtual void 
      setParams(const vector<int>& paramIDs, 
		const NOX::Abstract::MultiVector::DenseMatrix& vals);

      //@}

      /*! 
       * @name Implementation of LOCA::MultiContinuation::FiniteDifferenceGroup virtual methods.
       */
      //@{

      //! Assignment operator
      virtual LOCA::MultiContinuation::FiniteDifferenceGroup& 
      operator=(const LOCA::MultiContinuation::FiniteDifferenceGroup& source);

      //@}
      
    };
  } // namespace Abstract
} // namespace LOCA

#endif

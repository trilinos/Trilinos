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

#ifndef LOCA_MULTICONTINUATION_NATURALGROUP_H
#define LOCA_MULTICONTINUATION_NATURALGROUP_H

#include "LOCA_MultiContinuation_ExtendedGroup.H"       // base class

namespace LOCA { 

  namespace MultiContinuation { 

    
    class NaturalGroup : 
      public virtual LOCA::MultiContinuation::ExtendedGroup {

    public:
  
      //! Constructor with integer continuation parameter ids
      /*!
       * The argument \em grp is the underlying group which implements \f$F\f$,
       * \em paramID is the index of the continuation parameter, 
       * \em params
       * is the parameter list for the continuation groups described above.
       *
       * The \em grp argument is a non-const reference and may be modified.
       */
      NaturalGroup(LOCA::MultiContinuation::AbstractGroup& grp, 
		   const vector<int>& paramIDs,
		   NOX::Parameter::List& params);

      //! Constructor for single parameter continuation
      /*!
       * The argument \em grp is the underlying group which implements \f$F\f$,
       * \em paramID is the name of the continuation parameter, 
       * \em params
       * is the parameter list for the continuation groups described above.
       *
       * The \em grp argument is a non-const reference and may be modified.
       */
      NaturalGroup(LOCA::MultiContinuation::AbstractGroup& grp, 
		   const string& paramID,
		   NOX::Parameter::List& params);

      //! Copy constructor
      NaturalGroup(const NaturalGroup& source, 
		   NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~NaturalGroup();

      //! Assignment operator
      virtual NaturalGroup& 
      operator=(const NaturalGroup& source);

      /*! 
       * @name Implementation of LOCA::MultiContinuation::ExtendedGroup virtual methods
       */
      //@{

      //! Assignment operator
      virtual LOCA::MultiContinuation::ExtendedGroup& 
      operator=(const LOCA::MultiContinuation::ExtendedGroup& source);

      //! Compute continuation constraint equations
      virtual NOX::Abstract::Group::ReturnType
      computeConstraints();

      //! Compute continuation constraint equations
      virtual NOX::Abstract::Group::ReturnType
      computeConstraintDerivatives();

      //! Return \c true if constraint residuals are valid
      virtual bool isConstraints() const;

      //! Return \c true if constraint derivatives are valid
      virtual bool isConstraintDerivatives() const;

      //! Return constraint residuals
      virtual const NOX::Abstract::MultiVector::DenseMatrix&
      getConstraints() const;

      //! Return solution component of constraint derivatives
      virtual const NOX::Abstract::MultiVector*
      getConstraintDerivativesX() const;

      //! Return parameter component of constraint derivatives
      virtual const NOX::Abstract::MultiVector::DenseMatrix*
      getConstraintDerivativesP() const;

      /*! 
       * \brief Return \c true if solution component of constraint 
       * derivatives is zero
       */
      virtual bool isConstraintDerivativesXZero() const;

      /*! 
       * \brief Return \c true if solution component of constraint 
       * derivatives is zero
       */
      virtual bool isConstraintDerivativesPZero() const;

      //@}

      /*! 
       * @name Implementation of LOCA::Extended::AbstractGroup virtual methods 
       */
      //@{

      //! Assignment operator
      virtual LOCA::Extended::AbstractGroup& 
      operator=(const LOCA::Extended::AbstractGroup& source);

      //@}

      /*! 
       * @name Implementation of NOX::Abstract::Group virtual methods 
       */
      //@{

      //! Assignment operator
      virtual NOX::Abstract::Group& 
      operator=(const NOX::Abstract::Group& source);

      //! Clone function
      virtual NOX::Abstract::Group* 
      clone(NOX::CopyType type = NOX::DeepCopy) const;
      
      //@}

      /*! 
       * @name Implementation of LOCA::Stepper::AbstractGroup virtual methods 
       */
      //@{

      //! Set the previous solution vector y
      virtual void setPrevX(const NOX::Abstract::Vector& y);

      //! Set step size for continuation constraint equation \em i
      virtual void setStepSize(double deltaS, int i = 0);

      //@}

    protected:

      //! Resets all isValid flags to false
      virtual void resetIsValid();

    protected:

      //! Stores constraint residuals
      NOX::Abstract::MultiVector::DenseMatrix constraints;

      //! Stores derivative of constraints w.r.t. parameters
      NOX::Abstract::MultiVector::DenseMatrix dgdp;

      //! Flag indicating whether constraints are valid
      bool isValidConstraints;
    };

  } // namespace MultiContinuation
} // namespace LOCA

#endif

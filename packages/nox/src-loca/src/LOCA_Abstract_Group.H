// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

#ifndef LOCA_ABSTRACT_GROUP_H
#define LOCA_ABSTRACT_GROUP_H

#include "Teuchos_RefCountPtr.hpp"

#include "LOCA_Homotopy_AbstractGroup.H"
#include "LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.H"
#include "LOCA_Pitchfork_MooreSpence_AbstractGroup.H"
#include "LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.H"

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
      public virtual LOCA::Homotopy::AbstractGroup,
      public virtual LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup,
      public virtual LOCA::Pitchfork::MooreSpence::AbstractGroup,
      public virtual LOCA::Hopf::MooreSpence::FiniteDifferenceGroup
    {

    public:

      //! Constructor
      Group(const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data);

      //! Constructor
      Group(const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	    const Teuchos::RefCountPtr<LOCA::DerivUtils>& deriv);
 
      //! Copy constructor
      Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~Group();

      /*! 
       * @name Implementation of LOCA::Homotopy::AbstractGroup virtual methods.
       */
      //@{

      /*! 
       * \brief Replace Jacobian \f$J\f$ by \f$pJ+(1-p)I\f$ where \f$I\f$ is 
       * the identity matrix and \f$p\f$ is a scalar.
       */
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType 
      augmentJacobianForHomotopy(double conParamValue);

      //@}

      /*! 
       * @name Implementation of LOCA::TimeDependent::AbstractGroup virtual methods.
       */
      //@{

      //! Compute the shifted matrix
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeShiftedMatrix(double alpha, double beta);

      //! Multiply the shifted matrix by a vector.
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrix(const NOX::Abstract::Vector& input,
			 NOX::Abstract::Vector& result) const;

      //! Multiply the shifted matrix by a multi-vector.
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrixMultiVector(
				const NOX::Abstract::MultiVector& input,
				NOX::Abstract::MultiVector& result) const;
      
      /*! 
       * \brief Apply the inverse of the shifted matrix by a multi-vector, as 
       * needed by the shift-and-invert and generalized Cayley transformations.
       */
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrixInverseMultiVector(
			        Teuchos::ParameterList& params, 
				const NOX::Abstract::MultiVector& input,
				NOX::Abstract::MultiVector& result) const;

      //@}

      /*! 
       * @name Implementation of LOCA::Hopf::Moorespence::AbstractGroup virtual methods.
       */
      //@{

      //! Is  \f$J+i\omega B\f$ valid
      /*!
       * The implementation here always returns false.
       */
      virtual bool isComplex() const;

      //! Compute \f$J+i\omega B\f$
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeComplex(double frequency);

      //! Compute \f$(J+i\omega B)(y+iz)\f$
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyComplex(const NOX::Abstract::Vector& input_real,
		   const NOX::Abstract::Vector& input_imag,
		   NOX::Abstract::Vector& result_real,
		   NOX::Abstract::Vector& result_imag) const;

      //! Compute \f$(J+i\omega B)(y+iz)\f$
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyComplexMultiVector(const NOX::Abstract::MultiVector& input_real,
			      const NOX::Abstract::MultiVector& input_imag,
			      NOX::Abstract::MultiVector& result_real,
			      NOX::Abstract::MultiVector& result_imag) const;

      //! Solve \f$(J+i\omega B)(y+iz) = a+ib\f$
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyComplexInverseMultiVector(
			       Teuchos::ParameterList& params,
			       const NOX::Abstract::MultiVector& input_real,
			       const NOX::Abstract::MultiVector& input_imag,
			       NOX::Abstract::MultiVector& result_real,
			       NOX::Abstract::MultiVector& result_imag) const;

      //@}

      /*! 
       * @name Implementation of LOCA::MultiContinuation::AbstractGroup virtual methods.
       */
      //@{

      //! Assignment operator
      virtual void copy(const NOX::Abstract::Group& source);

      //! Set parameters indexed by (integer) paramIDs
      virtual void 
      setParamsMulti(const vector<int>& paramIDs, 
		     const NOX::Abstract::MultiVector::DenseMatrix& vals);

      //! Notify group that the continuation step is completed
      //! The default implementation here is to do nothing.
      virtual void notifyCompletedStep();

      //@}

      /*! 
       * @name Implementation of NOX::Abstract::Group virtual methods.
       */
      //@{

      virtual NOX::Abstract::Group& 
      operator=(const NOX::Abstract::Group& source);

      //@}

    private:

      //! Prohibit generation and use of operator=()
      Group& operator=(const Group& source);

    protected:

      //! Global data
      Teuchos::RefCountPtr<LOCA::GlobalData> globalData;
      
    };
  } // namespace Abstract
} // namespace LOCA

#endif

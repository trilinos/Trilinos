// $Id: LOCA_Abstract_Group.H,v 1.42 2007/06/21 16:22:52 rhoope Exp $ 
// $Source: /space/CVS/Trilinos/packages/nox/src-loca/src/LOCA_Abstract_Group.H,v $ 

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

#ifndef LOCA_ABSTRACT_GROUP_H
#define LOCA_ABSTRACT_GROUP_H

#include "Teuchos_RCP.hpp"

#include "LOCA_Homotopy_AbstractGroup.H"
#include "LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.H"
#include "LOCA_Pitchfork_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_FiniteDifferenceGroup.H"
#include "LOCA_PhaseTransition_AbstractGroup.H"
#include <vector>

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
      public virtual LOCA::Pitchfork::MinimallyAugmented::AbstractGroup,
      public virtual LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup,
      public virtual LOCA::PhaseTransition::AbstractGroup
    {

    public:

      //! Constructor
      Group(const Teuchos::RCP<LOCA::GlobalData>& global_data);

      //! Constructor
      Group(const Teuchos::RCP<LOCA::GlobalData>& global_data,
	    const Teuchos::RCP<LOCA::DerivUtils>& deriv);
 
      //! Copy constructor
      Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~Group();

      /*! 
       * @name Implementation of LOCA::Homotopy::AbstractGroup virtual methods.
       */
      //@{

      /*! 
       * \brief Replace Jacobian \f$J\f$ by \f$aJ+bI\f$ where \f$I\f$ is 
       * the identity matrix and \f$p\f$ is a scalar.
       */
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType 
      augmentJacobianForHomotopy(double a, double b);

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

      //! Compute the second shifted matrix. Can avoid recomputing if two are stored.
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeSecondShiftedMatrix(double alpha, double beta);

      //! Multiply the shifted matrix by a vector.
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType
      applySecondShiftedMatrix(const NOX::Abstract::Vector& input,
			 NOX::Abstract::Vector& result) const;

      //! Multiply the shifted matrix by a multi-vector.
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType
      applySecondShiftedMatrixMultiVector(
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
       * @name Implementation of LOCA::Hopf::MinimallyAugmented::AbstractGroup virtual methods.
       */
      //@{

      /*!
       * Computes conjugate-tranpose matrix vector product 
       * \f$ (J+i\omega B)^H (x + iy) \f$.
       */
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
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
       * @name Implementation of LOCA::MultiContinuation::AbstractGroup virtual methods.
       */
      //@{

      //! Assignment operator
      virtual void copy(const NOX::Abstract::Group& source);

      //! Set parameters indexed by (integer) paramIDs
      virtual void 
      setParamsMulti(const std::vector<int>& paramIDs, 
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

      /*! 
       * @name Implementation of LOCA::PhaseTransition::AbstractGroup virtual methods.
       */
      //@{

      virtual double computeFreeEnergy();

      //@}


    private:

      //! Prohibit generation and use of operator=()
      Group& operator=(const Group& source);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;
      
    };
  } // namespace Abstract
} // namespace LOCA

#endif

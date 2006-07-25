// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_ABSTRACT_GROUP_H
#define LOCA_ABSTRACT_GROUP_H

#include "Teuchos_RefCountPtr.hpp"

// #include "LOCA_Bifurcation_HopfBord_FiniteDifferenceGroup.H"  // base class
// #include "LOCA_Bifurcation_TPBord_SingularSolveGroup.H"       // base class
#include "LOCA_Homotopy_AbstractGroup.H"                      // base class
#include "LOCA_TimeDependent_AbstractGroup.H"
#include "LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.H"// base class
#include "LOCA_Pitchfork_MooreSpence_AbstractGroup.H"// base class

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
//       public LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup,
//       public LOCA::Bifurcation::TPBord::SingularSolveGroup,
      public virtual LOCA::Homotopy::AbstractGroup,
      public LOCA::TimeDependent::AbstractGroup,
      public LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup,
      public LOCA::Pitchfork::MooreSpence::AbstractGroup
    {

    public:

      //! Default constructor.
      Group(const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data);

      //! Default constructor.
      Group(const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	    const Teuchos::RefCountPtr<DerivUtils>& deriv);
  
//       /*! 
//        * \brief Constructor with parameter list 
//        * (for LOCA::Bifurcation::TPBord::SingularSolveGroup).
//        */
//       Group(Teuchos::ParameterList& params,
// 	    const DerivUtils& deriv = LOCA::DerivUtils());
 
      //! Copy constructor
      Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~Group();

      /*! 
       * @name Implementation of LOCA::TimeDependent::AbstractGroup virtual methods.
       */
      //@{

      /*! 
       * \brief Apply the inverse of the shifted matrix by a vector, as 
       * needed by the shift-and-invert and generalized Cayley transformations.
       */
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrixInverse(Teuchos::ParameterList& params, 
				const NOX::Abstract::Vector& input,
                                NOX::Abstract::Vector& result,
                                double shift);

      //@}

      /*! 
       * @name Implementation of LOCA::Bifurcation::HopfBord::AbstractGroup virtual methods.
       */
      //@{

      //! Solve \f$(J+i\omega B)(y+iz) = a+ib\f$
      /*!
       * Implementation here prints an error message and returns 
       * NOX::Abstract::Group::NotDefined.
       */
 //      virtual NOX::Abstract::Group::ReturnType
//       applyComplexInverse(Teuchos::ParameterList& params,
// 			  const NOX::Abstract::Vector& input_real,
// 			  const NOX::Abstract::Vector& input_imag,
// 			  double frequency,
// 			  NOX::Abstract::Vector& result_real,
// 			  NOX::Abstract::Vector& result_imag) const;

      //@}

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

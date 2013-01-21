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

#ifndef LOCA_HOPF_MINIMALLYAUGMENTED_ABSTRACTGROUP_H
#define LOCA_HOPF_MINIMALLYAUGMENTED_ABSTRACTGROUP_H

#include "LOCA_Hopf_MooreSpence_AbstractGroup.H"                // base class

namespace LOCA { 
  
  namespace Hopf { 

    /*! 
     * \brief Groups and vectors for locating Hopf bifurcations
     * using the minimally augmented Hopf formulation.
     */
    namespace MinimallyAugmented {

      /*! 
       * \brief Interface to underlying groups for Hopf calculations
       * using the minimally augmented formulation.
       */
      /*!
       * This abstract class provides the required interface for underlying
       * groups to locate Hopfs using the minimally augmented Hopf (see 
       * LOCA::Hopf::MinimallyAugmented::ExtendedGroup for a 
       * description of the governing equations).  
       *
       * This class is derived from the 
       * LOCA::Hopf::MooreSpence::AbstractGroup 
       * and declares several pure virtual methods to compute various 
       * derivatives of \f$w^H C e\f$ for a given \f$w = w_1 + i w_2\f$ and 
       * \f$e = y + i z\f$ where \f$C = J + i \omega M\f$.  Default 
       * implementations for the derivatives using finite differencing are 
       * implemented in the 
       * LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup.
       */
      class AbstractGroup : 
	public virtual LOCA::Hopf::MooreSpence::AbstractGroup {

      public:
  
	//! Default constructor.
	AbstractGroup() {}

	//! Destructor
	virtual ~AbstractGroup() {}

	/*! 
	 * @name Pure virtual methods 
	 * These methods must be defined by any concrete implementation
	 */
	//@{

	/*!
	 * Computes conjugate-tranpose matrix vector product 
	 * \f$ (J+i\omega M)^H (x + iy) \f$.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	applyComplexTranspose(const NOX::Abstract::Vector& input_real,
			      const NOX::Abstract::Vector& input_imag,
			      NOX::Abstract::Vector& result_real,
			      NOX::Abstract::Vector& result_imag) const = 0;

	/*!
	 * Computes conjugate-tranpose matrix vector product 
	 * \f$ (J+i\omega M)^H (x + iy) \f$.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	applyComplexTransposeMultiVector(
			   const NOX::Abstract::MultiVector& input_real,
			   const NOX::Abstract::MultiVector& input_imag,
			   NOX::Abstract::MultiVector& result_real,
			   NOX::Abstract::MultiVector& result_imag) const = 0;

	//! Solve \f$(J+i\omega M)^H (x + iy) = a+ib\f$
	virtual NOX::Abstract::Group::ReturnType
	applyComplexTransposeInverseMultiVector(
			    Teuchos::ParameterList& params,
			    const NOX::Abstract::MultiVector& input_real,
			    const NOX::Abstract::MultiVector& input_imag,
			    NOX::Abstract::MultiVector& result_real,
			    NOX::Abstract::MultiVector& result_imag) const = 0;

	/*! 
	 * \brief Computes the derivative \f$\partial w^H C e/\partial p\f$,
	 * \f$w = w_1 + i w_2\f$, \f$e = y + i z\f$, \f$C = J + i \omega M\f$.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDwtCeDp(const std::vector<int>& paramIDs, 
		       const NOX::Abstract::Vector& w1,
		       const NOX::Abstract::Vector& w2,
		       const NOX::Abstract::Vector& y,
		       const NOX::Abstract::Vector& z,
		       double omega,
		       NOX::Abstract::MultiVector::DenseMatrix& result_real,
		       NOX::Abstract::MultiVector::DenseMatrix& result_imag,
		       bool isValid) = 0;

	/*! 
	 * \brief Computes the derivative 
	 * \f$\frac{\partial w^H C e}{\partial x}\f$,
	 * \f$w = w_1 + i w_2\f$, \f$e = y + i z\f$, \f$C = J + i \omega M\f$.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	computeDwtCeDx(const NOX::Abstract::Vector& w1,
		       const NOX::Abstract::Vector& w2,
		       const NOX::Abstract::Vector& y,
		       const NOX::Abstract::Vector& z,
		       double omega, 
		       NOX::Abstract::Vector& result_real,
		       NOX::Abstract::Vector& result_imag) = 0;

	//@}

	//! Bring NOX::Abstract::Group::operator=() into scope
	using NOX::Abstract::Group::operator=;

      }; // class AbstractGroup

    } // namespace MinimallyAugmented

  } // namespace Hopf

} // namespace LOCA

#endif

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

#ifndef LOCA_MULTICONTINUATION_FINITEDIFFERENCEGROUP_H
#define LOCA_MULTICONTINUATION_FINITEDIFFERENCEGROUP_H

#include "LOCA_MultiContinuation_AbstractGroup.H"    // base class
#include "LOCA_DerivUtils.H"                    // class data element

namespace LOCA { 

  namespace MultiContinuation { 

    /*!  
     * \brief Concrete class that provides a concrete implementation of the
     * computeDfDp() method of the LOCA::Continuation::AbstractGroup using
     * first-order finite differencing.
    */
    /*! 
     * The finite-differencing calculations are actually implemented by
     * the LOCA::DerivUtils class, and a custom DerivUtils object can be
     * set by the setDerivUtils() method.  However, in the future
     * the calculations encapsulated in the DerivUtils class may be 
     * incorporated directly into this class and other finite-differencing
     * child classes.
     */
    class FiniteDifferenceGroup : 
      public virtual LOCA::MultiContinuation::AbstractGroup {

    public:
  
      //! Constructor.
      FiniteDifferenceGroup();
 
      //! Copy constructor
      FiniteDifferenceGroup(const FiniteDifferenceGroup& source, 
			    NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~FiniteDifferenceGroup();

      //! Copy
      virtual void copy(const NOX::Abstract::Group& source);

      //! Assignment
      virtual NOX::Abstract::Group& 
      operator=(const NOX::Abstract::Group& source);

      //! Set the LOCA::DerivUtils object
      virtual void setDerivUtils(
			 const Teuchos::RCP<LOCA::DerivUtils>& deriv);

      /*! 
       * Compute \f$\partial F/\partial p\f$ for each parameter \f$p\f$ 
       * indexed by paramIDs.  The first column of \em dfdp holds F,
       * which is valid if \em isValidF is true.  Otherwise F must be
       * computed.
       */
      /*! 
       * The calculation is implemented by calling the corresponding
       * LOCA::DerivUtils::computeDfDp() method of the passed 
       * LOCA::DerivUtils object.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeDfDpMulti(const std::vector<int>& paramIDs, 
		       NOX::Abstract::MultiVector& dfdp, 
		       bool isValidF);

    private:

      //! Prohibit generation and use of operator=()
      FiniteDifferenceGroup& operator=(const FiniteDifferenceGroup&);

    protected:

      //! Pointer to current DerivUtils derivative computation object
      Teuchos::RCP<DerivUtils> derivPtr;

    };
  } // namespace MultiContinuation
} // namespace LOCA

#endif

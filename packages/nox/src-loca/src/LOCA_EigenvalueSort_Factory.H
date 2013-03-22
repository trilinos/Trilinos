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

#ifndef LOCA_EIGENVALUESORT_FACTORY_H
#define LOCA_EIGENVALUESORT_FACTORY_H

#include "Teuchos_RCP.hpp"

// forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace EigenvalueSort {
    class AbstractStrategy;
  }
}

namespace LOCA {

  namespace EigenvalueSort {

    //! Factory for creating %EigenvalueSort strategy objects
    /*!
     * The parameters passed to the create() through the \c eigenParams
     * argument method should specify the sorting method "Sorting Method" 
     * as described below, as well as any additional parameters
     * for the particular strategy.
     * <ul>
     * <li> "Sorting Order" - Name of the sorting method. Valid choices are
     *   <ul>
     *   <li> "LM" (LOCA::EigenvalueSort::LargestMagnitude) [Default]
     *   <li> "LR" (LOCA::EigenvalueSort::LargestReal)
     *   <li> "LI" (LOCA::EigenvalueSort::LargestImaginary)
     *   <li> "SM" (LOCA::EigenvalueSort::SmallestMagnitude)
     *   <li> "SR" (LOCA::EigenvalueSort::SmallestReal)
     *   <li> "SI" (LOCA::EigenvalueSort::SmallestImaginary)
     *   <li> "CA" (LOCA::EigenvalueSort::LargestRealInverseCayley)
     *   </ul>
     * </ul>
     */
    class Factory {

    public:

      //! Constructor
      Factory(const Teuchos::RCP<LOCA::GlobalData>& global_data);

      //! Destructor
      virtual ~Factory();

      //! Create sorting strategy
      /*!
       * \param topParams [in] Parsed top-level parameter list.
       * \param eigenParams [in] %Eigensolver parameters as described above
       */
      Teuchos::RCP<LOCA::EigenvalueSort::AbstractStrategy>
      create(
	const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RCP<Teuchos::ParameterList>& eigenParams);

      //! Return strategy name given by \c eigenParams
      const std::string& strategyName(Teuchos::ParameterList& eigenParams) const;

    private:

      //! Private to prohibit copying
      Factory(const Factory&);

      //! Private to prohibit copying
      Factory& operator = (const Factory&);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

    }; // Class Factory

  } // Namespace EigenvalueSort

} // Namespace LOCA

#endif

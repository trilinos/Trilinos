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

#ifndef ANASAZI_LOCA_SORT_HPP
#define ANASAZI_LOCA_SORT_HPP

#include "Teuchos_RCP.hpp"

#include "AnasaziSortManager.hpp"

// Forward declarations
namespace NOX {
  namespace Abstract {
    class MultiVector;
  }
}
namespace LOCA {
  class GlobalData;
  namespace AnasaziOperator {
    class AbstractStrategy;
  }
  namespace EigenvalueSort {
    class AbstractStrategy;
  }
}

namespace Anasazi {

  /*!
   * \brief An implementation of the Anasazi::SortManager that performs a
   * collection of common sorting techniques.
   */
  class LOCASort :
    public Anasazi::SortManager< double > {

    //! Anasazi::Multivec typedef
    typedef NOX::Abstract::MultiVector MV;

    //! Anasazi::Operator typedef
    typedef LOCA::AnasaziOperator::AbstractStrategy OP;

  public:

    //! Constructor
    /*!
     * @param global_data [in] Global data object
     * @param strategy_ [in] Ref-count ptr to LOCA sorting strategy
     */
    LOCASort(const Teuchos::RCP<LOCA::GlobalData>& global_data,
         const Teuchos::RCP<LOCA::EigenvalueSort::AbstractStrategy>& strategy_);

    //! Destructor
    virtual ~LOCASort();

    //! Sort the vector of eigenvalues with respect to the chosen sorting type, optionally returning the permutation vector.
    /*!
     * \param solver [in] Eigensolver that is calling the sorting routine
     * \param n [in] Size of the array
     * \param evals [in/out] Array of length n containing the eigenvalues to be sorted
     * \param perm [out] Vector of length n to store the permutation (optional)
     * \return Returns the status of the sorting routine
     */
    virtual void
    sort(std::vector<double>& evals, Teuchos::RCP<std::vector<int> > perm = Teuchos::null, int n = -1) const;

    /*!
     * \brief Sort the vectors of eigenpairs with respect to the chosen
     * sorting type, optionally returning the permutation vector.
     */
    /*!
     * \param solver [in] Eigensolver that is calling the sorting routine
     * \param n [in] Size of the array
     * \param r_evals [in/out] Array of length n containing the real part of
     * the eigenvalues to be sorted
     * \param i_evals [in/out] Array of length n containing the imaginary
     * part of the eigenvalues to be sorted
     * \param perm [out] Vector of length n to store the permutation (optional)
     * \return Returns the status of the sorting routine
     */
    virtual void
    sort(std::vector<double>& r_evals, std::vector<double>& i_evals, Teuchos::RCP<std::vector<int> > perm = Teuchos::null, int n = -1) const;

  protected:

    //! Global data
    Teuchos::RCP<LOCA::GlobalData> globalData;

    //! EigenvalueSort strategy
    Teuchos::RCP<LOCA::EigenvalueSort::AbstractStrategy> strategy;

  }; // Class LOCASort

} // Namespace Anasazi

#endif // ANASAZI_LOCA_SORT_H


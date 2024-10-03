// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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


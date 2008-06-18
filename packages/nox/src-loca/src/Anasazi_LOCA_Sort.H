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


// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_ORTHOGONALIZATION_FACTORY_HPP
#define STOKHOS_ORTHOGONALIZATION_FACTORY_HPP

#include <string>
#include "Teuchos_Array.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Stokhos_SDMUtils.hpp"

namespace Stokhos {

  /*! 
   * \brief Encapsulate various orthogonalization (ie QR) methods
   */
  template <typename ordinal_type, typename value_type>
  class OrthogonalizationFactory {
  public:

    //! Constructor
    /*!
     * \param params Parameter dictating choice of reduction method
     */
    OrthogonalizationFactory() {}

    //! Destructor
    virtual ~OrthogonalizationFactory() {}

    typedef Teuchos::SerialDenseMatrix<ordinal_type,value_type> SDM;

    //! Create orthogonal basis via the method indicated by \c method
    static ordinal_type
    createOrthogonalBasis(const std::string& method, value_type threshold,
			  bool verbose, const SDM& A,
			  const Teuchos::Array<value_type>& w,
			  SDM& Q, SDM& R,
			  Teuchos::Array<ordinal_type>& piv) {

      ordinal_type m = A.numRows();
      ordinal_type n = A.numCols();
      ordinal_type rank = std::min(m,n);

      if (method == "SVD") { // A = U*diag(sigma)*V^T, Q = U, R = sigma*V^T
	Teuchos::Array<value_type> sigma;
	SDM Vt;
	rank = svd_threshold(threshold, A, sigma, Q, Vt);
	R.reshape(rank, Vt.numCols());
	for (ordinal_type j=0; j<Vt.numCols(); j++)
	  for (ordinal_type i=0; i<rank; i++)
	    R(i,j) = sigma[i]*Vt(i,j);
	piv.resize(n);
	for (int i=0; i<n; i++)
	  piv[i] = i;

	if (verbose) {
	  std::cout << "diag(sigma) = [ ";
	  for (ordinal_type i=0; i<rank; i++)
	    std::cout << sigma[i] << " ";
	  std::cout << "]" << std::endl;
	  
	  std::cout << "rank = " << rank << std::endl;
	}
      }

      else { // All QR-based methods

	if (method == "Householder")
	  rank = CPQR_Householder_threshold(threshold, A, w, Q, R, piv);

	else if (method == "Householder without Pivoting") {
	  QR_Householder(rank, A, w, Q, R);
	  piv.resize(n);
	  for (int i=0; i<n; i++)
	    piv[i] = i;
	}

	else if (method == "Modified Gram-Schmidt")
	  rank = CPQR_MGS_threshold(threshold, A, w, Q, R, piv);

	else if (method == "Modified Gram-Schmidt with Reorthogonalization")
	  rank = CPQR_MGS_reorthog_threshold(threshold, A, w, Q, R, piv);

	else if (method == "Modified Gram-Schmidt without Pivoting") {
	  QR_MGS(rank, A, w, Q, R);
	  piv.resize(n);
	  for (int i=0; i<n; i++)
	    piv[i] = i;
	}

	else if (method == "Modified Gram-Schmidt without Pivoting with Reorthogonalization") {
	  QR_MGS2(rank, A, w, Q, R);
	  piv.resize(n);
	  for (int i=0; i<n; i++)
	    piv[i] = i;
	}

	else
	  TEUCHOS_TEST_FOR_EXCEPTION(
	  true, std::logic_error, 
	  "Invalid orthogonalization method " << method);

	if (verbose) {
	  std::cout << "piv = [";
	  for (ordinal_type i=0; i<rank; i++)
	    std::cout << piv[i] << " ";
	  std::cout << "]" << std::endl;
    
	  std::cout << "diag(R) = [ ";
	  for (ordinal_type i=0; i<rank; i++)
	    std::cout << R(i,i) << " ";
	  std::cout << "]" << std::endl;
	  
	  std::cout << "rank = " << rank << std::endl;

	  // Check A*P = Q*R
	  std::cout << "||A*P-Q*R||_infty = " 
		    << Stokhos::residualCPQRError(A,Q,R,piv) << std::endl;
      
	  // Check Q^T*diag(w)*Q = I
	  std::cout << "||I - Q^T*diag(w)**Q||_infty = " 
		    << weightedQROrthogonalizationError(Q,w) << std::endl;
	}
      }

      return rank;
    }

  private:

    // Prohibit copying
    OrthogonalizationFactory(const OrthogonalizationFactory&);

    // Prohibit Assignment
    OrthogonalizationFactory& operator=(const OrthogonalizationFactory&);

  }; // class OrthogonalizationFactory

} // Namespace Stokhos

#endif

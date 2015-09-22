// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_EPETRA_OPERATOR_ORTHOG_POLY_HPP
#define STOKHOS_EPETRA_OPERATOR_ORTHOG_POLY_HPP

#include "Stokhos_VectorOrthogPoly.hpp"
#include "Stokhos_ProductEpetraOperator.hpp"

namespace Stokhos {

  /*! 
   * \brief A container class storing an orthogonal polynomial whose
   * coefficients are vectors, operators, or in general any type that 
   * would have an expensive copy constructor.  
   */
  class EpetraOperatorOrthogPoly : 
    public virtual VectorOrthogPoly<Epetra_Operator>,
    public virtual ProductEpetraOperator {
  public:

    //! Typename of values
    typedef double value_type;

    //! Typename of ordinals
    typedef int ordinal_type;

    /*! 
     * \brief Create a polynomial for basis \c basis where each coefficient is 
     * generated from the supplied map.
     */
    EpetraOperatorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_Map>& domain_base_map,
      const Teuchos::RCP<const Epetra_Map>& range_base_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm);

    /*! 
     * \brief Create a polynomial for basis \c basis where each coefficient is 
     * given by a created block vector
     */
    EpetraOperatorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_Map>& domain_base_map,
      const Teuchos::RCP<const Epetra_Map>& range_base_map,
      const Teuchos::RCP<const Epetra_Map>& range_product_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm);
    
    //! Copy constructor
    /*!
     * NOTE:  This is a shallow copy
     */
    EpetraOperatorOrthogPoly(const EpetraOperatorOrthogPoly& v);

    //! Destructor
    virtual ~EpetraOperatorOrthogPoly();

    //! Assignment
    /*!
     * NOTE:  This is a shallow copy
     */
    EpetraOperatorOrthogPoly& operator=(const EpetraOperatorOrthogPoly& v);

  protected:

    //! Protected constructor to allow 2-stage derived setup
    EpetraOperatorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm);

    //! Second stage of setup
    void setup(const Teuchos::RCP<const Epetra_Map>& domain_base_map,
	       const Teuchos::RCP<const Epetra_Map>& range_base_map);

  }; // class EpetraOperatorOrthogPoly

} // end namespace Stokhos

#endif  // STOKHOS_EPETRA_OPERATOR_ORTHOG_POLY_HPP

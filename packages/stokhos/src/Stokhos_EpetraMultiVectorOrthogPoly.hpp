// $Id$ 
// $Source$ 
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

#ifndef STOKHOS_EPETRAMULTIVECTORORTHOGPOLY_HPP
#define STOKHOS_EPETRAMULTIVECTORORTHOGPOLY_HPP

#include "Stokhos_VectorOrthogPoly.hpp"
#include "Stokhos_VectorOrthogPolyTraitsEpetra.hpp"
#include "EpetraExt_BlockMultiVector.h"

namespace Stokhos {

  /*! 
   * \brief A container class storing an orthogonal polynomial whose
   * coefficients are vectors, operators, or in general any type that 
   * would have an expensive copy constructor.  
   */
  class EpetraMultiVectorOrthogPoly : 
    public VectorOrthogPoly<Epetra_MultiVector> {
  public:

    //! Typename of values
    typedef double value_type;

    //! Typename of ordinals
    typedef int ordinal_type;

    //! Constructor with no basis
    /*!
     * Use with care!  Generally you will want to call reset() before using
     * any of the methods on this class.
     */
    EpetraMultiVectorOrthogPoly();

    /*! 
     * \brief Create a polynomial for basis \c basis with empty 
     * coefficients
     */
    EpetraMultiVectorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis);

    /*! 
     * \brief Create a polynomial for basis \c basis with empty 
     * coefficients of size sz
     */
    EpetraMultiVectorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      int sz);

    /*! 
     * \brief Create a polynomial for basis \c basis where each coefficient is 
     * generated from the supplied map.
     */
    EpetraMultiVectorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Epetra_BlockMap& coeff_map,
      int num_vectors);

    /*! 
     * \brief Create a polynomial for basis \c basis where each coefficient is 
     * generated from the supplied map.
     */
    EpetraMultiVectorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      const Epetra_BlockMap& coeff_map,
      int num_vectors,
      int sz);

    /*! 
     * \brief Create a polynomial for basis \c basis where each coefficient is 
     * given by a created block vector
     */
    EpetraMultiVectorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      Epetra_DataAccess CV,
      const Epetra_BlockMap& coeff_map,
      const Epetra_BlockMap& block_map,
      int num_vectors);

    /*! 
     * \brief Create a polynomial for basis \c basis where each coefficient is 
     * given by the supplied block vector.
     */
    EpetraMultiVectorOrthogPoly(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& basis,
      Epetra_DataAccess CV,
      const Epetra_BlockMap& coeff_map,
      const Epetra_MultiVector& block_vector);
    
    //! Copy constructor
    /*!
     * NOTE:  This is a shallow copy
     */
    EpetraMultiVectorOrthogPoly(const EpetraMultiVectorOrthogPoly& v);

    //! Destructor
    ~EpetraMultiVectorOrthogPoly();

    //! Assignment
    /*!
     * NOTE:  This is a shallow copy
     */
    EpetraMultiVectorOrthogPoly& operator=(const EpetraMultiVectorOrthogPoly& v);

    //! Assignment
    EpetraMultiVectorOrthogPoly& operator=(const Epetra_MultiVector& v);

    //! Assignment
    void assignToBlockMultiVector(Epetra_MultiVector& v) const;

    //! Assignment
    void assignFromBlockMultiVector(const Epetra_MultiVector& v);
      

    //! Reset to a new basis
    /*!
     * This resizes array to fit new basis.
     */
    void reset(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> >& new_basis,
      const Epetra_BlockMap& coeff_map,
      int num_vectors);

    //! Reset vector cofficients
    void resetCoefficients(Epetra_DataAccess CV,
			   const Epetra_BlockMap& coeff_map,
			   const Epetra_MultiVector& block_vector);

    //! Get block vector
    Teuchos::RCP<EpetraExt::BlockMultiVector> getBlockMultiVector();

    //! Get block vector
    Teuchos::RCP<const EpetraExt::BlockMultiVector> getBlockMultiVector() const;

    //! Set block vector
    void setBlockMultiVector(
      const Teuchos::RCP<EpetraExt::BlockMultiVector>& block_vec);

  protected:

    //! Block vector storing coefficients
    Teuchos::RCP<EpetraExt::BlockMultiVector> bv;

  }; // class EpetraMultiVectorOrthogPoly

} // end namespace Stokhos

#endif  // STOKHOS_EPETRAMULTIVECTORORTHOGPOLY_HPP

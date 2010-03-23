// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_TPETRA_VECTOR_DECL_HPP
#define THYRA_TPETRA_VECTOR_DECL_HPP


#include "Thyra_SpmdVectorBase_decl.hpp"
#include "Thyra_TpetraVectorSpace_decl.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"


namespace Thyra {


/** \brief Concrete Thyra::SpmdVectorBase using Tpetra::Vector.
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class TpetraVector : virtual public SpmdVectorBase<Scalar> {
public:

  /** @name Constructors/initializers */
  //@{

  /** \brief Construct to uninitialized. */
  TpetraVector();

  /** \brief Initialize. */
  void initialize(
    const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
    const RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVector
    );

  /** \brief Initialize. */
  void constInitialize(
    const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
    const RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVector
    );

  /** \brief Get the embedded non-const Tpetra::Vector. */
  RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getTpetraVector();

  /** \brief Get the embedded non-const Tpetra::Vector. */
  RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getConstTpetraVector() const;

  //@}

  /** @name Overridden from SpmdVectorBase */
  //@{

  /** \brief . */
  RCP<const SpmdVectorSpaceBase<Scalar> > spmdSpace() const;
  /** \brief . */
  void getNonconstLocalDataImpl(const Ptr<ArrayRCP<Scalar> > &localValues);
  /** \brief . */
  void getLocalDataImpl(const Ptr<ArrayRCP<const Scalar> > &localValues) const;

  //@}

private:

  // ///////////////////////////////////////
  // Private data members
  
  RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  tpetraVectorSpace_;
  
  Teuchos::ConstNonconstObjectContainer<Tpetra::Vector<Scalar, LocalOrdinal,GlobalOrdinal,Node> >
  tpetraVector_;

  typedef Tpetra::MultiVector<Scalar, LocalOrdinal,GlobalOrdinal,Node> TpetraMultiVector_t;

  // ////////////////////////////////////
  // Private member functions

  template<class TpetraVector_t>
  void initializeImpl(
    const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
    const RCP<TpetraVector_t> &tpetraVector
    );

};


/** \brief Nonmember constructor for TpetraVector.
 *
 * \relates TpetraVector
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
inline
RCP<TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
tpetraVector(
  const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
  const RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVector
  )
{
  RCP<TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > v =
    Teuchos::rcp(new TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>);
  v->initialize(tpetraVectorSpace, tpetraVector);
  return v;
}


/** \brief Nonmember constructor for TpetraVector.
 *
 * \relates TpetraVector
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
inline
RCP<const TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
constTpetraVector(
  const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
  const RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVector
  )
{
  RCP<TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > v =
    Teuchos::rcp(new TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>);
  v->constInitialize(tpetraVectorSpace, tpetraVector);
  return v;
}


} // end namespace Thyra


#endif // THYRA_TPETRA_VECTOR_DECL_HPP

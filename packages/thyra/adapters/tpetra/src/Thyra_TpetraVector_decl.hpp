// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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

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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov)
//
// ***********************************************************************
// @HEADER


#ifndef THYRA_TPETRA_VECTOR_SPACE_DECL_HPP
#define THYRA_TPETRA_VECTOR_SPACE_DECL_HPP


#include "Thyra_SpmdVectorSpaceDefaultBase.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"


namespace Thyra {


/** \brief Concrete implementation of an SPMD vector space for Tpetra.
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class TpetraVectorSpace : public SpmdVectorSpaceDefaultBase<Scalar>
{
public:

  /** \brief . */
  typedef TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> this_t;

  /** @name Constructors and initializers */
  //@{

  /** \brief Create with weak ownership to self. */
  static RCP<TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > create();

  /** \brief Initialize a serial space. */
  void initialize(
    const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > &tpetraMap
    );

  //@}

  /** @name Public overridden from VectorSpaceBase */
  //@{
  /** \brief Returns true if all the elements in <tt>rng</tt> are in this
   * process.
   */
  bool hasInCoreView(
    const Range1D& rng, const EViewType viewType, const EStrideType strideType
    ) const;
  /** \brief . */
  RCP< const VectorSpaceBase<Scalar> > clone() const;

  /** \brief Get the embedded Tpetra::Map. */
  RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
  getTpetraMap() const;

  //@}

protected:

  /** @name Protected overridden from VectorSpaceBase */
  //@{

  /** \brief . */
  RCP<VectorBase<Scalar> >
  createMember() const;
  /** \brief . */
  RCP<MultiVectorBase<Scalar> >
  createMembers(int numMembers) const;

public:

  RCP<MultiVectorBase<Scalar> >
  createCachedMembersView( const RTOpPack::SubMultiVectorView<Scalar> &raw_mv ) const;

  RCP<const MultiVectorBase<Scalar> >
  createCachedMembersView( const RTOpPack::ConstSubMultiVectorView<Scalar> &raw_mv ) const;

  //@}

public:

  /** @name Public overridden from SpmdVectorSpaceDefaultBase */
  //@{

  /** \brief . */
  RCP<const Teuchos::Comm<Ordinal> > getComm() const;
  /** \brief . */
  Ordinal localSubDim() const;

  //@}

private:

  // //////////////////////////////////////
  // Private data members

  RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > tpetraMap_;
  mutable RCP<TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tpetraDomainSpace_;
  mutable RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tpetraMV_;
  // The only reason Thyra needs this comm_ object is because Thyra
  // uses Ordinal as the Comm template parameter, while Tpetra uses
  // int.  Ordinal is some 64-bit type, which doesn't make any sense,
  // given that MPI implementations currently only allow 32-bit
  // process ranks.  This is why Thyra does not just use the Map's
  // stored communicator.
  RCP<const Teuchos::Comm<Ordinal> > comm_;
  RCP<this_t> weakSelfPtr_;

  // /////////////////////////////////////
  // Private member functions

  TpetraVectorSpace();

}; // end class TpetraVectorSpace


/** \brief Nonmember constructor that creats a serial vector space.
 *
 * \relates TpetraVectorSpace
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
tpetraVectorSpace(
  const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > &tpetraMap
  )
{
  RCP<TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > vs =
    TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node>::create();
  vs->initialize(tpetraMap);
  return vs;
}


} // end namespace Thyra


#endif // THYRA_TPETRA_VECTOR_SPACE_DECL_HPP

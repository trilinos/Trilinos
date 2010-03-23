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


#ifndef THYRA_TPETRA_VECTOR_SPACE_DECL_HPP
#define THYRA_TPETRA_VECTOR_SPACE_DECL_HPP


#include "Thyra_SpmdVectorSpaceDefaultBase.hpp"
#include "Tpetra_Map.hpp"


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
  RCP<const Teuchos::Comm<Ordinal> > comm_;
  Ordinal localSubDim_;
  int numProc_;
  int procRank_;
  RCP<this_t> weakSelfPtr_;

  // /////////////////////////////////////
  // Private member functions
 
  TpetraVectorSpace();
 
}; // end class TpetraVectorSpace


/** \brief Nonmember consturctor that creats a serial vector space.
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

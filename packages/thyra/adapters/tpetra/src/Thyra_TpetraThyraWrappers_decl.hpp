// @HEADER
// ***********************************************************************
// 
//               Thyra: Trilinos Solver Framework Core
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

#ifndef THYRA_TPETRA_THYRA_WRAPPERS_DECL_HPP
#define THYRA_TPETRA_THYRA_WRAPPERS_DECL_HPP


#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_TpetraVectorSpace_decl.hpp"
#include "Thyra_TpetraVector_decl.hpp"
#include "Thyra_TpetraMultiVector_decl.hpp"
#include "Thyra_TpetraLinearOp_decl.hpp"


namespace Thyra {


/** \brief Given an Tpetra <tt>Teuchos::Comm<int></tt> object, return an
 * equivalent <tt>Teuchos::Comm<Ordinal></tt> object.
 *
 * Will throw if conversion is not successful.
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
RCP<const Teuchos::Comm<Ordinal> >
convertTpetraToThyraComm( const RCP<const Teuchos::Comm<int> > &tpetraComm );


/** \brief Create a Thyra::VectorSpaceBase object given a Tpetra::Map.
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const VectorSpaceBase<Scalar> >
createVectorSpace(const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > &tpetraMap);


/** \brief .
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<VectorBase<Scalar> >
createVector(
  const RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVector,
  const RCP<const VectorSpaceBase<Scalar> > space = Teuchos::null
  );


/** \brief .
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const VectorBase<Scalar> >
createConstVector(
  const RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVector,
  const RCP<const VectorSpaceBase<Scalar> > space = Teuchos::null
  );


/** \brief .
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MultiVectorBase<Scalar> >
createMultiVector(
  const RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraMultiVector,
  const RCP<const VectorSpaceBase<Scalar> > rangeSpace = Teuchos::null,
  const RCP<const VectorSpaceBase<Scalar> > domainSpace = Teuchos::null
  );


/** \brief .
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const MultiVectorBase<Scalar> >
createConstMultiVector(
  const RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraMultiVector,
  const RCP<const VectorSpaceBase<Scalar> > rangeSpace = Teuchos::null,
  const RCP<const VectorSpaceBase<Scalar> > domainSpace = Teuchos::null
  );


/** \brief .
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<LinearOpBase<Scalar> >
createLinearOp(
  const RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraOperator,
  const RCP<const VectorSpaceBase<Scalar> > rangeSpace = Teuchos::null,
  const RCP<const VectorSpaceBase<Scalar> > domainSpace = Teuchos::null
  );


/** \brief .
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const LinearOpBase<Scalar> >
createConstLinearOp(
  const RCP<const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraOperator,
  const RCP<const VectorSpaceBase<Scalar> > rangeSpace = Teuchos::null,
  const RCP<const VectorSpaceBase<Scalar> > domainSpace = Teuchos::null
  );


/** \brief Traits class that enables the extraction of Tpetra operator/vector
 * objects wrapped in Thyra operator/vector objects.
 *
 * Example usage:
 
 \code

  typedef Thyra::TpetraObjectExtraction<Scalar,LO,GO,Node> TOE;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> TpetraMultiVector_t;
  
  RCP<TpetraMultiVector_t> tpetraMv = TOE::getTpetraMultiVector(thyraMv);
  RCP<TpetraVector_t> tpetraV = TOE::getTpetraVector(thyraV);

 \endcode

 *
 * \todo Finish documentation
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal=LocalOrdinal,
  class Node=Kokkos::DefaultNode::DefaultNodeType >
class TpetraOperatorVectorExtraction {
public:

  // ToDo: Get a Tpetra::Map from a Thyra::VectorSpaceBase?

  /** \brief Get a non-const Tpetra::Vector from a non-const
   * Thyra::VectorBase object.
   */
  static RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getTpetraVector(const RCP<VectorBase<Scalar> > &v);

  /** \brief Get a const Tpetra::Vector from a const
   * Thyra::VectorBase object.
   */
  static RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getConstTpetraVector(const RCP<const VectorBase<Scalar> > &v);

  /** \brief Get a non-const Tpetra::MultiVector from a non-const
   * Thyra::MultiVectorBase object.
   */
  static RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getTpetraMultiVector(const RCP<MultiVectorBase<Scalar> > &mv);

  /** \brief Get a const Tpetra::MultiVector from a const
   * Thyra::MultiVectorBase object.
   */
  static RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getConstTpetraMultiVector(const RCP<const MultiVectorBase<Scalar> > &mv);

  /** \brief Get a non-const Tpetra::Operator from a non-const
   * Thyra::LinearOpBase object.
   */
  static RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getTpetraOperator(const RCP<LinearOpBase<Scalar> > &op);

  /** \brief Get a const Tpetra::Operator from a const
   * Thyra::LinearOpBase object.
   */
  static RCP<const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getConstTpetraOperator(const RCP<const LinearOpBase<Scalar> > &op);

};


} // namespace Thyra


#endif // THYRA_TPETRA_THYRA_WRAPPERS_DECL_HPP

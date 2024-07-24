// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
template <class Scalar = Tpetra::Vector<>::scalar_type,
          class LocalOrdinal = Tpetra::Vector<>::local_ordinal_type,
          class GlobalOrdinal = Tpetra::Vector<>::global_ordinal_type,
          class Node = Tpetra::Vector<>::node_type>
class TpetraOperatorVectorExtraction {
public:

  /** \brief Get a const Tpetra::Map from a const Thyra::VectorSpaceBase object.
   */
  static RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
  getTpetraMap(const RCP<const VectorSpaceBase<Scalar> > &vs);

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

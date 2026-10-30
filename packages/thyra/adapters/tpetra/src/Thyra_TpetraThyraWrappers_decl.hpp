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


#include <type_traits>
#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_TpetraVectorSpace_decl.hpp"
#include "Thyra_TpetraVector_decl.hpp"
#include "Thyra_TpetraMultiVector_decl.hpp"
#include "Thyra_TpetraLinearOp_decl.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"

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

/** \brief Given an Tpetra <tt>Teuchos::Comm<Ordinal></tt> object, return an
 * equivalent <tt>Teuchos::Comm<int></tt> object.
 *
 * Will throw if conversion is not successful.
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
RCP<const Teuchos::Comm<int> >
convertThyraToTpetraComm( const RCP<const Teuchos::Comm<Ordinal> > &thyraComm );


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

/** \brief Initialize a preconditioner from a forward linear operator.
 *
*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Thyra::PreconditionerBase<Scalar>>
initializePrec(
  const PreconditionerFactoryBase<Scalar> &precFactory,
  const RCP<const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &tpetraFwdOp,
  const Teuchos::RCP<PreconditionerBase<Scalar>> &prec = Teuchos::null,
  const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED);

/** \brief Initialize a preconditioner from a forward linear operator.
 *
*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Thyra::PreconditionerBase<Scalar>>
initializePrec(
  const PreconditionerFactoryBase<Scalar> &precFactory,
  const RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &tpetraFwdOp,
  const Teuchos::RCP<PreconditionerBase<Scalar>> &prec = Teuchos::null,
  const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED);

/** \brief Create and initialize a <tt>LinearOpWithSolveBase</tt> object from
 * a <tt>Tpetra::Operator</tt> object using a
 * <tt>LinearOpWithSolveFactoryBase</tt> strategy object.
 * \relates LinearOpWithSolveFactoryBase */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<LinearOpWithSolveBase<Scalar>>
linearOpWithSolve(
  const LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const RCP<const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &tpetraFwdOp,
  const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED);


/** \brief Create and initialize a <tt>LinearOpWithSolveBase</tt> object from
 * a <tt>Tpetra::CrsMatrix</tt> object using a
 * <tt>LinearOpWithSolveFactoryBase</tt> strategy object.
 * \relates LinearOpWithSolveFactoryBase */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<LinearOpWithSolveBase<Scalar>>
linearOpWithSolve(
  const LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &tpetraFwdOp,
  const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED);

/** \brief Initialize a preconditioned LOWSB object given an external
 * preconditioner.
 *
 * \relates LinearOpWithSolveFactoryBase
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<LinearOpWithSolveBase<Scalar>>
initializePreconditionedOp(
  const LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const RCP<const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &tpetraFwdOp,
  const RCP<PreconditionerBase<Scalar>> &prec,
  const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED);

/** \brief Initialize a preconditioned LOWSB object given an external
 * preconditioner.
 *
 * \relates LinearOpWithSolveFactoryBase
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<LinearOpWithSolveBase<Scalar>>
initializePreconditionedOp(
  const LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &tpetraFwdOp,
  const RCP<PreconditionerBase<Scalar>> &prec,
  const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED);

/** \brief Call <tt>solve()</tt> as a non-member function
 *
 * \relates LinearOpWithSolveBase
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
SolveStatus<Scalar>
solve(
  const LinearOpWithSolveBase<Scalar> &A, const EOpTransp A_trans,
  const RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &tpetraB,
  const RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &tpetraX,
  const Ptr<const SolveCriteria<Scalar>> solveCriteria = Teuchos::null);

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

  /** \brief Get a non-const Tpetra::Vector from a non-const
   * Thyra::VectorBase object that is a Thyra::ProductVectorBase.
   */
  static RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getTpetraVectorBlock(const RCP<VectorBase<Scalar> > &v, const int k);

  /** \brief Get a const Tpetra::Vector from a const
   * Thyra::VectorBase object that is a Thyra::ProductVectorBase.
   */
  static RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getConstTpetraVectorBlock(const RCP<const VectorBase<Scalar> > &v, const int k);

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

  /** \brief Get a non-const Tpetra::MultiVector from a non-const
   * Thyra::MultiVectorBase object that is a Thyra::ProductMultiVectorBase.
   */
  static RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getTpetraMultiVectorBlock(const RCP<MultiVectorBase<Scalar> > &mv, const int k);

  /** \brief Get a const Tpetra::MultiVector from a const
   * Thyra::MultiVectorBase object that is a Thyra::ProductMultiVectorBase.
   */
  static RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getConstTpetraMultiVectorBlock(const RCP<const MultiVectorBase<Scalar> > &mv, const int k);

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


namespace {

template <typename Source, typename Target>
struct match_const {
  using type = std::remove_const_t<Target>;
};

template <typename Source, typename Target>
struct match_const<const Source, Target> {
  using type = const Target;
};

template <typename Source, typename Target>
using match_const_t = typename match_const<Source, Target>::type;

} // namespace


/** \brief Get a const Tpetra::Map from a const Thyra::VectorSpaceBase object.
 */
template <class Scalar = Tpetra::Vector<>::scalar_type,
          class LocalOrdinal = Tpetra::Vector<>::local_ordinal_type,
          class GlobalOrdinal = Tpetra::Vector<>::global_ordinal_type,
          class Node = Tpetra::Vector<>::node_type>
RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
toTpetra(const RCP<const Thyra::VectorSpaceBase<Scalar>> &vs);


/** \brief Get a Tpetra::Vector or Tpetra::MultiVector from a
 * Thyra::VectorBase or Thyra::MultiVectorBase object.
 */
template <class Scalar = Tpetra::Vector<>::scalar_type,
          class LocalOrdinal = Tpetra::Vector<>::local_ordinal_type,
          class GlobalOrdinal = Tpetra::Vector<>::global_ordinal_type,
          class Node = Tpetra::Vector<>::node_type,
          class thy_vec_type,
          class tp_vec_type = std::conditional_t<
              std::is_same_v<std::remove_const_t<thy_vec_type>,
                             MultiVectorBase<Scalar>>,
              Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>,
              Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>,
          class ret_type = match_const_t<thy_vec_type, tp_vec_type>>
std::enable_if_t<
    std::is_same_v<std::remove_const_t<thy_vec_type>, MultiVectorBase<Scalar >> ||
    std::is_same_v<std::remove_const_t<thy_vec_type>, VectorBase<Scalar>>,
    RCP<ret_type>>
toTpetra(const RCP<thy_vec_type> &mv);


/** \brief Get a Tpetra::Vector or Tpetra::MultiVector from a
 * Thyra::ProductVectorBase or Thyra::ProductMultiVectorBase object.
 */
template <class Scalar = Tpetra::Vector<>::scalar_type,
          class LocalOrdinal = Tpetra::Vector<>::local_ordinal_type,
          class GlobalOrdinal = Tpetra::Vector<>::global_ordinal_type,
          class Node = Tpetra::Vector<>::node_type,
          class thy_vec_type,
          class tp_vec_type = std::conditional_t<
              std::is_same_v<std::remove_const_t<thy_vec_type>,
                             MultiVectorBase<Scalar>>,
              Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>,
              Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>,
          class ret_type = match_const_t<thy_vec_type, tp_vec_type>>
RCP<ret_type>
toTpetra(const RCP<thy_vec_type> &mv, const int k);


/** \brief Get a Tpetra::Operator from a Thyra::LinearOpBase object.
 */
template <class Scalar = Tpetra::Vector<>::scalar_type,
          class LocalOrdinal = Tpetra::Vector<>::local_ordinal_type,
          class GlobalOrdinal = Tpetra::Vector<>::global_ordinal_type,
          class Node = Tpetra::Vector<>::node_type,
          class thy_op_type,
          class ret_type = match_const_t<
              thy_op_type,
              Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>
std::enable_if_t<std::is_same_v<std::remove_const_t<thy_op_type>,
                                LinearOpBase<Scalar>>, RCP<ret_type>>
toTpetra(const RCP<thy_op_type> &op);

} // namespace Thyra


#endif // THYRA_TPETRA_THYRA_WRAPPERS_DECL_HPP

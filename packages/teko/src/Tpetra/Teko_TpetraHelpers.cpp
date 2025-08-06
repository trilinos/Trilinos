// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_TpetraHelpers.hpp"
#include "Teko_ConfigDefs.hpp"

// Thyra Includes
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"

#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_ScalarProdVectorSpaceBase.hpp"

// Teko includes
#include "Teko_Utilities.hpp"

// Tpetra
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraMultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "TpetraExt_MatrixMatrix.hpp"

using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcpFromRef;

namespace Teko {
namespace TpetraHelpers {

/** \brief Convert a Tpetra_Vector into a diagonal linear operator.
 *
 * Convert a Tpetra_Vector into a diagonal linear operator.
 *
 * \param[in] tv  Tpetra_Vector to use as the diagonal
 * \param[in] map Map related to the Tpetra_Vector
 * \param[in] lbl String to easily label the operator
 *
 * \returns A diagonal linear operator using the vector
 */
const Teuchos::RCP<const Thyra::LinearOpBase<ST> > thyraDiagOp(
    const RCP<const Tpetra::Vector<ST, LO, GO, NT> >& tv, const Tpetra::Map<LO, GO, NT>& map,
    const std::string& lbl) {
  const RCP<const Thyra::VectorBase<ST> > thyraVec  // need a Thyra::VectorBase object
      = Thyra::createConstVector<ST, LO, GO, NT>(
          tv, Thyra::createVectorSpace<ST, LO, GO, NT>(rcpFromRef(map)));
  Teuchos::RCP<Thyra::LinearOpBase<ST> > op =
      Teuchos::rcp(new Thyra::DefaultDiagonalLinearOp<ST>(thyraVec));
  op->setObjectLabel(lbl);
  return op;
}

/** \brief Convert a Tpetra_Vector into a diagonal linear operator.
 *
 * Convert a Tpetra_Vector into a diagonal linear operator.
 *
 * \param[in] tv  Tpetra_Vector to use as the diagonal
 * \param[in] map Map related to the Tpetra_Vector
 * \param[in] lbl String to easily label the operator
 *
 * \returns A diagonal linear operator using the vector
 */
const Teuchos::RCP<Thyra::LinearOpBase<ST> > thyraDiagOp(
    const RCP<Tpetra::Vector<ST, LO, GO, NT> >& tv, const Tpetra::Map<LO, GO, NT>& map,
    const std::string& lbl) {
  const RCP<Thyra::VectorBase<ST> > thyraVec  // need a Thyra::VectorBase object
      = Thyra::createVector<ST, LO, GO, NT>(
          tv, Thyra::createVectorSpace<ST, LO, GO, NT>(rcpFromRef(map)));
  Teuchos::RCP<Thyra::LinearOpBase<ST> > op =
      Teuchos::rcp(new Thyra::DefaultDiagonalLinearOp<ST>(thyraVec));
  op->setObjectLabel(lbl);
  return op;
}

/** \brief Fill a Thyra vector with the contents of a tpetra vector. This prevents the
 *
 * Fill a Thyra vector with the contents of a tpetra vector. This prevents the need
 * to reallocate memory using a create_MultiVector routine. It also allows an aritrary
 * Thyra vector to be filled.
 *
 * \param[in,out] spmdMV Multi-vector to be filled.
 * \param[in]     mv     Tpetra multi-vector to be used in filling the Thyra vector.
 */
void fillDefaultSpmdMultiVector(Teuchos::RCP<Thyra::TpetraMultiVector<ST, LO, GO, NT> >& spmdMV,
                                Teuchos::RCP<Tpetra::MultiVector<ST, LO, GO, NT> >& tpetraMV) {
  // first get desired range and domain
  // const RCP<const Thyra::SpmdVectorSpaceBase<ST> > range  = spmdMV->spmdSpace();
  const RCP<Thyra::TpetraVectorSpace<ST, LO, GO, NT> > range =
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraMV->getMap());
  const RCP<const Thyra::ScalarProdVectorSpaceBase<ST> > domain =
      rcp_dynamic_cast<const Thyra::ScalarProdVectorSpaceBase<ST> >(spmdMV->domain());

  TEUCHOS_ASSERT((size_t)domain->dim() == tpetraMV->getNumVectors());

  // New local view of raw data
  if (!tpetraMV->isConstantStride())
    TEUCHOS_TEST_FOR_EXCEPT(true);  // ToDo: Implement views of non-contiguous mult-vectors!

  // Build the MultiVector
  spmdMV->initialize(range, domain, tpetraMV);

  // make sure the Epetra_MultiVector doesn't disappear prematurely
  Teuchos::set_extra_data<RCP<Tpetra::MultiVector<ST, LO, GO, NT> > >(
      tpetraMV, "Tpetra::MultiVector", Teuchos::outArg(spmdMV));
}

/** \brief Build a vector of the dirchlet row indices.
 *
 * Build a vector of the dirchlet row indices. That is, record the global
 * index of any row that is all zeros except for $1$ on the diagonal.
 *
 * \param[in]     rowMap   Map specifying which global indices this process examines
 * \param[in]     mat      Matrix to be examined
 * \param[in,out] indices Output list of indices corresponding to dirchlet rows (GIDs).
 */
void identityRowIndices(const Tpetra::Map<LO, GO, NT>& rowMap,
                        const Tpetra::CrsMatrix<ST, LO, GO, NT>& mat, std::vector<GO>& outIndices) {
  // loop over elements owned by this processor
  for (size_t i = 0; i < rowMap.getLocalNumElements(); i++) {
    bool rowIsIdentity = true;
    GO rowGID          = rowMap.getGlobalElement(i);

    size_t numEntries = mat.getNumEntriesInGlobalRow(i);
    auto indices = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_global_inds_host_view_type(
        Kokkos::ViewAllocateWithoutInitializing("rowIndices"), numEntries);
    auto values = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_values_host_view_type(
        Kokkos::ViewAllocateWithoutInitializing("rowIndices"), numEntries);

    mat.getGlobalRowCopy(rowGID, indices, values, numEntries);

    // loop over the columns of this row
    for (size_t j = 0; j < numEntries; j++) {
      GO colGID = indices(j);

      // look at row entries
      if (colGID == rowGID)
        rowIsIdentity &= values(j) == 1.0;
      else
        rowIsIdentity &= values(j) == 0.0;

      // not a dirchlet row...quit
      if (not rowIsIdentity) break;
    }

    // save a row that is dirchlet
    if (rowIsIdentity) outIndices.push_back(rowGID);
  }
}

/** \brief Zero out the value of a vector on the specified
 *        set of global indices.
 *
 * Zero out the value of a vector on the specified set of global
 * indices. The indices here are assumed to belong to the calling
 * process (i.e. zeroIndices $\in$ mv.Map()).
 *
 * \param[in,out] mv           Vector whose entries will be zeroed
 * \param[in]     zeroIndices Indices local to this process that need to be zeroed
 */
void zeroMultiVectorRowIndices(Tpetra::MultiVector<ST, LO, GO, NT>& mv,
                               const std::vector<GO>& zeroIndices) {
  LO colCnt = mv.getNumVectors();
  std::vector<GO>::const_iterator itr;

  // loop over the indices to zero
  for (itr = zeroIndices.begin(); itr != zeroIndices.end(); ++itr) {
    // loop over columns
    for (int j = 0; j < colCnt; j++) mv.replaceGlobalValue(*itr, j, 0.0);
  }
}

/** \brief Constructor for a ZeroedOperator.
 *
 * Build a ZeroedOperator based on a particular Epetra_Operator and
 * a set of indices to zero out. These indices must be local to this
 * processor as specified by RowMap().
 *
 * \param[in] zeroIndices Set of indices to zero out (must be local).
 * \param[in] op           Underlying epetra operator to use.
 */
ZeroedOperator::ZeroedOperator(const std::vector<GO>& zeroIndices,
                               const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> >& op)
    : zeroIndices_(zeroIndices), tpetraOp_(op) {}

//! Perform a matrix-vector product with certain rows zeroed out
void ZeroedOperator::apply(const Tpetra::MultiVector<ST, LO, GO, NT>& X,
                           Tpetra::MultiVector<ST, LO, GO, NT>& Y, Teuchos::ETransp mode, ST alpha,
                           ST beta) const {
  /*
     Epetra_MultiVector temp(X);
     zeroMultiVectorRowIndices(temp,zeroIndices_);
     int result = epetraOp_->Apply(temp,Y);
  */

  tpetraOp_->apply(X, Y, mode, alpha, beta);

  // zero a few of the rows
  zeroMultiVectorRowIndices(Y, zeroIndices_);
}

bool isTpetraLinearOp(const LinearOp& op) {
  // See if the operator is a TpetraLinearOp
  RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tOp =
      rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(op);
  if (!tOp.is_null()) return true;

  // See if the operator is a wrapped TpetraLinearOp
  ST scalar               = 0.0;
  Thyra::EOpTransp transp = Thyra::NOTRANS;
  RCP<const Thyra::LinearOpBase<ST> > wrapped_op;
  Thyra::unwrap(op, &scalar, &transp, &wrapped_op);
  tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(wrapped_op);
  if (!tOp.is_null()) return true;

  return false;
}

RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > getTpetraCrsMatrix(const LinearOp& op, ST* scalar,
                                                                 bool* transp) {
  // If the operator is a TpetraLinearOp
  RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tOp =
      rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(op);
  if (!tOp.is_null()) {
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > matrix =
        rcp_dynamic_cast<const Tpetra::CrsMatrix<ST, LO, GO, NT> >(tOp->getConstTpetraOperator(),
                                                                   true);
    *scalar = 1.0;
    *transp = false;
    return matrix;
  }

  // If the operator is a wrapped TpetraLinearOp
  RCP<const Thyra::LinearOpBase<ST> > wrapped_op;
  Thyra::EOpTransp eTransp = Thyra::NOTRANS;
  Thyra::unwrap(op, scalar, &eTransp, &wrapped_op);
  tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(wrapped_op, true);
  if (!tOp.is_null()) {
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > matrix =
        rcp_dynamic_cast<const Tpetra::CrsMatrix<ST, LO, GO, NT> >(tOp->getConstTpetraOperator(),
                                                                   true);
    *transp = true;
    if (eTransp == Thyra::NOTRANS) *transp = false;
    return matrix;
  }

  return Teuchos::null;
}

}  // end namespace TpetraHelpers
}  // end namespace Teko

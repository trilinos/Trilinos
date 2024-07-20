// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_Config.h"
#include "Teko_Utilities.hpp"

// Thyra includes
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_ZeroLinearOpBase.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include <utility>

#ifdef TEKO_HAVE_EPETRA
#include "Thyra_EpetraExtDiagScaledMatProdTransformer.hpp"
#include "Thyra_EpetraExtDiagScalingTransformer.hpp"
#include "Thyra_EpetraExtAddTransformer.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraOperatorWrapper.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#endif

// Teuchos includes
#include "Teuchos_Array.hpp"

// Epetra includes
#ifdef TEKO_HAVE_EPETRA
#include "Epetra_Operator.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "EpetraExt_Transpose_RowMatrix.h"
#include "EpetraExt_MatrixMatrix.h"
#include <EpetraExt_BlockMapOut.h>
#include <EpetraExt_RowMatrixOut.h>

#include "Teko_EpetraHelpers.hpp"
#include "Teko_EpetraOperatorWrapper.hpp"
#endif

// Anasazi includes
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziThyraAdapter.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBlockKrylovSchur.hpp"
#include "AnasaziStatusTestMaxIters.hpp"

// Isorropia includes
#ifdef Teko_ENABLE_Isorropia
#include "Isorropia_EpetraProber.hpp"
#endif

// Teko includes
#include "Teko_TpetraHelpers.hpp"
#include "Teko_TpetraOperatorWrapper.hpp"

// Tpetra
#include "Thyra_TpetraLinearOp.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "TpetraExt_MatrixMatrix.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include "MatrixMarket_Tpetra.hpp"

#include <cmath>

namespace Teko {

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;
#ifdef Teko_ENABLE_Isorropia
using Isorropia::Epetra::Prober;
#endif

const Teuchos::RCP<Teuchos::FancyOStream> getOutputStream() {
  Teuchos::RCP<Teuchos::FancyOStream> os = Teuchos::VerboseObjectBase::getDefaultOStream();

  // os->setShowProcRank(true);
  // os->setOutputToRootOnly(-1);
  return os;
}

// distance function...not parallel...entirely internal to this cpp file
inline double dist(int dim, double *coords, int row, int col) {
  double value = 0.0;
  for (int i = 0; i < dim; i++)
    value += std::pow(coords[dim * row + i] - coords[dim * col + i], 2.0);

  // the distance between the two
  return std::sqrt(value);
}

// distance function...not parallel...entirely internal to this cpp file
inline double dist(double *x, double *y, double *z, int stride, int row, int col) {
  double value = 0.0;
  if (x != 0) value += std::pow(x[stride * row] - x[stride * col], 2.0);
  if (y != 0) value += std::pow(y[stride * row] - y[stride * col], 2.0);
  if (z != 0) value += std::pow(z[stride * row] - z[stride * col], 2.0);

  // the distance between the two
  return std::sqrt(value);
}

/** \brief Build a graph Laplacian stenciled on a Epetra_CrsMatrix.
 *
 * This function builds a graph Laplacian given a (locally complete)
 * vector of coordinates and a stencil Epetra_CrsMatrix (could this be
 * a graph of Epetra_RowMatrix instead?). The resulting matrix will have
 * the negative of the inverse distance on off diagonals. And the sum
 * of the positive inverse distance of the off diagonals on the diagonal.
 * If there are no off diagonal entries in the stencil, the diagonal is
 * set to 0.
 *
 * \param[in]     dim     Number of physical dimensions (2D or 3D?).
 * \param[in]     coords  A vector containing the coordinates, with the <code>i</code>-th
 *                        coordinate beginning at <code>coords[i*dim]</code>.
 * \param[in]     stencil The stencil matrix used to describe the connectivity
 *                        of the graph Laplacian matrix.
 *
 * \returns The graph Laplacian matrix to be filled according to the <code>stencil</code> matrix.
 */
#ifdef TEKO_HAVE_EPETRA
RCP<Epetra_CrsMatrix> buildGraphLaplacian(int dim, double *coords,
                                          const Epetra_CrsMatrix &stencil) {
  // allocate a new matrix with storage for the laplacian...in case of diagonals add one extra
  // storage
  RCP<Epetra_CrsMatrix> gl = rcp(
      new Epetra_CrsMatrix(Copy, stencil.RowMap(), stencil.ColMap(), stencil.MaxNumEntries() + 1),
      true);

  // allocate an additional value for the diagonal, if neccessary
  std::vector<double> rowData(stencil.GlobalMaxNumEntries() + 1);
  std::vector<int> rowInd(stencil.GlobalMaxNumEntries() + 1);

  // loop over all the rows
  for (int j = 0; j < gl->NumMyRows(); j++) {
    int row          = gl->GRID(j);
    double diagValue = 0.0;
    int diagInd      = -1;
    int rowSz        = 0;

    // extract a copy of this row...put it in rowData, rowIndicies
    stencil.ExtractGlobalRowCopy(row, stencil.MaxNumEntries(), rowSz, &rowData[0], &rowInd[0]);

    // loop over elements of row
    for (int i = 0; i < rowSz; i++) {
      int col = rowInd[i];

      // is this a 0 entry masquerading as some thing else?
      double value = rowData[i];
      if (value == 0) continue;

      // for nondiagonal entries
      if (row != col) {
        double d   = dist(dim, coords, row, col);
        rowData[i] = -1.0 / d;
        diagValue += rowData[i];
      } else
        diagInd = i;
    }

    // handle diagonal entry
    if (diagInd < 0) {  // diagonal not in row
      rowData[rowSz] = -diagValue;
      rowInd[rowSz]  = row;
      rowSz++;
    } else {  // diagonal in row
      rowData[diagInd] = -diagValue;
      rowInd[diagInd]  = row;
    }

    // insert row data into graph Laplacian matrix
    TEUCHOS_TEST_FOR_EXCEPT(gl->InsertGlobalValues(row, rowSz, &rowData[0], &rowInd[0]));
  }

  gl->FillComplete();

  return gl;
}
#endif

RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> buildGraphLaplacian(
    int dim, ST *coords, const Tpetra::CrsMatrix<ST, LO, GO, NT> &stencil) {
  // allocate a new matrix with storage for the laplacian...in case of diagonals add one extra
  // storage
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> gl = rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(
      stencil.getRowMap(), stencil.getColMap(), stencil.getGlobalMaxNumRowEntries() + 1));

  // allocate an additional value for the diagonal, if neccessary
  auto rowInd = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_global_inds_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"),
      stencil.getGlobalMaxNumRowEntries() + 1);
  auto rowData = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_values_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"),
      stencil.getGlobalMaxNumRowEntries() + 1);

  // loop over all the rows
  for (LO j = 0; j < (LO)gl->getLocalNumRows(); j++) {
    GO row       = gl->getRowMap()->getGlobalElement(j);
    ST diagValue = 0.0;
    GO diagInd   = -1;
    size_t rowSz = 0;

    // extract a copy of this row...put it in rowData, rowIndicies
    stencil.getGlobalRowCopy(row, rowInd, rowData, rowSz);

    // loop over elements of row
    for (size_t i = 0; i < rowSz; i++) {
      GO col = rowInd(i);

      // is this a 0 entry masquerading as some thing else?
      ST value = rowData[i];
      if (value == 0) continue;

      // for nondiagonal entries
      if (row != col) {
        ST d       = dist(dim, coords, row, col);
        rowData[i] = -1.0 / d;
        diagValue += rowData(i);
      } else
        diagInd = i;
    }

    // handle diagonal entry
    if (diagInd < 0) {  // diagonal not in row
      rowData(rowSz) = -diagValue;
      rowInd(rowSz)  = row;
      rowSz++;
    } else {  // diagonal in row
      rowData(diagInd) = -diagValue;
      rowInd(diagInd)  = row;
    }

    // insert row data into graph Laplacian matrix
    gl->replaceGlobalValues(row, rowInd, rowData);
  }

  gl->fillComplete();

  return gl;
}

/** \brief Build a graph Laplacian stenciled on a Epetra_CrsMatrix.
 *
 * This function builds a graph Laplacian given a (locally complete)
 * vector of coordinates and a stencil Epetra_CrsMatrix (could this be
 * a graph of Epetra_RowMatrix instead?). The resulting matrix will have
 * the negative of the inverse distance on off diagonals. And the sum
 * of the positive inverse distance of the off diagonals on the diagonal.
 * If there are no off diagonal entries in the stencil, the diagonal is
 * set to 0.
 *
 * \param[in]     x       A vector containing the x-coordinates, with the <code>i</code>-th
 *                        coordinate beginning at <code>coords[i*stride]</code>.
 * \param[in]     y       A vector containing the y-coordinates, with the <code>i</code>-th
 *                        coordinate beginning at <code>coords[i*stride]</code>.
 * \param[in]     z       A vector containing the z-coordinates, with the <code>i</code>-th
 *                        coordinate beginning at <code>coords[i*stride]</code>.
 * \param[in]     stride  Stride between entries in the (x,y,z) coordinate array
 * \param[in]     stencil The stencil matrix used to describe the connectivity
 *                        of the graph Laplacian matrix.
 *
 * \returns The graph Laplacian matrix to be filled according to the <code>stencil</code> matrix.
 */
#ifdef TEKO_HAVE_EPETRA
RCP<Epetra_CrsMatrix> buildGraphLaplacian(double *x, double *y, double *z, int stride,
                                          const Epetra_CrsMatrix &stencil) {
  // allocate a new matrix with storage for the laplacian...in case of diagonals add one extra
  // storage
  RCP<Epetra_CrsMatrix> gl = rcp(
      new Epetra_CrsMatrix(Copy, stencil.RowMap(), stencil.ColMap(), stencil.MaxNumEntries() + 1),
      true);

  // allocate an additional value for the diagonal, if neccessary
  std::vector<double> rowData(stencil.GlobalMaxNumEntries() + 1);
  std::vector<int> rowInd(stencil.GlobalMaxNumEntries() + 1);

  // loop over all the rows
  for (int j = 0; j < gl->NumMyRows(); j++) {
    int row          = gl->GRID(j);
    double diagValue = 0.0;
    int diagInd      = -1;
    int rowSz        = 0;

    // extract a copy of this row...put it in rowData, rowIndicies
    stencil.ExtractGlobalRowCopy(row, stencil.MaxNumEntries(), rowSz, &rowData[0], &rowInd[0]);

    // loop over elements of row
    for (int i = 0; i < rowSz; i++) {
      int col = rowInd[i];

      // is this a 0 entry masquerading as some thing else?
      double value = rowData[i];
      if (value == 0) continue;

      // for nondiagonal entries
      if (row != col) {
        double d   = dist(x, y, z, stride, row, col);
        rowData[i] = -1.0 / d;
        diagValue += rowData[i];
      } else
        diagInd = i;
    }

    // handle diagonal entry
    if (diagInd < 0) {  // diagonal not in row
      rowData[rowSz] = -diagValue;
      rowInd[rowSz]  = row;
      rowSz++;
    } else {  // diagonal in row
      rowData[diagInd] = -diagValue;
      rowInd[diagInd]  = row;
    }

    // insert row data into graph Laplacian matrix
    TEUCHOS_TEST_FOR_EXCEPT(gl->InsertGlobalValues(row, rowSz, &rowData[0], &rowInd[0]));
  }

  gl->FillComplete();

  return gl;
}
#endif

RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> buildGraphLaplacian(
    ST *x, ST *y, ST *z, GO stride, const Tpetra::CrsMatrix<ST, LO, GO, NT> &stencil) {
  // allocate a new matrix with storage for the laplacian...in case of diagonals add one extra
  // storage
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> gl =
      rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(stencil.getRowMap(), stencil.getColMap(),
                                                stencil.getGlobalMaxNumRowEntries() + 1),
          true);

  // allocate an additional value for the diagonal, if neccessary
  auto rowInd = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_global_inds_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"),
      stencil.getGlobalMaxNumRowEntries() + 1);
  auto rowData = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_values_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"),
      stencil.getGlobalMaxNumRowEntries() + 1);

  // loop over all the rows
  for (LO j = 0; j < (LO)gl->getLocalNumRows(); j++) {
    GO row       = gl->getRowMap()->getGlobalElement(j);
    ST diagValue = 0.0;
    GO diagInd   = -1;
    size_t rowSz = 0;

    // extract a copy of this row...put it in rowData, rowIndicies
    stencil.getGlobalRowCopy(row, rowInd, rowData, rowSz);

    // loop over elements of row
    for (size_t i = 0; i < rowSz; i++) {
      GO col = rowInd(i);

      // is this a 0 entry masquerading as some thing else?
      ST value = rowData[i];
      if (value == 0) continue;

      // for nondiagonal entries
      if (row != col) {
        ST d       = dist(x, y, z, stride, row, col);
        rowData[i] = -1.0 / d;
        diagValue += rowData(i);
      } else
        diagInd = i;
    }

    // handle diagonal entry
    if (diagInd < 0) {  // diagonal not in row
      rowData(rowSz) = -diagValue;
      rowInd(rowSz)  = row;
      rowSz++;
    } else {  // diagonal in row
      rowData(diagInd) = -diagValue;
      rowInd(diagInd)  = row;
    }

    // insert row data into graph Laplacian matrix
    gl->replaceGlobalValues(row, rowInd, rowData);
  }

  gl->fillComplete();

  return gl;
}

/** \brief Apply a linear operator to a multivector (think of this as a matrix
 *        vector multiply).
 *
 * Apply a linear operator to a multivector. This also permits arbitrary scaling
 * and addition of the result. This function gives
 *
 *    \f$ y = \alpha A x + \beta y \f$
 *
 * \param[in]     A
 * \param[in]     x
 * \param[in,out] y
 * \param[in]     \alpha
 * \param[in]     \beta
 *
 */
void applyOp(const LinearOp &A, const MultiVector &x, MultiVector &y, double alpha, double beta) {
  Thyra::apply(*A, Thyra::NOTRANS, *x, y.ptr(), alpha, beta);
}

/** \brief Apply a transposed linear operator to a multivector (think of this as a matrix
 *        vector multiply).
 *
 * Apply a transposed linear operator to a multivector. This also permits arbitrary scaling
 * and addition of the result. This function gives
 *
 *    \f$ y = \alpha A^T x + \beta y \f$
 *
 * \param[in]     A
 * \param[in]     x
 * \param[in,out] y
 * \param[in]     \alpha
 * \param[in]     \beta
 *
 */
void applyTransposeOp(const LinearOp &A, const MultiVector &x, MultiVector &y, double alpha,
                      double beta) {
  Thyra::apply(*A, Thyra::TRANS, *x, y.ptr(), alpha, beta);
}

/** \brief Update the <code>y</code> vector so that \f$y = \alpha x+\beta y\f$
 */
void update(double alpha, const MultiVector &x, double beta, MultiVector &y) {
  Teuchos::Array<double> scale;
  Teuchos::Array<Teuchos::Ptr<const Thyra::MultiVectorBase<double>>> vec;

  // build arrays needed for linear combo
  scale.push_back(alpha);
  vec.push_back(x.ptr());

  // compute linear combination
  Thyra::linear_combination<double>(scale, vec, beta, y.ptr());
}

//! Get the strictly upper triangular portion of the matrix
BlockedLinearOp getUpperTriBlocks(const BlockedLinearOp &blo, bool callEndBlockFill) {
  int rows = blockRowCount(blo);

  TEUCHOS_ASSERT(rows == blockColCount(blo));

  RCP<const Thyra::ProductVectorSpaceBase<double>> range  = blo->productRange();
  RCP<const Thyra::ProductVectorSpaceBase<double>> domain = blo->productDomain();

  // allocate new operator
  BlockedLinearOp upper = createBlockedOp();

  // build new operator
  upper->beginBlockFill(rows, rows);

  for (int i = 0; i < rows; i++) {
    // put zero operators on the diagonal
    // this gurantees the vector space of
    // the new operator are fully defined
    RCP<const Thyra::LinearOpBase<double>> zed =
        Thyra::zero<double>(range->getBlock(i), domain->getBlock(i));
    upper->setBlock(i, i, zed);

    for (int j = i + 1; j < rows; j++) {
      // get block i,j
      LinearOp uij = blo->getBlock(i, j);

      // stuff it in U
      if (uij != Teuchos::null) upper->setBlock(i, j, uij);
    }
  }
  if (callEndBlockFill) upper->endBlockFill();

  return upper;
}

//! Get the strictly lower triangular portion of the matrix
BlockedLinearOp getLowerTriBlocks(const BlockedLinearOp &blo, bool callEndBlockFill) {
  int rows = blockRowCount(blo);

  TEUCHOS_ASSERT(rows == blockColCount(blo));

  RCP<const Thyra::ProductVectorSpaceBase<double>> range  = blo->productRange();
  RCP<const Thyra::ProductVectorSpaceBase<double>> domain = blo->productDomain();

  // allocate new operator
  BlockedLinearOp lower = createBlockedOp();

  // build new operator
  lower->beginBlockFill(rows, rows);

  for (int i = 0; i < rows; i++) {
    // put zero operators on the diagonal
    // this gurantees the vector space of
    // the new operator are fully defined
    RCP<const Thyra::LinearOpBase<double>> zed =
        Thyra::zero<double>(range->getBlock(i), domain->getBlock(i));
    lower->setBlock(i, i, zed);

    for (int j = 0; j < i; j++) {
      // get block i,j
      LinearOp uij = blo->getBlock(i, j);

      // stuff it in U
      if (uij != Teuchos::null) lower->setBlock(i, j, uij);
    }
  }
  if (callEndBlockFill) lower->endBlockFill();

  return lower;
}

/** \brief Build a zero operator mimicing the block structure
 *        of the passed in matrix.
 *
 * Build a zero operator mimicing the block structure
 * of the passed in matrix. Currently this function assumes
 * that the operator is "block" square. Also, this function
 * calls <code>beginBlockFill</code> but does not call
 * <code>endBlockFill</code>.  This is so that the user
 * can fill the matrix as they wish once created.
 *
 * \param[in] blo Blocked operator with desired structure.
 *
 * \returns A zero operator with the same block structure as
 *          the argument <code>blo</code>.
 *
 * \notes The caller is responsible for calling
 *        <code>endBlockFill</code> on the returned blocked
 *        operator.
 */
BlockedLinearOp zeroBlockedOp(const BlockedLinearOp &blo) {
  int rows = blockRowCount(blo);

  TEUCHOS_ASSERT(rows == blockColCount(blo));  // assert that matrix is square

  RCP<const Thyra::ProductVectorSpaceBase<double>> range  = blo->productRange();
  RCP<const Thyra::ProductVectorSpaceBase<double>> domain = blo->productDomain();

  // allocate new operator
  BlockedLinearOp zeroOp = createBlockedOp();

  // build new operator
  zeroOp->beginBlockFill(rows, rows);

  for (int i = 0; i < rows; i++) {
    // put zero operators on the diagonal
    // this gurantees the vector space of
    // the new operator are fully defined
    RCP<const Thyra::LinearOpBase<double>> zed =
        Thyra::zero<double>(range->getBlock(i), domain->getBlock(i));
    zeroOp->setBlock(i, i, zed);
  }

  return zeroOp;
}

//! Figure out if this operator is the zero operator (or null!)
bool isZeroOp(const LinearOp op) {
  // if operator is null...then its zero!
  if (op == Teuchos::null) return true;

  // try to cast it to a zero linear operator
  LinearOp test = rcp_dynamic_cast<const Thyra::ZeroLinearOpBase<double>>(op);

  // if it works...then its zero...otherwise its null
  if (test != Teuchos::null) return true;

  // See if the operator is a wrapped zero op
  ST scalar               = 0.0;
  Thyra::EOpTransp transp = Thyra::NOTRANS;
  RCP<const Thyra::LinearOpBase<ST>> wrapped_op;
  Thyra::unwrap(op, &scalar, &transp, &wrapped_op);
  test = rcp_dynamic_cast<const Thyra::ZeroLinearOpBase<double>>(wrapped_op);
  return test != Teuchos::null;
}

std::pair<ModifiableLinearOp, bool> getAbsRowSumMatrixEpetra(const LinearOp &op) {
#ifndef TEKO_HAVE_EPETRA
  return std::make_pair(ModifiableLinearOp{}, false);
#else
  RCP<const Epetra_CrsMatrix> eCrsOp;

  const auto eOp = rcp_dynamic_cast<const Thyra::EpetraLinearOp>(op);

  if (!eOp) {
    return std::make_pair(ModifiableLinearOp{}, false);
  }

  eCrsOp = rcp_dynamic_cast<const Epetra_CrsMatrix>(eOp->epetra_op(), true);

  // extract diagonal
  const auto ptrDiag  = rcp(new Epetra_Vector(eCrsOp->RowMap()));
  Epetra_Vector &diag = *ptrDiag;

  // compute absolute value row sum
  diag.PutScalar(0.0);
  for (int i = 0; i < eCrsOp->NumMyRows(); i++) {
    double *values = 0;
    int numEntries;
    eCrsOp->ExtractMyRowView(i, numEntries, values);

    // build abs value row sum
    for (int j = 0; j < numEntries; j++) diag[i] += std::abs(values[j]);
  }

  // build Thyra diagonal operator
  return std::make_pair(Teko::Epetra::thyraDiagOp(ptrDiag, eCrsOp->RowMap(),
                                                  "absRowSum( " + op->getObjectLabel() + " )"),
                        true);
#endif
}

std::pair<ModifiableLinearOp, bool> getAbsRowSumMatrixTpetra(const LinearOp &op) {
  RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOp;

  const auto tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT>>(op);

  tCrsOp = rcp_dynamic_cast<const Tpetra::CrsMatrix<ST, LO, GO, NT>>(tOp->getConstTpetraOperator(),
                                                                     true);

  // extract diagonal
  const auto ptrDiag = Tpetra::createVector<ST, LO, GO, NT>(tCrsOp->getRowMap());
  auto &diag         = *ptrDiag;

  // compute absolute value row sum
  diag.putScalar(0.0);
  for (LO i = 0; i < (LO)tCrsOp->getLocalNumRows(); i++) {
    auto numEntries = tCrsOp->getNumEntriesInLocalRow(i);
    typename Tpetra::CrsMatrix<ST, LO, GO, NT>::local_inds_host_view_type indices;
    typename Tpetra::CrsMatrix<ST, LO, GO, NT>::values_host_view_type values;
    tCrsOp->getLocalRowView(i, indices, values);

    // build abs value row sum
    for (size_t j = 0; j < numEntries; j++) diag.sumIntoLocalValue(i, std::abs(values(j)));
  }

  // build Thyra diagonal operator
  return std::make_pair(
      Teko::TpetraHelpers::thyraDiagOp(ptrDiag, *tCrsOp->getRowMap(),
                                       "absRowSum( " + op->getObjectLabel() + " ))"),
      true);
}

/** \brief Compute absolute row sum matrix.
 *
 * Compute the absolute row sum matrix. That is
 * a diagonal operator composed of the absolute value of the
 * row sum.
 *
 * \returns A diagonal operator.
 */
ModifiableLinearOp getAbsRowSumMatrix(const LinearOp &op) {
  try {
    auto eResult = getAbsRowSumMatrixEpetra(op);
    if (eResult.second) {
      return eResult.first;
    }

    auto tResult = getAbsRowSumMatrixTpetra(op);
    if (tResult.second) {
      return tResult.first;
    } else {
      throw std::logic_error("Neither Epetra nor Tpetra");
    }
  } catch (std::exception &e) {
    auto out = Teuchos::VerboseObjectBase::getDefaultOStream();

    *out << "Teko: getAbsRowSumMatrix requires an Epetra_CrsMatrix or a "
            "Tpetra::CrsMatrix\n";
    *out << "    Could not extract an Epetra_Operator or a Tpetra_Operator "
            "from a \""
         << op->description() << std::endl;
    *out << "           OR\n";
    *out << "    Could not cast an Epetra_Operator to a Epetra_CrsMatrix or "
            "a Tpetra_Operator to a Tpetra::CrsMatrix\n";
    *out << std::endl;
    *out << "*** THROWN EXCEPTION ***\n";
    *out << e.what() << std::endl;
    *out << "************************\n";

    throw e;
  }
}

/** \brief Compute inverse of the absolute row sum matrix.
 *
 * Compute the inverse of the absolute row sum matrix. That is
 * a diagonal operator composed of the inverse of the absolute value
 * of the row sum.
 *
 * \returns A diagonal operator.
 */
ModifiableLinearOp getAbsRowSumInvMatrix(const LinearOp &op) {
  // if this is a blocked operator, extract diagonals block by block
  // FIXME: this does not add in values from off-diagonal blocks
  RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_op =
      rcp_dynamic_cast<const Thyra::PhysicallyBlockedLinearOpBase<double>>(op);
  if (blocked_op != Teuchos::null) {
    int numRows = blocked_op->productRange()->numBlocks();
    TEUCHOS_ASSERT(blocked_op->productDomain()->numBlocks() == numRows);
    RCP<Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_diag =
        Thyra::defaultBlockedLinearOp<double>();
    blocked_diag->beginBlockFill(numRows, numRows);
    for (int r = 0; r < numRows; ++r) {
      for (int c = 0; c < numRows; ++c) {
        if (r == c)
          blocked_diag->setNonconstBlock(r, c, getAbsRowSumInvMatrix(blocked_op->getBlock(r, c)));
        else
          blocked_diag->setBlock(r, c,
                                 Thyra::zero<double>(blocked_op->getBlock(r, c)->range(),
                                                     blocked_op->getBlock(r, c)->domain()));
      }
    }
    blocked_diag->endBlockFill();
    return blocked_diag;
  }

  if (Teko::TpetraHelpers::isTpetraLinearOp(op)) {
    ST scalar   = 0.0;
    bool transp = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOp =
        Teko::TpetraHelpers::getTpetraCrsMatrix(op, &scalar, &transp);

    // extract diagonal
    const RCP<Tpetra::Vector<ST, LO, GO, NT>> ptrDiag =
        Tpetra::createVector<ST, LO, GO, NT>(tCrsOp->getRowMap());
    Tpetra::Vector<ST, LO, GO, NT> &diag = *ptrDiag;

    // compute absolute value row sum
    diag.putScalar(0.0);
    for (LO i = 0; i < (LO)tCrsOp->getLocalNumRows(); i++) {
      LO numEntries = tCrsOp->getNumEntriesInLocalRow(i);
      typename Tpetra::CrsMatrix<ST, LO, GO, NT>::local_inds_host_view_type indices;
      typename Tpetra::CrsMatrix<ST, LO, GO, NT>::values_host_view_type values;
      tCrsOp->getLocalRowView(i, indices, values);

      // build abs value row sum
      for (LO j = 0; j < numEntries; j++) diag.sumIntoLocalValue(i, std::abs(values(j)));
    }
    diag.scale(scalar);
    diag.reciprocal(diag);  // invert entries

    // build Thyra diagonal operator
    return Teko::TpetraHelpers::thyraDiagOp(ptrDiag, *tCrsOp->getRowMap(),
                                            "absRowSum( " + op->getObjectLabel() + " ))");

  } else {
#ifdef TEKO_HAVE_EPETRA
    RCP<const Thyra::EpetraLinearOp> eOp = rcp_dynamic_cast<const Thyra::EpetraLinearOp>(op, true);
    RCP<const Epetra_CrsMatrix> eCrsOp =
        rcp_dynamic_cast<const Epetra_CrsMatrix>(eOp->epetra_op(), true);

    // extract diagonal
    const RCP<Epetra_Vector> ptrDiag = rcp(new Epetra_Vector(eCrsOp->RowMap()));
    Epetra_Vector &diag              = *ptrDiag;

    // compute absolute value row sum
    diag.PutScalar(0.0);
    for (int i = 0; i < eCrsOp->NumMyRows(); i++) {
      double *values = 0;
      int numEntries;
      eCrsOp->ExtractMyRowView(i, numEntries, values);

      // build abs value row sum
      for (int j = 0; j < numEntries; j++) diag[i] += std::abs(values[j]);
    }
    diag.Reciprocal(diag);  // invert entries

    // build Thyra diagonal operator
    return Teko::Epetra::thyraDiagOp(ptrDiag, eCrsOp->RowMap(),
                                     "absRowSum( " + op->getObjectLabel() + " )");
#else
    throw std::logic_error(
        "getAbsRowSumInvMatrix is trying to use Epetra "
        "code, but TEKO_HAVE_EPETRA is disabled!");
#endif
  }
}

/** \brief Compute the lumped version of this matrix.
 *
 * Compute the lumped version of this matrix. That is
 * a diagonal operator composed of the row sum.
 *
 * \returns A diagonal operator.
 */
ModifiableLinearOp getLumpedMatrix(const LinearOp &op) {
  RCP<Thyra::VectorBase<ST>> ones = Thyra::createMember(op->domain());
  RCP<Thyra::VectorBase<ST>> diag = Thyra::createMember(op->range());

  // set to all ones
  Thyra::assign(ones.ptr(), 1.0);

  // compute lumped diagonal
  // Thyra::apply(*op,Thyra::NONCONJ_ELE,*ones,&*diag);
  Thyra::apply(*op, Thyra::NOTRANS, *ones, diag.ptr());

  return rcp(new Thyra::DefaultDiagonalLinearOp<ST>(diag));
}

/** \brief Compute the inverse of the lumped version of
 *        this matrix.
 *
 * Compute the inverse of the lumped version of this matrix.
 * That is a diagonal operator composed of the row sum.
 *
 * \returns A diagonal operator.
 */
ModifiableLinearOp getInvLumpedMatrix(const LinearOp &op) {
  RCP<Thyra::VectorBase<ST>> ones = Thyra::createMember(op->domain());
  RCP<Thyra::VectorBase<ST>> diag = Thyra::createMember(op->range());

  // set to all ones
  Thyra::assign(ones.ptr(), 1.0);

  // compute lumped diagonal
  Thyra::apply(*op, Thyra::NOTRANS, *ones, diag.ptr());
  Thyra::reciprocal(*diag, diag.ptr());

  return rcp(new Thyra::DefaultDiagonalLinearOp<ST>(diag));
}

const std::pair<ModifiableLinearOp, bool> getDiagonalOpEpetra(const LinearOp &op) {
#ifndef TEKO_HAVE_EPETRA
  return std::make_pair(ModifiableLinearOp{}, false);
#else
  RCP<const Epetra_CrsMatrix> eCrsOp;

  const auto eOp = rcp_dynamic_cast<const Thyra::EpetraLinearOp>(op);
  if (!eOp) {
    return std::make_pair(ModifiableLinearOp{}, false);
  }

  eCrsOp = rcp_dynamic_cast<const Epetra_CrsMatrix>(eOp->epetra_op(), true);

  // extract diagonal
  const auto diag = rcp(new Epetra_Vector(eCrsOp->RowMap()));
  TEUCHOS_TEST_FOR_EXCEPT(eCrsOp->ExtractDiagonalCopy(*diag));

  // build Thyra diagonal operator
  return std::make_pair(Teko::Epetra::thyraDiagOp(diag, eCrsOp->RowMap(),
                                                  "inv(diag( " + op->getObjectLabel() + " ))"),
                        true);
#endif
}

const std::pair<ModifiableLinearOp, bool> getDiagonalOpTpetra(const LinearOp &op) {
  RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOp;

  const auto tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT>>(op);
  if (!tOp) {
    return std::make_pair(ModifiableLinearOp{}, false);
  }

  tCrsOp = rcp_dynamic_cast<const Tpetra::CrsMatrix<ST, LO, GO, NT>>(tOp->getConstTpetraOperator(),
                                                                     true);

  // extract diagonal
  const auto diag = Tpetra::createVector<ST, LO, GO, NT>(tCrsOp->getRowMap());
  tCrsOp->getLocalDiagCopy(*diag);

  // build Thyra diagonal operator
  return std::make_pair(
      Teko::TpetraHelpers::thyraDiagOp(diag, *tCrsOp->getRowMap(),
                                       "inv(diag( " + op->getObjectLabel() + " ))"),
      true);
}

/** \brief Get the diaonal of a linear operator
 *
 * Get the diagonal of a linear operator. Currently
 * it is assumed that the underlying operator is
 * an Epetra_RowMatrix.
 *
 * \param[in] op The operator whose diagonal is to be
 *               extracted.
 *
 * \returns An diagonal operator.
 */
const ModifiableLinearOp getDiagonalOp(const LinearOp &op) {
  try {
    // get Epetra or Tpetra Operator
    const auto eDiagOp = getDiagonalOpEpetra(op);

    if (eDiagOp.second) {
      return eDiagOp.first;
    }

    const auto tDiagOp = getDiagonalOpTpetra(op);
    if (tDiagOp.second) {
      return tDiagOp.first;
    } else
      throw std::logic_error("Neither Epetra nor Tpetra");
  } catch (std::exception &e) {
    RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

    *out << "Teko: getDiagonalOp requires an Epetra_CrsMatrix or a Tpetra::CrsMatrix\n";
    *out << "    Could not extract an Epetra_Operator or a Tpetra_Operator from a \""
         << op->description() << std::endl;
    *out << "           OR\n";
    *out << "    Could not cast an Epetra_Operator to a Epetra_CrsMatrix or a Tpetra_Operator to a "
            "Tpetra::CrsMatrix\n";
    *out << std::endl;
    *out << "*** THROWN EXCEPTION ***\n";
    *out << e.what() << std::endl;
    *out << "************************\n";

    throw e;
  }
}

const MultiVector getDiagonal(const LinearOp &op) {
  try {
    // get Epetra or Tpetra Operator
    auto diagOp = getDiagonalOpEpetra(op);

    if (!diagOp.second) {
      diagOp = getDiagonalOpTpetra(op);

      if (!diagOp.second) {
        throw std::logic_error("Neither Epetra nor Tpetra");
      }
    }

    Teuchos::RCP<const Thyra::MultiVectorBase<double>> v =
        Teuchos::rcp_dynamic_cast<const Thyra::DiagonalLinearOpBase<double>>(diagOp.first)
            ->getDiag();
    return Teuchos::rcp_const_cast<Thyra::MultiVectorBase<double>>(v);
  } catch (std::exception &e) {
    RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

    *out << "Teko: getDiagonal requires an Epetra_CrsMatrix or a Tpetra::CrsMatrix\n";
    *out << "    Could not extract an Epetra_Operator or a Tpetra_Operator from a \""
         << op->description() << std::endl;
    *out << "           OR\n";
    *out << "    Could not cast an Epetra_Operator to a Epetra_CrsMatrix or a Tpetra_Operator to a "
            "Tpetra::CrsMatrix\n";
    *out << std::endl;
    *out << "*** THROWN EXCEPTION ***\n";
    *out << e.what() << std::endl;
    *out << "************************\n";

    throw e;
  }
}

const MultiVector getDiagonal(const Teko::LinearOp &A, const DiagonalType &dt) {
  LinearOp diagOp = Teko::getDiagonalOp(A, dt);

  Teuchos::RCP<const Thyra::MultiVectorBase<double>> v =
      Teuchos::rcp_dynamic_cast<const Thyra::DiagonalLinearOpBase<double>>(diagOp)->getDiag();
  return Teuchos::rcp_const_cast<Thyra::MultiVectorBase<double>>(v);
}

/** \brief Get the diaonal of a linear operator
 *
 * Get the inverse of the diagonal of a linear operator.
 * Currently it is assumed that the underlying operator is
 * an Epetra_RowMatrix.
 *
 * \param[in] op The operator whose diagonal is to be
 *               extracted and inverted
 *
 * \returns An diagonal operator.
 */
const ModifiableLinearOp getInvDiagonalOp(const LinearOp &op) {
  // if this is a diagonal linear op already, just take the reciprocal
  auto diagonal_op = rcp_dynamic_cast<const Thyra::DiagonalLinearOpBase<double>>(op);
  if (diagonal_op != Teuchos::null) {
    auto diag     = diagonal_op->getDiag();
    auto inv_diag = diag->clone_v();
    Thyra::reciprocal(*diag, inv_diag.ptr());
    return rcp(new Thyra::DefaultDiagonalLinearOp<double>(inv_diag));
  }

  // if this is a blocked operator, extract diagonals block by block
  RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_op =
      rcp_dynamic_cast<const Thyra::PhysicallyBlockedLinearOpBase<double>>(op);
  if (blocked_op != Teuchos::null) {
    int numRows = blocked_op->productRange()->numBlocks();
    TEUCHOS_ASSERT(blocked_op->productDomain()->numBlocks() == numRows);
    RCP<Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_diag =
        Thyra::defaultBlockedLinearOp<double>();
    blocked_diag->beginBlockFill(numRows, numRows);
    for (int r = 0; r < numRows; ++r) {
      for (int c = 0; c < numRows; ++c) {
        if (r == c)
          blocked_diag->setNonconstBlock(r, c, getInvDiagonalOp(blocked_op->getBlock(r, c)));
        else
          blocked_diag->setBlock(r, c,
                                 Thyra::zero<double>(blocked_op->getBlock(r, c)->range(),
                                                     blocked_op->getBlock(r, c)->domain()));
      }
    }
    blocked_diag->endBlockFill();
    return blocked_diag;
  }

  if (Teko::TpetraHelpers::isTpetraLinearOp(op)) {
    ST scalar   = 0.0;
    bool transp = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOp =
        Teko::TpetraHelpers::getTpetraCrsMatrix(op, &scalar, &transp);

    // extract diagonal
    const RCP<Tpetra::Vector<ST, LO, GO, NT>> diag =
        Tpetra::createVector<ST, LO, GO, NT>(tCrsOp->getRowMap());
    diag->scale(scalar);
    tCrsOp->getLocalDiagCopy(*diag);
    diag->reciprocal(*diag);

    // build Thyra diagonal operator
    return Teko::TpetraHelpers::thyraDiagOp(diag, *tCrsOp->getRowMap(),
                                            "inv(diag( " + op->getObjectLabel() + " ))");

  } else {
#ifdef TEKO_HAVE_EPETRA
    RCP<const Thyra::EpetraLinearOp> eOp = rcp_dynamic_cast<const Thyra::EpetraLinearOp>(op, true);
    RCP<const Epetra_CrsMatrix> eCrsOp =
        rcp_dynamic_cast<const Epetra_CrsMatrix>(eOp->epetra_op(), true);

    // extract diagonal
    const RCP<Epetra_Vector> diag = rcp(new Epetra_Vector(eCrsOp->RowMap()));
    TEUCHOS_TEST_FOR_EXCEPT(eCrsOp->ExtractDiagonalCopy(*diag));
    diag->Reciprocal(*diag);

    // build Thyra diagonal operator
    return Teko::Epetra::thyraDiagOp(diag, eCrsOp->RowMap(),
                                     "inv(diag( " + op->getObjectLabel() + " ))");
#else
    throw std::logic_error(
        "getInvDiagonalOp is trying to use Epetra "
        "code, but TEKO_HAVE_EPETRA is disabled!");
#endif
  }
}

/** \brief Multiply three linear operators.
 *
 * Multiply three linear operators. This currently assumes
 * that the underlying implementation uses Epetra_CrsMatrix.
 * The exception is that opm is allowed to be an diagonal matrix.
 *
 * \param[in] opl Left operator (assumed to be a Epetra_CrsMatrix)
 * \param[in] opm Middle operator (assumed to be a diagonal matrix)
 * \param[in] opr Right operator (assumed to be a Epetra_CrsMatrix)
 *
 * \returns Matrix product with a Epetra_CrsMatrix implementation
 */
const LinearOp explicitMultiply(const LinearOp &opl, const LinearOp &opm, const LinearOp &opr) {
  // if this is a blocked operator, multiply block by block
  // it is possible that not every factor in the product is blocked and these situations are handled
  // separately

  bool isBlockedL = isPhysicallyBlockedLinearOp(opl);
  bool isBlockedM = isPhysicallyBlockedLinearOp(opm);
  bool isBlockedR = isPhysicallyBlockedLinearOp(opr);

  // all factors blocked
  if ((isBlockedL && isBlockedM && isBlockedR)) {
    double scalarl = 0.0;
    bool transpl   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opl =
        getPhysicallyBlockedLinearOp(opl, &scalarl, &transpl);
    double scalarm = 0.0;
    bool transpm   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opm =
        getPhysicallyBlockedLinearOp(opm, &scalarm, &transpm);
    double scalarr = 0.0;
    bool transpr   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opr =
        getPhysicallyBlockedLinearOp(opr, &scalarr, &transpr);
    double scalar = scalarl * scalarm * scalarr;

    int numRows   = blocked_opl->productRange()->numBlocks();
    int numCols   = blocked_opr->productDomain()->numBlocks();
    int numMiddle = blocked_opm->productRange()->numBlocks();

    // Assume that the middle block is block nxn and that it's diagonal. Otherwise use the two
    // argument explicitMultiply twice
    TEUCHOS_ASSERT(blocked_opm->productDomain()->numBlocks() == numMiddle);
    TEUCHOS_ASSERT(blocked_opl->productDomain()->numBlocks() == numMiddle);
    TEUCHOS_ASSERT(blocked_opr->productRange()->numBlocks() == numMiddle);

    RCP<Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_product =
        Thyra::defaultBlockedLinearOp<double>();
    blocked_product->beginBlockFill(numRows, numCols);
    for (int r = 0; r < numRows; ++r) {
      for (int c = 0; c < numCols; ++c) {
        LinearOp product_rc = explicitMultiply(
            blocked_opl->getBlock(r, 0), blocked_opm->getBlock(0, 0), blocked_opr->getBlock(0, c));
        for (int m = 1; m < numMiddle; ++m) {
          LinearOp product_m =
              explicitMultiply(blocked_opl->getBlock(r, m), blocked_opm->getBlock(m, m),
                               blocked_opr->getBlock(m, c));
          product_rc = explicitAdd(product_rc, product_m);
        }
        blocked_product->setBlock(r, c, product_rc);
      }
    }
    blocked_product->endBlockFill();
    return Thyra::scale<double>(scalar, blocked_product.getConst());
  }

  // left and right factors blocked
  if (isBlockedL && !isBlockedM && isBlockedR) {
    double scalarl = 0.0;
    bool transpl   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opl =
        getPhysicallyBlockedLinearOp(opl, &scalarl, &transpl);
    double scalarr = 0.0;
    bool transpr   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opr =
        getPhysicallyBlockedLinearOp(opr, &scalarr, &transpr);
    double scalar = scalarl * scalarr;

    int numRows   = blocked_opl->productRange()->numBlocks();
    int numCols   = blocked_opr->productDomain()->numBlocks();
    int numMiddle = 1;

    // Assume that the middle block is 1x1 diagonal. Left must be rx1, right 1xc
    TEUCHOS_ASSERT(blocked_opl->productDomain()->numBlocks() == numMiddle);
    TEUCHOS_ASSERT(blocked_opr->productRange()->numBlocks() == numMiddle);

    RCP<Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_product =
        Thyra::defaultBlockedLinearOp<double>();
    blocked_product->beginBlockFill(numRows, numCols);
    for (int r = 0; r < numRows; ++r) {
      for (int c = 0; c < numCols; ++c) {
        LinearOp product_rc =
            explicitMultiply(blocked_opl->getBlock(r, 0), opm, blocked_opr->getBlock(0, c));
        blocked_product->setBlock(r, c, product_rc);
      }
    }
    blocked_product->endBlockFill();
    return Thyra::scale<double>(scalar, blocked_product.getConst());
  }

  // only right factor blocked
  if (!isBlockedL && !isBlockedM && isBlockedR) {
    double scalarr = 0.0;
    bool transpr   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opr =
        getPhysicallyBlockedLinearOp(opr, &scalarr, &transpr);
    double scalar = scalarr;

    int numRows   = 1;
    int numCols   = blocked_opr->productDomain()->numBlocks();
    int numMiddle = 1;

    // Assume that the middle block is 1x1 diagonal, left is 1x1. Right must be 1xc
    TEUCHOS_ASSERT(blocked_opr->productRange()->numBlocks() == numMiddle);

    RCP<Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_product =
        Thyra::defaultBlockedLinearOp<double>();
    blocked_product->beginBlockFill(numRows, numCols);
    for (int c = 0; c < numCols; ++c) {
      LinearOp product_c = explicitMultiply(opl, opm, blocked_opr->getBlock(0, c));
      blocked_product->setBlock(0, c, product_c);
    }
    blocked_product->endBlockFill();
    return Thyra::scale<double>(scalar, blocked_product.getConst());
  }

  // TODO: three more cases (only non-blocked - blocked - non-blocked not possible)

  bool isTpetral = Teko::TpetraHelpers::isTpetraLinearOp(opl);
  bool isTpetram = Teko::TpetraHelpers::isTpetraLinearOp(opm);
  bool isTpetrar = Teko::TpetraHelpers::isTpetraLinearOp(opr);

  if (isTpetral && isTpetram &&
      isTpetrar) {  // Both operators are Tpetra matrices so explicitly multiply them

    // Get left and right Tpetra crs operators
    ST scalarl   = 0.0;
    bool transpl = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpl =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opl, &scalarl, &transpl);
    ST scalarm   = 0.0;
    bool transpm = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpm =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opl, &scalarm, &transpm);
    ST scalarr   = 0.0;
    bool transpr = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpr =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opr, &scalarr, &transpr);

    // Build output operator
    RCP<Thyra::LinearOpBase<ST>> explicitOp = rcp(new Thyra::TpetraLinearOp<ST, LO, GO, NT>());
    RCP<Thyra::TpetraLinearOp<ST, LO, GO, NT>> tExplicitOp =
        rcp_dynamic_cast<Thyra::TpetraLinearOp<ST, LO, GO, NT>>(explicitOp);

    // Do explicit matrix-matrix multiply
    RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOplm =
        Tpetra::createCrsMatrix<ST, LO, GO, NT>(tCrsOpl->getRowMap());
    RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> explicitCrsOp =
        Tpetra::createCrsMatrix<ST, LO, GO, NT>(tCrsOpl->getRowMap());
    Tpetra::MatrixMatrix::Multiply<ST, LO, GO, NT>(*tCrsOpl, transpl, *tCrsOpm, transpm, *tCrsOplm);
    Tpetra::MatrixMatrix::Multiply<ST, LO, GO, NT>(*tCrsOplm, false, *tCrsOpr, transpr,
                                                   *explicitCrsOp);
    explicitCrsOp->resumeFill();
    explicitCrsOp->scale(scalarl * scalarm * scalarr);
    explicitCrsOp->fillComplete(tCrsOpr->getDomainMap(), tCrsOpl->getRangeMap());
    tExplicitOp->initialize(Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getRangeMap()),
                            Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getDomainMap()),
                            explicitCrsOp);
    return tExplicitOp;

  } else if (isTpetral && !isTpetram && isTpetrar) {  // Assume that the middle operator is diagonal

    // Get left and right Tpetra crs operators
    ST scalarl   = 0.0;
    bool transpl = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpl =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opl, &scalarl, &transpl);
    ST scalarr   = 0.0;
    bool transpr = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpr =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opr, &scalarr, &transpr);

    RCP<const Tpetra::Vector<ST, LO, GO, NT>> diagPtr;

    // Cast middle operator as DiagonalLinearOp and extract diagonal as Vector
    RCP<const Thyra::DiagonalLinearOpBase<ST>> dOpm =
        rcp_dynamic_cast<const Thyra::DiagonalLinearOpBase<ST>>(opm);
    if (dOpm != Teuchos::null) {
      RCP<const Thyra::TpetraVector<ST, LO, GO, NT>> tPtr =
          rcp_dynamic_cast<const Thyra::TpetraVector<ST, LO, GO, NT>>(dOpm->getDiag(), true);
      diagPtr = rcp_dynamic_cast<const Tpetra::Vector<ST, LO, GO, NT>>(tPtr->getConstTpetraVector(),
                                                                       true);
    }
    // If it's not diagonal, maybe it's zero
    else if (rcp_dynamic_cast<const Thyra::ZeroLinearOpBase<ST>>(opm) != Teuchos::null) {
      diagPtr = rcp(new Tpetra::Vector<ST, LO, GO, NT>(tCrsOpl->getDomainMap()));
    } else
      TEUCHOS_ASSERT(false);

    RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOplm =
        Tpetra::importAndFillCompleteCrsMatrix<Tpetra::CrsMatrix<ST, LO, GO, NT>>(
            tCrsOpl, Tpetra::Import<LO, GO, NT>(tCrsOpl->getRowMap(), tCrsOpl->getRowMap()));

    // Do the diagonal scaling
    tCrsOplm->rightScale(*diagPtr);

    // Build output operator
    RCP<Thyra::LinearOpBase<ST>> explicitOp = rcp(new Thyra::TpetraLinearOp<ST, LO, GO, NT>());
    RCP<Thyra::TpetraLinearOp<ST, LO, GO, NT>> tExplicitOp =
        rcp_dynamic_cast<Thyra::TpetraLinearOp<ST, LO, GO, NT>>(explicitOp);

    // Do explicit matrix-matrix multiply
    RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> explicitCrsOp =
        Tpetra::createCrsMatrix<ST, LO, GO, NT>(tCrsOpl->getRowMap());
    Tpetra::MatrixMatrix::Multiply<ST, LO, GO, NT>(*tCrsOplm, false, *tCrsOpr, transpr,
                                                   *explicitCrsOp);
    explicitCrsOp->resumeFill();
    explicitCrsOp->scale(scalarl * scalarr);
    explicitCrsOp->fillComplete(tCrsOpr->getDomainMap(), tCrsOpl->getRangeMap());
    tExplicitOp->initialize(Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getRangeMap()),
                            Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getDomainMap()),
                            explicitCrsOp);
    return tExplicitOp;

  } else {  // Assume Epetra and we can use transformers
#ifdef TEKO_HAVE_EPETRA
    // build implicit multiply
    const LinearOp implicitOp = Thyra::multiply(opl, opm, opr);

    // build transformer
    const RCP<Thyra::LinearOpTransformerBase<double>> prodTrans =
        Thyra::epetraExtDiagScaledMatProdTransformer();

    // build operator and multiply
    const RCP<Thyra::LinearOpBase<double>> explicitOp = prodTrans->createOutputOp();
    prodTrans->transform(*implicitOp, explicitOp.ptr());
    explicitOp->setObjectLabel("explicit( " + opl->getObjectLabel() + " * " +
                               opm->getObjectLabel() + " * " + opr->getObjectLabel() + " )");

    return explicitOp;
#else
    throw std::logic_error(
        "explicitMultiply is trying to use Epetra "
        "code, but TEKO_HAVE_EPETRA is disabled!");
#endif
  }
}

/** \brief Multiply three linear operators.
 *
 * Multiply three linear operators. This currently assumes
 * that the underlying implementation uses Epetra_CrsMatrix.
 * The exception is that opm is allowed to be an diagonal matrix.
 *
 * \param[in] opl Left operator (assumed to be a Epetra_CrsMatrix)
 * \param[in] opm Middle operator (assumed to be a Epetra_CrsMatrix or a diagonal matrix)
 * \param[in] opr Right operator (assumed to be a Epetra_CrsMatrix)
 * \param[in,out] destOp The operator to be used as the destination operator,
 *                       if this is null this function creates a new operator
 *
 * \returns Matrix product with a Epetra_CrsMatrix implementation
 */
const ModifiableLinearOp explicitMultiply(const LinearOp &opl, const LinearOp &opm,
                                          const LinearOp &opr, const ModifiableLinearOp &destOp) {
  bool isTpetral = Teko::TpetraHelpers::isTpetraLinearOp(opl);
  bool isTpetram = Teko::TpetraHelpers::isTpetraLinearOp(opm);
  bool isTpetrar = Teko::TpetraHelpers::isTpetraLinearOp(opr);

  if (isTpetral && isTpetram &&
      isTpetrar) {  // Both operators are Tpetra matrices so explicitly multiply them

    // Get left and right Tpetra crs operators
    ST scalarl   = 0.0;
    bool transpl = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpl =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opl, &scalarl, &transpl);
    ST scalarm   = 0.0;
    bool transpm = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpm =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opm, &scalarm, &transpm);
    ST scalarr   = 0.0;
    bool transpr = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpr =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opr, &scalarr, &transpr);

    // Build output operator
    auto tExplicitOp = rcp_dynamic_cast<Thyra::TpetraLinearOp<ST, LO, GO, NT>>(destOp);
    if (tExplicitOp.is_null()) tExplicitOp = rcp(new Thyra::TpetraLinearOp<ST, LO, GO, NT>());

    // Do explicit matrix-matrix multiply
    RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOplm =
        Tpetra::createCrsMatrix<ST, LO, GO, NT>(tCrsOpl->getRowMap());
    RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> explicitCrsOp =
        Tpetra::createCrsMatrix<ST, LO, GO, NT>(tCrsOpl->getRowMap());
    Tpetra::MatrixMatrix::Multiply<ST, LO, GO, NT>(*tCrsOpl, transpl, *tCrsOpm, transpm, *tCrsOplm);
    Tpetra::MatrixMatrix::Multiply<ST, LO, GO, NT>(*tCrsOplm, false, *tCrsOpr, transpr,
                                                   *explicitCrsOp);
    explicitCrsOp->resumeFill();
    explicitCrsOp->scale(scalarl * scalarm * scalarr);
    explicitCrsOp->fillComplete(tCrsOpr->getDomainMap(), tCrsOpl->getRangeMap());
    tExplicitOp->initialize(Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getRangeMap()),
                            Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getDomainMap()),
                            explicitCrsOp);
    return tExplicitOp;

  } else if (isTpetral && !isTpetram && isTpetrar) {  // Assume that the middle operator is diagonal

    // Get left and right Tpetra crs operators
    ST scalarl   = 0.0;
    bool transpl = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpl =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opl, &scalarl, &transpl);
    ST scalarr   = 0.0;
    bool transpr = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpr =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opr, &scalarr, &transpr);

    // Cast middle operator as DiagonalLinearOp and extract diagonal as Vector
    RCP<const Thyra::DiagonalLinearOpBase<ST>> dOpm =
        rcp_dynamic_cast<const Thyra::DiagonalLinearOpBase<ST>>(opm, true);
    RCP<const Thyra::TpetraVector<ST, LO, GO, NT>> tPtr =
        rcp_dynamic_cast<const Thyra::TpetraVector<ST, LO, GO, NT>>(dOpm->getDiag(), true);
    RCP<const Tpetra::Vector<ST, LO, GO, NT>> diagPtr =
        rcp_dynamic_cast<const Tpetra::Vector<ST, LO, GO, NT>>(tPtr->getConstTpetraVector(), true);
    RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOplm =
        Tpetra::importAndFillCompleteCrsMatrix<Tpetra::CrsMatrix<ST, LO, GO, NT>>(
            tCrsOpl, Tpetra::Import<LO, GO, NT>(tCrsOpl->getRowMap(), tCrsOpl->getRowMap()));

    // Do the diagonal scaling
    tCrsOplm->rightScale(*diagPtr);

    // Build output operator
    RCP<Thyra::LinearOpBase<ST>> explicitOp = rcp(new Thyra::TpetraLinearOp<ST, LO, GO, NT>());
    RCP<Thyra::TpetraLinearOp<ST, LO, GO, NT>> tExplicitOp =
        rcp_dynamic_cast<Thyra::TpetraLinearOp<ST, LO, GO, NT>>(explicitOp);

    // Do explicit matrix-matrix multiply
    RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> explicitCrsOp =
        Tpetra::createCrsMatrix<ST, LO, GO, NT>(tCrsOpl->getRowMap());
    Tpetra::MatrixMatrix::Multiply<ST, LO, GO, NT>(*tCrsOplm, false, *tCrsOpr, transpr,
                                                   *explicitCrsOp);
    explicitCrsOp->resumeFill();
    explicitCrsOp->scale(scalarl * scalarr);
    explicitCrsOp->fillComplete(tCrsOpr->getDomainMap(), tCrsOpl->getRangeMap());
    tExplicitOp->initialize(Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getRangeMap()),
                            Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getDomainMap()),
                            explicitCrsOp);
    return tExplicitOp;

  } else {  // Assume Epetra and we can use transformers
#ifdef TEKO_HAVE_EPETRA
    // build implicit multiply
    const LinearOp implicitOp = Thyra::multiply(opl, opm, opr);

    // build transformer
    const RCP<Thyra::LinearOpTransformerBase<double>> prodTrans =
        Thyra::epetraExtDiagScaledMatProdTransformer();

    // build operator destination operator
    ModifiableLinearOp explicitOp;

    // if neccessary build a operator to put the explicit multiply into
    if (destOp == Teuchos::null)
      explicitOp = prodTrans->createOutputOp();
    else
      explicitOp = destOp;

    // perform multiplication
    prodTrans->transform(*implicitOp, explicitOp.ptr());

    // label it
    explicitOp->setObjectLabel("explicit( " + opl->getObjectLabel() + " * " +
                               opm->getObjectLabel() + " * " + opr->getObjectLabel() + " )");

    return explicitOp;
#else
    throw std::logic_error(
        "explicitMultiply is trying to use Epetra "
        "code, but TEKO_HAVE_EPETRA is disabled!");
#endif
  }
}

/** \brief Multiply two linear operators.
 *
 * Multiply two linear operators. This currently assumes
 * that the underlying implementation uses Epetra_CrsMatrix.
 *
 * \param[in] opl Left operator (assumed to be a Epetra_CrsMatrix)
 * \param[in] opr Right operator (assumed to be a Epetra_CrsMatrix)
 *
 * \returns Matrix product with a Epetra_CrsMatrix implementation
 */
const LinearOp explicitMultiply(const LinearOp &opl, const LinearOp &opr) {
  // if this is a blocked operator, multiply block by block
  // it is possible that not every factor in the product is blocked and these situations are handled
  // separately

  bool isBlockedL = isPhysicallyBlockedLinearOp(opl);
  bool isBlockedR = isPhysicallyBlockedLinearOp(opr);

  // both factors blocked
  if ((isBlockedL && isBlockedR)) {
    double scalarl = 0.0;
    bool transpl   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opl =
        getPhysicallyBlockedLinearOp(opl, &scalarl, &transpl);
    double scalarr = 0.0;
    bool transpr   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opr =
        getPhysicallyBlockedLinearOp(opr, &scalarr, &transpr);
    double scalar = scalarl * scalarr;

    int numRows   = blocked_opl->productRange()->numBlocks();
    int numCols   = blocked_opr->productDomain()->numBlocks();
    int numMiddle = blocked_opl->productDomain()->numBlocks();

    TEUCHOS_ASSERT(blocked_opr->productRange()->numBlocks() == numMiddle);

    RCP<Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_product =
        Thyra::defaultBlockedLinearOp<double>();
    blocked_product->beginBlockFill(numRows, numCols);
    for (int r = 0; r < numRows; ++r) {
      for (int c = 0; c < numCols; ++c) {
        LinearOp product_rc =
            explicitMultiply(blocked_opl->getBlock(r, 0), blocked_opr->getBlock(0, c));
        for (int m = 1; m < numMiddle; ++m) {
          LinearOp product_m =
              explicitMultiply(blocked_opl->getBlock(r, m), blocked_opr->getBlock(m, c));
          product_rc = explicitAdd(product_rc, product_m);
        }
        blocked_product->setBlock(r, c, Thyra::scale(scalar, product_rc));
      }
    }
    blocked_product->endBlockFill();
    return blocked_product;
  }

  // only left factor blocked
  if ((isBlockedL && !isBlockedR)) {
    double scalarl = 0.0;
    bool transpl   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opl =
        getPhysicallyBlockedLinearOp(opl, &scalarl, &transpl);
    double scalar = scalarl;

    int numRows   = blocked_opl->productRange()->numBlocks();
    int numCols   = 1;
    int numMiddle = 1;

    TEUCHOS_ASSERT(blocked_opl->productDomain()->numBlocks() == numMiddle);

    RCP<Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_product =
        Thyra::defaultBlockedLinearOp<double>();
    blocked_product->beginBlockFill(numRows, numCols);
    for (int r = 0; r < numRows; ++r) {
      LinearOp product_r = explicitMultiply(blocked_opl->getBlock(r, 0), opr);
      blocked_product->setBlock(r, 0, Thyra::scale(scalar, product_r));
    }
    blocked_product->endBlockFill();
    return blocked_product;
  }

  // only right factor blocked
  if ((!isBlockedL && isBlockedR)) {
    double scalarr = 0.0;
    bool transpr   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opr =
        getPhysicallyBlockedLinearOp(opr, &scalarr, &transpr);
    double scalar = scalarr;

    int numRows   = 1;
    int numCols   = blocked_opr->productDomain()->numBlocks();
    int numMiddle = 1;

    TEUCHOS_ASSERT(blocked_opr->productRange()->numBlocks() == numMiddle);

    RCP<Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_product =
        Thyra::defaultBlockedLinearOp<double>();
    blocked_product->beginBlockFill(numRows, numCols);
    for (int c = 0; c < numCols; ++c) {
      LinearOp product_c = explicitMultiply(opl, blocked_opr->getBlock(0, c));
      blocked_product->setBlock(0, c, Thyra::scale(scalar, product_c));
    }
    blocked_product->endBlockFill();
    return blocked_product;
  }

  bool isTpetral = Teko::TpetraHelpers::isTpetraLinearOp(opl);
  bool isTpetrar = Teko::TpetraHelpers::isTpetraLinearOp(opr);

  if (isTpetral && isTpetrar) {  // Both operators are Tpetra matrices so explicitly multiply them
    // Get left and right Tpetra crs operators
    ST scalarl   = 0.0;
    bool transpl = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpl =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opl, &scalarl, &transpl);
    ST scalarr   = 0.0;
    bool transpr = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpr =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opr, &scalarr, &transpr);

    // Build output operator
    RCP<Thyra::LinearOpBase<ST>> explicitOp = rcp(new Thyra::TpetraLinearOp<ST, LO, GO, NT>());
    RCP<Thyra::TpetraLinearOp<ST, LO, GO, NT>> tExplicitOp =
        rcp_dynamic_cast<Thyra::TpetraLinearOp<ST, LO, GO, NT>>(explicitOp);

    // Do explicit matrix-matrix multiply
    RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> explicitCrsOp =
        Tpetra::createCrsMatrix<ST, LO, GO, NT>(tCrsOpl->getRowMap());
    Tpetra::MatrixMatrix::Multiply<ST, LO, GO, NT>(*tCrsOpl, transpl, *tCrsOpr, transpr,
                                                   *explicitCrsOp);
    explicitCrsOp->resumeFill();
    explicitCrsOp->scale(scalarl * scalarr);
    explicitCrsOp->fillComplete(tCrsOpr->getDomainMap(), tCrsOpl->getRangeMap());
    tExplicitOp->initialize(Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getRangeMap()),
                            Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getDomainMap()),
                            explicitCrsOp);
    return tExplicitOp;

  } else if (isTpetral && !isTpetrar) {  // Assume that the right operator is diagonal

    // Get left Tpetra crs operator
    ST scalarl   = 0.0;
    bool transpl = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpl =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opl, &scalarl, &transpl);

    // Cast right operator as DiagonalLinearOp and extract diagonal as Vector
    RCP<const Thyra::DiagonalLinearOpBase<ST>> dOpr =
        rcp_dynamic_cast<const Thyra::DiagonalLinearOpBase<ST>>(opr, true);
    RCP<const Thyra::TpetraVector<ST, LO, GO, NT>> tPtr =
        rcp_dynamic_cast<const Thyra::TpetraVector<ST, LO, GO, NT>>(dOpr->getDiag(), true);
    RCP<const Tpetra::Vector<ST, LO, GO, NT>> diagPtr =
        rcp_dynamic_cast<const Tpetra::Vector<ST, LO, GO, NT>>(tPtr->getConstTpetraVector(), true);
    RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> explicitCrsOp =
        Tpetra::importAndFillCompleteCrsMatrix<Tpetra::CrsMatrix<ST, LO, GO, NT>>(
            tCrsOpl, Tpetra::Import<LO, GO, NT>(tCrsOpl->getRowMap(), tCrsOpl->getRowMap()));

    explicitCrsOp->rightScale(*diagPtr);
    explicitCrsOp->resumeFill();
    explicitCrsOp->scale(scalarl);
    explicitCrsOp->fillComplete(tCrsOpl->getDomainMap(), tCrsOpl->getRangeMap());

    return Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getRangeMap()),
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getDomainMap()), explicitCrsOp);

  } else if (!isTpetral && isTpetrar) {  // Assume that the left operator is diagonal

    // Get right Tpetra crs operator
    ST scalarr   = 0.0;
    bool transpr = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpr =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opr, &scalarr, &transpr);

    RCP<const Tpetra::Vector<ST, LO, GO, NT>> diagPtr;

    // Cast left operator as DiagonalLinearOp and extract diagonal as Vector
    RCP<const Thyra::DiagonalLinearOpBase<ST>> dOpl =
        rcp_dynamic_cast<const Thyra::DiagonalLinearOpBase<ST>>(opl);
    if (dOpl != Teuchos::null) {
      RCP<const Thyra::TpetraVector<ST, LO, GO, NT>> tPtr =
          rcp_dynamic_cast<const Thyra::TpetraVector<ST, LO, GO, NT>>(dOpl->getDiag(), true);
      diagPtr = rcp_dynamic_cast<const Tpetra::Vector<ST, LO, GO, NT>>(tPtr->getConstTpetraVector(),
                                                                       true);
    }
    // If it's not diagonal, maybe it's zero
    else if (rcp_dynamic_cast<const Thyra::ZeroLinearOpBase<ST>>(opl) != Teuchos::null) {
      diagPtr = rcp(new Tpetra::Vector<ST, LO, GO, NT>(tCrsOpr->getRangeMap()));
    } else
      TEUCHOS_ASSERT(false);

    RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> explicitCrsOp =
        Tpetra::importAndFillCompleteCrsMatrix<Tpetra::CrsMatrix<ST, LO, GO, NT>>(
            tCrsOpr, Tpetra::Import<LO, GO, NT>(tCrsOpr->getRowMap(), tCrsOpr->getRowMap()));

    explicitCrsOp->leftScale(*diagPtr);
    explicitCrsOp->resumeFill();
    explicitCrsOp->scale(scalarr);
    explicitCrsOp->fillComplete(tCrsOpr->getDomainMap(), tCrsOpr->getRangeMap());

    return Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getRangeMap()),
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getDomainMap()), explicitCrsOp);

  } else {  // Assume Epetra and we can use transformers
#ifdef TEKO_HAVE_EPETRA
    // build implicit multiply
    const LinearOp implicitOp = Thyra::multiply(opl, opr);

    // build a scaling transformer
    RCP<Thyra::LinearOpTransformerBase<double>> prodTrans =
        Thyra::epetraExtDiagScalingTransformer();

    // check to see if a scaling transformer works: if not use the
    // DiagScaledMatrixProduct transformer
    if (not prodTrans->isCompatible(*implicitOp))
      prodTrans = Thyra::epetraExtDiagScaledMatProdTransformer();

    // build operator and multiply
    const RCP<Thyra::LinearOpBase<double>> explicitOp = prodTrans->createOutputOp();
    prodTrans->transform(*implicitOp, explicitOp.ptr());
    explicitOp->setObjectLabel("explicit( " + opl->getObjectLabel() + " * " +
                               opr->getObjectLabel() + " )");

    return explicitOp;
#else
    throw std::logic_error(
        "explicitMultiply is trying to use Epetra "
        "code, but TEKO_HAVE_EPETRA is disabled!");
#endif
  }
}

/** \brief Multiply two linear operators.
 *
 * Multiply two linear operators. This currently assumes
 * that the underlying implementation uses Epetra_CrsMatrix.
 * The exception is that opm is allowed to be an diagonal matrix.
 *
 * \param[in] opl Left operator (assumed to be a Epetra_CrsMatrix)
 * \param[in] opr Right operator (assumed to be a Epetra_CrsMatrix)
 * \param[in,out] destOp The operator to be used as the destination operator,
 *                       if this is null this function creates a new operator
 *
 * \returns Matrix product with a Epetra_CrsMatrix implementation
 */
const ModifiableLinearOp explicitMultiply(const LinearOp &opl, const LinearOp &opr,
                                          const ModifiableLinearOp &destOp) {
  // if this is a blocked operator, multiply block by block
  // it is possible that not every factor in the product is blocked and these situations are handled
  // separately

  bool isBlockedL = isPhysicallyBlockedLinearOp(opl);
  bool isBlockedR = isPhysicallyBlockedLinearOp(opr);

  // both factors blocked
  if ((isBlockedL && isBlockedR)) {
    double scalarl = 0.0;
    bool transpl   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opl =
        getPhysicallyBlockedLinearOp(opl, &scalarl, &transpl);
    double scalarr = 0.0;
    bool transpr   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opr =
        getPhysicallyBlockedLinearOp(opr, &scalarr, &transpr);
    double scalar = scalarl * scalarr;

    int numRows   = blocked_opl->productRange()->numBlocks();
    int numCols   = blocked_opr->productDomain()->numBlocks();
    int numMiddle = blocked_opl->productDomain()->numBlocks();

    TEUCHOS_ASSERT(blocked_opr->productRange()->numBlocks() == numMiddle);

    RCP<Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_product =
        Thyra::defaultBlockedLinearOp<double>();
    blocked_product->beginBlockFill(numRows, numCols);
    for (int r = 0; r < numRows; ++r) {
      for (int c = 0; c < numCols; ++c) {
        LinearOp product_rc =
            explicitMultiply(blocked_opl->getBlock(r, 0), blocked_opr->getBlock(0, c));
        for (int m = 1; m < numMiddle; ++m) {
          LinearOp product_m =
              explicitMultiply(blocked_opl->getBlock(r, m), blocked_opr->getBlock(m, c));
          product_rc = explicitAdd(product_rc, product_m);
        }
        blocked_product->setBlock(r, c, Thyra::scale(scalar, product_rc));
      }
    }
    blocked_product->endBlockFill();
    return blocked_product;
  }

  // only left factor blocked
  if ((isBlockedL && !isBlockedR)) {
    double scalarl = 0.0;
    bool transpl   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opl =
        getPhysicallyBlockedLinearOp(opl, &scalarl, &transpl);
    double scalar = scalarl;

    int numRows   = blocked_opl->productRange()->numBlocks();
    int numCols   = 1;
    int numMiddle = 1;

    TEUCHOS_ASSERT(blocked_opl->productDomain()->numBlocks() == numMiddle);

    RCP<Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_product =
        Thyra::defaultBlockedLinearOp<double>();
    blocked_product->beginBlockFill(numRows, numCols);
    for (int r = 0; r < numRows; ++r) {
      LinearOp product_r = explicitMultiply(blocked_opl->getBlock(r, 0), opr);
      blocked_product->setBlock(r, 0, Thyra::scale(scalar, product_r));
    }
    blocked_product->endBlockFill();
    return blocked_product;
  }

  // only right factor blocked
  if ((!isBlockedL && isBlockedR)) {
    double scalarr = 0.0;
    bool transpr   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opr =
        getPhysicallyBlockedLinearOp(opr, &scalarr, &transpr);
    double scalar = scalarr;

    int numRows   = 1;
    int numCols   = blocked_opr->productDomain()->numBlocks();
    int numMiddle = 1;

    TEUCHOS_ASSERT(blocked_opr->productRange()->numBlocks() == numMiddle);

    RCP<Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_product =
        Thyra::defaultBlockedLinearOp<double>();
    blocked_product->beginBlockFill(numRows, numCols);
    for (int c = 0; c < numCols; ++c) {
      LinearOp product_c = explicitMultiply(opl, blocked_opr->getBlock(0, c));
      blocked_product->setBlock(0, c, Thyra::scale(scalar, product_c));
    }
    blocked_product->endBlockFill();
    return blocked_product;
  }

  bool isTpetral = Teko::TpetraHelpers::isTpetraLinearOp(opl);
  bool isTpetrar = Teko::TpetraHelpers::isTpetraLinearOp(opr);

  if (isTpetral && isTpetrar) {  // Both operators are Tpetra matrices so use the explicit Tpetra
                                 // matrix-matrix multiply

    // Get left and right Tpetra crs operators
    ST scalarl   = 0.0;
    bool transpl = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpl =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opl, &scalarl, &transpl);
    ST scalarr   = 0.0;
    bool transpr = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpr =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opr, &scalarr, &transpr);

    // Build output operator
    RCP<Thyra::LinearOpBase<ST>> explicitOp;
    if (destOp != Teuchos::null)
      explicitOp = destOp;
    else
      explicitOp = rcp(new Thyra::TpetraLinearOp<ST, LO, GO, NT>());
    RCP<Thyra::TpetraLinearOp<ST, LO, GO, NT>> tExplicitOp =
        rcp_dynamic_cast<Thyra::TpetraLinearOp<ST, LO, GO, NT>>(explicitOp);

    // Do explicit matrix-matrix multiply
    RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> explicitCrsOp =
        Tpetra::createCrsMatrix<ST, LO, GO, NT>(tCrsOpl->getRowMap());
    Tpetra::MatrixMatrix::Multiply<ST, LO, GO, NT>(*tCrsOpl, transpl, *tCrsOpr, transpr,
                                                   *explicitCrsOp);
    explicitCrsOp->resumeFill();
    explicitCrsOp->scale(scalarl * scalarr);
    explicitCrsOp->fillComplete(tCrsOpr->getDomainMap(), tCrsOpl->getRangeMap());
    tExplicitOp->initialize(Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getRangeMap()),
                            Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getDomainMap()),
                            explicitCrsOp);
    return tExplicitOp;

  } else if (isTpetral && !isTpetrar) {  // Assume that the right operator is diagonal

    // Get left Tpetra crs operator
    ST scalarl   = 0.0;
    bool transpl = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpl =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opl, &scalarl, &transpl);

    // Cast right operator as DiagonalLinearOp and extract diagonal as Vector
    RCP<const Thyra::DiagonalLinearOpBase<ST>> dOpr =
        rcp_dynamic_cast<const Thyra::DiagonalLinearOpBase<ST>>(opr);
    RCP<const Thyra::TpetraVector<ST, LO, GO, NT>> tPtr =
        rcp_dynamic_cast<const Thyra::TpetraVector<ST, LO, GO, NT>>(dOpr->getDiag(), true);
    RCP<const Tpetra::Vector<ST, LO, GO, NT>> diagPtr =
        rcp_dynamic_cast<const Tpetra::Vector<ST, LO, GO, NT>>(tPtr->getConstTpetraVector(), true);

    // Scale by the diagonal operator
    RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> explicitCrsOp =
        Tpetra::importAndFillCompleteCrsMatrix<Tpetra::CrsMatrix<ST, LO, GO, NT>>(
            tCrsOpl, Tpetra::Import<LO, GO, NT>(tCrsOpl->getRowMap(), tCrsOpl->getRowMap()));
    explicitCrsOp->rightScale(*diagPtr);
    explicitCrsOp->resumeFill();
    explicitCrsOp->scale(scalarl);
    explicitCrsOp->fillComplete(tCrsOpl->getDomainMap(), tCrsOpl->getRangeMap());
    return Thyra::tpetraLinearOp<ST, LO, GO, NT>(
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getRangeMap()),
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getDomainMap()), explicitCrsOp);

  } else if (!isTpetral && isTpetrar) {  // Assume that the left operator is diagonal

    // Get right Tpetra crs operator
    ST scalarr   = 0.0;
    bool transpr = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpr =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opr, &scalarr, &transpr);

    // Cast leftt operator as DiagonalLinearOp and extract diagonal as Vector
    RCP<const Thyra::DiagonalLinearOpBase<ST>> dOpl =
        rcp_dynamic_cast<const Thyra::DiagonalLinearOpBase<ST>>(opl, true);
    RCP<const Thyra::TpetraVector<ST, LO, GO, NT>> tPtr =
        rcp_dynamic_cast<const Thyra::TpetraVector<ST, LO, GO, NT>>(dOpl->getDiag(), true);
    RCP<const Tpetra::Vector<ST, LO, GO, NT>> diagPtr =
        rcp_dynamic_cast<const Tpetra::Vector<ST, LO, GO, NT>>(tPtr->getConstTpetraVector(), true);

    // Scale by the diagonal operator
    RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> explicitCrsOp =
        Tpetra::importAndFillCompleteCrsMatrix<Tpetra::CrsMatrix<ST, LO, GO, NT>>(
            tCrsOpr, Tpetra::Import<LO, GO, NT>(tCrsOpr->getRowMap(), tCrsOpr->getRowMap()));
    explicitCrsOp->leftScale(*diagPtr);
    explicitCrsOp->resumeFill();
    explicitCrsOp->scale(scalarr);
    explicitCrsOp->fillComplete(tCrsOpr->getDomainMap(), tCrsOpr->getRangeMap());
    return Thyra::tpetraLinearOp<ST, LO, GO, NT>(
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getRangeMap()),
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getDomainMap()), explicitCrsOp);

  } else {  // Assume Epetra and we can use transformers
#ifdef TEKO_HAVE_EPETRA
    // build implicit multiply
    const LinearOp implicitOp = Thyra::multiply(opl, opr);

    // build a scaling transformer

    RCP<Thyra::LinearOpTransformerBase<double>> prodTrans =
        Thyra::epetraExtDiagScalingTransformer();

    // check to see if a scaling transformer works: if not use the
    // DiagScaledMatrixProduct transformer
    if (not prodTrans->isCompatible(*implicitOp))
      prodTrans = Thyra::epetraExtDiagScaledMatProdTransformer();

    // build operator destination operator
    ModifiableLinearOp explicitOp;

    // if neccessary build a operator to put the explicit multiply into
    if (destOp == Teuchos::null)
      explicitOp = prodTrans->createOutputOp();
    else
      explicitOp = destOp;

    // perform multiplication
    prodTrans->transform(*implicitOp, explicitOp.ptr());

    // label it
    explicitOp->setObjectLabel("explicit( " + opl->getObjectLabel() + " * " +
                               opr->getObjectLabel() + " )");

    return explicitOp;
#else
    throw std::logic_error(
        "explicitMultiply is trying to use Epetra "
        "code, but TEKO_HAVE_EPETRA is disabled!");
#endif
  }
}

/** \brief Add two linear operators.
 *
 * Add two linear operators. This currently assumes
 * that the underlying implementation uses Epetra_CrsMatrix.
 *
 * \param[in] opl Left operator (assumed to be a Epetra_CrsMatrix)
 * \param[in] opr Right operator (assumed to be a Epetra_CrsMatrix)
 *
 * \returns Matrix sum with a Epetra_CrsMatrix implementation
 */
const LinearOp explicitAdd(const LinearOp &opl_in, const LinearOp &opr_in) {
  // if both blocked, add block by block
  if (isPhysicallyBlockedLinearOp(opl_in) && isPhysicallyBlockedLinearOp(opr_in)) {
    double scalarl = 0.0;
    bool transpl   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opl =
        getPhysicallyBlockedLinearOp(opl_in, &scalarl, &transpl);

    double scalarr = 0.0;
    bool transpr   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opr =
        getPhysicallyBlockedLinearOp(opr_in, &scalarr, &transpr);

    int numRows = blocked_opl->productRange()->numBlocks();
    int numCols = blocked_opl->productDomain()->numBlocks();
    TEUCHOS_ASSERT(blocked_opr->productRange()->numBlocks() == numRows);
    TEUCHOS_ASSERT(blocked_opr->productDomain()->numBlocks() == numCols);

    RCP<Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_sum =
        Thyra::defaultBlockedLinearOp<double>();
    blocked_sum->beginBlockFill(numRows, numCols);
    for (int r = 0; r < numRows; ++r)
      for (int c = 0; c < numCols; ++c)
        blocked_sum->setBlock(r, c,
                              explicitAdd(Thyra::scale(scalarl, blocked_opl->getBlock(r, c)),
                                          Thyra::scale(scalarr, blocked_opr->getBlock(r, c))));
    blocked_sum->endBlockFill();
    return blocked_sum;
  }

  // if only one is blocked, it must be 1x1
  LinearOp opl = opl_in;
  LinearOp opr = opr_in;
  if (isPhysicallyBlockedLinearOp(opl_in)) {
    double scalarl = 0.0;
    bool transpl   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opl =
        getPhysicallyBlockedLinearOp(opl_in, &scalarl, &transpl);
    TEUCHOS_ASSERT(blocked_opl->productRange()->numBlocks() == 1);
    TEUCHOS_ASSERT(blocked_opl->productDomain()->numBlocks() == 1);
    opl = Thyra::scale(scalarl, blocked_opl->getBlock(0, 0));
  }
  if (isPhysicallyBlockedLinearOp(opr_in)) {
    double scalarr = 0.0;
    bool transpr   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opr =
        getPhysicallyBlockedLinearOp(opr_in, &scalarr, &transpr);
    TEUCHOS_ASSERT(blocked_opr->productRange()->numBlocks() == 1);
    TEUCHOS_ASSERT(blocked_opr->productDomain()->numBlocks() == 1);
    opr = Thyra::scale(scalarr, blocked_opr->getBlock(0, 0));
  }

  bool isTpetral = Teko::TpetraHelpers::isTpetraLinearOp(opl);
  bool isTpetrar = Teko::TpetraHelpers::isTpetraLinearOp(opr);

  // if one of the operators in the sum is a thyra zero op
  if (isZeroOp(opl)) {
    if (isZeroOp(opr)) return opr;  // return a zero op if both are zero
    if (isTpetrar) {                // if other op is tpetra, replace this with a zero crs matrix
      ST scalar     = 0.0;
      bool transp   = false;
      auto crs_op   = Teko::TpetraHelpers::getTpetraCrsMatrix(opr, &scalar, &transp);
      auto zero_crs = Tpetra::createCrsMatrix<ST, LO, GO, NT>(crs_op->getRowMap());
      zero_crs->fillComplete();
      opl = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
          Thyra::tpetraVectorSpace<ST, LO, GO, NT>(crs_op->getRangeMap()),
          Thyra::tpetraVectorSpace<ST, LO, GO, NT>(crs_op->getDomainMap()), zero_crs);
      isTpetral = true;
    } else
      return opr->clone();
  }
  if (isZeroOp(opr)) {
    if (isTpetral) {  // if other op is tpetra, replace this with a zero crs matrix
      ST scalar     = 0.0;
      bool transp   = false;
      auto crs_op   = Teko::TpetraHelpers::getTpetraCrsMatrix(opr, &scalar, &transp);
      auto zero_crs = Tpetra::createCrsMatrix<ST, LO, GO, NT>(crs_op->getRowMap());
      zero_crs->fillComplete();
      opr = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
          Thyra::tpetraVectorSpace<ST, LO, GO, NT>(crs_op->getRangeMap()),
          Thyra::tpetraVectorSpace<ST, LO, GO, NT>(crs_op->getDomainMap()), zero_crs);
      isTpetrar = true;
    } else
      return opl->clone();
  }

  if (isTpetral && isTpetrar) {  // Both operators are Tpetra matrices so use the explicit Tpetra
                                 // matrix-matrix add

    // Get left and right Tpetra crs operators
    ST scalarl   = 0.0;
    bool transpl = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpl =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opl, &scalarl, &transpl);
    ST scalarr   = 0.0;
    bool transpr = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpr =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opr, &scalarr, &transpr);

    // Build output operator
    RCP<Thyra::LinearOpBase<ST>> explicitOp = rcp(new Thyra::TpetraLinearOp<ST, LO, GO, NT>());
    RCP<Thyra::TpetraLinearOp<ST, LO, GO, NT>> tExplicitOp =
        rcp_dynamic_cast<Thyra::TpetraLinearOp<ST, LO, GO, NT>>(explicitOp);

    // Do explicit matrix-matrix add
    RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> explicitCrsOp =
        Tpetra::MatrixMatrix::add<ST, LO, GO, NT>(scalarl, transpl, *tCrsOpl, scalarr, transpr,
                                                  *tCrsOpr);
    tExplicitOp->initialize(Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getRangeMap()),
                            Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getDomainMap()),
                            explicitCrsOp);
    return tExplicitOp;

  } else {  // Assume Epetra
#ifdef TEKO_HAVE_EPETRA
    // build implicit add
    const LinearOp implicitOp = Thyra::add(opl, opr);

    // build transformer
    const RCP<Thyra::LinearOpTransformerBase<double>> prodTrans = Thyra::epetraExtAddTransformer();

    // build operator and add
    const RCP<Thyra::LinearOpBase<double>> explicitOp = prodTrans->createOutputOp();
    prodTrans->transform(*implicitOp, explicitOp.ptr());
    explicitOp->setObjectLabel("explicit( " + opl->getObjectLabel() + " + " +
                               opr->getObjectLabel() + " )");

    return explicitOp;
#else
    throw std::logic_error(
        "explicitAdd is trying to use Epetra "
        "code, but TEKO_HAVE_EPETRA is disabled!");
#endif
  }
}

/** \brief Add two linear operators.
 *
 * Add two linear operators. This currently assumes
 * that the underlying implementation uses Epetra_CrsMatrix.
 *
 * \param[in] opl Left operator (assumed to be a Epetra_CrsMatrix)
 * \param[in] opr Right operator (assumed to be a Epetra_CrsMatrix)
 * \param[in,out] destOp The operator to be used as the destination operator,
 *                       if this is null this function creates a new operator
 *
 * \returns Matrix sum with a Epetra_CrsMatrix implementation
 */
const ModifiableLinearOp explicitAdd(const LinearOp &opl_in, const LinearOp &opr_in,
                                     const ModifiableLinearOp &destOp) {
  // if blocked, add block by block
  if (isPhysicallyBlockedLinearOp(opl_in) && isPhysicallyBlockedLinearOp(opr_in)) {
    double scalarl = 0.0;
    bool transpl   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opl =
        getPhysicallyBlockedLinearOp(opl_in, &scalarl, &transpl);

    double scalarr = 0.0;
    bool transpr   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opr =
        getPhysicallyBlockedLinearOp(opr_in, &scalarr, &transpr);

    int numRows = blocked_opl->productRange()->numBlocks();
    int numCols = blocked_opl->productDomain()->numBlocks();
    TEUCHOS_ASSERT(blocked_opr->productRange()->numBlocks() == numRows);
    TEUCHOS_ASSERT(blocked_opr->productDomain()->numBlocks() == numCols);

    RCP<Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_sum =
        Teuchos::rcp_dynamic_cast<Thyra::DefaultBlockedLinearOp<double>>(destOp);
    if (blocked_sum.is_null()) {
      // take care of the null case, this means we must alllocate memory
      blocked_sum = Thyra::defaultBlockedLinearOp<double>();

      blocked_sum->beginBlockFill(numRows, numCols);
      for (int r = 0; r < numRows; ++r) {
        for (int c = 0; c < numCols; ++c) {
          auto block =
              explicitAdd(Thyra::scale(scalarl, blocked_opl->getBlock(r, c)),
                          Thyra::scale(scalarr, blocked_opr->getBlock(r, c)), Teuchos::null);
          blocked_sum->setNonconstBlock(r, c, block);
        }
      }
      blocked_sum->endBlockFill();

    } else {
      // in this case memory can be reused
      for (int r = 0; r < numRows; ++r)
        for (int c = 0; c < numCols; ++c)
          explicitAdd(Thyra::scale(scalarl, blocked_opl->getBlock(r, c)),
                      Thyra::scale(scalarr, blocked_opr->getBlock(r, c)),
                      blocked_sum->getNonconstBlock(r, c));
    }

    return blocked_sum;
  }

  LinearOp opl = opl_in;
  LinearOp opr = opr_in;
  // if only one is blocked, it must be 1x1
  if (isPhysicallyBlockedLinearOp(opl)) {
    double scalarl = 0.0;
    bool transpl   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opl =
        getPhysicallyBlockedLinearOp(opl, &scalarl, &transpl);
    TEUCHOS_ASSERT(blocked_opl->productRange()->numBlocks() == 1);
    TEUCHOS_ASSERT(blocked_opl->productDomain()->numBlocks() == 1);
    opl = Thyra::scale(scalarl, blocked_opl->getBlock(0, 0));
  }
  if (isPhysicallyBlockedLinearOp(opr)) {
    double scalarr = 0.0;
    bool transpr   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_opr =
        getPhysicallyBlockedLinearOp(opr, &scalarr, &transpr);
    TEUCHOS_ASSERT(blocked_opr->productRange()->numBlocks() == 1);
    TEUCHOS_ASSERT(blocked_opr->productDomain()->numBlocks() == 1);
    opr = Thyra::scale(scalarr, blocked_opr->getBlock(0, 0));
  }

  bool isTpetral = Teko::TpetraHelpers::isTpetraLinearOp(opl);
  bool isTpetrar = Teko::TpetraHelpers::isTpetraLinearOp(opr);

  if (isTpetral && isTpetrar) {  // Both operators are Tpetra matrices so use the explicit
                                 // Tpetra matrix-matrix add

    // Get left and right Tpetra crs operators
    ST scalarl   = 0.0;
    bool transpl = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpl =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opl, &scalarl, &transpl);
    ST scalarr   = 0.0;
    bool transpr = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOpr =
        Teko::TpetraHelpers::getTpetraCrsMatrix(opr, &scalarr, &transpr);

    // Build output operator
    RCP<Thyra::LinearOpBase<ST>> explicitOp;
    RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> explicitCrsOp;
    if (!destOp.is_null()) {
      explicitOp = destOp;
      RCP<Thyra::TpetraLinearOp<ST, LO, GO, NT>> tOp =
          rcp_dynamic_cast<Thyra::TpetraLinearOp<ST, LO, GO, NT>>(destOp);
      if (!tOp.is_null())
        explicitCrsOp =
            rcp_dynamic_cast<Tpetra::CrsMatrix<ST, LO, GO, NT>>(tOp->getTpetraOperator());
      bool needNewTpetraMatrix =
          (explicitCrsOp.is_null()) || (tCrsOpl == explicitCrsOp) || (tCrsOpr == explicitCrsOp);
      if (!needNewTpetraMatrix) {
        try {
          // try to reuse matrix sparsity with Add. If it fails, build new operator with add
          Tpetra::MatrixMatrix::Add<ST, LO, GO, NT>(*tCrsOpl, transpl, scalarl, *tCrsOpr, transpr,
                                                    scalarr, explicitCrsOp);
        } catch (std::logic_error &e) {
          RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
          *out << "*** THROWN EXCEPTION ***\n";
          *out << e.what() << std::endl;
          *out << "************************\n";
          *out << "Teko: explicitAdd unable to reuse existing operator. Creating new operator.\n"
               << std::endl;
          needNewTpetraMatrix = true;
        }
      }
      if (needNewTpetraMatrix)
        // Do explicit matrix-matrix add
        explicitCrsOp = Tpetra::MatrixMatrix::add<ST, LO, GO, NT>(scalarl, transpl, *tCrsOpl,
                                                                  scalarr, transpr, *tCrsOpr);
    } else {
      explicitOp = rcp(new Thyra::TpetraLinearOp<ST, LO, GO, NT>());
      // Do explicit matrix-matrix add
      explicitCrsOp = Tpetra::MatrixMatrix::add<ST, LO, GO, NT>(scalarl, transpl, *tCrsOpl, scalarr,
                                                                transpr, *tCrsOpr);
    }
    RCP<Thyra::TpetraLinearOp<ST, LO, GO, NT>> tExplicitOp =
        rcp_dynamic_cast<Thyra::TpetraLinearOp<ST, LO, GO, NT>>(explicitOp);

    tExplicitOp->initialize(Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getRangeMap()),
                            Thyra::tpetraVectorSpace<ST, LO, GO, NT>(explicitCrsOp->getDomainMap()),
                            explicitCrsOp);
    return tExplicitOp;

  } else {  // Assume Epetra
#ifdef TEKO_HAVE_EPETRA
    // build implicit add
    const LinearOp implicitOp = Thyra::add(opl, opr);

    // build transformer
    const RCP<Thyra::LinearOpTransformerBase<double>> prodTrans = Thyra::epetraExtAddTransformer();

    // build or reuse destination operator
    RCP<Thyra::LinearOpBase<double>> explicitOp;
    if (destOp != Teuchos::null)
      explicitOp = destOp;
    else
      explicitOp = prodTrans->createOutputOp();

    // perform add
    prodTrans->transform(*implicitOp, explicitOp.ptr());
    explicitOp->setObjectLabel("explicit( " + opl->getObjectLabel() + " + " +
                               opr->getObjectLabel() + " )");

    return explicitOp;
#else
    throw std::logic_error(
        "explicitAdd is trying to use Epetra "
        "code, but TEKO_HAVE_EPETRA is disabled!");
#endif
  }
}

/** \brief Sum an operator
 *
 * \returns Matrix sum with a Epetra_CrsMatrix implementation
 */
const ModifiableLinearOp explicitSum(const LinearOp &op, const ModifiableLinearOp &destOp) {
#ifdef TEKO_HAVE_EPETRA
  // convert operators to Epetra_CrsMatrix
  const RCP<const Epetra_CrsMatrix> epetraOp =
      rcp_dynamic_cast<const Epetra_CrsMatrix>(get_Epetra_Operator(*op), true);

  if (destOp == Teuchos::null) {
    Teuchos::RCP<Epetra_Operator> epetraDest = Teuchos::rcp(new Epetra_CrsMatrix(*epetraOp));

    return Thyra::nonconstEpetraLinearOp(epetraDest);
  }

  const RCP<Epetra_CrsMatrix> epetraDest =
      rcp_dynamic_cast<Epetra_CrsMatrix>(get_Epetra_Operator(*destOp), true);

  EpetraExt::MatrixMatrix::Add(*epetraOp, false, 1.0, *epetraDest, 1.0);

  return destOp;
#else
  throw std::logic_error(
      "explicitSum is trying to use Epetra "
      "code, but TEKO_HAVE_EPETRA is disabled!");
#endif
}

const LinearOp explicitTranspose(const LinearOp &op) {
  if (Teko::TpetraHelpers::isTpetraLinearOp(op)) {
    RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT>> tOp =
        rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT>>(op, true);
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOp =
        rcp_dynamic_cast<const Tpetra::CrsMatrix<ST, LO, GO, NT>>(tOp->getConstTpetraOperator(),
                                                                  true);

    Tpetra::RowMatrixTransposer<ST, LO, GO, NT> transposer(tCrsOp);
    RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> transOp = transposer.createTranspose();

    return Thyra::tpetraLinearOp<ST, LO, GO, NT>(
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(transOp->getRangeMap()),
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(transOp->getDomainMap()), transOp);

  } else {
#ifdef TEKO_HAVE_EPETRA
    RCP<const Epetra_Operator> eOp = Thyra::get_Epetra_Operator(*op);
    TEUCHOS_TEST_FOR_EXCEPTION(eOp == Teuchos::null, std::logic_error,
                               "Teko::explicitTranspose Not an Epetra_Operator");
    RCP<const Epetra_RowMatrix> eRowMatrixOp =
        Teuchos::rcp_dynamic_cast<const Epetra_RowMatrix>(eOp);
    TEUCHOS_TEST_FOR_EXCEPTION(eRowMatrixOp == Teuchos::null, std::logic_error,
                               "Teko::explicitTranspose Not an Epetra_RowMatrix");

    // we now have a delete transpose operator
    EpetraExt::RowMatrix_Transpose tranposeOp;
    Epetra_RowMatrix &eMat = tranposeOp(const_cast<Epetra_RowMatrix &>(*eRowMatrixOp));

    // this copy is because of a poor implementation of the EpetraExt::Transform
    // implementation
    Teuchos::RCP<Epetra_CrsMatrix> crsMat =
        Teuchos::rcp(new Epetra_CrsMatrix(dynamic_cast<Epetra_CrsMatrix &>(eMat)));

    return Thyra::epetraLinearOp(crsMat);
#else
    throw std::logic_error(
        "explicitTranspose is trying to use Epetra "
        "code, but TEKO_HAVE_EPETRA is disabled!");
#endif
  }
}

const LinearOp explicitScale(double scalar, const LinearOp &op) {
  if (Teko::TpetraHelpers::isTpetraLinearOp(op)) {
    RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT>> tOp =
        rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT>>(op, true);
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOp =
        rcp_dynamic_cast<const Tpetra::CrsMatrix<ST, LO, GO, NT>>(tOp->getConstTpetraOperator(),
                                                                  true);
    auto crsOpNew = rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(*tCrsOp, Teuchos::Copy));
    crsOpNew->scale(scalar);
    return Thyra::tpetraLinearOp<ST, LO, GO, NT>(
        Thyra::createVectorSpace<ST, LO, GO, NT>(crsOpNew->getRangeMap()),
        Thyra::createVectorSpace<ST, LO, GO, NT>(crsOpNew->getDomainMap()), crsOpNew);
  } else {
#ifdef TEKO_HAVE_EPETRA
    RCP<const Thyra::EpetraLinearOp> eOp = rcp_dynamic_cast<const Thyra::EpetraLinearOp>(op, true);
    RCP<const Epetra_CrsMatrix> eCrsOp =
        rcp_dynamic_cast<const Epetra_CrsMatrix>(eOp->epetra_op(), true);
    Teuchos::RCP<Epetra_CrsMatrix> crsMat = Teuchos::rcp(new Epetra_CrsMatrix(*eCrsOp));

    crsMat->Scale(scalar);

    return Thyra::epetraLinearOp(crsMat);
#else
    throw std::logic_error(
        "explicitScale is trying to use Epetra "
        "code, but TEKO_HAVE_EPETRA is disabled!");
#endif
  }
}

double frobeniusNorm(const LinearOp &op_in) {
  LinearOp op;
  double scalar = 1.0;

  // if blocked, must be 1x1
  if (isPhysicallyBlockedLinearOp(op_in)) {
    bool transp = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_op =
        getPhysicallyBlockedLinearOp(op_in, &scalar, &transp);
    TEUCHOS_ASSERT(blocked_op->productRange()->numBlocks() == 1);
    TEUCHOS_ASSERT(blocked_op->productDomain()->numBlocks() == 1);
    op = blocked_op->getBlock(0, 0);
  } else
    op = op_in;

  if (Teko::TpetraHelpers::isTpetraLinearOp(op)) {
    const RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT>> tOp =
        rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT>>(op);
    const RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> crsOp =
        rcp_dynamic_cast<const Tpetra::CrsMatrix<ST, LO, GO, NT>>(tOp->getConstTpetraOperator(),
                                                                  true);
    return crsOp->getFrobeniusNorm();
  } else {
#ifdef TEKO_HAVE_EPETRA
    const RCP<const Epetra_Operator> epOp   = Thyra::get_Epetra_Operator(*op);
    const RCP<const Epetra_CrsMatrix> crsOp = rcp_dynamic_cast<const Epetra_CrsMatrix>(epOp, true);
    return crsOp->NormFrobenius();
#else
    throw std::logic_error(
        "frobeniusNorm is trying to use Epetra "
        "code, but TEKO_HAVE_EPETRA is disabled!");
#endif
  }
}

double oneNorm(const LinearOp &op) {
  if (Teko::TpetraHelpers::isTpetraLinearOp(op)) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "One norm not currently implemented for Tpetra matrices");

  } else {
#ifdef TEKO_HAVE_EPETRA
    const RCP<const Epetra_Operator> epOp   = Thyra::get_Epetra_Operator(*op);
    const RCP<const Epetra_CrsMatrix> crsOp = rcp_dynamic_cast<const Epetra_CrsMatrix>(epOp, true);
    return crsOp->NormOne();
#else
    throw std::logic_error(
        "oneNorm is trying to use Epetra "
        "code, but TEKO_HAVE_EPETRA is disabled!");
#endif
  }
}

double infNorm(const LinearOp &op) {
  if (Teko::TpetraHelpers::isTpetraLinearOp(op)) {
    ST scalar   = 0.0;
    bool transp = false;
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tCrsOp =
        Teko::TpetraHelpers::getTpetraCrsMatrix(op, &scalar, &transp);

    // extract diagonal
    const RCP<Tpetra::Vector<ST, LO, GO, NT>> ptrDiag =
        Tpetra::createVector<ST, LO, GO, NT>(tCrsOp->getRowMap());
    Tpetra::Vector<ST, LO, GO, NT> &diag = *ptrDiag;

    // compute absolute value row sum
    diag.putScalar(0.0);
    for (LO i = 0; i < (LO)tCrsOp->getLocalNumRows(); i++) {
      LO numEntries = tCrsOp->getNumEntriesInLocalRow(i);
      typename Tpetra::CrsMatrix<ST, LO, GO, NT>::local_inds_host_view_type indices;
      typename Tpetra::CrsMatrix<ST, LO, GO, NT>::values_host_view_type values;
      tCrsOp->getLocalRowView(i, indices, values);

      // build abs value row sum
      for (LO j = 0; j < numEntries; j++) diag.sumIntoLocalValue(i, std::abs(values(j)));
    }
    return diag.normInf() * scalar;

  } else {
#ifdef TEKO_HAVE_EPETRA
    const RCP<const Epetra_Operator> epOp   = Thyra::get_Epetra_Operator(*op);
    const RCP<const Epetra_CrsMatrix> crsOp = rcp_dynamic_cast<const Epetra_CrsMatrix>(epOp, true);
    return crsOp->NormInf();
#else
    throw std::logic_error(
        "infNorm is trying to use Epetra "
        "code, but TEKO_HAVE_EPETRA is disabled!");
#endif
  }
}

const LinearOp buildDiagonal(const MultiVector &src, const std::string &lbl) {
  RCP<Thyra::VectorBase<double>> dst = Thyra::createMember(src->range());
  Thyra::copy(*src->col(0), dst.ptr());

  return Thyra::diagonal<double>(dst, lbl);
}

const LinearOp buildInvDiagonal(const MultiVector &src, const std::string &lbl) {
  const RCP<const Thyra::VectorBase<double>> srcV = src->col(0);
  RCP<Thyra::VectorBase<double>> dst              = Thyra::createMember(srcV->range());
  Thyra::reciprocal<double>(*srcV, dst.ptr());

  return Thyra::diagonal<double>(dst, lbl);
}

//! build a BlockedMultiVector from a vector of MultiVectors
BlockedMultiVector buildBlockedMultiVector(const std::vector<MultiVector> &mvv) {
  Teuchos::Array<MultiVector> mvA;
  Teuchos::Array<VectorSpace> vsA;

  // build arrays of multi vectors and vector spaces
  std::vector<MultiVector>::const_iterator itr;
  for (itr = mvv.begin(); itr != mvv.end(); ++itr) {
    mvA.push_back(*itr);
    vsA.push_back((*itr)->range());
  }

  // construct the product vector space
  const RCP<const Thyra::DefaultProductVectorSpace<double>> vs =
      Thyra::productVectorSpace<double>(vsA);

  return Thyra::defaultProductMultiVector<double>(vs, mvA);
}

/** Construct an indicator vector specified by a vector of indices to
 * be set to ``on''.
 *
 * \param[in] indices Vector of indicies to turn on
 * \param[in] vs Vector space to construct the vector from
 * \param[in] onValue Value to set in the vector to on
 * \param[in] offValue Value to set in the vector to off
 *
 * \return Vector of on and off values.
 */
Teuchos::RCP<Thyra::VectorBase<double>> indicatorVector(const std::vector<int> &indices,
                                                        const VectorSpace &vs, double onValue,
                                                        double offValue)

{
  using Teuchos::RCP;

  // create a new vector
  RCP<Thyra::VectorBase<double>> v = Thyra::createMember<double>(vs);
  Thyra::put_scalar<double>(offValue, v.ptr());  // fill it with "off" values

  // set on values
  for (std::size_t i = 0; i < indices.size(); i++)
    Thyra::set_ele<double>(indices[i], onValue, v.ptr());

  return v;
}

/** \brief Compute the spectral radius of a matrix
 *
 * Compute the spectral radius of matrix A.  This utilizes the
 * Trilinos-Anasazi BlockKrylovShcur method for computing eigenvalues.
 * It attempts to compute the largest (in magnitude) eigenvalue to a given
 * level of tolerance.
 *
 * \param[in] A   matrix whose spectral radius is needed
 * \param[in] tol The <em>most</em> accuracy needed (the algorithm will run until
 *            it reaches this level of accuracy and then it will quit).
 *            If this level is not reached it will return something to indicate
 *            it has not converged.
 * \param[out] eigVec Eigenvectors
 * \param[in] isHermitian Is the matrix Hermitian
 * \param[in] numBlocks The size of the orthogonal basis built (like in GMRES) before
 *                  restarting.  Increase the memory usage by O(restart*n). At least
 *                  restart=3 is required.
 * \param[in] restart How many restarts are permitted
 * \param[in] verbosity See the Anasazi documentation
 *
 * \return The spectral radius of the matrix.  If the algorithm didn't converge the
 *         number is the negative of the ritz-values. If a <code>NaN</code> is returned
 *         there was a problem constructing the Anasazi problem
 */
double computeSpectralRad(const RCP<const Thyra::LinearOpBase<double>> &A, double tol,
                          bool isHermitian, int numBlocks, int restart, int verbosity) {
  typedef Thyra::LinearOpBase<double> OP;
  typedef Thyra::MultiVectorBase<double> MV;

  int startVectors = 1;

  // construct an initial guess
  const RCP<MV> ivec = Thyra::createMember(A->domain());
  Thyra::randomize(-1.0, 1.0, ivec.ptr());

  RCP<Anasazi::BasicEigenproblem<double, MV, OP>> eigProb =
      rcp(new Anasazi::BasicEigenproblem<double, MV, OP>(A, ivec));
  eigProb->setNEV(1);
  eigProb->setHermitian(isHermitian);

  // set the problem up
  if (not eigProb->setProblem()) {
    // big time failure!
    return Teuchos::ScalarTraits<double>::nan();
  }

  // we want largert magnitude eigenvalue
  std::string which("LM");  // largest magnitude

  // Create the parameter list for the eigensolver
  // verbosity+=Anasazi::TimingDetails;
  Teuchos::ParameterList MyPL;
  MyPL.set("Verbosity", verbosity);
  MyPL.set("Which", which);
  MyPL.set("Block Size", startVectors);
  MyPL.set("Num Blocks", numBlocks);
  MyPL.set("Maximum Restarts", restart);
  MyPL.set("Convergence Tolerance", tol);

  // build status test manager
  // RCP<Anasazi::StatusTestMaxIters<double,MV,OP> > statTest
  //       = rcp(new Anasazi::StatusTestMaxIters<double,MV,OP>(numBlocks*(restart+1)));

  // Create the Block Krylov Schur solver
  // This takes as inputs the eigenvalue problem and the solver parameters
  Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MyBlockKrylovSchur(eigProb, MyPL);

  // Solve the eigenvalue problem, and save the return code
  Anasazi::ReturnType solverreturn = MyBlockKrylovSchur.solve();

  if (solverreturn == Anasazi::Unconverged) {
    double real = MyBlockKrylovSchur.getRitzValues().begin()->realpart;
    double comp = MyBlockKrylovSchur.getRitzValues().begin()->imagpart;

    return -std::abs(std::sqrt(real * real + comp * comp));

    // cout << "Anasazi::BlockKrylovSchur::solve() did not converge!" << std::endl;
    // return -std::abs(MyBlockKrylovSchur.getRitzValues().begin()->realpart);
  } else {  // solverreturn==Anasazi::Converged
    double real = eigProb->getSolution().Evals.begin()->realpart;
    double comp = eigProb->getSolution().Evals.begin()->imagpart;

    return std::abs(std::sqrt(real * real + comp * comp));

    // cout << "Anasazi::BlockKrylovSchur::solve() converged!" << endl;
    // return std::abs(eigProb->getSolution().Evals.begin()->realpart);
  }
}

/** \brief Compute the smallest eigenvalue of an operator
 *
 * Compute the smallest eigenvalue of matrix A.  This utilizes the
 * Trilinos-Anasazi BlockKrylovShcur method for computing eigenvalues.
 * It attempts to compute the smallest (in magnitude) eigenvalue to a given
 * level of tolerance.
 *
 * \param[in] A   matrix whose spectral radius is needed
 * \param[in] tol The <em>most</em> accuracy needed (the algorithm will run until
 *            it reaches this level of accuracy and then it will quit).
 *            If this level is not reached it will return something to indicate
 *            it has not converged.
 * \param[in] isHermitian Is the matrix Hermitian
 * \param[in] numBlocks The size of the orthogonal basis built (like in GMRES) before
 *                  restarting.  Increase the memory usage by O(restart*n). At least
 *                  restart=3 is required.
 * \param[in] restart How many restarts are permitted
 * \param[in] verbosity See the Anasazi documentation
 *
 * \return The smallest magnitude eigenvalue of the matrix.  If the algorithm didn't converge
 * the number is the negative of the ritz-values. If a <code>NaN</code> is returned there was a
 * problem constructing the Anasazi problem
 */
double computeSmallestMagEig(const RCP<const Thyra::LinearOpBase<double>> &A, double tol,
                             bool isHermitian, int numBlocks, int restart, int verbosity) {
  typedef Thyra::LinearOpBase<double> OP;
  typedef Thyra::MultiVectorBase<double> MV;

  int startVectors = 1;

  // construct an initial guess
  const RCP<MV> ivec = Thyra::createMember(A->domain());
  Thyra::randomize(-1.0, 1.0, ivec.ptr());

  RCP<Anasazi::BasicEigenproblem<double, MV, OP>> eigProb =
      rcp(new Anasazi::BasicEigenproblem<double, MV, OP>(A, ivec));
  eigProb->setNEV(1);
  eigProb->setHermitian(isHermitian);

  // set the problem up
  if (not eigProb->setProblem()) {
    // big time failure!
    return Teuchos::ScalarTraits<double>::nan();
  }

  // we want largert magnitude eigenvalue
  std::string which("SM");  // smallest magnitude

  // Create the parameter list for the eigensolver
  Teuchos::ParameterList MyPL;
  MyPL.set("Verbosity", verbosity);
  MyPL.set("Which", which);
  MyPL.set("Block Size", startVectors);
  MyPL.set("Num Blocks", numBlocks);
  MyPL.set("Maximum Restarts", restart);
  MyPL.set("Convergence Tolerance", tol);

  // build status test manager
  // RCP<Anasazi::StatusTestMaxIters<double,MV,OP> > statTest
  //       = rcp(new Anasazi::StatusTestMaxIters<double,MV,OP>(10));

  // Create the Block Krylov Schur solver
  // This takes as inputs the eigenvalue problem and the solver parameters
  Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MyBlockKrylovSchur(eigProb, MyPL);

  // Solve the eigenvalue problem, and save the return code
  Anasazi::ReturnType solverreturn = MyBlockKrylovSchur.solve();

  if (solverreturn == Anasazi::Unconverged) {
    // cout << "Anasazi::BlockKrylovSchur::solve() did not converge! " << std::endl;
    return -std::abs(MyBlockKrylovSchur.getRitzValues().begin()->realpart);
  } else {  // solverreturn==Anasazi::Converged
    // cout << "Anasazi::BlockKrylovSchur::solve() converged!" << endl;
    return std::abs(eigProb->getSolution().Evals.begin()->realpart);
  }
}

/** Get a diagonal operator from a matrix. The mechanism for computing
 * the diagonal is specified by a <code>DiagonalType</code> arugment.
 *
 * \param[in] A <code>Epetra_CrsMatrix</code> to extract the diagonal from.
 * \param[in] dt Specifies the type of diagonal that is desired.
 *
 * \returns A diagonal operator.
 */
ModifiableLinearOp getDiagonalOp(const Teko::LinearOp &A, const DiagonalType &dt) {
  switch (dt) {
    case Diagonal: return getDiagonalOp(A);
    case Lumped: return getLumpedMatrix(A);
    case AbsRowSum: return getAbsRowSumMatrix(A);
    case NotDiag:
    default: TEUCHOS_TEST_FOR_EXCEPT(true);
  };

  return Teuchos::null;
}

/** Get the inverse of a diagonal operator from a matrix. The mechanism for computing
 * the diagonal is specified by a <code>DiagonalType</code> arugment.
 *
 * \param[in] A <code>Epetra_CrsMatrix</code> to extract the diagonal from.
 * \param[in] dt Specifies the type of diagonal that is desired.
 *
 * \returns A inverse of a diagonal operator.
 */
ModifiableLinearOp getInvDiagonalOp(const Teko::LinearOp &A, const Teko::DiagonalType &dt) {
  switch (dt) {
    case Diagonal: return getInvDiagonalOp(A);
    case Lumped: return getInvLumpedMatrix(A);
    case AbsRowSum: return getAbsRowSumInvMatrix(A);
    case NotDiag:
    default: TEUCHOS_TEST_FOR_EXCEPT(true);
  };

  return Teuchos::null;
}

/** Get a string corresponding to the type of diagonal specified.
 *
 * \param[in] dt The type of diagonal.
 *
 * \returns A string name representing this diagonal type.
 */
std::string getDiagonalName(const DiagonalType &dt) {
  switch (dt) {
    case Diagonal: return "Diagonal";
    case Lumped: return "Lumped";
    case AbsRowSum: return "AbsRowSum";
    case NotDiag: return "NotDiag";
    case BlkDiag: return "BlkDiag";
  };

  return "<error>";
}

/** Get a type corresponding to the name of a diagonal specified.
 *
 * \param[in] name String representing the diagonal type
 *
 * \returns The type representation of the string, if the
 *          string is not recognized this function returns
 *          a <code>NotDiag</code>
 */
DiagonalType getDiagonalType(std::string name) {
  if (name == "Diagonal") return Diagonal;
  if (name == "Lumped") return Lumped;
  if (name == "AbsRowSum") return AbsRowSum;
  if (name == "BlkDiag") return BlkDiag;

  return NotDiag;
}

#ifdef TEKO_HAVE_EPETRA
LinearOp probe(Teuchos::RCP<const Epetra_CrsGraph> &G, const LinearOp &Op) {
#ifdef Teko_ENABLE_Isorropia
  Teuchos::ParameterList probeList;
  Prober prober(G, probeList, true);
  Teuchos::RCP<Epetra_CrsMatrix> Mat = rcp(new Epetra_CrsMatrix(Copy, *G));
  Teko::Epetra::EpetraOperatorWrapper Mwrap(Op);
  prober.probe(Mwrap, *Mat);
  return Thyra::epetraLinearOp(Mat);
#else
  (void)G;
  (void)Op;
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Probe requires Isorropia");
#endif
}
#endif

double norm_1(const MultiVector &v, std::size_t col) {
  Teuchos::Array<double> n(v->domain()->dim());
  Thyra::norms_1<double>(*v, n);

  return n[col];
}

double norm_2(const MultiVector &v, std::size_t col) {
  Teuchos::Array<double> n(v->domain()->dim());
  Thyra::norms_2<double>(*v, n);

  return n[col];
}

#ifdef TEKO_HAVE_EPETRA
void putScalar(const ModifiableLinearOp &op, double scalar) {
  try {
    // get Epetra_Operator
    RCP<Epetra_Operator> eOp = Thyra::get_Epetra_Operator(*op);

    // cast it to a CrsMatrix
    RCP<Epetra_CrsMatrix> eCrsOp = rcp_dynamic_cast<Epetra_CrsMatrix>(eOp, true);

    eCrsOp->PutScalar(scalar);
  } catch (std::exception &e) {
    RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

    *out << "Teko: putScalar requires an Epetra_CrsMatrix\n";
    *out << "    Could not extract an Epetra_Operator from a \"" << op->description() << std::endl;
    *out << "           OR\n";
    *out << "    Could not cast an Epetra_Operator to a Epetra_CrsMatrix\n";
    *out << std::endl;
    *out << "*** THROWN EXCEPTION ***\n";
    *out << e.what() << std::endl;
    *out << "************************\n";

    throw e;
  }
}
#endif

void clipLower(MultiVector &v, double lowerBound) {
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // cast so entries are accessible
  // RCP<Thyra::SpmdMultiVectorBase<double> > spmdMVec
  //    = rcp_dynamic_cast<Thyra::DefaultSpmdMultiVector<double> >(v);

  for (Thyra::Ordinal i = 0; i < v->domain()->dim(); i++) {
    RCP<Thyra::SpmdVectorBase<double>> spmdVec =
        rcp_dynamic_cast<Thyra::SpmdVectorBase<double>>(v->col(i), true);

    Teuchos::ArrayRCP<double> values;
    // spmdMVec->getNonconstLocalData(Teuchos::ptrFromRef(values),Teuchos::ptrFromRef(i));
    spmdVec->getNonconstLocalData(Teuchos::ptrFromRef(values));
    for (Teuchos::ArrayRCP<double>::size_type j = 0; j < values.size(); j++) {
      if (values[j] < lowerBound) values[j] = lowerBound;
    }
  }
}

void clipUpper(MultiVector &v, double upperBound) {
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // cast so entries are accessible
  // RCP<Thyra::SpmdMultiVectorBase<double> > spmdMVec
  //   = rcp_dynamic_cast<Thyra::DefaultSpmdMultiVector<double> >(v);
  for (Thyra::Ordinal i = 0; i < v->domain()->dim(); i++) {
    RCP<Thyra::SpmdVectorBase<double>> spmdVec =
        rcp_dynamic_cast<Thyra::SpmdVectorBase<double>>(v->col(i), true);

    Teuchos::ArrayRCP<double> values;
    // spmdMVec->getNonconstLocalData(Teuchos::ptrFromRef(values),Teuchos::ptrFromRef(i));
    spmdVec->getNonconstLocalData(Teuchos::ptrFromRef(values));
    for (Teuchos::ArrayRCP<double>::size_type j = 0; j < values.size(); j++) {
      if (values[j] > upperBound) values[j] = upperBound;
    }
  }
}

void replaceValue(MultiVector &v, double currentValue, double newValue) {
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // cast so entries are accessible
  // RCP<Thyra::SpmdMultiVectorBase<double> > spmdMVec
  //    = rcp_dynamic_cast<Thyra::SpmdMultiVectorBase<double> >(v,true);
  for (Thyra::Ordinal i = 0; i < v->domain()->dim(); i++) {
    RCP<Thyra::SpmdVectorBase<double>> spmdVec =
        rcp_dynamic_cast<Thyra::SpmdVectorBase<double>>(v->col(i), true);

    Teuchos::ArrayRCP<double> values;
    // spmdMVec->getNonconstLocalData(Teuchos::ptrFromRef(values),Teuchos::ptrFromRef(i));
    spmdVec->getNonconstLocalData(Teuchos::ptrFromRef(values));
    for (Teuchos::ArrayRCP<double>::size_type j = 0; j < values.size(); j++) {
      if (values[j] == currentValue) values[j] = newValue;
    }
  }
}

void columnAverages(const MultiVector &v, std::vector<double> &averages) {
  averages.resize(v->domain()->dim());

  // sum over each column
  Thyra::sums<double>(*v, averages);

  // build averages
  Thyra::Ordinal rows = v->range()->dim();
  for (std::size_t i = 0; i < averages.size(); i++) averages[i] = averages[i] / rows;
}

double average(const MultiVector &v) {
  Thyra::Ordinal rows = v->range()->dim();
  Thyra::Ordinal cols = v->domain()->dim();

  std::vector<double> averages;
  columnAverages(v, averages);

  double sum = 0.0;
  for (std::size_t i = 0; i < averages.size(); i++) sum += averages[i] * rows;

  return sum / (rows * cols);
}

bool isPhysicallyBlockedLinearOp(const LinearOp &op) {
  // See if the operator is a PBLO
  RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> pblo =
      rcp_dynamic_cast<const Thyra::PhysicallyBlockedLinearOpBase<double>>(op);
  if (!pblo.is_null()) return true;

  // See if the operator is a wrapped PBLO
  ST scalar               = 0.0;
  Thyra::EOpTransp transp = Thyra::NOTRANS;
  RCP<const Thyra::LinearOpBase<ST>> wrapped_op;
  Thyra::unwrap(op, &scalar, &transp, &wrapped_op);
  pblo = rcp_dynamic_cast<const Thyra::PhysicallyBlockedLinearOpBase<double>>(wrapped_op);
  if (!pblo.is_null()) return true;

  return false;
}

RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> getPhysicallyBlockedLinearOp(
    const LinearOp &op, ST *scalar, bool *transp) {
  // If the operator is a TpetraLinearOp
  RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> pblo =
      rcp_dynamic_cast<const Thyra::PhysicallyBlockedLinearOpBase<double>>(op);
  if (!pblo.is_null()) {
    *scalar = 1.0;
    *transp = false;
    return pblo;
  }

  // If the operator is a wrapped TpetraLinearOp
  RCP<const Thyra::LinearOpBase<ST>> wrapped_op;
  Thyra::EOpTransp eTransp = Thyra::NOTRANS;
  Thyra::unwrap(op, scalar, &eTransp, &wrapped_op);
  pblo = rcp_dynamic_cast<const Thyra::PhysicallyBlockedLinearOpBase<double>>(wrapped_op, true);
  if (!pblo.is_null()) {
    *transp = true;
    if (eTransp == Thyra::NOTRANS) *transp = false;
    return pblo;
  }

  return Teuchos::null;
}

std::string formatBlockName(const std::string &prefix, int i, int j, int nrow) {
  unsigned digits = 0;
  auto blockId    = nrow - 1;
  do {
    blockId /= 10;
    digits++;
  } while (blockId);

  std::ostringstream ss;
  ss << prefix << "_";
  ss << std::setfill('0') << std::setw(digits) << i;
  ss << "_";
  ss << std::setfill('0') << std::setw(digits) << j;
  ss << ".mm";
  return ss.str();
}

void writeMatrix(const std::string &filename, const Teko::LinearOp &op) {
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
#ifdef TEKO_HAVE_EPETRA
  const RCP<const Thyra::EpetraLinearOp> eOp = rcp_dynamic_cast<const Thyra::EpetraLinearOp>(op);
#endif

  if (Teko::TpetraHelpers::isTpetraLinearOp(op)) {
    const RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT>> tOp =
        rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT>>(op);
    const RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> crsOp =
        rcp_dynamic_cast<const Tpetra::CrsMatrix<ST, LO, GO, NT>>(tOp->getConstTpetraOperator(),
                                                                  true);
    using Writer = Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<ST, LO, GO, NT>>;
    Writer::writeMapFile(("rowmap_" + filename).c_str(), *(crsOp->getRowMap()));
    Writer::writeMapFile(("colmap_" + filename).c_str(), *(crsOp->getColMap()));
    Writer::writeMapFile(("domainmap_" + filename).c_str(), *(crsOp->getDomainMap()));
    Writer::writeMapFile(("rangemap_" + filename).c_str(), *(crsOp->getRangeMap()));
    Writer::writeSparseFile(filename.c_str(), crsOp);
  }
#ifdef TEKO_HAVE_EPETRA
  else if (eOp != Teuchos::null) {
    const RCP<const Epetra_CrsMatrix> crsOp =
        rcp_dynamic_cast<const Epetra_CrsMatrix>(eOp->epetra_op(), true);
    EpetraExt::BlockMapToMatrixMarketFile(("rowmap_" + filename).c_str(), crsOp->RowMap());
    EpetraExt::BlockMapToMatrixMarketFile(("colmap_" + filename).c_str(), crsOp->ColMap());
    EpetraExt::BlockMapToMatrixMarketFile(("domainmap_" + filename).c_str(), crsOp->DomainMap());
    EpetraExt::BlockMapToMatrixMarketFile(("rangemap_" + filename).c_str(), crsOp->RangeMap());
    EpetraExt::RowMatrixToMatrixMarketFile(filename.c_str(), *crsOp);
  }
#endif
  else if (isPhysicallyBlockedLinearOp(op)) {
    double scalar = 0.0;
    bool transp   = false;
    RCP<const Thyra::PhysicallyBlockedLinearOpBase<double>> blocked_op =
        getPhysicallyBlockedLinearOp(op, &scalar, &transp);

    int numRows = blocked_op->productRange()->numBlocks();
    int numCols = blocked_op->productDomain()->numBlocks();

    for (int r = 0; r < numRows; ++r)
      for (int c = 0; c < numCols; ++c) {
        auto block = Teko::explicitScale(scalar, blocked_op->getBlock(r, c));
        if (transp) block = Teko::explicitTranspose(block);
        writeMatrix(formatBlockName(filename, r, c, numRows), block);
      }
  } else {
    TEUCHOS_ASSERT(false);
  }
}

}  // namespace Teko

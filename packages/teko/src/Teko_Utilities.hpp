/*
// @HEADER
//
// ***********************************************************************
//
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

*/

/** \file Teko_Utilities.hpp
 *
 * This file contains a number of useful functions and classes
 * used in Teko. They are distinct from the core functionality of
 * the preconditioner factory, however, the functions are critical
 * to construction of the preconditioners themselves.
 */

#ifndef __Teko_Utilities_hpp__
#define __Teko_Utilities_hpp__

#include "Teko_ConfigDefs.hpp"

#ifdef TEKO_HAVE_EPETRA
#include "Epetra_CrsMatrix.h"
#endif

#include "Tpetra_CrsMatrix.hpp"

// Teuchos includes
#include "Teuchos_VerboseObject.hpp"

// Thyra includes
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"

#ifdef _MSC_VER
#ifndef _MSC_EXTENSIONS
#define _MSC_EXTENSIONS
#define TEKO_DEFINED_MSC_EXTENSIONS
#endif
#include <iso646.h>  // For C alternative tokens
#endif

// #define Teko_DEBUG_OFF
#define Teko_DEBUG_INT 5

namespace Teko {

using Thyra::add;
using Thyra::block1x2;
using Thyra::block2x1;
using Thyra::block2x2;
using Thyra::identity;
using Thyra::multiply;
using Thyra::scale;
using Thyra::zero;  // make it to take one argument (square matrix)

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
Teuchos::RCP<Epetra_CrsMatrix> buildGraphLaplacian(int dim, double *coords,
                                                   const Epetra_CrsMatrix &stencil);
#endif

Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > buildGraphLaplacian(
    int dim, ST *coords, const Tpetra::CrsMatrix<ST, LO, GO, NT> &stencil);

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
Teuchos::RCP<Epetra_CrsMatrix> buildGraphLaplacian(double *x, double *y, double *z, int stride,
                                                   const Epetra_CrsMatrix &stencil);
#endif

Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > buildGraphLaplacian(
    ST *x, ST *y, ST *z, GO stride, const Tpetra::CrsMatrix<ST, LO, GO, NT> &stencil);

/** \brief Function used internally by Teko to find the output stream.
 *
 * Function used internally by Teko to find the output stream.
 *
 * \returns An output stream to use for printing
 */
const Teuchos::RCP<Teuchos::FancyOStream> getOutputStream();
// inline const Teuchos::RCP<Teuchos::FancyOStream> getOutputStream();
// { return Teuchos::VerboseObjectBase::getDefaultOStream(); }

#ifndef Teko_DEBUG_OFF
//#if 0
#define Teko_DEBUG_EXPR(str) str
#define Teko_DEBUG_MSG(str, level)                                     \
  if (level <= Teko_DEBUG_INT) {                                       \
    Teuchos::RCP<Teuchos::FancyOStream> out = Teko::getOutputStream(); \
    *out << "Teko: " << str << std::endl;                              \
  }
#define Teko_DEBUG_MSG_BEGIN(level)                        \
  if (level <= Teko_DEBUG_INT) {                           \
    Teko::getOutputStream()->pushTab(3);                   \
    *Teko::getOutputStream() << "Teko: Begin debug MSG\n"; \
    std::ostream &DEBUG_STREAM = *Teko::getOutputStream(); \
    Teko::getOutputStream()->pushTab(3);
#define Teko_DEBUG_MSG_END()                           \
  Teko::getOutputStream()->popTab();                   \
  *Teko::getOutputStream() << "Teko: End debug MSG\n"; \
  Teko::getOutputStream()->popTab();                   \
  }
#define Teko_DEBUG_PUSHTAB() Teko::getOutputStream()->pushTab(3)
#define Teko_DEBUG_POPTAB() Teko::getOutputStream()->popTab()
#define Teko_DEBUG_SCOPE(str, level)

//   struct __DebugScope__ {
//      __DebugScope__(const std::string & str,int level)
//         : str_(str), level_(level)
//      { Teko_DEBUG_MSG("BEGIN "+str_,level_); Teko_DEBUG_PUSHTAB(); }
//      ~__DebugScope__()
//      { Teko_DEBUG_POPTAB(); Teko_DEBUG_MSG("END "+str_,level_); }
//      std::string str_; int level_; };
//   #define Teko_DEBUG_SCOPE(str,level) __DebugScope__ __dbgScope__(str,level);
#else
#define Teko_DEBUG_EXPR(str)
#define Teko_DEBUG_MSG(str, level)
#define Teko_DEBUG_MSG_BEGIN(level) \
  if (false) {                      \
    std::ostream &DEBUG_STREAM = *Teko::getOutputStream();
#define Teko_DEBUG_MSG_END() }
#define Teko_DEBUG_PUSHTAB()
#define Teko_DEBUG_POPTAB()
#define Teko_DEBUG_SCOPE(str, level)
#endif

// typedefs for increased simplicity
typedef Teuchos::RCP<const Thyra::VectorSpaceBase<double> > VectorSpace;

// ----------------------------------------------------------------------------

//! @name MultiVector utilities
//@{

typedef Teuchos::RCP<Thyra::ProductMultiVectorBase<double> > BlockedMultiVector;
typedef Teuchos::RCP<Thyra::MultiVectorBase<double> > MultiVector;

//! Convert to a MultiVector from a BlockedMultiVector
inline MultiVector toMultiVector(BlockedMultiVector &bmv) { return bmv; }

//! Convert to a MultiVector from a BlockedMultiVector
inline const MultiVector toMultiVector(const BlockedMultiVector &bmv) { return bmv; }

//! Convert to a BlockedMultiVector from a MultiVector
inline const BlockedMultiVector toBlockedMultiVector(const MultiVector &bmv) {
  return Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<double> >(bmv);
}

//! Get the column count in a block linear operator
inline int blockCount(const BlockedMultiVector &bmv) { return bmv->productSpace()->numBlocks(); }

//! Get the <code>i</code>th block from a BlockedMultiVector object
inline MultiVector getBlock(int i, const BlockedMultiVector &bmv) {
  return Teuchos::rcp_const_cast<Thyra::MultiVectorBase<double> >(bmv->getMultiVectorBlock(i));
}

//! Perform a deep copy of the vector
inline MultiVector deepcopy(const MultiVector &v) { return v->clone_mv(); }

//! Perform a deep copy of the vector
inline MultiVector copyAndInit(const MultiVector &v, double scalar) {
  MultiVector mv = v->clone_mv();
  Thyra::assign(mv.ptr(), scalar);
  return mv;
}

//! Perform a deep copy of the blocked vector
inline BlockedMultiVector deepcopy(const BlockedMultiVector &v) {
  return toBlockedMultiVector(v->clone_mv());
}

/** \brief Copy the contents of a multivector to a destination vector.
 *
 * Copy the contents of a multivector to a new vector. If the destination
 * vector is null, a deep copy of the source multivector is made to a newly allocated
 * vector. Also, if the destination and the source do not match, a new destination
 * object is allocated and returned to the user.
 *
 * \param[in] src Source multivector to be copied.
 * \param[in] dst Destination multivector.  If null a new multivector will be allocated.
 *
 * \returns A copy of the source multivector. If dst is not null a pointer to this object
 *          is returned. Otherwise a new multivector is returned.
 */
inline MultiVector datacopy(const MultiVector &src, MultiVector &dst) {
  if (dst == Teuchos::null) return deepcopy(src);

  bool rangeCompat  = src->range()->isCompatible(*dst->range());
  bool domainCompat = src->domain()->isCompatible(*dst->domain());

  if (not(rangeCompat && domainCompat)) return deepcopy(src);

  // perform data copy
  Thyra::assign<double>(dst.ptr(), *src);
  return dst;
}

/** \brief Copy the contents of a blocked multivector to a destination vector.
 *
 * Copy the contents of a blocked multivector to a new vector. If the destination
 * vector is null, a deep copy of the source multivector is made to a newly allocated
 * vector. Also, if the destination and the source do not match, a new destination
 * object is allocated and returned to the user.
 *
 * \param[in] src Source multivector to be copied.
 * \param[in] dst Destination multivector.  If null a new multivector will be allocated.
 *
 * \returns A copy of the source multivector. If dst is not null a pointer to this object
 *          is returned. Otherwise a new multivector is returned.
 */
inline BlockedMultiVector datacopy(const BlockedMultiVector &src, BlockedMultiVector &dst) {
  if (dst == Teuchos::null) return deepcopy(src);

  bool rangeCompat  = src->range()->isCompatible(*dst->range());
  bool domainCompat = src->domain()->isCompatible(*dst->domain());

  if (not(rangeCompat && domainCompat)) return deepcopy(src);

  // perform data copy
  Thyra::assign<double>(dst.ptr(), *src);
  return dst;
}

//! build a BlockedMultiVector from a vector of MultiVectors
BlockedMultiVector buildBlockedMultiVector(const std::vector<MultiVector> &mvs);

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
Teuchos::RCP<Thyra::VectorBase<double> > indicatorVector(const std::vector<int> &indices,
                                                         const VectorSpace &vs,
                                                         double onValue  = 1.0,
                                                         double offValue = 0.0);

//@}

// ----------------------------------------------------------------------------

//! @name LinearOp utilities
//@{
typedef Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<ST> > BlockedLinearOp;
typedef Teuchos::RCP<const Thyra::LinearOpBase<ST> > LinearOp;
typedef Teuchos::RCP<Thyra::LinearOpBase<ST> > InverseLinearOp;
typedef Teuchos::RCP<Thyra::LinearOpBase<ST> > ModifiableLinearOp;

//! Build a square zero operator from a single vector space
inline LinearOp zero(const VectorSpace &vs) { return Thyra::zero<ST>(vs, vs); }

inline LinearOp zero(const VectorSpace &range, const VectorSpace &domain) {
  return Thyra::zero<ST>(range, domain);
}

//! Replace nonzeros with a scalar value, used to zero out an operator
#ifdef TEKO_HAVE_EPETRA
void putScalar(const ModifiableLinearOp &op, double scalar);
#endif

//! Get the range space of a linear operator
inline VectorSpace rangeSpace(const LinearOp &lo) { return lo->range(); }

//! Get the domain space of a linear operator
inline VectorSpace domainSpace(const LinearOp &lo) { return lo->domain(); }

//! Converts a LinearOp to a BlockedLinearOp
inline BlockedLinearOp toBlockedLinearOp(LinearOp &clo) {
  Teuchos::RCP<Thyra::LinearOpBase<double> > lo =
      Teuchos::rcp_const_cast<Thyra::LinearOpBase<double> >(clo);
  return Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<double> >(lo);
}

//! Converts a LinearOp to a BlockedLinearOp
inline const BlockedLinearOp toBlockedLinearOp(const LinearOp &clo) {
  Teuchos::RCP<Thyra::LinearOpBase<double> > lo =
      Teuchos::rcp_const_cast<Thyra::LinearOpBase<double> >(clo);
  return Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<double> >(lo);
}

//! Convert to a LinearOp from a BlockedLinearOp
inline LinearOp toLinearOp(BlockedLinearOp &blo) { return blo; }

//! Convert to a LinearOp from a BlockedLinearOp
inline const LinearOp toLinearOp(const BlockedLinearOp &blo) { return blo; }

//! Convert to a LinearOp from a BlockedLinearOp
inline LinearOp toLinearOp(ModifiableLinearOp &blo) { return blo; }

//! Convert to a LinearOp from a BlockedLinearOp
inline const LinearOp toLinearOp(const ModifiableLinearOp &blo) { return blo; }

//! Get the row count in a block linear operator
inline int blockRowCount(const BlockedLinearOp &blo) { return blo->productRange()->numBlocks(); }

//! Get the column count in a block linear operator
inline int blockColCount(const BlockedLinearOp &blo) { return blo->productDomain()->numBlocks(); }

//! Get the <code>i,j</code> block in a BlockedLinearOp object
inline LinearOp getBlock(int i, int j, const BlockedLinearOp &blo) { return blo->getBlock(i, j); }

//! Set the <code>i,j</code> block in a BlockedLinearOp object
inline void setBlock(int i, int j, BlockedLinearOp &blo, const LinearOp &lo) {
  return blo->setBlock(i, j, lo);
}

//! Build a new blocked linear operator
inline BlockedLinearOp createBlockedOp() {
  return rcp(new Thyra::DefaultBlockedLinearOp<double>());
}

/** \brief Let the blocked operator know that you are going to
 *        set the sub blocks.
 *
 * Let the blocked operator know that you are going to
 * set the sub blocks. This is a simple wrapper around the
 * member function of the same name in Thyra.
 *
 * \param[in,out] blo Blocked operator to have its fill stage activated
 * \param[in]  rowCnt Number of block rows in this operator
 * \param[in]  colCnt Number of block columns in this operator
 */
inline void beginBlockFill(BlockedLinearOp &blo, int rowCnt, int colCnt) {
  blo->beginBlockFill(rowCnt, colCnt);
}

/** \brief Let the blocked operator know that you are going to
 *        set the sub blocks.
 *
 * Let the blocked operator know that you are going to
 * set the sub blocks. This is a simple wrapper around the
 * member function of the same name in Thyra.
 *
 * \param[in,out] blo Blocked operator to have its fill stage activated
 */
inline void beginBlockFill(BlockedLinearOp &blo) { blo->beginBlockFill(); }

//! Notify the blocked operator that the fill stage is completed.
inline void endBlockFill(BlockedLinearOp &blo) { blo->endBlockFill(); }

//! Get the strictly upper triangular portion of the matrix
BlockedLinearOp getUpperTriBlocks(const BlockedLinearOp &blo, bool callEndBlockFill = true);

//! Get the strictly lower triangular portion of the matrix
BlockedLinearOp getLowerTriBlocks(const BlockedLinearOp &blo, bool callEndBlockFill = true);

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
 * \note The caller is responsible for calling
 *       <code>endBlockFill</code> on the returned blocked
 *       operator.
 */
BlockedLinearOp zeroBlockedOp(const BlockedLinearOp &blo);

//! Figure out if this operator is the zero operator (or null!)
bool isZeroOp(const LinearOp op);

/** \brief Compute absolute row sum matrix.
 *
 * Compute the absolute row sum matrix. That is
 * a diagonal operator composed of the absolute value of the
 * row sum.
 *
 * \returns A diagonal operator.
 */
ModifiableLinearOp getAbsRowSumMatrix(const LinearOp &op);

/** \brief Compute inverse of the absolute row sum matrix.
 *
 * Compute the inverse of the absolute row sum matrix. That is
 * a diagonal operator composed of the inverse of the absolute value
 * of the row sum.
 *
 * \returns A diagonal operator.
 */
ModifiableLinearOp getAbsRowSumInvMatrix(const LinearOp &op);

/** \brief Compute the lumped version of this matrix.
 *
 * Compute the lumped version of this matrix. That is
 * a diagonal operator composed of the row sum.
 *
 * \returns A diagonal operator.
 */
ModifiableLinearOp getLumpedMatrix(const LinearOp &op);

/** \brief Compute the inverse of the lumped version of
 *        this matrix.
 *
 * Compute the inverse of the lumped version of this matrix.
 * That is a diagonal operator composed of the row sum.
 *
 * \returns A diagonal operator.
 */
ModifiableLinearOp getInvLumpedMatrix(const LinearOp &op);

//@}

//! @name Mathematical functions
//@{

/** \brief Apply a linear operator to a multivector (think of this as a matrix
 *        vector multiply).
 *
 * Apply a linear operator to a multivector. This also permits arbitrary scaling
 * and addition of the result. This function gives
 *
 *    \f$ y = \alpha A x + \beta y \f$
 *
 * It is required that the range space of <code>A</code> is compatible with <code>y</code> and the
 * domain space of <code>A</code> is compatible with <code>x</code>.
 *
 * \param[in]     A
 * \param[in]     x
 * \param[in,out] y
 * \param[in]     alpha
 * \param[in]     beta
 *
 */
void applyOp(const LinearOp &A, const MultiVector &x, MultiVector &y, double alpha = 1.0,
             double beta = 0.0);

/** \brief Apply a transposed linear operator to a multivector (think of this as a matrix
 *        vector multiply).
 *
 * Apply a transposed linear operator to a multivector. This also permits arbitrary scaling
 * and addition of the result. This function gives
 *
 *    \f$ y = \alpha A^T x + \beta y \f$
 *
 * It is required that the domain space of <code>A</code> is compatible with <code>y</code> and the
 * range space of <code>A</code> is compatible with <code>x</code>.
 *
 * \param[in]     A
 * \param[in]     x
 * \param[in,out] y
 * \param[in]     alpha
 * \param[in]     beta
 *
 */
void applyTransposeOp(const LinearOp &A, const MultiVector &x, MultiVector &y, double alpha = 1.0,
                      double beta = 0.0);

/** \brief Apply a linear operator to a blocked multivector (think of this as a matrix
 *        vector multiply).
 *
 * Apply a linear operator to a blocked multivector. This also permits arbitrary scaling
 * and addition of the result. This function gives
 *
 *    \f$ y = \alpha A x + \beta y \f$
 *
 * It is required that the range space of <code>A</code> is compatible with <code>y</code> and the
 * domain space of <code>A</code> is compatible with <code>x</code>.
 *
 * \param[in]     A
 * \param[in]     x
 * \param[in,out] y
 * \param[in]     alpha
 * \param[in]     beta
 *
 */
inline void applyOp(const LinearOp &A, const BlockedMultiVector &x, BlockedMultiVector &y,
                    double alpha = 1.0, double beta = 0.0) {
  const MultiVector x_mv = toMultiVector(x);
  MultiVector y_mv       = toMultiVector(y);
  applyOp(A, x_mv, y_mv, alpha, beta);
}

/** \brief Apply a transposed linear operator to a blocked multivector (think of this as a matrix
 *        vector multiply).
 *
 * Apply a transposed linear operator to a blocked multivector. This also permits arbitrary scaling
 * and addition of the result. This function gives
 *
 *    \f$ y = \alpha A^T x + \beta y \f$
 *
 * It is required that the domain space of <code>A</code> is compatible with <code>y</code> and the
 * range space of <code>A</code> is compatible with <code>x</code>.
 *
 * \param[in]     A
 * \param[in]     x
 * \param[in,out] y
 * \param[in]     alpha
 * \param[in]     beta
 *
 */
inline void applyTransposeOp(const LinearOp &A, const BlockedMultiVector &x, BlockedMultiVector &y,
                             double alpha = 1.0, double beta = 0.0) {
  const MultiVector x_mv = toMultiVector(x);
  MultiVector y_mv       = toMultiVector(y);
  applyTransposeOp(A, x_mv, y_mv, alpha, beta);
}

/** \brief Update the <code>y</code> vector so that \f$y = \alpha x+\beta y\f$
 *
 * Compute the linear combination \f$y=\alpha x + \beta y\f$.
 *
 * \param[in]     alpha
 * \param[in]     x
 * \param[in]     beta
 * \param[in,out] y
 */
void update(double alpha, const MultiVector &x, double beta, MultiVector &y);

//! \brief Update for a BlockedMultiVector
inline void update(double alpha, const BlockedMultiVector &x, double beta, BlockedMultiVector &y) {
  MultiVector x_mv = toMultiVector(x);
  MultiVector y_mv = toMultiVector(y);
  update(alpha, x_mv, beta, y_mv);
}

//! Scale a multivector by a constant
inline void scale(const double alpha, MultiVector &x) { Thyra::scale<double>(alpha, x.ptr()); }

//! Scale a multivector by a constant
inline void scale(const double alpha, BlockedMultiVector &x) {
  MultiVector x_mv = toMultiVector(x);
  scale(alpha, x_mv);
}

//! Scale a modifiable linear op by a constant
inline LinearOp scale(const double alpha, ModifiableLinearOp &a) {
  return Thyra::nonconstScale(alpha, a);
}

//! Construct an implicit adjoint of the linear operators
inline LinearOp adjoint(ModifiableLinearOp &a) { return Thyra::nonconstAdjoint(a); }

//@}

//! \name Epetra_Operator specific functions
//@{

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
const ModifiableLinearOp getDiagonalOp(const LinearOp &op);

/** \brief Get the diagonal of a linear operator
 *
 * Get the diagonal of a linear operator, putting it
 * in the first column of a multivector.
 */
const MultiVector getDiagonal(const LinearOp &op);

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
const ModifiableLinearOp getInvDiagonalOp(const LinearOp &op);

/** \brief Multiply three linear operators.
 *
 * Multiply three linear operators. This currently assumes
 * that the underlying implementation uses Epetra_CrsMatrix.
 * The exception is that opm is allowed to be an diagonal matrix.
 *
 * \param[in] opl Left operator (assumed to be a Epetra_CrsMatrix)
 * \param[in] opm Middle operator (assumed to be a Epetra_CrsMatrix or a diagonal matrix)
 * \param[in] opr Right operator (assumed to be a Epetra_CrsMatrix)
 *
 * \returns Matrix product with a Epetra_CrsMatrix implementation
 */
const LinearOp explicitMultiply(const LinearOp &opl, const LinearOp &opm, const LinearOp &opr);

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
                                          const LinearOp &opr, const ModifiableLinearOp &destOp);

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
const LinearOp explicitMultiply(const LinearOp &opl, const LinearOp &opr);

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
                                          const ModifiableLinearOp &destOp);

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
const LinearOp explicitAdd(const LinearOp &opl, const LinearOp &opr);

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
const ModifiableLinearOp explicitAdd(const LinearOp &opl, const LinearOp &opr,
                                     const ModifiableLinearOp &destOp);

/** Sum into the modifiable linear op.
 */
const ModifiableLinearOp explicitSum(const LinearOp &opl, const ModifiableLinearOp &destOp);

/** Build an explicit transpose of a linear operator. (Concrete data
 * underneath.
 */
const LinearOp explicitTranspose(const LinearOp &op);

/** Explicitely scale a linear operator.
 */
const LinearOp explicitScale(double scalar, const LinearOp &op);

/** Rturn the frobenius norm of a linear operator
 */
double frobeniusNorm(const LinearOp &op);
double oneNorm(const LinearOp &op);
double infNorm(const LinearOp &op);

/** \brief Take the first column of a multivector and build a
 *        diagonal linear operator
 */
const LinearOp buildDiagonal(const MultiVector &v, const std::string &lbl = "ANYM");

/** \brief Using the first column of a multivector, take the elementwise build a
 *        inverse and build the inverse diagonal operator.
 */
const LinearOp buildInvDiagonal(const MultiVector &v, const std::string &lbl = "ANYM");

//@}

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
double computeSpectralRad(const Teuchos::RCP<const Thyra::LinearOpBase<double> > &A, double tol,
                          bool isHermitian = false, int numBlocks = 5, int restart = 0,
                          int verbosity = 0);

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
 * \return The smallest magnitude eigenvalue of the matrix.  If the algorithm didn't converge the
 *         number is the negative of the ritz-values. If a <code>NaN</code> is returned
 *         there was a problem constructing the Anasazi problem
 */
double computeSmallestMagEig(const Teuchos::RCP<const Thyra::LinearOpBase<double> > &A, double tol,
                             bool isHermitian = false, int numBlocks = 5, int restart = 0,
                             int verbosity = 0);

//! Type describing the type of diagonal to construct.
typedef enum {
  Diagonal  //! Specifies that just the diagonal is used
  ,
  Lumped  //! Specifies that row sum is used to form a diagonal
  ,
  AbsRowSum  //! Specifies that the \f$i^{th}\f$ diagonal entry is \f$\sum_j |A_{ij}|\f$
  ,
  BlkDiag  //! Specifies that a block diagonal approximation is to be used
  ,
  NotDiag  //! For user convenience, if Teko recieves this value, exceptions will be thrown
} DiagonalType;

/** Get a diagonal operator from a matrix. The mechanism for computing
 * the diagonal is specified by a <code>DiagonalType</code> arugment.
 *
 * \param[in] A <code>Epetra_CrsMatrix</code> to extract the diagonal from.
 * \param[in] dt Specifies the type of diagonal that is desired.
 *
 * \returns A diagonal operator.
 */
ModifiableLinearOp getDiagonalOp(const Teko::LinearOp &A, const DiagonalType &dt);

/** Get the inverse of a diagonal operator from a matrix. The mechanism for computing
 * the diagonal is specified by a <code>DiagonalType</code> arugment.
 *
 * \param[in] A <code>Epetra_CrsMatrix</code> to extract the diagonal from.
 * \param[in] dt Specifies the type of diagonal that is desired.
 *
 * \returns A inverse of a diagonal operator.
 */
ModifiableLinearOp getInvDiagonalOp(const Teko::LinearOp &A, const DiagonalType &dt);

/** \brief Get the diagonal of a sparse linear operator
 *
 * \param[in] Op Sparse linear operator to get diagonal of
 * \param[in] dt Type of diagonal operator required.
 */
const MultiVector getDiagonal(const LinearOp &op, const DiagonalType &dt);

/** Get a string corresponding to the type of digaonal specified.
 *
 * \param[in] dt The type of diagonal.
 *
 * \returns A string name representing this diagonal type.
 */
std::string getDiagonalName(const DiagonalType &dt);

/** Get a type corresponding to the name of a diagonal specified.
 *
 * \param[in] name String representing the diagonal type
 *
 * \returns The type representation of the string, if the
 *          string is not recognized this function returns
 *          a <code>NotDiag</code>
 */
DiagonalType getDiagonalType(std::string name);

#ifdef TEKO_HAVE_EPETRA
LinearOp probe(Teuchos::RCP<const Epetra_CrsGraph> &G, const LinearOp &Op);
#endif

/** Get the one norm of the vector
 */
double norm_1(const MultiVector &v, std::size_t col);

/** Get the two norm of the vector
 */
double norm_2(const MultiVector &v, std::size_t col);

/** This replaces entries of a vector falling below a particular
 * bound. Thus a an entry will be greater than or equal to \code{lowerBound}.
 */
void clipLower(MultiVector &v, double lowerBound);

/** This replaces entries of a vector above a particular
 * bound. Thus a an entry will be less than or equal to \code{upperBound}.
 */
void clipUpper(MultiVector &v, double upperBound);

/** This replaces entries of a vector equal to a particular value
 * with a new value.
 */
void replaceValue(MultiVector &v, double currentValue, double newValue);

/** Compute the averages of each column of the multivector.
 */
void columnAverages(const MultiVector &v, std::vector<double> &averages);

/** Compute the average of the solution.
 */
double average(const MultiVector &v);

/** Is this operator a physically blocked linear op?
 */
bool isPhysicallyBlockedLinearOp(const LinearOp &op);

/** Return a physically blocked linear op and whether it is scaled or transpose in its wrapper
 */
Teuchos::RCP<const Thyra::PhysicallyBlockedLinearOpBase<double> > getPhysicallyBlockedLinearOp(
    const LinearOp &op, ST *scalar, bool *transp);

//! Construct filename string for writing blocks to matrix-market format
std::string formatBlockName(const std::string &prefix, int i, int j, int nrow);

//! Write a matrix to file
void writeMatrix(const std::string &filename, const Teko::LinearOp &op);

}  // end namespace Teko

#ifdef _MSC_VER
#ifdef TEKO_DEFINED_MSC_EXTENSIONS
#undef _MSC_EXTENSIONS
#endif
#endif

#endif

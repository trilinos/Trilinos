/** \file PB_Utilities.hpp
  *
  * This file contains a number of useful functions and classes
  * used in PB. They are distinct from the core functionality of
  * the preconditioner factory, however, the functions are critical
  * to construction of the preconditioners themselves.
  */

#ifndef __PB_Utilities_hpp__
#define __PB_Utilities_hpp__

#include "Epetra_CrsMatrix.h"

// Thyra includes
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"

namespace PB {

using Thyra::multiply;
using Thyra::scale;
using Thyra::add;
using Thyra::identity;
using Thyra::zero; // make it to take one argument (square matrix)
using Thyra::block2x2;
using Thyra::block2x1;
using Thyra::block1x2;

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
Teuchos::RCP<Epetra_CrsMatrix> buildGraphLaplacian(int dim,double * coords,const Epetra_CrsMatrix & stencil);


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
Teuchos::RCP<Epetra_CrsMatrix> buildGraphLaplacian(double * x,double * y,double * z,int stride,const Epetra_CrsMatrix & stencil);

// typedefs for increased simplicity
typedef Teuchos::RCP<const Thyra::VectorSpaceBase<double> > VectorSpace;

// ----------------------------------------------------------------------------

//! @name MultiVector utilities
//@{

typedef Teuchos::RCP<Thyra::ProductMultiVectorBase<double> > BlockedMultiVector;
typedef Teuchos::RCP<Thyra::MultiVectorBase<double> > MultiVector;

//! Convert to a MultiVector from a BlockedMultiVector
inline MultiVector toMultiVector(BlockedMultiVector & bmv) { return bmv; }

//! Convert to a MultiVector from a BlockedMultiVector
inline const MultiVector toMultiVector(const BlockedMultiVector & bmv) { return bmv; }

//! Convert to a BlockedMultiVector from a MultiVector
inline const BlockedMultiVector toBlockedMultiVector(const MultiVector & bmv) 
{ return Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<double> >(bmv); }

//! Get the column count in a block linear operator
inline int blockCount(const BlockedMultiVector & bmv)
{ return bmv->productSpace()->numBlocks(); }

//! Get the <code>i</code>th block from a BlockedMultiVector object
inline MultiVector getBlock(int i,const BlockedMultiVector & bmv)
{ return Teuchos::rcp_const_cast<Thyra::MultiVectorBase<double> >(bmv->getMultiVectorBlock(i)); }

//! Performa deep copy of the vector
inline MultiVector deepcopy(const MultiVector & v)
{ return v->clone_mv(); }

//! Performa deep copy of the blocked vector
inline BlockedMultiVector deepcopy(const BlockedMultiVector & v)
{ return toBlockedMultiVector(v->clone_mv()); }

//! build a BlockedMultiVector from a vector of MultiVectors
BlockedMultiVector buildBlockedMultiVector(const std::vector<MultiVector> & mvs);

//@}

// ----------------------------------------------------------------------------

//! @name LinearOp utilities
//@{
typedef Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > BlockedLinearOp;
typedef Teuchos::RCP<const Thyra::LinearOpBase<double> > LinearOp;
typedef Teuchos::RCP<Thyra::LinearOpBase<double> > InverseLinearOp;

//! Build a square zero operator from a single vector space
inline LinearOp zero(const VectorSpace & vs)
{ return Thyra::zero<double>(vs,vs); }

//! Get the range space of a linear operator
inline VectorSpace rangeSpace(const LinearOp & lo)
{ return lo->range(); }

//! Get the domain space of a linear operator
inline VectorSpace domainSpace(const LinearOp & lo)
{ return lo->domain(); }

//! Converts a LinearOp to a BlockedLinearOp
inline BlockedLinearOp toBlockedLinearOp(LinearOp & clo)
{
   Teuchos::RCP<Thyra::LinearOpBase<double> > lo = Teuchos::rcp_const_cast<Thyra::LinearOpBase<double> >(clo);
   return Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<double> > (lo);
}

//! Converts a LinearOp to a BlockedLinearOp
inline const BlockedLinearOp toBlockedLinearOp(const LinearOp & clo)
{
   Teuchos::RCP<Thyra::LinearOpBase<double> > lo = Teuchos::rcp_const_cast<Thyra::LinearOpBase<double> >(clo);
   return Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<double> > (lo);
}

//! Convert to a LinearOp from a BlockedLinearOp
inline LinearOp toLinearOp(BlockedLinearOp & blo) { return blo; }

//! Convert to a LinearOp from a BlockedLinearOp
inline const LinearOp toLinearOp(const BlockedLinearOp & blo) { return blo; }

//! Get the row count in a block linear operator
inline int blockRowCount(const BlockedLinearOp & blo)
{ return blo->productRange()->numBlocks(); }

//! Get the column count in a block linear operator
inline int blockColCount(const BlockedLinearOp & blo)
{ return blo->productDomain()->numBlocks(); }

//! Get the <code>i,j</code> block in a BlockedLinearOp object
inline LinearOp getBlock(int i,int j,const BlockedLinearOp & blo)
{ return blo->getBlock(i,j); }

//! Set the <code>i,j</code> block in a BlockedLinearOp object
inline void setBlock(int i,int j,BlockedLinearOp & blo, const LinearOp & lo)
{ return blo->setBlock(i,j,lo); }

//! Build a new blocked linear operator
inline BlockedLinearOp createBlockedOp()
{ return rcp(new Thyra::DefaultBlockedLinearOp<double>()); }

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
inline void beginBlockFill(BlockedLinearOp & blo,int rowCnt,int colCnt)
{ blo->beginBlockFill(rowCnt,colCnt); }

/** \brief Let the blocked operator know that you are going to 
  *        set the sub blocks.
  *
  * Let the blocked operator know that you are going to 
  * set the sub blocks. This is a simple wrapper around the
  * member function of the same name in Thyra.
  *
  * \param[in,out] blo Blocked operator to have its fill stage activated
  */
inline void beginBlockFill(BlockedLinearOp & blo)
{ blo->beginBlockFill(); }

//! Notify the blocked operator that the fill stage is completed.
inline void endBlockFill(BlockedLinearOp & blo)
{ blo->endBlockFill(); }

//! Get the strictly upper triangular portion of the matrix
BlockedLinearOp getUpperTriBlocks(const BlockedLinearOp & blo);

//! Get the strictly lower triangular portion of the matrix
BlockedLinearOp getLowerTriBlocks(const BlockedLinearOp & blo);

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
BlockedLinearOp zeroBlockedOp(const BlockedLinearOp & blo);

//! Figure out if this operator is the zero operator (or null!)
bool isZeroOp(const LinearOp op);

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
  * It is required that the range space of <code>A</code> is compatible with <code>y</code> and the domain space
  * of <code>A</code> is compatible with <code>x</code>.
  *
  * \param[in]     A
  * \param[in]     x
  * \param[in,out] y
  * \param[in]     alpha
  * \param[in]     beta
  *
  */
void applyOp(const LinearOp & A,const MultiVector & x,MultiVector & y,double alpha=1.0,double beta=0.0);

/** \brief Apply a linear operator to a blocked multivector (think of this as a matrix
  *        vector multiply).
  *
  * Apply a linear operator to a blocked multivector. This also permits arbitrary scaling
  * and addition of the result. This function gives
  *     
  *    \f$ y = \alpha A x + \beta y \f$
  *
  * It is required that the range space of <code>A</code> is compatible with <code>y</code> and the domain space
  * of <code>A</code> is compatible with <code>x</code>.
  *
  * \param[in]     A
  * \param[in]     x
  * \param[in,out] y
  * \param[in]     alpha
  * \param[in]     beta
  *
  */
inline void applyOp(const LinearOp & A,const BlockedMultiVector & x,BlockedMultiVector & y,double alpha=1.0,double beta=0.0)
{ const MultiVector x_mv = toMultiVector(x); MultiVector y_mv = toMultiVector(y);
  applyOp(A,x_mv,y_mv,alpha,beta); }

/** \brief Update the <code>y</code> vector so that \f$y = \alpha x+\beta y\f$
  *
  * Compute the linear combination \f$y=\alpha x + \beta y\f$. 
  *
  * \param[in]     alpha
  * \param[in]     x 
  * \param[in]     beta 
  * \param[in,out] y
  */
void update(double alpha,const MultiVector & x,double beta,MultiVector & y);

//! \brief Update for a BlockedMultiVector
inline void update(double alpha,const BlockedMultiVector & x,double beta,BlockedMultiVector & y)
{ MultiVector x_mv = toMultiVector(x); MultiVector y_mv = toMultiVector(y);
  update(alpha,x_mv,beta,y_mv); }

//! Scale a multivector by a constant
inline void scale(const double alpha,MultiVector & x) { Thyra::scale<double>(alpha,x.ptr()); }

//! Scale a multivector by a constant
inline void scale(const double alpha,BlockedMultiVector & x) 
{  MultiVector x_mv = toMultiVector(x); scale(alpha,x_mv); }

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
const LinearOp getDiagonalOp(const LinearOp & op);

/** \brief Get the diagonal of a linear operator
  *
  * Get the diagonal of a linear operator, putting it
  * in the first column of a multivector.
  */
const MultiVector getDiagonal(const LinearOp & op);

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
const LinearOp getInvDiagonalOp(const LinearOp & op);

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
const LinearOp explicitMultiply(const LinearOp & opl,const LinearOp & opm,const LinearOp & opr);

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
const LinearOp explicitMultiply(const LinearOp & opl,const LinearOp & opr);

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
const LinearOp explicitAdd(const LinearOp & opl,const LinearOp & opr);

/** \brief Take the first column of a multivector and build a
  *        diagonal linear operator
  */
const LinearOp buildDiagonal(const MultiVector & v,const std::string & lbl="ANYM");

/** \brief Using the first column of a multivector, take the elementwise build a
  *        inverse and build the inverse diagonal operator.
  */
const LinearOp buildInvDiagonal(const MultiVector & v,const std::string & lbl="ANYM");

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
double computeSpectralRad(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A,double tol,
                          bool isHermitian=false,int numBlocks=5,int restart=0,int verbosity=0);

} // end namespace PB

#endif

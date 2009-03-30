#ifndef __PB_Utilities_hpp__
#define __PB_Utilities_hpp__

#include "Epetra_CrsMatrix.h"

// Thyra includes
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

namespace PB {

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
  * \param[in,out] gl      The graph Laplacian matrix to be filled according
  *                        to the <code>stencil</code> matrix.
  *
  * \pre Assumes the <code>gl</code> argument is constructed to have a row map
  *      equivalent in size to <code>stencil</code>
  */
// void buildGraphLaplacian(int dim,double * x,double * y,double * z,const Epetra_CrsMatrix & stencil,Epetra_CrsMatrix & gl);
void buildGraphLaplacian(int dim,double * coords,const Epetra_CrsMatrix & stencil,Epetra_CrsMatrix & gl);

// Function that might be useful
// multiply
// add
// buildInverse(InverseFactory,LinearOp) 


// typedefs for increased simplicity
typedef Teuchos::RCP<const Thyra::VectorSpaceBase<double> > VectorSpace;

//! @name MultiVector utilities
//@{

typedef Teuchos::RCP<Thyra::ProductMultiVectorBase<double> > BlockedMultiVector;
typedef Teuchos::RCP<Thyra::MultiVectorBase<double> > MultiVector;

//! Convert to a MultiVector from a BlockedMultiVector
inline MultiVector toMultiVector(BlockedMultiVector & bmv) { return bmv; }

//! Convert to a MultiVector from a BlockedMultiVector
inline const MultiVector toMultiVector(const BlockedMultiVector & bmv) { return bmv; }

//! Get the <code>i</code>th block from a BlockedMultiVector object
inline MultiVector getBlock(int i,const BlockedMultiVector & bmv)
{ return Teuchos::rcp_const_cast<Thyra::MultiVectorBase<double> >(bmv->getMultiVectorBlock(i)); }

inline MultiVector deepcopy(const MultiVector & v)
{ return v->clone_mv(); }

//@}

//! @name LinearOp utilities
//@{
typedef Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > BlockedLinearOp;
typedef Teuchos::RCP<const Thyra::LinearOpBase<double> > LinearOp;

//! Convert to a LinearOp from a BlockedLinearOp
inline LinearOp toLinearOp(BlockedLinearOp & bmv) { return bmv; }

//! Convert to a LinearOp from a BlockedLinearOp
inline const LinearOp toLinearOp(const BlockedLinearOp & bmv) { return bmv; }

//! Get the <code>i,j</code> block in a BlockedLinearOp object
inline LinearOp getBlock(int i,int j,const BlockedLinearOp & blo)
{ return blo->getBlock(i,j); }

//! Set the <code>i,j</code> block in a BlockedLinearOp object
inline void setBlock(int i,int j,BlockedLinearOp & blo, const LinearOp & lo)
{ return blo->setBlock(i,j,lo); }

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
  * \param[in]     A
  * \param[in]     x
  * \param[in,out] y
  * \param[in]     alpha
  * \param[in]     beta
  *
  */
void applyOp(const LinearOp & A,const MultiVector & x,MultiVector & y,double alpha=1.0,double beta=0.0);

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

//! Scale a multivector by a constant
inline void scale(const double alpha,MultiVector & x) { Thyra::scale<double>(alpha,x.ptr()); }

//! Scale a multivector by a constant
inline void scale(const double alpha,BlockedMultiVector & x) { scale(alpha,x); }

//@}

} // end namespace PB

#endif

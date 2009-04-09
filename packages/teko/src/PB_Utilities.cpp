#include "PB_Utilities.hpp"

#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_ZeroLinearOpBase.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_EpetraExtDiagScaledMatProdTransformer.hpp"
#include "Thyra_EpetraExtAddTransformer.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"

#include "Teuchos_Array.hpp"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"

// PB includes
#include "Epetra/PB_EpetraHelpers.hpp"

#include <cmath>

namespace PB {

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;

// distance function...not parallel...entirely internal to this cpp file
inline double dist(int dim,double * coords,int row,int col)
{
   double value = 0.0;
   for(int i=0;i<dim;i++)
      value += std::pow(coords[dim*row+i]-coords[dim*col+i],2.0);

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
  * \param[in,out] gl      The graph Laplacian matrix to be filled according
  *                        to the <code>stencil</code> matrix.
  *
  * \pre Assumes the <code>gl</code> argument is constructed to have a row map
  *      equivalent in size to <code>stencil</code>
  */
void buildGraphLaplacian(int dim,double * coords,const Epetra_CrsMatrix & stencil,Epetra_CrsMatrix & gl)
{
   // allocate an additional value for the diagonal, if neccessary
   double rowData[stencil.GlobalMaxNumEntries()+1];
   int rowInd[stencil.GlobalMaxNumEntries()+1];

   // loop over all the rows
   for(int j=0;j<gl.NumMyRows();j++) {
      int row = gl.GRID(j);
      double diagValue = 0.0;
      int diagInd = -1;
      int rowSz = 0;

      // extract a copy of this row...put it in rowData, rowIndicies
      stencil.ExtractGlobalRowCopy(row,stencil.MaxNumEntries(),rowSz,rowData,rowInd);
 
      // loop over elements of row
      for(int i=0;i<rowSz;i++) {
         int col = rowInd[i];

         // for nondiagonal entries
         if(row!=col) {
            double d = dist(dim,coords,row,col);
            rowData[i] = -1.0/d;
            diagValue += rowData[i];
         }
         else 
            diagInd = i;
      }
    
      // handle diagonal entry
      if(diagInd<0) { // diagonal not in row
         rowData[rowSz] = -diagValue;
         rowInd[rowSz] = row;
         rowSz++;
      }
      else { // diagonal in row
         rowData[diagInd] = -diagValue;
         rowInd[diagInd] = row;
      }

      // insert row data into graph Laplacian matrix
      gl.InsertGlobalValues(row,rowSz,rowData,rowInd);
   }

   gl.FillComplete();
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
void applyOp(const LinearOp & A,const MultiVector & x,MultiVector & y,double alpha,double beta)
{
   Thyra::apply(*A,Thyra::NONCONJ_ELE,*x,&*y,alpha,beta);
}

/** \brief Update the <code>y</code> vector so that \f$y = \alpha x+\beta y\f$
  */
void update(double alpha,const MultiVector & x,double beta,MultiVector & y)
{
   Teuchos::Array<double> scale;
   Teuchos::Array<Teuchos::Ptr<const Thyra::MultiVectorBase<double> > >  vec;

   // build arrays needed for linear combo
   scale.push_back(alpha);
   vec.push_back(x.ptr());

   // compute linear combination
   Thyra::linear_combination<double>(scale,vec,beta,y.ptr());
}

//! Get the strictly upper triangular portion of the matrix
BlockedLinearOp getUpperTriBlocks(const BlockedLinearOp & blo)
{
   int rows = blockRowCount(blo);

   TEUCHOS_ASSERT(rows==blockColCount(blo));

   RCP<const Thyra::ProductVectorSpaceBase<double> > range = blo->productRange();
   RCP<const Thyra::ProductVectorSpaceBase<double> > domain = blo->productDomain();

   // allocate new operator
   BlockedLinearOp upper = createNewBlockedOp();
 
   // build new operator 
   upper->beginBlockFill(rows,rows);

   for(int i=0;i<rows;i++) {
      // put zero operators on the diagonal
      // this gurantees the vector space of
      // the new operator are fully defined
      RCP<const Thyra::LinearOpBase<double> > zed 
            = Thyra::zero<double>(range->getBlock(i),domain->getBlock(i));
      upper->setBlock(i,i,zed);

      for(int j=i+1;j<rows;j++) {
         // get block i,j
         LinearOp uij = blo->getBlock(i,j);

         // stuff it in U
         if(uij!=Teuchos::null)
            upper->setBlock(i,j,uij);
      }
   }
   upper->endBlockFill();

   return upper;
}

//! Get the strictly lower triangular portion of the matrix
BlockedLinearOp getLowerTriBlocks(const BlockedLinearOp & blo)
{
   int rows = blockRowCount(blo);

   TEUCHOS_ASSERT(rows==blockColCount(blo));

   RCP<const Thyra::ProductVectorSpaceBase<double> > range = blo->productRange();
   RCP<const Thyra::ProductVectorSpaceBase<double> > domain = blo->productDomain();

   // allocate new operator
   BlockedLinearOp lower = createNewBlockedOp();
 
   // build new operator 
   lower->beginBlockFill(rows,rows);

   for(int i=0;i<rows;i++) {
      // put zero operators on the diagonal
      // this gurantees the vector space of
      // the new operator are fully defined
      RCP<const Thyra::LinearOpBase<double> > zed 
            = Thyra::zero<double>(range->getBlock(i),domain->getBlock(i));
      lower->setBlock(i,i,zed);

      for(int j=0;j<i;j++) {
         // get block i,j
         LinearOp uij = blo->getBlock(i,j);

         // stuff it in U
         if(uij!=Teuchos::null)
            lower->setBlock(i,j,uij);
      }
   }
   lower->endBlockFill();

   return lower;
}

//! Figure out if this operator is the zero operator (or null!)
bool isZeroOp(const LinearOp op)
{
   // if operator is null...then its zero!
   if(op==Teuchos::null) return true;

   // try to cast it to a zero linear operator
   const LinearOp test = rcp_dynamic_cast<const Thyra::ZeroLinearOpBase<double> >(op);

   // if it works...then its zero...otherwise its null
   return test!=Teuchos::null;
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
const LinearOp getDiagonalOp(const LinearOp & op)
{
   // get Epetra_CrsMatrix
   RCP<const Epetra_CrsMatrix> eOp = rcp_dynamic_cast<const Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*op));
   TEST_FOR_EXCEPTION(eOp==Teuchos::null,std::runtime_error,"getDiagonalOp requires an Epetra_CrsMatrix");

   // extract diagonal
   const RCP<Epetra_Vector> diag = rcp(new Epetra_Vector(eOp->OperatorRangeMap()));
   eOp->ExtractDiagonalCopy(*diag);

   // build Thyra diagonal operator
   return PB::Epetra::thyraDiagOp(diag,eOp->OperatorRangeMap());
}

const MultiVector getDiagonal(const LinearOp & op)
{
   // get Epetra_CrsMatrix
   RCP<const Epetra_CrsMatrix> eOp = rcp_dynamic_cast<const Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*op));
   TEST_FOR_EXCEPTION(eOp==Teuchos::null,std::runtime_error,"getDiagonalOp requires an Epetra_CrsMatrix");

   // extract diagonal
   const RCP<Epetra_Vector> diag = rcp(new Epetra_Vector(eOp->RowMap()));
   eOp->ExtractDiagonalCopy(*diag);

   return Thyra::create_Vector(diag,Thyra::create_VectorSpace(Teuchos::rcpFromRef(eOp->RowMap())));
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
const LinearOp getInvDiagonalOp(const LinearOp & op)
{
   // get Epetra_CrsMatrix
   RCP<const Epetra_CrsMatrix> eOp = rcp_dynamic_cast<const Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*op));
   TEST_FOR_EXCEPTION(eOp==Teuchos::null,std::runtime_error,"getDiagonalOp requires an Epetra_CrsMatrix");

   // extract diagonal
   const RCP<Epetra_Vector> diag = rcp(new Epetra_Vector(eOp->RowMap()));
   eOp->ExtractDiagonalCopy(*diag);
   diag->Reciprocal(*diag);

   // build Thyra diagonal operator
   return PB::Epetra::thyraDiagOp(diag,eOp->RowMap());
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
const LinearOp explicitMultiply(const LinearOp & opl,const LinearOp & opm,const LinearOp & opr)
{
   // build implicit multiply
   const LinearOp implicitOp = Thyra::multiply(opl,opm,opr);

   // build transformer
   const RCP<Thyra::LinearOpTransformerBase<double> > prodTrans =
       Thyra::epetraExtDiagScaledMatProdTransformer();

   // build operator and multiply
   const RCP<Thyra::LinearOpBase<double> > explicitOp = prodTrans->createOutputOp();
   prodTrans->transform(*implicitOp,explicitOp.ptr());

   return explicitOp;
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
const LinearOp explicitMultiply(const LinearOp & opl,const LinearOp & opr)
{
   // build implicit multiply
   const LinearOp implicitOp = Thyra::multiply(opl,opr);
 
   // build transformer
   const RCP<Thyra::LinearOpTransformerBase<double> > prodTrans =
       Thyra::epetraExtDiagScaledMatProdTransformer();

   // build operator and multiply
   const RCP<Thyra::LinearOpBase<double> > explicitOp = prodTrans->createOutputOp();
   prodTrans->transform(*implicitOp,explicitOp.ptr());

   return explicitOp;
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
const LinearOp explicitAdd(const LinearOp & opl,const LinearOp & opr)
{
   // build implicit multiply
   const LinearOp implicitOp = Thyra::add(opl,opr);
 
   // build transformer
   const RCP<Thyra::LinearOpTransformerBase<double> > prodTrans =
       Thyra::epetraExtAddTransformer();

   // build operator and multiply
   const RCP<Thyra::LinearOpBase<double> > explicitOp = prodTrans->createOutputOp();
   prodTrans->transform(*implicitOp,explicitOp.ptr());

   return explicitOp;
}

const LinearOp buildDiagonal(const MultiVector & v)
{
   return Thyra::diagonal<double>(v->col(0));
}

const LinearOp buildInvDiagonal(const MultiVector & src)
{
   MultiVector dst = deepcopy(src); 
   Thyra::reciprocal<double>(dst->col(0).ptr(),*src->col(0));

   return Thyra::diagonal<double>(dst->col(0));
}

//! build a BlockedMultiVector from a vector of MultiVectors
BlockedMultiVector buildBlockedMultiVector(const std::vector<MultiVector> & mvv)
{
   Teuchos::Array<MultiVector> mvA;
   Teuchos::Array<VectorSpace> vsA;

   // build arrays of multi vectors and vector spaces
   std::vector<MultiVector>::const_iterator itr;
   for(itr=mvv.begin();itr!=mvv.end();++itr) {
      mvA.push_back(*itr);
      vsA.push_back((*itr)->range());
   }

   // construct the product vector space
   const RCP<const Thyra::DefaultProductVectorSpace<double> > vs
         = Thyra::productVectorSpace<double>(vsA);

   return Thyra::defaultProductMultiVector<double>(vs,mvA);
}

}

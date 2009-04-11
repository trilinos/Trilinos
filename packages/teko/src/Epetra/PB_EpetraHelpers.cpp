#include "PB_EpetraHelpers.hpp"

// Thyra Includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"

// Epetra includes
#include "Epetra_Vector.h"

// EpetraExt includes
#include "EpetraExt_ProductOperator.h"
#include "EpetraExt_MatrixMatrix.h"

// PB includes
#include "PB_EpetraOperatorWrapper.hpp"
#include "PB_Helpers.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::rcp_dynamic_cast;
using Teuchos::null;

namespace PB {
namespace Epetra {

// build an epetra operator from a block 2x2 matrix
Epetra_Operator * block2x2(const Epetra_Operator * A,const Epetra_Operator * B,
                         const Epetra_Operator * C,const Epetra_Operator * D,const std::string & str)
{
   using Teuchos::rcpFromRef;
   using Thyra::zero;

   // not all things can be zero
   assert(not (A==0 && B==0));
   assert(not (A==0 && C==0));
   assert(not (B==0 && D==0));
   assert(not (C==0 && D==0));


   RCP<const Thyra::LinearOpBase<double> > ptrA = A!=0 ? Thyra::epetraLinearOp(rcp(A,false),str+"_00"): null;
   RCP<const Thyra::LinearOpBase<double> > ptrB = B!=0 ? Thyra::epetraLinearOp(rcp(B,false),str+"_01"): null;
   RCP<const Thyra::LinearOpBase<double> > ptrC = C!=0 ? Thyra::epetraLinearOp(rcp(C,false),str+"_10"): null;
   RCP<const Thyra::LinearOpBase<double> > ptrD = D!=0 ? Thyra::epetraLinearOp(rcp(D,false),str+"_11"): null;

   // get some ranges so the null matrices can be set to zero
   const RCP<const Thyra::VectorSpaceBase<double> > row0Range  = A!=0 ? ptrA->range()  : ptrB->range();
   const RCP<const Thyra::VectorSpaceBase<double> > row1Range  = C!=0 ? ptrC->range()  : ptrD->range();
   const RCP<const Thyra::VectorSpaceBase<double> > col0Domain = A!=0 ? ptrA->domain() : ptrC->domain();
   const RCP<const Thyra::VectorSpaceBase<double> > col1Domain = B!=0 ? ptrB->domain() : ptrD->domain();

   // set previously null pointers to 0
   if(A==0) ptrA = Thyra::zero(row0Range,col0Domain);
   if(B==0) ptrB = Thyra::zero(row0Range,col1Domain);
   if(C==0) ptrC = Thyra::zero(row1Range,col0Domain);
   if(D==0) ptrD = Thyra::zero(row1Range,col1Domain);

   const RCP<const Thyra::LinearOpBase<double> > mat = Thyra::block2x2(ptrA,ptrB,ptrC,ptrD,str);

   return new PB::Epetra::EpetraOperatorWrapper(mat);
}

// Convert an Epetra "inverse" operator to a "forward" operator
Epetra_Operator * mechanicalInverse(const Epetra_Operator * inverse)
{
   // convert to a forward operator
   Teuchos::RCP<const Epetra_Operator> opAr[] = {rcp(inverse,false)};
   Teuchos::ETransp epetraOpsTransp[] = { Teuchos::NO_TRANS };
   EpetraExt::ProductOperator::EApplyMode epetraOpsApplyMode[] = { EpetraExt::ProductOperator::APPLY_MODE_APPLY_INVERSE };

   return new EpetraExt::ProductOperator(1,opAr,epetraOpsTransp,epetraOpsApplyMode);
}

const Teuchos::RCP<const Thyra::LinearOpBase<double> > thyraDiagOp(const RCP<const Epetra_Vector> & ev,const Epetra_Map & map,
                                                                   const std::string & lbl)
{
   const RCP<const Thyra::VectorBase<double> > thyraVec  // need a Thyra::VectorBase object
         = Thyra::create_Vector(ev,Thyra::create_VectorSpace(rcpFromRef(map)));
   Teuchos::RCP<Thyra::LinearOpBase<double> > op 
         = Teuchos::rcp(new Thyra::DefaultDiagonalLinearOp<double>(thyraVec));
   op->setObjectLabel(lbl);
   return op;
}

} // end namespace Epetra
} // end namespace PB

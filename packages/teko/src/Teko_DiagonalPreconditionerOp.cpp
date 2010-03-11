#include "Teko_DiagonalPreconditionerOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "EpetraExt_PointToBlockDiagPermute.h"
#include "Epetra_MultiVector.h"

using Teuchos::rcpFromRef;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcp_const_cast;
using Teuchos::RCP;

using Thyra::MultiVectorBase;

namespace Teko {

DiagonalPreconditionerOp::DiagonalPreconditionerOp(Teuchos::RCP<EpetraExt_PointToBlockDiagPermute> BDP, const VectorSpace range, const VectorSpace domain):
  BDP_(BDP),
  range_(range),
  domain_(domain){}


void DiagonalPreconditionerOp::implicitApply(const MultiVector & x, MultiVector & y,
					      const double alpha, const double beta) const{

  // Get the Multivectors into Epetra land
  // NTS: Thyra inexplicably wants maps, even when they are completely unecessary.
  const Epetra_Map & rangemap_=BDP_->OperatorRangeMap();
  const Epetra_Map & domainmap_=BDP_->OperatorDomainMap();

  RCP<const Epetra_MultiVector> x_=Thyra::get_Epetra_MultiVector(domainmap_,x);
  RCP<Epetra_MultiVector> y_=Thyra::get_Epetra_MultiVector(rangemap_,y); 
  TEUCHOS_ASSERT(x_!=Teuchos::null);
  TEUCHOS_ASSERT(y_!=Teuchos::null);


  // y = \alpha M x + \beta y $
  if(beta==0.0){
    BDP_->ApplyInverse(*x_,*y_);
    scale(alpha,y);
  }
  else{
    MultiVector y0=deepcopy(y);
    BDP_->ApplyInverse(*x_,*y_);
    update(alpha,y,beta,y0);  
  }
}

bool DiagonalPreconditionerOp::opSupportedImpl(const Thyra::EOpTransp M_trans) const
{
  return (M_trans == Thyra::NOTRANS);
}

void DiagonalPreconditionerOp::applyImpl(
  const Thyra::EOpTransp M_trans,
  const Thyra::MultiVectorBase<double> & x,
  const Teuchos::Ptr<Thyra::MultiVectorBase<double> > & y,
  const double alpha,
  const double beta
  ) const
{
   TEST_FOR_EXCEPTION(M_trans!=Thyra::NOTRANS, std::runtime_error,
     "Linear operators of inherited type ImplicitLinearOp "
     "cannot handle conjugation (yet!)");

   // cast source vector
   RCP<const MultiVectorBase<double> > src =
     rcp_dynamic_cast<const MultiVectorBase<double> >(rcpFromRef(x));
   MultiVector srcX = rcp_const_cast<MultiVectorBase<double> >(src);

   // cast destination vector
   MultiVector destY =
     Teuchos::rcp_dynamic_cast<MultiVectorBase<double> >(rcpFromPtr(y));

   // call apply
   implicitApply(srcX,destY,alpha,beta);
}



 void DiagonalPreconditionerOp::describe(Teuchos::FancyOStream & out_arg,
                                      const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::OSTab;

  switch(verbLevel) {
     case Teuchos::VERB_DEFAULT:
     case Teuchos::VERB_LOW:
        out_arg << this->description() << std::endl;
        break;
     case Teuchos::VERB_MEDIUM:
     case Teuchos::VERB_HIGH:
     case Teuchos::VERB_EXTREME:
       if(BDP_!=Teuchos::null)
	 BDP_->Print(out_arg);
       break;
     default:
        TEST_FOR_EXCEPT(true); // Should never get here!
  }
}

} // end namespace Teko


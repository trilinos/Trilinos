#include "Teko_ReorderedLinearOp.hpp"

namespace Teko {

ReorderedLinearOp::ReorderedLinearOp(const Teuchos::RCP<const BlockReorderManager> & mgr,
                                     const Teuchos::RCP<Thyra::LinearOpBase<double> > & blockedOp)
   : mgr_(mgr), blockedOp_(blockedOp)
{
   
   range_ = buildFlatVectorSpace(*mgr_,blockedOp_->range());
   domain_ = buildFlatVectorSpace(*mgr_,blockedOp_->domain());
}

VectorSpace ReorderedLinearOp::range() const
{
   return range_;
}

VectorSpace ReorderedLinearOp::domain() const
{
   return domain_;
}

void ReorderedLinearOp::implicitApply(const MultiVector & x, MultiVector & y,
                                      const double alpha, const double beta) const
{
   using Teuchos::rcp_dynamic_cast;

   Teuchos::RCP<const Thyra::MultiVectorBase<double> > reorderX 
      = Teko::buildReorderedMultiVector(*mgr_,rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<double> >(x));
   MultiVector reorderY = Teko::buildReorderedMultiVector(*mgr_,rcp_dynamic_cast<Thyra::ProductMultiVectorBase<double> >(y));

   // this will automatically fill the right data
   Thyra::apply(*blockedOp_,Thyra::NOTRANS,*reorderX,reorderY.ptr(),alpha,beta);
}

void ReorderedLinearOp::describe(Teuchos::FancyOStream & out_arg,
                                 const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::RCP;
  using Teuchos::OSTab;

  RCP<Teuchos::FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  switch(verbLevel) {
     case Teuchos::VERB_DEFAULT:
     case Teuchos::VERB_LOW:
        *out << this->description() << std::endl;
        break;
     case Teuchos::VERB_MEDIUM:
     case Teuchos::VERB_HIGH:
     case Teuchos::VERB_EXTREME:
        {
           *out << Teuchos::Describable::description() << "{"
                << "rangeDim=" << this->range()->dim()
                << ",domainDim=" << this->domain()->dim()
                << "}\n";
           {
              OSTab tab(out);
              *out << "[Blocked Op] = ";
              *out << Teuchos::describe(*blockedOp_,verbLevel);
           }
           {
              OSTab tab(out);
              *out << "[Blocked Manager] = ";
              *out << mgr_->toString() << std::endl;
           }
           break;
        }
     default:
        TEUCHOS_TEST_FOR_EXCEPT(true); // Should never get here!
  }
}

} // end namespace Teko

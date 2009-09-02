#include "Thyra_ModelEvaluatorBase_decl.hpp"


namespace Thyra {


std::string
ModelEvaluatorBase::DerivativeSupport::description() const
{
  std::ostringstream oss;
  oss << "DerivativeSupport{";
  if (none()) {
    oss << "none";
  }
  else {
    bool wroteOutput = false;
    if (supportsLinearOp_) {
      oss << "DERIV_LINEAR_OP";
      wroteOutput = true;
    }
    if (supportsMVByCol_) {
      oss << (wroteOutput?",":"") << toString(DERIV_MV_BY_COL);
      wroteOutput = true;
    }
    if (supportsTransMVByRow_) {
      oss << (wroteOutput?",":"") << toString(DERIV_TRANS_MV_BY_ROW);
      wroteOutput = true;
    }
  }
  oss << "}";
  return oss.str();
}


} // namespace Thyra


#ifdef HAVE_THYRA_EXPLICIT_INSTANTIATION


#include "Thyra_ModelEvaluatorBase_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Thyra {

TEUCHOS_MACRO_TEMPLATE_INSTANT_SCALAR_TYPES(THYRA_MODEL_EVALUATOR_BASE_INSTANT)

} // namespace Thyra

#endif // HAVE_THYRA_EXPLICIT_INSTANTIATION

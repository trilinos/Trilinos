#include "Tempus_ResidualModelEvaluator.hpp"
#include "Tempus_ResidualModelEvaluator_impl.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION

namespace Tempus {

template class ResidualModelEvaluator<double>;

template <typename double>
ResidualModelEvaluator<double>::
ResidualModelEvaluator(
  const Teuchos::RCP<Thyra::ModelEvaluator<double> >& transientModel);

} // namespace Tempus

#endif

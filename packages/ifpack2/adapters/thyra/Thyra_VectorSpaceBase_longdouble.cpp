#include "Teuchos_ConfigDefs.hpp"
#include "Thyra_ConfigDefs.hpp"

#if defined(HAVE_THYRA_EXPLICIT_INSTANTIATION) && defined(HAVE_TEUCHOS_LONG_DOUBLE)

#include "Thyra_VectorSpaceBase_decl.hpp"
#include "Thyra_VectorSpaceBase_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Thyra {

  // mfh 12 Sep 2015: Leaving this in results in duplicate
  // instantiation of VectorSpaceBase.
  //
  //TEUCHOS_CLASS_TEMPLATE_INSTANT_LONG_DOUBLE(VectorSpaceBase)

TEUCHOS_MACRO_TEMPLATE_INSTANT_LONG_DOUBLE(THYRA_VECTOR_SPACE_BASE_INSTANT)

} // namespace Thyra

#endif // HAVE_THYRA_EXPLICIT_INSTANTIATION && HAVE_TEUCHOS_LONG_DOUBLE

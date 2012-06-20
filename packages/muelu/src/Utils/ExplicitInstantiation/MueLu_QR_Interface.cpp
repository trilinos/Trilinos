#include "MueLu_ExplicitInstantiation.hpp"

#include "MueLu_QR_Interface_def.hpp"

#if defined(HAVE_MUELU_INST_DOUBLE_INT_INT) || defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT)
template class MueLu::QR_Interface<double, int>;
#endif

#ifdef HAVE_MUELU_INST_COMPLEX_INT_INT
# ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
template class MueLu::QR_Interface<std::complex<double>, int>;
# else
# warning To compile MueLu with 'complex' support, please turn on Teuchos_ENABLE_COMPLEX
# endif
#endif

#ifdef HAVE_MUELU_INST_PCE_INT_INT
#include "Stokhos_Sacado.hpp"
typedef Stokhos::StandardStorage<int,double> Storage;
typedef Sacado::PCE::OrthogPoly<double,Storage> pce_type;
template class MueLu::QR_Interface<pce_type, int>;
#endif

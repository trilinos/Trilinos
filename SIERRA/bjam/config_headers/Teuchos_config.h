#include "TrilinosSierraConfig.h"
#include "Teuchos_config_fcs.h"
#if ! (defined (__IBMCPP__) || defined(__IBMC__))
# define HAVE_TEUCHOS_BLASFLOAT_DOUBLE_RETURN
# define HAVE_SLAPY2_PROBLEM
# define HAVE_SLAPY2_DOUBLE_RETURN
# define HAVE_FIXABLE_COMPLEX_BLAS_PROBLEM
#endif

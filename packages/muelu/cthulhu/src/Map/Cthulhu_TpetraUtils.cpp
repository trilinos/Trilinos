#include "Cthulhu_TpetraUtils.hpp"
#include "Cthulhu_Exceptions.hpp"

#ifdef HAVE_CTHULHU_TPETRA

namespace Cthulhu {

  Tpetra::LocalGlobal toTpetra(LocalGlobal lg) {

    if (lg == Cthulhu::LocallyReplicated)
      return Tpetra::LocallyReplicated;
    if (lg == Cthulhu::GloballyDistributed)
      return Tpetra::GloballyDistributed;
    
    TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::RuntimeError, "Unknown LocalGlobal"); 

  }

}

#endif // HAVE_CTHULHU_TPETRA

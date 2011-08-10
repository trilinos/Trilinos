#include "Cthulhu_LookupStatus.hpp"
#include "Cthulhu_Exceptions.hpp"

#ifdef HAVE_CTHULHU_TPETRA

namespace Cthulhu {

  const Tpetra::LookupStatus toTpetra(Cthulhu::LookupStatus LS) {
    
    if (LS == Cthulhu::AllIDsPresent)
      return Tpetra::AllIDsPresent;
    if (LS == Cthulhu::IDNotPresent)
      return Tpetra::IDNotPresent;
    
    TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::RuntimeError, "Unknown LookupStatus"); 
    
  }

  const Cthulhu::LookupStatus toCthulhu(Tpetra::LookupStatus LS) {
    
    if (LS == Tpetra::AllIDsPresent)
      return Cthulhu::AllIDsPresent;
    if (LS == Tpetra::IDNotPresent)
      return Cthulhu::IDNotPresent;
    
    TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::RuntimeError, "Unknown LookupStatus"); 
    
  }
  
}

#endif // HAVE_CTHULHU_TPETRA

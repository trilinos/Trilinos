#include "Cthulhu_Utils.hpp"
#include "Cthulhu_Exceptions.hpp"

namespace Cthulhu {

#ifdef HAVE_CTHULHU_TPETRA

  Cthulhu::LookupStatus toCthulhu(Tpetra::LookupStatus ls) {
    
    if (ls == Tpetra::AllIDsPresent)
      return Cthulhu::AllIDsPresent;
    if (ls == Tpetra::IDNotPresent)
      return Cthulhu::IDNotPresent;
    
    TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::RuntimeError, "Unknown LookupStatus"); 
    
  }

  Tpetra::ProfileType toTpetra(Cthulhu::ProfileType pt) {

    if (pt == Cthulhu::StaticProfile)
      return Tpetra::StaticProfile;
    if (pt == Cthulhu::DynamicProfile)
      return Tpetra::DynamicProfile;
    
    TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::RuntimeError, "Unknown ProfileType"); 

  }

  Tpetra::OptimizeOption toTpetra(Cthulhu::OptimizeOption os) {

    if (os == Cthulhu::DoOptimizeStorage)
      return Tpetra::DoOptimizeStorage;
    if (os == Cthulhu::DoNotOptimizeStorage)
      return Tpetra::DoNotOptimizeStorage;
    
    TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::RuntimeError, "Unknown OptimizeOption"); 

  }

  Tpetra::CombineMode toTpetra(Cthulhu::CombineMode cm) { 
  
    if (cm == Cthulhu::ADD)
      return Tpetra::ADD;
  
    if (cm == Cthulhu::INSERT)
      return Tpetra::INSERT;
  
    if (cm == Cthulhu::ABSMAX) {
      return Tpetra::ABSMAX;
    }
  
    TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::RuntimeError, "Cannot convert Cthulhu::CombineMode to Tpetra::CombineMode: unsupported CombineMode."); 

  }

  Tpetra::LocalGlobal toTpetra(LocalGlobal lg) {

    if (lg == Cthulhu::LocallyReplicated)
      return Tpetra::LocallyReplicated;
    if (lg == Cthulhu::GloballyDistributed)
      return Tpetra::GloballyDistributed;
    
    TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::RuntimeError, "Unknown LocalGlobal"); 

  }

#endif // HAVE_CTHULHU_TPETRA

#ifdef HAVE_CTHULHU_EPETRA

  Cthulhu::LookupStatus toCthulhu(int) {
    return Cthulhu::AllIDsPresent; // TODO: manage error of EpetraMap RemoteIDList (return -1) + return the correct LookupStatus
  }

  bool toEpetra(Cthulhu::ProfileType pt) {
    
    if (pt == Cthulhu::StaticProfile)
      return true;
    if (pt == Cthulhu::DynamicProfile)
      return false;
    
    TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::RuntimeError, "Unknown ProfileType"); 
  }

  bool toEpetra(Cthulhu::OptimizeOption os) {

    if (os == Cthulhu::DoOptimizeStorage)
      return true;
    if (os == Cthulhu::DoNotOptimizeStorage)
      return false;
    
    TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::RuntimeError, "Unknown OptimizeOption"); 

  }
  
  Epetra_CombineMode toEpetra(Cthulhu::CombineMode cm) { 
    // Note: all the CombineMode are not supported.
    // According to Chris B., the behavior in Tpetra is the same as Epetra but I prefer to limit my tests for now.
    // See also the discussion of March 22 on the Tpetra developers mailing list.

    if (cm == Cthulhu::ADD)
      return Add;
    if (cm == Cthulhu::INSERT)
      return Insert;
    if (cm == Cthulhu::ABSMAX)
      return AbsMax;
  
    TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::RuntimeError, "Cannot convert Cthulhu::CombineMode to Epetra_CombineMode: unsupported CombineMode."); 

  }

#endif // HAVE_CTHULHU_EPETRA

} // namespace cthulhu

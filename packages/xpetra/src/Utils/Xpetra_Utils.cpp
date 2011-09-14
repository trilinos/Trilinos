#include "Xpetra_Utils.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

  std::string toString(Xpetra::UnderlyingLib lib) {
    if (lib == Xpetra::UseTpetra) {
      return "Tpetra";
    } else if (lib == Xpetra::UseEpetra) {
      return "Epetra";
    } else {
      TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "lib != UseTpetra && lib != UseEpetra");
    }
  }

#ifdef HAVE_XPETRA_TPETRA

  Xpetra::LookupStatus toXpetra(Tpetra::LookupStatus ls) {
    
    if (ls == Tpetra::AllIDsPresent)
      return Xpetra::AllIDsPresent;
    if (ls == Tpetra::IDNotPresent)
      return Xpetra::IDNotPresent;
    
    TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::RuntimeError, "Unknown LookupStatus"); 
    
  }

  Tpetra::ProfileType toTpetra(Xpetra::ProfileType pt) {

    if (pt == Xpetra::StaticProfile)
      return Tpetra::StaticProfile;
    if (pt == Xpetra::DynamicProfile)
      return Tpetra::DynamicProfile;
    
    TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::RuntimeError, "Unknown ProfileType"); 

  }

  Tpetra::OptimizeOption toTpetra(Xpetra::OptimizeOption os) {

    if (os == Xpetra::DoOptimizeStorage)
      return Tpetra::DoOptimizeStorage;
    if (os == Xpetra::DoNotOptimizeStorage)
      return Tpetra::DoNotOptimizeStorage;
    
    TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::RuntimeError, "Unknown OptimizeOption"); 

  }

  Tpetra::CombineMode toTpetra(Xpetra::CombineMode cm) { 
  
    if (cm == Xpetra::ADD)
      return Tpetra::ADD;
  
    if (cm == Xpetra::INSERT)
      return Tpetra::INSERT;
  
    if (cm == Xpetra::ABSMAX) {
      return Tpetra::ABSMAX;
    }
  
    TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::RuntimeError, "Cannot convert Xpetra::CombineMode to Tpetra::CombineMode: unsupported CombineMode."); 

  }

  Tpetra::LocalGlobal toTpetra(LocalGlobal lg) {

    if (lg == Xpetra::LocallyReplicated)
      return Tpetra::LocallyReplicated;
    if (lg == Xpetra::GloballyDistributed)
      return Tpetra::GloballyDistributed;
    
    TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::RuntimeError, "Unknown LocalGlobal"); 

  }

#endif // HAVE_XPETRA_TPETRA

#ifdef HAVE_XPETRA_EPETRA

  Xpetra::LookupStatus toXpetra(int) {
    return Xpetra::AllIDsPresent; // TODO: manage error of EpetraMap RemoteIDList (return -1) + return the correct LookupStatus
  }

  bool toEpetra(Xpetra::ProfileType pt) {
    
    if (pt == Xpetra::StaticProfile)
      return true;
    if (pt == Xpetra::DynamicProfile)
      return false;
    
    TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::RuntimeError, "Unknown ProfileType"); 
  }

  bool toEpetra(Xpetra::OptimizeOption os) {

    if (os == Xpetra::DoOptimizeStorage)
      return true;
    if (os == Xpetra::DoNotOptimizeStorage)
      return false;
    
    TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::RuntimeError, "Unknown OptimizeOption"); 

  }
  
  Epetra_CombineMode toEpetra(Xpetra::CombineMode cm) { 
    // Note: all the CombineMode are not supported.
    // According to Chris B., the behavior in Tpetra is the same as Epetra but I prefer to limit my tests for now.
    // See also the discussion of March 22 on the Tpetra developers mailing list.

    if (cm == Xpetra::ADD)
      return Add;
    if (cm == Xpetra::INSERT)
      return Insert;
    if (cm == Xpetra::ABSMAX)
      return AbsMax;
  
    TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::RuntimeError, "Cannot convert Xpetra::CombineMode to Epetra_CombineMode: unsupported CombineMode."); 

  }

#endif // HAVE_XPETRA_EPETRA

} // namespace xpetra

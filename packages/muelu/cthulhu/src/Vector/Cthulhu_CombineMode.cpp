#include "Cthulhu_CombineMode.hpp"
#include "Cthulhu_Exceptions.hpp"

//JG: Cthulhu currently supports only INSERT and ABSMAX for Vector/MultiVector.
//CombineMode semantics for a DistObject depends on the underlying class. According to Chris B., the behavior in Tpetra is the same as Epetra but I prefer to limit my tests for now.
//See also the discussion of March 22 on the Tpetra developers mailing list.

#ifdef HAVE_CTHULHU_EPETRA

namespace Cthulhu {

  //! Convert a Cthulhu Combine Mode to an Epetra Combine Mode.
  Epetra_CombineMode toEpetra(Cthulhu::CombineMode CM) { 
  
    if (CM == Cthulhu::ADD)
      return Add;
    if (CM == Cthulhu::INSERT)
      return Insert;
    if (CM == Cthulhu::ABSMAX)
      return AbsMax;
  
    TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::RuntimeError, "Cannot convert Cthulhu::CombineMode to Epetra_CombineMode: unsupported CombineMode."); 

  }

}

#endif

#ifdef HAVE_CTHULHU_TPETRA

namespace Cthulhu {

  //! Convert a Cthulhu Combine Mode to a Tpetra Combine Mode.
  Tpetra::CombineMode toTpetra(Cthulhu::CombineMode CM) { 
  
    if (CM == Cthulhu::ADD)
      return Tpetra::ADD;
  
    if (CM == Cthulhu::INSERT)
      return Tpetra::INSERT;
  
    if (CM == Cthulhu::ABSMAX) {
      //std::cerr << "toTpetra ERROR !!!" << std::endl;
      //return Tpetra::INSERT;
      return Tpetra::ABSMAX; //TODO
    }
  
    TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::RuntimeError, "Cannot convert Cthulhu::CombineMode to Tpetra::CombineMode: unsupported CombineMode."); 

  }

}
#endif

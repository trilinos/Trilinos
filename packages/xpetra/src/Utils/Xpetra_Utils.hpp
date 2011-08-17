#ifndef CTHULHU_LOOKUPSTATUS_HPP
#define CTHULHU_LOOKUPSTATUS_HPP

#include "Cthulhu_ConfigDefs.hpp"
#ifdef HAVE_CTHULHU_TPETRA  
#include "Tpetra_ConfigDefs.hpp"
#endif

#ifdef HAVE_CTHULHU_EPETRA  
#include "Epetra_CombineMode.h"
#endif

namespace Cthulhu {

#ifdef HAVE_CTHULHU_TPETRA  
  
  //! Convert a Tpetra::LookupStatus to a Cthulhu::LookupStatus.
  Cthulhu::LookupStatus  toCthulhu(Tpetra::LookupStatus);

  //! Convert a Cthulhu::OptimizeOption to a Tpetra::OptimizeOption.
  Tpetra::ProfileType    toTpetra(Cthulhu::ProfileType);

  //! Convert a Cthulhu::OptimizeOption to a Tpetra::OptimizeOption.
  Tpetra::OptimizeOption toTpetra(Cthulhu::OptimizeOption);

  //! Convert a Cthulhu::CombineMode to a Tpetra::CombineMode.
  Tpetra::CombineMode    toTpetra(Cthulhu::CombineMode CM);

  //! Convert a Cthulhu::LocalGlobal to a Tpetra::LocalGlobal.  
  Tpetra::LocalGlobal    toTpetra(LocalGlobal lg);

#endif // HAVE_CTHULHU_TPETRA

#ifdef HAVE_CTHULHU_EPETRA  
  
  //! Convert a Epetra return value to a Cthulhu::LookupStatus.
  Cthulhu::LookupStatus  toCthulhu(int);
  
  //! Convert a Cthulhu::ProfileType to an Epetra StaticProfil boolean
  bool                   toEpetra(Cthulhu::ProfileType);
  
  //! Convert a Cthulhu::OptimizeOption to an Epetra OptimizeDataStorage boolean
  bool                   toEpetra(Cthulhu::OptimizeOption);
  
  //! Convert a Cthulhu::CombineMode to an Epetra_CombineMode.
  Epetra_CombineMode     toEpetra(Cthulhu::CombineMode CM);

#endif // HAVE_CTHULHU_EPETRA

}

#endif // CTHULHU_LOOKUPSTATUS_HPP

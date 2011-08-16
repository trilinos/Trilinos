#ifndef CTHULHU_LOOKUPSTATUS_HPP
#define CTHULHU_LOOKUPSTATUS_HPP

#include "Cthulhu_ConfigDefs.hpp"
#ifdef HAVE_CTHULHU_TPETRA  
#include "Tpetra_ConfigDefs.hpp"
#endif

namespace Cthulhu {

#ifdef HAVE_CTHULHU_TPETRA  
  const Cthulhu::LookupStatus  toCthulhu(Tpetra::LookupStatus);
  const Tpetra::ProfileType    toTpetra(Cthulhu::ProfileType);
  const Tpetra::OptimizeOption toTpetra(Cthulhu::OptimizeOption);
#endif // HAVE_CTHULHU_TPETRA

#ifdef HAVE_CTHULHU_EPETRA  
  const Cthulhu::LookupStatus  toCthulhu(int);
  const bool                   toEpetra(Cthulhu::ProfileType);
  const bool                   toEpetra(Cthulhu::OptimizeOption);
#endif // HAVE_CTHULHU_TPETRA

}

#endif // CTHULHU_LOOKUPSTATUS_HPP

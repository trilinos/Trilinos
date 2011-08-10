#ifndef CTHULHU_LOOKUPSTATUS_HPP
#define CTHULHU_LOOKUPSTATUS_HPP

#include "Cthulhu_ConfigDefs.hpp"
#ifdef HAVE_CTHULHU_TPETRA  
#include "Tpetra_ConfigDefs.hpp"
#endif

namespace Cthulhu {

#ifdef HAVE_CTHULHU_TPETRA  
  const Tpetra::LookupStatus toTpetra(Cthulhu::LookupStatus);
  const Cthulhu::LookupStatus toCthulhu(Tpetra::LookupStatus);
  
#endif // HAVE_CTHULHU_TPETRA

}

#endif // CTHULHU_LOOKUPSTATUS_HPP

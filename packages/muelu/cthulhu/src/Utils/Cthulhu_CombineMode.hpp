#ifndef CTHULHU_COMBINEMODE_HPP
#define CTHULHU_COMBINEMODE_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifdef HAVE_CTHULHU_EPETRA
#include <Epetra_CombineMode.h>

namespace Cthulhu {

  //! Convert a Cthulhu Combine Mode to an Epetra Combine Mode.
  Epetra_CombineMode toEpetra(Cthulhu::CombineMode CM);

}

#endif

#ifdef HAVE_CTHULHU_TPETRA

#include <Tpetra_ConfigDefs.hpp>

namespace Cthulhu {

  //! Convert a Cthulhu Combine Mode to a Tpetra Combine Mode.
  Tpetra::CombineMode toTpetra(Cthulhu::CombineMode CM);

}

#endif

#endif

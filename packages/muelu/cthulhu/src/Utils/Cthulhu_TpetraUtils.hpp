#ifndef CTHULHU_TPETRAUTILS_HPP
#define CTHULHU_TPETRAUTILS_HPP

#include "Cthulhu_TpetraConfigDefs.hpp"
#include "Tpetra_ConfigDefs.hpp"

namespace Cthulhu {
  Tpetra::LocalGlobal toTpetra(LocalGlobal lg);
}

#endif // CTHULHU_TPETRAUTILS_HPP

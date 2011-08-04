#ifndef CTHULHU_COMBINEMODE_HPP
#define CTHULHU_COMBINEMODE_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifdef HAVE_CTHULHU_EPETRA

#include <Epetra_CombineMode.h>

//! Convert a Cthulhu Combine Mode to an Epetra Combine Mode.
const Epetra_CombineMode Cthulhu2Epetra_CombineMode(const Cthulhu::CombineMode& CM);

#endif

#ifdef HAVE_CTHULHU_TPETRA

#include <Tpetra_ConfigDefs.hpp>

//! Convert a Cthulhu Combine Mode to a Tpetra Combine Mode.
const Tpetra::CombineMode Cthulhu2Tpetra_CombineMode(const Cthulhu::CombineMode& CM);

#endif

#endif

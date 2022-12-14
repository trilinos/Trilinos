#ifndef MUELU_ETI_4ARGUMENT_HPP
#define MUELU_ETI_4ARGUMENT_HPP

// The macro "MUELU_ETI_GROUP" must be defined prior to including this file.

// We need to define these typedefs as it is not possible to properly expand
// macros with colons in them
# include <TpetraCore_config.h>
# include <TpetraCore_ETIHelperMacros.h>
TPETRA_ETI_MANGLING_TYPEDEFS()

// Epetra = on, Tpetra = off

// Epetra = on, Tpetra = on

// Epetra = off, Tpetra = on
  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(MUELU_ETI_GROUP)

#endif //ifndef MUELU_ETI_4ARGUMENT_HPP

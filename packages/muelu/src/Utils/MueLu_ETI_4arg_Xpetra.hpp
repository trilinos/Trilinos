#ifndef MUELU_ETI_5ARGUMENT_HPP
#define MUELU_ETI_5ARGUMENT_HPP

// The macros "MUELU_ETI_GROUP_EPETRA" and "MUELU_ETI_GROUP_TPETRA" must be defined prior to including this file.

// We need to define these typedefs as it is not possible to properly expand
// macros with colons in them
#if defined(HAVE_MUELU_TPETRA)
# include <TpetraCore_config.h>
# include <TpetraCore_ETIHelperMacros.h>
TPETRA_ETI_MANGLING_TYPEDEFS()
#endif
#if defined(HAVE_MUELU_EPETRA)
# include <Epetra_config.h>
#endif

#if   (defined(HAVE_MUELU_EPETRA) &&  defined(EPETRA_HAVE_OMP) && (!defined(HAVE_MUELU_TPETRA) || !defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT)))
  // Epetra is enabled with OpenMP node, but Tpetra is a) not enabled, or b) is not instantiated on OpenMP, or c) is not instantiated on OpenMP with <double,int,int>
  typedef Kokkos::Compat::KokkosOpenMPWrapperNode EpetraNode;
#elif (defined(HAVE_MUELU_EPETRA) && !defined(EPETRA_HAVE_OMP) && (!defined(HAVE_MUELU_TPETRA) || !defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT)))
  // Epetra is enabled with Serial node, but Tpetra is a) not enabled, or b) is not instantiated on Serial, or c) is not instantiated on Serial with <double,int,int>
  typedef Kokkos::Compat::KokkosSerialWrapperNode EpetraNode;
#endif

// Epetra = on, Tpetra = off
#if defined(HAVE_MUELU_EPETRA) && !defined(HAVE_MUELU_TPETRA)
  MUELU_ETI_GROUP_EPETRA(double,int,int,EpetraNode)
#endif



// Epetra = on, Tpetra = on
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_TPETRA)
  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(MUELU_ETI_GROUP_TPETRA)
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
  MUELU_ETI_GROUP_EPETRA(double,int,int,EpetraNode)
# endif

#endif

// Epetra = off, Tpetra = on
#if !defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_TPETRA)
  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(MUELU_ETI_GROUP_TPETRA)
#endif

#endif //ifndef MUELU_ETI_4ARGUMENT_HPP

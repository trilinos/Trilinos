#ifndef XPETRA_ETI_4ARGUMENT_SLGN_HPP
#define XPETRA_ETI_4ARGUMENT_SLGN_HPP

// The macro "XPETRA_ETI_GROUP" must be defined prior to including this file.

// We need to define these typedefs as it is not possible to properly expand
// macros with colons in them
#if defined(HAVE_XPETRA_TPETRA)
# include <TpetraCore_config.h>
# include <TpetraCore_ETIHelperMacros.h>
TPETRA_ETI_MANGLING_TYPEDEFS()
#endif
#if defined(HAVE_XPETRA_EPETRA)
# include <Epetra_config.h>
#endif

// This is a kludge relying on the fact that LO=int always and everywhere
#define INT_LO_GO_NO_GROUP(LO,GO,NO) XPETRA_ETI_GROUP(int,LO,GO,NO)
#define INT_LONG_NO_GROUP(NO) XPETRA_ETI_GROUP(int,int,long,NO)
#define INT_LONG_LONG_NO_GROUP(NO) XPETRA_ETI_GROUP(int,int, long long,NO)

#if   (defined(HAVE_XPETRA_EPETRA) &&  defined(EPETRA_HAVE_OMP) && (!defined(HAVE_XPETRA_TPETRA) || !defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT)))
  // Epetra is enabled with OpenMP node, but Tpetra is a) not enabled, or b) is not instantiated on OpenMP, or c) is not instantiated on OpenMP with <double,int,int>
  typedef Kokkos::Compat::KokkosOpenMPWrapperNode EpetraNode;
#elif (defined(HAVE_XPETRA_EPETRA) && !defined(EPETRA_HAVE_OMP) && (!defined(HAVE_XPETRA_TPETRA) || !defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT)))
  // Epetra is enabled with Serial node, but Tpetra is a) not enabled, or b) is not instantiated on Serial, or c) is not instantiated on Serial with <double,int,int>
  typedef Kokkos::Compat::KokkosSerialWrapperNode EpetraNode;
#endif

// Epetra = on, Tpetra = off
#if defined(HAVE_XPETRA_EPETRA) && !defined(HAVE_XPETRA_TPETRA)
  XPETRA_ETI_GROUP(double,int,int,EpetraNode)
  XPETRA_ETI_GROUP(int,int,int,EpetraNode)
#endif

// Epetra = on, Tpetra = on
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_TPETRA)
  TPETRA_INSTANTIATE_SLGN(XPETRA_ETI_GROUP)

#if defined(HAVE_TPETRA_INST_INT_INT)
  #if defined(HAVE_TPETRA_INST_INT_LONG)
  TPETRA_INSTANTIATE_N(INT_LONG_NO_GROUP)
  #endif
  #if defined(HAVE_TPETRA_INST_INT_LONG_LONG)
  TPETRA_INSTANTIATE_N(INT_LONG_LONG_NO_GROUP)
  #endif
#else
  TPETRA_INSTANTIATE_LGN(INT_LO_GO_NO_GROUP)
#endif


#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
  XPETRA_ETI_GROUP(double,int,int,EpetraNode)
  XPETRA_ETI_GROUP(int,int,int,EpetraNode)
#endif

#endif

// Epetra = off, Tpetra = on
#if !defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_TPETRA)
  TPETRA_INSTANTIATE_SLGN(XPETRA_ETI_GROUP)

#if defined(HAVE_TPETRA_INST_INT_INT)
  #if defined(HAVE_TPETRA_INST_INT_LONG)
  TPETRA_INSTANTIATE_N(INT_LONG_NO_GROUP)
  #endif
  #if defined(HAVE_TPETRA_INST_INT_LONG_LONG)
  TPETRA_INSTANTIATE_N(INT_LONG_LONG_NO_GROUP)
  #endif
#else
  TPETRA_INSTANTIATE_LGN(INT_LO_GO_NO_GROUP)
#endif
#endif

#endif //ifndef XPETRA_ETI_4ARGUMENT_SLGN_HPP

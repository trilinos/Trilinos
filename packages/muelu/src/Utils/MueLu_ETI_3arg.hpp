#ifndef MUELU_ETI_3ARGUMENT_HPP
#define MUELU_ETI_3ARGUMENT_HPP

#if defined(HAVE_MUELU_EPETRA)  && !defined(HAVE_MUELU_TPETRA)
  MUELU_ETI_GROUP(double,int,int,Kokkos::Serial);
#endif

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_TPETRA)
# include <TpetraCore_config.h>
# include <TpetraCore_ETIHelperMacros.h>

  TPETRA_ETI_MANGLING_TYPEDEFS()
  TPETRA_INSTANTIATE_LGN(MUELU_ETI_GROUP) 
#if !defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT)
    MUELU_ETI_GROUP(double,int,int,Kokkos::Serial);
# endif

#endif

#if !defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_TPETRA)
# include <TpetraCore_config.h>
# include <TpetraCore_ETIHelperMacros.h>

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN(MUELU_ETI_GROUP)
#endif

#endif //ifndef MUELU_ETI_3ARGUMENT_HPP

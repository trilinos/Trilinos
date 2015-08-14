#ifndef IFPACK2_ETIHELPERMACROS_H_
#define IFPACK2_ETIHELPERMACROS_H_

#include <Ifpack2_ConfigDefs.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>

#define IFPACK2_INSTANTIATE_SL(INSTMACRO)\
	INSTMACRO( double , int )


#define IFPACK2_INSTANTIATE_SL_REAL(INSTMACRO)\
	INSTMACRO( double , int , int )\
	INSTMACRO( double , int , long )


#define IFPACK2_INSTANTIATE_L(INSTMACRO)\
	INSTMACRO( int )


#define IFPACK2_INSTANTIATE_SLG(INSTMACRO)\
	INSTMACRO( double , int , int )\
	INSTMACRO( double , int , long )


#define IFPACK2_INSTANTIATE_SLG_REAL(INSTMACRO)\
	INSTMACRO( double , int , int )\
	INSTMACRO( double , int , long )


#define IFPACK2_INSTANTIATE_LG(INSTMACRO)\
	INSTMACRO( int , int )\
	INSTMACRO( int , long )


#define IFPACK2_INSTANTIATE_SLGN(INSTMACRO)\
	INSTMACRO( double , int , int , Kokkos_Compat_KokkosSerialWrapperNode )\
	INSTMACRO( double , int , long , Kokkos_Compat_KokkosSerialWrapperNode )


#define IFPACK2_INSTANTIATE_SLGN_REAL(INSTMACRO)\
	INSTMACRO( double , int , int , Kokkos_Compat_KokkosSerialWrapperNode )\
	INSTMACRO( double , int , long , Kokkos_Compat_KokkosSerialWrapperNode )


#define IFPACK2_INSTANTIATE_LGN(INSTMACRO)\
	INSTMACRO( double , int , int , Kokkos_Compat_KokkosSerialWrapperNode )\
	INSTMACRO( double , int , long , Kokkos_Compat_KokkosSerialWrapperNode )


#define IFPACK2_ETI_MANGLING_TYPEDEFS()  \
	typedef Kokkos::Compat::KokkosSerialWrapperNode Kokkos_Compat_KokkosSerialWrapperNode;

#endif // IFPACK2_ETIHELPERMACROS_H_

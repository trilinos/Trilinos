#ifndef IFPACK2_ETIHELPERMACROS_H_
#define IFPACK2_ETIHELPERMACROS_H_

#include <Ifpack2_ConfigDefs.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>

#define IFPACK2_INSTANTIATE_SLG(INSTMACRO )\
   INSTMACRO( double , int , int )\
   INSTMACRO( double , int , long )\
   INSTMACRO( std_complex0double0 , int , int )\
   INSTMACRO( std_complex0double0 , int , long )


#define IFPACK2_INSTANTIATE_SLG_call(INSTMACRO, CALL_MACRO )\
   INSTMACRO( double , int , int, CALL_MACRO )\
   INSTMACRO( double , int , long, CALL_MACRO )\
   INSTMACRO( std_complex0double0 , int , int, CALL_MACRO )\
   INSTMACRO( std_complex0double0 , int , long, CALL_MACRO )


#define IFPACK2_INSTANTIATE_SLG_REAL(INSTMACRO )\
   INSTMACRO( double , int , int )\
   INSTMACRO( double , int , long )\
   INSTMACRO( std_complex0double0 , int , int )\
   INSTMACRO( std_complex0double0 , int , long )

#define IFPACK2_INSTANTIATE_SLG_REAL_call(INSTMACRO, CALL_MACRO )\
   INSTMACRO( double , int , int, CALL_MACRO )\
   INSTMACRO( double , int , long, CALL_MACRO )\
   INSTMACRO( std_complex0double0 , int , int, CALL_MACRO )\
   INSTMACRO( std_complex0double0 , int , long, CALL_MACRO )


#define IFPACK2_INSTANTIATE_LG(INSTMACRO)\
	INSTMACRO( int , int )\
	INSTMACRO( int , long )


#define IFPACK2_INSTANTIATE_LG_call(INSTMACRO, CALL_MACRO) \
   INSTMACRO( int , int, CALL_MACRO )                                 \
   INSTMACRO( int , long, CALL_MACRO )


#ifdef HAVE_TPETRA_INST_SERIAL
   #define Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNode( SCALAR, LO, GO, INSTAMACRO ) \
      INSTAMACRO( SCALAR, LO, GO, Kokkos_Compat_KokkosSerialWrapperNode )
   #define Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNodeNS( LO, GO, INSTAMACRO ) \
      INSTAMACRO( LO, GO, Kokkos_Compat_KokkosSerialWrapperNode )
   #define Inst_Kokkos_Compat_KokkosSerialWrapperNode( INSTAMACRO ) \
      INSTAMACRO( Kokkos_Compat_KokkosSerialWrapperNode )
   #define IFPACK2_ETI_SERIAL_TYPEDEF\
      typedef Kokkos::Compat::KokkosSerialWrapperNode Kokkos_Compat_KokkosSerialWrapperNode;
#else
   #define Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNode( SCALAR, LO, GO, INSTAMACRO )
   #define Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNodeNS( LO, GO, INSTAMACRO )
   #define Inst_Kokkos_Compat_KokkosSerialWrapperNode( INSTAMACRO )
   #define IFPACK2_ETI_SERIAL_TYPEDEF
#endif

#ifdef HAVE_TPETRA_INST_PTHREAD
   #define Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNode( SCALAR, LO, GO, INSTAMACRO ) \
      INSTAMACRO( SCALAR, LO, GO, Kokkos_Compat_KokkosThreadsWrapperNode )
   #define Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNodeNS( LO, GO, INSTAMACRO ) \
      INSTAMACRO( LO, GO, Kokkos_Compat_KokkosThreadsWrapperNode )
   #define Inst_Kokkos_Compat_KokkosThreadsWrapperNode( INSTAMACRO ) \
      INSTAMACRO( Kokkos_Compat_KokkosThreadsWrapperNode )
   #define IFPACK2_ETI_THREADS_TYPEDEF\
      typedef Kokkos::Compat::KokkosThreadsWrapperNode Kokkos_Compat_KokkosThreadsWrapperNode;
#else
   #define Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNode( SCALAR, LO, GO, INSTAMACRO )
   #define Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNodeNS( LO, GO, INSTAMACRO )
   #define Inst_Kokkos_Compat_KokkosThreadsWrapperNode( INSTAMACRO )
   #define IFPACK2_ETI_THREADS_TYPEDEF
#endif

#ifdef HAVE_TPETRA_INST_OPENMP
   #define Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNode( SCALAR, LO, GO, INSTAMACRO ) \
      INSTAMACRO( SCALAR, LO, GO, Kokkos_Compat_KokkosOpenMPWrapperNode )
   #define Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNodeNS( LO, GO, INSTAMACRO ) \
      INSTAMACRO( LO, GO, Kokkos_Compat_KokkosOpenMPWrapperNode )
   #define Inst_Kokkos_Compat_KokkosOpenMPWrapperNode( INSTAMACRO ) \
      INSTAMACRO( Kokkos_Compat_KokkosOpenMPWrapperNode )
   #define IFPACK2_ETI_OPENMP_TYPEDEF\
	typedef Kokkos::Compat::KokkosOpenMPWrapperNode Kokkos_Compat_KokkosOpenMPWrapperNode;
#else
   #define Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNode( SCALAR, LO, GO, INSTAMACRO )
   #define Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNodeNS( LO, GO, INSTAMACRO )
   #define Inst_Kokkos_Compat_KokkosOpenMPWrapperNode( INSTAMACRO )
   #define IFPACK2_ETI_OPENMP_TYPEDEF
#endif

#ifdef HAVE_TPETRA_INST_CUDA
   #define Exp_Macro_Kokkos_Compat_KokkosCudaWrapperNode( SCALAR, LO, GO, INSTAMACRO ) \
      INSTAMACRO( SCALAR, LO, GO, Kokkos_Compat_KokkosCudaWrapperNode )
   #define Exp_Macro_Kokkos_Compat_KokkosCudaWrapperNodeNS( LO, GO, INSTAMACRO ) \
      INSTAMACRO( LO, GO, Kokkos_Compat_KokkosCudaWrapperNode )
   #define Inst_Kokkos_Compat_KokkosCudaWrapperNode( INSTAMACRO ) \
      INSTAMACRO( Kokkos_Compat_KokkosCudaWrapperNode )
   #define IFPACK2_ETI_CUDA_TYPEDEF
      typedef Kokkos::Compat::KokkosOpenMPWrapperNode Kokkos_Compat_KokkosCudaWrapperNode;
#else
   #define Exp_Macro_Kokkos_Compat_KokkosCudaWrapperNode( SCALAR, LO, GO, INSTAMACRO )
   #define Exp_Macro_Kokkos_Compat_KokkosCudaWrapperNodeNS( LO, GO, INSTAMACRO )
   #define Inst_Kokkos_Compat_KokkosCudaWrapperNode( INSTAMACRO )
   #define IFPACK2_ETI_CUDA_TYPEDEF
#endif


#define IFPACK2_INSTANTIATE_SLGN(INSTMACRO)\
   IFPACK2_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNode, INSTMACRO )\
   IFPACK2_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNode, INSTMACRO )\
   IFPACK2_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNode, INSTMACRO )\
   IFPACK2_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosCudaWrapperNode, INSTMACRO )


#define IFPACK2_INSTANTIATE_SLGN_REAL(INSTMACRO)\
   IFPACK2_INSTANTIATE_SLG_REAL_call( Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNode, INSTMACRO )\
   IFPACK2_INSTANTIATE_SLG_REAL_call( Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNode, INSTMACRO )\
   IFPACK2_INSTANTIATE_SLG_REAL_call( Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNode, INSTMACRO )\
   IFPACK2_INSTANTIATE_SLG_REAL_call( Exp_Macro_Kokkos_Compat_KokkosCudaWrapperNode, INSTMACRO )


#define IFPACK2_INSTANTIATE_LGN(INSTMACRO)\
   IFPACK2_INSTANTIATE_LG_call( Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNodeNS, INSTMACRO )\
   IFPACK2_INSTANTIATE_LG_call( Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNodeNS, INSTMACRO )\
   IFPACK2_INSTANTIATE_LG_call( Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNodeNS, INSTMACRO )\
   IFPACK2_INSTANTIATE_LG_call( Exp_Macro_Kokkos_Compat_KokkosCudaWrapperNodeNS, INSTMACRO )


#define IFPACK2_ETI_MANGLING_TYPEDEFS()\
	IFPACK2_ETI_SERIAL_TYPEDEF \
	IFPACK2_ETI_THREADS_TYPEDEF \
	IFPACK2_ETI_OPENMP_TYPEDEF \
	IFPACK2_ETI_CUDA_TYPEDEF \
	typedef std::complex<double> std_complex0double0; 

#endif // IFPACK2_ETIHELPERMACROS_H_

#ifndef TPETRACORE_ETIHELPERMACROS_H
#define TPETRACORE_ETIHELPERMACROS_H

#include <Tpetra_ConfigDefs.hpp>

/* Tpetra provides official support for the following nodes */
#include <Kokkos_DefaultNode.hpp>

/* Tpetra provides official support for dd_real and qd_real */
#if defined(HAVE_TPETRA_QD)
#include <qd/qd_real.h>
#endif

#define TPETRA_INSTANTIATE_SLG(INSTMACRO )\
   INSTMACRO( int , int , int )\
   INSTMACRO( int , int , long )\
   INSTMACRO( int , long , long )\
   INSTMACRO( long , int , long )\
   INSTMACRO( long , long , long )\
   INSTMACRO( double , int , int )\
   INSTMACRO( double , int , long )\
   INSTMACRO( double , long , long )\
   INSTMACRO( std_complex0double0 , int , int )\
   INSTMACRO( std_complex0double0 , int , long )\
   INSTMACRO( std_complex0double0 , long , long )

#define TPETRA_INSTANTIATE_SLG_call(INSTMACRO, CALL_MACRO )\
   INSTMACRO( int , int , int, CALL_MACRO )\
   INSTMACRO( int , int , long, CALL_MACRO )\
   INSTMACRO( int , long , long, CALL_MACRO )\
   INSTMACRO( long , int , long, CALL_MACRO )\
   INSTMACRO( long , long , long, CALL_MACRO )\
   INSTMACRO( double , int , int, CALL_MACRO )\
   INSTMACRO( double , int , long, CALL_MACRO )\
   INSTMACRO( double , long , long, CALL_MACRO )\
   INSTMACRO( std_complex0double0 , int , int, CALL_MACRO )\
   INSTMACRO( std_complex0double0 , int , long, CALL_MACRO )\
   INSTMACRO( std_complex0double0 , long , long, CALL_MACRO )

#ifdef HAVE_TPETRA_INST_SERIAL
   #define Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNode( SCALAR, LO, GO, INSTAMACRO ) \
      INSTAMACRO( SCALAR, LO, GO, Kokkos_Compat_KokkosSerialWrapperNode )
   #define Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNodeNS( LO, GO, INSTAMACRO ) \
      INSTAMACRO( LO, GO, Kokkos_Compat_KokkosSerialWrapperNode )
   #define Inst_Kokkos_Compat_KokkosSerialWrapperNode( INSTAMACRO ) \
      INSTAMACRO( Kokkos_Compat_KokkosSerialWrapperNode )
   #define TPETRA_ETI_SERIAL_TYPEDEF\
      typedef Kokkos::Compat::KokkosSerialWrapperNode Kokkos_Compat_KokkosSerialWrapperNode;
#else
   #define Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNode( SCALAR, LO, GO, INSTAMACRO )
   #define Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNodeNS( LO, GO, INSTAMACRO )
   #define Inst_Kokkos_Compat_KokkosSerialWrapperNode( INSTAMACRO )
   #define TPETRA_ETI_SERIAL_TYPEDEF
#endif

#ifdef HAVE_TPETRA_INST_PTHREAD
   #define Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNode( SCALAR, LO, GO, INSTAMACRO ) \
      INSTAMACRO( SCALAR, LO, GO, Kokkos_Compat_KokkosThreadsWrapperNode )
   #define Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNodeNS( LO, GO, INSTAMACRO ) \
      INSTAMACRO( LO, GO, Kokkos_Compat_KokkosThreadsWrapperNode )
   #define Inst_Kokkos_Compat_KokkosThreadsWrapperNode( INSTAMACRO ) \
      INSTAMACRO( Kokkos_Compat_KokkosThreadsWrapperNode )
   #define TPETRA_ETI_THREADS_TYPEDEF\
      typedef Kokkos::Compat::KokkosThreadsWrapperNode Kokkos_Compat_KokkosThreadsWrapperNode;
#else
   #define Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNode( SCALAR, LO, GO, INSTAMACRO )
   #define Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNodeNS( LO, GO, INSTAMACRO )
   #define Inst_Kokkos_Compat_KokkosThreadsWrapperNode( INSTAMACRO )
   #define TPETRA_ETI_THREADS_TYPEDEF
#endif

#ifdef HAVE_TPETRA_INST_OPENMP
   #define Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNode( SCALAR, LO, GO, INSTAMACRO ) \
      INSTAMACRO( SCALAR, LO, GO, Kokkos_Compat_KokkosOpenMPWrapperNode )
   #define Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNodeNS( LO, GO, INSTAMACRO ) \
      INSTAMACRO( LO, GO, Kokkos_Compat_KokkosOpenMPWrapperNode )
   #define Inst_Kokkos_Compat_KokkosOpenMPWrapperNode( INSTAMACRO ) \
      INSTAMACRO( Kokkos_Compat_KokkosOpenMPWrapperNode )
   #define TPETRA_ETI_OPENMP_TYPEDEF\
	typedef Kokkos::Compat::KokkosOpenMPWrapperNode Kokkos_Compat_KokkosOpenMPWrapperNode;
#else
   #define Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNode( SCALAR, LO, GO, INSTAMACRO )
   #define Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNodeNS( LO, GO, INSTAMACRO )
   #define Inst_Kokkos_Compat_KokkosOpenMPWrapperNode( INSTAMACRO )
   #define TPETRA_ETI_OPENMP_TYPEDEF
#endif

#ifdef HAVE_TPETRA_INST_CUDA
   #define Exp_Macro_Kokkos_Compat_KokkosCudaWrapperNode( SCALAR, LO, GO, INSTAMACRO ) \
      INSTAMACRO( SCALAR, LO, GO, Kokkos_Compat_KokkosCudaWrapperNode )
   #define Exp_Macro_Kokkos_Compat_KokkosCudaWrapperNodeNS( LO, GO, INSTAMACRO ) \
      INSTAMACRO( LO, GO, Kokkos_Compat_KokkosCudaWrapperNode )
   #define Inst_Kokkos_Compat_KokkosCudaWrapperNode( INSTAMACRO ) \
      INSTAMACRO( Kokkos_Compat_KokkosCudaWrapperNode )
   #define TPETRA_ETI_CUDA_TYPEDEF
      typedef Kokkos::Compat::KokkosOpenMPWrapperNode Kokkos_Compat_KokkosCudaWrapperNode;
#else
   #define Exp_Macro_Kokkos_Compat_KokkosCudaWrapperNode( SCALAR, LO, GO, INSTAMACRO )
   #define Exp_Macro_Kokkos_Compat_KokkosCudaWrapperNodeNS( LO, GO, INSTAMACRO )
   #define Inst_Kokkos_Compat_KokkosCudaWrapperNode( INSTAMACRO )
   #define TPETRA_ETI_CUDA_TYPEDEF
#endif

#define TPETRA_INSTANTIATE_LG(INSTMACRO)\
   INSTMACRO( int , int )\
   INSTMACRO( int , long )\
   INSTMACRO( long , long )

#define TPETRA_INSTANTIATE_LG_call(INSTMACRO, CALL_MACRO)\
   INSTMACRO( int , int, CALL_MACRO )\
   INSTMACRO( int , long, CALL_MACRO )\
   INSTMACRO( long , long, CALL_MACRO )

#define TPETRA_INSTANTIATE_TSLGN(INSTMACRO)


#define TPETRA_INSTANTIATE_CONVERT(INSTMACRO)


#define TPETRA_INSTANTIATE_VECTOR(INSTMACRO)\
   TPETRA_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNode, INSTMACRO ) \
   TPETRA_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNode, INSTMACRO )\
   TPETRA_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNode, INSTMACRO )\
   TPETRA_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosCudaWrapperNode, INSTMACRO )


#define TPETRA_INSTANTIATE_SLGN(INSTMACRO)\
   TPETRA_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNode, INSTMACRO )\
   TPETRA_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNode, INSTMACRO )\
   TPETRA_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNode, INSTMACRO )\
   TPETRA_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosCudaWrapperNode, INSTMACRO )


#define TPETRA_INSTANTIATE_PLGN(INSTMACRO)\
   TPETRA_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNode, INSTMACRO )\
   TPETRA_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNode, INSTMACRO )\
   TPETRA_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNode, INSTMACRO )\
   TPETRA_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosCudaWrapperNode, INSTMACRO )


#define TPETRA_INSTANTIATE_LGN(INSTMACRO)\
   TPETRA_INSTANTIATE_LG_call( Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNodeNS, INSTMACRO )\
   TPETRA_INSTANTIATE_LG_call( Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNodeNS, INSTMACRO )\
   TPETRA_INSTANTIATE_LG_call( Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNodeNS, INSTMACRO )\
   TPETRA_INSTANTIATE_LG_call( Exp_Macro_Kokkos_Compat_KokkosCudaWrapperNodeNS, INSTMACRO )


#define TPETRA_INSTANTIATE_TSLG(INSTMACRO)


#define TPETRA_INSTANTIATE_SL(INSTMACRO)\
   INSTMACRO( double , int )\
   INSTMACRO( double , int )\
   INSTMACRO( std_complex0double0 , int )\
   INSTMACRO( std_complex0double0 , int )


#define TPETRA_INSTANTIATE_N(INSTMACRO)\
   Inst_Kokkos_Compat_KokkosSerialWrapperNode(INSTMACRO)\
   Inst_Kokkos_Compat_KokkosThreadsWrapperNode(INSTMACRO)\
   Inst_Kokkos_Compat_KokkosOpenMPWrapperNode(INSTMACRO)\
   Inst_Kokkos_Compat_KokkosCudaWrapperNode(INSTMACRO)


#define TPETRA_INSTANTIATE_SLGN_NOGPU(INSTMACRO)\
   TPETRA_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNode, INSTMACRO )\
   TPETRA_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNode, INSTMACRO )\
   TPETRA_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNode, INSTMACRO )


#define TPETRA_INSTANTIATE_LGN_NOGPU(INSTMACRO)\
   TPETRA_INSTANTIATE_LG_call( Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNodeNS, INSTMACRO )\
   TPETRA_INSTANTIATE_LG_call( Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNodeNS, INSTMACRO )\
   TPETRA_INSTANTIATE_LG_call( Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNodeNS, INSTMACRO )


#define TPETRA_INSTANTIATE_SLGN_NOGPU(INSTMACRO)\
   TPETRA_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNode, INSTMACRO )\
   TPETRA_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNode, INSTMACRO )\
   TPETRA_INSTANTIATE_SLG_call( Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNode, INSTMACRO )


#define TPETRA_INSTANTIATE_LGN_NOGPU(INSTMACRO)\
   TPETRA_INSTANTIATE_LG_call( Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNodeNS, INSTMACRO )\
   TPETRA_INSTANTIATE_LG_call( Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNodeNS, INSTMACRO )\
   TPETRA_INSTANTIATE_LG_call( Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNodeNS, INSTMACRO )


#define TPETRA_INSTANTIATE_SLG_NOGPU(INSTMACRO)\
   TPETRA_INSTANTIATE_SLG(INSTMACRO)


#define TPETRA_INSTANTIATE_LG_NOGPU(INSTMACRO)\
   TPETRA_INSTANTIATE_LG(INSTMACRO)


#define TPETRA_INSTANTIATE_N_NOGPU(INSTMACRO)\
   Inst_Kokkos_Compat_KokkosSerialWrapperNode(INSTMACRO)\
   Inst_Kokkos_Compat_KokkosThreadsWrapperNode(INSTMACRO)\
   Inst_Kokkos_Compat_KokkosOpenMPWrapperNode(INSTMACRO)


#define TPETRA_INSTANTIATE_TSLGN_NOGPU(INSTMACRO)


#define TPETRA_INSTANTIATE_TSLG_NOGPU(INSTMACRO)


#define TPETRA_INSTANTIATE_CONVERT_NOGPU(INSTMACRO)


#define TPETRA_INSTANTIATE_CONVERT_NOGPU_SSL(INSTMACRO)


#define TPETRA_INSTANTIATE_TESTMV(INSTMACRO)\
	INSTMACRO( double , int , int , Kokkos_Compat_KokkosSerialWrapperNode )\
	INSTMACRO( double , int , long , Kokkos_Compat_KokkosSerialWrapperNode )\
	INSTMACRO( std_complex0double0 , int , int , Kokkos_Compat_KokkosSerialWrapperNode )\
	INSTMACRO( std_complex0double0 , int , long , Kokkos_Compat_KokkosSerialWrapperNode )


#define TPETRA_INSTANTIATE_TESTMV_NOGPU(INSTMACRO)


#define TPETRA_INSTANTIATE_DOUBLE_INT_INT_N(INSTMACRO)\
   Exp_Macro_Kokkos_Compat_KokkosSerialWrapperNode(double , int , int , INSTMACRO )\
   Exp_Macro_Kokkos_Compat_KokkosThreadsWrapperNode( double , int , int , INSTMACRO )\
   Exp_Macro_Kokkos_Compat_KokkosOpenMPWrapperNode( double , int , int , INSTMACRO )\
   Exp_Macro_Kokkos_Compat_KokkosCudaWrapperNode( double , int , int , INSTMACRO )


#define TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(INSTMACRO)\
	INSTMACRO( double , int , int , Kokkos_Compat_KokkosSerialWrapperNode )\
	INSTMACRO( double , int , long , Kokkos_Compat_KokkosSerialWrapperNode )\
	INSTMACRO( std_complex0double0 , int , int , Kokkos_Compat_KokkosSerialWrapperNode )\
	INSTMACRO( std_complex0double0 , int , long , Kokkos_Compat_KokkosSerialWrapperNode )


#define TPETRA_INSTANTIATE_SLG_NO_ORDINAL_SCALAR(INSTMACRO)\
   INSTMACRO( double , int , int )\
   INSTMACRO( double , int , long )\
   INSTMACRO( std_complex0double0 , int , int )\
   INSTMACRO( std_complex0double0 , int , long )


#define TPETRA_ETI_MANGLING_TYPEDEFS()\
	TPETRA_ETI_SERIAL_TYPEDEF \
	TPETRA_ETI_THREADS_TYPEDEF \
	TPETRA_ETI_OPENMP_TYPEDEF \
	TPETRA_ETI_CUDA_TYPEDEF \
	typedef std::complex<double> std_complex0double0; 

#endif // TPETRACORE_ETIHELPERMACROS_H

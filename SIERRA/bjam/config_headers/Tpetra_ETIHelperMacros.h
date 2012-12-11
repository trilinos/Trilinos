#ifndef TPETRA_ETIHELPERMACROS_H_
#define TPETRA_ETIHELPERMACROS_H_

#include <Tpetra_ConfigDefs.hpp>

/* Tpetra provides official support for the following nodes */
#include <Kokkos_SerialNode.hpp>
#if defined(HAVE_KOKKOSCLASSIC_TBB)
#  include <Kokkos_TBBNode.hpp>
#endif
#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
#  include <Kokkos_TPINode.hpp>
#endif
#if defined(HAVE_KOKKOSCLASSIC_OPENMP)
#  include <Kokkos_OpenMPNode.hpp>
#endif
#if defined(HAVE_KOKKOSCLASSIC_THRUST)
#  include <Kokkos_ThrustGPUNode.hpp>
#endif

/* Tpetra provides official support for dd_real and qd_real */
#if defined(HAVE_TPETRA_QD)
#include <qd/qd_real.h>
#endif

#define TPETRA_INSTANTIATE_TSLGN(INSTMACRO)


#define TPETRA_INSTANTIATE_CONVERT(INSTMACRO)


#define TPETRA_INSTANTIATE_SLGN(INSTMACRO)


#define TPETRA_INSTANTIATE_LGN(INSTMACRO)


#define TPETRA_INSTANTIATE_TSLG(INSTMACRO)


#define TPETRA_INSTANTIATE_SLG(INSTMACRO)


#define TPETRA_INSTANTIATE_LG(INSTMACRO)


#define TPETRA_INSTANTIATE_N(INSTMACRO)


#define TPETRA_INSTANTIATE_SLGN_NOGPU(INSTMACRO)


#define TPETRA_INSTANTIATE_LGN_NOGPU(INSTMACRO)
  

#define TPETRA_INSTANTIATE_SLG_NOGPU(INSTMACRO)
  

#define TPETRA_INSTANTIATE_LG_NOGPU(INSTMACRO)
   

#define TPETRA_INSTANTIATE_N_NOGPU(INSTMACRO)
    

#define TPETRA_INSTANTIATE_TSLGN_NOGPU(INSTMACRO)


#define TPETRA_INSTANTIATE_TSLG_NOGPU(INSTMACRO)
 

#define TPETRA_INSTANTIATE_CONVERT_NOGPU(INSTMACRO)


#define TPETRA_INSTANTIATE_TESTMV(INSTMACRO)


#define TPETRA_INSTANTIATE_TESTMV_NOGPU(INSTMACRO)


#define TPETRA_INSTANTIATE_SLGNN(INSTMACRO)\
	INSTMACRO( double , int , int , Kokkos_SerialNode , Kokkos_TPINode )\
	INSTMACRO( double , int , int , Kokkos_TPINode , Kokkos_SerialNode )


#define TPETRA_INSTANTIATE_LGNN(INSTMACRO)\
	INSTMACRO( int , int , Kokkos_SerialNode , Kokkos_TPINode )\
	INSTMACRO( int , int , Kokkos_TPINode , Kokkos_SerialNode )


#define TPETRA_ETI_MANGLING_TYPEDEFS()  \
	typedef Kokkos::SerialNode Kokkos_SerialNode; \
	typedef Kokkos::TPINode Kokkos_TPINode;

#endif // TPETRA_ETIHELPERMACROS_H_

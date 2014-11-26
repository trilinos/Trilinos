#ifndef TPETRA_ETIHELPERMACROS_H_
#define TPETRA_ETIHELPERMACROS_H_

#include <Tpetra_ConfigDefs.hpp>

/* Tpetra provides official support for the following nodes */
#include <Kokkos_DefaultNode.hpp>

/* Tpetra provides official support for dd_real and qd_real */
#if defined(HAVE_TPETRA_QD)
#include <qd/qd_real.h>
#endif

#define TPETRA_INSTANTIATE_TSLGN(INSTMACRO)


#define TPETRA_INSTANTIATE_CONVERT(INSTMACRO)


#define TPETRA_INSTANTIATE_VECTOR(INSTMACRO)\
	INSTMACRO( double , int , int , KokkosClassic_SerialNode )\
	INSTMACRO( double , int , int , KokkosClassic_TPINode )


#define TPETRA_INSTANTIATE_SLGN(INSTMACRO)\
	INSTMACRO( double , int , int , KokkosClassic_SerialNode )\
	INSTMACRO( double , int , int , KokkosClassic_TPINode )


#define TPETRA_INSTANTIATE_PLGN(INSTMACRO)\
	INSTMACRO( double , int , int , KokkosClassic_SerialNode )\
	INSTMACRO( double , int , int , KokkosClassic_TPINode )\
	INSTMACRO( int , int , int , KokkosClassic_SerialNode )\
	INSTMACRO( int , int , int , KokkosClassic_TPINode )


#define TPETRA_INSTANTIATE_LGN(INSTMACRO)\
	INSTMACRO( int , int , KokkosClassic_SerialNode )\
	INSTMACRO( int , int , KokkosClassic_TPINode )


#define TPETRA_INSTANTIATE_TSLG(INSTMACRO)


#define TPETRA_INSTANTIATE_SLG(INSTMACRO)\
	INSTMACRO( double , int , int )


#define TPETRA_INSTANTIATE_LG(INSTMACRO)\
	INSTMACRO( int , int )


#define TPETRA_INSTANTIATE_SL(INSTMACRO)\
	INSTMACRO( double , int )


#define TPETRA_INSTANTIATE_N(INSTMACRO)\
	INSTMACRO( KokkosClassic_SerialNode )\
	INSTMACRO( KokkosClassic_TPINode )


#define TPETRA_INSTANTIATE_SLGN_NOGPU(INSTMACRO)\
	INSTMACRO( double , int , int , KokkosClassic_SerialNode )\
	INSTMACRO( double , int , int , KokkosClassic_TPINode )


#define TPETRA_INSTANTIATE_LGN_NOGPU(INSTMACRO)\
	INSTMACRO( int , int , KokkosClassic_SerialNode )\
	INSTMACRO( int , int , KokkosClassic_TPINode )


#define TPETRA_INSTANTIATE_SLG_NOGPU(INSTMACRO)\
	INSTMACRO( double , int , int )


#define TPETRA_INSTANTIATE_LG_NOGPU(INSTMACRO)\
	INSTMACRO( int , int )


#define TPETRA_INSTANTIATE_N_NOGPU(INSTMACRO)\
	INSTMACRO( KokkosClassic_SerialNode )\
	INSTMACRO( KokkosClassic_TPINode )


#define TPETRA_INSTANTIATE_TSLGN_NOGPU(INSTMACRO)


#define TPETRA_INSTANTIATE_TSLG_NOGPU(INSTMACRO)


#define TPETRA_INSTANTIATE_CONVERT_NOGPU(INSTMACRO)


#define TPETRA_INSTANTIATE_TESTMV(INSTMACRO)\
	INSTMACRO( double , int , int , KokkosClassic_SerialNode )\
	INSTMACRO( double , int , int , KokkosClassic_TPINode )


#define TPETRA_INSTANTIATE_TESTMV_NOGPU(INSTMACRO)\
	INSTMACRO( double , int , int , KokkosClassic_SerialNode )\
	INSTMACRO( double , int , int , KokkosClassic_TPINode )


#define TPETRA_INSTANTIATE_DOUBLE_INT_INT_N(INSTMACRO)\
	INSTMACRO( double , int , int , KokkosClassic_SerialNode )\
	INSTMACRO( double , int , int , KokkosClassic_TPINode )


#define TPETRA_ETI_MANGLING_TYPEDEFS()  \
	typedef KokkosClassic::SerialNode KokkosClassic_SerialNode; \
	typedef KokkosClassic::TPINode KokkosClassic_TPINode;

#endif // TPETRA_ETIHELPERMACROS_H_

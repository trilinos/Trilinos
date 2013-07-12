#ifndef IFPACK2_ETIHELPERMACROS_H_
#define IFPACK2_ETIHELPERMACROS_H_

#include <Ifpack2_ConfigDefs.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>

#define IFPACK2_INSTANTIATE_SLG(INSTMACRO)\
	INSTMACRO( double , int , int )\
	INSTMACRO( double , int , long )


#define IFPACK2_INSTANTIATE_SLG_REAL(INSTMACRO)\
	INSTMACRO( double , int , int )\
	INSTMACRO( double , int , long )


#define IFPACK2_INSTANTIATE_LG(INSTMACRO)\
	INSTMACRO( int , int )\
	INSTMACRO( int , long )


#define IFPACK2_ETI_MANGLING_TYPEDEFS() 

#endif // IFPACK2_ETIHELPERMACROS_H_

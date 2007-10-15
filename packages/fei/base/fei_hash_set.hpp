#ifndef _fei_hash_set_hpp_
#define _fei_hash_set_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

//
//Some fei implementation code can use a hash-set. That code refers
//to the hash-set type using the macro 'FEI_HASH_SET'. That macro
//should be defined in this header to the platform-specific type that
//provides a hash-set implementation.
//

////////////GNUC case//////////////////
#ifdef __GNUC__
#include <ext/hash_set>

#if !defined(__INTEL_COMPILER)
#define FEI_HASH_SET __gnu_cxx::hash_set
#endif

#endif
///////end of GNUC case///////////////////

////////sgi case///////////////////////
//#ifdef __sgi
//#include <hash_set>
//#define FEI_HASH_SET std::hash_set
//#endif
/////////end of sgi case//////////////
#endif

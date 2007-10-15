#ifndef _fei_hash_map_hpp_
#define _fei_hash_map_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

//
//Some fei implementation code can use a hash-map. That code refers
//to the hash-map type using the macro 'FEI_HASH_MAP'. That macro
//should be defined in this header to the platform-specific type that
//provides a hash-map implementation.
//

//////////////////////////////
//GNUC case covers g++ as well as intel compiler...

#if defined(__GNUC__)
#include <ext/hash_map>

#if !defined(__INTEL_COMPILER)
#define FEI_HASH_MAP __gnu_cxx::hash_map
#endif

#endif

//end of GNUC case
////////////////////////////


////////////////////////////
//sgi case
//#ifdef __sgi
//#include <hash_map>
//#define FEI_HASH_MAP std::hash_map
//#endif

//end of sgi case
/////////////////////////////


#endif

#ifndef __ESI_basicTypes_h
#define __ESI_basicTypes_h

/** This header is compatible with ANSI C or C++ compilers.
    This defines some basic infrastructure such as:

- 'bool' if ESI_SIMULATE_BOOL is defined
- esi::ErrorCode which is the return-value of all ESI interface functions

*/

namespace esi {

/** ESI return codes

    \verbatim
    All ESI functions returning an int error code follow this
    convention:
    
       N == 0: all ok.
       N < 0: some error condition logically fatal has occured.
       N > 0: some informational condition.

    The non-zero values of N are (for ESI compliance purposes)
    significant only in the sign bit. Individual implementations may 
    assign significance to bits other than the sign bit, for instance
    to carry an enum.
    \endverbatim
*/

typedef int ErrorCode;

#ifdef __cplusplus

#ifdef ESI_SIMULATE_BOOL
//
// Following are some #defines that will provide bool support
// when using a compiler that doesn't have a native 'bool' type.
//

#ifdef bool
#undef bool
#endif
#ifdef true
#undef true
#endif
#ifdef false
#undef false
#endif

#define bool int
#define true 1
#define false 0

#endif /* ESI_SIMULATE_BOOL */

#endif /* __cplusplus */

};     /* esi namespace */

#endif /* __ESI_basicTypes_h */


/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_util_util_NamedPair_hpp
#define stk_util_util_NamedPair_hpp

#define NAMED_PAIR( STRUCT_NAME , FIRST_TYPE , FIRST_NAME , SECOND_TYPE , SECOND_NAME )	\
struct STRUCT_NAME {	\
  typedef FIRST_TYPE first_type; \
  typedef SECOND_TYPE second_type; \
  \
  first_type  FIRST_NAME ; \
  second_type SECOND_NAME ; \
  \
  STRUCT_NAME ( const first_type & arg_ ## FIRST_NAME , \
                const second_type & arg_ ## SECOND_NAME ) \
   : FIRST_NAME ( arg_ ## FIRST_NAME ) , \
     SECOND_NAME ( arg_ ## SECOND_NAME ) {} \
  \
  ~ STRUCT_NAME () {} \
  STRUCT_NAME () {} \
  STRUCT_NAME ( const STRUCT_NAME & arg_rhs ) \
   : FIRST_NAME ( arg_rhs. FIRST_NAME ) , \
     SECOND_NAME ( arg_rhs. SECOND_NAME ) {} \
  STRUCT_NAME & operator = ( const STRUCT_NAME & arg_rhs ) \
   { FIRST_NAME  = arg_rhs. FIRST_NAME ; \
     SECOND_NAME = arg_rhs. SECOND_NAME ; \
     return *this ; } \
}; \
\
inline \
bool operator == ( const STRUCT_NAME & arg_lhs , \
                   const STRUCT_NAME & arg_rhs ) \
  { return arg_lhs. FIRST_NAME  == arg_rhs. FIRST_NAME && \
           arg_lhs. SECOND_NAME == arg_rhs. SECOND_NAME ; } \
\
inline \
bool operator != ( const STRUCT_NAME & arg_lhs , \
                   const STRUCT_NAME & arg_rhs ) \
  { return arg_lhs. FIRST_NAME  != arg_rhs. FIRST_NAME || \
           arg_lhs. SECOND_NAME != arg_rhs. SECOND_NAME ; } \
\
inline \
bool operator < ( const STRUCT_NAME & arg_lhs , \
                  const STRUCT_NAME & arg_rhs ) \
  { return arg_lhs. FIRST_NAME < arg_rhs. FIRST_NAME || \
           ( ! ( arg_rhs. FIRST_NAME < arg_lhs. FIRST_NAME ) && \
           arg_lhs. SECOND_NAME < arg_rhs. SECOND_NAME ) ; } \
\
inline \
bool operator > ( const STRUCT_NAME & arg_lhs , \
                  const STRUCT_NAME & arg_rhs ) \
  { return operator < ( arg_rhs , arg_lhs ); } \

#endif


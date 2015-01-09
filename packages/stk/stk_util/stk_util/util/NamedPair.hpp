// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

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


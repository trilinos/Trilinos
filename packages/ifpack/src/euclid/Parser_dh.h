/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef PARSER_DH_DH
#define PARSER_DH_DH

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "euclid_common.h"

#ifdef __cplusplus
extern "C"
{
#endif

  extern void Parser_dhCreate (Parser_dh * p);
  extern void Parser_dhDestroy (Parser_dh p);

  extern bool Parser_dhHasSwitch (Parser_dh p, char *in);
  extern bool Parser_dhReadString (Parser_dh p, char *in, char **out);
  extern bool Parser_dhReadInt (Parser_dh p, char *in, int *out);
  extern bool Parser_dhReadDouble (Parser_dh p, char *in, double *out);
  /* if the flag (char *in) is found, these four return 
     true and set "out" accordingly.  If not found, they return 
     false, and "out" is unaltered.
   */

  extern void Parser_dhPrint (Parser_dh p, FILE * fp, bool allPrint);
  /* Outputs all <flag,value> pairs.  "bool allPrint" is
   * only meaningful when Euclid is compiled in MPI mode
   */

  extern void Parser_dhInsert (Parser_dh p, char *name, char *value);
  /* For inserting a new <flag,value> pair, or altering
   * the value of an existing pair from within user apps.
   */

  extern void Parser_dhUpdateFromFile (Parser_dh p, char *name);

  extern void Parser_dhInit (Parser_dh p, int argc, char *argv[]);
  /* Init enters <flag,value> pairs in its internal database in
     the following order:

     (1)   $PCPACK_DIR/options_database  
     (2)   "database" in local directory, if the file exists
     (3)   "pathname/foo" if argv[] contains a pair of entries:
     -db_filename pathname/foo
     (4)   flag,value pairs from the command line (ie, argv)

     If a flag already exists, its value is updated if it is
     encountered a second time.  

     WARNING! to enter a negative value, you must use two dashes, e.g:
     -myvalue  --0.1
     otherwise, if you code "-myvalue -0.1" you will have entered
     the pair of entries <-myvalue, 1>, <-0.1, 1>.  Yuck!@#
     But this works, since Euclid doesn't use negative numbers much.

     If the 2nd entry is missing, a value of "1" is assumed (this only
     works on the command line; for files, you must explicitly code a
     value.  See $PCPACK_DIR/options_database for samples.

     The following will cause Parser_dhHasSwitch(p, "-phoo") to return false:
     -phoo 0
     -phoo false
     -phoo False
     -phoo FALSE
     any other value, including something silly like -phoo 0.0
     will return true.
   */

#ifdef __cplusplus
}
#endif
#endif

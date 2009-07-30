/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef PARSER_DH_DH
#define PARSER_DH_DH

#include "euclid_common.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void Parser_dhCreate(Parser_dh *p);
extern void Parser_dhDestroy(Parser_dh p);

extern bool Parser_dhHasSwitch(Parser_dh p, char *in);
extern bool Parser_dhReadString(Parser_dh p, char *in, char **out);
extern bool Parser_dhReadInt(Parser_dh p, char *in, int *out);
extern bool Parser_dhReadDouble(Parser_dh p, char *in, double *out);
  /* if the flag (char *in) is found, these four return 
     true and set "out" accordingly.  If not found, they return 
     false, and "out" is unaltered.
   */

extern void Parser_dhPrint(Parser_dh p, FILE *fp, bool allPrint);
  /* Outputs all <flag,value> pairs.  "bool allPrint" is
   * only meaningful when Euclid is compiled in MPI mode
   */

extern void Parser_dhInsert(Parser_dh p, char *name, char *value);
  /* For inserting a new <flag,value> pair, or altering
   * the value of an existing pair from within user apps.
   */

extern void Parser_dhUpdateFromFile(Parser_dh p, char *name);

extern void Parser_dhInit(Parser_dh p, int argc, char *argv[]);
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


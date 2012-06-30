/*
 * Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
/* ifdefx.h - ifdef file for cgi shell routines
 * This file is used to define the system dependent ways C is
 * called from FORTRAN.  Underscores are used by default.
 *
 * SUN DEC/ULTRIX ALLIANT : C routines must have underscores
 * SGI CONVEX             : C routines must have underscores
 *
 * VAX HP IBM/aix         : C routines do not have underscores
 *
 * CRAY/UNICOS            : C routines must be capitalized,
 *                            and no underscores
 *
 * This include file is used by MDCGI.C and all CGISxxx.C files.
 */

#if defined(ADDC_)
#endif
#if !defined(CRA) && !defined(ADDC_) && !defined(COUGAR)
#define xcoon_ xcoon
#define xcooff_ xcooff
#define xcact_ xcact
#define xcdact_ xcdact
#define xcsol_ xcsol
#define cgia60_ cgia60
#define cgi16c_ cgi16c
#define cgi35c_ cgi35c
#define cgi810_ cgi810
#define cgi48l_ cgi48l
#define cgi24l_ cgi24l
#define cgi35a_ cgi35a
#define cgimet_ cgimet
#define cginws_ cginws
#define cgipst_ cgipst
#define cgiqms_ cgiqms
#define cgit05_ cgit05
#define cgit07_ cgit07
#define cgit15_ cgit15
#define cgitk4_ cgitk4
#define cgiv34_ cgiv34
#define cgix11_ cgix11
#define cgihcb_ cgihcb
#define cgif8t_ cgif8t
#define cgif3c_ cgif3c
#define cgifsq_ cgifsq
#endif
#if defined(CRA)
#define xcoon_ XCOON
#define xcooff_ XCOOFF
#define xcact_ XCACT
#define xcdact_ XCDACT
#define xcsol_ XCSOL
#define cgi16c_ CGI16C
#define cgi35c_ CGI35C
#define cgi48l_ CGI48L
#define cgi24l_ CGI24L
#define cgi35a_ CGI35A
#define cgi810_ CGI810
#define cgia60_ CGIA60
#define cgimet_ CGIMET
#define cginws_ CGINWS
#define cgipst_ CGIPST
#define cgiqms_ CGIQMS
#define cgit05_ CGIT05
#define cgit07_ CGIT07
#define cgit15_ CGIT15
#define cgitk4_ CGITK4
#define cgiv34_ CGIV34
#define cgix11_ CGIX11
#define cgihcb_ CGIHCB
#define cgif8t_ CGIF8T
#define cgif3c_ CGIF3C
#define cgifsq_ CGIFSQ
#endif

/* end ifdefx.h */

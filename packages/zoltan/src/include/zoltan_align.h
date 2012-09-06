/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */

#ifndef __ZOLTAN_ALIGN_H
#define __ZOLTAN_ALIGN_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* 
 * Plauger alignment algorithm, The Standard C Library.
 * Forces malloc'ed variable size struct alignment.
 * ZOLTAN_ALIGN_VAL is defined in Zoltan/include/zoltan_align.h;
 * values are 0,1,3,7U depending upon machine.
 * (E.g., 7U == 8 byte alignment.)
 */

#ifndef ZOLTAN_ALIGN_VAL
#define ZOLTAN_ALIGN_VAL 7U
#endif

extern int Zoltan_Align(int);
extern size_t Zoltan_Align_size_t(size_t);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif

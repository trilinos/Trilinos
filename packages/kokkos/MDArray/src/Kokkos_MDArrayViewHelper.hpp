/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#ifndef KOKKOS_MDARRAYVIEWHELPER
#define KOKKOS_MDARRAYVIEWHELPER


#define MULTI_INDEX_LEFT_8( i0 , i1 , i2 , i3 , i4 , i5 , i6 , i7 , dim ) \
  ( i0 + dim[0] * ( i1 + dim[1] * ( i2 + dim[2] * ( i3 + dim[3] * ( i4 + dim[4] * ( i5 + dim[5] * ( i6 * dim[6] * ( i7 ))))))))

#define MULTI_INDEX_LEFT_7( i0 , i1 , i2 , i3 , i4 , i5 , i6 , dim ) \
  ( i0 + dim[0] * ( i1 + dim[1] * ( i2 + dim[2] * ( i3 + dim[3] * ( i4 + dim[4] * ( i5 + dim[5] * ( i6 )))))))

#define MULTI_INDEX_LEFT_6( i0 , i1 , i2 , i3 , i4 , i5 , dim ) \
  ( i0 + dim[0] * ( i1 + dim[1] * ( i2 + dim[2] * ( i3 + dim[3] * ( i4 + dim[4] * ( i5 ))))))

#define MULTI_INDEX_LEFT_5( i0 , i1 , i2 , i3 , i4 , dim ) \
  ( i0 + dim[0] * ( i1 + dim[1] * ( i2 + dim[2] * ( i3 + dim[3] * ( i4 )))))

#define MULTI_INDEX_LEFT_4( i0 , i1 , i2 , i3 , dim ) \
  ( i0 + dim[0] * ( i1 + dim[1] * ( i2 + dim[2] * ( i3 ))))

#define MULTI_INDEX_LEFT_3( i0 , i1 , i2 , dim ) \
  ( i0 + dim[0] * ( i1 + dim[1] * ( i2 )))

#define MULTI_INDEX_LEFT_2( i0 , i1 , dim ) \
  ( i0 + dim[0] * ( i1 ))

#define MULTI_INDEX_LEFT_1( i0 , dim ) \
  ( i0 )

#endif /* KOKKOS_MDARRAYVIEWHELPER */


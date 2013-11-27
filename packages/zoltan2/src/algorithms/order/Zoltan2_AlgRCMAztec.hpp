// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef _ZOLTAN2_ALGRCMAZTEC_HPP_
#define _ZOLTAN2_ALGRCMAZTEC_HPP_

#if 0 // DOESN'T COMPILE FOR NOW!
#include "az_aztec.h"


////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgRCMAztec.hpp
//! \brief Interface to Aztec's RCM, which is based on George & Liu.

namespace Zoltan2{

  void RCMAztec(int &perm, int N)
  // TODO: Pass CRS graph
  // TODO: What about inverse perm?
  {
   int i, root, total, ccsize;
   int mask[N];

   // Convert matrix to Fortran 1-base
   for (i = N+1 ; i < bindx2[N]; i++ ) bindx2[i]++;
   for (i = 0 ; i <= N ; i++ )         bindx2[i] -= N;

   // Call RCM for each connected component
   for (i = 0 ; i < N ; i++ ) mask[i] = 1;
   root = 1;
   while (total != N ) {
      AZ_FNROOT_F77(&root,bindx2,&(bindx2[N+1]),mask, &nlvl,
              &((*ordering)[total]), *inv_ordering);
      AZ_RCM_F77(&root,bindx2,&(bindx2[N+1]),mask,&((*ordering)[total]),
           &ccsize, *inv_ordering);

      if ( ccsize != N) {
         for (i = 0 ; i < ccsize ; i++) mask[(*ordering)[total+i]-1] = 0;
         for (i = 0 ; i < N ; i++ ) {
            if ( mask[i] == 1) break;
         }
         root = i+1;
      }
      total += ccsize;
      if (ccsize == 0) {
         // Error; throw exception
      }
   }

   // Convert perm back to 0-base
   for (i = 0 ; i < N ; i++ ) 
     (*ordering)[i]--;

  }

}
#endif
#endif

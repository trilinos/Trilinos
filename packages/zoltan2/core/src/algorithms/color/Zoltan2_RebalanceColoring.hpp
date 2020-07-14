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
#ifndef _ZOLTAN2_REBALANCECOLORING_HPP_
#define _ZOLTAN2_REBALANCECOLORING_HPP_

#include <Zoltan2_Standards.hpp>

// Given a valid coloring of a graph, rebalance the colors 
// such that either:
// (a) color classes are as balanced as possible, or
// (b) every color class has at least minSize vertices.
// This function can be called as a post-processing after any initial
// coloring.
// This is a greedy heuristic so there is no guarantee of success,
// though in practice we believe it will work well.
   int rebalanceColoring(
     const lno_t nVtx,
     ArrayView<const lno_t> edgeIds,
     ArrayView<const offset_t> offsets,
     ArrayRCP<int> colors,
     const int balanceColors,
     const lno_t minSize
     )
   {
#ifndef INCLUDE_ZOLTAN2_EXPERIMENTAL
    Z2_THROW_EXPERIMENTAL("Zoltan2 rebalanceColoring is experimental and not tested");
#else
   // Experimental: Not tested yet!
   //
   // Count size of each color class.
   // No vertex weights, may add in future.
   lno_t maxColor = 0;
   Teuchos::Array<lno_t> colorSize(nVtx,0);
   for (lno_t i=0; i < nVtx; i++){
     if (colors[i] > maxColor) 
       maxColor = colors[i];
   }
   for (lno_t i=0; i < nVtx; i++){
     colorSize[colors[i]]++;
   }
   lno_t targetSize = 0;
   if (balanceColors > 0)
     targetSize = nVtx/maxColor;
   else
     targetSize = minSize;

   // Vertex-centric version: Loop over vertices
   // and move it to different color if needed.
   Teuchos::Array<int> forbidden(maxDegree+2, 0);
   for (lno_t i=0; i < nVtx; i++){
     if (colorSize[colors[i]] > targetSize){
       // Find first available color that is not full yet
       for (offset_t j=offsets[i]; j<offsets[i+1]; j++){
         lno_t nbor = edgeIds[j];
         if (colors[nbor] > 0){
           // Neighbors' colors are forbidden
           forbidden[colors[nbor]] = i;
         }
       }
       // Pick first (smallest) underfull color > 0.
       // If no such color, just keep colors[i].
       int newcolor = colors[i];
       for (int c=1; c <= maxColor; c++){
         if ((forbidden[c] != i) && (colorSize[c]<targetSize)){ 
           newcolor = c;
           break;
         }
       }
       
       // Move vertex i to underfull color.
       colorSize[colors[i]]--;
       colorSize[newcolor]++;
       colors[i] = newcolor;
     } 
   }

#endif
   }

#endif

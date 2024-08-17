// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

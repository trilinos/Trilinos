// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file   Intrepid_CellToolsDef.hpp
    \brief  Definition file for the Intrepid2::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/
#ifndef __INTREPID2_CELLTOOLS_DEF_HPP__
#define __INTREPID2_CELLTOOLS_DEF_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  //============================================================================================//
  //                                                                                            //
  //                                           Debug                                            //
  //                                                                                            //
  //============================================================================================//


  template<class Scalar>
  void CellTools<Scalar>::printSubcellVertices(const int subcellDim,
                                               const int subcellOrd,
                                               const shards::CellTopology & parentCell){
  
    // Get number of vertices for the specified subcell and parent cell dimension
    int subcVertexCount = parentCell.getVertexCount(subcellDim, subcellOrd);
    int cellDim         = parentCell.getDimension();
  
    // Allocate space for the subcell vertex coordinates
    FieldContainer<double> subcellVertices(subcVertexCount, cellDim);
  
    // Retrieve the vertex coordinates
    getReferenceSubcellVertices(subcellVertices,
                                subcellDim,
                                subcellOrd,
                                parentCell);
  
    // Print the vertices
    std::cout 
      << " Subcell " << std::setw(2) << subcellOrd 
      <<  " is " << parentCell.getName(subcellDim, subcellOrd) << " with vertices = {";
  
    // Loop over subcell vertices
    for(int subcVertOrd = 0; subcVertOrd < subcVertexCount; subcVertOrd++){
      std::cout<< "(";
    
      // Loop over vertex Cartesian coordinates
      for(int dim = 0; dim < (int)parentCell.getDimension(); dim++){
        std::cout << subcellVertices(subcVertOrd, dim);
        if(dim < (int)parentCell.getDimension()-1 ) { std::cout << ","; }
      }
      std::cout<< ")";
      if(subcVertOrd < subcVertexCount - 1) { std::cout << ", "; }
    }
    std::cout << "}\n";
  }
  

  template<class Scalar>
  template<class ArrayCell>
  void CellTools<Scalar>::printWorksetSubcell(const ArrayCell &             cellWorkset,
                                              const shards::CellTopology &  parentCell,
                                              const int&                    pCellOrd,
                                              const int&                    subcellDim,
                                              const int&                    subcellOrd,
                                              const int&                    fieldWidth){
  
    // Get the ordinals, relative to reference cell, of subcell cellWorkset
    int subcNodeCount = parentCell.getNodeCount(subcellDim, subcellOrd);
    int pCellDim      = parentCell.getDimension();
    std::vector<int> subcNodeOrdinals(subcNodeCount);
  
    for(int i = 0; i < subcNodeCount; i++){
      subcNodeOrdinals[i] = parentCell.getNodeMap(subcellDim, subcellOrd, i);
    }
  
    // Loop over parent cells and print subcell cellWorkset
  
    std::cout 
      << " Subcell " << subcellOrd << " on parent cell " << pCellOrd << " is " 
      << parentCell.getName(subcellDim, subcellOrd) << " with node(s) \n ({";
  
    for(int i = 0; i < subcNodeCount; i++){
    
      // print Cartesian coordinates of the node
      for(int dim = 0; dim < pCellDim; dim++){
        std::cout
          << std::setw(fieldWidth) << std::right << cellWorkset(pCellOrd, subcNodeOrdinals[i], dim); 
        if(dim < pCellDim - 1){ std::cout << ","; }
      }
      std::cout << "}";
      if(i < subcNodeCount - 1){ std::cout <<", {"; }
    }
    std::cout << ")\n\n";
  }

}

#endif

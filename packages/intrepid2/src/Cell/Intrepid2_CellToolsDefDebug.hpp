// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file   Intrepid2_CellToolsDefDebug.hpp
    \brief  Definition file for the debug functions of the Intrepid2::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/
#ifndef __INTREPID2_CELLTOOLS_DEBUG_DEF_HPP__
#define __INTREPID2_CELLTOOLS_DEBUG_DEF_HPP__

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
  void CellTools<Scalar>::printSubcellVertices(const ordinal_type subcellDim,
                                               const ordinal_type subcellOrd,
                                               const shards::CellTopology & parentCell){
  
    // Get number of vertices for the specified subcell and parent cell dimension
    ordinal_type subcVertexCount = parentCell.getVertexCount(subcellDim, subcellOrd);
    ordinal_type cellDim         = parentCell.getDimension();
  
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
    for(ordinal_type subcVertOrd = 0; subcVertOrd < subcVertexCount; subcVertOrd++){
      std::cout<< "(";
    
      // Loop over vertex Cartesian coordinates
      for(ordinal_type dim = 0; dim < (ordinal_type)parentCell.getDimension(); dim++){
        std::cout << subcellVertices(subcVertOrd, dim);
        if(dim < (ordinal_type)parentCell.getDimension()-1 ) { std::cout << ","; }
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
                                              const ordinal_type&                    pCellOrd,
                                              const ordinal_type&                    subcellDim,
                                              const ordinal_type&                    subcellOrd,
                                              const ordinal_type&                    fieldWidth){
  
    // Get the ordinals, relative to reference cell, of subcell cellWorkset
    ordinal_type subcNodeCount = parentCell.getNodeCount(subcellDim, subcellOrd);
    ordinal_type pCellDim      = parentCell.getDimension();
    std::vector<ordinal_type> subcNodeOrdinals(subcNodeCount);
  
    for(ordinal_type i = 0; i < subcNodeCount; i++){
      subcNodeOrdinals[i] = parentCell.getNodeMap(subcellDim, subcellOrd, i);
    }
  
    // Loop over parent cells and print subcell cellWorkset
  
    std::cout 
      << " Subcell " << subcellOrd << " on parent cell " << pCellOrd << " is " 
      << parentCell.getName(subcellDim, subcellOrd) << " with node(s) \n ({";
  
    for(ordinal_type i = 0; i < subcNodeCount; i++){
    
      // print Cartesian coordinates of the node
      for(ordinal_type dim = 0; dim < pCellDim; dim++){
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

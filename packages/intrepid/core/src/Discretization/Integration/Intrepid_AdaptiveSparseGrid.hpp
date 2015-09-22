// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_AdaptiveSparseGrid.hpp
    \brief  Header file for the Intrepid::AdaptiveSparseGrid class.
    \author Created by D. Kouri and D. Ridzal.
*/

#ifndef INTREPID_ADAPTIVESPARSEGRID_HPP
#define INTREPID_ADAPTIVESPARSEGRID_HPP

#include "Intrepid_AdaptiveSparseGridInterface.hpp"

namespace Intrepid {

/** \class Intrepid::AdaptiveSparseGrid
    \brief Builds general adaptive sparse grid rules (Gerstner and Griebel) 
           using the 1D cubature rules in the Intrepid::CubatureLineSorted 
	   class.
*/

template<class Scalar, class UserVector>
class AdaptiveSparseGrid {
public:
  /** \brief Update adaptive sparse grid.
      
      \param indexSet      [in/out] - Admissible set of multi-indices.
      \param integralValue [in/out] - Value of integral.
      \param problem_data  [in]     - User defined problem data.
  */
  static Scalar refine_grid(
              typename std::multimap<Scalar,std::vector<int> > & indexSet,      
	      UserVector & integralValue,
	      AdaptiveSparseGridInterface<Scalar,UserVector> & problem_data); 

  /** \brief Update adaptive sparse grid.
      
      \param activeIndex          [in/out] - Active Indices
      \param oldIndex             [in/out] - Inactive Indices
      \param integralValue        [in/out] - Value of Integral
      \param globalErrorIndicator [in/out] - Error Indicator
      \param problem_data         [in]     - User defined problem data.
  */
  static Scalar refine_grid(
	      typename std::multimap<Scalar,std::vector<int> > & activeIndex, 
	      std::set<std::vector<int> > & oldIndex, 
	      UserVector & integralValue,
	      Scalar globalErrorIndicator,
	      AdaptiveSparseGridInterface<Scalar,UserVector> & problem_data); 
 
  /** \brief Update adaptive sparse grid.
      
      \param activeIndex          [in/out] - Active Indices
      \param oldIndex             [in/out] - Inactive Indices
      \param integralValue        [in/out] - Value of Integral
      \param cubRule              [in/out] - Sparse Grid Points and Weights
      \param globalErrorIndicator [in/out] - Error Indicator
      \param problem_data         [in]     - User defined problem data.
  */ 
  static Scalar refine_grid(
	      typename std::multimap<Scalar,std::vector<int> > & activeIndex, 
	      std::set<std::vector<int> > & oldIndex, 
	      UserVector & integralValue,
	      CubatureTensorSorted<Scalar> & cubRule,
	      Scalar globalErrorIndicator,
	      AdaptiveSparseGridInterface<Scalar,UserVector> & problem_data); 
  /*
  static void coarsen_grid(
              typename std::multimap<Scalar,std::vector<int> > & indexSet,
	      int dimension,
	      std::vector<EIntrepidBurkardt> rule1D, 
	      std::vector<EIntrepidGrowth> growth1D);

  static void coarsen_grid(
              typename std::multimap<Scalar,std::vector<int> > & activeIndex, 
	      std::set<std::vector<int> > oldIndex,
	      int dimension;
	      std::vector<EIntrepidBurkardt> rule1D, 
	      std::vector<EIntrepidGrowth> growth1D);
  */

  /** \brief Given an index, build the corresponding differential cubature rule.
      
      \param outRule       [out] - Cubature nodes/weights for differential rule.
      \param index         [in]  - Multi-index of cubature levels.
      \param problem_data  [in]  - User defined problem data.
  */
  static void build_diffRule(
	      CubatureTensorSorted<Scalar> & outRule, 
	      std::vector<int> index, 
	      AdaptiveSparseGridInterface<Scalar,UserVector> & problem_data);

  /** \brief Given an index, build the corresponding differential cubature rule.
      
      \param outRule       [out] - Cubature nodes/weights for differential rule.
      \param index         [in]  - Multi-index of cubature levels.
      \param dimension     [in]  - Dimension of integration domain.
      \param rule1D        [in]  - 1D cubature rule names.
      \param growth1D      [in]  - 1D cubature growth rules.
      \param problem_data  [in]  - User defined problem data.
  */
  static void build_diffRule(
	      CubatureTensorSorted<Scalar> & outRule, 
	      std::vector<int> index, 
	      int dimension,
	      std::vector<EIntrepidBurkardt> rule1D, 
	      std::vector<EIntrepidGrowth> growth1D,
	      bool isNormalized);
 
   /** \brief Check admissibility of an index set, outputs true if admissible.
      
      \param index         [in]  - Multi-index of cubature levels.
      \param direction     [in]  - Search direction.
      \param inOldIndex    [in]  - Input index set.
      \param problem_data  [in]  - User defined problem data.
  */
  static bool isAdmissible(
	      std::vector<int> index,
	      int direction, 
	      std::set<std::vector<int> > inOldIndex,
	      AdaptiveSparseGridInterface<Scalar,UserVector> & problem_data);
 
  /** \brief Build a classic isotropic sparse grid.
      
      \param output        [out] - Cubature points/weights for sparse grid
      \param dimension     [in]  - Dimension of integration domain, dimension = 2,3,4,5.
      \param maxlevel      [in]  - Maximum level of sparse grid.
      \param rule1D        [in]  - 1D cubature rule names.
      \param growth1D      [in]  - 1D cubature growth rules.
      \param problem_data  [in]  - User defined problem data.
  */
  static void buildSparseGrid(
	      CubatureTensorSorted<Scalar> & output,
	      int dimension, int maxlevel,
	      std::vector<EIntrepidBurkardt> rule1D, 
	      std::vector<EIntrepidGrowth> growth1D,
	      bool isNormalized);
					 
 
  
};

} // End Intrepid namespace

// include templated definitions
#include <Intrepid_AdaptiveSparseGridDef.hpp>

#endif

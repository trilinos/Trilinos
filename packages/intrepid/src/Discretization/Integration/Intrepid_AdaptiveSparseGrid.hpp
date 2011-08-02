// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
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

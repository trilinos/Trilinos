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

/** \file   Intrepid_AdaptiveSparseGridDef.hpp
    \brief  Definition file for the Intrepid::AdaptiveSparseGrid class.
    \author Created by D. Kouri and D. Ridzal.
*/

namespace Intrepid {

template<class Scalar, class UserVector>
bool AdaptiveSparseGrid<Scalar,UserVector>::isAdmissible(
               std::vector<int> index, 
	       int direction, 
	       std::set<std::vector<int> > inOldIndex,
	       AdaptiveSparseGridInterface<Scalar,UserVector> & problem_data) {
  /*
    Check if inOldIndex remains admissible if index is added, i.e.
       index-ek in inOldIndex for all k=1,...,dim.
  */
  int dimension = problem_data.getDimension();
  for (int i=0; i<dimension; i++) {
    if (index[i]>1 && i!=direction) {
      index[i]--;
      if (!inOldIndex.count(index)) {
	return false;
      }
      index[i]++;
    }
  }
  return true;
}

template<class Scalar, class UserVector>
void AdaptiveSparseGrid<Scalar,UserVector>::build_diffRule(
	       CubatureTensorSorted<Scalar> & outRule, 
	       std::vector<int> index,
	       AdaptiveSparseGridInterface<Scalar,UserVector> & problem_data) {
	       
  int numPoints     = 0;
  int dimension     = problem_data.getDimension();
  bool isNormalized = problem_data.isNormalized();
  std::vector<EIntrepidBurkardt> rule1D; problem_data.getRule(rule1D);
  std::vector<EIntrepidGrowth> growth1D; problem_data.getGrowth(growth1D);

  for (int i=0; i<dimension; i++) {
    // Compute 1D rules
    numPoints = growthRule1D(index[i],growth1D[i],rule1D[i]);
    CubatureLineSorted<Scalar> diffRule1(rule1D[i],numPoints,isNormalized);
 
    if (numPoints!=1) { // Compute differential rule
      numPoints = growthRule1D(index[i]-1,growth1D[i],rule1D[i]);
      CubatureLineSorted<Scalar> rule1(rule1D[i],numPoints,isNormalized);
      diffRule1.update(-1.0,rule1,1.0);
    }
    // Build Tensor Rule
    outRule = kron_prod<Scalar>(outRule,diffRule1); 
  }
}

template<class Scalar, class UserVector>
void AdaptiveSparseGrid<Scalar,UserVector>::build_diffRule(
	       CubatureTensorSorted<Scalar> & outRule, 
	       std::vector<int> index, int dimension,
	       std::vector<EIntrepidBurkardt> rule1D, 
	       std::vector<EIntrepidGrowth> growth1D,
	       bool isNormalized) {
	       
  int numPoints = 0;
  for (int i=0; i<dimension; i++) {
    // Compute 1D rules
    numPoints = growthRule1D(index[i],growth1D[i],rule1D[i]);
    CubatureLineSorted<Scalar> diffRule1(rule1D[i],numPoints,isNormalized);
 
    if (numPoints!=1) { // Differential Rule
      numPoints = growthRule1D(index[i]-1,growth1D[i],rule1D[i]);
      CubatureLineSorted<Scalar> rule1(rule1D[i],numPoints,isNormalized);
      diffRule1.update(-1.0,rule1,1.0);
    }
    // Build Tensor Rule
    outRule = kron_prod<Scalar>(outRule,diffRule1); 
  }
}
 
// Update Index Set - no knowledge of active or old indices
template<class Scalar, class UserVector>
Scalar AdaptiveSparseGrid<Scalar,UserVector>::refine_grid(
   	         typename std::multimap<Scalar,std::vector<int> >  & indexSet, 
		 UserVector & integralValue,
		 AdaptiveSparseGridInterface<Scalar,UserVector> & problem_data) {
  
  int dimension = problem_data.getDimension();
  std::vector<EIntrepidBurkardt> rule1D; problem_data.getRule(rule1D);
  std::vector<EIntrepidGrowth> growth1D; problem_data.getGrowth(growth1D);
  
  // Copy Multimap into a Set for ease of use
  typename std::multimap<Scalar,std::vector<int> >::iterator it;
  std::set<std::vector<int> > oldSet;
  std::set<std::vector<int> >::iterator it1(oldSet.begin());
  for (it=indexSet.begin(); it!=indexSet.end(); it++) {
    oldSet.insert(it1,it->second);
    it1++;
  }
  indexSet.clear();
  
  // Find Possible Active Points
  int flag = 1;
  std::vector<int> index(dimension,0);
  typename std::multimap<Scalar,std::vector<int> > activeSet;
  for (it1=oldSet.begin(); it1!=oldSet.end(); it1++) {
    index = *it1;
    for (int i=0; i<dimension; i++) {
      index[i]++;
      flag = (int)(!oldSet.count(index));
      index[i]--;
      if (flag) {
	activeSet.insert(std::pair<Scalar,std::vector<int> >(1.0,index));
	oldSet.erase(it1);
	break;
      }
    }
  }

  // Compute local and global error indicators for active set
  typename std::multimap<Scalar,std::vector<int> >::iterator it2;
  Scalar eta = 0.0;
  Scalar G   = 0.0;
  Teuchos::RCP<UserVector> s = integralValue.Create();
  for (it2=activeSet.begin(); it2!=activeSet.end(); it2++) {
    // Build Differential Quarature Rule
    index = it2->second;
    CubatureTensorSorted<Scalar> diffRule(0,dimension);
    build_diffRule(diffRule,index,problem_data);
    
    // Apply Rule to function
    problem_data.eval_cubature(*s,diffRule);
    
    // Update local error indicator and index set
    G  = problem_data.error_indicator(*s);
    activeSet.erase(it2);
    activeSet.insert(it2,std::pair<Scalar,std::vector<int> >(G,index));
    eta += G;
  }

  // Refine Sparse Grid
  eta = refine_grid(activeSet,oldSet,integralValue,eta,
		    dimension,rule1D,growth1D);

  // Insert New Active and Old Index sets into indexSet
  indexSet.insert(activeSet.begin(),activeSet.end());
  for (it1=oldSet.begin(); it1!=oldSet.end(); it1++) {
    index = *it1;
    indexSet.insert(std::pair<Scalar,std::vector<int> >(-1.0,index));
  }
    
  return eta;
}

// Update index set and output integral
template<class Scalar, class UserVector> 
Scalar AdaptiveSparseGrid<Scalar,UserVector>::refine_grid(
          typename std::multimap<Scalar,std::vector<int> > & activeIndex, 
      	  std::set<std::vector<int> > & oldIndex, 
	  UserVector & integralValue,
	  Scalar globalErrorIndicator,
	  AdaptiveSparseGridInterface<Scalar,UserVector> & problem_data) {
 
 TEST_FOR_EXCEPTION((activeIndex.empty()),std::out_of_range,
              ">>> ERROR (AdaptiveSparseGrid): Active Index set is empty.");  

  int dimension = problem_data.getDimension();
  std::vector<EIntrepidBurkardt> rule1D; problem_data.getRule(rule1D);
  std::vector<EIntrepidGrowth> growth1D; problem_data.getGrowth(growth1D);

  // Initialize Flags
  bool maxLevelFlag     = true;
  bool isAdmissibleFlag = true;

  // Initialize Cubature Rule
  Teuchos::RCP<UserVector> s = integralValue.Create();
  // Initialize iterator at end of inOldIndex
  std::set<std::vector<int> >::iterator it1(oldIndex.end()); 
  
  if (oldIndex.end()!=oldIndex.begin()) 
    it1--;  
 
  // Initialize iterator at end of inActiveIndex
  typename std::multimap<Scalar,std::vector<int> >::iterator it;

  // Obtain Global Error Indicator as sum of key values of inActiveIndex
  Scalar eta = globalErrorIndicator;

  // Select Index to refine
  it = activeIndex.end(); it--;        // Decrement to position of final value
  Scalar G               = it->first;  // Largest Error Indicator is at End 
  eta                   -= G;          // Update global error indicator
  std::vector<int> index = it->second; // Get Corresponding index
  activeIndex.erase(it);               // Erase Index from active index set
  oldIndex.insert(it1,index); it1++;   // Insert Index into old index set

  // Refinement process
  for (int k=0; k<dimension; k++) {
    index[k]++; // index + ek
    // Check Max Level
    maxLevelFlag = problem_data.max_level(index);
    if (maxLevelFlag) {
      // Check Admissibility
      isAdmissibleFlag = isAdmissible(index,k,oldIndex,problem_data);
      if (isAdmissibleFlag) { // If admissible
	// Build Differential Quarature Rule
	CubatureTensorSorted<Scalar> diffRule(0,dimension);
	build_diffRule(diffRule,index,problem_data);	

	// Apply Rule to function
	problem_data.eval_cubature(*s,diffRule);
	
	// Update integral value
	integralValue.Update(*s);

	// Update local error indicator and index set
	G  = problem_data.error_indicator(*s); 
	if (activeIndex.end()!=activeIndex.begin()) 
	  activeIndex.insert(activeIndex.end()--,
			   std::pair<Scalar,std::vector<int> >(G,index));
	else
	  activeIndex.insert(std::pair<Scalar,std::vector<int> >(G,index));

	// Update global error indicators
	eta += G;
      }
    }
    else { // Max Level Exceeded 
      //std::cout << "Maximum Level Exceeded" << std::endl;
    }
    index[k]--;
  }
  return eta;
}

// Update index set and output integral/sparse grid
template<class Scalar, class UserVector> 
Scalar AdaptiveSparseGrid<Scalar,UserVector>::refine_grid(
	typename std::multimap<Scalar,std::vector<int> > & activeIndex, 
	std::set<std::vector<int> > & oldIndex, 
	UserVector & integralValue,
	CubatureTensorSorted<Scalar> & cubRule,
	Scalar globalErrorIndicator,
	AdaptiveSparseGridInterface<Scalar,UserVector> & problem_data) {

  TEST_FOR_EXCEPTION((activeIndex.empty()),std::out_of_range,
              ">>> ERROR (AdaptiveSparseGrid): Active Index set is empty.");  

  int dimension = problem_data.getDimension();
  std::vector<EIntrepidBurkardt> rule1D; problem_data.getRule(rule1D);
  std::vector<EIntrepidGrowth> growth1D; problem_data.getGrowth(growth1D);

  // Initialize Flags
  bool maxLevelFlag     = true;
  bool isAdmissibleFlag = true;

  // Initialize Cubature Rule
  Teuchos::RCP<UserVector> s = integralValue.Create();

  // Initialize iterator at end of inOldIndex
  std::set<std::vector<int> >::iterator it1(oldIndex.end());  
  if (oldIndex.end()!=oldIndex.begin()) 
    it1--;
  // Initialize iterator at end of inActiveIndex
  typename std::multimap<Scalar,std::vector<int> >::iterator it;

  // Obtain Global Error Indicator as sum of key values of inActiveIndex
  Scalar eta = globalErrorIndicator;

  // Select Index to refine
  it = activeIndex.end(); it--;        // Decrement to position of final value 
  Scalar G               = it->first;  // Largest Error Indicator is at End
  eta                   -= G;          // Update global error indicator
  std::vector<int> index = it->second; // Get Corresponding index
  activeIndex.erase(it);               // Erase Index from active index set
  oldIndex.insert(it1,index); it1++;   // Insert Index into old index set
  
  // Refinement process
  for (int k=0; k<dimension; k++) {
    index[k]++; // index + ek
    // Check Max Level
    maxLevelFlag = problem_data.max_level(index);
    if (maxLevelFlag) {
      // Check Admissibility
      isAdmissibleFlag = isAdmissible(index,k,oldIndex,problem_data);
      if (isAdmissibleFlag) { // If admissible
	// Build Differential Quarature Rule
	CubatureTensorSorted<Scalar> diffRule(0,dimension);
	build_diffRule(diffRule,index,problem_data);
	
	// Apply Rule to function
	problem_data.eval_cubature(*s,diffRule);
	
	// Update integral value
	integralValue.Update(*s);
	
	// Update local error indicator and index set
	G  = problem_data.error_indicator(*s); 	
	if (activeIndex.end()!=activeIndex.begin()) 
	  activeIndex.insert(activeIndex.end()--,
			   std::pair<Scalar,std::vector<int> >(G,index));
	else
	  activeIndex.insert(std::pair<Scalar,std::vector<int> >(G,index));
		
	// Update global error indicators
	eta += G;

	// Update adapted quadrature rule nodes and weights
	cubRule.update(1.0,diffRule,1.0);
      }
    }
    else { // Max Level Exceeded 
      //std::cout << "Maximum Level Exceeded" << std::endl;
    }
    index[k]--;
  }
  return eta;
}

template<class Scalar, class UserVector>
void AdaptiveSparseGrid<Scalar,UserVector>::buildSparseGrid(
	      CubatureTensorSorted<Scalar> & output,
	      int dimension, int maxlevel,
	      std::vector<EIntrepidBurkardt> rule1D, 
	      std::vector<EIntrepidGrowth> growth1D,
	      bool isNormalized) {

  if (dimension == 2) {
    std::vector<int> index(dimension,0);
    for (int i=0; i<maxlevel; i++) {
      for (int j=0; j<maxlevel; j++) {
	if(i+j+dimension <= maxlevel+dimension-1) {
	  index[0] = i+1; index[1] = j+1; 
	  CubatureTensorSorted<Scalar> diffRule(0,dimension);
	  build_diffRule(diffRule,index,dimension,rule1D,growth1D,isNormalized);
	  output.update(1.0,diffRule,1.0);
	}
      }
    } 
  }
  else if (dimension == 3) {    
    std::vector<int> index(dimension,0);
    for (int i=0; i<maxlevel; i++) {
      for (int j=0; j<maxlevel; j++) {
	for (int k=0; k<maxlevel; k++) {
	  if(i+j+k+dimension <= maxlevel+dimension-1) {
	    index[0] = i+1; index[1] = j+1; index[2] = k+1;
	    CubatureTensorSorted<Scalar> diffRule(0,dimension);
	    build_diffRule(diffRule,index,dimension,rule1D,
			  growth1D,isNormalized);
	    output.update(1.0,diffRule,1.0);
	  }
	}
      }
    } 
  }
  else if (dimension == 4) {
    std::vector<int> index(dimension,0);
    for (int i=0; i<maxlevel; i++) {
      for (int j=0; j<maxlevel; j++) {
	for (int k=0; k<maxlevel; k++) {
	  for (int l=0; l<maxlevel; l++) {
	    if(i+j+k+l+dimension <= maxlevel+dimension-1) {
	      index[0] = i+1; index[1] = j+1; index[2] = k+1; index[3] = l+1;
	      CubatureTensorSorted<Scalar> diffRule(0,dimension);
	      build_diffRule(diffRule,index,dimension,rule1D,
			    growth1D,isNormalized);
	      output.update(1.0,diffRule,1.0);
	    }
	  }
	}
      } 
    }
  }
  else if (dimension == 5) {
    std::vector<int> index(dimension,0);
    for (int i=0; i<maxlevel; i++) {
      for (int j=0; j<maxlevel; j++) {
	for (int k=0; k<maxlevel; k++) {
	  for (int l=0; l<maxlevel; l++) {
	    for (int m=0; m<maxlevel; m++) {
	      if(i+j+k+l+m+dimension <= maxlevel+dimension-1) {
		index[0] = i+1; index[1] = j+1; index[2] = k+1; 
		index[3] = l+1; index[4] = m+1;
		CubatureTensorSorted<Scalar> diffRule(0,dimension);
		build_diffRule(diffRule,index,dimension,rule1D,
			      growth1D,isNormalized);
		output.update(1.0,diffRule,1.0);
	      }
	    }
	  }
	} 
      }
    }
  }
  else 
    std::cout << "Dimension Must Be Less Than 5\n";
}

} // end Intrepid namespace

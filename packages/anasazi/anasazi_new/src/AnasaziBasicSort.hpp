// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef ANASAZI_BASIC_SORT_HPP
#define ANASAZI_BASIC_SORT_HPP

/*!    \class Anasazi::BasicSort
       \brief An implementation of the Anasazi::SortManager that performs a collection
       of common sorting techniques.

       \author Heidi Thornquist
*/

#include "AnasaziSortManager.hpp"
#include "Teuchos_LAPACK.hpp"

namespace Anasazi {

  template<class TYPE>
  class BasicSort : public SortManager<TYPE> {
    
  public:
    
    //! Constructor
    /**
       @param which [in] The eigenvalues of interest for this eigenproblem.
       <ul>
       <li> "LM" - Largest Magnitude [ default ]
       <li> "SM" - Smallest Magnitude
       <li> "LR" - Largest Real 
       <li> "SR" - Smallest Real 
       <li> "LI" - Largest Imaginary 
       <li> "SI" - Smallest Imaginary 
       </ul>
    */
    BasicSort( const string which = "LM" ) { _which = which; };

    //! Destructor
    virtual ~BasicSort() {};


    //! Sort the vector of eigenvalues with respect to the chosen sorting type, optionally returning the permutation vector.
    /**
       @param n [in] Size of the array

       @param evals [in/out] Array of length n containing the eigenvalues to be sorted

       @param perm [out] Array of length n to store the permutation (optional)

       @return Returns the status of the sorting routine [ Undefined by default ] 
    */
    ReturnType sort(int n, TYPE *evals, int *perm = 0) const;
    
    //! Sort the vectors of eigenpairs with respect to the chosen sorting type, optionally returning the permutation vector.
    /**
       @param n [in] Size of the array

       @param r_evals [in/out] Array of length n containing the real part of the eigenvalues to be sorted 

       @param i_evals [in/out] Array of length n containing the imaginary part of the eigenvalues to be sorted 

       @param perm [out] Array of length n to store the permutation (optional)

       @return Returns the status of the sorting routine [ Undefined by default ] 
    */
    ReturnType sort(int n, TYPE *r_evals, TYPE *i_evals, int *perm = 0) const;
    
  protected: 
    
    string _which;

  };

  template<class TYPE>
  ReturnType BasicSort<TYPE>::sort(int n, TYPE *evals, int *perm) const 
  {
    int i, j, tempord;
    TYPE temp, temp2;
    Teuchos::LAPACK<int,TYPE> lapack;
    //
    // Reset the permutation if it is required.
    //		
    if (perm) {
      for (i=0; i < n; i++) {
	perm[i] = i;
      }
    }
    //
    // These methods use an insertion sort method to circument recursive calls.
    //---------------------------------------------------------------
    // Sort eigenvalues in increasing order of magnitude
    //---------------------------------------------------------------
    if (!_which.compare("SM")) {
      for (j=1; j < n; ++j) {
	temp = evals[j]; 
	if (perm)
	  tempord = perm[j];
	temp2 = evals[j]*evals[j];
	for (i=j-1; i>=0 && (evals[i]*evals[i])>temp2; --i) {
	  evals[i+1]=evals[i];
	  if (perm)
	    perm[i+1]=perm[i];
	}
	evals[i+1] = temp; 
	if (perm) 
	  perm[i+1] = tempord;	
      }
      return Ok;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in increasing order of real part
    //---------------------------------------------------------------
    if (!_which.compare("SR")) {
      for (j=1; j < n; ++j) {
	temp = evals[j]; 
	if (perm)
	  tempord = perm[j];
	for (i=j-1; i>=0 && evals[i]>temp; --i) {
	  evals[i+1]=evals[i];
	  if (perm)
	    perm[i+1]=perm[i];
	}
	evals[i+1] = temp; 
	if (perm)
	  perm[i+1] = tempord;	
      }
      return Ok;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in increasing order of imaginary part
    // NOTE:  There is no implementation for this since this sorting
    // method assumes only real eigenvalues.
    //---------------------------------------------------------------
    if (!_which.compare("SI")) {
      return Undefined;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in decreasing order of magnitude
    //---------------------------------------------------------------
    if (!_which.compare("LM")) {
      for (j=1; j < n; ++j) {
	temp = evals[j]; 
	if (perm)
	  tempord = perm[j];
	temp2 = evals[j]*evals[j];
	for (i=j-1; i>=0 && (evals[i]*evals[i])<temp2; --i) {
	  evals[i+1]=evals[i];
	  if (perm)
	    perm[i+1]=perm[i];
	}
	evals[i+1] = temp; 
	if (perm)
	  perm[i+1] = tempord;	
      }
      return Ok;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in decreasing order of real part
    //---------------------------------------------------------------
    if (!_which.compare("LR")) {
      for (j=1; j < n; ++j) {
	temp = evals[j]; 
	if (perm)
	  tempord = perm[j];
	for (i=j-1; i>=0 && evals[i]<temp; --i) {
	  evals[i+1]=evals[i];
	  if (perm)
	    perm[i+1]=perm[i];
	}
	evals[i+1] = temp; 
	if (perm)
	  perm[i+1] = tempord;	
      }
      return Ok;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in decreasing order of imaginary part
    // NOTE:  There is no implementation for this since this sorting
    // method assumes only real eigenvalues.
    //---------------------------------------------------------------
    if (!_which.compare("LI")) {
      return Undefined;
    }
    
    // The character string held by this class is not valid.  
    return Undefined;    
  }
  

  template<class TYPE>
  ReturnType BasicSort<TYPE>::sort(int n, TYPE *r_evals, TYPE *i_evals, int *perm) const {
    int i, j, tempord;
    TYPE temp, tempr, tempi;
    Teuchos::LAPACK<int,TYPE> lapack;
    //
    // Reset the index
    //		
    if (perm) {
      for (i=0; i < n; i++) {
	perm[i] = i;
      }
    }
    //
    // These methods use an insertion sort method to circument recursive calls.
    //---------------------------------------------------------------
    // Sort eigenvalues in increasing order of magnitude
    //---------------------------------------------------------------
    if (!_which.compare("SM")) {
      for (j=1; j < n; ++j) {
	tempr = r_evals[j]; tempi = i_evals[j]; 
	if (perm)
	  tempord = perm[j];
	temp=lapack.LAPY2(r_evals[j],i_evals[j]);
	for (i=j-1; i>=0 && lapack.LAPY2(r_evals[i],i_evals[i])>temp; --i) {
	  r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
	  if (perm)
	    perm[i+1]=perm[i];
	}
	r_evals[i+1] = tempr; i_evals[i+1] = tempi; 
	if (perm)
	  perm[i+1] = tempord;	
      }	
      return Ok;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in increasing order of real part
    //---------------------------------------------------------------
    if (!_which.compare("SR")) {
      for (j=1; j < n; ++j) {
	tempr = r_evals[j]; tempi = i_evals[j]; 
	if (perm)
	  tempord = perm[j];
	for (i=j-1; i>=0 && r_evals[i]>tempr; --i) {
	  r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
	  if (perm)
	    perm[i+1]=perm[i];
	}
	r_evals[i+1] = tempr; i_evals[i+1] = tempi; 
	if (perm)
	  perm[i+1] = tempord;	
      }	
      return Ok;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in increasing order of imaginary part
    //---------------------------------------------------------------
    if (!_which.compare("SI")) {
      for (j=1; j < n; ++j) {
	tempr = r_evals[j]; tempi = i_evals[j]; 
	if (perm)
	  tempord = perm[j];
	for (i=j-1; i>=0 && i_evals[i]>tempi; --i) {
	  r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
	  if (perm)
	    perm[i+1]=perm[i];
	}
	r_evals[i+1] = tempr; i_evals[i+1] = tempi; 
	if (perm)
	  perm[i+1] = tempord;	
      }
      return Ok;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in decreasing order of magnitude
    //---------------------------------------------------------------
    if (!_which.compare("LM")) {
      for (j=1; j < n; ++j) {
	tempr = r_evals[j]; tempi = i_evals[j]; 
	if (perm)
	  tempord = perm[j];
	temp=lapack.LAPY2(r_evals[j],i_evals[j]);
	for (i=j-1; i>=0 && lapack.LAPY2(r_evals[i],i_evals[i])<temp; --i) {
	  r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
	  if (perm)
	    perm[i+1]=perm[i];
	}
	r_evals[i+1] = tempr; i_evals[i+1] = tempi; 
	if (perm)
	  perm[i+1] = tempord;	
      }	
      return Ok;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in decreasing order of real part
    //---------------------------------------------------------------
    if (!_which.compare("LR")) {
      for (j=1; j < n; ++j) {
	tempr = r_evals[j]; tempi = i_evals[j]; 
	if (perm)
	  tempord = perm[j];
	for (i=j-1; i>=0 && r_evals[i]<tempr; --i) {
	  r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
	  if (perm)
	    perm[i+1]=perm[i];
	}
	r_evals[i+1] = tempr; i_evals[i+1] = tempi; 
	if (perm)
	  perm[i+1] = tempord;	
      }	
      return Ok;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in decreasing order of imaginary part
    //---------------------------------------------------------------
    if (!_which.compare("LI")) {
      for (j=1; j < n; ++j) {
	tempr = r_evals[j]; tempi = i_evals[j]; 
	if (perm)
	  tempord = perm[j];
	for (i=j-1; i>=0 && i_evals[i]<tempi; --i) {
	  r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
	  if (perm)
	    perm[i+1]=perm[i];
	}
	r_evals[i+1] = tempr; i_evals[i+1] = tempi; 
	if (perm)
	  perm[i+1] = tempord;	
      }
      return Ok;
    }

    // The character string held by this class is not valid.  
    return Undefined;      
  }
  
} // namespace Anasazi

#endif // ANASAZI_BASIC_SORT_HPP


//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef ANASAZI_LOCA_SORT_HPP
#define ANASAZI_LOCA_SORT_HPP

/*!    \class Anasazi::LOCA::Sort
       \brief An implementation of the Anasazi::SortManager that performs a collection
       of common sorting techniques.
*/

#include "AnasaziSortManager.hpp"
#include "Teuchos_LAPACK.hpp"

namespace Anasazi {
  
  namespace LOCA {
    
    template<class ScalarType, class MV, class OP>
    class Sort : public SortManager<ScalarType,MV,OP> {
      
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
	 <li> "CA" - Cayley Sorting
	 </ul>
      */
      Sort( const string which = "LM", 
	    ScalarType cayleyPole = Teuchos::ScalarTraits<ScalarType>::zero(),
	    ScalarType cayleyZero = Teuchos::ScalarTraits<ScalarType>::zero()
	    ) 
      { _which = which; _sigma = cayleyPole; _mu = cayleyZero; };
      
      //! Destructor
      virtual ~Sort() {};
      
      
      //! Sort the vector of eigenvalues with respect to the chosen sorting type, optionally returning the permutation vector.
      /**
	 @param solver [in] Eigensolver that is calling the sorting routine
	 
	 @param n [in] Size of the array
	 
	 @param evals [in/out] Array of length n containing the eigenvalues to be sorted
	 
	 @param perm [out] Vector of length n to store the permutation (optional)
	 
	 @return Returns the status of the sorting routine [ Undefined by default ] 
      */
      ReturnType sort(Eigensolver<ScalarType,MV,OP>* solver, int n, ScalarType *evals, std::vector<int> *perm = 0) const;
      
      //! Sort the vectors of eigenpairs with respect to the chosen sorting type, optionally returning the permutation vector.
      /**
	 @param solver [in] Eigensolver that is calling the sorting routine
	 
	 @param n [in] Size of the array
	 
	 @param r_evals [in/out] Array of length n containing the real part of the eigenvalues to be sorted 
	 
	 @param i_evals [in/out] Array of length n containing the imaginary part of the eigenvalues to be sorted 
	 
	 @param perm [out] Vector of length n to store the permutation (optional)
	 
	 @return Returns the status of the sorting routine [ Undefined by default ] 
      */
      ReturnType sort(Eigensolver<ScalarType,MV,OP>* solver, int n, ScalarType *r_evals, ScalarType *i_evals, std::vector<int> *perm = 0) const;
      
    protected: 
      
      string _which;
      ScalarType _sigma, _mu;
      
    private:
      
      // Helper method for sorting cayley transform
      ScalarType realLambda(const ScalarType er, const ScalarType ei) const;      
      
    };
    
    template<class ScalarType, class MV, class OP>
    ReturnType Sort<ScalarType,MV,OP>::sort(Eigensolver<ScalarType,MV,OP>* solver, int n, ScalarType *evals, std::vector<int> *perm) const 
    {
      int i, j, tempord;
      ScalarType temp, temp2;
      Teuchos::LAPACK<int,ScalarType> lapack;
      //
      // Reset the permutation if it is required.
      //		
      if (perm) {
	for (i=0; i < n; i++) {
	  (*perm)[i] = i;
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
	    tempord = (*perm)[j];
	  temp2 = evals[j]*evals[j];
	  for (i=j-1; i>=0 && (evals[i]*evals[i])>temp2; --i) {
	    evals[i+1]=evals[i];
	    if (perm)
	      (*perm)[i+1]=(*perm)[i];
	  }
	  evals[i+1] = temp; 
	  if (perm) 
	    (*perm)[i+1] = tempord;	
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
	    tempord = (*perm)[j];
	  for (i=j-1; i>=0 && evals[i]>temp; --i) {
	    evals[i+1]=evals[i];
	    if (perm)
	      (*perm)[i+1]=(*perm)[i];
	  }
	  evals[i+1] = temp; 
	  if (perm)
	    (*perm)[i+1] = tempord;	
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
	    tempord = (*perm)[j];
	  temp2 = evals[j]*evals[j];
	  for (i=j-1; i>=0 && (evals[i]*evals[i])<temp2; --i) {
	    evals[i+1]=evals[i];
	    if (perm)
	      (*perm)[i+1]=(*perm)[i];
	  }
	  evals[i+1] = temp; 
	  if (perm)
	    (*perm)[i+1] = tempord;	
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
	    tempord = (*perm)[j];
	  for (i=j-1; i>=0 && evals[i]<temp; --i) {
	    evals[i+1]=evals[i];
	    if (perm)
	      (*perm)[i+1]=(*perm)[i];
	  }
	  evals[i+1] = temp; 
	  if (perm)
	    (*perm)[i+1] = tempord;	
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
    //---------------------------------------------------------------
    // Sort eigenvalues however Andy wants them ( which is always correct )
    //---------------------------------------------------------------
    if (!_which.compare("CA")) {
      // Sigma and mu are available here      
      // use _evalr and _evali and the ordering goes in _order
      // remember to actually sort the eigenvalues yourself
      //
      // LM code to start with....
      ScalarType templambda;
	for (j=1; j < n; ++j) {
	  temp = evals[j]; 
	  tempord = (*perm)[j];
	  templambda=realLambda(evals[j],0);
	  for (i=j-1; i>=0 && realLambda(evals[i],0)<templambda; --i) {
	    evals[i+1]=evals[i]; 
	    (*perm)[i+1]=(*perm)[i];
	  }
	  evals[i+1] = temp; (*perm)[i+1] = tempord;	
	}	
	return Ok;
    }
    
    // The character string held by this class is not valid.  
    return Undefined;    
    }
    
    
    template<class ScalarType, class MV, class OP>
    ReturnType Sort<ScalarType,MV,OP>::sort(Eigensolver<ScalarType,MV,OP>* solver, int n, ScalarType *r_evals, ScalarType *i_evals, std::vector<int> *perm) const {
      int i, j, tempord;
      ScalarType temp, tempr, tempi;
      Teuchos::LAPACK<int,ScalarType> lapack;
      //
      // Reset the index
      //		
      if (perm) {
	for (i=0; i < n; i++) {
	  (*perm)[i] = i;
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
	    tempord = (*perm)[j];
	  temp=lapack.LAPY2(r_evals[j],i_evals[j]);
	  for (i=j-1; i>=0 && lapack.LAPY2(r_evals[i],i_evals[i])>temp; --i) {
	    r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
	    if (perm)
	      (*perm)[i+1]=(*perm)[i];
	  }
	  r_evals[i+1] = tempr; i_evals[i+1] = tempi; 
	  if (perm)
	    (*perm)[i+1] = tempord;	
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
	    tempord = (*perm)[j];
	  for (i=j-1; i>=0 && r_evals[i]>tempr; --i) {
	    r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
	    if (perm)
	      (*perm)[i+1]=(*perm)[i];
	  }
	  r_evals[i+1] = tempr; i_evals[i+1] = tempi; 
	  if (perm)
	    (*perm)[i+1] = tempord;	
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
	    tempord = (*perm)[j];
	  for (i=j-1; i>=0 && i_evals[i]>tempi; --i) {
	    r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
	    if (perm)
	      (*perm)[i+1]=(*perm)[i];
	  }
	  r_evals[i+1] = tempr; i_evals[i+1] = tempi; 
	  if (perm)
	    (*perm)[i+1] = tempord;	
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
	    tempord = (*perm)[j];
	  temp=lapack.LAPY2(r_evals[j],i_evals[j]);
	  for (i=j-1; i>=0 && lapack.LAPY2(r_evals[i],i_evals[i])<temp; --i) {
	    r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
	    if (perm)
	      (*perm)[i+1]=(*perm)[i];
	  }
	  r_evals[i+1] = tempr; i_evals[i+1] = tempi; 
	  if (perm)
	    (*perm)[i+1] = tempord;	
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
	    tempord = (*perm)[j];
	  for (i=j-1; i>=0 && r_evals[i]<tempr; --i) {
	    r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
	    if (perm)
	      (*perm)[i+1]=(*perm)[i];
	  }
	  r_evals[i+1] = tempr; i_evals[i+1] = tempi; 
	  if (perm)
	    (*perm)[i+1] = tempord;	
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
	    tempord = (*perm)[j];
	  for (i=j-1; i>=0 && i_evals[i]<tempi; --i) {
	    r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
	    if (perm)
	      (*perm)[i+1]=(*perm)[i];
	  }
	  r_evals[i+1] = tempr; i_evals[i+1] = tempi; 
	  if (perm)
	    (*perm)[i+1] = tempord;	
	}
	return Ok;
      }
    //---------------------------------------------------------------
    // Sort eigenvalues however Andy wants them ( which is always correct )
    //---------------------------------------------------------------
    if (!_which.compare("CA")) {
      // Sigma and mu are available here      
      // use _evalr and _evali and the ordering goes in _order
      // remember to actually sort the eigenvalues yourself
      //
      // LM code to start with....
      for (j=1; j < n; ++j) {
	tempr = r_evals[j]; tempi = i_evals[j]; 
	tempord = (*perm)[j];
	temp=realLambda(r_evals[j],i_evals[j]);
	for (i=j-1; i>=0 && realLambda(r_evals[i],i_evals[i])<temp; --i) {
	  r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
	  (*perm)[i+1]=(*perm)[i];
	}
	r_evals[i+1] = tempr; i_evals[i+1] = tempi; (*perm)[i+1] = tempord;	
      }
      return Ok;
    }      
    // The character string held by this class is not valid.  
    return Undefined;      
    }
    
    template<class ScalarType, class MV, class OP>
    ScalarType Sort<ScalarType, MV, OP>::realLambda(const ScalarType er, const ScalarType ei) const {
      // Utility that returns the real part of the un-Cayley-transformed eigenvalue for sorting
      // Reject if it is to the right of sigma --- these are junk
      // This is temporary, to be replaced by sorting class
      ScalarType reLambda =  (_sigma*(er*er+ei*ei) - (_sigma+_mu)*er + _mu)/ ( (er-1.0)*(er-1.0) + ei*ei);
      if (reLambda > _sigma) return -1.0e6;
      else                   return reLambda;
    }
    
  }// namespace LOCA
  
} // namespace Anasazi

#endif // ANASAZI_LOCA_SORT_HPP


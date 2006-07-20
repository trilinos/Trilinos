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

/*! \file AnasaziBasicSort.hpp
  \brief Basic implementation of the Anasazi::SortManager class
*/

#ifndef ANASAZI_BASIC_SORT_HPP
#define ANASAZI_BASIC_SORT_HPP

/*!    \class Anasazi::BasicSort
       \brief An implementation of the Anasazi::SortManager that performs a collection
       of common sorting techniques.

       \author Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziSortManager.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Anasazi {

  template<class ScalarType, class MV, class OP>
  class BasicSort : public SortManager<ScalarType,MV,OP> {
    
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

    //! Set sort type
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
    void SetSortType( const string which ) { _which = which; };
    
    //! Sort the vector of eigenvalues with respect to the chosen sorting type, optionally returning the permutation vector.
    /**
       @param solver [in] Eigensolver that is calling the sorting routine

       @param n [in] Size of the array

       @param evals [in/out] Array of length n containing the eigenvalues to be sorted

       @param perm [out] Vector of length n to store the permutation (optional)
    */
    void sort(Eigensolver<ScalarType,MV,OP>* solver, int n, ScalarType *evals, std::vector<int> *perm = 0) const;
    
    //! Sort the vectors of eigenpairs with respect to the chosen sorting type, optionally returning the permutation vector.
    /**
       @param solver [in] Eigensolver that is calling the sorting routine

       @param n [in] Size of the array

       @param r_evals [in/out] Array of length n containing the real part of the eigenvalues to be sorted 

       @param i_evals [in/out] Array of length n containing the imaginary part of the eigenvalues to be sorted 

       @param perm [out] Vector of length n to store the permutation (optional)
    */
    void sort(Eigensolver<ScalarType,MV,OP>* solver, int n, ScalarType *r_evals, ScalarType *i_evals, std::vector<int> *perm = 0) const;
    
  protected: 
    
    //! Sorting type
    /*! \note Sorting choices:
       <ul>
       <li> "LM" - Largest Magnitude [ default ]
       <li> "SM" - Smallest Magnitude
       <li> "LR" - Largest Real 
       <li> "SR" - Smallest Real 
       <li> "LI" - Largest Imaginary 
       <li> "SI" - Smallest Imaginary 
       </ul>
    */
    string _which;

  };

  template<class ScalarType, class MV, class OP>
  void BasicSort<ScalarType,MV,OP>::sort(Eigensolver<ScalarType,MV,OP>* solver, int n, ScalarType *evals, std::vector<int> *perm) const 
  {
    int i=0, j=0;

    // Temp integer for swapping the index of the permutation, used in all sorting types.
    int tempord=0;

    // Temp variable for the magnitude of the ScalarType used in sorting "LM" and "SM".
    typename Teuchos::ScalarTraits<ScalarType>::magnitudeType temp2;

    // Temp variable for swapping the eigenvalue used in all sorting types.
    ScalarType temp;

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
        if (perm) {
          tempord = (*perm)[j];
        }
        temp2 = Teuchos::ScalarTraits<ScalarType>::magnitude(evals[j]);
        for (i=j-1; i>=0 && Teuchos::ScalarTraits<ScalarType>::magnitude(evals[i])>temp2; --i) {
          evals[i+1]=evals[i];
          if (perm)
            (*perm)[i+1]=(*perm)[i];
        }
        evals[i+1] = temp; 
        if (perm) 
          (*perm)[i+1] = tempord;
      }
      return;
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
      return;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in increasing order of imaginary part
    // NOTE:  There is no implementation for this since this sorting
    // method assumes only real eigenvalues.
    //---------------------------------------------------------------
    TEST_FOR_EXCEPTION(!_which.compare("SI"), SortManagerError, 
                       "Anasazi::BasicSort::sort() assumes real eigenvalues");
    //---------------------------------------------------------------
    // Sort eigenvalues in decreasing order of magnitude
    //---------------------------------------------------------------
    if (!_which.compare("LM")) {
      for (j=1; j < n; ++j) {
        temp = evals[j]; 
        if (perm)
          tempord = (*perm)[j];
        temp2 = Teuchos::ScalarTraits<ScalarType>::magnitude(evals[j]);
        for (i=j-1; i>=0 && Teuchos::ScalarTraits<ScalarType>::magnitude(evals[i])<temp2; --i) {
          evals[i+1]=evals[i];
          if (perm)
            (*perm)[i+1]=(*perm)[i];
        }
        evals[i+1] = temp; 
        if (perm)
          (*perm)[i+1] = tempord;
        }
      return;
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
      return;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in decreasing order of imaginary part
    // NOTE:  There is no implementation for this since this templating
    // assumes only real eigenvalues.
    //---------------------------------------------------------------
    TEST_FOR_EXCEPTION(!_which.compare("LI"), SortManagerError, 
                       "Anasazi::BasicSort::sort() assumes real eigenvalues");
    
    // The character string held by this class is not valid.  
    TEST_FOR_EXCEPTION(true, SortManagerError, 
                       "Anasazi::BasicSort::sort(): sorting order is not valid");
  }


  template<class ScalarType, class MV, class OP>
  void BasicSort<ScalarType,MV,OP>::sort(Eigensolver<ScalarType,MV,OP>* solver, int n, ScalarType *r_evals, ScalarType *i_evals, std::vector<int> *perm) const {
    int i=0, j=0, tempord=0;
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
      return;
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
      return;
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
      return;
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
      return;
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
      return;
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
      return;
    }

    TEST_FOR_EXCEPTION(true, SortManagerError, 
                       "Anasazi::BasicSort::sort(): sorting order is not valid");
  }
  
#if ( defined(HAVE_COMPLEX) || defined(HAVE_COMPLEX_H) ) && defined(HAVE_TEUCHOS_COMPLEX)


  // Define a macro for which complex class we have available
  
#if defined(HAVE_COMPLEX)
#define ANSZI_CPLX_CLASS std::complex
#elif  defined(HAVE_COMPLEX_H)
#define ANSZI_CPLX_CLASS ::complex
#endif

  // ----------------------------------------------------------------------------
  //  Template specialization for complex<> scalar types
  // ----------------------------------------------------------------------------
  template<class ScalarType, class MV, class OP>
  class BasicSort<ANSZI_CPLX_CLASS<ScalarType>, MV, OP> : public SortManager<ANSZI_CPLX_CLASS<ScalarType>,MV,OP> {
    
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

    //! Set sort type
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
    void SetSortType( const string which ) { _which = which; };
    
    //! Sort the vector of eigenvalues with respect to the chosen sorting type, optionally returning the permutation vector.
    /**
       @param solver [in] Eigensolver that is calling the sorting routine

       @param n [in] Size of the array

       @param evals [in/out] Array of length n containing the eigenvalues to be sorted

       @param perm [out] Vector of length n to store the permutation (optional)
    */
    void sort(Eigensolver<ANSZI_CPLX_CLASS<ScalarType>,MV,OP>* solver, int n, ANSZI_CPLX_CLASS<ScalarType> *evals, std::vector<int> *perm = 0) const;
    
    //! Sort the vectors of eigenpairs with respect to the chosen sorting type, optionally returning the permutation vector.
    /**
       @param solver [in] Eigensolver that is calling the sorting routine

       @param n [in] Size of the array

       @param r_evals [in/out] Array of length n containing the real part of the eigenvalues to be sorted 

       @param i_evals [in/out] Array of length n containing the imaginary part of the eigenvalues to be sorted 

       @param perm [out] Vector of length n to store the permutation (optional)
    */
    void sort(Eigensolver<ANSZI_CPLX_CLASS<ScalarType>,MV,OP>* solver, int n, ANSZI_CPLX_CLASS<ScalarType> *r_evals, 
		    ANSZI_CPLX_CLASS<ScalarType> *i_evals, std::vector<int> *perm = 0) const;
    
  protected: 
    
    //! Sorting type
    /*! \note Sorting choices:
       <ul>
       <li> "LM" - Largest Magnitude [ default ]
       <li> "SM" - Smallest Magnitude
       <li> "LR" - Largest Real 
       <li> "SR" - Smallest Real 
       <li> "LI" - Largest Imaginary 
       <li> "SI" - Smallest Imaginary 
       </ul>
    */
    string _which;

  };

  // Partial specialization for complex numbers templated on real type ScalarType
  template<class ScalarType, class MV, class OP>
  void BasicSort<ANSZI_CPLX_CLASS<ScalarType>,MV,OP>::sort(Eigensolver<ANSZI_CPLX_CLASS<ScalarType>,MV,OP>* solver, 
								 int n, ANSZI_CPLX_CLASS<ScalarType> *evals, std::vector<int> *perm) const 
  {
    int i=0, j=0;
    
    // Temp integer for swapping the index of the permutation, used in all sorting types.
    int tempord=0;
    
    // Temp variable for the magnitude of the ScalarType used in sorting "LM" and "SM".
    typename Teuchos::ScalarTraits<ANSZI_CPLX_CLASS<ScalarType> >::magnitudeType temp2;
    
    // Temp variable for swapping the eigenvalue used in all sorting types.
    ANSZI_CPLX_CLASS<ScalarType> temp;

    Teuchos::LAPACK<int,ANSZI_CPLX_CLASS<ScalarType> > lapack;

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
        if (perm) {
          tempord = (*perm)[j];
        }
        temp2 = Teuchos::ScalarTraits<complex<ScalarType> >::magnitude(evals[j]);
        for (i=j-1; i>=0 && Teuchos::ScalarTraits<complex<ScalarType> >::magnitude(evals[i])>temp2; --i) {
          evals[i+1]=evals[i];
          if (perm)
            (*perm)[i+1]=(*perm)[i];
        }
        evals[i+1] = temp; 
        if (perm) 
          (*perm)[i+1] = tempord;
      }
      return;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in increasing order of real part
    //---------------------------------------------------------------
    if (!_which.compare("SR")) {
      for (j=1; j < n; ++j) {
        temp = evals[j]; 
        if (perm)
          tempord = (*perm)[j];
        for (i=j-1; i>=0 && evals[i].real()>temp.real(); --i) {
          evals[i+1]=evals[i];
          if (perm)
            (*perm)[i+1]=(*perm)[i];
        }
        evals[i+1] = temp; 
        if (perm)
          (*perm)[i+1] = tempord;
      }
      return;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in increasing order of imagninary part
    //---------------------------------------------------------------
    if (!_which.compare("SI")) {
      for (j=1; j < n; ++j) {
        temp = evals[j]; 
        if (perm)
          tempord = (*perm)[j];
        for (i=j-1; i>=0 && evals[i].imag()>temp.imag(); --i) {
          evals[i+1]=evals[i];
          if (perm)
            (*perm)[i+1]=(*perm)[i];
        }
        evals[i+1] = temp; 
        if (perm)
          (*perm)[i+1] = tempord;
      }
      return;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in decreasing order of magnitude
    //---------------------------------------------------------------
    if (!_which.compare("LM")) {
      for (j=1; j < n; ++j) {
        temp = evals[j]; 
        if (perm)
          tempord = (*perm)[j];
        temp2 = Teuchos::ScalarTraits<ANSZI_CPLX_CLASS<ScalarType> >::magnitude(evals[j]);
        for (i=j-1; i>=0 && Teuchos::ScalarTraits<ANSZI_CPLX_CLASS<ScalarType> >::magnitude(evals[i])<temp2; --i) {
          evals[i+1]=evals[i];
          if (perm)
            (*perm)[i+1]=(*perm)[i];
        }
        evals[i+1] = temp; 
        if (perm)
          (*perm)[i+1] = tempord;
        }
      return;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in decreasing order of real part
    //---------------------------------------------------------------
    if (!_which.compare("LR")) {
      for (j=1; j < n; ++j) {
        temp = evals[j]; 
        if (perm)
        tempord = (*perm)[j];
        for (i=j-1; i>=0 && evals[i].real()<temp.real(); --i) {
          evals[i+1]=evals[i];
          if (perm)
            (*perm)[i+1]=(*perm)[i];
        }
        evals[i+1] = temp; 
        if (perm)
          (*perm)[i+1] = tempord;
      }
      return;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in decreasing order of imaginary part
    //---------------------------------------------------------------
    if (!_which.compare("LI")) {
      for (j=1; j < n; ++j) {
        temp = evals[j]; 
        if (perm)
        tempord = (*perm)[j];
        for (i=j-1; i>=0 && evals[i].real()<temp.real(); --i) {
          evals[i+1]=evals[i];
          if (perm)
            (*perm)[i+1]=(*perm)[i];
        }
        evals[i+1] = temp; 
        if (perm)
          (*perm)[i+1] = tempord;
      }
      return;
    }

    // The character string held by this class is not valid.  
    TEST_FOR_EXCEPTION(true, SortManagerError, 
                       "Anasazi::BasicSort::sort(): sorting order is not valid");
  }

  // Partial specialization for complex numbers templated on real type ScalarType
  // NOTE:  This implementation doesn't perform any operations on i_evals because r_evals is a vector of complex 
  //        scalar type and should hold all the approximate eigenvalues that need to be sorted.
  template<class ScalarType, class MV, class OP>
  void BasicSort<ANSZI_CPLX_CLASS<ScalarType>,MV,OP>::sort(Eigensolver<ANSZI_CPLX_CLASS<ScalarType>,MV,OP>* solver, 
								 int n, ANSZI_CPLX_CLASS<ScalarType> *r_evals, 
								 ANSZI_CPLX_CLASS<ScalarType> *i_evals,
								 std::vector<int> *perm) const 
  {
    int i=0, j=0;
    
    // Temp integer for swapping the index of the permutation, used in all sorting types.
    int tempord=0;
    
    // Temp variable for the magnitude of the ScalarType used in sorting "LM" and "SM".
    typename Teuchos::ScalarTraits<ANSZI_CPLX_CLASS<ScalarType> >::magnitudeType temp2;
    
    // Temp variable for swapping the eigenvalue used in all sorting types.
    ANSZI_CPLX_CLASS<ScalarType> temp;

    Teuchos::LAPACK<int,ANSZI_CPLX_CLASS<ScalarType> > lapack;

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
        temp = r_evals[j]; 
        if (perm) {
          tempord = (*perm)[j];
        }
        temp2 = Teuchos::ScalarTraits<complex<ScalarType> >::magnitude(r_evals[j]);
        for (i=j-1; i>=0 && Teuchos::ScalarTraits<complex<ScalarType> >::magnitude(r_evals[i])>temp2; --i) {
          r_evals[i+1]=r_evals[i];
          if (perm)
            (*perm)[i+1]=(*perm)[i];
        }
        r_evals[i+1] = temp; 
        if (perm) 
          (*perm)[i+1] = tempord;
      }
      return;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in increasing order of real part
    //---------------------------------------------------------------
    if (!_which.compare("SR")) {
      for (j=1; j < n; ++j) {
        temp = r_evals[j]; 
        if (perm)
          tempord = (*perm)[j];
        for (i=j-1; i>=0 && r_evals[i].real()>temp.real(); --i) {
          r_evals[i+1]=r_evals[i];
          if (perm)
            (*perm)[i+1]=(*perm)[i];
        }
        r_evals[i+1] = temp; 
        if (perm)
          (*perm)[i+1] = tempord;
      }
      return;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in increasing order of imagninary part
    //---------------------------------------------------------------
    if (!_which.compare("SI")) {
      for (j=1; j < n; ++j) {
        temp = r_evals[j]; 
        if (perm)
          tempord = (*perm)[j];
        for (i=j-1; i>=0 && r_evals[i].imag()>temp.imag(); --i) {
          r_evals[i+1]=r_evals[i];
          if (perm)
            (*perm)[i+1]=(*perm)[i];
        }
        r_evals[i+1] = temp; 
        if (perm)
          (*perm)[i+1] = tempord;
      }
      return;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in decreasing order of magnitude
    //---------------------------------------------------------------
    if (!_which.compare("LM")) {
      for (j=1; j < n; ++j) {
        temp = r_evals[j]; 
        if (perm)
          tempord = (*perm)[j];
        temp2 = Teuchos::ScalarTraits<ANSZI_CPLX_CLASS<ScalarType> >::magnitude(r_evals[j]);
        for (i=j-1; i>=0 && Teuchos::ScalarTraits<ANSZI_CPLX_CLASS<ScalarType> >::magnitude(r_evals[i])<temp2; --i) {
          r_evals[i+1]=r_evals[i];
          if (perm)
            (*perm)[i+1]=(*perm)[i];
        }
        r_evals[i+1] = temp; 
        if (perm)
          (*perm)[i+1] = tempord;
        }
      return;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in decreasing order of real part
    //---------------------------------------------------------------
    if (!_which.compare("LR")) {
      for (j=1; j < n; ++j) {
        temp = r_evals[j]; 
        if (perm)
        tempord = (*perm)[j];
        for (i=j-1; i>=0 && r_evals[i].real()<temp.real(); --i) {
          r_evals[i+1]=r_evals[i];
          if (perm)
            (*perm)[i+1]=(*perm)[i];
        }
        r_evals[i+1] = temp; 
        if (perm)
          (*perm)[i+1] = tempord;
      }
      return;
    }
    //---------------------------------------------------------------
    // Sort eigenvalues in decreasing order of imaginary part
    //---------------------------------------------------------------
    if (!_which.compare("LI")) {
      for (j=1; j < n; ++j) {
        temp = r_evals[j]; 
        if (perm)
        tempord = (*perm)[j];
        for (i=j-1; i>=0 && r_evals[i].real()<temp.real(); --i) {
          r_evals[i+1]=r_evals[i];
          if (perm)
            (*perm)[i+1]=(*perm)[i];
        }
        r_evals[i+1] = temp; 
        if (perm)
          (*perm)[i+1] = tempord;
      }
      return;
    }

    // The character string held by this class is not valid.  
    TEST_FOR_EXCEPTION(true, SortManagerError, 
                       "Anasazi::BasicSort::sort(): sorting order is not valid");
  }

#endif // ( defined(HAVE_COMPLEX) || defined(HAVE_COMPLEX_H) ) && defined(HAVE_TEUCHOS_COMPLEX)

  
} // namespace Anasazi

#endif // ANASAZI_BASIC_SORT_HPP


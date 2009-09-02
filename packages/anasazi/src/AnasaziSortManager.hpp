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

#ifndef ANASAZI_SORTMANAGER_HPP
#define ANASAZI_SORTMANAGER_HPP

/*!     \file AnasaziSortManager.hpp
        \brief Virtual base class which defines the interface between an eigensolver and a class whose
	job is the sorting of the computed eigenvalues
*/

/*!    \class Anasazi::SortManager

       \brief Anasazi's templated pure virtual class for managing the sorting
       of approximate eigenvalues computed by the eigensolver. A concrete
       implementation of this class is necessary. 

       \author Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_ParameterList.hpp"



namespace Anasazi {

  //! @name LOBPCG Exceptions
  //@{ 
  /** \brief SortManagerError is thrown when the Anasazi::SortManager is unable to sort the numbers,
   *  due to some failure of the sort method or error in calling it.
   */
  class SortManagerError : public AnasaziError
  {public: SortManagerError(const std::string& what_arg) : AnasaziError(what_arg) {}};

  //@}

  template<class MagnitudeType>
  class SortManager {
    
  public:

    //! Default constructor.
    SortManager() {};
    
    //! Constructor accepting a Teuchos::ParameterList. This is the default mode for instantiating a SortManager.
    SortManager(Teuchos::ParameterList &pl) {};

    //! Destructor
    virtual ~SortManager() {};

    /*! \brief Sort real eigenvalues, optionally returning the permutation vector.

       @param evals [in/out] Vector of length at least \c n containing the eigenvalues to be sorted.  <br>
                     On output, the first \c n eigenvalues will be sorted. The rest will be unchanged.

       @param perm [out] Vector of length at least \c n to store the permutation index (optional).  <br>
       If specified, on output the first \c n eigenvalues will contain the permutation indices, in the range <tt>[0,n-1]</tt>, such that <tt>evals_out[i] = evals_in[perm[i]]</tt>

       @param n [in] Number of values in evals to be sorted. If <tt>n == -1</tt>, all values will be sorted.
    */
    virtual void sort(std::vector<MagnitudeType> &evals, Teuchos::RCP<std::vector<int> > perm = Teuchos::null, int n = -1) const = 0;

    /*! \brief Sort complex eigenvalues, optionally returning the permutation vector.

       This routine takes two vectors, one for each part of a complex
       eigenvalue. This is helpful for solving real, non-symmetric eigenvalue
       problems.

       @param r_evals [in/out] Vector of length at least \c n containing the real part of the eigenvalues to be sorted.  <br>
                     On output, the first \c n eigenvalues will be sorted. The rest will be unchanged.

       @param i_evals [in/out] Vector of length at least \c n containing the imaginary part of the eigenvalues to be sorted.  <br>
                     On output, the first \c n eigenvalues will be sorted. The rest will be unchanged.

       @param perm [out] Vector of length at least \c n to store the permutation index (optional).  <br>
       If specified, on output the first \c n eigenvalues will contain the permutation indices, in the range <tt>[0,n-1]</tt>, such that <tt>r_evals_out[i] = r_evals_in[perm[i]]</tt>
       and similarly for \c i_evals.

       @param n [in] Number of values in \c r_evals, \c i_evals to be sorted. If <tt>n == -1</tt>, all values will be sorted.
    */
    virtual void sort(std::vector<MagnitudeType> &r_evals, 
                      std::vector<MagnitudeType> &i_evals, 
                      Teuchos::RCP<std::vector<int> > perm = Teuchos::null,
                      int n = -1) const = 0;
    
  };
  
}

#endif // ANASAZI_SORTMANAGER_HPP


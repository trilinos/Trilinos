// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_FAD_FADVECTOR_HPP
#define SACADO_FAD_FADVECTOR_HPP

#include <vector>
#include "Sacado_Fad_DVFad.hpp"

namespace Sacado {

  namespace Fad {

    /*!
     * \brief A class for storing a contiguously allocated array of Fad
     * objects where the values and derivative arrays for each Fad object
     * are stored in contiguous memory.  To preserve this structure, many
     * vector operations aren't supported (like resizing).
     */
    template <typename ValueT, 
	      typename ScalarT = typename ScalarValueType<ValueT>::type >
    class FadVector {
    public:

      //! Synonym for Fad type
      typedef Sacado::Fad::DVFad<ValueT,ScalarT> FadType;

      //! Constructor 
      FadVector(int vec_size, int deriv_size);

      //! Copy constructor
      FadVector(const FadVector& fv);

      //! Destructor
      ~FadVector();

      //! Assignment
      FadVector& operator=(const FadVector& fv);

      //! Vector size
      int size() const { return vec_.size(); }

      //! Derivative size
      int deriv_size() const { return deriv_size_; }

      //! Array access
      FadType& operator[] (int i) { return vec_[i]; }

      //! Array access
      const FadType& operator[](int i) const { return vec_[i]; }

      //! Pointer to values
      ValueT* vals();

      //! Pointer to values
      const ValueT* vals() const;

      //! Pointer to derivatives
      ValueT* dx();

      //! Pointer to values
      const ValueT* dx() const;

    protected:

      //! Size of derivative array
      int deriv_size_;

      //! Vector of Fad's
      std::vector<FadType> vec_;

    }; // class FadVector

  } // namespace Fad

} // namespace Sacado

#include "Sacado_Fad_FadVectorImp.hpp"

#endif // SACADO_FAD_FADVECTOR_HPP

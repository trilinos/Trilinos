// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_VECTOR_HPP
#define SACADO_FAD_VECTOR_HPP

#include <vector>
#include "Sacado_Fad_DVFad.hpp"

namespace Sacado {

  namespace Fad {

    enum VectorDerivOrientation {
      Row,      //! Derivatives are stored row-wise (strided)
      Column    //! Derivatives ar stored column-wise (unit stride)
    };

    /*!
     * \brief A class for storing a contiguously allocated array of Fad
     * objects.  This is a general definition that will work for all Fad types,
     * and is merely a wrapper around std::vector.  A specialization for
     * Sacado::Fad::DVFad providing contiguous allocation of values and 
     * derivatives is below.
     */
    template <typename OrdinalType, typename FadType >
    class Vector {
    public:

      //! Typename of values
      typedef typename Sacado::ValueType<FadType>::type ValueType;

      //! Constructor 
      Vector(OrdinalType vec_size, OrdinalType deriv_sz,
	     VectorDerivOrientation orient = Row) :
	deriv_size_(deriv_sz), vec_(vec_size) {
	for (OrdinalType i=0; i<vec_size; i++)
	  vec_[i] = FadType(deriv_size_, ValueType(0.0));
      }

      //! Copy constructor
      Vector(const Vector& fv) : deriv_size_(fv.deriv_size_), vec_(fv.vec_) {}

      //! Destructor
      ~Vector() {}

      //! Assignment
      Vector& operator=(const Vector& fv) {
	deriv_size_ = fv.deriv_size_;
	vec_ = fv.vec_;
	return *this;
      }

      //! Vector size
      OrdinalType size() const { return vec_.size(); }

      //! Derivative size
      OrdinalType deriv_size() const { return deriv_size_; }

      //! Derivative array stride
      OrdinalType deriv_stride() const { return 1; }

      //! Derivative array orientation
      VectorDerivOrientation deriv_orientation() const { return Column; }

      //! Array access
      FadType& operator[] (OrdinalType i) { return vec_[i]; }

      //! Array access
      const FadType& operator[](OrdinalType i) const { return vec_[i]; }

    protected:

      //! Size of derivative array
      OrdinalType deriv_size_;

      //! Vector of Fad's
      std::vector<FadType> vec_;

    }; // class Vector

    /*!
     * \brief A class for storing a contiguously allocated array of Fad
     * objects where the values and derivative arrays for each Fad object
     * are stored in contiguous memory.  To preserve this structure, many
     * vector operations aren't supported (like resizing).
     */
    template <typename OrdinalType, typename ValueType>
    class Vector< OrdinalType, Sacado::Fad::DVFad<ValueType> > {
    public:

      //! Synonym for Fad type
      typedef Sacado::Fad::DVFad<ValueType> FadType;

      //! Constructor 
      Vector(OrdinalType vec_size, OrdinalType deriv_size, 
	     VectorDerivOrientation orient = Row);

      //! Copy constructor
      Vector(const Vector& fv);

      //! Destructor
      ~Vector();

      //! Assignment
      Vector& operator=(const Vector& fv);

      //! Vector size
      OrdinalType size() const { return vec_.size(); }

      //! Derivative size
      OrdinalType deriv_size() const { return deriv_size_; }

      //! Derivative array stride
      OrdinalType deriv_stride() const { return stride_; }

      //! Derivative array orientation
      VectorDerivOrientation deriv_orientation() const { return orient_; }

      //! Array access
      FadType& operator[] (OrdinalType i) { return vec_[i]; }

      //! Array access
      const FadType& operator[](OrdinalType i) const { return vec_[i]; }

      //! Pointer to values
      ValueType* vals();

      //! Pointer to values
      const ValueType* vals() const;

      //! Pointer to derivatives
      ValueType* dx();

      //! Pointer to values
      const ValueType* dx() const;

    protected:

      //! Size of derivative array
      OrdinalType deriv_size_;

      //! Derivative array orientation
      VectorDerivOrientation orient_;

      //! Derivative array stride
      OrdinalType stride_;

      //! Vector of Fad's
      std::vector<FadType> vec_;

    }; // class Vector

  } // namespace Fad

} // namespace Sacado

#include "Sacado_Fad_VectorImp.hpp"

#endif // SACADO_FAD_VECTOR_HPP

// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#include "Thyra_VectorBase.hpp"

#ifndef THYRA_EXPLICIT_VECTOR_VIEW_HPP
#define THYRA_EXPLICIT_VECTOR_VIEW_HPP

namespace Thyra {

/** \brief Create an explicit non-mutable (const) view of a <tt>VectorBase</tt> object.
 *
 * This utility class makes it easy to explicitly access a contiguous
 * subset of elements in any <tt>const VectorBase</tt> object.
 *
 * <b>Warning!</b> Creating an explicit view of an arbitrary
 * <tt>%VectorBase</tt> object may be a very expensive operation (such as
 * with distributed-memory vectors) and should only be done in special
 * cases (such as when the vector is an in-core vector).  There
 * several specialized use cases where creating these types of
 * explicit views are necessary but in most cases this should not be
 * done.
 *
 * The following functions show four different ways to extract a range
 * of elements from in any <tt>const VectorBase</tt> object and copy
 * then into a raw C++ array:
 *
 \code

 //
 // Copy using unit-strided pointer access
 //
 template<class Scalar>
 void copyToArrayPointerUnit(
   const Thyra::VectorBase<Scalar>   &v
   ,const Thyra::Range1D             &rng
   ,Scalar                           x[]   // Size == rng.size()
   )
 {
    Thyra::ExplicitVectorView<Scalar> v_ev(v,rng,true);  // Force unit stride
    const Scalar *v_ptr = v_ev.values();                 // Get pointer to unit-stride data
    for( int k = 0; k < n; ++k )                         // For each element in view:
      x[k] = v_ptr[k];                                   //   Copy elements
    // When this function returns then v_ev will be destroyed and the view will be freed
 }

 //
 // Copy using non-unit-strided pointer access
 //
 template<class Scalar>
 void copyToArrayPointerNonunit(
   const Thyra::VectorBase<Scalar>   &v
   ,const Thyra::Range1D             &rng
   ,Scalar                           x[]   // Size == rng.size()
   )
 {
    Thyra::ExplicitVectorView<Scalar> v_ev(v,rng);           // Allow non-unit stride
    const Scalar        *v_ptr    = v_ev.values();           // Get pointer to non-unit-stride data
    const Thyra::Index  *v_stride = v_ev.stride();           // Get stride between vector data
    for( int k = 0; k < rng.size(); ++k, v_ptr += v_stride ) // For each element in view:
      x[k] = *v_ptr;                                         //   Copy elements
    // When this function returns then v_ev will be destroyed and the view will be freed
 }

 //
 // Copy using possibly non-unit stride with zero-based overloaded function operator[]()
 //
 template<class Scalar>
 void copyToArrayZeroBasedOperator(
   const Thyra::VectorBase<Scalar>   &v
   ,const Thyra::Range1D             &rng
   ,Scalar                           x[]   // Size == rng.size()
   )
 {
    Thyra::ExplicitVectorView<Scalar> v_ev(v,rng);       // Allow non-unit stride 
    for( int k = 0; k < rng.size(); ++k )                // For each element in view:
      x[k] = v_ev[k];                                    //   Copy elements using operator[]()
    // When this function returns then v_ev will be destroyed and the view will be freed
 }

 //
 // Copy using possibly non-unit stride with one-based overloaded function operator()()
 //
 template<class Scalar>
 void copyToArrayOneBasedOperator(
   const Thyra::VectorBase<Scalar>   &v
   ,const Thyra::Range1D             &rng
   ,Scalar                           x[]   // Size == rng.size()
   )
 {
    Thyra::ExplicitVectorView<Scalar> v_ev(v,rng);       // Allow non-unit stride 
    for( int k = 0; k < rng.size(); ++k )                // For each element in view:
      x[k] = v_ev(k+1);                                  //   Copy elements using operator()()
    // When this function returns then v_ev will be destroyed and the view will be freed
 }

 \endcode
 *
 * If one wants to modify the elements in a <tt>VectorBase</tt> object
 * then one should use the utility class
 * <tt>ExplicitMutableVectorView</tt>.
 *
 * The default constructor, copy constructor and assignment operators
 * are not allowed.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class ExplicitVectorView {
public:
  /** \brief Construct an explicit non-mutable (const) view of a subset of elements.
   *
   * @param  v     [in] The vector that a view will be taken.  This object must be maintained
   *               until <tt>*this</tt> is destroyed.
   * @param  rng   [in] the range of element indices that the explicit view will be taken.
   * @param  forceUnitStride
   *               [in] If <tt>true</tt> then the view will have unit stride.
   *
   * Preconditions:<ul>
   * <li>[<tt>rng.full_range()==false</tt>] <tt>rng.ubound() <= v.space()->dim()</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li><tt>this->sv()</tt> returns the created view
   * <li><tt>this->globalOffset()==rng.lbound()-1</tt>
   * <li><tt>this->subDim()==rng.size()</tt>
   * <li><tt>this->values()</tt> returns a pointer to a <tt>Scalar</tt> array
   * <li><tt>this->stride()</tt> returns the stride between the elements pointed it in <tt>this->values()</tt>
   * <li>[<tt>forceUnitStride==true</tt>] <tt>this->stride()==1</tt>
   * </ul>
   */
  ExplicitVectorView( const VectorBase<Scalar>& v, const Range1D &rng = Range1D(), const bool forceUnitStride = false )
    :v_(v)
    {
      v_.getSubVector(rng,&sv_s_);
      if( forceUnitStride && sv_s_.stride() != 1 ) {
        Scalar *values = new Scalar[sv_s_.subDim()];
        RTOp_index_type i; const Scalar *sv_v;
        for( sv_v = sv_s_.values(), i=0; i < sv_s_.subDim(); ++i, sv_v += sv_s_.stride() )
          values[i] = *sv_v;
        sv_.initialize(sv_s_.globalOffset(),sv_s_.subDim(),values,1);
      }
      else {
        sv_ = sv_s_;
      }
    }
  /// Free the explicit view on the <tt>VectorBase</tt> object <tt>v</tt> passed to <tt>ExplicitVectorView()</tt>
  ~ExplicitVectorView()
    {
      if( sv_s_.stride() != sv_.stride() )
        delete [] const_cast<Scalar*>(sv_.values());
      v_.freeSubVector(&sv_s_);
    }
  /// Returns the explict view as an <tt>RTOpPack::SubVectorT<Scalar></tt> object
  const RTOpPack::SubVectorT<Scalar>& sv() const { return sv_; }
  /// Returns the global offset for the explicit view
  RTOp_index_type   globalOffset() const { return sv_.globalOffset(); }
  /// Returns the dimension of the explicit view
  RTOp_index_type   subDim()       const { return sv_.subDim();  }
  /// Return a pointer to a <tt>Scalar</tt> array containing the explicit view
  const Scalar*     values()       const { return sv_.values();  }
  /// Return the stride between elements in the array returned from <tt>this->values()</tt>
  ptrdiff_t         stride()       const { return sv_.stride();  }
  /// Zero-based indexing: Preconditions: <tt>values()!=NULL && (0 <= i <= subDim()-1)</tt>
  const Scalar& operator[](RTOp_index_type i) const { return sv_[i]; }
  /// One-based indexing: Preconditions: <tt>values()!=NULL && (1 <= i <= subDim())</tt>
  const Scalar& operator()(RTOp_index_type i) const { return sv_(i); }
private:
  const VectorBase<Scalar>          &v_;
  RTOpPack::SubVectorT<Scalar>  sv_s_;
  RTOpPack::SubVectorT<Scalar>  sv_;
  // Not defined and not to be called
  ExplicitVectorView();
  ExplicitVectorView(const ExplicitVectorView<Scalar>&);
  ExplicitVectorView<Scalar>& operator==(const ExplicitVectorView<Scalar>&);
};
 
/** \brief Create an explicit mutable (non-const) view of a <tt>VectorBase</tt> object.
 *
 * This utility class makes it easy to explicitly access a contiguous
 * subset of elements in any <tt>VectorBase</tt> object and change the
 * <tt>VectorBase object</tt>.
 *
 * <b>Warning!</b> Creating an explicit view of an arbitrary
 * <tt>%VectorBase</tt> object may be a very expensive operation (such as
 * with distributed-memory vectors) and should only be done in special
 * cases (such as when the vector is an in-core vector).  There
 * several specialized use cases where creating these types of
 * explicit views are necessary but in most cases this should not be
 * done.
 *
 * The following functions show four different ways to access a range
 * of elements from in any <tt>VectorBase</tt> object and then add to them
 * the values form a raw C++ array.
 *
 \code

 //
 // Add-to using unit-strided pointer access
 //
 template<class Scalar>
 void addToArrayPointerUnit(
   const Thyra::Range1D            &rng
   ,const Scalar                   x[]   // Size == rng.size()
   ,Thyra::VectorBase<Scalar>      *v
   )
 {
    Thyra::ExplicitMutableVectorView<Scalar> v_ev(rng,*v,true);  // Force unit stride
    const Scalar *v_ptr = v_ev.values();                         // Get pointer to unit-stride data
    for( int k = 0; k < n; ++k )                                 // For each element in view:
      v_ptr[k] += x[k];                                          //   Add-to elements
    // When this function returns then v_ev will be destroyed and the view will be commited back and *v modified
 }

 //
 // Add-to using non-unit-strided pointer access
 //
 template<class Scalar>
 void addToArrayPointerNonunit(
   const Thyra::Range1D          &rng
   ,const Scalar                 x[]   // Size == rng.size()
   ,Thyra::VectorBase<Scalar>    *v
   )
 {
    Thyra::ExplicitMutableVectorView<Scalar> v_ev(rng,*v);   // Allow non-unit stride
    const Scalar          *v_ptr    = v_ev.values();         // Get pointer to non-unit-stride data
    const Thyra::Index  *v_stride = v_ev.stride();           // Get stride between vector data
    for( int k = 0; k < rng.size(); ++k, v_ptr += v_stride ) // For each element in view:
      *v_ptr + x[k];                                         //   Add-to elements
    // When this function returns then v_ev will be destroyed and the view will be commited back and *v modified
 }

 //
 // Add-to using possibly non-unit stride with zero-based overloaded function operator[]()
 //
 template<class Scalar>
 void addToArrayZeroBasedOperator(
   const Thyra::Range1D          &rng
   ,const Scalar                 x[]   // Size == rng.size()
   ,Thyra::VectorBase<Scalar>    *v
   )
 {
    Thyra::ExplicitMutableVectorView<Scalar> v_ev(rng,*v);   // Allow non-unit stride 
    for( int k = 0; k < rng.size(); ++k )                    // For each element in view:
      v_ev[k] += x[k];                                       //   Add-to elements using operator[]()
    // When this function returns then v_ev will be destroyed and the view will be commited back and *v modified
 }

 //
 // Add-to using possibly non-unit stride with one-based overloaded function operator()()
 //
 template<class Scalar>
 void addToArrayOneBasedOperator(
   const Thyra::Range1D          &rng
   ,const Scalar                 x[]   // Size == rng.size()
   ,Thyra::VectorBase<Scalar>    *v
   )
 {
    Thyra::ExplicitMutableVectorView<Scalar> v_ev(rng,*v); // Allow non-unit stride 
    for( int k = 0; k < rng.size(); ++k )                  // For each element in view:
      v_ev(k+1) += x[k];                                   //   Add-to elements using operator()()
    // When this function returns then v_ev will be destroyed and the view will be commited back and *v modified
 }

 \endcode
 *
 * If one wants to only read the elements in a <tt>VectorBase</tt> object
 * then one should use the utility class <tt>ExplicitVectorView</tt>.
 *
 * The default constructor, copy constructor and assignment operators
 * are not allowed.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class ExplicitMutableVectorView {
public:
  /** \brief Construct an explicit mutable (non-const) view of a subset of elements.
   *
   * @param  v     [in] The vector that a view will be taken.  This object must be maintained
   *               until <tt>*this</tt> is destroyed.  The elements in <tt>v</tt> are not
   *               guaranteed to be updated until <tt>*this</tt> is destroyed.
   * @param  rng   [in] the range of element indices that the explicit view will be taken.
   * @param  forceUnitStride
   *               [in] If <tt>true</tt> then the view will have unit stride.
   *
   * Preconditions:<ul>
   * <li>[<tt>rng.full_range()==false</tt>] <tt>rng.ubound() <= v.space()->dim()</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li><tt>this->sv()</tt> returns the created view
   * <li><tt>this->globalOffset()==rng.lbound()-1</tt>
   * <li><tt>this->subDim()==rng.size()</tt>
   * <li><tt>this->values()</tt> returns a pointer to a <tt>Scalar</tt> array
   * <li><tt>this->stride()</tt> returns the stride between the elements pointed it in <tt>this->values()</tt>
   * <li>[<tt>forceUnitStride==true</tt>] <tt>this->stride()==1</tt>
   * </ul>
   */
  ExplicitMutableVectorView( VectorBase<Scalar>& v, const Range1D &rng = Range1D(), const bool forceUnitStride = false )
    :v_(v)
    {
      v_.getSubVector(rng,&sv_s_);
      if( forceUnitStride && sv_s_.stride() != 1 ) {
        Scalar *values = new Scalar[sv_s_.subDim()];
        RTOp_index_type i; const Scalar *sv_v;
        for( sv_v = sv_s_.values(), i=0; i < sv_s_.subDim(); ++i, sv_v += sv_s_.stride() )
          values[i] = *sv_v;
        sv_.initialize(sv_s_.globalOffset(),sv_s_.subDim(),values,1);
      }
      else {
        sv_ = sv_s_;
      }
    }
  /// Commits back the the explicit view on the <tt>VectorBase</tt> object <tt>v</tt> passed to <tt>ExplicitMutableVectorView()</tt>
  ~ExplicitMutableVectorView()
    {
      if( sv_s_.stride() != sv_.stride() ) {
        RTOp_index_type i; Scalar *sv_v; const Scalar *values;
        for( sv_v = sv_s_.values(), values = sv_.values(), i=0; i < sv_s_.subDim(); ++i, sv_v += sv_s_.stride() )
          *sv_v = *values++;
        delete [] sv_.values();
      }
      v_.commitSubVector(&sv_s_);
    }
  /// Returns the explict view as an <tt>RTOpPack::SubVectorT<Scalar></tt> object
  const RTOpPack::MutableSubVectorT<Scalar>& sv() const { return sv_; }
  /// Returns the global offset for the explicit view
  RTOp_index_type   globalOffset() const { return sv_.globalOffset(); }
  /// Returns the dimension of the explicit view
  RTOp_index_type   subDim()       const { return sv_.subDim();  }
  /// Return a pointer to a <tt>Scalar</tt> array containing the explicit view
  Scalar*           values()       const { return sv_.values();  }
  /// Return the stride between elements in the array returned from <tt>this->values()</tt>
  ptrdiff_t         stride()       const { return sv_.stride();  }
  /// Zero-based indexing: Preconditions: <tt>values()!=NULL && (0 <= i <= subDim()-1)</tt>
  Scalar& operator[](RTOp_index_type i) const { return sv_[i]; }
  /// One-based indexing: Preconditions: <tt>values()!=NULL && (1 <= i <= subDim())</tt>
  Scalar& operator()(RTOp_index_type i) const { return sv_(i); }
private:
  VectorBase<Scalar>                       &v_;
  RTOpPack::MutableSubVectorT<Scalar>  sv_s_;
  RTOpPack::MutableSubVectorT<Scalar>  sv_;
  // Not defined and not to be called
  ExplicitMutableVectorView();
  ExplicitMutableVectorView(const ExplicitMutableVectorView<Scalar>&);
  ExplicitMutableVectorView<Scalar>& operator==(const ExplicitMutableVectorView<Scalar>&);
};

} // namespace Thyra

#endif // THYRA_EXPLICIT_VECTOR_VIEW_HPP

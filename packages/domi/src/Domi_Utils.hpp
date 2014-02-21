// @HEADER
// ***********************************************************************
//
//            Domi: Multidimensional Datastructures Package
//                 Copyright (2013) Sandia Corporation
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef DOMI_UTILS_HPP
#define DOMI_UTILS_HPP

// Domi includes
#include "Domi_ConfigDefs.hpp"

// Teuchos includes
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"

#ifdef HAVE_MPI
// MPI include
#include "mpi.h"
#endif

namespace Domi
{

/** \brief Used to specify the order in which elements of an MDArray,
 *         MDArrayView, or MDArrayRCP are stored.
 *
 *  Note that there are only two orderings supported, but multiple
 *  ways to refer to them.
 */
enum ELayout
{
  /** \brief C order, in which the last index varies fastest */
  C_ORDER             = 0,
  /** \brief Fortran order, in which the first index varies fastest */
  FORTRAN_ORDER       = 1,
  /** \brief Row-major order, in which the last index varies fastest */
  ROW_MAJOR           = 0,
  /** \brief Column-major order, in which the first index varies fastest */
  COLUMN_MAJOR        = 1,
  /** \brief Last index fastest order */
  LAST_INDEX_FASTEST  = 0,
  /** \brief First index fastest order */
  FIRST_INDEX_FASTEST = 1,
  /** \brief Default order, currently first index fastest */
  DEFAULT_ORDER       = 1
};

////////////////////////////////////////////////////////////////////////

/** \brief Provide capability to declare a variable as non-const, even
 *         if template parameter is const
 *
 * Note that this capability is provided in the C++11 standard via
 * std::remove_const<>, but that we cannot asssume that we are using a
 * C++11 compiler.
 */
template< class T >
struct remove_const
{
  /** \brief Typedef for the template parameter
   */
  typedef T type;
};

/** \brief Specialization of remove_const when template parameter is
 *         already const
 */
template< class T >
struct remove_const< const T >
{
  /** \brief Typedef for the non-const version of the template
   *         parameter
   */
  typedef T type;
};

////////////////////////////////////////////////////////////////////////

/** \brief Compute the strides of an <tt>MDArray</tt>,
 *         <tt>MDArrayView</tt>, or <tt>MDArrayRCP</tt>, given the
 *         dimensions (as an Array) and the storage order.
 */
template< typename T >
Teuchos::Array< T > computeStrides(const Teuchos::ArrayView< T > & dimensions,
                                   const ELayout layout)
{
  typedef typename Teuchos::Array< T >::size_type size_type;
  size_type n = dimensions.size();
  Teuchos::Array< T > strides(n);
  if (layout == FIRST_INDEX_FASTEST)
  {
    strides[0] = 1;
    for (size_type axis = 1; axis < n; ++axis)
      strides[axis] = strides[axis-1] * dimensions[axis-1];
  }
  else
  {
    strides[n-1] = 1;
    for (size_type axis = n-2; axis >= 0; --axis)
      strides[axis] = strides[axis+1] * dimensions[axis+1];
  }
  return strides;
}

////////////////////////////////////////////////////////////////////////

/** \brief Compute the strides of an <tt>MDArray</tt>,
 *         <tt>MDArrayView</tt>, or <tt>MDArrayRCP</tt>, given its
 *         dimensions as an ArrayView.
 */
template< typename T >
Teuchos::Array< T >
computeStrides(const Teuchos::Array< T > & dimensions,
               const ELayout layout)
{
  // In the MDArray<T>(const MDArrayView<T> &) constructor, I try to
  // pass the MDArrayView dimensions to computeStrides(), but they
  // come in as ArrayView<const T> (for reasons I can't determine) and
  // cause all sorts of const-correctness problems.  So I copy them
  // into a new Array<T> and pass its view to the main
  // computeStrides() function.  Fortunately, the array of dimensions
  // is small.
  Teuchos::Array< T > nonConstDims(0);
  nonConstDims.insert(nonConstDims.begin(),
                      dimensions.begin(),
                      dimensions.end());
  return computeStrides(nonConstDims(), layout);
}

////////////////////////////////////////////////////////////////////////

/** \brief Compute the minimum size required for an <tt>MDArray</tt>,
 *         <tt>MDArrayView</tt>, or <tt>MDArrayRCP</tt>, given its
 *         dimensions as an Arrayview.
 */
template< typename T >
T computeSize(const Teuchos::ArrayView< T > & dimensions)
{
  typedef typename Teuchos::ArrayView< T >::size_type size_type;
  T result = 1;
  for (size_type axis = 0; axis < dimensions.size(); ++axis)
    result *= dimensions[axis];
  return result;
}

////////////////////////////////////////////////////////////////////////

/** \brief Compute the minimum size required for an <tt>MDArray</tt>,
 *         <tt>MDArrayView</tt>, or <tt>MDArrayRCP</tt>, given its
 *         dimensions as an Array.
 */
template< typename T >
T computeSize(const Teuchos::Array< T > & dimensions)
{
  // In the MDArray<T>(const MDArrayView<T> &) constructor, I try to
  // pass the MDArrayView dimensions to computeSize(), but they come
  // in as ArrayView<const T> (for reasons I can't determine) and
  // cause all sorts of const-correctness problems.  So I copy them
  // into a new Array<T> and pass its view to the main computeSize()
  // function.  Fortunately, the array of dimensions is small.
  Teuchos::Array< T > nonConstDims(0);
  nonConstDims.insert(nonConstDims.begin(),
                      dimensions.begin(),
                      dimensions.end());
  return computeSize(nonConstDims());
}

////////////////////////////////////////////////////////////////////////

/** \brief Compute the minimum size required for an <tt>MDArray</tt>,
 *         <tt>MDArrayView</tt>, or <tt>MDArrayRCP</tt>, given its
 *         dimensions and strides.
 */
template< typename T >
T computeSize(const Teuchos::ArrayView< T > & dimensions,
	      const Teuchos::ArrayView< T > & strides)
{
  typedef typename Teuchos::ArrayView< T >::size_type size_type;
  // T might be a const type, but we need result to be non-const
  typename remove_const< T >::type result = 1;
  for (size_type axis = 0; axis < dimensions.size(); ++axis)
    result += (dimensions[axis]-1) * strides[axis];
  return result;
}

////////////////////////////////////////////////////////////////////////

/** \brief Compute a valid axisCommSizes array, given the number of
 *         processors, the number of dimensions, and a candidate
 *         axisCommSizes array.
 *
 *  The candidate array can have fewer entries than the number of
 *  dimensions -- this function will fill in the extra ones.  It can
 *  also have negative values -- this function will substitute valid
 *  processor counts in place of negative values.  If a valid array is
 *  not possible, i.e. the number of processors is not evenly
 *  divisible into the number of dimensions or the given concrete axis
 *  sizes, then an exception will be raised.
 *
 *  The current algorithm is not terribly smart -- the first negative
 *  axis size encountered will result in all the available processors
 *  being assigned to that axis and any further negative axis sizes
 *  will be given a single process.
 */
Teuchos::Array< int >
regularizeAxisSizes(int numProcs,
                    int numDims,
                    const Teuchos::ArrayView< int > & axisCommSizes);

////////////////////////////////////////////////////////////////////////

/** \brief Compute the axis ranks for a given processor rank, given
 *         the number of processors along each axis.
 */
Teuchos::Array< int >
computeAxisRanks(int rank,
                 const Teuchos::ArrayView< int > & axisCommSizes);

////////////////////////////////////////////////////////////////////////

/** \brief Compute the axis ranks for a given processor rank, given
 *         the offset and the processors strides along each axis.
 */
Teuchos::Array< int >
computeAxisRanks(int rank,
                 int offset,
                 const Teuchos::ArrayView< int > & axisStrides);

////////////////////////////////////////////////////////////////////////

/** \brief Compute the array of periodic flags, given the user input.
 *
 * If the input array is shorter than the number of dimensions, fill
 * the remainder of the result with zeros.
 */
Teuchos::Array< int >
computePeriodic(int numDims,
                const Teuchos::ArrayView< int > & periodic);

////////////////////////////////////////////////////////////////////////

/** \brief Given a std::string which contains comma-separated integers,
 *         return an array of ints.
 */
void splitStringOfIntsWithCommas(std::string data,
                                 Teuchos::Array< int > & result);

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_MPI

/** \brief Return an MPI_Datatype, given its corresponding C/C++ type,
           specified as template parameter T
 */
template< class T > MPI_Datatype mpiType();

#endif

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_MPI

/** \brief Return an MPI flag for data layout, given its corresponding
           Domi enumeration
 */
int mpiOrder(ELayout layout);

#endif

}       // End Domi namespace

#endif	// DOMI_MDARRAY_UTILS_HPP

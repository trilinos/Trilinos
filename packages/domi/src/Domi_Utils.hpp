// @HEADER
// ***********************************************************************
//
//     Domi: Multi-dimensional Distributed Linear Algebra Services
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
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
#include "Teuchos_ParameterList.hpp"

#ifdef HAVE_MPI
// MPI include
#include "mpi.h"
#endif

namespace Domi
{

/** \name Domi types */
//@{

////////////////////////////////////////////////////////////////////////

/** \brief The Domi size_type has the capacity to describe the entire
 *         size of the allocated buffer
 *
 * Data buffers in Domi are stored within MDArray, MDArrayView, and
 * MDArrayRCP objects.  These in turn use std::vector to allocate
 * their data.  Thus, by default, the Domi::size_type is the same type
 * as std:vector<T>::size_type, although it can be changed at
 * configuration time by setting Teuchos_ORDINAL_TYPE.
 */
typedef Teuchos::Ordinal size_type;

/** \brief The Domi dim_type is the ordinal type used by dimensions
 *         and indexes, both local and global
 *
 * The default type for dim_type is int.  This means that the maximum
 * size of a multi-dimensional Domi data structure is (2**31)**n =
 * ~2.1B**n, where n is the number of dimensions (or the maximum value
 * of size_type, whichever is smaller).  This should be sufficient for
 * the largest problems under consideration for the forseeable future.
 * MPI arrays also use ints for the same reason.  The type specified
 * by dim_type can be altered at configuration time by setting
 * Domi_ORDINAL_TYPE.
 */
typedef Ordinal dim_type;

/** \brief The Domi difference_type is the same as the dim_type, not
 *         the size_type 
*/
typedef Ordinal difference_type;

/** \brief The Domi padding type simply encapsulates a lower padding
 *         and upper padding along a single axis
 *
 * This is simply a Teuchos Tuple of length 2, with element 0
 * referencing the lower pad and element 1 referencing the upper pad.
 */
typedef Teuchos::Tuple< int, 2 > padding_type;

////////////////////////////////////////////////////////////////////////

/** \brief Layout enumeration, used to specify the order in which
 *         elements of an MDArray, MDArrayView, or MDArrayRCP are
 *         stored.
 *
 *  Note that there are only two orderings supported, but multiple
 *  ways to refer to them.
 */
enum Layout
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

//@}

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
 *
 * \param dimensions [in] an array of dimensions
 *
 * \param layout [in] the memory storage order
 */
template< class SIZE_TYPE, class DIM_TYPE >
Teuchos::Array< SIZE_TYPE >
computeStrides(const Teuchos::Array< DIM_TYPE > & dimensions,
               const Layout layout)
{
  int n = dimensions.size();
  Teuchos::Array< SIZE_TYPE > strides(n);
  if (n == 0) return strides;

  if (layout == FIRST_INDEX_FASTEST)
  {
    strides[0] = 1;
    for (int axis = 1; axis < n; ++axis)
      strides[axis] = strides[axis-1] * dimensions[axis-1];
  }
  else
  {
    strides[n-1] = 1;
    for (int axis = n-2; axis >= 0; --axis)
      strides[axis] = strides[axis+1] * dimensions[axis+1];
  }
  return strides;
}

////////////////////////////////////////////////////////////////////////

/** \brief Compute the strides of an <tt>MDArray</tt>,
 *         <tt>MDArrayView</tt>, or <tt>MDArrayRCP</tt>, given its
 *         dimensions as an ArrayView.
 *
 * \param dimensions [in] an array of dimensions
 *
 * \param layout [in] the memory storage order
 */
template< class SIZE_TYPE, class DIM_TYPE >
Teuchos::Array< SIZE_TYPE >
computeStrides(const Teuchos::ArrayView< DIM_TYPE > & dimensions,
               const Layout layout)
{
  // In the MDArray<T>(const MDArrayView<T> &) constructor, I try to
  // pass the MDArrayView dimensions to computeStrides(), but they
  // come in as ArrayView<const T> (for reasons I can't determine) and
  // cause all sorts of const-correctness problems.  So I copy them
  // into a new Array<T> and pass its reference to the main
  // computeStrides() function.  Fortunately, the array of dimensions
  // is small.
  Teuchos::Array< DIM_TYPE > nonConstDims(0);
  nonConstDims.insert(nonConstDims.begin(),
                      dimensions.begin(),
                      dimensions.end());
  return computeStrides< SIZE_TYPE, DIM_TYPE >(nonConstDims, layout);
}

////////////////////////////////////////////////////////////////////////

/** \brief Compute the minimum size required for an <tt>MDArray</tt>,
 *         <tt>MDArrayView</tt>, or <tt>MDArrayRCP</tt>, given its
 *         dimensions as an Arrayview.
 *
 * \param dimensions [in] an array of dimensions
 */
template< class DIM_TYPE >
size_type computeSize(const Teuchos::ArrayView< DIM_TYPE > & dimensions)
{
  size_type result = 1;
  for (int axis = 0; axis < dimensions.size(); ++axis)
    result *= dimensions[axis];
  return result;
}

////////////////////////////////////////////////////////////////////////

/** \brief Compute the minimum size required for an <tt>MDArray</tt>,
 *         <tt>MDArrayView</tt>, or <tt>MDArrayRCP</tt>, given its
 *         dimensions as an Array.
 *
 * \param dimensions [in] an array of dimensions
 */
template< class DIM_TYPE >
size_type computeSize(const Teuchos::Array< DIM_TYPE > & dimensions)
{
  // In the MDArray<T>(const MDArrayView<T> &) constructor, I try to
  // pass the MDArrayView dimensions to computeSize(), but they come
  // in as ArrayView<const T> (for reasons I can't determine) and
  // cause all sorts of const-correctness problems.  So I copy them
  // into a new Array<T> and pass its view to the main computeSize()
  // function.  Fortunately, the array of dimensions is small.
  Teuchos::Array< DIM_TYPE > nonConstDims(0);
  nonConstDims.insert(nonConstDims.begin(),
                      dimensions.begin(),
                      dimensions.end());
  return computeSize(nonConstDims());
}

////////////////////////////////////////////////////////////////////////

/** \brief Compute the minimum size required for an <tt>MDArray</tt>,
 *         <tt>MDArrayView</tt>, or <tt>MDArrayRCP</tt>, given its
 *         dimensions and strides.
 *
 * \param dimensions [in] an array of dimensions
 *
 * \param strides [in] an array of strides
 */
template< class SIZE_TYPE, class DIM_TYPE >
SIZE_TYPE computeSize(const Teuchos::ArrayView< DIM_TYPE > & dimensions,
                      const Teuchos::ArrayView< SIZE_TYPE > & strides)
{
  // SIZE_TYPE might be a const type, but we need result to be non-const
  typename remove_const< SIZE_TYPE >::type result = 1;
  for (int axis = 0; axis < dimensions.size(); ++axis)
    result += (dimensions[axis]-1) * strides[axis];
  return result;
}

////////////////////////////////////////////////////////////////////////

/** \brief Return and array of integers that represent the prime
 *         factors of the input argument
 *
 * \param n [in] integer to be factored
 */
Teuchos::Array< int >
factor(int n);

////////////////////////////////////////////////////////////////////////

/** \brief Return the index of the maximum value in the provided
 *         sequence of floating point values
 *
 * \param seq [in] a sequence of floating point values
 */
int indexOfMax(const Teuchos::ArrayView< const float > & seq);

////////////////////////////////////////////////////////////////////////

/** \brief Return a decomposition of processors, given on the total
 *         number of processors and the dimensions of the field that
 *         will be decomposed. The size of the dimensions input
 *         argument determines the size of the output
 *         decomposition. The product of the output decomposition will
 *         equal the total number of processors. The decomposition
 *         values will be as close to the relative sizes of the
 *         dimensions as the algorithm can get.
 *
 * \param nprocs [in] total number of processors
 *
 * \param dimensions [in] an array of dimensions of the field that
 *                        will be decomposed
 */
Teuchos::Array< int >
decomposeProcs(int nprocs,
               const Teuchos::ArrayView< dim_type > & dimensions);

////////////////////////////////////////////////////////////////////////

/** \brief Compute a valid commDims array, given the number of
 *         processors, the number of dimensions, and a candidate
 *         commDims array.
 *
 * \param numProcs [in] the number of processors
 *
 * \param numDims [in] the number of dimensions
 *
 * \param commDims [in] the communicor dimensions to be regularized
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
regularizeCommDims(int numProcs,
                   int numDims,
                   const Teuchos::ArrayView< const int > & commDims);

////////////////////////////////////////////////////////////////////////

/** \brief Given a Domi ParameterList, replace the "comm dimensions"
 *         parameter with a valid commDims array, given the number of
 *         processors and data from the ParameterList: the number of
 *         dimensions, the provided commDim array, and the dimensions
 *         of the field to be decomposed.
 *
 * \param numDims [in] the number of dimensions
 *
 * \param plist [in/out] ParameterList with construction information
 *        \htmlonly
 *        <iframe src="domi.xml" width="100%" scrolling="no" frameborder="0">
 *        </iframe>
 *        <hr />
 *        \endhtmlonly
 *
 * This function assumes that the ParameterList has already been
 * validated prior to being called.
 */
Teuchos::Array< int >
regularizeCommDims(int numProcs,
                   Teuchos::ParameterList & plist);

////////////////////////////////////////////////////////////////////////

/** \brief Compute the axis ranks for a given processor rank, given
 *         the communicator stride sizes along each axis.
 *
 * \param rank [in] the rank of the given processor
 *
 * \param commStrides [in] an array of the communicator stride sizes
 *        along each axis
 */
Teuchos::Array< int >
computeCommIndexes(int rank,
                   const Teuchos::ArrayView< int > & commStrides);

////////////////////////////////////////////////////////////////////////

/** \brief Compute an array of integers, given the user input.
 *
 * If the input array is shorter than the number of dimensions, fill
 * the remainder of the result with zeros.
 *
 * \param numDims [in] the number of dimensions
 *
 * \param source [in] the source array of ints
 */
Teuchos::Array< int >
createArrayOfInts(int numDims,
                  const Teuchos::ArrayView< const int > & source);

////////////////////////////////////////////////////////////////////////

/** \brief Given a std::string which contains comma-separated integers,
 *         return an array of ints.
 *
 * \param string [in] a string of comma-separated integers
 */
Teuchos::Array< int >
splitStringOfIntsWithCommas(std::string data);

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_MPI

/** \brief Return an MPI_Datatype, given its corresponding C/C++ type,
 *         specified as template parameter T
 */
template< class T > MPI_Datatype mpiType();

#endif

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_MPI

/** \brief Return an MPI flag for data layout, given its corresponding
 *         Domi enumeration
 *
 * \param layout [in] the Domi memory storage order
 */
int mpiOrder(Layout layout);

#endif

}       // End Domi namespace

#endif	// DOMI_MDARRAY_UTILS_HPP

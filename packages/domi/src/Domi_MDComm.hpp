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

#ifndef DOMI_MDCOMM_HPP
#define DOMI_MDCOMM_HPP

// Teuchos includes
#include "Teuchos_Comm.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"

// Domi includes
#include "Domi_ConfigDefs.hpp"
#include "Domi_Slice.hpp"

namespace Domi
{

/** \brief Provide a simple typedef for the Teuchos communicator
 */
typedef Teuchos::RCP< const Teuchos::Comm<int> > TeuchosComm;

/** \brief Multi-dimensional communicator object
 *
 * The <tt>MDComm</tt> is a relatively simple class that contains a
 * <tt>Teuchos::RCP<const Teuchos::Comm<int> ></tt> and arranges the
 * processor ranks of that communicator in a multi-dimensional grid
 * suitable for decomposing a multi-dimensional data structure.
 *
 * To construct an <tt>MDComm</tt>, the user provides the number of
 * dimensions, an array of sizes along each axis, or both.  The user
 * can request the constructor to compute some of the decomposition by
 * making the axis sizes list shorter than the number of dimensions,
 * or by providing axis sizes that are negative.  Currently, the
 * algorithm assigns all remaining processors to the first unspecified
 * axis size, and assigns all the rest to be one.
 *
 * An <tt>MDComm</tt> can also be constructed that is a slice of a
 * parent <tt>MDComm</tt>.  Such an object is constructed by providing
 * an RCP of the parent and an array of Slices that define the
 * sub-communicator.
 *
 * The <tt>MDComm</tt> class provides two methods that return the
 * processor rank of neighboring processors in the grid of arrays:
 *
 *     \code
 *     int getLowerNeighbor(int axis) const;
 *     int getUpperNeighbor(int axis) const;
 *     \code
 *
 * This should be sufficient to provide the necessary processor ranks
 * for conductiong halo updates.
 */
class MDComm
{
public:

  /** \name Teuchos::Array typedefs */
  //@{

  typedef Teuchos::Array< int >::size_type size_type;

  //@}

  /** \name Constructors and Destructor */
  //@{

  /** \brief Constructor with axis sizes
   *
   * \param teuchosComm [in] The Teuchos Communicator
   *
   * \param axisSizes [in] An array containing the sizes of the
   *        <tt>MDComm</tt> along each axis.  The size of
   *        <tt>axisSizes</tt> determines the number of dimensions.
   *        Negative values will be converted to positive such that
   *        the product of the resulting axis sizes will equal the
   *        number of processors in the Teuchos communicator.
   */
  MDComm(const TeuchosComm teuchosComm,
         const Teuchos::ArrayView< int > & axisSizes);

  /** \brief Constructor with number of dimensions
   *
   * \param teuchosComm [in] The Teuchos Communicator
   *
   * \param numDims [in] The number of dimensions in the
   *        <tt>MDComm</tt>.  Currently, all of the processors are
   *        allocated to the first axis, and the other axes are assign
   *        a size of one.
   */
  MDComm(const TeuchosComm teuchosComm,
         int numDims);

  /** \brief Constructor with number of dimensions and axis sizes
   *
   * \param teuchosComm [in] The Teuchos Communicator
   *
   * \param numDims [in] The number of dimensions in the
   *        <tt>MDComm</tt>.
   *
   * \param axisSizes [in] An array containing the sizes of the
   *        <tt>MDComm</tt> along each axis.  If the size of
   *        <tt>axisSizes</tt> is less than <tt>numDims</tt>, then the
   *        missing values are treated as unspecified.  Negative
   *        values will also be treated as unspecified.  Unspecified
   *        vlaues will be converted to positive such that the product
   *        of the resulting axis sizes will equal the number of
   *        processors in the Teuchos communicator.
   */
  MDComm(const TeuchosComm teuchosComm,
         int numDims,
         const Teuchos::ArrayView< int > & axisSizes);

  /** \brief Sub-communicator constructor
   *
   * \param parent [in] A <tt>Teuchos::RCP</tt> of the parent
   *        <tt>MDComm</tt> object.
   *
   * \param slices [in] An array of <tt>Slice</tt> objects that
   *        defines what portions of the parent will be translated to
   *        the sub-communicator along each axis.
   */
  MDComm(const Teuchos::RCP< MDComm > parent,
         const Teuchos::ArrayView< Slice > & slices);

  /** Destructor
   */
  ~MDComm();

  //@}

  /** \name Accessor methods */
  //@{

  /** \brief Get the Teuchos communicator
   */
  TeuchosComm getTeuchosComm() const;

  /** \brief Get the number of dimensions
   */
  int getNumDims() const;

  /** \brief Get the size along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   */
  int getAxisSize(int axis) const;

  /** \brief Get the axis-rank along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   */
  int getAxisRank(int axis) const;

  /** \brief Get the rank of the lower neighbor
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   */
  int getLowerNeighbor(int axis) const;

  /** \brief Get the rank of the upper neighbor
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   */
  int getUpperNeighbor(int axis) const;

  //@}

protected:

  // Not implemented
  MDComm();

private:

  // The Teuchos communicator
  TeuchosComm           _teuchosComm;

  // An array of the sizes along each axis
  Teuchos::Array< int > _axisSizes;

  // An array of the axis-ranks for this processor along each axis
  Teuchos::Array< int > _axisRanks;

  // An array of the strides between processor ranks along each axis.
  Teuchos::Array< int > _axisStrides;

};

}  // namespace Domi

#endif

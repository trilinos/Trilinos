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

#ifndef DOMI_MDCOMM_HPP
#define DOMI_MDCOMM_HPP

// Teuchos includes
#include "Teuchos_Comm.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_ParameterList.hpp"

// Domi includes
#include "Domi_ConfigDefs.hpp"
#include "Domi_Slice.hpp"

#ifdef HAVE_EPETRA
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#endif

namespace Domi
{

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
 * sub-communicator.  When you construct a sub-communicator, you run
 * the risk of calling methods of the <tt>MDComm</tt> object or of
 * objects built on top of the <tt>MDComm</tt> object from processors
 * that are not a part of the sub-communicator.  You can always check
 * whether this is safe with the
 *
 *     \code
 *     bool onSubcommunicator() const;
 *     \endcode
 *
 * method.
 *
 * The <tt>MDComm</tt> class provides two methods that return the
 * processor rank of neighboring processors in the grid of arrays:
 *
 *     \code
 *     int getLowerNeighbor(int axis) const;
 *     int getUpperNeighbor(int axis) const;
 *     \endcode
 *
 * This should be sufficient to provide the necessary processor ranks
 * for conducting communication padding updates.
 */
class MDComm
{
public:

  /** \name Teuchos::Array typedefs */
  //@{

  typedef Teuchos::Array< int >::size_type size_type;

  //@}

  /** \name MDComm layout */
  //@{

  static const Layout commLayout;

  //@}

  /** \name Constructors and Destructor */
  //@{

  /** \brief Constructor with default Teuchos comm and axis sizes
   *
   * \param commDims [in] An array containing the sizes of the
   *        <tt>MDComm</tt> along each axis.  The size of
   *        <tt>commDims</tt> determines the number of dimensions.
   *        Negative values will be converted to positive such that
   *        the product of the resulting axis sizes will equal the
   *        number of processors in the Teuchos communicator.
   *
   * \param periodic [in] An array of ints which are simple flags
   *        denoting whether each axis is periodic.  If this array is
   *        shorter than the size of commDims, the unspecified axes
   *        are assumed to be zero (false).
   *
   * This constructor uses the Teuchos::DefaultComm communicator
   */
  MDComm(const Teuchos::ArrayView< int > & commDims,
         const Teuchos::ArrayView< int > & periodic =
           Teuchos::ArrayView< int >());

  /** \brief Constructor with Teuchos Comm and axis sizes
   *
   * \param teuchosComm [in] The Teuchos Communicator
   *
   * \param commDims [in] An array containing the sizes of the
   *        <tt>MDComm</tt> along each axis.  The size of
   *        <tt>commDims</tt> determines the number of dimensions.
   *        Negative values will be converted to positive such that
   *        the product of the resulting axis sizes will equal the
   *        number of processors in the Teuchos communicator.
   *
   * \param periodic [in] An array of ints which are simple flags
   *        denoting whether each axis is periodic.  If this array is
   *        shorter than the size of commDims, the unspecified axes
   *        are assumed to be zero (false).
   */
  MDComm(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
         const Teuchos::ArrayView< int > & commDims,
         const Teuchos::ArrayView< int > & periodic =
           Teuchos::ArrayView< int >());

  /** \brief Constructor with ParameterList
   * 
   * \param plist [in] ParameterList with construction information
   *        \htmlonly
   *        <iframe src="domi.xml" width="100%" scrolling="no" frameborder="0">
   *        </iframe>
   *        <hr />
   *        \endhtmlonly
   *
   * This constructor uses the Teuchos::DefaultComm
   */
  MDComm(Teuchos::ParameterList & plist);

  /** \brief Constructor with Teuchos Comm and ParameterList
   * 
   * \param teuchosComm [in] The Teuchos Communicator
   *
   * \param plist [in] ParameterList with construction information
   *        \htmlonly
   *        <iframe src="domi.xml" width="100%" scrolling="no" frameborder="0">
   *        </iframe>
   *        <hr />
   *        \endhtmlonly
   */
  MDComm(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
         Teuchos::ParameterList & plist);

  /** \brief Constructor with number of dimensions
   *
   * \param numDims [in] The number of dimensions in the
   *        <tt>MDComm</tt>.  Currently, all of the processors are
   *        allocated to the first axis, and the other axes are
   *        assigned a size of one.
   *
   * This constructor uses the Teuchos::DefaultComm
   */
  MDComm(int numDims);

  /** \brief Constructor with number of dimensions
   *
   * \param teuchosComm [in] The Teuchos Communicator
   *
   * \param numDims [in] The number of dimensions in the
   *        <tt>MDComm</tt>.  Currently, all of the processors are
   *        allocated to the first axis, and the other axes are
   *        assigned a size of one.
   */
  MDComm(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
         int numDims);

  /** \brief Constructor with number of dimensions and axis sizes
   *
   * \param numDims [in] The number of dimensions in the
   *        <tt>MDComm</tt>.
   *
   * \param commDims [in] An array containing the sizes of the
   *        <tt>MDComm</tt> along each axis.  If the size of
   *        <tt>commDims</tt> is less than <tt>numDims</tt>, then the
   *        missing values are treated as unspecified.  Negative
   *        values will also be treated as unspecified.  Unspecified
   *        values will be converted to positive such that the product
   *        of the resulting axis sizes will equal the number of
   *        processors in the Teuchos communicator.
   *
   * \param periodic [in] An array of ints which are simple flags
   *        denoting whether each axis is periodic.  If this array is
   *        shorter than numDims, the unspecified axes are assumed to
   *        be zero (false).
   *
   * This constructor uses the Teuchos::DefaultComm
   */
  MDComm(int numDims,
         const Teuchos::ArrayView< int > & commDims,
         const Teuchos::ArrayView< int > & periodic =
           Teuchos::ArrayView< int >());

  /** \brief Constructor with Teuchos Comm, number of dimensions and axis sizes
   *
   * \param teuchosComm [in] The Teuchos Communicator
   *
   * \param numDims [in] The number of dimensions in the
   *        <tt>MDComm</tt>.
   *
   * \param commDims [in] An array containing the sizes of the
   *        <tt>MDComm</tt> along each axis.  If the size of
   *        <tt>commDims</tt> is less than <tt>numDims</tt>, then the
   *        missing values are treated as unspecified.  Negative
   *        values will also be treated as unspecified.  Unspecified
   *        values will be converted to positive such that the product
   *        of the resulting axis sizes will equal the number of
   *        processors in the Teuchos communicator.
   *
   * \param periodic [in] An array of ints which are simple flags
   *        denoting whether each axis is periodic.  If this array is
   *        shorter than numDims, the unspecified axes are assumed to
   *        be zero (false).
   */
  MDComm(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
         int numDims,
         const Teuchos::ArrayView< int > & commDims,
         const Teuchos::ArrayView< int > & periodic =
           Teuchos::ArrayView< int >());

  /** \brief Axis rank sub-communicator constructor
   *
   * \param parent [in] The parent <tt>MDComm</tt> object.
   *
   * \param axis [in] The axis to which the axisRank argument applies
   *
   * \param axisRank [in] The value of the rank along the given axis
   *
   * This constructor will return an MDComm that is one dimension less
   * than the parent MDComm (unless the parent MDComm is already
   * one-dimensional), equivalent to a slice of the parent at the
   * given axis rank value.
   */
   MDComm(const MDComm & parent,
         int axis,
         int axisRank);

  /** \brief Slice sub-communicator constructor
   *
   * \param parent [in] The parent <tt>MDComm</tt> object.
   *
   * \param axis [in] The axis to which the slice argument applies
   *
   * \param slice [in] A <tt>Slice</tt> object that defines what
   *        portion of the parent will be translated to the
   *        sub-communicator along the given axis axis.
   */
  MDComm(const MDComm & parent,
         int axis,
         const Slice & slice);

  /** \brief Array of Slices sub-communicator constructor
   *
   * \param parent [in] The parent <tt>MDComm</tt> object.
   *
   * \param slices [in] An array of <tt>Slice</tt> objects that
   *        defines what portions of the parent will be translated to
   *        the sub-communicator along each axis.
   */
  MDComm(const MDComm & parent,
         const Teuchos::ArrayView< Slice > & slices);

  /** \brief Copy constructor
   *
   * \param source [in] Source MDComm to be copied
   */
  MDComm(const MDComm & source);

  /** Destructor
   */
  ~MDComm();

  /** \brief Assignment operator
   *
   * \param source [in] MDComm to be copied
   */
  MDComm & operator=(const MDComm & source);

  //@}

  /** \name Accessor methods */
  //@{

  /** \brief Query whether this processor is on the sub-communicator.
   *
   * Sub-communicators are formed when a parent MDComm is sliced by
   * using the (parent,slices) constructor.  For a full communicator,
   * this method will always return true.
   */
  bool onSubcommunicator() const;

  /** \brief Get the Teuchos communicator
   *
   * Note that if the communicator is not a full communicator, i.e. a
   * sub-communicator, that the underlying Comm pointer may be NULL,
   * depending on this processor's rank.
   */
  Teuchos::RCP< const Teuchos::Comm< int > > getTeuchosComm() const;

#ifdef HAVE_EPETRA
  /** \brief Get an equivalent Epetra communicator
   *
   * Note that if the communicator is not a full communicator, i.e. a
   * sub-communicator, that the underlying Comm pointer may be NULL,
   * depending on this processor's rank.
   */
  Teuchos::RCP< const Epetra_Comm > getEpetraComm() const;
#endif

  /** \brief Get the number of dimensions
   *
   * This method will return 0 if the communicator is a
   * sub-communicator and this processor does not belong to the
   * sub-communicator.
   */
  int numDims() const;

  /** \brief Get the communicator size along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This method will throw a Domi::SubcommunicatorError if the
   * communicator is a sub-communicator and this processor does not
   * belong to the sub-communicator.
   */
  int getCommDim(int axis) const;

  /** \brief Return the periodic flag for the given axis.
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This method will throw a Domi::SubcommunicatorError if the
   * communicator is a sub-communicator and this processor does not
   * belong to the sub-communicator.
   */
  bool isPeriodic(int axis) const;

  /** \brief Get the comm index along the given axis
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This method will throw a Domi::SubcommunicatorError if the
   * communicator is a sub-communicator and this processor does not
   * belong to the sub-communicator.
   */
  int getCommIndex(int axis) const;

  /** \brief Get the rank of the lower neighbor
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This method will throw a Domi::SubcommunicatorError if the
   * communicator is a sub-communicator and this processor does not
   * belong to the sub-communicator.
   *
   * If the periodic flag for the given axis is off, and the axis rank
   * of the calling processor is zero, then this method returns -1.
   *
   * If the periodic flag for the given axis is on, and the axis rank
   * of the calling processor is zero, then the returned lower
   * neighbor will be the highest axis rank processor along this axis.
   */
  int getLowerNeighbor(int axis) const;

  /** \brief Get the rank of the upper neighbor
   *
   * \param axis [in] the index of the axis (from zero to the number
   *        of dimensions - 1)
   *
   * This method will throw a Domi::SubcommunicatorError if the
   * communicator is a sub-communicator and this processor does not
   * belong to the sub-communicator.
   *
   * If the periodic flag for the given axis is off, and the axis rank
   * of the calling processor is the highest axis rank along the axis,
   * then this method returns -1.
   *
   * If the periodic flag for the given axis is on, and the axis rank
   * of the calling processor is the highest axis rank along the axis,
   * then this method will return axis rank zero.
   */
  int getUpperNeighbor(int axis) const;

  //@}

protected:

  // Not implemented
  MDComm();

private:

  // The Teuchos communicator
  Teuchos::RCP< const Teuchos::Comm< int > > _teuchosComm;

#ifdef HAVE_EPETRA
  // An equivalent Epetra communicator.  This is mutable because we
  // only compute it if requested by a get method that is const.
  mutable Teuchos::RCP< const Epetra_Comm > _epetraComm;
#endif

  // An array of the sizes of the communicator along each axis
  Teuchos::Array< int > _commDims;

  // An array of the strides between processor ranks along each axis.
  Teuchos::Array< int > _commStrides;

  // The comm index for this processor along each axis
  Teuchos::Array< int > _commIndex;

  // An array of flags denoting periodic axes
  Teuchos::Array< int > _periodic;

};

}  // namespace Domi

#endif

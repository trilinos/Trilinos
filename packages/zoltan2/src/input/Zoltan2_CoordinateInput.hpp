// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_CoordinateInput.hpp
    \brief Defines the CoordinateInput interface.     
*/

#ifndef _ZOLTAN2_COORDINATEINPUT_HPP_
#define _ZOLTAN2_COORDINATEINPUT_HPP_

#include <Zoltan2_InputAdapter.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

#include <string>

namespace Zoltan2 {

/*!  \brief CoordinateInput defines the interface for input geometric
                coordinates.

    InputAdapter objects provide access for Zoltan2 to the user's data.
    Many built-in adapters are already defined for common data structures,
    such as Tpetra and Epetra objects and C-language pointers to arrays.

    Data types:
    \li \c scalar_t coordinate values and weights
    \li \c lno_t    local indices and local counts
    \li \c gno_t    global indices and global counts
    \li \c gid_t    application global Ids
    \li \c node_t is a sub class of Kokkos::StandardNodeMemoryModel

    See IdentifierTraits to understand why the user's global ID type (\c gid_t)
    may differ from that used by Zoltan2 (\c gno_t).

    The Kokkos node type can be safely ignored.

    The template parameter \c User is a user-defined data type
    which, through a traits mechanism, provides the actual data types
    with which the Zoltan2 library will be compiled.
    \c User may be the actual class or structure used by application to
    represent a vector, or it may be the helper class BasicUserTypes.
    See InputTraits for more information.


    The \c scalar_t type, representing use data such as matrix values, is
    used by Zoltan2 for weights, coordinates, and quality metrics.
    Some User types (like Tpetra::CrsMatrix) have an inherent scalar type,
    and some
    (like Tpetra::CrsGraph) do not.  For such objects, the scalar type is
    set by Zoltan2 to \c float.  If you wish to change it to double, set
    the second template parameter to \c double.


    Input adapters provide access for Zoltan2 to the user's data.  The
    methods in the interface must be defined by users.  Many built-in
    adapters are already defined for common data structures, such as
    Tpetra and Epetra objects and C-language pointers to arrays.

    \todo We don't really need global Ids.  They should be optional
    \todo Migration doesn't move weights.  Should it?
*/

template <typename User>
  class CoordinateInput : public InputAdapter<User> {

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef User user_t;
#endif

  /*! \brief Destructor
   */
  virtual ~CoordinateInput() {};

  ////////////////////////////////////////////////////
  // The InputAdapter interface.
  ////////////////////////////////////////////////////

  enum InputAdapterType inputAdapterType() const {return CoordinateAdapterType;}

  ////////////////////////////////////////////////////
  // My interface.
  ////////////////////////////////////////////////////

  /*! \brief Return dimension of the coordinates.
   */
  virtual int getCoordinateDimension() const = 0;

  /*! \brief Return the number of weights per coordinate.
   *   \return the count of weights, zero or more per coordinate.
   */
  virtual int getNumberOfWeights() const = 0;

  /*! \brief Return the number of coordinates on this process.
   *   \return  the count of coordinates on the local process.
   */
  virtual size_t getLocalNumberOfCoordinates() const = 0;

  /*! \brief Provide a pointer to one dimension of this process' coordinates.
      \param coordDim  is a value from 0 to one less than 
         getLocalNumberOfCoordinates() specifying which dimension is
         being provided in the coords list.
      \param coords  points to a list of coordinate values for the dimension.
      \param stride  describes the layout of the coordinate values in
              the coords list.  If stride is one, then the ith coordinate
              value is coords[i], but if stride is two, then the
              ith coordinate value is coords[2*i].

       \return The length of the \c coords list.  This may be more than
              getLocalNumberOfCoordinates() because the \c stride
              may be more than one.

      Zoltan2 does not copy your data.  The data pointed to coords
      must remain valid for the lifetime of this InputAdapter.
   */

  virtual size_t getCoordinates(int coordDim, const gid_t *&gids, 
    const scalar_t *&coords, int &stride) const = 0;

  /*! \brief  Provide a pointer to the weights, if any, corresponding 
       to the coordinates returned in getCoordinates(). 

      \param weightDim ranges from zero to one less than getNumberOfWeights()
      \param weights is the list of weights of the given dimension for
           the coordinates listed in getCoordinates().  If weights for
           this dimension are to be uniform for all coordinates in the
           global problem, the \c weights should be a NULL pointer.
       \param stride The k'th weight is located at weights[stride*k]
       \return The number of weights listed, which should be at least
                  the local number of coordinates times the stride for
                  non-uniform weights, zero otherwise.
   */

  virtual size_t getCoordinateWeights(int weightDim,
     const scalar_t *&weights, int &stride) const = 0;

  /*! \brief Apply a PartitioningSolution to an input.
   *
   *  This is not a required part of the CoordinateInput interface. However
   *  if the Caller calls a Problem method to redistribute data, it needs
   *  this method to perform the redistribution.
   *
   *  \param in  An input object with a structure and assignment of
   *           of global Ids to processes that matches that of the input
   *           data that instantiated this InputAdapter.
   *  \param out On return this should point to a newly created object
   *            with the specified partitioning.
   *  \param solution  The Solution object created by a Problem should
   *      be supplied as the third argument.  It must have been templated
   *      on user data that has the same global ID distribution as this
   *      user data.
   *  \return   Returns the number of local Ids in the new partitioning.
   */

  template <typename Adapter>
    size_t applyPartitioningSolution(User &in, User *&out,
         const PartitioningSolution<Adapter> &solution) const
  {
    return 0;
  } 

private:
};
  
  
}  //namespace Zoltan2
  
#endif

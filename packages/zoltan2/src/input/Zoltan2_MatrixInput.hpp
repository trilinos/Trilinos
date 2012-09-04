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

/*! \file Zoltan2_MatrixInput.hpp
    \brief Defines the MatrixInput adapter interface.
*/

#ifndef _ZOLTAN2_MATRIXINPUT_HPP_
#define _ZOLTAN2_MATRIXINPUT_HPP_

#include <Zoltan2_InputAdapter.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

namespace Zoltan2 {

/*!  \brief MatrixInput defines the interface for input adapters for matrices.

    InputAdapter objects provide access for Zoltan2 to the user's data.
    Many built-in adapters are already defined for common data structures,
    such as Tpetra and Epetra objects and C-language pointers to arrays.

    Data types:
    \li \c scalar_t row, column or non-zero weights
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
    used by Zoltan2 for weights, coordinates, part sizes and
    quality metrics.
    Some User types (like Tpetra::CrsMatrix) have an inherent scalar type,
    and some
    (like Tpetra::CrsGraph) do not.  For such objects, the scalar type is
    set by Zoltan2 to \c float.  If you wish to change it to double, set
    the second template parameter to \c double.

     \todo Create BasicCrsMatrixInput subclass
     \todo Do we want to require input adapters to give us the global
               number of rows, columns etc?  We can figure that out.
      \todo  This is a row-oriented matrix.  Do we need a column-oriented
              matrix?  In particular - we assumed coordinates are for rows.
      \todo  If the user can tell us there are no diagonal entries
        in a square matrix, it can save time if we have to remove
        them for the algorithm.  Should we have a set method in 
        subclasses for setMatrixHasDiagonalEntries yes, no and maybe?
*/

template <typename User>
  class MatrixInput : public InputAdapter<User> {
private:

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef User user_t;
#endif

  enum InputAdapterType inputAdapterType() const {return MatrixAdapterType;}

  /*! \brief Destructor
   */
  virtual ~MatrixInput(){};

  /*! \brief Returns the number rows on this process.
   */
  virtual size_t getLocalNumRows() const = 0;

  /*! \brief Returns the number columns on this process.
   */
  virtual size_t getLocalNumColumns() const = 0;

  /*! \brief Sets pointers to this process' matrix entries.
      If this optional call is defined in the adapter, it can save a memory
      copy of application data.
      \param rowIds will on return a pointer to row global Ids
      \param offsets is an array of size numRows + 1.  The column Ids for
          rowId[i] begin at colIds[offsets[i]].  The last element of offsets
          is the size of the colIds array.
      \param colIds on return will point to the global column Ids for
         the non-zeros for each row.
       \return The number of ids in the rowIds list.
   */

  virtual size_t getRowListView(const gid_t *&rowIds, 
    const lno_t *&offsets, const gid_t *& colIds) const = 0;

  /*! \brief Returns the dimension (0 or greater) of row weights.

      Row weights may be used when partitioning matrix rows.
      If no weights are provided, and rowWeightIsNumberOfNonZeros()
      is \c false, then it is assumed rows are equally weighted.
   */
  virtual int getRowWeightDimension() const = 0;

  /*! \brief  Provide a pointer to the row weights, if any.

      \param weightDim ranges from zero to one less than
                   getRowWeightDimension().
      \param weights is the list of weights of the given dimension for
           the rows returned in getRowListView().  If weights for
           this dimension are to be uniform for all rows in the
           global problem, the \c weights should be a NULL pointer.
       \param stride The k'th weight is located at weights[stride*k]
      \return The number of weights listed, which should be at least
                  the local number of rows times the stride for
                  non-uniform weights, zero otherwise.

      Zoltan2 does not copy your data.  The data pointed to by weights
      must remain valid for the lifetime of this InputAdapter.
   */

  virtual size_t getRowWeights(int weightDim,
     const scalar_t *&weights, int &stride) const = 0;

  /*! \brief Is the row weight for a dimension the number of row non-zeros?
   *   \param dim a value between zero and getRowWeightDimension() that
   *     specifies the weight dimension that will be equal to the
   *      number of row non-zeros.
   *
   *  If you wish for Zoltan2 to automatically assign for one of
   *  the weight dimensions for each row 
   *  a weight that is equal to the number of non-zeros in the row,
   *  then return true here.  Otherwise return false.
   */

  virtual bool getRowWeightIsNumberOfNonZeros(int dim) const = 0;

  /*! \brief Returns the dimension of the geometry, if any.
   *
   *  Some algorithms can use geometric row or column coordinate
   *    information if it is present.  Given the problem parameters
   *    supplied by the user, it may make sense to use row coordinates
   *    or it may make sense to use column coordinates.
   */

  virtual int getCoordinateDimension() const = 0;

  /*! \brief Provide a pointer to one dimension of row coordinates.
      \param coordDim [input] is a value from 0 to one less than
         getCoordinateDimension() specifying which dimension is
         being provided in the coords list.
      \param coords  [output] points to a list of coordinate values for the dimension.
             The order of \c coords should correspond to the order of \c rowIds
             in getRowListView().
      \param stride  [output] describes the layout of the coordinate values in
              the coords list.  If stride is one, then the ith coordinate
              value is coords[i], but if stride is two, then the
              ith coordinate value is coords[2*i].

       \return The length of the \c coords list.  This may be more than
              getLocalNumberOfVertices() because the \c stride
              may be more than one.  It may also be zero if column 
              coordinates but not row coordinates are supplied.

      Zoltan2 does not copy your data.  The data pointed to by coords
      must remain valid for the lifetime of this InputAdapter.
   */

  virtual size_t getRowCoordinates(int coordDim,
    const scalar_t *&coords, int &stride) const = 0;

  /*! \brief Apply the solution to a partitioning problem to an input.  
   *
   *  This is not a required part of the MatrixInput interface.  
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
   *  \return   Returns the number of Ids in the new partitioning.
   */

  template <typename Adapter>
    size_t applyPartitioningSolution(const User &in, User *&out,
         const PartitioningSolution<Adapter> &solution) const
  {
    return 0;
  }

};
  
}  //namespace Zoltan2
  
#endif

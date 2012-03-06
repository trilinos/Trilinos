// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_MatrixInput.hpp
    \brief Defines the MatrixInput adapter interface.
*/

#ifndef _ZOLTAN2_MATRIXINPUT_HPP_
#define _ZOLTAN2_MATRIXINPUT_HPP_

#include <Zoltan2_InputAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>
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

     \todo Do we want to require input adapters to give us the global
               number of rows, columns etc?  We can figure that out.
     \todo Do we want to add the ability for the user to supply row
             or column weights, or is this something the algorithm
             will add?
      \todo  This is a row-oriented matrix.  Do we need a column-oriented
              matrix?
*/

template <typename User>
  class MatrixInput : public InputAdapter {
private:

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t scalar_t;
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

  /*! \brief Returns the global number rows.
   */
  virtual global_size_t getGlobalNumRows() const = 0;

  /*! \brief Returns the number columns on this process.
   */
  virtual size_t getLocalNumColumns() const = 0;

  /*! \brief Returns the global number columns.
   */
  virtual global_size_t getGlobalNumColumns() const = 0;

  /*! \brief Return true if the sparse square matrix may globally have
   *  diagonal entries.  Return false otherwise.
   */
  virtual bool diagonalEntriesMayBePresent() const = 0;

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

  /*! \brief Returns the dimension of the geometry, if any.
   *
   *  Some algorithms can use geometric row or column coordinate
   *    information if it is present.  Given the problem parameters
   *    supplied by the user, it may make sense to use row coordinates
   *    or it may make sense to use column coordinates.
   */

  virtual int getCoordinateDimension() const = 0;

  /*! \brief Provide a pointer to one dimension of row coordinates.
      \param coordDim  is a value from 0 to one less than
         getCoordinateDimension() specifying which dimension is
         being provided in the coords list.
      \param coords  points to a list of coordinate values for the dimension.
             The order of \c coords should correspond to the order of \c rowIds
             in getRowListView().
      \param stride  describes the layout of the coordinate values in
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

  template <typename User2>
    size_t applyPartitioningSolution(const User &in, User *&out,
         const PartitioningSolution<User2> &solution) const
  {
    return 0;
  }

};
  
}  //namespace Zoltan2
  
#endif

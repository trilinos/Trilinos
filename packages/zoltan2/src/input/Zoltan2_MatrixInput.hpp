// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_MatrixInput.hpp

    \brief The abstract interface for a graph input adapter.
*/

#ifndef _ZOLTAN2_MATRIXINPUT_HPP_
#define _ZOLTAN2_MATRIXINPUT_HPP_

#include <Zoltan2_InputAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

namespace Zoltan2 {

/*! Zoltan2::MatrixInput
    \brief MatrixInput defines the interface for input adapters for
            matrices.

    The Matrix accessor methods defined here mimic those of 
    Tpetra::CrsMatrix and Tpetra::Map.

    scalar_t: This data type is used for matrix non-zeros
    gid_t: the type for the application's global Ids
    lno_t: the integral type that Zoltan2 will use for local counters.
    gno_t: the integral type that Zoltan2 will use for the global 
           counts and identifiers.  It needs to be large enough for the
           problem's number of objects.
*/

template <typename User>
  class MatrixInput : public InputAdapter {
private:

public:

  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef User user_t;

  // adapterType == MatrixAdapterType
  // Function must return one of Zoltan2's enumerated types in InputAdapter
  // User should not rewrite this function.
  enum InputAdapterType inputAdapterType() const {return MatrixAdapterType;}

  /*! Pure virtual Destructor
   */
  virtual ~MatrixInput(){};

  /*! Returns the number rows on this process.
   */
  virtual size_t getLocalNumRows() const = 0;

  /*! Returns the global number rows.
   */
  virtual global_size_t getGlobalNumRows() const = 0;

  /*! Returns the number columns on this process.
   */
  virtual size_t getLocalNumColumns() const = 0;

  /*! Returns the global number columns.
   */
  virtual global_size_t getGlobalNumColumns() const = 0;

  /*! Return true if the sparse square matrix may globally have
   *  diagonal entries.  Return false otherwise.
   */
  virtual bool diagonalEntriesMayBePresent() const = 0;

  /*! Sets pointers to this process' matrix entries.
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

  /*! Apply the solution to a partitioning problem to an input.  
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
         const PartitioningSolution<User2> &solution)
  {
    return 0;
  }

};
  
}  //namespace Zoltan2
  
#endif

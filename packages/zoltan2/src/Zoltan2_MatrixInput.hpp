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

#include <string>
#include <Zoltan2_InputAdapter.hpp>

namespace Zoltan2 {

/*! Zoltan2::MatrixInput
    \brief MatrixInput defines the interface for input adapters for
            matrices.

    The Matrix accessor methods defined here mimic those of 
    Tpetra::CrsMatrix and Tpetra::Map.

    scalar_t: This data type is used for matrix non-zeros
    lid_t: the type for the application's local Ids
    gid_t: the type for the application's global Ids
    lno_t: the integral type that Zoltan2 will use for local counters.
    gno_t: the integral type that Zoltan2 will use for the global 
           counts and identifiers.  It needs to be large enough for the
           problem's number of objects.
*/

template <typename User>
  class MatrixInput : public InputAdapter<User> {
private:

public:

  typedef typename InputAdapter<User>::scalar_t scalar_t;
  typedef typename InputAdapter<User>::lno_t    lno_t;
  typedef typename InputAdapter<User>::gno_t    gno_t;
  typedef typename InputAdapter<User>::lid_t    lid_t;
  typedef typename InputAdapter<User>::gid_t    gid_t;
  typedef typename InputAdapter<User>::node_t   node_t;

  // adapterType == MatrixAdapterType
  // Function must return one of Zoltan2's enumerated types in InputAdapter
  // User should not rewrite this function.
  enum InputAdapterType inputAdapterType() {return MatrixAdapterType;}

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

  /*! Returns list of this process' matrix entries.
      \param rowIds will on return a list of row global Ids
      \param localIds can, optionally, on return hold a list of locally
        relevant values that the process will use to refer to the objects
        listed in the first list.  If localIds are omitted and
        haveConsecutiveLocalIds is true, it is assumed that the
        global Ids are in local Id order.
      \param rowSize on return will list the number of non-zeros
        for each row.
      \param colIds on return will list the global column Ids for
         the non-zeros for each row.
   */

  virtual void getRowListCopy(std::vector<gid_t> &rowIds, 
    std::vector<lid_t> &localIds, std::vector<lno_t> &rowSize,
    std::vector<gid_t> &colIds) const = 0;

  /*! Sets pointers to this process' matrix entries.
      If this optional call is defined in the adapter, it can save a memory
      copy of application data.
      \param rowIds will on return a pointer to row global Ids
      \param localIds can, optionally, on return hold a list of locally
        relevant values that the process will use to refer to the objects
        listed in the first list.  If localIds are omitted and
        haveConsecutiveLocalIds is true, it is assumed that the
        global Ids are in local Id order.
      \param offsets is an array of size numRows + 1.  The column Ids for
          rowId[i] begin at colIds[offsets[i]].  The last element of offsets
          is the size of the colIds array.
      \param colIds on return will point to the global column Ids for
         the non-zeros for each row.
       \return The number of ids in the rowIds list.
   */

  virtual size_t getRowListView(const gid_t *&rowIds, const lid_t *&localIds, 
    const lno_t *&offsets, const gid_t *& colIds) const = 0;
};
  
  
}  //namespace Zoltan2
  
#endif

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

    Scalar: This data type is used for matrix non-zeros
    LID: the type for the application's local Ids
    GID: the type for the application's global Ids
    LNO: the integral type that Zoltan2 will use for local counters.
    GNO: the integral type that Zoltan2 will use for the global 
      counts and identifiers.  It needs to be large enough for the
      problem's number of objects.
*/

CONSISTENT_CLASS_TEMPLATE_LINE
  class MatrixInput : public InputAdapter {
private:

public:

  typedef Scalar scalarType;
  typedef LID lidType;
  typedef GID gidType;
  typedef LNO lnoType;
  typedef GNO gnoType;
  typedef Node  nodeType;

  // adapterType == MatrixAdapterType
  // Function must return one of Zoltan2's enumerated types in InputAdapter
  // User should not rewrite this function.
  enum InputAdapterType inputAdapterType() {return MatrixAdapterType;}

  // TODO - what to do with destructor
  //virtual ~MatrixInput();

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

  virtual void getRowListCopy(std::vector<GID> &rowIds, 
    std::vector<LID> &localIds, std::vector<LNO> &rowSize,
    std::vector<GID> &colIds) const = 0;

  /*! Sets pointers to this process' matrix entries.
      If this optional call is defined in the adapter, it can save a memory
      copy of application data.
      \param rowIds will on return a pointer to row global Ids
      \param localIds can, optionally, on return hold a list of locally
        relevant values that the process will use to refer to the objects
        listed in the first list.  If localIds are omitted and
        haveConsecutiveLocalIds is true, it is assumed that the
        global Ids are in local Id order.
      \param rowSize on return will point to the number of non-zeros
        for each row.
      \param colIds on return will point to the global column Ids for
         the non-zeros for each row.
       \return The number of ids in the rowIds list.
   */

  LNO getRowListView(GID *&rowIds, LID *&localIds, 
    LNO *&rowSize, GID *& colIds)
  {
    rowIds = NULL;
    localIds = NULL;
    rowSize = NULL;
    colIds = NULL;
    return 0;
  }
};
  
  
}  //namespace Zoltan2
  
#endif

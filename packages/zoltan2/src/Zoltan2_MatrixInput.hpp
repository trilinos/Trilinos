// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_MatrixInput.hpp

    \brief The abstract interface for a matrix input adapter.

    \author Siva Rajamanickam
*/


#ifndef _ZOLTAN2_MATRIXINPUT_HPP_
#define _ZOLTAN2_MATRIXINPUT_HPP_

#include <Zoltan2_TemplateMacros.hpp>
#include <Zoltan2_InputAdapter.hpp>

namespace Zoltan2 {

/*! Zoltan2::MatrixInput
    \brief The MatrixInput is the abstract base class for matrix input adapters.

    The Matrix accessor methods defined here mimic those of Tpetra::CrsMatrix
    These public methods define the graph adapter interface to Zoltan2 models.
    TODO: It may be necessary to put a migration interface at this level.
*/

CONSISTENT_CLASS_TEMPLATE_LINE
class MatrixInput : public InputAdapter{
private:

public:

  // adapterType == MatrixAdapterType
  // Function must return one of Zoltan2's enumerated types in InputAdapter
  // User should not rewrite this function.
  enum InputAdapterType adapterType() {return MatrixAdapterType;}

  /*! Destructor TODO don't know what to do about destructor */
  virtual ~MatrixInput(){};

  /*! Returns the number rows on this process.
   */
  virtual LNO getLocalNumRows() const = 0;

  /*! Returns true if input adapter uses local Ids.
   */
  virtual bool haveLocalIds() const = 0;

  /*! Return true if local Ids are consecutive integral
   *   values and supply the base.  Providing this information
   *   can save memory, making local Id lists unneccesary.
   */
  virtual bool haveConsecutiveLocalIds(LID &base) const = 0;

  /*! Returns the number columns used by rows on this process
   */
  virtual LNO getLocalNumColumns() const = 0;

  /*! Return the total number of non-zero entries on this process.
  */
  virtual LNO getLocalNumNonZeros() const = 0;

  /*! Return the maximum number of non-zero entries in any local row.
  */
  virtual LNO getLocalMaxNumNonZeros() const = 0;

  /*! Return the local row information
      \param Ids will on return hold a list of the global Ids for
        each row on this process.
      \param localIds can, optionally, on return hold a list of locally
        relevant values that the process will use to refer to the objects
        listed in the first list.  If localIds are omitted and
        haveConsecutiveLocalIds is true, it is assumed that the
        global Ids are in local Id order.
      \param nnz will on return hold the number of non zeros in the
         cooresponding row.
  */
  virtual void getRowListCopy(std::vector<GID> &Ids,
    std::vector<LID> &localIds, std::vector<LNO> &nnz);

  /*! Sets pointers to this process' row Ids and non-zero count.
      If this optional call is defined in the adapter, it can save a memory
      copy of application data.
      \param Ids will on return point to the list of the global Ids for
        each row on this process.
      \param localIds can, optionally, on return point to a list of locally
        relevant values that the process will use to refer to the objects
        listed in the first list. If localIds is NULL and
        haveConsecutiveLocalIds is true, it is assumed that the
        global Ids are in local ID order.
      \param nnz will on return point to a list of the number of non-zeros
         in the corresponding row.
       \return The number of ids in the Ids list.
   */

  LNO getRowListView(GID *&Ids, LID *&localIds, LNO *nnz)
  {
    Ids = NULL;
    localIds = NULL;
    nnz = NULL;
    return 0;
  }

  /*! Return the column Ids of the non-zeros for the given row.
      \param Id  global Id for a row on this process
      \param localId  app's local Id, if any, associated with this row
      \param columnId on return will contain the list of global column Ids
   */
  virtual void getRowNonZeroCopy(GID Id, LID localId,
    std::vector<GID> &columnId);

  /*! Obtain a read-only view, if possible, of the column Ids of the
      input row.
      \param Id  global Id for a row on this process
      \param localId  if input adapter supplied local Ids, this
         is that localId
      \param columnId on return will point a list of global column global Ids.
      \return The number of ids in the columnId list.
   */
  LNO getRowNonZeroView(GID Id, LID localId, GID *&columnId) const
  {
    columnId = NULL;
    return 0;
  }

  /*! Return true of matrix is globally lower triangular.
   */
  virtual bool isLowerTriangular() const;

  /*! Return true of matrix is globally upper triangular.
   */
  virtual bool isUpperTriangular() const;

  /*! Return true of matrix globally has any diagonal entries.
   */
  virtual bool hasDiagonalEntries() const;
};
  
  
}  //namespace Zoltan2
  
#endif

// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_XpetraCrsMatrixInput.hpp

    \brief An input adapter for a Xpetra::CrsMatrix.
*/

#ifndef _ZOLTAN2_XPETRACRSMATRIXINPUT_HPP_
#define _ZOLTAN2_XPETRACRSMATRIXINPUT_HPP_

#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Zoltan2_MatrixInput.hpp>
#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_Util.hpp>

namespace Zoltan2 {

//////////////////////////////////////////////////////////////////////////////
/*! Zoltan2::XpetraCrsMatrixInput
    \brief Provides access for Zoltan2 to Xpetra::CrsMatrix data.

    TODO: we assume FillComplete has been called.  We should support
                objects that are not FillCompleted.

    The template parameter is the user's input object - an Epetra
    matrix or a templated Tpetra matrix 
    or a templated Xpetra::CrsMatrix.
*/

template <typename User>
class XpetraCrsMatrixInput : public MatrixInput<User> {
public:

  typedef typename InputAdapter<User>::scalar_t scalar_t;
  typedef typename InputAdapter<User>::lno_t    lno_t;
  typedef typename InputAdapter<User>::gno_t    gno_t;
  typedef typename InputAdapter<User>::lid_t    lid_t;
  typedef typename InputAdapter<User>::gid_t    gid_t;
  typedef typename InputAdapter<User>::node_t   node_t;

  typedef Epetra_CrsMatrix ematrix_t;
  typedef Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> tmatrix_t;
  typedef Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> xmatrix_t;
  typedef Xpetra::TpetraCrsMatrix<scalar_t, lno_t, gno_t, node_t> xtmatrix_t;
  typedef Xpetra::EpetraCrsMatrix xematrix_t;

  // TODO - add this value to the traits.
  enum InputAdapterType inputAdapterType() {return XpetraCrsMatrixAdapterType;}

  /*! Name of input adapter type
   */
  std::string inputAdapterName()const {return std::string("XpetraCrsMatrix");}

  /*! Destructor
   */
  ~XpetraCrsMatrixInput() { }

  /*! Constructor   
   */
  // Constructor 
  XpetraCrsMatrixInput(const RCP<const User> &inmatrix):
    inmatrix_(inmatrix), 
    matrix_(),
    rowMap_(),
    colMap_(),
    base_(),
    offset_(),
    columnIds_()
  {
    matrix_ = XpetraTraits<User>::convertToXpetra(inmatrix);
    rowMap_ = matrix_->getRowMap();
    colMap_ = matrix_->getColMap();
    base_ = rowMap_->getIndexBase();

    size_t nrows = matrix_->getNodeNumRows();
    size_t nnz = matrix_->getNodeNumEntries();
 
    offset_.resize(nrows+1, lid_t(0));
    columnIds_.resize(nnz);
    ArrayView<const lid_t> indices;
    ArrayView<const scalar_t> nzs;
    lid_t next = 0;
    for (size_t i=0; i < nrows; i++){
      lid_t row = i + base_;
      lid_t nnz = matrix_->getNumEntriesInLocalRow(row);
      matrix_->getLocalRowView(row, indices, nzs);
      for (lid_t j=0; j < nnz; j++){
        // TODO - this will be slow
        //   Is it possible that global columns ids might be stored in order?
        columnIds_[next++] = colMap_->getGlobalElement(indices[j]);
      }
      offset_[i+1] = offset_[i] + nnz;
    } 
  };

  ////////////////////////////////////////////////////
  // The MatrixInput interface.
  ////////////////////////////////////////////////////

  /*! Returns the number rows on this process.
   */
  size_t getLocalNumRows() const { 
    return matrix_->getNodeNumRows();
  }

  /*! Returns the number rows in the entire matrix.
   */
  global_size_t getGlobalNumRows() const { 
    return matrix_->getGlobalNumRows();
  }

  /*! Return whether input adapter wants to use local IDs.
   */

  bool haveLocalIds() const {return true;}

  /*! Return whether local ids are consecutive and if so the base.
   */

  bool haveConsecutiveLocalIds (size_t &base) const
  {
    base = static_cast<size_t>(base_);
    return true;
  }

  /*! Returns the number columns on this process.
   */
  size_t getLocalNumColumns() const { 
    return matrix_->getNodeNumCols();
  }

  /*! Returns the number columns on this entire matrix.
   *    what about directional columns, count twice?
   */
  global_size_t getGlobalNumColumns() const { 
    return matrix_->getGlobalNumCols();
  }

  /*! Get copy of matrix entries on local process
   */
  void getRowListCopy(std::vector<gid_t> &rowIds,
    std::vector<lid_t> &localIds, std::vector<lno_t> &offsets,
    std::vector<gid_t> &colIds) const
  {
    size_t nrows = getLocalNumRows();
    size_t nnz = matrix_->getNodeNumEntries();
    size_t maxrow = matrix_->getNodeMaxNumRowEntries();
    size_t next = 0;

    rowIds.resize(nrows);
    offsets.resize(nrows+1);
    colIds.resize(nnz);
    localIds.clear();   // consecutive integer IDs implied

    Teuchos::Array<lno_t> indices(maxrow);
    Teuchos::Array<scalar_t> nzs(maxrow);

    offsets[0] = 0;

    for (size_t i=0; i < nrows; i++){
      lno_t row = i + base_;
      lno_t nnz = matrix_->getNumEntriesInLocalRow(row);
      size_t n;
      matrix_->getLocalRowCopy(row, indices.view(0,nnz), nzs.view(0,nnz), n);
      for (lno_t j=0; j < nnz; j++){
        colIds[next++] = colMap_->getGlobalElement(indices[j]);
      }
      rowIds[i] = rowMap_->getGlobalElement(row);
      offsets[i+1] = offsets[i] + nnz;
    }
  } 

  /*! Access to xpetra matrix
   */

  const RCP<const xmatrix_t> &getMatrix() const
  {
    return matrix_;
  }

  /*! Return a read only view of the data.
     \param rowIds  Global row ids.  The memory for the global 
          row IDs persists until the underlying Xpetra::CrsMatrix is deleted.
     \param localIds on return is NULL, signifying that local IDs are
           contiguous integers starting at 0.
     \param offsets The columns for rowIds[i] begin at colIds[offsets[i]].  There are
           numRows+1 offsets.  The last offset is the length of the colIds array.
           The memory pointed to by offsets persists in the GraphModel is deleted.
     \param colIds The global column Ids. The memory pointed to by colIds 
          persists until the GraphModel is deleted.

     \return  The number rows in the rowIds list is returned.
   */
  //lno_t getRowListView(gid_t *&rowIds, lid_t *&localIds,
  //  lno_t *&rowSize, gid_t *& colIds) const

  size_t getRowListView(const gid_t *&rowIds, const lid_t *&localIds,
    const lid_t *&offsets, const gid_t *& colIds) const
  {
    size_t nrows = getLocalNumRows();

    ArrayView<const gid_t> rowView = rowMap_->getNodeElementList();
    rowIds = rowView.getRawPtr();
   
    localIds = NULL;   // Implies consecutive integers

    offsets = offset_.getRawPtr();
    colIds = columnIds_.getRawPtr();
    return nrows;
  }

#if 0
  /*! Apply a partitioning solution to the matrix.
   *   Return a new matrix reflecting the solution.
   */
  User *applyPartitioningSolution(const User &in,
    lno_t numIds, lno_t numParts, 
    const gid_t *gid, const lid_t *lid, const lno_t *partition)
  {
    RCP<User> inptr = rcp(&in, false);
    std::string userType(InputTraits<User>::name());

    bool isEpetra = (userType == std::string("Epetra_CrsMatrix"));
    bool isTpetra = (userType == std::string("Tpetra::CrsMatrix"));
    bool isXpetra = (userType == std::string("Xpetra::CrsMatrix"));

    RCP<xmatrix_t> inMatrix;

    if (isEpetra){
      RCP<xematrix_t> ematrix = rcp(new xematrix_t(inptr));
      inMatrix = rcp_implicit_cast<xmatrix_t>(ematrix);
    }
    else if (isTpetra){
      RCP<xtmatrix_t > tmatrix = rcp(new xtmatrix_t(inptr));
      inMatrix = rcp_implicit_cast<xmatrix_t>(tmatrix);
    }
    else if (isXpetra){
      inMatrix = inptr;
    }
    else{
      throw std::logic_error("should be epetra, tpetra or xpetra");
    }

    RCP<const Comm<int> > comm = inMatrix->getRowMap()->getComm();

    size_t numLocalRows = inMatrix->getNodeNumRows();
    size_t numGlobalRows = inMatrix->getGlobalNumRows();
    gno_t base = inMatrix->getRowMap()->getIndexBase();
    size_t numNewRows = 0;

    ArrayRCP<lno_t> rowSizes(numIds);
    for (lno_t i=0; i < numIds; i++){
      rowSizes[i] = inMatrix->getNumEntriesInLocalRow(lid[i]);
    }

    ArrayRCP<lno_t> partArray(
      const_cast<lno_t *>(partition), 0, numLocalRows, false);
    ArrayRCP<gid_t> gidArray(const_cast<gid_t *>(gid), 0, numLocalRows, false);
    ArrayRCP<gid_t> newGidArray;
    ArrayRCP<size_t> newRowSizes;   // Epetra uses int, Tpetra uses size_t.
  
    try{
      numNewRows = convertPartitionListToImportList(comm,
        partArray, gidArray, rowSizes, newGidArray, newRowSizes);
    }
    catch (std::exception &e){
      Z2_THROW_ZOLTAN2_ERROR(env, e);
    }

    // We work with the underlying Tpetra or Epetra object.  This
    // is because when we create xpetra maps and importers and
    // a new matrix, they must be of a concrete type, not of type Xpetra.

    if (isXpetra){
      Xpetra::UnderlyingLib l = inMatrix->getRowMap()->lib();
      if (l == Xpetra::UseEpetra)
        isEpetra = true;
      else
        isTpetra = true;
    }

    if (isEpetra){
      const RCP<const Epetra_Comm> ecomm = Xpetra::toEpetra(comm);
    
      // source map
      RCP<xematrix_t> xematrix = rcp_dynamic_cast<xematrix_t>(inMatrix);
      RCP<const Epetra_CrsMatrix> ematrix = xematrix->getEpetra_CrsMatrix();
      const Epetra_Map &sMap = ematrix->RowMap();

      // target map
      // TODO - doesn't compile if GIDs are not ints.
      Epetra_Map tMap(numGlobalRows, numNewRows, newGidArray.getRawPtr(),
        1, base, ecomm);

      // target matrix
    
      Array<int> intRowSizes;
      int *nnz=NULL;
    
      if (sizeof(size_t) != sizeof(int)){
        intRowSizes = Array<int>(numNewRows);
        for (lno_t i=0; i < numNewRows; i++){
          intRowSizes[i] = static_cast<int>(newRowSizes[i]);
        }
        nnz = intRowSizes.getRawPtr();
      }
      else{ 
        nnz = static_cast<int *>(newRowSizes.getRawPtr());
      }

      Epetra_CrsMatrix *M = new Epetra_CrsMatrix(Copy, tMap, nnz, true);
      Epetra_Import importer(tMap, sMap);

      M->Import(*ematrix, importer, Insert);

      return M;
    }
    else{
      typedef typename Tpetra::Map<lno_t,gno_t,node_t> tmap_t;

      // source map
      RCP<xtmatrix_t> xtmatrix = rcp_dynamic_cast<xtmatrix_t>(inMatrix);
      RCP<const tmatrix_t> tmatrix = xtmatrix->getTpetra_CrsMatrix();
      const RCP<const tmap_t> &sMap = tmatrix->getRowMap();

      // target map
      RCP<const tmap_t> tMap =
        rcp(new tmap_t(numGlobalRows, newGidArray.view(0, numNewRows), base, comm));

      // target matrix

      tmatrix_t *M = new tmatrix_t(tMap, newRowSizes, Tpetra::StaticProfile);

      Tpetra::Import<lno_t, gno_t, node_t> importer(sMap, tMap);

      M->doImport(*tmatrix, importer, Tpetra::INSERT);

      return M;
    }
  }
#endif


private:

  RCP<const User> inmatrix_;
  RCP<const xmatrix_t> matrix_;
  RCP<const Xpetra::Map<lno_t, gno_t, node_t> > rowMap_;
  RCP<const Xpetra::Map<lno_t, gno_t, node_t> > colMap_;
  lno_t base_;
  ArrayRCP<lno_t> offset_;
  ArrayRCP<gno_t> columnIds_;

};
  
}  //namespace Zoltan2
  
#endif

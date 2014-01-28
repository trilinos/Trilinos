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

/*! \file Zoltan2_MatrixAdapter.hpp
    \brief Defines the MatrixAdapter interface.
*/

#ifndef _ZOLTAN2_MATRIXADAPTER_HPP_
#define _ZOLTAN2_MATRIXADAPTER_HPP_

#include <Zoltan2_Adapter.hpp>

namespace Zoltan2 {

enum MatrixEntityType {
  MATRIX_ROW,
  MATRIX_COLUMN,
  MATRIX_NONZERO
};

/*!  \brief MatrixAdapter defines the adapter interface for matrices.

    Adapter objects provide access for Zoltan2 to the user's data.
    Many built-in adapters are already defined for common data structures,
    such as Tpetra and Epetra objects and C-language pointers to arrays.

    Data types:
    \li \c scalar_t row, column or non-zero weights
    \li \c lno_t    local indices and local counts
    \li \c gno_t    global indices and global counts
    \li \c gid_t    application global Ids
    \li \c node_t is a sub class of KokkosClassic::StandardNodeMemoryModel

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

     \todo Create BasicCrsMatrixAdapter subclass
     \todo Do we want to require adapters to give us the global
               number of rows, columns etc?  We can figure that out.
      \todo  This is a row-oriented matrix.  Do we need a column-oriented
              matrix?  In particular - we assumed coordinates are for rows.
      \todo  If the user can tell us there are no diagonal entries
        in a square matrix, it can save time if we have to remove
        them for the algorithm.  Should we have a set method in 
        subclasses for setMatrixHasDiagonalEntries yes, no and maybe?
*/

template <typename User>
  class MatrixAdapter : public BaseAdapter<User> {
private:
  enum MatrixEntityType primaryEntityType;

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef User user_t;
#endif

  enum BaseAdapterType adapterType() const {return MatrixAdapterType;}

  // Constructor; sets default primaryEntityType to MATRIX_ROW.
  MatrixAdapter() : primaryEntityType(MATRIX_ROW) {}

  /*! \brief Destructor
   */
  virtual ~MatrixAdapter(){}

  /*! \brief Returns the number of rows on this process.
   */
  virtual size_t getLocalNumRows() const = 0;

  /*! \brief Returns the number of columns on this process.
   */
  virtual size_t getLocalNumColumns() const = 0;

  /*! \brief Returns the number of nonzeros on this process.
   */
  virtual size_t getLocalNumEntries() const = 0;



  /*! \brief Indicates whether the MatrixAdapter implements a view of the 
             matrix in compressed sparse row (CRS) format.
             All matrix adapters must implement either getCRSView or 
             getCCSView, but implementation of both is not required.
   */
  virtual bool CRSViewAvailable() const { return false; }

  /*! \brief Sets pointer to this process' rows' global IDs.
      \param rowIds will on return a pointer to row global Ids
   */
  virtual void getRowIDsView(const gid_t *&rowIds) const
  {
    rowIds = NULL;
    Z2_THROW_NOT_IMPLEMENTED_ERROR
  }

  /*! \brief Sets pointers to this process' matrix entries using
             compressed sparse row (CRS) format.
             All matrix adapters must implement either getCRSView or 
             getCCSView, but implementation of both is not required.
      \param offsets is an array of size numRows + 1.  The column Ids for
          rowIds[i] (returned by getRowIDsView) 
          begin at colIds[offsets[i]].  The last element of offsets
          is the size of the colIds array.
      \param colIds on return will point to the global column Ids for
         the non-zeros for each row.
   */
  virtual void getCRSView(const lno_t *&offsets,
                          const gid_t *&colIds) const 
  {
    // Default implementation; no CRS view provided.
    offsets = NULL;
    colIds = NULL;
    Z2_THROW_NOT_IMPLEMENTED_ERROR
  }

  /*! \brief Sets pointers to this process' matrix entries 
             and their values using
             compressed sparse row (CRS) format.
             All matrix adapters must implement either getCRSView or 
             getCCSView, but implementation of both is not required.
      \param offsets is an array of size numRows + 1.  The column Ids for
          rowIds[i] (returned by getRowIDsView) 
          begin at colIds[offsets[i]].  The last element of offsets
          is the size of the colIds array.
      \param colIds on return will point to the global column Ids for
         the non-zeros for each row.
      \param values on return will point to the values stored in the 
         non-zeros for each row.
   */
  virtual void getCRSView(const lno_t *&offsets,
                          const gid_t *& colIds,
                          const scalar_t *&values) const 
  {
    // Default implementation; no CRS view provided.
    offsets = NULL;
    colIds = NULL;
    values = NULL;
    Z2_THROW_NOT_IMPLEMENTED_ERROR
  }

  /*! \brief Returns the number of weights per row (0 or greater).
      Row weights may be used when partitioning matrix rows.
   */
  virtual int getNumWeightsPerRow() const { return 0;}

  /*! \brief  Provide a pointer to the row weights, if any.
      \param weights is the list of weights with a given index for
           the rows returned in getRowIDsView().  If weights for
           this index are to be uniform for all rows in the
           global problem, the \c weights should be a NULL pointer.
      \param stride The k'th weight is located at weights[stride*k]
      \param idx ranges from zero to one less than getNumWeightsPerRow().
   */
  virtual void getRowWeightsView(const scalar_t *&weights, int &stride,
                                 int idx = 0) const
  {
    // Default implementation
    weights = NULL;
    stride = 0;
    Z2_THROW_NOT_IMPLEMENTED_ERROR
  }

  /*! \brief Indicate whether row weight with index idx should be the
   *         global number of nonzeros in the row.
   */
  virtual bool useNumNonzerosAsRowWeight(int idx) const { return 0; }



  /*! \brief Indicates whether the MatrixAdapter implements a view of the 
             matrix in compressed sparse column (CCS) format.
             All matrix adapters must implement either getCRSView or 
             getCCSView, but implementation of both is not required.
   */
  virtual bool CCSViewAvailable() const { return false; }

  /*! \brief Sets pointer to this process' columns' global IDs.
      \param colIds will on return a pointer to column global Ids
   */
  virtual void getColumnIDsView(const gid_t *&colIds) const
  {
    colIds = NULL;
    Z2_THROW_NOT_IMPLEMENTED_ERROR
  }

  /*! \brief Sets pointers to this process' matrix entries using
             compressed sparse column (CCS) format.
             All matrix adapters must implement either getCRSView or 
             getCCSView, but implementation of both is not required.
      \param offsets is an array of size numCols + 1.  The row Ids for
          colIds[i] (returned by getColumnIDsView) 
          begin at rowIds[offsets[i]].  The last element of offsets
          is the size of the rowIds array.
      \param rowIds on return will point to the global row Ids for
         the non-zeros for each column.
   */
  virtual void getCCSView(const lno_t *&offsets,
                          const gid_t *&rowIds) const 
  {
    // Default implementation; no CCS view provided.
    offsets = NULL;
    rowIds = NULL;
    Z2_THROW_NOT_IMPLEMENTED_ERROR
  }

  /*! \brief Sets pointers to this process' matrix entries 
             and their values using
             compressed sparse column (CCS) format.
             All matrix adapters must implement either getCRSView or 
             getCCSView, but implementation of both is not required.
      \param offsets is an array of size numCols + 1.  The row Ids for
          colIds[i] (returned by getColumnIDsView) 
          begin at rowIds[offsets[i]].  The last element of offsets
          is the size of the rowIds array.
      \param rowIds on return will point to the global row Ids for
         the non-zeros for each column.
      \param values on return will point to the values stored in the
         non-zeros for each column.
   */
  virtual void getCCSView(const lno_t *&offsets,
                          const gid_t *&rowIds,
                          const scalar_t *&values) const
  {
    // Default implementation; no CCS view provided.
    offsets = NULL;
    rowIds = NULL;
    values = NULL;
    Z2_THROW_NOT_IMPLEMENTED_ERROR
  }

  /*! \brief Returns the number of weights per column (0 or greater).
      Column weights may be used when partitioning matrix columns.
   */
  virtual int getNumWeightsPerColumn() const { return 0; }

  /*! \brief  Provide a pointer to the column weights, if any.
      \param weights is the list of weights with a given index for
           the columns returned in getColumnIDsView().  If weights for
           this index are to be uniform for all columns in the
           global problem, the \c weights should be a NULL pointer.
      \param stride The k'th weight is located at weights[stride*k]
      \param idx ranges from zero to one less than getNumWeightsPerColumn().
   */
  virtual void getColumnWeightsView(const scalar_t *&weights, int &stride,
                                    int idx = 0) const 
  {
    // Default implementation
    weights = NULL;
    stride = 0;
    Z2_THROW_NOT_IMPLEMENTED_ERROR
  }

  /*! \brief Indicate whether column weight with index idx should be the
   *         global number of nonzeros in the column.
   */
  virtual bool useNumNonzerosAsColumnWeight(int idx) const { return 0; }

#ifdef FUTURE_FEATURE
  /*! method saying whether the matrix is using symmetric storage; that is,
   *  for symmetric matrices, is only the upper or lower triangular matrix
   *  stored, or is the entire matrix stored?
   */
  virtual bool symmetricStorage() const {return false;}
#endif

  /*! \brief Returns the dimension of the geometry, if any.
   *
   *  Some algorithms can use geometric row or column coordinate
   *    information if it is present.  Given the problem parameters
   *    supplied by the user, it may make sense to use row coordinates
   *    or it may make sense to use column coordinates.
   */
//KDDDEC  Instead of having coordinate-based functions, 
//KDDDEC  see if can have a method SetCoordinateInput that
//KDDDEC  allows user to provide a coordinate input adapter,
//KDDDEC  and just save it.

  virtual int getCoordinateDimension() const { return 0; }

  /*! \brief Provide a pointer to one dimension of row coordinates.
      \param coords  [output] points to a list of coordinate values for the 
             dimension.
             The order of \c coords should correspond to the order of \c rowIds
             in getRowIDsView().
      \param stride  [output] describes the layout of the coordinate values in
              the coords list.  If stride is one, then the ith coordinate
              value is coords[i], but if stride is two, then the
              ith coordinate value is coords[2*i].
      \param dim [input] is a value from 0 to one less than
         getCoordinateDimension() specifying which dimension is
         being provided in the coords list.
   */

  virtual void getRowCoordinatesView(const scalar_t *&coords, int &stride,
                                     int dim) const 
  {
    coords = NULL;
    stride = 0;
    Z2_THROW_NOT_IMPLEMENTED_ERROR
  }

  /*! \brief Provide a pointer to one dimension of column coordinates.
      \param coords  [output] points to a list of coordinate values for the
             dimension.
             The order of \c coords should correspond to the order of \c colIds
             in getColumnIDsView().
      \param stride  [output] describes the layout of the coordinate values in
              the coords list.  If stride is one, then the ith coordinate
              value is coords[i], but if stride is two, then the
              ith coordinate value is coords[2*i].
      \param dim [input] is a value from 0 to one less than
         getDimension() specifying which dimension is
         being provided in the coords list.
   */

  virtual void getColumnCoordinatesView(const scalar_t *&coords, int &stride,
                                        int dim) const 
  {
    coords = NULL;
    stride = 0;
    Z2_THROW_NOT_IMPLEMENTED_ERROR
  }


  ////////////////////////////////////////////////////////////////////////////
  // Implementations of base-class methods and other methods shared by all

  /*! \brief Returns the entity to be partitioned, ordered, colored, etc.
   *  Valid values are MATRIX_ROW, MATRIX_COLUMN, MATRIX_NONZERO
   */
  inline enum MatrixEntityType getPrimaryEntityType() const {
    return this->primaryEntityType;
  }

  /*! \brief Sets the primary entity type.  Called by algorithm based on
   *  parameter value in parameter list from application.
   *  Also sets to adjacencyEntityType to something reasonable:  opposite of
   *  primaryEntityType.
   */
  void setPrimaryEntityType(string typestr) {
    if (typestr == "row") {
      this->primaryEntityType = MATRIX_ROW;
    }
    else if (typestr == "column") {
      this->primaryEntityType = MATRIX_COLUMN;
    }
    else if (typestr == "nonzero") {
      this->primaryEntityType = MATRIX_NONZERO;
    }
    else {
      std::ostringstream emsg;
      emsg << __FILE__ << "," << __LINE__
           << " error:  Invalid MatrixEntityType " << typestr << std::endl;
      emsg << "Valid values are 'row', 'column' and 'nonzero'." << std::endl;
      throw std::runtime_error(emsg.str());
    }
  }

  // Functions from the BaseAdapter interface
  size_t getLocalNumIDs() const {
    switch (getPrimaryEntityType()) {
    case MATRIX_ROW:
      return getLocalNumRows();
    case MATRIX_COLUMN:
      return getLocalNumColumns();
    case MATRIX_NONZERO:
      return getLocalNumEntries();
    default:   // Shouldn't reach default; just making compiler happy
      return 0;
    }
  }

  void getIDsView(const gid_t *&Ids) const {
    switch (getPrimaryEntityType()) {
    case MATRIX_ROW:
      getRowIDsView(Ids);
      break;
    case MATRIX_COLUMN:
      getColumnIDsView(Ids);
      break;
    case MATRIX_NONZERO: {
      // TODO:  Need getNonzeroIDsView?  What is a Nonzero ID?  
      // TODO:  std::pair<gid_t, gid_t>?
      std::ostringstream emsg;
      emsg << __FILE__ << "," << __LINE__
           << " error:  getIDsView not yet supported for matrix nonzeros." 
           << std::endl;
      throw std::runtime_error(emsg.str());
      break;
      }
    default:   // Shouldn't reach default; just making compiler happy
      break;
    }
  }

  int getNumWeightsPerID() const {
    switch (getPrimaryEntityType()) {
    case MATRIX_ROW:
      return getNumWeightsPerRow();
    case MATRIX_COLUMN:
      return getNumWeightsPerColumn();
    case MATRIX_NONZERO:
      return 0;  //TODO: weights not yet supported for nonzeros
    default:   // Shouldn't reach default; just making compiler happy
      return 0;
    }
  }

  void getWeightsView(const scalar_t *&wgt, int &stride, int idx = 0) const {
    switch (getPrimaryEntityType()) {
    case MATRIX_ROW:
      getRowWeightsView(wgt, stride, idx);
      break;
    case MATRIX_COLUMN:
      getColumnWeightsView(wgt, stride, idx);
      break;
    case MATRIX_NONZERO:
      {
      // TODO:  Need getNonzeroWeightsView with Nonzeros as primary object? 
      // TODO:  That is, get Nonzeros' weights based on some nonzero ID?
      std::ostringstream emsg;
      emsg << __FILE__ << "," << __LINE__
           << " error:  getWeightsView not yet supported for matrix nonzeros." 
           << std::endl;
      throw std::runtime_error(emsg.str());
      break;
      }
    default:   // Shouldn't reach default; just making compiler happy
      break;
    }
  }

  void getCoordinatesView(const scalar_t *&coords, int &stride, int dim) const
  {
    switch (getPrimaryEntityType()) {
    case MATRIX_ROW:
      getRowCoordinatesView(coords, stride, dim);
      break;
    case MATRIX_COLUMN:
      getColumnCoordinatesView(coords, stride, dim);
      break;
    case MATRIX_NONZERO: {
      // TODO:  Need getCoordinatesView with Nonzeros as primary object? 
      // TODO:  Could return (i,j), but need an ordering for nonzeros.
      std::ostringstream emsg;
      emsg << __FILE__ << "," << __LINE__
           << " error:  getCoordinatesView not supported for matrix nonzeros." 
           << std::endl;
      throw std::runtime_error(emsg.str());
      break;
      }
    default:   // Shouldn't reach default; just making compiler happy
      break;
    }
  }
};
  
}  //namespace Zoltan2
  
#endif

// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_MatrixAdapter.hpp
    \brief Defines the MatrixAdapter interface.
*/

#ifndef _ZOLTAN2_MATRIXADAPTER_HPP_
#define _ZOLTAN2_MATRIXADAPTER_HPP_

#include <Zoltan2_Adapter.hpp>
#include <Zoltan2_VectorAdapter.hpp>

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
    \li \c node_t   is a Kokkos Node type

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

template <typename User, typename UserCoord=User>
  class MatrixAdapter : public AdapterWithCoordsWrapper<User, UserCoord> {
private:
  enum MatrixEntityType primaryEntityType_;
  VectorAdapter<UserCoord> *coordinateInput_;
  bool haveCoordinateInput_;

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  using scalar_t = typename InputTraits<User>::scalar_t;
  using lno_t = typename InputTraits<User>::lno_t;
  using gno_t = typename InputTraits<User>::gno_t;
  using part_t = typename InputTraits<User>::part_t;
  using node_t = typename InputTraits<User>::node_t;
  using offset_t = typename InputTraits<User>::offset_t;
  using user_t = User;
  using userCoord_t = UserCoord;
  using base_adapter_t = MatrixAdapter<User,UserCoord>;
  using device_t = typename node_t::device_type;
#endif

  enum BaseAdapterType adapterType() const override {return MatrixAdapterType;}

  // Constructor; sets default primaryEntityType to MATRIX_ROW.
  MatrixAdapter() : primaryEntityType_(MATRIX_ROW),
                    coordinateInput_(),
                    haveCoordinateInput_(false) {}

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
  virtual void getRowIDsView(const gno_t *&rowIds) const
  {
    rowIds = NULL;
    Z2_THROW_NOT_IMPLEMENTED
  }

  virtual void getRowIDsHostView(typename BaseAdapter<User>::ConstIdsHostView& rowIds) const
  {
    Z2_THROW_NOT_IMPLEMENTED
  }

  virtual void getRowIDsDeviceView(typename BaseAdapter<User>::ConstIdsDeviceView& rowIds) const
  {
    Z2_THROW_NOT_IMPLEMENTED
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
  virtual void getCRSView(ArrayRCP<const offset_t> &offsets, ArrayRCP<const gno_t> &colIds) const
  {
    // Default implementation; no CRS view provided.
    offsets = ArrayRCP<const offset_t>();
    colIds = ArrayRCP<const gno_t>();
    Z2_THROW_NOT_IMPLEMENTED
  }

  virtual void getCRSHostView(typename BaseAdapter<User>::ConstOffsetsHostView& offsets,
                              typename BaseAdapter<User>::ConstIdsHostView& colIds) const
  {
      Z2_THROW_NOT_IMPLEMENTED
  }

  virtual void getCRSDeviceView(typename BaseAdapter<User>::ConstOffsetsDeviceView& offsets,
                                typename BaseAdapter<User>::ConstIdsDeviceView& colIds) const
  {
      Z2_THROW_NOT_IMPLEMENTED
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
  virtual void getCRSView(ArrayRCP<const offset_t> &offsets,
                          ArrayRCP<const gno_t> &colIds,
                          ArrayRCP<const scalar_t> &values) const
  {
    // Default implementation; no CRS view provided.
    offsets = ArrayRCP<const offset_t>();
    colIds = ArrayRCP<const gno_t>();
    values = ArrayRCP<const scalar_t>();
    Z2_THROW_NOT_IMPLEMENTED
  }

  virtual void getCRSHostView(typename BaseAdapter<User>::ConstOffsetsHostView& offsets,
                              typename BaseAdapter<User>::ConstIdsHostView& colIds,
                              typename BaseAdapter<User>::ConstScalarsHostView& values) const
  {
      Z2_THROW_NOT_IMPLEMENTED
  }

  virtual void getCRSDeviceView(typename BaseAdapter<User>::ConstOffsetsDeviceView& offsets,
                                typename BaseAdapter<User>::ConstIdsDeviceView& colIds,
                                typename BaseAdapter<User>::ConstScalarsDeviceView& values) const
  {
      Z2_THROW_NOT_IMPLEMENTED
  }


  /*! \brief Returns the number of weights per row (0 or greater).
      Row weights may be used when partitioning matrix rows.
   */
  virtual int getNumWeightsPerRow() const { return 0;}

  /*! \brief  Provide a pointer to the row weights, if any.
      \param weights is the list of weights with a given index for
           the rows returned in getRowIDsView().
      \param stride The k'th weight is located at weights[stride*k]
      \param idx ranges from zero to one less than getNumWeightsPerRow().
   */
  virtual void getRowWeightsView(const scalar_t *&weights, int &stride,
                                 int idx = 0) const
  {
    // Default implementation
    weights = NULL;
    stride = 0;
    Z2_THROW_NOT_IMPLEMENTED
  }

  virtual void getRowWeightsHostView(typename BaseAdapter<User>::WeightsHostView1D& weights,
                                     int /* idx */ = 0) const {
      Z2_THROW_NOT_IMPLEMENTED
  }

  virtual void getRowWeightsHostView(typename BaseAdapter<User>::WeightsHostView &weights) const {
    Z2_THROW_NOT_IMPLEMENTED
  }

  virtual void getRowWeightsDeviceView(typename BaseAdapter<User>::WeightsDeviceView1D& weights,
                                       int /* idx */ = 0) const {
      Z2_THROW_NOT_IMPLEMENTED
  }

  virtual void getRowWeightsDeviceView(typename BaseAdapter<User>::WeightsDeviceView &weights) const {
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief Indicate whether row weight with index idx should be the
   *         global number of nonzeros in the row.
   */
  virtual bool useNumNonzerosAsRowWeight(int idx) const
  {
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief Indicates whether the MatrixAdapter implements a view of the
             matrix in compressed sparse column (CCS) format.
             All matrix adapters must implement either getCRSView or
             getCCSView, but implementation of both is not required.
   */
  virtual bool CCSViewAvailable() const { return false; }

  /*! \brief Sets pointer to this process' columns' global IDs.
      \param colIds will on return a pointer to column global Ids
   */
  virtual void getColumnIDsView(const gno_t *&colIds) const
  {
    colIds = NULL;
    Z2_THROW_NOT_IMPLEMENTED
  }

  virtual void getColumnIDsHostView(typename BaseAdapter<User>::ConstIdsHostView& colIds) const
  {
    Z2_THROW_NOT_IMPLEMENTED
  }

  virtual void getColumnIDsDeviceView(typename BaseAdapter<User>::ConstIdsDeviceView& colIds) const
  {
    Z2_THROW_NOT_IMPLEMENTED
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
  virtual void getCCSView(ArrayRCP<const offset_t> &offsets,
                          ArrayRCP<const gno_t> &rowIds) const
  {
    // Default implementation; no CCS view provided.
    offsets = ArrayRCP<const offset_t>();
    rowIds = ArrayRCP<const gno_t>();
    Z2_THROW_NOT_IMPLEMENTED
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
  virtual void getCCSView(ArrayRCP<const offset_t> &offsets,
                          ArrayRCP<const gno_t> &rowIds,
                          ArrayRCP<const scalar_t> &values) const
  {
    // Default implementation; no CCS view provided.
    offsets = ArrayRCP<const offset_t>();
    rowIds = ArrayRCP<const gno_t>();
    values = ArrayRCP<const scalar_t>();
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief Returns the number of weights per column (0 or greater).
      Column weights may be used when partitioning matrix columns.
   */
  virtual int getNumWeightsPerColumn() const { return 0; }

  /*! \brief  Provide a pointer to the column weights, if any.
      \param weights is the list of weights with a given index for
           the columns returned in getColumnIDsView().
      \param stride The k'th weight is located at weights[stride*k]
      \param idx ranges from zero to one less than getNumWeightsPerColumn().
   */
  virtual void getColumnWeightsView(const scalar_t *&weights, int &stride,
                                    int idx = 0) const
  {
    // Default implementation
    weights = NULL;
    stride = 0;
    Z2_THROW_NOT_IMPLEMENTED
  }

  virtual void getColumnWeightsHostView(typename BaseAdapter<User>::WeightsHostView1D& weights,
                                     int /* idx */ = 0) const {
      Z2_THROW_NOT_IMPLEMENTED
  }

  virtual void getColumnWeightsHostView(typename BaseAdapter<User>::WeightsHostView &weights) const {
    Z2_THROW_NOT_IMPLEMENTED
  }

  virtual void getColumnWeightsDeviceView(typename BaseAdapter<User>::WeightsDeviceView1D& weights,
                                       int /* idx */ = 0) const {
      Z2_THROW_NOT_IMPLEMENTED
  }

  virtual void getColumnWeightsDeviceView(typename BaseAdapter<User>::WeightsDeviceView &weights) const {
    Z2_THROW_NOT_IMPLEMENTED
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

  /*! \brief Allow user to provide additional data that contains coordinate
   *         info associated with the MatrixAdapter's primaryEntityType.
   *         Associated data must have the same parallel distribution and
   *         ordering of entries as the primaryEntityType.
   *
   *  \param coordData is a pointer to a VectorAdapter with the user's
   *         coordinate data.
   */
  void setCoordinateInput(VectorAdapter<UserCoord> *coordData) override
  {
    coordinateInput_ = coordData;
    haveCoordinateInput_ = true;
  }

  /*! \brief Indicate whether coordinate information has been set for this
   *         MatrixAdapter
   */
  bool coordinatesAvailable() const { return haveCoordinateInput_; }

  /*! \brief Obtain the coordinate data registered by the user.
   *  \return pointer a VectorAdapter with the user's coordinate data.
   */
  VectorAdapter<UserCoord> *getCoordinateInput() const override
  {
    return coordinateInput_;
  }

  ////////////////////////////////////////////////////////////////////////////
  // Implementations of base-class methods and other methods shared by all

  /*! \brief Returns the entity to be partitioned, ordered, colored, etc.
   *  Valid values are MATRIX_ROW, MATRIX_COLUMN, MATRIX_NONZERO
   */
  inline enum MatrixEntityType getPrimaryEntityType() const
  {
    return this->primaryEntityType_;
  }

  /*! \brief Sets the primary entity type.  Called by algorithm based on
   *  parameter value in parameter list from application.
   *  Also sets to adjacencyEntityType to something reasonable:  opposite of
   *  primaryEntityType.
   */
  void setPrimaryEntityType(std::string typestr)
  {
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
  size_t getLocalNumIDs() const override
  {
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

  void getIDsView(const gno_t *&Ids) const override
  {
    switch (getPrimaryEntityType()) {
    case MATRIX_ROW:
      getRowIDsView(Ids);
      break;
    case MATRIX_COLUMN:
      getColumnIDsView(Ids);
      break;
    case MATRIX_NONZERO: {
      // TODO:  Need getNonzeroIDsView?  What is a Nonzero ID?
      // TODO:  std::pair<gno_t, gno_t>?
      std::ostringstream emsg;
      emsg << __FILE__ << "," << __LINE__
           << " error:  getIDsView not yet supported for matrix nonzeros."
           << std::endl;
      throw std::runtime_error(emsg.str());
      }
    default:   // Shouldn't reach default; just making compiler happy
      break;
    }
  }

  void getIDsHostView(typename BaseAdapter<User>::ConstIdsHostView& ids) const override {
    switch (getPrimaryEntityType()) {
    case MATRIX_ROW:
      getRowIDsHostView(ids);
      break;
    case MATRIX_COLUMN:
      getColumnIDsHostView(ids);
      break;
    case MATRIX_NONZERO: {
      // TODO:  Need getNonzeroIDsHostView?  What is a Nonzero ID?
      // TODO:  std::pair<gno_t, gno_t>?
      std::ostringstream emsg;
      emsg << __FILE__ << "," << __LINE__
           << " error:  getIDsView not yet supported for matrix nonzeros."
           << std::endl;
      throw std::runtime_error(emsg.str());
      }
    default:   // Shouldn't reach default; just making compiler happy
      break;
    }
  }

  void getIDsDeviceView(typename BaseAdapter<User>::ConstIdsDeviceView& ids) const override {
    switch (getPrimaryEntityType()) {
    case MATRIX_ROW:
      getRowIDsDeviceView(ids);
      break;
    case MATRIX_COLUMN:
      getColumnIDsDeviceView(ids);
      break;
    case MATRIX_NONZERO: {
      // TODO:  Need getNonzeroIDsDeviceView?  What is a Nonzero ID?
      // TODO:  std::pair<gno_t, gno_t>?
      std::ostringstream emsg;
      emsg << __FILE__ << "," << __LINE__
           << " error:  getIDsView not yet supported for matrix nonzeros."
           << std::endl;
      throw std::runtime_error(emsg.str());
      }
    default:   // Shouldn't reach default; just making compiler happy
      break;
    }
  }

  int getNumWeightsPerID() const override
  {
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

  void getWeightsView(const scalar_t *&wgt, int &stride,
                      int idx = 0) const override
  {
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
      }
    default:   // Shouldn't reach default; just making compiler happy
      break;
    }
  }

  void getWeightsHostView(typename BaseAdapter<User>::WeightsHostView &hostWgts) const override {
      switch (getPrimaryEntityType()) {
      case MATRIX_ROW:
        getRowWeightsHostView(hostWgts);
        break;
      case MATRIX_COLUMN:
        getColumnWeightsHostView(hostWgts);
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
        }
      default:   // Shouldn't reach default; just making compiler happy
        break;
      }  }

  void getWeightsHostView(typename BaseAdapter<User>::WeightsHostView1D &hostWgts,
                                  int idx = 0) const override {
      switch (getPrimaryEntityType()) {
      case MATRIX_ROW:
        getRowWeightsHostView(hostWgts, idx);
        break;
      case MATRIX_COLUMN:
        getColumnWeightsHostView(hostWgts, idx);
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
        }
      default:   // Shouldn't reach default; just making compiler happy
        break;
      }  }

  void getWeightsDeviceView(typename BaseAdapter<User>::WeightsDeviceView& deviceWgts) const override {
      switch (getPrimaryEntityType()) {
      case MATRIX_ROW:
        getRowWeightsDeviceView(deviceWgts);
        break;
      case MATRIX_COLUMN:
        getColumnWeightsDeviceView(deviceWgts);
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
        }
      default:   // Shouldn't reach default; just making compiler happy
        break;
      }  }

  void getWeightsDeviceView(typename BaseAdapter<User>::WeightsDeviceView1D& deviceWgts,
                                    int idx = 0) const override {
      switch (getPrimaryEntityType()) {
      case MATRIX_ROW:
        getRowWeightsDeviceView(deviceWgts, idx);
        break;
      case MATRIX_COLUMN:
        getColumnWeightsDeviceView(deviceWgts, idx);
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
        }
      default:   // Shouldn't reach default; just making compiler happy
        break;
      }  }

  bool useDegreeAsWeight(int idx) const
  {
    if (this->getPrimaryEntityType() == MATRIX_ROW)
      return useNumNonzerosAsRowWeight(idx);
    else {
      std::ostringstream emsg;
      emsg << __FILE__ << "," << __LINE__
           << " error:  useDegreeAsWeight is currently supported only for rows"
           << std::endl;
      throw std::runtime_error(emsg.str());
    }
  }
};

}  //namespace Zoltan2

#endif

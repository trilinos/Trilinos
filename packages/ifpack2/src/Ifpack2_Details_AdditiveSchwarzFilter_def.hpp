// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_ADDITIVE_SCHWARZ_FILTER_DEF_HPP
#define IFPACK2_ADDITIVE_SCHWARZ_FILTER_DEF_HPP

#include "Ifpack2_Details_AdditiveSchwarzFilter_decl.hpp"
#include "KokkosKernels_Sorting.hpp"
#include "KokkosSparse_SortCrs.hpp"
#include "Kokkos_Bitset.hpp"

namespace Ifpack2
{
namespace Details
{

template<class MatrixType>
AdditiveSchwarzFilter<MatrixType>::
AdditiveSchwarzFilter(const Teuchos::RCP<const row_matrix_type>& A_unfiltered,
    const Teuchos::ArrayRCP<local_ordinal_type>& perm,
    const Teuchos::ArrayRCP<local_ordinal_type>& reverseperm,
    bool filterSingletons)
{
  setup(A_unfiltered, perm, reverseperm, filterSingletons);
}

template<class MatrixType>
void AdditiveSchwarzFilter<MatrixType>::
setup(const Teuchos::RCP<const row_matrix_type>& A_unfiltered,
    const Teuchos::ArrayRCP<local_ordinal_type>& perm,
    const Teuchos::ArrayRCP<local_ordinal_type>& reverseperm,
    bool filterSingletons)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  //Check that A is an instance of allowed types, either:
  //  - Tpetra::CrsMatrix
  //  - Ifpack2::OverlappingRowMatrix (which always consists of two CrsMatrices, A_ and ExtMatrix_)
  TEUCHOS_TEST_FOR_EXCEPTION(
    !rcp_dynamic_cast<const crs_matrix_type>(A_unfiltered) &&
    !rcp_dynamic_cast<const overlapping_matrix_type>(A_unfiltered),
    std::invalid_argument,
    "Ifpack2::Details::AdditiveSchwarzFilter: The input matrix must be a Tpetra::CrsMatrix or an Ifpack2::OverlappingRowMatrix");
  A_unfiltered_ = A_unfiltered;
  local_matrix_type A_local, A_halo;
  bool haveHalo = !rcp_dynamic_cast<const overlapping_matrix_type>(A_unfiltered_).is_null();
  if(haveHalo)
  {
    auto A_overlapping = rcp_dynamic_cast<const overlapping_matrix_type>(A_unfiltered_); 
    A_local = A_overlapping->getUnderlyingMatrix()->getLocalMatrixDevice();
    A_halo = A_overlapping->getExtMatrix()->getLocalMatrixDevice();
  }
  else
  {
    auto A_crs = rcp_dynamic_cast<const crs_matrix_type>(A_unfiltered_); 
    A_local = A_crs->getLocalMatrixDevice();
    //leave A_halo empty in this case
  }
  //Check that perm and reversePerm are the same length and match the number of local rows in A
  TEUCHOS_TEST_FOR_EXCEPTION(
    perm.size() != reverseperm.size(),
    std::invalid_argument,
    "Ifpack2::Details::AdditiveSchwarzFilter: perm and reverseperm should be the same length");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (size_t) perm.size() != (size_t) A_unfiltered_->getLocalNumRows(),
    std::invalid_argument,
    "Ifpack2::Details::AdditiveSchwarzFilter: length of perm and reverseperm should be the same as the number of local rows in unfiltered A");
  //First, compute the permutation tables (that exclude singletons, if requested)
  FilterSingletons_ = filterSingletons;
  local_ordinal_type numLocalRows;
  local_ordinal_type totalLocalRows = A_local.numRows() + A_halo.numRows();
  if(FilterSingletons_)
  {
    //Fill singletons and singletonDiagonals (allocate them to the upper bound at first, then shrink it to size)
    singletons_ = Kokkos::DualView<local_ordinal_type*, device_type>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "singletons_"), totalLocalRows);
    singletonDiagonals_ = Kokkos::DualView<impl_scalar_type*, device_type>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "singletonDiagonals_"), totalLocalRows);
    auto singletonsDevice = singletons_.view_device();
    singletons_.modify_device();
    auto singletonDiagonalsDevice = singletonDiagonals_.view_device();
    singletonDiagonals_.modify_device();
    local_ordinal_type numSingletons;
    Kokkos::Bitset<device_type> isSingletonBitset(totalLocalRows);
    Kokkos::parallel_scan(policy_type(0, totalLocalRows),
        KOKKOS_LAMBDA(local_ordinal_type i, local_ordinal_type& lnumSingletons, bool finalPass)
        {
          bool isSingleton = true;
          impl_scalar_type singletonValue = Kokkos::ArithTraits<impl_scalar_type>::zero();
          if(i < A_local.numRows())
          {
            //row i is in original local matrix
            size_type rowBegin = A_local.graph.row_map(i);
            size_type rowEnd = A_local.graph.row_map(i + 1);
            for(size_type j = rowBegin; j < rowEnd; j++)
            {
              local_ordinal_type entry = A_local.graph.entries(j);
              if(entry < totalLocalRows && entry != i)
              {
                isSingleton = false;
                break;
              }
              if(finalPass && entry == i)
              {
                //note: using a sum to compute the diagonal value, in case entries are not compressed.
                singletonValue += A_local.values(j);
              }
            }
          }
          else
          {
            //row i is in halo
            local_ordinal_type row = i - A_local.numRows();
            size_type rowBegin = A_halo.graph.row_map(row);
            size_type rowEnd = A_halo.graph.row_map(row + 1);
            for(size_type j = rowBegin; j < rowEnd; j++)
            {
              local_ordinal_type entry = A_halo.graph.entries(j);
              if(entry < totalLocalRows && entry != i)
              {
                isSingleton = false;
                break;
              }
              if(finalPass && entry == i)
              {
                singletonValue += A_halo.values(j);
              }
            }
          }
          if(isSingleton)
          {
            if(finalPass)
            {
              isSingletonBitset.set(i);
              singletonsDevice(lnumSingletons) = i;
              singletonDiagonalsDevice(lnumSingletons) = singletonValue;
            }
            lnumSingletons++;
          }
        }, numSingletons);
    numSingletons_ = numSingletons;
    //Each local row in A_unfiltered is either a singleton or a row in the filtered matrix.
    numLocalRows = totalLocalRows - numSingletons_;
    //Using the list of singletons, create the reduced permutations.
    perm_ = Kokkos::DualView<local_ordinal_type*, device_type>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "perm_"), totalLocalRows);
    perm_.modify_device();
    reverseperm_ = Kokkos::DualView<local_ordinal_type*, device_type>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "perm_"), numLocalRows);
    reverseperm_.modify_device();
    auto permView = perm_.view_device();
    auto reversepermView = reverseperm_.view_device();
    //Get the original inverse permutation on device
    Kokkos::View<local_ordinal_type*, device_type> origpermView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "input reverse perm"), totalLocalRows);
    Kokkos::View<local_ordinal_type*, Kokkos::HostSpace> origpermHost(reverseperm.get(), totalLocalRows);
    Kokkos::deep_copy(execution_space(), origpermView, origpermHost);
    //reverseperm is a list of local rows in A_unfiltered, so filter out the elements of reverseperm where isSingleton is true. Then regenerate the forward permutation.
    Kokkos::parallel_scan(policy_type(0, totalLocalRows),
        KOKKOS_LAMBDA(local_ordinal_type i, local_ordinal_type& lindex, bool finalPass)
        {
          local_ordinal_type origRow = origpermView(i);
          if(!isSingletonBitset.test(origRow))
          {
            if(finalPass)
            {
              //mark the mapping in both directions between origRow and lindex (a row in filtered A)
              reversepermView(lindex) = origRow;
              permView(origRow) = lindex;
            }
            lindex++;
          }
          else
          {
            //origRow is a singleton, so mark a -1 in the new forward permutation to show that it has no corresponding row in filtered A.
            if(finalPass)
              permView(origRow) = local_ordinal_type(-1);
          }
        });
    perm_.sync_host();
    reverseperm_.sync_host();
  }
  else
  {
    //We will keep all the local rows, so use the permutation as is
    perm_ = Kokkos::DualView<local_ordinal_type*, device_type>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "perm_"), perm.size());
    perm_.modify_host();
    auto permHost = perm_.view_host();
    for(size_t i = 0; i < (size_t) perm.size(); i++)
      permHost(i) = perm[i];
    perm_.sync_device();
    reverseperm_ = Kokkos::DualView<local_ordinal_type*, device_type>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "reverseperm_"), reverseperm.size());
    reverseperm_.modify_host();
    auto reversepermHost = reverseperm_.view_host();
    for(size_t i = 0; i < (size_t) reverseperm.size(); i++)
      reversepermHost(i) = reverseperm[i];
    reverseperm_.sync_device();
    numSingletons_ = 0;
    numLocalRows = totalLocalRows;
  }
  auto permView = perm_.view_device();
  auto reversepermView = reverseperm_.view_device();
  //Fill the local matrix of A_ (filtered and permuted)
  //First, construct the rowmap by counting the entries in each row
  row_map_type localRowptrs(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Filtered rowptrs"), numLocalRows + 1);
  Kokkos::parallel_for(policy_type(0, numLocalRows + 1),
      KOKKOS_LAMBDA(local_ordinal_type i)
      {
        if(i == numLocalRows)
        {
          localRowptrs(i) = 0;
          return;
        }
        //Count the entries that will be in filtered row i.
        //This means entries which correspond to local, non-singleton rows/columns.
        local_ordinal_type numEntries = 0;
        local_ordinal_type origRow = reversepermView(i);
        if(origRow < A_local.numRows())
        {
          //row i is in original local matrix
          size_type rowBegin = A_local.graph.row_map(origRow);
          size_type rowEnd = A_local.graph.row_map(origRow + 1);
          for(size_type j = rowBegin; j < rowEnd; j++)
          {
            local_ordinal_type entry = A_local.graph.entries(j);
            if(entry < totalLocalRows && permView(entry) != -1)
              numEntries++;
          }
        }
        else
        {
          //row i is in halo
          local_ordinal_type row = origRow - A_local.numRows();
          size_type rowBegin = A_halo.graph.row_map(row);
          size_type rowEnd = A_halo.graph.row_map(row + 1);
          for(size_type j = rowBegin; j < rowEnd; j++)
          {
            local_ordinal_type entry = A_halo.graph.entries(j);
            if(entry < totalLocalRows && permView(entry) != -1)
              numEntries++;
          }
        }
        localRowptrs(i) = numEntries;
      });
  //Prefix sum to finish computing the rowptrs
  size_type numLocalEntries;
  Kokkos::parallel_scan(policy_type(0, numLocalRows + 1),
      KOKKOS_LAMBDA(local_ordinal_type i, size_type& lnumEntries, bool finalPass)
      {
        size_type numEnt = localRowptrs(i);
        if(finalPass)
          localRowptrs(i) = lnumEntries;
        lnumEntries += numEnt;
      }, numLocalEntries);
  //Allocate and fill local entries and values
  entries_type localEntries(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Filtered entries"), numLocalEntries);
  values_type localValues(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Filtered values"), numLocalEntries);
  //Create the local matrix with the uninitialized entries/values, then fill it.
  local_matrix_type localMatrix("AdditiveSchwarzFilter", numLocalRows, numLocalRows, numLocalEntries, localValues, localRowptrs, localEntries);
  fillLocalMatrix(localMatrix);
  //Create a serial comm and the map for the final filtered CrsMatrix (each process uses its own local map)
#ifdef HAVE_IFPACK2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! mapPairsAreFitted (*A_unfiltered_), std::invalid_argument, "Ifpack2::LocalFilter: "
    "A's Map pairs are not fitted to each other on Process "
    << A_->getRowMap ()->getComm ()->getRank () << " of the input matrix's "
    "communicator.  "
    "This means that LocalFilter does not currently know how to work with A.  "
    "This will change soon.  Please see discussion of Bug 5992.");
#endif // HAVE_IFPACK2_DEBUG

  // Build the local communicator (containing this process only).
  RCP<const Teuchos::Comm<int> > localComm;
#ifdef HAVE_MPI
  localComm = rcp (new Teuchos::MpiComm<int> (MPI_COMM_SELF));
#else
  localComm = rcp (new Teuchos::SerialComm<int> ());
#endif // HAVE_MPI
  //All 4 maps are the same for a local square operator.
  localMap_ = rcp (new map_type (numLocalRows, 0, localComm, Tpetra::GloballyDistributed));
  //Create the inner filtered matrix.
  auto crsParams = rcp(new Teuchos::ParameterList);
  crsParams->template set<bool>("sorted", true);
  //NOTE: this is the fastest way to construct A_ - it's created as fillComplete,
  //and no communication needs to be done since localMap_ uses a local comm.
  //It does need to copy the whole local matrix to host when DualViews are constructed
  //(only an issue with non-UVM GPU backends) but this is just unavoidable when creating a Tpetra::CrsMatrix.
  //It also needs to compute local constants (maxNumRowEntries) but this should be a
  //cheap operation relative to what this constructor already did.
  A_ = rcp(new crs_matrix_type(localMap_, localMap_, localMatrix, crsParams));
}

template<class MatrixType>
AdditiveSchwarzFilter<MatrixType>::~AdditiveSchwarzFilter () {}

template<class MatrixType>
void AdditiveSchwarzFilter<MatrixType>::updateMatrixValues()
{
  //Get the local matrix back from A_
  auto localMatrix = A_->getLocalMatrixDevice();
  //Copy new values from A_unfiltered to the local matrix, and then reconstruct A_.
  fillLocalMatrix(localMatrix);
  A_->setAllValues(localMatrix.graph.row_map, localMatrix.graph.entries, localMatrix.values);
}

template<class MatrixType>
Teuchos::RCP<const typename AdditiveSchwarzFilter<MatrixType>::crs_matrix_type>
AdditiveSchwarzFilter<MatrixType>::getFilteredMatrix() const
{
  return A_;
}

template<class MatrixType>
void AdditiveSchwarzFilter<MatrixType>::fillLocalMatrix(typename AdditiveSchwarzFilter<MatrixType>::local_matrix_type localMatrix)
{
  auto localRowptrs = localMatrix.graph.row_map;
  auto localEntries = localMatrix.graph.entries;
  auto localValues = localMatrix.values;
  auto permView = perm_.view_device();
  auto reversepermView = reverseperm_.view_device();
  local_matrix_type A_local, A_halo;
  bool haveHalo = !Teuchos::rcp_dynamic_cast<const overlapping_matrix_type>(A_unfiltered_).is_null();
  if(haveHalo)
  {
    auto A_overlapping = Teuchos::rcp_dynamic_cast<const overlapping_matrix_type>(A_unfiltered_); 
    A_local = A_overlapping->getUnderlyingMatrix()->getLocalMatrixDevice();
    A_halo = A_overlapping->getExtMatrix()->getLocalMatrixDevice();
  }
  else
  {
    auto A_crs = Teuchos::rcp_dynamic_cast<const crs_matrix_type>(A_unfiltered_); 
    A_local = A_crs->getLocalMatrixDevice();
    //leave A_halo empty in this case
  }
  local_ordinal_type totalLocalRows = A_local.numRows() + A_halo.numRows();
  //Fill entries and values
  Kokkos::parallel_for(policy_type(0, localMatrix.numRows()),
      KOKKOS_LAMBDA(local_ordinal_type i)
      {
        //Count the entries that will be in filtered row i.
        //This means entries which correspond to local, non-singleton rows/columns.
        size_type outRowStart = localRowptrs(i);
        local_ordinal_type insertPos = 0;
        local_ordinal_type origRow = reversepermView(i);
        if(origRow < A_local.numRows())
        {
          //row i is in original local matrix
          size_type rowBegin = A_local.graph.row_map(origRow);
          size_type rowEnd = A_local.graph.row_map(origRow + 1);
          for(size_type j = rowBegin; j < rowEnd; j++)
          {
            local_ordinal_type entry = A_local.graph.entries(j);
            if(entry < totalLocalRows)
            {
              auto newCol = permView(entry);
              if(newCol != -1)
              {
                localEntries(outRowStart + insertPos) = newCol;
                localValues(outRowStart + insertPos) = A_local.values(j);
                insertPos++;
              }
            }
          }
        }
        else
        {
          //row i is in halo
          local_ordinal_type row = origRow - A_local.numRows();
          size_type rowBegin = A_halo.graph.row_map(row);
          size_type rowEnd = A_halo.graph.row_map(row + 1);
          for(size_type j = rowBegin; j < rowEnd; j++)
          {
            local_ordinal_type entry = A_halo.graph.entries(j);
            if(entry < totalLocalRows)
            {
              auto newCol = permView(entry);
              if(newCol != -1)
              {
                localEntries(outRowStart + insertPos) = newCol;
                localValues(outRowStart + insertPos) = A_halo.values(j);
                insertPos++;
              }
            }
          }
        }
      });
  //Sort the matrix, since the reordering can shuffle the entries into any order.
  KokkosSparse::sort_crs_matrix(localMatrix);
}

template<class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> > AdditiveSchwarzFilter<MatrixType>::getComm() const
{
  return localMap_->getComm();
}

template<class MatrixType>
Teuchos::RCP<const typename AdditiveSchwarzFilter<MatrixType>::map_type>
AdditiveSchwarzFilter<MatrixType>::getRowMap() const
{
  return localMap_;
}

template<class MatrixType>
Teuchos::RCP<const typename AdditiveSchwarzFilter<MatrixType>::map_type>
AdditiveSchwarzFilter<MatrixType>::getColMap() const
{
  return localMap_;
}

template<class MatrixType>
Teuchos::RCP<const typename AdditiveSchwarzFilter<MatrixType>::map_type>
AdditiveSchwarzFilter<MatrixType>::getDomainMap() const
{
  return localMap_;
}

template<class MatrixType>
Teuchos::RCP<const typename AdditiveSchwarzFilter<MatrixType>::map_type>
AdditiveSchwarzFilter<MatrixType>::getRangeMap() const
{
  return localMap_;
}

template<class MatrixType>
Teuchos::RCP<const Tpetra::RowGraph<typename MatrixType::local_ordinal_type,
                                    typename MatrixType::global_ordinal_type,
                                    typename MatrixType::node_type> >
AdditiveSchwarzFilter<MatrixType>::getGraph() const
{
  //NOTE BMK 6-22: this is to maintain compatibilty with LocalFilter.
  //Situations like overlapping AdditiveSchwarz + BlockRelaxation
  //require the importer of the original distributed graph, even though the
  //BlockRelaxation is preconditioning a local matrix (A_).
  return A_unfiltered_->getGraph();
}

template<class MatrixType>
typename MatrixType::local_ordinal_type AdditiveSchwarzFilter<MatrixType>::getBlockSize() const
{
  return A_->getBlockSize();
}

template<class MatrixType>
global_size_t AdditiveSchwarzFilter<MatrixType>::getGlobalNumRows() const
{
  return A_->getGlobalNumRows();
}

template<class MatrixType>
global_size_t AdditiveSchwarzFilter<MatrixType>::getGlobalNumCols() const
{
  return A_->getGlobalNumCols();
}

template<class MatrixType>
size_t AdditiveSchwarzFilter<MatrixType>::getLocalNumRows() const
{
  return A_->getLocalNumRows();
}

template<class MatrixType>
size_t AdditiveSchwarzFilter<MatrixType>::getLocalNumCols() const
{
  return A_->getLocalNumCols();
}

template<class MatrixType>
typename MatrixType::global_ordinal_type AdditiveSchwarzFilter<MatrixType>::getIndexBase() const
{
  return A_->getIndexBase();
}

template<class MatrixType>
global_size_t AdditiveSchwarzFilter<MatrixType>::getGlobalNumEntries() const
{
  return getLocalNumEntries();
}

template<class MatrixType>
size_t AdditiveSchwarzFilter<MatrixType>::getLocalNumEntries() const
{
  return A_->getLocalNumEntries();
}

template<class MatrixType>
size_t AdditiveSchwarzFilter<MatrixType>::
getNumEntriesInGlobalRow (global_ordinal_type globalRow) const
{
  return A_->getNumEntriesInGlobalRow(globalRow);
}

template<class MatrixType>
size_t AdditiveSchwarzFilter<MatrixType>::
getNumEntriesInLocalRow (local_ordinal_type localRow) const
{
  return A_->getNumEntriesInLocalRow(localRow);
}

template<class MatrixType>
size_t AdditiveSchwarzFilter<MatrixType>::getGlobalMaxNumRowEntries() const
{
  //Use A_unfiltered_ to get a valid upper bound for this.
  //This lets us support this function without computing global constants on A_.
  return A_unfiltered_->getGlobalMaxNumRowEntries();
}

template<class MatrixType>
size_t AdditiveSchwarzFilter<MatrixType>::getLocalMaxNumRowEntries() const
{
  //Use A_unfiltered_ to get a valid upper bound for this
  return A_unfiltered_->getLocalMaxNumRowEntries();
}


template<class MatrixType>
bool AdditiveSchwarzFilter<MatrixType>::hasColMap() const
{
  return true;
}


template<class MatrixType>
bool AdditiveSchwarzFilter<MatrixType>::isLocallyIndexed() const
{
  return true;
}


template<class MatrixType>
bool AdditiveSchwarzFilter<MatrixType>::isGloballyIndexed() const
{
  return false;
}


template<class MatrixType>
bool AdditiveSchwarzFilter<MatrixType>::isFillComplete() const
{
  return true;
}


template<class MatrixType>
void AdditiveSchwarzFilter<MatrixType>::
 getGlobalRowCopy (global_ordinal_type globalRow,
                   nonconst_global_inds_host_view_type &globalInd,
                   nonconst_values_host_view_type &val,
                   size_t& numEntries) const
{
  throw std::runtime_error("Ifpack2::Details::AdditiveSchwarzFilter does not implement getGlobalRowCopy.");
}

template<class MatrixType>
void AdditiveSchwarzFilter<MatrixType>::
getLocalRowCopy (local_ordinal_type LocalRow,
    nonconst_local_inds_host_view_type &Indices,
    nonconst_values_host_view_type &Values,
    size_t& NumEntries) const

{
  A_->getLocalRowCopy(LocalRow, Indices, Values, NumEntries);
}

template<class MatrixType>
void AdditiveSchwarzFilter<MatrixType>::getGlobalRowView(global_ordinal_type /* GlobalRow */,
                                                  global_inds_host_view_type &/*indices*/,
                                                  values_host_view_type &/*values*/) const
{
  throw std::runtime_error("Ifpack2::AdditiveSchwarzFilter: does not support getGlobalRowView.");
}

template<class MatrixType>
void AdditiveSchwarzFilter<MatrixType>::getLocalRowView(local_ordinal_type LocalRow,
    local_inds_host_view_type & indices,
    values_host_view_type & values) const
{
  A_->getLocalRowView(LocalRow, indices, values);
}

template<class MatrixType>
void AdditiveSchwarzFilter<MatrixType>::
getLocalDiagCopy (Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &diag) const
{
  // This is somewhat dubious as to how the maps match.
  A_->getLocalDiagCopy(diag);
}

template<class MatrixType>
void AdditiveSchwarzFilter<MatrixType>::leftScale(const Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& /* x */)
{
  throw std::runtime_error("Ifpack2::AdditiveSchwarzFilter does not support leftScale.");
}

template<class MatrixType>
void AdditiveSchwarzFilter<MatrixType>::rightScale(const Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& /* x */)
{
  throw std::runtime_error("Ifpack2::AdditiveSchwarzFilter does not support rightScale.");
}

template<class MatrixType>
void AdditiveSchwarzFilter<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  A_->apply(X, Y, mode, alpha, beta);
}

template<class MatrixType>
void AdditiveSchwarzFilter<MatrixType>::
CreateReducedProblem(
    const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& OverlappingB,
    Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& OverlappingY,
    Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& ReducedReorderedB) const
{
  //NOTE: the three vectors here are always constructed within AdditiveSchwarz and are not subviews,
  //so they are all constant stride (this avoids a lot of boilerplate with whichVectors)
  TEUCHOS_TEST_FOR_EXCEPTION(!OverlappingB.isConstantStride() || !OverlappingY.isConstantStride() || !ReducedReorderedB.isConstantStride(),
      std::logic_error, "Ifpack2::AdditiveSchwarzFilter::CreateReducedProblem ERROR: One of the input MultiVectors is not constant stride.");
  local_ordinal_type numVecs = OverlappingB.getNumVectors();
  auto b = OverlappingB.getLocalViewDevice(Tpetra::Access::ReadOnly);
  auto reducedB = ReducedReorderedB.getLocalViewDevice(Tpetra::Access::OverwriteAll);
  auto singletons = singletons_.view_device();
  auto singletonDiags = singletonDiagonals_.view_device();
  auto revperm = reverseperm_.view_device();
  //First, solve the singletons.
  {
    auto y = OverlappingY.getLocalViewDevice(Tpetra::Access::ReadWrite);
    Kokkos::parallel_for(policy_2d_type({0, 0}, {(local_ordinal_type) numSingletons_, numVecs}),
      KOKKOS_LAMBDA(local_ordinal_type singletonIndex, local_ordinal_type col)
      {
        local_ordinal_type row = singletons(singletonIndex);
        y(row, col) = b(row, col) / singletonDiags(singletonIndex);
      });
  }
  //Next, permute OverlappingB to ReducedReorderedB.
  Kokkos::parallel_for(policy_2d_type({0, 0}, {(local_ordinal_type) reducedB.extent(0), numVecs}),
    KOKKOS_LAMBDA(local_ordinal_type row, local_ordinal_type col)
    {
      reducedB(row, col) = b(revperm(row), col);
    });
  //Finally, if there are singletons, eliminate the matrix entries which are in singleton columns ("eliminate" here means update reducedB like in row reduction)
  //This could be done more efficiently by storing a separate matrix of just the singleton columns and permuted non-singleton rows, but this adds a lot of complexity.
  //Instead, just apply() the unfiltered matrix to OverlappingY (which is 0, except for singletons), and then permute the result of that and subtract it from reducedB.
  if(numSingletons_)
  {
    mv_type singletonUpdates(A_unfiltered_->getRowMap(), numVecs, false);
    A_unfiltered_->apply(OverlappingY, singletonUpdates);
    auto updatesView = singletonUpdates.getLocalViewDevice(Tpetra::Access::ReadOnly);
    // auto revperm = reverseperm_.view_device();
    Kokkos::parallel_for(policy_2d_type({0, 0}, {(local_ordinal_type) reducedB.extent(0), numVecs}),
      KOKKOS_LAMBDA(local_ordinal_type row, local_ordinal_type col)
      {
        local_ordinal_type origRow = revperm(row);
        reducedB(row, col) -= updatesView(origRow, col);
      });
  }
}

template<class MatrixType>
void AdditiveSchwarzFilter<MatrixType>::UpdateLHS(
    const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& ReducedReorderedY,
    Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& OverlappingY) const
{
  //NOTE: both vectors here are always constructed within AdditiveSchwarz and are not subviews,
  //so they are all constant stride (this avoids a lot of boilerplate with whichVectors)
  TEUCHOS_TEST_FOR_EXCEPTION(!ReducedReorderedY.isConstantStride() || !OverlappingY.isConstantStride(),
      std::logic_error, "Ifpack2::AdditiveSchwarzFilter::UpdateLHS ERROR: One of the input MultiVectors is not constant stride.");
  auto reducedY = ReducedReorderedY.getLocalViewDevice(Tpetra::Access::ReadOnly);
  auto y = OverlappingY.getLocalViewDevice(Tpetra::Access::ReadWrite);
  auto revperm = reverseperm_.view_device();
  Kokkos::parallel_for(policy_2d_type({0, 0}, {(local_ordinal_type) reducedY.extent(0), (local_ordinal_type) reducedY.extent(1)}),
    KOKKOS_LAMBDA(local_ordinal_type row, local_ordinal_type col)
    {
      y(revperm(row), col) = reducedY(row, col);
    });
}

template<class MatrixType>
bool AdditiveSchwarzFilter<MatrixType>::hasTransposeApply() const
{
  return true;
}


template<class MatrixType>
bool AdditiveSchwarzFilter<MatrixType>::supportsRowViews() const
{
  return true;
}

template<class MatrixType>
typename AdditiveSchwarzFilter<MatrixType>::mag_type AdditiveSchwarzFilter<MatrixType>::getFrobeniusNorm() const
{
  // Reordering doesn't change the Frobenius norm.
  return A_->getFrobeniusNorm ();
}

template<class MatrixType>
bool
AdditiveSchwarzFilter<MatrixType>::
mapPairsAreFitted (const row_matrix_type& A)
{
  const map_type& rangeMap = * (A.getRangeMap ());
  const map_type& rowMap = * (A.getRowMap ());
  const bool rangeAndRowFitted = mapPairIsFitted (rowMap, rangeMap);

  const map_type& domainMap = * (A.getDomainMap ());
  const map_type& columnMap = * (A.getColMap ());
  const bool domainAndColumnFitted = mapPairIsFitted (columnMap, domainMap);

  //Note BMK 6-22: Map::isLocallyFitted is a local-only operation, not a collective.
  //This means that it can return different values on different ranks. This can cause MPI to hang,
  //even though it's supposed to terminate globally when any single rank does.
  //
  //This function doesn't need to be fast since it's debug-only code.
  int localSuccess = rangeAndRowFitted && domainAndColumnFitted;
  int globalSuccess;

  Teuchos::reduceAll<int, int> (*(A.getComm()), Teuchos::REDUCE_MIN, localSuccess, Teuchos::outArg (globalSuccess));

  return globalSuccess == 1;
}


template<class MatrixType>
bool
AdditiveSchwarzFilter<MatrixType>::
mapPairIsFitted (const map_type& map1, const map_type& map2)
{
  return map1.isLocallyFitted (map2);
}


}} // namespace Ifpack2::Details

#define IFPACK2_DETAILS_ADDITIVESCHWARZFILTER_INSTANT(S,LO,GO,N)                        \
  template class Ifpack2::Details::AdditiveSchwarzFilter< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif


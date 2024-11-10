// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// clang-format off

#ifndef TPETRA_BLOCKCRSMATRIX_DEF_HPP
#define TPETRA_BLOCKCRSMATRIX_DEF_HPP

/// \file Tpetra_BlockCrsMatrix_def.hpp
/// \brief Definition of Tpetra::BlockCrsMatrix

#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_PackTraits.hpp"
#include "Tpetra_Details_Profiling.hpp"
#include "Tpetra_BlockMultiVector.hpp"
#include "Tpetra_BlockView.hpp"

#include "Teuchos_TimeMonitor.hpp"
#ifdef HAVE_TPETRA_DEBUG
#  include <set>
#endif // HAVE_TPETRA_DEBUG

#include "KokkosSparse_spmv.hpp"

namespace Tpetra {

namespace Impl {

  template<typename T>
  struct BlockCrsRowStruct {
    T totalNumEntries, totalNumBytes, maxRowLength;
    KOKKOS_DEFAULTED_FUNCTION BlockCrsRowStruct() = default;
    KOKKOS_DEFAULTED_FUNCTION BlockCrsRowStruct(const BlockCrsRowStruct &b) = default;
    KOKKOS_INLINE_FUNCTION BlockCrsRowStruct(const T& numEnt, const T& numBytes, const T& rowLength)
      : totalNumEntries(numEnt), totalNumBytes(numBytes), maxRowLength(rowLength) {}
  };

  template<typename T>
  static
  KOKKOS_INLINE_FUNCTION
  void operator+=(BlockCrsRowStruct<T> &a, const BlockCrsRowStruct<T> &b) {
    a.totalNumEntries += b.totalNumEntries;
    a.totalNumBytes   += b.totalNumBytes;
    a.maxRowLength     = a.maxRowLength > b.maxRowLength ? a.maxRowLength : b.maxRowLength;
  }

  template<typename T, typename ExecSpace>
  struct BlockCrsReducer {
    typedef BlockCrsReducer reducer;
    typedef T value_type;
    typedef Kokkos::View<value_type,ExecSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;
    value_type *value;

    KOKKOS_INLINE_FUNCTION
    BlockCrsReducer(value_type &val) : value(&val) {}

    KOKKOS_INLINE_FUNCTION void join(value_type &dst, value_type &src) const { dst += src; }
    KOKKOS_INLINE_FUNCTION void join(value_type &dst, const value_type &src) const { dst += src; }
    KOKKOS_INLINE_FUNCTION void init(value_type &val) const { val = value_type(); }
    KOKKOS_INLINE_FUNCTION value_type& reference() { return *value; }
    KOKKOS_INLINE_FUNCTION result_view_type view() const { return result_view_type(value); }
  };

} // namespace Impl

namespace { // (anonymous)

// Implementation of BlockCrsMatrix::getLocalDiagCopy (non-deprecated
// version that takes two Kokkos::View arguments).
template<class Scalar, class LO, class GO, class Node>
class GetLocalDiagCopy {
public:
  typedef typename Node::device_type device_type;
  typedef size_t diag_offset_type;
  typedef Kokkos::View<const size_t*, device_type,
                       Kokkos::MemoryUnmanaged> diag_offsets_type;
  typedef typename ::Tpetra::CrsGraph<LO, GO, Node> global_graph_type;
  typedef typename global_graph_type::local_graph_device_type local_graph_device_type;
  typedef typename local_graph_device_type::row_map_type row_offsets_type;
  typedef typename ::Tpetra::BlockMultiVector<Scalar, LO, GO, Node>::impl_scalar_type IST;
  typedef Kokkos::View<IST***, device_type, Kokkos::MemoryUnmanaged> diag_type;
  typedef Kokkos::View<const IST*, device_type, Kokkos::MemoryUnmanaged> values_type;

  // Constructor
  GetLocalDiagCopy (const diag_type& diag,
                    const values_type& val,
                    const diag_offsets_type& diagOffsets,
                    const row_offsets_type& ptr,
                    const LO blockSize) :
    diag_ (diag),
    diagOffsets_ (diagOffsets),
    ptr_ (ptr),
    blockSize_ (blockSize),
    offsetPerBlock_ (blockSize_*blockSize_),
    val_(val)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator() (const LO& lclRowInd) const
  {
    using Kokkos::ALL;

    // Get row offset
    const size_t absOffset = ptr_[lclRowInd];

    // Get offset relative to start of row
    const size_t relOffset = diagOffsets_[lclRowInd];

    // Get the total offset
    const size_t pointOffset = (absOffset+relOffset)*offsetPerBlock_;

    // Get a view of the block.  BCRS currently uses LayoutRight
    // regardless of the device.
    typedef Kokkos::View<const IST**, Impl::BlockCrsMatrixLittleBlockArrayLayout,
      device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      const_little_block_type;
    const_little_block_type D_in (val_.data () + pointOffset,
                                  blockSize_, blockSize_);
    auto D_out = Kokkos::subview (diag_, lclRowInd, ALL (), ALL ());
    COPY (D_in, D_out);
  }

  private:
    diag_type diag_;
    diag_offsets_type diagOffsets_;
    row_offsets_type ptr_;
    LO blockSize_;
    LO offsetPerBlock_;
    values_type val_;
  };
} // namespace (anonymous)

  template<class Scalar, class LO, class GO, class Node>
  std::ostream&
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  markLocalErrorAndGetStream ()
  {
    * (this->localError_) = true;
    if ((*errs_).is_null ()) {
      *errs_ = Teuchos::rcp (new std::ostringstream ());
    }
    return **errs_;
  }

  template<class Scalar, class LO, class GO, class Node>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  BlockCrsMatrix () :
    dist_object_type (Teuchos::rcp (new map_type ())), // nonnull, so DistObject doesn't throw
    graph_ (Teuchos::rcp (new map_type ()), 0), // FIXME (mfh 16 May 2014) no empty ctor yet
    blockSize_ (static_cast<LO> (0)),
    X_colMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    Y_rowMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    pointImporter_ (new Teuchos::RCP<typename crs_graph_type::import_type> ()),
    offsetPerBlock_ (0),
    localError_ (new bool (false)),
    errs_ (new Teuchos::RCP<std::ostringstream> ()) // ptr to a null ptr
  {
  }

  template<class Scalar, class LO, class GO, class Node>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  BlockCrsMatrix (const crs_graph_type& graph,
                  const LO blockSize) :
    dist_object_type (graph.getMap ()),
    graph_ (graph),
    rowMeshMap_ (* (graph.getRowMap ())),
    blockSize_ (blockSize),
    X_colMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    Y_rowMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    pointImporter_ (new Teuchos::RCP<typename crs_graph_type::import_type> ()),
    offsetPerBlock_ (blockSize * blockSize),
    localError_ (new bool (false)),
    errs_ (new Teuchos::RCP<std::ostringstream> ()) // ptr to a null ptr
  {

    /// KK : additional check is needed that graph is fill complete.
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! graph_.isSorted (), std::invalid_argument, "Tpetra::"
      "BlockCrsMatrix constructor: The input CrsGraph does not have sorted "
      "rows (isSorted() is false).  This class assumes sorted rows.");

    graphRCP_ = Teuchos::rcpFromRef(graph_);

    // Trick to test whether LO is nonpositive, without a compiler
    // warning in case LO is unsigned (which is generally a bad idea
    // anyway).  I don't promise that the trick works, but it
    // generally does with gcc at least, in my experience.
    const bool blockSizeIsNonpositive = (blockSize + 1 <= 1);
    TEUCHOS_TEST_FOR_EXCEPTION(
      blockSizeIsNonpositive, std::invalid_argument, "Tpetra::"
      "BlockCrsMatrix constructor: The input blockSize = " << blockSize <<
      " <= 0.  The block size must be positive.");

    domainPointMap_ = BMV::makePointMap (* (graph.getDomainMap ()), blockSize);
    rangePointMap_ = BMV::makePointMap (* (graph.getRangeMap ()), blockSize);

    {
      // These are rcp
      const auto domainPointMap = getDomainMap();
      const auto colPointMap = Teuchos::rcp 
        (new typename BMV::map_type (BMV::makePointMap (*graph_.getColMap(), blockSize_)));
      *pointImporter_ = Teuchos::rcp 
        (new typename crs_graph_type::import_type (domainPointMap, colPointMap));
    }
    {
      auto local_graph_h = graph.getLocalGraphHost ();
      auto ptr_h = local_graph_h.row_map;
      ptrHost_ = decltype(ptrHost_)(Kokkos::ViewAllocateWithoutInitializing("graph row offset"), ptr_h.extent(0));
      Kokkos::deep_copy(ptrHost_, ptr_h);

      auto ind_h = local_graph_h.entries;
      indHost_ = decltype(indHost_)(Kokkos::ViewAllocateWithoutInitializing("graph column indices"), ind_h.extent(0));
      // DEEP_COPY REVIEW - HOST-TO-HOST
      Kokkos::deep_copy (indHost_, ind_h);

      const auto numValEnt = ind_h.extent(0) * offsetPerBlock ();
      val_ = decltype (val_) (impl_scalar_type_dualview("val", numValEnt));
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  BlockCrsMatrix (const crs_graph_type& graph,
                  const typename local_matrix_device_type::values_type& values,
                  const LO blockSize) :
    dist_object_type (graph.getMap ()),
    graph_ (graph),
    rowMeshMap_ (* (graph.getRowMap ())),
    blockSize_ (blockSize),
    X_colMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    Y_rowMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    pointImporter_ (new Teuchos::RCP<typename crs_graph_type::import_type> ()),
    offsetPerBlock_ (blockSize * blockSize),
    localError_ (new bool (false)),
    errs_ (new Teuchos::RCP<std::ostringstream> ()) // ptr to a null ptr
  {
    /// KK : additional check is needed that graph is fill complete.
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! graph_.isSorted (), std::invalid_argument, "Tpetra::"
      "BlockCrsMatrix constructor: The input CrsGraph does not have sorted "
      "rows (isSorted() is false).  This class assumes sorted rows.");

    graphRCP_ = Teuchos::rcpFromRef(graph_);

    // Trick to test whether LO is nonpositive, without a compiler
    // warning in case LO is unsigned (which is generally a bad idea
    // anyway).  I don't promise that the trick works, but it
    // generally does with gcc at least, in my experience.
    const bool blockSizeIsNonpositive = (blockSize + 1 <= 1);
    TEUCHOS_TEST_FOR_EXCEPTION(
      blockSizeIsNonpositive, std::invalid_argument, "Tpetra::"
      "BlockCrsMatrix constructor: The input blockSize = " << blockSize <<
      " <= 0.  The block size must be positive.");

    domainPointMap_ = BMV::makePointMap (* (graph.getDomainMap ()), blockSize);
    rangePointMap_ = BMV::makePointMap (* (graph.getRangeMap ()), blockSize);

    {
      // These are rcp
      const auto domainPointMap = getDomainMap();
      const auto colPointMap = Teuchos::rcp
        (new typename BMV::map_type (BMV::makePointMap (*graph_.getColMap(), blockSize_)));
      *pointImporter_ = Teuchos::rcp
        (new typename crs_graph_type::import_type (domainPointMap, colPointMap));
    }
    {
      auto local_graph_h = graph.getLocalGraphHost ();
      auto ptr_h = local_graph_h.row_map;
      ptrHost_ = decltype(ptrHost_)(Kokkos::ViewAllocateWithoutInitializing("graph row offset"), ptr_h.extent(0));
      Kokkos::deep_copy(ptrHost_, ptr_h);

      auto ind_h = local_graph_h.entries;
      indHost_ = decltype(indHost_)(Kokkos::ViewAllocateWithoutInitializing("graph column indices"), ind_h.extent(0));
      Kokkos::deep_copy (indHost_, ind_h);

      const auto numValEnt = ind_h.extent(0) * offsetPerBlock ();
      TEUCHOS_ASSERT_EQUALITY(numValEnt, values.size());
      val_ = decltype (val_) (values);
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  BlockCrsMatrix (const crs_graph_type& graph,
                  const map_type& domainPointMap,
                  const map_type& rangePointMap,
                  const LO blockSize) :
    dist_object_type (graph.getMap ()),
    graph_ (graph),
    rowMeshMap_ (* (graph.getRowMap ())),
    domainPointMap_ (domainPointMap),
    rangePointMap_ (rangePointMap),
    blockSize_ (blockSize),
    X_colMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    Y_rowMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    pointImporter_ (new Teuchos::RCP<typename crs_graph_type::import_type> ()),
    offsetPerBlock_ (blockSize * blockSize),
    localError_ (new bool (false)),
    errs_ (new Teuchos::RCP<std::ostringstream> ()) // ptr to a null ptr
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! graph_.isSorted (), std::invalid_argument, "Tpetra::"
      "BlockCrsMatrix constructor: The input CrsGraph does not have sorted "
      "rows (isSorted() is false).  This class assumes sorted rows.");

    graphRCP_ = Teuchos::rcpFromRef(graph_);

    // Trick to test whether LO is nonpositive, without a compiler
    // warning in case LO is unsigned (which is generally a bad idea
    // anyway).  I don't promise that the trick works, but it
    // generally does with gcc at least, in my experience.
    const bool blockSizeIsNonpositive = (blockSize + 1 <= 1);
    TEUCHOS_TEST_FOR_EXCEPTION(
      blockSizeIsNonpositive, std::invalid_argument, "Tpetra::"
      "BlockCrsMatrix constructor: The input blockSize = " << blockSize <<
      " <= 0.  The block size must be positive.");
    {
      // These are rcp
      const auto rcpDomainPointMap = getDomainMap();
      const auto colPointMap = Teuchos::rcp 
        (new typename BMV::map_type (BMV::makePointMap (*graph_.getColMap(), blockSize_)));
      *pointImporter_ = Teuchos::rcp 
        (new typename crs_graph_type::import_type (rcpDomainPointMap, colPointMap));
    }
    {
      auto local_graph_h = graph.getLocalGraphHost ();
      auto ptr_h = local_graph_h.row_map;
      ptrHost_ = decltype(ptrHost_)(Kokkos::ViewAllocateWithoutInitializing("graph row offset"), ptr_h.extent(0));
      // DEEP_COPY REVIEW - HOST-TO-HOST
      Kokkos::deep_copy(ptrHost_, ptr_h);

      auto ind_h = local_graph_h.entries;
      indHost_ = decltype(indHost_)(Kokkos::ViewAllocateWithoutInitializing("graph column indices"), ind_h.extent(0));
      // DEEP_COPY REVIEW - HOST-TO-HOST
      Kokkos::deep_copy (indHost_, ind_h);

      const auto numValEnt = ind_h.extent(0) * offsetPerBlock ();
      val_ = decltype (val_) (impl_scalar_type_dualview("val", numValEnt));

    }
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<const typename BlockCrsMatrix<Scalar, LO, GO, Node>::map_type>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getDomainMap () const
  { // Copy constructor of map_type does a shallow copy.
    // We're only returning by RCP for backwards compatibility.
    return Teuchos::rcp (new map_type (domainPointMap_));
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<const typename BlockCrsMatrix<Scalar, LO, GO, Node>::map_type>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getRangeMap () const
  { // Copy constructor of map_type does a shallow copy.
    // We're only returning by RCP for backwards compatibility.
    return Teuchos::rcp (new map_type (rangePointMap_));
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<const typename BlockCrsMatrix<Scalar, LO, GO, Node>::map_type>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getRowMap () const
  {
    return graph_.getRowMap();
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<const typename BlockCrsMatrix<Scalar, LO, GO, Node>::map_type>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getColMap () const
  {
    return graph_.getColMap();
  }

  template<class Scalar, class LO, class GO, class Node>
  global_size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalNumRows() const
  {
    return graph_.getGlobalNumRows();
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalNumRows() const
  {
    return graph_.getLocalNumRows();
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalMaxNumRowEntries() const
  {
    return graph_.getLocalMaxNumRowEntries();
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  apply (const mv_type& X,
         mv_type& Y,
         Teuchos::ETransp mode,
         Scalar alpha,
         Scalar beta) const
  {
    using this_BCRS_type = BlockCrsMatrix<Scalar, LO, GO, Node>;
    TEUCHOS_TEST_FOR_EXCEPTION(
      mode != Teuchos::NO_TRANS && mode != Teuchos::TRANS && mode != Teuchos::CONJ_TRANS,
      std::invalid_argument, "Tpetra::BlockCrsMatrix::apply: "
      "Invalid 'mode' argument.  Valid values are Teuchos::NO_TRANS, "
      "Teuchos::TRANS, and Teuchos::CONJ_TRANS.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      !X.isConstantStride() || !Y.isConstantStride(),
      std::invalid_argument, "Tpetra::BlockCrsMatrix::apply: "
      "X and Y must both be constant stride");

    BMV X_view;
    BMV Y_view;
    const LO blockSize = getBlockSize ();
    try {
      X_view = BMV (X, * (graph_.getDomainMap ()), blockSize);
      Y_view = BMV (Y, * (graph_.getRangeMap ()), blockSize);
    }
    catch (std::invalid_argument& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::BlockCrsMatrix::"
        "apply: Either the input MultiVector X or the output MultiVector Y "
        "cannot be viewed as a BlockMultiVector, given this BlockCrsMatrix's "
        "graph.  BlockMultiVector's constructor threw the following exception: "
        << e.what ());
    }

    try {
      // mfh 16 May 2014: Casting 'this' to nonconst is icky, but it's
      // either that or mark fields of this class as 'mutable'.  The
      // problem is that applyBlock wants to do lazy initialization of
      // temporary block multivectors.
      const_cast<this_BCRS_type&> (*this).applyBlock (X_view, Y_view, mode, alpha, beta);
    } catch (std::invalid_argument& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::BlockCrsMatrix::"
        "apply: The implementation method applyBlock complained about having "
        "an invalid argument.  It reported the following: " << e.what ());
    } catch (std::logic_error& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::BlockCrsMatrix::"
        "apply: The implementation method applyBlock complained about a "
        "possible bug in its implementation.  It reported the following: "
        << e.what ());
    } catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::BlockCrsMatrix::"
        "apply: The implementation method applyBlock threw an exception which "
        "is neither std::invalid_argument nor std::logic_error, but is a "
        "subclass of std::exception.  It reported the following: "
        << e.what ());
    } catch (...) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error, "Tpetra::BlockCrsMatrix::"
        "apply: The implementation method applyBlock threw an exception which "
        "is not an instance of a subclass of std::exception.  This probably "
        "indicates a bug.  Please report this to the Tpetra developers.");
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  applyBlock (const BlockMultiVector<Scalar, LO, GO, Node>& X,
              BlockMultiVector<Scalar, LO, GO, Node>& Y,
              Teuchos::ETransp mode,
              const Scalar alpha,
              const Scalar beta)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      X.getBlockSize () != Y.getBlockSize (), std::invalid_argument,
      "Tpetra::BlockCrsMatrix::applyBlock: "
      "X and Y have different block sizes.  X.getBlockSize() = "
      << X.getBlockSize () << " != Y.getBlockSize() = "
      << Y.getBlockSize () << ".");

    if (mode == Teuchos::NO_TRANS) {
      applyBlockNoTrans (X, Y, alpha, beta);
    } else if (mode == Teuchos::TRANS || mode == Teuchos::CONJ_TRANS) {
      applyBlockTrans (X, Y, mode, alpha, beta);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::BlockCrsMatrix::"
        "applyBlock: Invalid 'mode' argument.  Valid values are "
        "Teuchos::NO_TRANS, Teuchos::TRANS, and Teuchos::CONJ_TRANS.");
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  importAndFillComplete (Teuchos::RCP<BlockCrsMatrix<Scalar, LO, GO, Node> >& destMatrix,
                         const Import<LO, GO, Node>& importer) const
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using this_BCRS_type = BlockCrsMatrix<Scalar, LO, GO, Node>;

    // Right now, we make many assumptions...
    TEUCHOS_TEST_FOR_EXCEPTION(!destMatrix.is_null(), std::invalid_argument,
                               "destMatrix is required to be null.");
 
    // BlockCrsMatrix requires a complete graph at construction.
    // So first step is to import and fill complete the destGraph.
    RCP<crs_graph_type>  srcGraph = rcp (new  crs_graph_type(this->getCrsGraph()));
    RCP<crs_graph_type> destGraph = importAndFillCompleteCrsGraph<crs_graph_type>(srcGraph, importer,
                                                                                  srcGraph->getDomainMap(),
                                                                                  srcGraph->getRangeMap());


    // Final step, create and import the destMatrix.
    destMatrix = rcp (new this_BCRS_type (*destGraph, getBlockSize()));
    destMatrix->doImport(*this, importer, Tpetra::INSERT);
  }


  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  exportAndFillComplete (Teuchos::RCP<BlockCrsMatrix<Scalar, LO, GO, Node> >& destMatrix,
                         const Export<LO, GO, Node>& exporter) const
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using this_BCRS_type = BlockCrsMatrix<Scalar, LO, GO, Node>;

    // Right now, we make many assumptions...
    TEUCHOS_TEST_FOR_EXCEPTION(!destMatrix.is_null(), std::invalid_argument,
                               "destMatrix is required to be null.");
 
    // BlockCrsMatrix requires a complete graph at construction.
    // So first step is to import and fill complete the destGraph.
    RCP<crs_graph_type>  srcGraph = rcp (new  crs_graph_type(this->getCrsGraph()));
    RCP<crs_graph_type> destGraph = exportAndFillCompleteCrsGraph<crs_graph_type>(srcGraph, exporter,
                                                                                  srcGraph->getDomainMap(),
                                                                                  srcGraph->getRangeMap());


    // Final step, create and import the destMatrix.
    destMatrix = rcp (new this_BCRS_type (*destGraph, getBlockSize()));
    destMatrix->doExport(*this, exporter, Tpetra::INSERT);
  }


  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  setAllToScalar (const Scalar& alpha)
  {
    auto val_d = val_.getDeviceView(Access::OverwriteAll);
    Kokkos::deep_copy(val_d, alpha);
  }

  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  replaceLocalValues (const LO localRowInd,
                      const LO colInds[],
                      const Scalar vals[],
                      const LO numColInds) const
  {
    std::vector<ptrdiff_t> offsets(numColInds);
    const LO numOffsets = this->getLocalRowOffsets(localRowInd, offsets.data(), colInds, numColInds);
    const LO validCount = this->replaceLocalValuesByOffsets(localRowInd, offsets.data(), vals, numOffsets);
    return validCount;
  }

  template <class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar,LO,GO,Node>::
  getLocalDiagOffsets (const Kokkos::View<size_t*, device_type,
                         Kokkos::MemoryUnmanaged>& offsets) const
  {
    graph_.getLocalDiagOffsets (offsets);
  }

  template <class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar,LO,GO,Node>::
  getLocalDiagCopy (const Kokkos::View<impl_scalar_type***, device_type,
                                       Kokkos::MemoryUnmanaged>& diag,
                    const Kokkos::View<const size_t*, device_type,
                                       Kokkos::MemoryUnmanaged>& offsets) const
  {
    using Kokkos::parallel_for;
    const char prefix[] = "Tpetra::BlockCrsMatrix::getLocalDiagCopy (2-arg): ";

    const LO lclNumMeshRows = static_cast<LO> (rowMeshMap_.getLocalNumElements ());
    const LO blockSize = this->getBlockSize ();
    TEUCHOS_TEST_FOR_EXCEPTION
      (static_cast<LO> (diag.extent (0)) < lclNumMeshRows ||
       static_cast<LO> (diag.extent (1)) < blockSize ||
       static_cast<LO> (diag.extent (2)) < blockSize,
       std::invalid_argument, prefix <<
       "The input Kokkos::View is not big enough to hold all the data.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (static_cast<LO> (offsets.size ()) < lclNumMeshRows, std::invalid_argument,
       prefix << "offsets.size() = " << offsets.size () << " < local number of "
       "diagonal blocks " << lclNumMeshRows << ".");

    typedef Kokkos::RangePolicy<execution_space, LO> policy_type;
    typedef GetLocalDiagCopy<Scalar, LO, GO, Node> functor_type;

    // FIXME (mfh 26 May 2016) Not really OK to const_cast here, since
    // we reserve the right to do lazy allocation of device data.  (We
    // don't plan to do lazy allocation for host data; the host
    // version of the data always exists.)
    auto val_d = val_.getDeviceView(Access::ReadOnly);
    parallel_for (policy_type (0, lclNumMeshRows),
                  functor_type (diag, val_d, offsets,
                                graph_.getLocalGraphDevice ().row_map, blockSize_));
  }

  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  absMaxLocalValues (const LO localRowInd,
                     const LO colInds[],
                     const Scalar vals[],
                     const LO numColInds) const
  {
    std::vector<ptrdiff_t> offsets(numColInds);
    const LO numOffsets = this->getLocalRowOffsets(localRowInd, offsets.data(), colInds, numColInds);
    const LO validCount = this->absMaxLocalValuesByOffsets(localRowInd, offsets.data(), vals, numOffsets);
    return validCount;
  }


  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  sumIntoLocalValues (const LO localRowInd,
                      const LO colInds[],
                      const Scalar vals[],
                      const LO numColInds) const
  {
    std::vector<ptrdiff_t> offsets(numColInds);
    const LO numOffsets = this->getLocalRowOffsets(localRowInd, offsets.data(), colInds, numColInds);
    const LO validCount = this->sumIntoLocalValuesByOffsets(localRowInd, offsets.data(), vals, numOffsets);
    return validCount;
  }
  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalRowCopy (LO LocalRow,
                   nonconst_local_inds_host_view_type &Indices,
                   nonconst_values_host_view_type &Values,
                   size_t& NumEntries) const 
  {
    auto vals = getValuesHost(LocalRow);
    const LO numInds = ptrHost_(LocalRow+1) - ptrHost_(LocalRow);
    if (numInds > (LO)Indices.extent(0) || numInds*blockSize_*blockSize_ > (LO)Values.extent(0)) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                  "Tpetra::BlockCrsMatrix::getLocalRowCopy : Column and/or values array is not large enough to hold "
                  << numInds << " row entries");
    }
    const LO * colInds = indHost_.data() + ptrHost_(LocalRow);
    for (LO i=0; i<numInds; ++i) {
      Indices[i] = colInds[i];
    }
    for (LO i=0; i<numInds*blockSize_*blockSize_; ++i) {
      Values[i] = vals[i];
    }
    NumEntries = numInds;
  }

  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalRowOffsets (const LO localRowInd,
                      ptrdiff_t offsets[],
                      const LO colInds[],
                      const LO numColInds) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      // We got no offsets, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }

    const LO LINV = Teuchos::OrdinalTraits<LO>::invalid ();
    LO hint = 0; // Guess for the relative offset into the current row
    LO validCount = 0; // number of valid column indices in colInds

    for (LO k = 0; k < numColInds; ++k) {
      const LO relBlockOffset =
        this->findRelOffsetOfColumnIndex (localRowInd, colInds[k], hint);
      if (relBlockOffset != LINV) {
        offsets[k] = static_cast<ptrdiff_t> (relBlockOffset);
        hint = relBlockOffset + 1;
        ++validCount;
      }
    }
    return validCount;
  }


  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  replaceLocalValuesByOffsets (const LO localRowInd,
                               const ptrdiff_t offsets[],
                               const Scalar vals[],
                               const LO numOffsets) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      // We modified no values, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }
    const impl_scalar_type* const vIn = reinterpret_cast<const impl_scalar_type*> (vals);
    using this_BCRS_type = BlockCrsMatrix<Scalar, LO, GO, Node>;
    auto val_out = const_cast<this_BCRS_type&>(*this).getValuesHostNonConst(localRowInd);
    impl_scalar_type* vOut = val_out.data();

    const size_t perBlockSize = static_cast<LO> (offsetPerBlock ());
    const ptrdiff_t STINV = Teuchos::OrdinalTraits<ptrdiff_t>::invalid ();
    size_t pointOffset = 0; // Current offset into input values
    LO validCount = 0; // number of valid offsets

    for (LO k = 0; k < numOffsets; ++k, pointOffset += perBlockSize) {
      const size_t blockOffset = offsets[k]*perBlockSize;
      if (offsets[k] != STINV) {
        little_block_type A_old = getNonConstLocalBlockFromInput (vOut, blockOffset);
        const_little_block_type A_new = getConstLocalBlockFromInput (vIn, pointOffset);
        COPY (A_new, A_old);
        ++validCount;
      }
    }
    return validCount;
  }


  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  absMaxLocalValuesByOffsets (const LO localRowInd,
                              const ptrdiff_t offsets[],
                              const Scalar vals[],
                              const LO numOffsets) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      // We modified no values, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }
    const impl_scalar_type* const vIn = reinterpret_cast<const impl_scalar_type*> (vals);
    auto val_out = getValuesHost(localRowInd);
    impl_scalar_type* vOut = const_cast<impl_scalar_type*>(val_out.data());

    const size_t perBlockSize = static_cast<LO> (offsetPerBlock ());
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    size_t pointOffset = 0; // Current offset into input values
    LO validCount = 0; // number of valid offsets

    for (LO k = 0; k < numOffsets; ++k, pointOffset += perBlockSize) {
      const size_t blockOffset = offsets[k]*perBlockSize;
      if (blockOffset != STINV) {
        little_block_type A_old = getNonConstLocalBlockFromInput (vOut, blockOffset);
        const_little_block_type A_new = getConstLocalBlockFromInput (vIn, pointOffset);
        ::Tpetra::Impl::absMax (A_old, A_new);
        ++validCount;
      }
    }
    return validCount;
  }


  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  sumIntoLocalValuesByOffsets (const LO localRowInd,
                               const ptrdiff_t offsets[],
                               const Scalar vals[],
                               const LO numOffsets) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      // We modified no values, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }
    const impl_scalar_type ONE = static_cast<impl_scalar_type> (1.0);
    const impl_scalar_type* const vIn = reinterpret_cast<const impl_scalar_type*> (vals);
    typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_BCRS_type;
    auto val_out = const_cast<this_BCRS_type&>(*this).getValuesHostNonConst(localRowInd);
    impl_scalar_type* vOut = val_out.data();

    const size_t perBlockSize = static_cast<LO> (offsetPerBlock ());
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    size_t pointOffset = 0; // Current offset into input values
    LO validCount = 0; // number of valid offsets

    for (LO k = 0; k < numOffsets; ++k, pointOffset += perBlockSize) {
      const size_t blockOffset = offsets[k]*perBlockSize;
      if (blockOffset != STINV) {
        little_block_type A_old = getNonConstLocalBlockFromInput (vOut, blockOffset);
        const_little_block_type A_new = getConstLocalBlockFromInput (vIn, pointOffset);
        AXPY (ONE, A_new, A_old);
        ++validCount;
      }
    }
    return validCount;
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::
  impl_scalar_type_dualview::t_host::const_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getValuesHost () const 
  {
    return val_.getHostView(Access::ReadOnly);
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::
  impl_scalar_type_dualview::t_dev::const_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getValuesDevice () const 
  {
    return val_.getDeviceView(Access::ReadOnly);
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::
  impl_scalar_type_dualview::t_host
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getValuesHostNonConst () const 
  {
    return val_.getHostView(Access::ReadWrite);
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::
  impl_scalar_type_dualview::t_dev
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getValuesDeviceNonConst () const 
  {
    return val_.getDeviceView(Access::ReadWrite);
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::
  impl_scalar_type_dualview::t_host::const_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getValuesHost (const LO& lclRow) const 
  {
    const size_t perBlockSize = static_cast<LO> (offsetPerBlock ());
    auto val = val_.getHostView(Access::ReadOnly);
    auto r_val = Kokkos::subview(val, Kokkos::pair<LO,LO>(ptrHost_(lclRow)*perBlockSize, ptrHost_(lclRow+1)*perBlockSize)); 
    return r_val;
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::
  impl_scalar_type_dualview::t_dev::const_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getValuesDevice (const LO& lclRow) const 
  {
    const size_t perBlockSize = static_cast<LO> (offsetPerBlock ());
    auto val = val_.getDeviceView(Access::ReadOnly);
    auto r_val = Kokkos::subview(val, Kokkos::pair<LO,LO>(ptrHost_(lclRow)*perBlockSize, ptrHost_(lclRow+1)*perBlockSize)); 
    return r_val;
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::impl_scalar_type_dualview::t_host
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getValuesHostNonConst (const LO& lclRow) 
  {
    const size_t perBlockSize = static_cast<LO> (offsetPerBlock ());
    auto val = val_.getHostView(Access::ReadWrite);
    auto r_val = Kokkos::subview(val, Kokkos::pair<LO,LO>(ptrHost_(lclRow)*perBlockSize, ptrHost_(lclRow+1)*perBlockSize)); 
    return r_val;
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::impl_scalar_type_dualview::t_dev
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getValuesDeviceNonConst (const LO& lclRow) 
  {
    const size_t perBlockSize = static_cast<LO> (offsetPerBlock ());
    auto val = val_.getDeviceView(Access::ReadWrite);
    auto r_val = Kokkos::subview(val, Kokkos::pair<LO,LO>(ptrHost_(lclRow)*perBlockSize, ptrHost_(lclRow+1)*perBlockSize)); 
    return r_val;
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNumEntriesInLocalRow (const LO localRowInd) const
  {
    const size_t numEntInGraph = graph_.getNumEntriesInLocalRow (localRowInd);
    if (numEntInGraph == Teuchos::OrdinalTraits<size_t>::invalid ()) {
      return static_cast<LO> (0); // the calling process doesn't have that row
    } else {
      return static_cast<LO> (numEntInGraph);
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::local_matrix_device_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalMatrixDevice () const
  {
    auto numCols = this->graph_.getColMap()->getLocalNumElements();
    auto val = val_.getDeviceView(Access::ReadWrite);
    const LO blockSize = this->getBlockSize ();
    const auto graph = this->graph_.getLocalGraphDevice ();

    return local_matrix_device_type("Tpetra::BlockCrsMatrix::lclMatrixDevice",
                                    numCols,
                                    val,
                                    graph,
                                    blockSize);
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  applyBlockTrans (const BlockMultiVector<Scalar, LO, GO, Node>& X,
                   BlockMultiVector<Scalar, LO, GO, Node>& Y,
                   const Teuchos::ETransp mode,
                   const Scalar alpha,
                   const Scalar beta)
  {
    (void) X;
    (void) Y;
    (void) mode;
    (void) alpha;
    (void) beta;

    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::BlockCrsMatrix::apply: "
      "transpose and conjugate transpose modes are not yet implemented.");
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  applyBlockNoTrans (const BlockMultiVector<Scalar, LO, GO, Node>& X,
                     BlockMultiVector<Scalar, LO, GO, Node>& Y,
                     const Scalar alpha,
                     const Scalar beta)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef ::Tpetra::Import<LO, GO, Node> import_type;
    typedef ::Tpetra::Export<LO, GO, Node> export_type;
    const Scalar zero = STS::zero ();
    const Scalar one = STS::one ();
    RCP<const import_type> import = graph_.getImporter ();
    // "export" is a reserved C++ keyword, so we can't use it.
    RCP<const export_type> theExport = graph_.getExporter ();
    const char prefix[] = "Tpetra::BlockCrsMatrix::applyBlockNoTrans: ";

    if (alpha == zero) {
      if (beta == zero) {
        Y.putScalar (zero); // replace Inf or NaN (BLAS rules)
      }
      else if (beta != one) {
        Y.scale (beta);
      }
    }
    else { // alpha != 0
      const BMV* X_colMap = NULL;
      if (import.is_null ()) {
        try {
          X_colMap = &X;
        }
        catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION
            (true, std::logic_error, prefix << "Tpetra::MultiVector::"
             "operator= threw an exception: " << std::endl << e.what ());
        }
      }
      else {
        // X_colMap_ is a pointer to a pointer to BMV.  Ditto for
        // Y_rowMap_ below.  This lets us do lazy initialization
        // correctly with view semantics of BlockCrsMatrix.  All views
        // of this BlockCrsMatrix have the same outer pointer.  That
        // way, we can set the inner pointer in one view, and all
        // other views will see it.
        if ((*X_colMap_).is_null () ||
            (**X_colMap_).getNumVectors () != X.getNumVectors () ||
            (**X_colMap_).getBlockSize () != X.getBlockSize ()) {
          *X_colMap_ = rcp (new BMV (* (graph_.getColMap ()), getBlockSize (),
                                     static_cast<LO> (X.getNumVectors ())));
        }
        (*X_colMap_)->getMultiVectorView().doImport (X.getMultiVectorView (),
                                                     **pointImporter_,
                                                     ::Tpetra::REPLACE);
        try {
          X_colMap = &(**X_colMap_);
        }
        catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION
            (true, std::logic_error, prefix << "Tpetra::MultiVector::"
            "operator= threw an exception: " << std::endl << e.what ());
        }
      }

      BMV* Y_rowMap = NULL;
      if (theExport.is_null ()) {
        Y_rowMap = &Y;
      }
      else if ((*Y_rowMap_).is_null () ||
                 (**Y_rowMap_).getNumVectors () != Y.getNumVectors () ||
                 (**Y_rowMap_).getBlockSize () != Y.getBlockSize ()) {
        *Y_rowMap_ = rcp (new BMV (* (graph_.getRowMap ()), getBlockSize (),
                                   static_cast<LO> (X.getNumVectors ())));
        try {
          Y_rowMap = &(**Y_rowMap_);
        }
        catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION(
            true, std::logic_error, prefix << "Tpetra::MultiVector::"
            "operator= threw an exception: " << std::endl << e.what ());
        }
      }

      try {
        localApplyBlockNoTrans (*X_colMap, *Y_rowMap, alpha, beta);
      }
      catch (std::exception& e) {
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::runtime_error, prefix << "localApplyBlockNoTrans threw "
           "an exception: " << e.what ());
      }
      catch (...) {
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::runtime_error, prefix << "localApplyBlockNoTrans threw "
           "an exception not a subclass of std::exception.");
      }

      if (! theExport.is_null ()) {
        Y.doExport (*Y_rowMap, *theExport, ::Tpetra::REPLACE);
      }
    }
  }

// clang-format on
template <class Scalar, class LO, class GO, class Node>
void BlockCrsMatrix<Scalar, LO, GO, Node>::localApplyBlockNoTrans(
    const BlockMultiVector<Scalar, LO, GO, Node> &X,
    BlockMultiVector<Scalar, LO, GO, Node> &Y, const Scalar alpha,
    const Scalar beta) {
  ::Tpetra::Details::ProfilingRegion profile_region(
      "Tpetra::BlockCrsMatrix::localApplyBlockNoTrans");
  const impl_scalar_type alpha_impl = alpha;
  const auto graph = this->graph_.getLocalGraphDevice();

  mv_type X_mv = X.getMultiVectorView();
  mv_type Y_mv = Y.getMultiVectorView();
  auto X_lcl = X_mv.getLocalViewDevice(Access::ReadOnly);
  auto Y_lcl = Y_mv.getLocalViewDevice(Access::ReadWrite);

#if KOKKOSKERNELS_VERSION >= 40299
  auto A_lcl = getLocalMatrixDevice();
  if(!applyHelper.get()) {
    // The apply helper does not exist, so create it
    applyHelper = std::make_shared<ApplyHelper>(A_lcl.nnz(), A_lcl.graph.row_map);
  }
  if(applyHelper->shouldUseIntRowptrs())
  {
    auto A_lcl_int_rowptrs = applyHelper->getIntRowptrMatrix(A_lcl);
    KokkosSparse::spmv(
        &applyHelper->handle_int, KokkosSparse::NoTranspose,
        alpha_impl, A_lcl_int_rowptrs, X_lcl, beta, Y_lcl);
  }
  else
  {
    KokkosSparse::spmv(
        &applyHelper->handle, KokkosSparse::NoTranspose,
        alpha_impl, A_lcl, X_lcl, beta, Y_lcl);
  }
#else
  auto A_lcl = getLocalMatrixDevice();
  KokkosSparse::spmv(KokkosSparse::NoTranspose, alpha_impl, A_lcl, X_lcl, beta,
                     Y_lcl);
#endif
}
// clang-format off

  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  findRelOffsetOfColumnIndex (const LO localRowIndex,
                              const LO colIndexToFind,
                              const LO hint) const
  {
    const size_t absStartOffset = ptrHost_[localRowIndex];
    const size_t absEndOffset = ptrHost_[localRowIndex+1];
    const LO numEntriesInRow = static_cast<LO> (absEndOffset - absStartOffset);
    // Amortize pointer arithmetic over the search loop.
    const LO* const curInd = indHost_.data () + absStartOffset;

    // If the hint was correct, then the hint is the offset to return.
    if (hint < numEntriesInRow && curInd[hint] == colIndexToFind) {
      // Always return the offset relative to the current row.
      return hint;
    }

    // The hint was wrong, so we must search for the given column
    // index in the column indices for the given row.
    LO relOffset = Teuchos::OrdinalTraits<LO>::invalid ();

    // We require that the graph have sorted rows.  However, binary
    // search only pays if the current row is longer than a certain
    // amount.  We set this to 32, but you might want to tune this.
    const LO maxNumEntriesForLinearSearch = 32;
    if (numEntriesInRow > maxNumEntriesForLinearSearch) {
      // Use binary search.  It would probably be better for us to
      // roll this loop by hand.  If we wrote it right, a smart
      // compiler could perhaps use conditional loads and avoid
      // branches (according to Jed Brown on May 2014).
      const LO* beg = curInd;
      const LO* end = curInd + numEntriesInRow;
      std::pair<const LO*, const LO*> p =
        std::equal_range (beg, end, colIndexToFind);
      if (p.first != p.second) {
        // offset is relative to the current row
        relOffset = static_cast<LO> (p.first - beg);
      }
    }
    else { // use linear search
      for (LO k = 0; k < numEntriesInRow; ++k) {
        if (colIndexToFind == curInd[k]) {
          relOffset = k; // offset is relative to the current row
          break;
        }
      }
    }

    return relOffset;
  }

  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  offsetPerBlock () const
  {
    return offsetPerBlock_;
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::const_little_block_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getConstLocalBlockFromInput (const impl_scalar_type* val,
                               const size_t pointOffset) const
  {
    // Row major blocks
    const LO rowStride = blockSize_;
    return const_little_block_type (val + pointOffset, blockSize_, rowStride);
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::little_block_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNonConstLocalBlockFromInput (impl_scalar_type* val,
                                  const size_t pointOffset) const
  {
    // Row major blocks
    const LO rowStride = blockSize_;
    return little_block_type (val + pointOffset, blockSize_, rowStride);
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::little_block_host_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNonConstLocalBlockFromInputHost (impl_scalar_type* val,
                                  const size_t pointOffset) const
  {
    // Row major blocks
    const LO rowStride = blockSize_;
    const size_t bs2 = blockSize_ * blockSize_;
    return little_block_host_type (val + bs2 * pointOffset, blockSize_, rowStride);
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::little_block_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalBlockDeviceNonConst (const LO localRowInd, const LO localColInd) const
  {
    using this_BCRS_type = BlockCrsMatrix<Scalar, LO, GO, Node>;

    const size_t absRowBlockOffset = ptrHost_[localRowInd];
    const LO relBlockOffset = this->findRelOffsetOfColumnIndex (localRowInd, localColInd);
    if (relBlockOffset != Teuchos::OrdinalTraits<LO>::invalid ()) {
      const size_t absBlockOffset = absRowBlockOffset + relBlockOffset;
      auto vals = const_cast<this_BCRS_type&>(*this).getValuesDeviceNonConst();
      auto r_val = getNonConstLocalBlockFromInput (vals.data(), absBlockOffset);      
      return r_val; 
    }
    else {
      return little_block_type ();
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::little_block_host_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalBlockHostNonConst (const LO localRowInd, const LO localColInd) const
  {
    using this_BCRS_type = BlockCrsMatrix<Scalar, LO, GO, Node>;

    const size_t absRowBlockOffset = ptrHost_[localRowInd];
    const LO relBlockOffset = this->findRelOffsetOfColumnIndex (localRowInd, localColInd);
    if (relBlockOffset != Teuchos::OrdinalTraits<LO>::invalid ()) {
      const size_t absBlockOffset = absRowBlockOffset + relBlockOffset;
      auto vals = const_cast<this_BCRS_type&>(*this).getValuesHostNonConst();
      auto r_val = getNonConstLocalBlockFromInputHost (vals.data(), absBlockOffset);      
      return r_val; 
    }
    else {
      return little_block_host_type ();
    }
  }


  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  checkSizes (const ::Tpetra::SrcDistObject& source)
  {
    using std::endl;
    typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_BCRS_type;
    const this_BCRS_type* src = dynamic_cast<const this_BCRS_type* > (&source);

    if (src == NULL) {
      std::ostream& err = markLocalErrorAndGetStream ();
      err << "checkSizes: The source object of the Import or Export "
        "must be a BlockCrsMatrix with the same template parameters as the "
        "target object." << endl;
    }
    else {
      // Use a string of ifs, not if-elseifs, because we want to know
      // all the errors.
      if (src->getBlockSize () != this->getBlockSize ()) {
        std::ostream& err = markLocalErrorAndGetStream ();
        err << "checkSizes: The source and target objects of the Import or "
            << "Export must have the same block sizes.  The source's block "
            << "size = " << src->getBlockSize () << " != the target's block "
            << "size = " << this->getBlockSize () << "." << endl;
      }
      if (! src->graph_.isFillComplete ()) {
        std::ostream& err = markLocalErrorAndGetStream ();
        err << "checkSizes: The source object of the Import or Export is "
          "not fill complete.  Both source and target objects must be fill "
          "complete." << endl;
      }
      if (! this->graph_.isFillComplete ()) {
        std::ostream& err = markLocalErrorAndGetStream ();
        err << "checkSizes: The target object of the Import or Export is "
          "not fill complete.  Both source and target objects must be fill "
          "complete." << endl;
      }
      if (src->graph_.getColMap ().is_null ()) {
        std::ostream& err = markLocalErrorAndGetStream ();
        err << "checkSizes: The source object of the Import or Export does "
          "not have a column Map.  Both source and target objects must have "
          "column Maps." << endl;
      }
      if (this->graph_.getColMap ().is_null ()) {
        std::ostream& err = markLocalErrorAndGetStream ();
        err << "checkSizes: The target object of the Import or Export does "
          "not have a column Map.  Both source and target objects must have "
          "column Maps." << endl;
      }
    }

    return ! (* (this->localError_));
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  copyAndPermute
  (const ::Tpetra::SrcDistObject& source,
   const size_t numSameIDs,
   const Kokkos::DualView<const local_ordinal_type*,
     buffer_device_type>& permuteToLIDs,
   const Kokkos::DualView<const local_ordinal_type*,
     buffer_device_type>& permuteFromLIDs,
   const CombineMode /*CM*/)
  {
    using ::Tpetra::Details::Behavior;
    using ::Tpetra::Details::dualViewStatusToString;
    using ::Tpetra::Details::ProfilingRegion;
    using std::endl;
    using this_BCRS_type = BlockCrsMatrix<Scalar, LO, GO, Node>;

    ProfilingRegion profile_region("Tpetra::BlockCrsMatrix::copyAndPermute");
    const bool debug = Behavior::debug();
    const bool verbose = Behavior::verbose();

    // Define this function prefix
    std::string prefix;
    {
      std::ostringstream os;
      const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
      os << "Proc " << myRank << ": BlockCrsMatrix::copyAndPermute : " << endl;
      prefix = os.str();
    }

    // check if this already includes a local error
    if (* (this->localError_)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix
          << "The target object of the Import or Export is already in an error state."
          << endl;
      return;
    }

    //
    // Verbose input dual view status
    //
    if (verbose) {
      std::ostringstream os;
      os << prefix << endl
         << prefix << "  " << dualViewStatusToString (permuteToLIDs, "permuteToLIDs") << endl
         << prefix << "  " << dualViewStatusToString (permuteFromLIDs, "permuteFromLIDs") << endl;
      std::cerr << os.str ();
    }

    ///
    /// Check input valid
    ///
    if (permuteToLIDs.extent (0) != permuteFromLIDs.extent (0)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix
          << "permuteToLIDs.extent(0) = " << permuteToLIDs.extent (0)
          << " != permuteFromLIDs.extent(0) = " << permuteFromLIDs.extent(0)
          << "." << endl;
      return;
    }
    if (permuteToLIDs.need_sync_host () || permuteFromLIDs.need_sync_host ()) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix
          << "Both permuteToLIDs and permuteFromLIDs must be sync'd to host."
          << endl;
      return;
    }

    const this_BCRS_type* src = dynamic_cast<const this_BCRS_type* > (&source);
    if (src == NULL) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix
          << "The source (input) object of the Import or "
        "Export is either not a BlockCrsMatrix, or does not have the right "
        "template parameters.  checkSizes() should have caught this.  "
        "Please report this bug to the Tpetra developers." << endl;
      return;
    }

    bool lclErr = false;
#ifdef HAVE_TPETRA_DEBUG
    std::set<LO> invalidSrcCopyRows;
    std::set<LO> invalidDstCopyRows;
    std::set<LO> invalidDstCopyCols;
    std::set<LO> invalidDstPermuteCols;
    std::set<LO> invalidPermuteFromRows;
#endif // HAVE_TPETRA_DEBUG

    // Copy the initial sequence of rows that are the same.
    //
    // The two graphs might have different column Maps, so we need to
    // do this using global column indices.  This is purely local, so
    // we only need to check for local sameness of the two column
    // Maps.

#ifdef HAVE_TPETRA_DEBUG
    const map_type& srcRowMap = * (src->graph_.getRowMap ());
#endif // HAVE_TPETRA_DEBUG
    const map_type& dstRowMap = * (this->graph_.getRowMap ());
    const map_type& srcColMap = * (src->graph_.getColMap ());
    const map_type& dstColMap = * (this->graph_.getColMap ());
    const bool canUseLocalColumnIndices = srcColMap.locallySameAs (dstColMap);

    const size_t numPermute = static_cast<size_t> (permuteFromLIDs.extent(0));
    if (verbose) {
      std::ostringstream os;
      os << prefix
         << "canUseLocalColumnIndices: "
         << (canUseLocalColumnIndices ? "true" : "false")
         << ", numPermute: " << numPermute
         << endl;
      std::cerr << os.str ();
    }

    const auto permuteToLIDsHost = permuteToLIDs.view_host();
    const auto permuteFromLIDsHost = permuteFromLIDs.view_host();

    if (canUseLocalColumnIndices) {
      // Copy local rows that are the "same" in both source and target.
      for (LO localRow = 0; localRow < static_cast<LO> (numSameIDs); ++localRow) {
#ifdef HAVE_TPETRA_DEBUG
        if (! srcRowMap.isNodeLocalElement (localRow)) {
          lclErr = true;
          invalidSrcCopyRows.insert (localRow);
          continue; // skip invalid rows
        }
#endif // HAVE_TPETRA_DEBUG

        local_inds_host_view_type lclSrcCols;
        values_host_view_type vals;
        LO numEntries;
        // If this call fails, that means the mesh row local index is
        // invalid.  That means the Import or Export is invalid somehow.
        src->getLocalRowView (localRow, lclSrcCols, vals); numEntries = lclSrcCols.extent(0);
        if (numEntries > 0) {
          LO err = this->replaceLocalValues (localRow, lclSrcCols.data(), reinterpret_cast<const scalar_type*>(vals.data()), numEntries);
          if (err != numEntries) {
            lclErr = true;
            if (! dstRowMap.isNodeLocalElement (localRow)) {
#ifdef HAVE_TPETRA_DEBUG
              invalidDstCopyRows.insert (localRow);
#endif // HAVE_TPETRA_DEBUG
            }
            else {
              // Once there's an error, there's no sense in saving
              // time, so we check whether the column indices were
              // invalid.  However, only remember which ones were
              // invalid in a debug build, because that might take a
              // lot of space.
              for (LO k = 0; k < numEntries; ++k) {
                if (! dstColMap.isNodeLocalElement (lclSrcCols[k])) {
                  lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
                  (void) invalidDstCopyCols.insert (lclSrcCols[k]);
#endif // HAVE_TPETRA_DEBUG
                }
              }
            }
          }
        }
      } // for each "same" local row

      // Copy the "permute" local rows.
      for (size_t k = 0; k < numPermute; ++k) {
        const LO srcLclRow = static_cast<LO> (permuteFromLIDsHost(k));
        const LO dstLclRow = static_cast<LO> (permuteToLIDsHost(k));

        local_inds_host_view_type lclSrcCols;
        values_host_view_type vals;
        LO numEntries;
        src->getLocalRowView (srcLclRow, lclSrcCols, vals); numEntries = lclSrcCols.extent(0);
        if (numEntries > 0) {
          LO err = this->replaceLocalValues (dstLclRow, lclSrcCols.data(), reinterpret_cast<const scalar_type*>(vals.data()), numEntries);
          if (err != numEntries) {
            lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
            for (LO c = 0; c < numEntries; ++c) {
              if (! dstColMap.isNodeLocalElement (lclSrcCols[c])) {
                invalidDstPermuteCols.insert (lclSrcCols[c]);
              }
            }
#endif // HAVE_TPETRA_DEBUG
          }
        }
      }
    }
    else { // must convert column indices to global
      // Reserve space to store the destination matrix's local column indices.
      const size_t maxNumEnt = src->graph_.getLocalMaxNumRowEntries ();
      Teuchos::Array<LO> lclDstCols (maxNumEnt);

      // Copy local rows that are the "same" in both source and target.
      for (LO localRow = 0; localRow < static_cast<LO> (numSameIDs); ++localRow) {
        local_inds_host_view_type lclSrcCols;
        values_host_view_type vals;
        LO numEntries;

        // If this call fails, that means the mesh row local index is
        // invalid.  That means the Import or Export is invalid somehow.
        try {
          src->getLocalRowView (localRow, lclSrcCols, vals); numEntries = lclSrcCols.extent(0);
        } catch (std::exception& e) {
          if (debug) {
            std::ostringstream os;
            const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
            os << "Proc " << myRank << ": copyAndPermute: At \"same\" localRow "
               << localRow << ", src->getLocalRowView() threw an exception: "
               << e.what ();
            std::cerr << os.str ();
          }
          throw e;
        }

        if (numEntries > 0) {
          if (static_cast<size_t> (numEntries) > static_cast<size_t> (lclDstCols.size ())) {
            lclErr = true;
            if (debug) {
              std::ostringstream os;
              const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
              os << "Proc " << myRank << ": copyAndPermute: At \"same\" localRow "
                 << localRow << ", numEntries = " << numEntries << " > maxNumEnt = "
                 << maxNumEnt << endl;
              std::cerr << os.str ();
            }
          }
          else {
            // Convert the source matrix's local column indices to the
            // destination matrix's local column indices.
            Teuchos::ArrayView<LO> lclDstColsView = lclDstCols.view (0, numEntries);
            for (LO j = 0; j < numEntries; ++j) {
              lclDstColsView[j] = dstColMap.getLocalElement (srcColMap.getGlobalElement (lclSrcCols[j]));
              if (lclDstColsView[j] == Teuchos::OrdinalTraits<LO>::invalid ()) {
                lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
                invalidDstCopyCols.insert (lclDstColsView[j]);
#endif // HAVE_TPETRA_DEBUG
              }
            }
            LO err(0);
            try {
              err = this->replaceLocalValues (localRow, lclDstColsView.getRawPtr (), reinterpret_cast<const scalar_type*>(vals.data()), numEntries);
            } catch (std::exception& e) {
              if (debug) {
                std::ostringstream os;
                const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
                os << "Proc " << myRank << ": copyAndPermute: At \"same\" localRow "
                   << localRow << ", this->replaceLocalValues() threw an exception: "
                   << e.what ();
                std::cerr << os.str ();
              }
              throw e;
            }
            if (err != numEntries) {
              lclErr = true;
              if (debug) {
                std::ostringstream os;
                const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
                os << "Proc " << myRank << ": copyAndPermute: At \"same\" "
                  "localRow " << localRow << ", this->replaceLocalValues "
                  "returned " << err << " instead of numEntries = "
                   << numEntries << endl;
                std::cerr << os.str ();
              }
            }
          }
        }
      }

      // Copy the "permute" local rows.
      for (size_t k = 0; k < numPermute; ++k) {
        const LO srcLclRow = static_cast<LO> (permuteFromLIDsHost(k));
        const LO dstLclRow = static_cast<LO> (permuteToLIDsHost(k));

        local_inds_host_view_type lclSrcCols;
        values_host_view_type vals;
        LO numEntries;

        try {
          src->getLocalRowView (srcLclRow, lclSrcCols, vals); numEntries = lclSrcCols.extent(0);
        } catch (std::exception& e) {
          if (debug) {
            std::ostringstream os;
            const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
            os << "Proc " << myRank << ": copyAndPermute: At \"permute\" "
              "srcLclRow " << srcLclRow << " and dstLclRow " << dstLclRow
               << ", src->getLocalRowView() threw an exception: " << e.what ();
            std::cerr << os.str ();
          }
          throw e;
        }

        if (numEntries > 0) {
          if (static_cast<size_t> (numEntries) > static_cast<size_t> (lclDstCols.size ())) {
            lclErr = true;
          }
          else {
            // Convert the source matrix's local column indices to the
            // destination matrix's local column indices.
            Teuchos::ArrayView<LO> lclDstColsView = lclDstCols.view (0, numEntries);
            for (LO j = 0; j < numEntries; ++j) {
              lclDstColsView[j] = dstColMap.getLocalElement (srcColMap.getGlobalElement (lclSrcCols[j]));
              if (lclDstColsView[j] == Teuchos::OrdinalTraits<LO>::invalid ()) {
                lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
                invalidDstPermuteCols.insert (lclDstColsView[j]);
#endif // HAVE_TPETRA_DEBUG
              }
            }
            LO err = this->replaceLocalValues (dstLclRow, lclDstColsView.getRawPtr (), reinterpret_cast<const scalar_type*>(vals.data()), numEntries);
            if (err != numEntries) {
              lclErr = true;
            }
          }
        }
      }
    }

    if (lclErr) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
#ifdef HAVE_TPETRA_DEBUG
      err << "copyAndPermute: The graph structure of the source of the "
        "Import or Export must be a subset of the graph structure of the "
        "target.  ";
      err << "invalidSrcCopyRows = [";
      for (typename std::set<LO>::const_iterator it = invalidSrcCopyRows.begin ();
           it != invalidSrcCopyRows.end (); ++it) {
        err << *it;
        typename std::set<LO>::const_iterator itp1 = it;
        itp1++;
        if (itp1 != invalidSrcCopyRows.end ()) {
          err << ",";
        }
      }
      err << "], invalidDstCopyRows = [";
      for (typename std::set<LO>::const_iterator it = invalidDstCopyRows.begin ();
           it != invalidDstCopyRows.end (); ++it) {
        err << *it;
        typename std::set<LO>::const_iterator itp1 = it;
        itp1++;
        if (itp1 != invalidDstCopyRows.end ()) {
          err << ",";
        }
      }
      err << "], invalidDstCopyCols = [";
      for (typename std::set<LO>::const_iterator it = invalidDstCopyCols.begin ();
           it != invalidDstCopyCols.end (); ++it) {
        err << *it;
        typename std::set<LO>::const_iterator itp1 = it;
        itp1++;
        if (itp1 != invalidDstCopyCols.end ()) {
          err << ",";
        }
      }
      err << "], invalidDstPermuteCols = [";
      for (typename std::set<LO>::const_iterator it = invalidDstPermuteCols.begin ();
           it != invalidDstPermuteCols.end (); ++it) {
        err << *it;
        typename std::set<LO>::const_iterator itp1 = it;
        itp1++;
        if (itp1 != invalidDstPermuteCols.end ()) {
          err << ",";
        }
      }
      err << "], invalidPermuteFromRows = [";
      for (typename std::set<LO>::const_iterator it = invalidPermuteFromRows.begin ();
           it != invalidPermuteFromRows.end (); ++it) {
        err << *it;
        typename std::set<LO>::const_iterator itp1 = it;
        itp1++;
        if (itp1 != invalidPermuteFromRows.end ()) {
          err << ",";
        }
      }
      err << "]" << endl;
#else
      err << "copyAndPermute: The graph structure of the source of the "
        "Import or Export must be a subset of the graph structure of the "
        "target." << endl;
#endif // HAVE_TPETRA_DEBUG
    }

    if (debug) {
      std::ostringstream os;
      const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
      const bool lclSuccess = ! (* (this->localError_));
      os << "*** Proc " << myRank << ": copyAndPermute "
         << (lclSuccess ? "succeeded" : "FAILED");
      if (lclSuccess) {
        os << endl;
      } else {
        os << ": error messages: " << this->errorMessages (); // comes w/ endl
      }
      std::cerr << os.str ();
    }
  }

  namespace { // (anonymous)

    /// \brief Return the (maximum) number of bytes required to pack a
    ///   block row's entries.
    ///
    /// \param numEnt [in] Number of block entries in the row.
    ///
    /// \param numBytesPerValue [in] Maximum number of bytes per
    ///   scalar (not block!) entry (value) of the row.
    ///
    /// \param blkSize [in] Block size of the block sparse matrix.
    ///
    /// If \c Scalar (the type of entries in the matrix) is a plain
    /// old data (POD) type like \c float or \c double, or a struct of
    /// POD (like <tt>std::complex<double></tt>), then the second
    /// argument is just <tt>sizeof(Scalar)</tt>.  If a \c Scalar
    /// instance has a size determined at run time (e.g., when calling
    /// its constructor), then the second argument is the result of
    /// <tt>PackTraits<Scalar>::packValueCount</tt>, called on a
    /// <tt>Scalar</tt> value with the correct run-time size.
    template<class LO, class GO>
    size_t
    packRowCount (const size_t numEnt,
                  const size_t numBytesPerValue,
                  const size_t blkSize)
    {
      using ::Tpetra::Details::PackTraits;

      if (numEnt == 0) {
        // Empty rows always take zero bytes, to ensure sparsity.
        return 0;
      }
      else {
        // We store the number of entries as a local index (LO).
        LO numEntLO = 0; // packValueCount wants this.
        GO gid {};
        const size_t numEntLen = PackTraits<LO>::packValueCount (numEntLO);
        const size_t gidsLen = numEnt * PackTraits<GO>::packValueCount (gid);
        const size_t valsLen = numEnt * numBytesPerValue * blkSize * blkSize;
        return numEntLen + gidsLen + valsLen;
      }
    }

    /// \brief Unpack and return the number of (block) entries in the
    ///   packed row.
    ///
    /// \param imports [in] All the packed data.
    /// \param offset [in] Index of \c imports at which the row starts.
    /// \param numBytes [in] Number of bytes in the packed row.
    /// \param numBytesPerValue [in] Maximum number of bytes per
    ///   scalar (not block!) entry (value) of the row.
    ///
    /// \return Number of (block) entries in the packed row.
    template<class ST, class LO, class GO>
    size_t
    unpackRowCount (const typename ::Tpetra::Details::PackTraits<LO>::input_buffer_type& imports,
                    const size_t offset,
                    const size_t numBytes,
                    const size_t /* numBytesPerValue */)
    {
      using Kokkos::subview;
      using ::Tpetra::Details::PackTraits;

      if (numBytes == 0) {
        // Empty rows always take zero bytes, to ensure sparsity.
        return static_cast<size_t> (0);
      }
      else {
        LO numEntLO = 0;
        const size_t theNumBytes = PackTraits<LO>::packValueCount (numEntLO);
        TEUCHOS_TEST_FOR_EXCEPTION
          (theNumBytes > numBytes, std::logic_error, "unpackRowCount: "
           "theNumBytes = " << theNumBytes << " < numBytes = " << numBytes
           << ".");
        const char* const inBuf = imports.data () + offset;
        const size_t actualNumBytes = PackTraits<LO>::unpackValue (numEntLO, inBuf);
        TEUCHOS_TEST_FOR_EXCEPTION
          (actualNumBytes > numBytes, std::logic_error, "unpackRowCount: "
           "actualNumBytes = " << actualNumBytes << " < numBytes = " << numBytes
           << ".");
        return static_cast<size_t> (numEntLO);
      }
    }

    /// \brief Pack the block row (stored in the input arrays).
    ///
    /// \return The number of bytes packed.
    template<class ST, class LO, class GO>
    size_t
    packRowForBlockCrs (const typename ::Tpetra::Details::PackTraits<LO>::output_buffer_type exports,
                        const size_t offset,
                        const size_t numEnt,
                        const typename ::Tpetra::Details::PackTraits<GO>::input_array_type& gidsIn,
                        const typename ::Tpetra::Details::PackTraits<ST>::input_array_type& valsIn,
                        const size_t numBytesPerValue,
                        const size_t blockSize)
    {
      using Kokkos::subview;
      using ::Tpetra::Details::PackTraits;

      if (numEnt == 0) {
        // Empty rows always take zero bytes, to ensure sparsity.
        return 0;
      }
      const size_t numScalarEnt = numEnt * blockSize * blockSize;

      const GO gid = 0; // packValueCount wants this
      const LO numEntLO = static_cast<size_t> (numEnt);

      const size_t numEntBeg = offset;
      const size_t numEntLen = PackTraits<LO>::packValueCount (numEntLO);
      const size_t gidsBeg = numEntBeg + numEntLen;
      const size_t gidsLen = numEnt * PackTraits<GO>::packValueCount (gid);
      const size_t valsBeg = gidsBeg + gidsLen;
      const size_t valsLen = numScalarEnt * numBytesPerValue;

      char* const numEntOut = exports.data () + numEntBeg;
      char* const gidsOut = exports.data () + gidsBeg;
      char* const valsOut = exports.data () + valsBeg;

      size_t numBytesOut = 0;
      int errorCode = 0;
      numBytesOut += PackTraits<LO>::packValue (numEntOut, numEntLO);

      {
        Kokkos::pair<int, size_t> p;
        p = PackTraits<GO>::packArray (gidsOut, gidsIn.data (), numEnt);
        errorCode += p.first;
        numBytesOut += p.second;

        p = PackTraits<ST>::packArray (valsOut, valsIn.data (), numScalarEnt);
        errorCode += p.first;
        numBytesOut += p.second;
      }

      const size_t expectedNumBytes = numEntLen + gidsLen + valsLen;
      TEUCHOS_TEST_FOR_EXCEPTION
        (numBytesOut != expectedNumBytes, std::logic_error,
         "packRowForBlockCrs: numBytesOut = " << numBytesOut
         << " != expectedNumBytes = " << expectedNumBytes << ".");

      TEUCHOS_TEST_FOR_EXCEPTION
        (errorCode != 0, std::runtime_error, "packRowForBlockCrs: "
         "PackTraits::packArray returned a nonzero error code");

      return numBytesOut;
    }

    // Return the number of bytes actually read / used.
    template<class ST, class LO, class GO>
    size_t
    unpackRowForBlockCrs (const typename ::Tpetra::Details::PackTraits<GO>::output_array_type& gidsOut,
                          const typename ::Tpetra::Details::PackTraits<ST>::output_array_type& valsOut,
                          const typename ::Tpetra::Details::PackTraits<int>::input_buffer_type& imports,
                          const size_t offset,
                          const size_t numBytes,
                          const size_t numEnt,
                          const size_t numBytesPerValue,
                          const size_t blockSize)
    {
      using ::Tpetra::Details::PackTraits;

      if (numBytes == 0) {
        // Rows with zero bytes always have zero entries.
        return 0;
      }
      const size_t numScalarEnt = numEnt * blockSize * blockSize;
      TEUCHOS_TEST_FOR_EXCEPTION
        (static_cast<size_t> (imports.extent (0)) <= offset,
         std::logic_error, "unpackRowForBlockCrs: imports.extent(0) = "
         << imports.extent (0) << " <= offset = " << offset << ".");
      TEUCHOS_TEST_FOR_EXCEPTION
        (static_cast<size_t> (imports.extent (0)) < offset + numBytes,
         std::logic_error, "unpackRowForBlockCrs: imports.extent(0) = "
         << imports.extent (0) << " < offset + numBytes = "
         << (offset + numBytes) << ".");

      const GO gid = 0; // packValueCount wants this
      const LO lid = 0; // packValueCount wants this

      const size_t numEntBeg = offset;
      const size_t numEntLen = PackTraits<LO>::packValueCount (lid);
      const size_t gidsBeg = numEntBeg + numEntLen;
      const size_t gidsLen = numEnt * PackTraits<GO>::packValueCount (gid);
      const size_t valsBeg = gidsBeg + gidsLen;
      const size_t valsLen = numScalarEnt * numBytesPerValue;

      const char* const numEntIn = imports.data () + numEntBeg;
      const char* const gidsIn = imports.data () + gidsBeg;
      const char* const valsIn = imports.data () + valsBeg;

      size_t numBytesOut = 0;
      int errorCode = 0;
      LO numEntOut;
      numBytesOut += PackTraits<LO>::unpackValue (numEntOut, numEntIn);
      TEUCHOS_TEST_FOR_EXCEPTION
        (static_cast<size_t> (numEntOut) != numEnt, std::logic_error,
         "unpackRowForBlockCrs: Expected number of entries " << numEnt
         << " != actual number of entries " << numEntOut << ".");

      {
        Kokkos::pair<int, size_t> p;
        p = PackTraits<GO>::unpackArray (gidsOut.data (), gidsIn, numEnt);
        errorCode += p.first;
        numBytesOut += p.second;

        p = PackTraits<ST>::unpackArray (valsOut.data (), valsIn, numScalarEnt);
        errorCode += p.first;
        numBytesOut += p.second;
      }

      TEUCHOS_TEST_FOR_EXCEPTION
        (numBytesOut != numBytes, std::logic_error,
         "unpackRowForBlockCrs: numBytesOut = " << numBytesOut
         << " != numBytes = " << numBytes << ".");

      const size_t expectedNumBytes = numEntLen + gidsLen + valsLen;
      TEUCHOS_TEST_FOR_EXCEPTION
        (numBytesOut != expectedNumBytes, std::logic_error,
         "unpackRowForBlockCrs: numBytesOut = " << numBytesOut
         << " != expectedNumBytes = " << expectedNumBytes << ".");

      TEUCHOS_TEST_FOR_EXCEPTION
        (errorCode != 0, std::runtime_error, "unpackRowForBlockCrs: "
         "PackTraits::unpackArray returned a nonzero error code");

      return numBytesOut;
    }
  } // namespace (anonymous)

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  packAndPrepare
  (const ::Tpetra::SrcDistObject& source,
   const Kokkos::DualView<const local_ordinal_type*,
     buffer_device_type>& exportLIDs,
   Kokkos::DualView<packet_type*,
     buffer_device_type>& exports, // output
   Kokkos::DualView<size_t*,
     buffer_device_type> numPacketsPerLID, // output
   size_t& constantNumPackets)
  {
    using ::Tpetra::Details::Behavior;
    using ::Tpetra::Details::dualViewStatusToString;
    using ::Tpetra::Details::ProfilingRegion;
    using ::Tpetra::Details::PackTraits;

    typedef typename Kokkos::View<int*, device_type>::HostMirror::execution_space host_exec;

    typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_BCRS_type;

    ProfilingRegion profile_region("Tpetra::BlockCrsMatrix::packAndPrepare");

    const bool debug = Behavior::debug();
    const bool verbose = Behavior::verbose();

    // Define this function prefix
    std::string prefix;
    {
      std::ostringstream os;
      const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
      os << "Proc " << myRank << ": BlockCrsMatrix::packAndPrepare: " << std::endl;
      prefix = os.str();
    }

    // check if this already includes a local error
    if (* (this->localError_)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix
          << "The target object of the Import or Export is already in an error state."
          << std::endl;
      return;
    }

    //
    // Verbose input dual view status
    //
    if (verbose) {
      std::ostringstream os;
      os << prefix << std::endl
         << prefix << "  " << dualViewStatusToString (exportLIDs, "exportLIDs") << std::endl
         << prefix << "  " << dualViewStatusToString (exports, "exports") << std::endl
         << prefix << "  " << dualViewStatusToString (numPacketsPerLID, "numPacketsPerLID") << std::endl;
      std::cerr << os.str ();
    }

    ///
    /// Check input valid
    ///
    if (exportLIDs.extent (0) != numPacketsPerLID.extent (0)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix
          << "exportLIDs.extent(0) = " << exportLIDs.extent (0)
          << " != numPacketsPerLID.extent(0) = " << numPacketsPerLID.extent(0)
          << "." << std::endl;
      return;
    }
    if (exportLIDs.need_sync_host ()) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix << "exportLIDs be sync'd to host." << std::endl;
      return;
    }

    const this_BCRS_type* src = dynamic_cast<const this_BCRS_type* > (&source);
    if (src == NULL) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix
          << "The source (input) object of the Import or "
        "Export is either not a BlockCrsMatrix, or does not have the right "
        "template parameters.  checkSizes() should have caught this.  "
        "Please report this bug to the Tpetra developers." << std::endl;
      return;
    }

    // Graphs and matrices are allowed to have a variable number of
    // entries per row.  We could test whether all rows have the same
    // number of entries, but DistObject can only use this
    // optimization if all rows on _all_ processes have the same
    // number of entries.  Rather than do the all-reduce necessary to
    // test for this unlikely case, we tell DistObject (by setting
    // constantNumPackets to zero) to assume that different rows may
    // have different numbers of entries.
    constantNumPackets = 0;

    // const values
    const crs_graph_type& srcGraph = src->graph_;
    const size_t blockSize = static_cast<size_t> (src->getBlockSize ());
    const size_t numExportLIDs = exportLIDs.extent (0);
    size_t numBytesPerValue(0);
    {
      auto val_host = val_.getHostView(Access::ReadOnly);
      numBytesPerValue =
        PackTraits<impl_scalar_type>
        ::packValueCount(val_host.extent(0) ? val_host(0) : impl_scalar_type());
    }

    // Compute the number of bytes ("packets") per row to pack.  While
    // we're at it, compute the total # of block entries to send, and
    // the max # of block entries in any of the rows we're sending.

    Impl::BlockCrsRowStruct<size_t> rowReducerStruct;

    // Graph information is on host; let's do this on host parallel reduce
    auto exportLIDsHost = exportLIDs.view_host();
    auto numPacketsPerLIDHost = numPacketsPerLID.view_host(); // we will modify this.
    numPacketsPerLID.modify_host ();
    {
      rowReducerStruct = Impl::BlockCrsRowStruct<size_t>();
      for (size_t i = 0; i < numExportLIDs; ++i) {
          const LO lclRow = exportLIDsHost(i);
          size_t numEnt = srcGraph.getNumEntriesInLocalRow (lclRow);
          numEnt = (numEnt == Teuchos::OrdinalTraits<size_t>::invalid () ? 0 : numEnt);

          const size_t numBytes =
            packRowCount<LO, GO> (numEnt, numBytesPerValue, blockSize);
          numPacketsPerLIDHost(i) = numBytes;
          rowReducerStruct += Impl::BlockCrsRowStruct<size_t>(numEnt, numBytes, numEnt);
        }
    }

    // Compute the number of bytes ("packets") per row to pack.  While
    // we're at it, compute the total # of block entries to send, and
    // the max # of block entries in any of the rows we're sending.
    const size_t totalNumBytes   = rowReducerStruct.totalNumBytes;
    const size_t totalNumEntries = rowReducerStruct.totalNumEntries;
    const size_t maxRowLength    = rowReducerStruct.maxRowLength;

    if (verbose) {
      std::ostringstream os;
      os << prefix
         << "totalNumBytes = " << totalNumBytes << ", totalNumEntries = " << totalNumEntries
         << std::endl;
      std::cerr << os.str ();
    }

    // We use a "struct of arrays" approach to packing each row's
    // entries.  All the column indices (as global indices) go first,
    // then all their owning process ranks, and then the values.
    if(exports.extent(0) != totalNumBytes)
    {
      const std::string oldLabel = exports.d_view.label ();
      const std::string newLabel = (oldLabel == "") ? "exports" : oldLabel;
      exports = Kokkos::DualView<packet_type*, buffer_device_type>(newLabel, totalNumBytes);
    }
    if (totalNumEntries > 0) {
      // Current position (in bytes) in the 'exports' output array.
      Kokkos::View<size_t*, host_exec> offset("offset", numExportLIDs+1);
      {
        const auto policy = Kokkos::RangePolicy<host_exec>(size_t(0), numExportLIDs+1);
        Kokkos::parallel_scan
          (policy,
           [=](const size_t &i, size_t &update, const bool &final) {
            if (final) offset(i) = update;
            update += (i == numExportLIDs ? 0 : numPacketsPerLIDHost(i));
          });
      }
      if (offset(numExportLIDs) != totalNumBytes) {
        std::ostream& err = this->markLocalErrorAndGetStream ();
        err << prefix
            << "At end of method, the final offset (in bytes) "
            << offset(numExportLIDs) << " does not equal the total number of bytes packed "
            << totalNumBytes << ".  "
            << "Please report this bug to the Tpetra developers." << std::endl;
        return;
      }

      // For each block row of the matrix owned by the calling
      // process, pack that block row's column indices and values into
      // the exports array.

      // Source matrix's column Map.  We verified in checkSizes() that
      // the column Map exists (is not null).
      const map_type& srcColMap = * (srcGraph.getColMap ());

      // Pack the data for each row to send, into the 'exports' buffer.
      // exports will be modified on host.
      exports.modify_host();
      {
        typedef Kokkos::TeamPolicy<host_exec> policy_type;
        const auto policy =
          policy_type(numExportLIDs, 1, 1)
          .set_scratch_size(0, Kokkos::PerTeam(sizeof(GO)*maxRowLength));
        // The following parallel_for needs const access to the local values of src.
        // (the local graph is also accessed on host, but this does not use WDVs).
        getValuesHost();
        Details::disableWDVTracking();
        Kokkos::parallel_for
          (policy,
           [=](const typename policy_type::member_type &member) {
            const size_t i = member.league_rank();
            Kokkos::View<GO*, typename host_exec::scratch_memory_space>
              gblColInds(member.team_scratch(0), maxRowLength);

            const LO  lclRowInd = exportLIDsHost(i);
            local_inds_host_view_type lclColInds;
            values_host_view_type vals;

            // It's OK to ignore the return value, since if the calling
            // process doesn't own that local row, then the number of
            // entries in that row on the calling process is zero.
            src->getLocalRowView (lclRowInd, lclColInds, vals); 
            const size_t numEnt = lclColInds.extent(0);

            // Convert column indices from local to global.
            for (size_t j = 0; j < numEnt; ++j)
              gblColInds(j) = srcColMap.getGlobalElement (lclColInds(j));

            // Kyungjoo: additional wrapping scratch view is necessary
            //   the following function interface need the same execution space
            //   host scratch space somehow is not considered same as the host_exec
            // Copy the row's data into the current spot in the exports array.
            const size_t numBytes =
              packRowForBlockCrs<impl_scalar_type, LO, GO>
              (exports.view_host(),
               offset(i),
               numEnt,
               Kokkos::View<const GO*, host_exec>(gblColInds.data(), maxRowLength),
               vals,
               numBytesPerValue,
               blockSize);

            // numBytes should be same as the difference between offsets
            if (debug) {
              const size_t offsetDiff = offset(i+1) - offset(i);
              if (numBytes != offsetDiff) {
                std::ostringstream os;
                os << prefix
                   << "numBytes computed from packRowForBlockCrs is different from "
                   << "precomputed offset values, LID = " << i << std::endl;
                std::cerr << os.str ();
              }
            }
          }); // for each LID (of a row) to send
        Details::enableWDVTracking();
      }
    } // if totalNumEntries > 0

    if (debug) {
      std::ostringstream os;
      const bool lclSuccess = ! (* (this->localError_));
      os << prefix
         << (lclSuccess ? "succeeded" : "FAILED")
         << std::endl;
      std::cerr << os.str ();
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  unpackAndCombine
  (const Kokkos::DualView<const local_ordinal_type*,
     buffer_device_type>& importLIDs,
   Kokkos::DualView<packet_type*,
     buffer_device_type> imports,
   Kokkos::DualView<size_t*,
     buffer_device_type> numPacketsPerLID,
   const size_t /* constantNumPackets */,
   const CombineMode combineMode)
  {
    using ::Tpetra::Details::Behavior;
    using ::Tpetra::Details::dualViewStatusToString;
    using ::Tpetra::Details::ProfilingRegion;
    using ::Tpetra::Details::PackTraits;
    using std::endl;
    using host_exec =
      typename Kokkos::View<int*, device_type>::HostMirror::execution_space;

    ProfilingRegion profile_region("Tpetra::BlockCrsMatrix::unpackAndCombine");
    const bool verbose = Behavior::verbose ();

    // Define this function prefix
    std::string prefix;
    {
      std::ostringstream os;
      auto map = this->graph_.getRowMap ();
      auto comm = map.is_null () ? Teuchos::null : map->getComm ();
      const int myRank = comm.is_null () ? -1 : comm->getRank ();
      os << "Proc " << myRank << ": Tpetra::BlockCrsMatrix::unpackAndCombine: ";
      prefix = os.str ();
      if (verbose) {
        os << "Start" << endl;
        std::cerr << os.str ();
      }
    }

    // check if this already includes a local error
    if (* (this->localError_)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      std::ostringstream os;
      os << prefix << "{Im/Ex}port target is already in error." << endl;
      if (verbose) {
        std::cerr << os.str ();
      }
      err << os.str ();
      return;
    }

    ///
    /// Check input valid
    ///
    if (importLIDs.extent (0) != numPacketsPerLID.extent (0)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      std::ostringstream os;
      os << prefix << "importLIDs.extent(0) = " << importLIDs.extent (0)
         << " != numPacketsPerLID.extent(0) = " << numPacketsPerLID.extent(0)
         << "." << endl;
      if (verbose) {
        std::cerr << os.str ();
      }
      err << os.str ();
      return;
    }

    if (combineMode != ADD     && combineMode != INSERT &&
        combineMode != REPLACE && combineMode != ABSMAX &&
        combineMode != ZERO) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      std::ostringstream os;
      os << prefix
         << "Invalid CombineMode value " << combineMode << ".  Valid "
         << "values include ADD, INSERT, REPLACE, ABSMAX, and ZERO."
         << std::endl;
      if (verbose) {
        std::cerr << os.str ();
      }
      err << os.str ();
      return;
    }

    if (this->graph_.getColMap ().is_null ()) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      std::ostringstream os;
      os << prefix << "Target matrix's column Map is null." << endl;
      if (verbose) {
        std::cerr << os.str ();
      }
      err << os.str ();
      return;
    }

    // Target matrix's column Map.  Use to convert the global column
    // indices in the receive buffer to local indices.  We verified in
    // checkSizes() that the column Map exists (is not null).
    const map_type& tgtColMap = * (this->graph_.getColMap ());

    // Const values
    const size_t blockSize = this->getBlockSize ();
    const size_t numImportLIDs = importLIDs.extent(0);
    // FIXME (mfh 06 Feb 2019) For scalar types with run-time sizes, a
    // default-constructed instance could have a different size than
    // other instances.  (We assume all nominally constructed
    // instances have the same size; that's not the issue here.)  This
    // could be bad if the calling process has no entries, but other
    // processes have entries that they want to send to us.
    size_t numBytesPerValue(0);
    {
      auto val_host = val_.getHostView(Access::ReadOnly);
      numBytesPerValue =
        PackTraits<impl_scalar_type>::packValueCount
        (val_host.extent (0) ? val_host(0) : impl_scalar_type ());
    }
    const size_t maxRowNumEnt = graph_.getLocalMaxNumRowEntries ();
    const size_t maxRowNumScalarEnt = maxRowNumEnt * blockSize * blockSize;

    if (verbose) {
      std::ostringstream os;
      os << prefix << "combineMode: "
         << ::Tpetra::combineModeToString (combineMode)
         << ", blockSize: " << blockSize
         << ", numImportLIDs: " << numImportLIDs
         << ", numBytesPerValue: " << numBytesPerValue
         << ", maxRowNumEnt: " << maxRowNumEnt
         << ", maxRowNumScalarEnt: " << maxRowNumScalarEnt << endl;
      std::cerr << os.str ();
    }

    if (combineMode == ZERO || numImportLIDs == 0) {
      if (verbose) {
        std::ostringstream os;
        os << prefix << "Nothing to unpack. Done!" << std::endl;
        std::cerr << os.str ();
      }
      return;
    }

    // NOTE (mfh 07 Feb 2019) If we ever implement unpack on device,
    // we can remove this sync.
    if (imports.need_sync_host ()) {
      imports.sync_host ();
    }

    // NOTE (mfh 07 Feb 2019) DistObject::doTransferNew has already
    // sync'd numPacketsPerLID to host, since it needs to do that in
    // order to post MPI_Irecv messages with the correct lengths.  A
    // hypothetical device-specific boundary exchange implementation
    // could possibly receive data without sync'ing lengths to host,
    // but we don't need to design for that nonexistent thing yet.

    if (imports.need_sync_host () || numPacketsPerLID.need_sync_host () ||
        importLIDs.need_sync_host ()) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      std::ostringstream os;
      os << prefix << "All input DualViews must be sync'd to host by now. "
         << "imports_nc.need_sync_host()="
         << (imports.need_sync_host () ? "true" : "false")
         << ", numPacketsPerLID.need_sync_host()="
         << (numPacketsPerLID.need_sync_host () ? "true" : "false")
         << ", importLIDs.need_sync_host()="
         << (importLIDs.need_sync_host () ? "true" : "false")
         << "." << endl;
      if (verbose) {
        std::cerr << os.str ();
      }
      err << os.str ();
      return;
    }

    const auto importLIDsHost = importLIDs.view_host ();
    const auto numPacketsPerLIDHost = numPacketsPerLID.view_host ();

    // FIXME (mfh 06 Feb 2019) We could fuse the scan with the unpack
    // loop, by only unpacking on final==true (when we know the
    // current offset's value).

    Kokkos::View<size_t*, host_exec> offset ("offset", numImportLIDs+1);
    {
      const auto policy = Kokkos::RangePolicy<host_exec>(size_t(0), numImportLIDs+1);
      Kokkos::parallel_scan
        ("Tpetra::BlockCrsMatrix::unpackAndCombine: offsets", policy,
         [=] (const size_t &i, size_t &update, const bool &final) {
          if (final) offset(i) = update;
          update += (i == numImportLIDs ? 0 : numPacketsPerLIDHost(i));
        });
    }

    // this variable does not matter with a race condition (just error flag)
    //
    // NOTE (mfh 06 Feb 2019) CUDA doesn't necessarily like reductions
    // or atomics on bool, so we use int instead.  (I know we're not
    // launching a CUDA loop here, but we could in the future, and I'd
    // like to avoid potential trouble.)
    Kokkos::View<int, host_exec, Kokkos::MemoryTraits<Kokkos::Atomic> >
      errorDuringUnpack ("errorDuringUnpack");
    errorDuringUnpack () = 0;
    {
      using policy_type = Kokkos::TeamPolicy<host_exec>;
      size_t scratch_per_row = sizeof(GO) * maxRowNumEnt + sizeof (LO) * maxRowNumEnt + numBytesPerValue * maxRowNumScalarEnt
        + 2 * sizeof(GO); // Yeah, this is a fudge factor

      const auto policy = policy_type (numImportLIDs, 1, 1)     
        .set_scratch_size (0, Kokkos::PerTeam (scratch_per_row));
      using host_scratch_space = typename host_exec::scratch_memory_space;
      
      using pair_type = Kokkos::pair<size_t, size_t>;

      //The following parallel_for modifies values on host while unpacking.
      getValuesHostNonConst();
      Details::disableWDVTracking();
      Kokkos::parallel_for
        ("Tpetra::BlockCrsMatrix::unpackAndCombine: unpack", policy,
         [=] (const typename policy_type::member_type& member) {
          const size_t i = member.league_rank();
          Kokkos::View<GO*, host_scratch_space> gblColInds
            (member.team_scratch (0), maxRowNumEnt);
          Kokkos::View<LO*, host_scratch_space> lclColInds
            (member.team_scratch (0), maxRowNumEnt);
          Kokkos::View<impl_scalar_type*, host_scratch_space> vals
            (member.team_scratch (0), maxRowNumScalarEnt);
          

          const size_t offval = offset(i);
          const LO lclRow = importLIDsHost(i);
          const size_t numBytes = numPacketsPerLIDHost(i);
          const size_t numEnt =
            unpackRowCount<impl_scalar_type, LO, GO>
            (imports.view_host (), offval, numBytes, numBytesPerValue);

          if (numBytes > 0) {
            if (numEnt > maxRowNumEnt) {
              errorDuringUnpack() = 1;
              if (verbose) {
                std::ostringstream os;
                os << prefix
                   << "At i = " << i << ", numEnt = " << numEnt
                   << " > maxRowNumEnt = " << maxRowNumEnt
                   << std::endl;
                std::cerr << os.str ();
              }
              return;
            }
          }
          const size_t numScalarEnt = numEnt*blockSize*blockSize;
          auto gidsOut = Kokkos::subview(gblColInds, pair_type(0, numEnt));
          auto lidsOut = Kokkos::subview(lclColInds, pair_type(0, numEnt));
          auto valsOut = Kokkos::subview(vals,       pair_type(0, numScalarEnt));

          // Kyungjoo: additional wrapping scratch view is necessary
          //   the following function interface need the same execution space
          //   host scratch space somehow is not considered same as the host_exec
          size_t numBytesOut = 0;
          try {
            numBytesOut =
              unpackRowForBlockCrs<impl_scalar_type, LO, GO>
              (Kokkos::View<GO*,host_exec>(gidsOut.data(), numEnt),
               Kokkos::View<impl_scalar_type*,host_exec>(valsOut.data(), numScalarEnt),
               imports.view_host(),
               offval, numBytes, numEnt,
               numBytesPerValue, blockSize);
          }
          catch (std::exception& e) {
            errorDuringUnpack() = 1;
            if (verbose) {
              std::ostringstream os;
              os << prefix << "At i = " << i << ", unpackRowForBlockCrs threw: "
                 << e.what () << endl;
              std::cerr << os.str ();
            }
            return;
          }

          if (numBytes != numBytesOut) {
            errorDuringUnpack() = 1;
            if (verbose) {
              std::ostringstream os;
              os << prefix
                 << "At i = " << i << ", numBytes = " << numBytes
                 << " != numBytesOut = " << numBytesOut << "."
                 << std::endl;
              std::cerr << os.str();
            }
            return;
          }

          // Convert incoming global indices to local indices.
          for (size_t k = 0; k < numEnt; ++k) {
            lidsOut(k) = tgtColMap.getLocalElement (gidsOut(k));
            if (lidsOut(k) == Teuchos::OrdinalTraits<LO>::invalid ()) {
              if (verbose) {
                std::ostringstream os;
                os << prefix
                    << "At i = " << i << ", GID " << gidsOut(k)
                    << " is not owned by the calling process."
                    << std::endl;
                std::cerr << os.str();
              }
              return;
            }
          }

          // Combine the incoming data with the matrix's current data.
          LO numCombd = 0;
          const LO* const lidsRaw = const_cast<const LO*> (lidsOut.data ());
          const Scalar* const valsRaw = reinterpret_cast<const Scalar*>
            (const_cast<const impl_scalar_type*> (valsOut.data ()));
          if (combineMode == ADD) {
            numCombd = this->sumIntoLocalValues (lclRow, lidsRaw, valsRaw, numEnt);
          } else if (combineMode == INSERT || combineMode == REPLACE) {
            numCombd = this->replaceLocalValues (lclRow, lidsRaw, valsRaw, numEnt);
          } else if (combineMode == ABSMAX) {
            numCombd = this->absMaxLocalValues (lclRow, lidsRaw, valsRaw, numEnt);
          }

          if (static_cast<LO> (numEnt) != numCombd) {
            errorDuringUnpack() = 1;
            if (verbose) {
              std::ostringstream os;
              os << prefix << "At i = " << i << ", numEnt = " << numEnt
                 << " != numCombd = " << numCombd << "."
                 << endl;
              std::cerr << os.str ();
            }
            return;
          }
        }); // for each import LID i
      Details::enableWDVTracking();
    }

    if (errorDuringUnpack () != 0) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix << "Unpacking failed.";
      if (! verbose) {
        err << "  Please run again with the environment variable setting "
          "TPETRA_VERBOSE=1 to get more verbose diagnostic output.";
      }
      err << endl;
    }

    if (verbose) {
      std::ostringstream os;
      const bool lclSuccess = ! (* (this->localError_));
      os << prefix
         << (lclSuccess ? "succeeded" : "FAILED")
         << std::endl;
      std::cerr << os.str ();
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  std::string
  BlockCrsMatrix<Scalar, LO, GO, Node>::description () const
  {
    using Teuchos::TypeNameTraits;
    std::ostringstream os;
    os << "\"Tpetra::BlockCrsMatrix\": { "
       << "Template parameters: { "
       << "Scalar: " << TypeNameTraits<Scalar>::name ()
       << "LO: " << TypeNameTraits<LO>::name ()
       << "GO: " << TypeNameTraits<GO>::name ()
       << "Node: " << TypeNameTraits<Node>::name ()
       << " }"
       << ", Label: \"" << this->getObjectLabel () << "\""
       << ", Global dimensions: ["
       << graph_.getDomainMap ()->getGlobalNumElements () << ", "
       << graph_.getRangeMap ()->getGlobalNumElements () << "]"
       << ", Block size: " << getBlockSize ()
       << " }";
    return os.str ();
  }


  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::CommRequest;
    using Teuchos::FancyOStream;
    using Teuchos::getFancyOStream;
    using Teuchos::ireceive;
    using Teuchos::isend;
    using Teuchos::outArg;
    using Teuchos::TypeNameTraits;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    // using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    using Teuchos::RCP;
    using Teuchos::wait;
    using std::endl;

    // Set default verbosity if applicable.
    const Teuchos::EVerbosityLevel vl =
      (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

    if (vl == VERB_NONE) {
      return; // print nothing
    }

    // describe() always starts with a tab before it prints anything.
    Teuchos::OSTab tab0 (out);

    out << "\"Tpetra::BlockCrsMatrix\":" << endl;
    Teuchos::OSTab tab1 (out);

    out << "Template parameters:" << endl;
    {
      Teuchos::OSTab tab2 (out);
      out << "Scalar: " << TypeNameTraits<Scalar>::name () << endl
          << "LO: " << TypeNameTraits<LO>::name () << endl
          << "GO: " << TypeNameTraits<GO>::name () << endl
          << "Node: " << TypeNameTraits<Node>::name () << endl;
    }
    out << "Label: \"" << this->getObjectLabel () << "\"" << endl
        << "Global dimensions: ["
        << graph_.getDomainMap ()->getGlobalNumElements () << ", "
        << graph_.getRangeMap ()->getGlobalNumElements () << "]" << endl;

    const LO blockSize = getBlockSize ();
    out << "Block size: " << blockSize << endl;

    // constituent objects
    if (vl >= VERB_MEDIUM) {
      const Teuchos::Comm<int>& comm = * (graph_.getMap ()->getComm ());
      const int myRank = comm.getRank ();
      if (myRank == 0) {
        out << "Row Map:" << endl;
      }
      getRowMap()->describe (out, vl);

      if (! getColMap ().is_null ()) {
        if (getColMap() == getRowMap()) {
          if (myRank == 0) {
            out << "Column Map: same as row Map" << endl;
          }
        }
        else {
          if (myRank == 0) {
            out << "Column Map:" << endl;
          }
          getColMap ()->describe (out, vl);
        }
      }
      if (! getDomainMap ().is_null ()) {
        if (getDomainMap () == getRowMap ()) {
          if (myRank == 0) {
            out << "Domain Map: same as row Map" << endl;
          }
        }
        else if (getDomainMap () == getColMap ()) {
          if (myRank == 0) {
            out << "Domain Map: same as column Map" << endl;
          }
        }
        else {
          if (myRank == 0) {
            out << "Domain Map:" << endl;
          }
          getDomainMap ()->describe (out, vl);
        }
      }
      if (! getRangeMap ().is_null ()) {
        if (getRangeMap () == getDomainMap ()) {
          if (myRank == 0) {
            out << "Range Map: same as domain Map" << endl;
          }
        }
        else if (getRangeMap () == getRowMap ()) {
          if (myRank == 0) {
            out << "Range Map: same as row Map" << endl;
          }
        }
        else {
          if (myRank == 0) {
            out << "Range Map: " << endl;
          }
          getRangeMap ()->describe (out, vl);
        }
      }
    }

    if (vl >= VERB_EXTREME) {
      const Teuchos::Comm<int>& comm = * (graph_.getMap ()->getComm ());
      const int myRank = comm.getRank ();
      const int numProcs = comm.getSize ();

      // Print the calling process' data to the given output stream.
      RCP<std::ostringstream> lclOutStrPtr (new std::ostringstream ());
      RCP<FancyOStream> osPtr = getFancyOStream (lclOutStrPtr);
      FancyOStream& os = *osPtr;
      os << "Process " << myRank << ":" << endl;
      Teuchos::OSTab tab2 (os);

      const map_type& meshRowMap = * (graph_.getRowMap ());
      const map_type& meshColMap = * (graph_.getColMap ());
      for (LO meshLclRow = meshRowMap.getMinLocalIndex ();
           meshLclRow <= meshRowMap.getMaxLocalIndex ();
           ++meshLclRow) {
        const GO meshGblRow = meshRowMap.getGlobalElement (meshLclRow);
        os << "Row " << meshGblRow << ": {";

        local_inds_host_view_type lclColInds;
        values_host_view_type vals;
        LO numInds = 0;
        this->getLocalRowView (meshLclRow, lclColInds, vals); numInds = lclColInds.extent(0);

        for (LO k = 0; k < numInds; ++k) {
          const GO gblCol = meshColMap.getGlobalElement (lclColInds[k]);

          os << "Col " << gblCol << ": [";
          for (LO i = 0; i < blockSize; ++i) {
            for (LO j = 0; j < blockSize; ++j) {
              os << vals[blockSize*blockSize*k + i*blockSize + j];
              if (j + 1 < blockSize) {
                os << ", ";
              }
            }
            if (i + 1 < blockSize) {
              os << "; ";
            }
          }
          os << "]";
          if (k + 1 < numInds) {
            os << ", ";
          }
        }
        os << "}" << endl;
      }

      // Print data on Process 0.  This will automatically respect the
      // current indentation level.
      if (myRank == 0) {
        out << lclOutStrPtr->str ();
        lclOutStrPtr = Teuchos::null; // clear it to save space
      }

      const int sizeTag = 1337;
      const int dataTag = 1338;

      ArrayRCP<char> recvDataBuf; // only used on Process 0

      // Send string sizes and data from each process in turn to
      // Process 0, and print on that process.
      for (int p = 1; p < numProcs; ++p) {
        if (myRank == 0) {
          // Receive the incoming string's length.
          ArrayRCP<size_t> recvSize (1);
          recvSize[0] = 0;
          RCP<CommRequest<int> > recvSizeReq =
            ireceive<int, size_t> (recvSize, p, sizeTag, comm);
          wait<int> (comm, outArg (recvSizeReq));
          const size_t numCharsToRecv = recvSize[0];

          // Allocate space for the string to receive.  Reuse receive
          // buffer space if possible.  We can do this because in the
          // current implementation, we only have one receive in
          // flight at a time.  Leave space for the '\0' at the end,
          // in case the sender doesn't send it.
          if (static_cast<size_t>(recvDataBuf.size()) < numCharsToRecv + 1) {
            recvDataBuf.resize (numCharsToRecv + 1);
          }
          ArrayRCP<char> recvData = recvDataBuf.persistingView (0, numCharsToRecv);
          // Post the receive of the actual string data.
          RCP<CommRequest<int> > recvDataReq =
            ireceive<int, char> (recvData, p, dataTag, comm);
          wait<int> (comm, outArg (recvDataReq));

          // Print the received data.  This will respect the current
          // indentation level.  Make sure that the string is
          // null-terminated.
          recvDataBuf[numCharsToRecv] = '\0';
          out << recvDataBuf.getRawPtr ();
        }
        else if (myRank == p) { // if I am not Process 0, and my rank is p
          // This deep-copies the string at most twice, depending on
          // whether std::string reference counts internally (it
          // generally does, so this won't deep-copy at all).
          const std::string stringToSend = lclOutStrPtr->str ();
          lclOutStrPtr = Teuchos::null; // clear original to save space

          // Send the string's length to Process 0.
          const size_t numCharsToSend = stringToSend.size ();
          ArrayRCP<size_t> sendSize (1);
          sendSize[0] = numCharsToSend;
          RCP<CommRequest<int> > sendSizeReq =
            isend<int, size_t> (sendSize, 0, sizeTag, comm);
          wait<int> (comm, outArg (sendSizeReq));

          // Send the actual string to Process 0.  We know that the
          // string has length > 0, so it's save to take the address
          // of the first entry.  Make a nonowning ArrayRCP to hold
          // the string.  Process 0 will add a null termination
          // character at the end of the string, after it receives the
          // message.
          ArrayRCP<const char> sendData (&stringToSend[0], 0, numCharsToSend, false);
          RCP<CommRequest<int> > sendDataReq =
            isend<int, char> (sendData, 0, dataTag, comm);
          wait<int> (comm, outArg (sendDataReq));
        }
      } // for each process rank p other than 0
    } // extreme verbosity level (print the whole matrix)
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<const Teuchos::Comm<int> >
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getComm() const
  {
    return graph_.getComm();
  }


  template<class Scalar, class LO, class GO, class Node>
  global_size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalNumCols() const
  {
    return graph_.getGlobalNumCols();
  }


  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalNumCols() const
  {
    return graph_.getLocalNumCols();
  }

  template<class Scalar, class LO, class GO, class Node>
  GO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getIndexBase() const
  {
    return graph_.getIndexBase();
  }

  template<class Scalar, class LO, class GO, class Node>
  global_size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalNumEntries() const
  {
    return graph_.getGlobalNumEntries();
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalNumEntries() const
  {
    return graph_.getLocalNumEntries();
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNumEntriesInGlobalRow (GO globalRow) const
  {
    return graph_.getNumEntriesInGlobalRow(globalRow);
  }


  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalMaxNumRowEntries() const
  {
    return graph_.getGlobalMaxNumRowEntries();
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  hasColMap() const
  {
    return graph_.hasColMap();
  }


  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  isLocallyIndexed() const
  {
    return graph_.isLocallyIndexed();
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  isGloballyIndexed() const
  {
    return graph_.isGloballyIndexed();
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  isFillComplete() const
  {
    return graph_.isFillComplete ();
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  supportsRowViews() const
  {
    return false;
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalRowCopy (GO /*GlobalRow*/,
                    nonconst_global_inds_host_view_type &/*Indices*/,
                    nonconst_values_host_view_type &/*Values*/,
                    size_t& /*NumEntries*/) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::BlockCrsMatrix::getGlobalRowCopy: "
      "This class doesn't support global matrix indexing.");

  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalRowView (GO /* GlobalRow */,
                    global_inds_host_view_type &/* indices */,
                    values_host_view_type &/* values */) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::BlockCrsMatrix::getGlobalRowView: "
      "This class doesn't support global matrix indexing.");

  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalRowView (LO localRowInd,
                   local_inds_host_view_type &colInds,
                   values_host_view_type &vals) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      colInds = local_inds_host_view_type();
      vals = values_host_view_type();
    }
    else {
      const size_t absBlockOffsetStart = ptrHost_[localRowInd];
      const size_t numInds = ptrHost_[localRowInd + 1] - absBlockOffsetStart;
      colInds = local_inds_host_view_type(indHost_.data()+absBlockOffsetStart, numInds);

      vals = getValuesHost (localRowInd);
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalRowViewNonConst (LO localRowInd,
                           local_inds_host_view_type &colInds,
                           nonconst_values_host_view_type &vals) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      colInds = nonconst_local_inds_host_view_type();
      vals = nonconst_values_host_view_type();
    }
    else {
      const size_t absBlockOffsetStart = ptrHost_[localRowInd];
      const size_t numInds = ptrHost_[localRowInd + 1] - absBlockOffsetStart;
      colInds = local_inds_host_view_type(indHost_.data()+absBlockOffsetStart, numInds);

      using this_BCRS_type = BlockCrsMatrix<Scalar, LO, GO, Node>;
      vals = const_cast<this_BCRS_type&>(*this).getValuesHostNonConst(localRowInd);
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalDiagCopy (::Tpetra::Vector<Scalar,LO,GO,Node>& diag) const
  {
    const size_t lclNumMeshRows = graph_.getLocalNumRows ();

    Kokkos::View<size_t*, device_type> diagOffsets ("diagOffsets", lclNumMeshRows);
    graph_.getLocalDiagOffsets (diagOffsets);

    // The code below works on host, so use a host View.
    auto diagOffsetsHost = Kokkos::create_mirror_view (diagOffsets);
    // DEEP_COPY REVIEW - DEVICE-TO-HOSTMIRROR
    Kokkos::deep_copy (execution_space(), diagOffsetsHost, diagOffsets);

    auto vals_host_out = val_.getHostView(Access::ReadOnly);
    const impl_scalar_type* vals_host_out_raw = vals_host_out.data();

    // TODO amk: This is a temporary measure to make the code run with Ifpack2
    size_t rowOffset = 0;
    size_t offset = 0;
    LO bs = getBlockSize();
    for(size_t r=0; r<getLocalNumRows(); r++)
    {
      // move pointer to start of diagonal block
      offset = rowOffset + diagOffsetsHost(r)*bs*bs;
      for(int b=0; b<bs; b++)
      {
        diag.replaceLocalValue(r*bs+b, vals_host_out_raw[offset+b*(bs+1)]);
      }
      // move pointer to start of next block row
      rowOffset += getNumEntriesInLocalRow(r)*bs*bs;
    }

    //diag.template sync<memory_space> (); // sync vec of diag entries back to dev
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  leftScale (const ::Tpetra::Vector<Scalar, LO, GO, Node>& /* x */)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::BlockCrsMatrix::leftScale: "
      "not implemented.");

  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  rightScale (const ::Tpetra::Vector<Scalar, LO, GO, Node>& /* x */)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::BlockCrsMatrix::rightScale: "
      "not implemented.");

  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<const ::Tpetra::RowGraph<LO, GO, Node> >
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGraph() const
  {
    return graphRCP_;
  }

  template<class Scalar, class LO, class GO, class Node>
  typename ::Tpetra::RowMatrix<Scalar, LO, GO, Node>::mag_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getFrobeniusNorm () const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::BlockCrsMatrix::getFrobeniusNorm: "
      "not implemented.");
  }

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_BLOCKCRSMATRIX_INSTANT(S,LO,GO,NODE) \
  template class BlockCrsMatrix< S, LO, GO, NODE >;

#endif // TPETRA_BLOCKCRSMATRIX_DEF_HPP

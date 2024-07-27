// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_CRSMATRIX_DEF_HPP
#define TPETRA_CRSMATRIX_DEF_HPP

/// \file Tpetra_CrsMatrix_def.hpp
/// \brief Definition of the Tpetra::CrsMatrix class
///
/// If you want to use Tpetra::CrsMatrix, include
/// "Tpetra_CrsMatrix.hpp" (a file which CMake generates and installs
/// for you).  If you only want the declaration of Tpetra::CrsMatrix,
/// include "Tpetra_CrsMatrix_decl.hpp".

#include "Tpetra_Import_Util.hpp"
#include "Tpetra_Import_Util2.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_LocalCrsMatrixOperator.hpp"

#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_castAwayConstDualView.hpp"
#include "Tpetra_Details_computeOffsets.hpp"
#include "Tpetra_Details_copyOffsets.hpp"
#include "Tpetra_Details_createMirrorView.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Details_getDiagCopyWithoutOffsets.hpp"
#include "Tpetra_Details_leftScaleLocalCrsMatrix.hpp"
#include "Tpetra_Details_Profiling.hpp"
#include "Tpetra_Details_rightScaleLocalCrsMatrix.hpp"
#include "Tpetra_Details_ScalarViewTraits.hpp"
#include "Tpetra_Details_copyConvert.hpp"
#include "Tpetra_Details_iallreduce.hpp"
#include "Tpetra_Details_getEntryOnHost.hpp"
#include "Tpetra_Details_packCrsMatrix.hpp"
#include "Tpetra_Details_unpackCrsMatrixAndCombine.hpp"
#include "Tpetra_Details_crsUtils.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_DataAccess.hpp"
#include "Teuchos_SerialDenseMatrix.hpp" // unused here, could delete
#include "KokkosBlas1_scal.hpp"
#include "KokkosSparse_getDiagCopy.hpp"
#include "KokkosSparse_spmv.hpp"

#include <memory>
#include <sstream>
#include <typeinfo>
#include <utility>
#include <vector>

namespace Tpetra {

namespace { // (anonymous)

  template<class T, class BinaryFunction>
  T atomic_binary_function_update (volatile T* const dest,
                                   const T& inputVal,
                                   BinaryFunction f)
  {
    T oldVal = *dest;
    T assume;

    // NOTE (mfh 30 Nov 2015) I do NOT need a fence here for IBM
    // POWER architectures, because 'newval' depends on 'assume',
    // which depends on 'oldVal', which depends on '*dest'.  This
    // sets up a chain of read dependencies that should ensure
    // correct behavior given a sane memory model.
    do {
      assume = oldVal;
      T newVal = f (assume, inputVal);
      oldVal = Kokkos::atomic_compare_exchange (dest, assume, newVal);
    } while (assume != oldVal);

    return oldVal;
  }
} // namespace (anonymous)

//
// Users must never rely on anything in the Details namespace.
//
namespace Details {

/// \struct AbsMax
/// \brief Functor for the the ABSMAX CombineMode of Import and Export operations.
/// \tparam Scalar Same as the Scalar template parameter of CrsMatrix.
///
/// \warning This is an implementation detail of CrsMatrix.  Users
///   must not rely on this class.  It may disappear or its
///   interface may change at any time.
///
/// \tparam Scalar Same as the Scalar template parameter of CrsMatrix.
template<class Scalar>
struct AbsMax {
  //! Return the maximum of the magnitudes (absolute values) of x and y.
  Scalar operator() (const Scalar& x, const Scalar& y) {
    typedef Teuchos::ScalarTraits<Scalar> STS;
    return std::max (STS::magnitude (x), STS::magnitude (y));
  }
};

} // namespace Details
} // namespace Tpetra

namespace Tpetra {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
             size_t maxNumEntriesPerRow,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
  {
    const char tfecfFuncName[] = "CrsMatrix(RCP<const Map>, size_t "
      "[, RCP<ParameterList>]): ";
    Teuchos::RCP<crs_graph_type> graph;
    try {
      graph = Teuchos::rcp (new crs_graph_type (rowMap, maxNumEntriesPerRow,
                                                params));
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "CrsGraph constructor (RCP<const Map>, "
         "size_t [, RCP<ParameterList>]) threw an exception: "
         << e.what ());
    }
    // myGraph_ not null means that the matrix owns the graph.  That's
    // different than the const CrsGraph constructor, where the matrix
    // does _not_ own the graph.
    myGraph_ = graph;
    staticGraph_ = myGraph_;
    resumeFill (params);
    checkInternalState ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
             const Teuchos::ArrayView<const size_t>& numEntPerRowToAlloc,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
  {
    const char tfecfFuncName[] = "CrsMatrix(RCP<const Map>, "
      "ArrayView<const size_t>[, RCP<ParameterList>]): ";
    Teuchos::RCP<crs_graph_type> graph;
    try {
      using Teuchos::rcp;
      graph = rcp(new crs_graph_type(rowMap, numEntPerRowToAlloc,
                                     params));
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "CrsGraph constructor "
         "(RCP<const Map>, ArrayView<const size_t>"
         "[, RCP<ParameterList>]) threw an exception: "
         << e.what ());
    }
    // myGraph_ not null means that the matrix owns the graph.  That's
    // different than the const CrsGraph constructor, where the matrix
    // does _not_ own the graph.
    myGraph_ = graph;
    staticGraph_ = graph;
    resumeFill (params);
    checkInternalState ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
             const Teuchos::RCP<const map_type>& colMap,
             const size_t maxNumEntPerRow,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
  {
    const char tfecfFuncName[] = "CrsMatrix(RCP<const Map>, "
      "RCP<const Map>, size_t[, RCP<ParameterList>]): ";
    const char suffix[] =
      "  Please report this bug to the Tpetra developers.";

    // An artifact of debugging something a while back.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! staticGraph_.is_null (), std::logic_error,
       "staticGraph_ is not null at the beginning of the constructor."
       << suffix);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! myGraph_.is_null (), std::logic_error,
       "myGraph_ is not null at the beginning of the constructor."
       << suffix);
    Teuchos::RCP<crs_graph_type> graph;
    try {
      graph = Teuchos::rcp (new crs_graph_type (rowMap, colMap,
                                                maxNumEntPerRow,
                                                params));
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "CrsGraph constructor (RCP<const Map>, "
         "RCP<const Map>, size_t[, RCP<ParameterList>]) threw an "
         "exception: " << e.what ());
    }
    // myGraph_ not null means that the matrix owns the graph.  That's
    // different than the const CrsGraph constructor, where the matrix
    // does _not_ own the graph.
    myGraph_ = graph;
    staticGraph_ = myGraph_;
    resumeFill (params);
    checkInternalState ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
             const Teuchos::RCP<const map_type>& colMap,
             const Teuchos::ArrayView<const size_t>& numEntPerRowToAlloc,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
  {
    const char tfecfFuncName[] =
      "CrsMatrix(RCP<const Map>, RCP<const Map>, "
      "ArrayView<const size_t>[, RCP<ParameterList>]): ";
    Teuchos::RCP<crs_graph_type> graph;
    try {
      graph = Teuchos::rcp (new crs_graph_type (rowMap, colMap,
                                                numEntPerRowToAlloc,
                                                params));
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "CrsGraph constructor (RCP<const Map>, "
         "RCP<const Map>, ArrayView<const size_t>[, "
         "RCP<ParameterList>]) threw an exception: " << e.what ());
    }
    // myGraph_ not null means that the matrix owns the graph.  That's
    // different than the const CrsGraph constructor, where the matrix
    // does _not_ own the graph.
    myGraph_ = graph;
    staticGraph_ = graph;
    resumeFill (params);
    checkInternalState ();
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  CrsMatrix (const Teuchos::RCP<const crs_graph_type>& graph,
             const Teuchos::RCP<Teuchos::ParameterList>& /* params */) :
    dist_object_type (graph->getRowMap ()),
    staticGraph_ (graph),
    storageStatus_ (Details::STORAGE_1D_PACKED)
  {
    using std::endl;
    typedef typename local_matrix_device_type::values_type values_type;
    const char tfecfFuncName[] = "CrsMatrix(RCP<const CrsGraph>[, "
      "RCP<ParameterList>]): ";
    const bool verbose = Details::Behavior::verbose("CrsMatrix");

    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsMatrix", "CrsMatrix(graph,params)");
      std::ostringstream os;
      os << *prefix << "Start" << endl;
      std::cerr << os.str ();
    }

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (graph.is_null (), std::runtime_error, "Input graph is null.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! graph->isFillComplete (), std::runtime_error, "Input graph "
       "is not fill complete. You must call fillComplete on the "
       "graph before using it to construct a CrsMatrix.  Note that "
       "calling resumeFill on the graph makes it not fill complete, "
       "even if you had previously called fillComplete.  In that "
       "case, you must call fillComplete on the graph again.");

    // The graph is fill complete, so it is locally indexed and has a
    // fixed structure.  This means we can allocate the (1-D) array of
    // values and build the local matrix right now.  Note that the
    // local matrix's number of columns comes from the column Map, not
    // the domain Map.

    const size_t numEnt = graph->lclIndsPacked_wdv.extent (0);
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Allocate values: " << numEnt << endl;
      std::cerr << os.str ();
    }

    values_type val ("Tpetra::CrsMatrix::values", numEnt);
    valuesPacked_wdv = values_wdv_type(val);
    valuesUnpacked_wdv = valuesPacked_wdv;

    checkInternalState ();

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str ();
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  CrsMatrix(CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& matrix,
            const Teuchos::RCP<const crs_graph_type>& graph,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (graph->getRowMap ()),
    staticGraph_ (graph),
    storageStatus_ (matrix.storageStatus_)
  {
    const char tfecfFuncName[] = "CrsMatrix(RCP<const CrsGraph>, "
      "local_matrix_device_type::values_type, "
      "[,RCP<ParameterList>]): ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (graph.is_null (), std::runtime_error, "Input graph is null.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! graph->isFillComplete (), std::runtime_error, "Input graph "
       "is not fill complete. You must call fillComplete on the "
       "graph before using it to construct a CrsMatrix.  Note that "
       "calling resumeFill on the graph makes it not fill complete, "
       "even if you had previously called fillComplete.  In that "
       "case, you must call fillComplete on the graph again.");

    size_t numValuesPacked = graph->lclIndsPacked_wdv.extent(0);
    valuesPacked_wdv = values_wdv_type(matrix.valuesPacked_wdv, 0, numValuesPacked);

    size_t numValuesUnpacked = graph->lclIndsUnpacked_wdv.extent(0);
    valuesUnpacked_wdv = values_wdv_type(matrix.valuesUnpacked_wdv, 0, numValuesUnpacked);

    checkInternalState();
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  CrsMatrix (const Teuchos::RCP<const crs_graph_type>& graph,
             const typename local_matrix_device_type::values_type& values,
             const Teuchos::RCP<Teuchos::ParameterList>& /* params */) :
    dist_object_type (graph->getRowMap ()),
    staticGraph_ (graph),
    storageStatus_ (Details::STORAGE_1D_PACKED)
  {
    const char tfecfFuncName[] = "CrsMatrix(RCP<const CrsGraph>, "
      "local_matrix_device_type::values_type, "
      "[,RCP<ParameterList>]): ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (graph.is_null (), std::runtime_error, "Input graph is null.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! graph->isFillComplete (), std::runtime_error, "Input graph "
       "is not fill complete. You must call fillComplete on the "
       "graph before using it to construct a CrsMatrix.  Note that "
       "calling resumeFill on the graph makes it not fill complete, "
       "even if you had previously called fillComplete.  In that "
       "case, you must call fillComplete on the graph again.");

    // The graph is fill complete, so it is locally indexed and has a
    // fixed structure.  This means we can allocate the (1-D) array of
    // values and build the local matrix right now.  Note that the
    // local matrix's number of columns comes from the column Map, not
    // the domain Map.

    valuesPacked_wdv = values_wdv_type(values);
    valuesUnpacked_wdv = valuesPacked_wdv;

    // FIXME (22 Jun 2016) I would very much like to get rid of
    // k_values1D_ at some point.  I find it confusing to have all
    // these extra references lying around.
    // KDDKDD ALMOST THERE, MARK!
//    k_values1D_ = valuesUnpacked_wdv.getDeviceView(Access::ReadWrite);

    checkInternalState ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
             const Teuchos::RCP<const map_type>& colMap,
             const typename local_graph_device_type::row_map_type& rowPointers,
             const typename local_graph_device_type::entries_type::non_const_type& columnIndices,
             const typename local_matrix_device_type::values_type& values,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap),
    storageStatus_ (Details::STORAGE_1D_PACKED)
  {
    using Details::getEntryOnHost;
    using Teuchos::RCP;
    using std::endl;
    const char tfecfFuncName[] = "Tpetra::CrsMatrix(RCP<const Map>, "
      "RCP<const Map>, ptr, ind, val[, params]): ";
    const char suffix[] =
      ".  Please report this bug to the Tpetra developers.";
    const bool debug = Details::Behavior::debug("CrsMatrix");
    const bool verbose = Details::Behavior::verbose("CrsMatrix");

    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix(
        "CrsMatrix", "CrsMatrix(rowMap,colMap,ptr,ind,val[,params])");
      std::ostringstream os;
      os << *prefix << "Start" << endl;
      std::cerr << os.str ();
    }

    // Check the user's input.  Note that this might throw only on
    // some processes but not others, causing deadlock.  We prefer
    // deadlock due to exceptions to segfaults, because users can
    // catch exceptions.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (values.extent(0) != columnIndices.extent(0),
       std::invalid_argument, "values.extent(0)=" << values.extent(0)
       << " != columnIndices.extent(0) = " << columnIndices.extent(0)
       << ".");
    if (debug && rowPointers.extent(0) != 0) {
      const size_t numEnt =
        getEntryOnHost(rowPointers, rowPointers.extent(0) - 1);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (numEnt != size_t(columnIndices.extent(0)) ||
         numEnt != size_t(values.extent(0)),
         std::invalid_argument, "Last entry of rowPointers says that "
         "the matrix has " << numEnt << " entr"
         << (numEnt != 1 ? "ies" : "y") << ", but the dimensions of "
         "columnIndices and values don't match this.  "
         "columnIndices.extent(0)=" << columnIndices.extent (0)
         << " and values.extent(0)=" << values.extent (0) << ".");
    }

    RCP<crs_graph_type> graph;
    try {
      graph = Teuchos::rcp (new crs_graph_type (rowMap, colMap, rowPointers,
                                                columnIndices, params));
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "CrsGraph constructor (RCP<const Map>, "
         "RCP<const Map>, ptr, ind[, params]) threw an exception: "
         << e.what ());
    }
    // The newly created CrsGraph _must_ have a local graph at this
    // point.  We don't really care whether CrsGraph's constructor
    // deep-copies or shallow-copies the input, but the dimensions
    // have to be right.  That's how we tell whether the CrsGraph has
    // a local graph.
    auto lclGraph = graph->getLocalGraphDevice ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (lclGraph.row_map.extent (0) != rowPointers.extent (0) ||
       lclGraph.entries.extent (0) != columnIndices.extent (0),
       std::logic_error, "CrsGraph's constructor (rowMap, colMap, ptr, "
       "ind[, params]) did not set the local graph correctly." << suffix);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (lclGraph.entries.extent (0) != values.extent (0),
       std::logic_error, "CrsGraph's constructor (rowMap, colMap, ptr, ind[, "
       "params]) did not set the local graph correctly.  "
       "lclGraph.entries.extent(0) = " << lclGraph.entries.extent (0)
       << " != values.extent(0) = " << values.extent (0) << suffix);

    // myGraph_ not null means that the matrix owns the graph.  This
    // is true because the column indices come in as nonconst,
    // implying shared ownership.
    myGraph_ = graph;
    staticGraph_ = graph;

    // The graph may not be fill complete yet.  However, it is locally
    // indexed (since we have a column Map) and has a fixed structure
    // (due to the input arrays).  This means we can allocate the
    // (1-D) array of values and build the local matrix right now.
    // Note that the local matrix's number of columns comes from the
    // column Map, not the domain Map.

    valuesPacked_wdv = values_wdv_type(values);
    valuesUnpacked_wdv = valuesPacked_wdv;

    // FIXME (22 Jun 2016) I would very much like to get rid of
    // k_values1D_ at some point.  I find it confusing to have all
    // these extra references lying around.
//    this->k_values1D_ = valuesPacked_wdv.getDeviceView(Access::ReadWrite);

    checkInternalState ();
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
             const Teuchos::RCP<const map_type>& colMap,
             const Teuchos::ArrayRCP<size_t>& ptr,
             const Teuchos::ArrayRCP<LocalOrdinal>& ind,
             const Teuchos::ArrayRCP<Scalar>& val,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap),
    storageStatus_ (Details::STORAGE_1D_PACKED)
  {
    using Kokkos::Compat::getKokkosViewDeepCopy;
    using Teuchos::av_reinterpret_cast;
    using Teuchos::RCP;
    using values_type = typename local_matrix_device_type::values_type;
    using IST = impl_scalar_type;
    const char tfecfFuncName[] = "Tpetra::CrsMatrix(RCP<const Map>, "
      "RCP<const Map>, ptr, ind, val[, params]): ";

    RCP<crs_graph_type> graph;
    try {
      graph = Teuchos::rcp (new crs_graph_type (rowMap, colMap, ptr,
                                                ind, params));
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "CrsGraph constructor (RCP<const Map>, "
         "RCP<const Map>, ArrayRCP<size_t>, ArrayRCP<LocalOrdinal>[, "
         "RCP<ParameterList>]) threw an exception: " << e.what ());
    }
    // myGraph_ not null means that the matrix owns the graph.  This
    // is true because the column indices come in as nonconst,
    // implying shared ownership.
    myGraph_ = graph;
    staticGraph_ = graph;

    // The graph may not be fill complete yet.  However, it is locally
    // indexed (since we have a column Map) and has a fixed structure
    // (due to the input arrays).  This means we can allocate the
    // (1-D) array of values and build the local matrix right now.
    // Note that the local matrix's number of columns comes from the
    // column Map, not the domain Map.

    // The graph _must_ have a local graph at this point.  We don't
    // really care whether CrsGraph's constructor deep-copies or
    // shallow-copies the input, but the dimensions have to be right.
    // That's how we tell whether the CrsGraph has a local graph.
    auto lclGraph = staticGraph_->getLocalGraphDevice ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (size_t (lclGraph.row_map.extent (0)) != size_t (ptr.size ()) ||
       size_t (lclGraph.entries.extent (0)) != size_t (ind.size ()),
       std::logic_error, "CrsGraph's constructor (rowMap, colMap, "
       "ptr, ind[, params]) did not set the local graph correctly.  "
       "Please report this bug to the Tpetra developers.");

    values_type valIn =
      getKokkosViewDeepCopy<device_type> (av_reinterpret_cast<IST> (val ()));
    valuesPacked_wdv = values_wdv_type(valIn);
    valuesUnpacked_wdv = valuesPacked_wdv;

    // FIXME (22 Jun 2016) I would very much like to get rid of
    // k_values1D_ at some point.  I find it confusing to have all
    // these extra references lying around.
//    this->k_values1D_ = valuesPacked_wdv.getDeviceView(Access::ReadWrite);

    checkInternalState ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
             const Teuchos::RCP<const map_type>& colMap,
             const local_matrix_device_type& lclMatrix,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap),
    storageStatus_ (Details::STORAGE_1D_PACKED),
    fillComplete_ (true)
  {
    const char tfecfFuncName[] = "Tpetra::CrsMatrix(RCP<const Map>, "
      "RCP<const Map>, local_matrix_device_type[, RCP<ParameterList>]): ";
    const char suffix[] =
      "  Please report this bug to the Tpetra developers.";

    Teuchos::RCP<crs_graph_type> graph;
    try {
      graph = Teuchos::rcp (new crs_graph_type (rowMap, colMap,
                                                lclMatrix.graph, params));
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "CrsGraph constructor (RCP<const Map>, "
         "RCP<const Map>, local_graph_device_type[, RCP<ParameterList>]) threw an "
         "exception: " << e.what ());
    }
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (!graph->isFillComplete (), std::logic_error, "CrsGraph constructor (RCP"
       "<const Map>, RCP<const Map>, local_graph_device_type[, RCP<ParameterList>]) "
       "did not produce a fill-complete graph.  Please report this bug to the "
       "Tpetra developers.");
    // myGraph_ not null means that the matrix owns the graph.  This
    // is true because the column indices come in as nonconst through
    // the matrix, implying shared ownership.
    myGraph_ = graph;
    staticGraph_ = graph;

    valuesPacked_wdv = values_wdv_type(lclMatrix.values);
    valuesUnpacked_wdv = valuesPacked_wdv;

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (isFillActive (), std::logic_error,
       "At the end of a CrsMatrix constructor that should produce "
       "a fillComplete matrix, isFillActive() is true." << suffix);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! isFillComplete (), std::logic_error, "At the end of a "
       "CrsMatrix constructor that should produce a fillComplete "
       "matrix, isFillComplete() is false." << suffix);
    checkInternalState ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  CrsMatrix (const local_matrix_device_type& lclMatrix,
             const Teuchos::RCP<const map_type>& rowMap,
             const Teuchos::RCP<const map_type>& colMap,
             const Teuchos::RCP<const map_type>& domainMap,
             const Teuchos::RCP<const map_type>& rangeMap,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap),
    storageStatus_ (Details::STORAGE_1D_PACKED),
    fillComplete_ (true)
  {
    const char tfecfFuncName[] = "Tpetra::CrsMatrix(RCP<const Map>, "
      "RCP<const Map>, RCP<const Map>, RCP<const Map>, "
      "local_matrix_device_type[, RCP<ParameterList>]): ";
    const char suffix[] =
      "  Please report this bug to the Tpetra developers.";

    Teuchos::RCP<crs_graph_type> graph;
    try {
      graph = Teuchos::rcp (new crs_graph_type (lclMatrix.graph, rowMap, colMap,
                                                domainMap, rangeMap, params));
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "CrsGraph constructor (RCP<const Map>, "
         "RCP<const Map>, RCP<const Map>, RCP<const Map>, local_graph_device_type[, "
         "RCP<ParameterList>]) threw an exception: " << e.what ());
    }
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! graph->isFillComplete (), std::logic_error, "CrsGraph "
       "constructor (RCP<const Map>, RCP<const Map>, RCP<const Map>, "
       "RCP<const Map>, local_graph_device_type[, RCP<ParameterList>]) did "
       "not produce a fillComplete graph." << suffix);
    // myGraph_ not null means that the matrix owns the graph.  This
    // is true because the column indices come in as nonconst through
    // the matrix, implying shared ownership.
    myGraph_ = graph;
    staticGraph_ = graph;

    valuesPacked_wdv = values_wdv_type(lclMatrix.values);
    valuesUnpacked_wdv = valuesPacked_wdv;

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (isFillActive (), std::logic_error,
       "At the end of a CrsMatrix constructor that should produce "
       "a fillComplete matrix, isFillActive() is true." << suffix);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! isFillComplete (), std::logic_error, "At the end of a "
       "CrsMatrix constructor that should produce a fillComplete "
       "matrix, isFillComplete() is false." << suffix);
    checkInternalState ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  CrsMatrix (const local_matrix_device_type& lclMatrix,
             const Teuchos::RCP<const map_type>& rowMap,
             const Teuchos::RCP<const map_type>& colMap,
             const Teuchos::RCP<const map_type>& domainMap,
             const Teuchos::RCP<const map_type>& rangeMap,
             const Teuchos::RCP<const import_type>& importer,
             const Teuchos::RCP<const export_type>& exporter,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap),
    storageStatus_ (Details::STORAGE_1D_PACKED),
    fillComplete_ (true)
  {
    using Teuchos::rcp;
    const char tfecfFuncName[] = "Tpetra::CrsMatrix"
      "(lclMat,Map,Map,Map,Map,Import,Export,params): ";
    const char suffix[] =
      "  Please report this bug to the Tpetra developers.";

    Teuchos::RCP<crs_graph_type> graph;
    try {
      graph = rcp (new crs_graph_type (lclMatrix.graph, rowMap, colMap,
                                       domainMap, rangeMap, importer,
                                       exporter, params));
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "CrsGraph constructor "
         "(local_graph_device_type, Map, Map, Map, Map, Import, Export, "
         "params) threw: " << e.what ());
    }
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (!graph->isFillComplete (), std::logic_error, "CrsGraph "
       "constructor (local_graph_device_type, Map, Map, Map, Map, Import, "
       "Export, params) did not produce a fill-complete graph.  "
       "Please report this bug to the Tpetra developers.");
    // myGraph_ not null means that the matrix owns the graph.  This
    // is true because the column indices come in as nonconst through
    // the matrix, implying shared ownership.
    myGraph_ = graph;
    staticGraph_ = graph;

    valuesPacked_wdv = values_wdv_type(lclMatrix.values);
    valuesUnpacked_wdv = valuesPacked_wdv;

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (isFillActive (), std::logic_error,
       "At the end of a CrsMatrix constructor that should produce "
       "a fillComplete matrix, isFillActive() is true." << suffix);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! isFillComplete (), std::logic_error, "At the end of a "
       "CrsMatrix constructor that should produce a fillComplete "
       "matrix, isFillComplete() is false." << suffix);
    checkInternalState ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  CrsMatrix (const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source,
             const Teuchos::DataAccess copyOrView):
    dist_object_type (source.getCrsGraph()->getRowMap ()),
    staticGraph_ (source.getCrsGraph()),
    storageStatus_ (source.storageStatus_)
  {
    const char tfecfFuncName[] = "Tpetra::CrsMatrix("
      "const CrsMatrix&, const Teuchos::DataAccess): ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! source.isFillComplete (), std::invalid_argument,
       "Source graph must be fillComplete().");

    if (copyOrView == Teuchos::Copy) {
      using values_type = typename local_matrix_device_type::values_type;
      auto vals = source.getLocalValuesDevice (Access::ReadOnly);
      using Kokkos::view_alloc;
      using Kokkos::WithoutInitializing;
      values_type newvals (view_alloc ("val", WithoutInitializing),
                           vals.extent (0));
          // DEEP_COPY REVIEW - DEVICE-TO_DEVICE
      Kokkos::deep_copy (newvals, vals);
      valuesPacked_wdv = values_wdv_type(newvals);
      valuesUnpacked_wdv = valuesPacked_wdv;
      fillComplete (source.getDomainMap (), source.getRangeMap ());
    }
    else if (copyOrView == Teuchos::View) {
      valuesPacked_wdv = values_wdv_type(source.valuesPacked_wdv);
      valuesUnpacked_wdv = values_wdv_type(source.valuesUnpacked_wdv);
      fillComplete (source.getDomainMap (), source.getRangeMap ());
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::invalid_argument, "Second argument 'copyOrView' "
         "has an invalid value " << copyOrView << ".  Valid values "
         "include Teuchos::Copy = " << Teuchos::Copy << " and "
         "Teuchos::View = " << Teuchos::View << ".");
    }
    checkInternalState();
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  swap(CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & crs_matrix)
  {
    std::swap(crs_matrix.importMV_,      this->importMV_);
    std::swap(crs_matrix.exportMV_,      this->exportMV_);
    std::swap(crs_matrix.staticGraph_,   this->staticGraph_);
    std::swap(crs_matrix.myGraph_,       this->myGraph_);
    std::swap(crs_matrix.valuesPacked_wdv, this->valuesPacked_wdv);
    std::swap(crs_matrix.valuesUnpacked_wdv, this->valuesUnpacked_wdv);
    std::swap(crs_matrix.storageStatus_, this->storageStatus_);
    std::swap(crs_matrix.fillComplete_,  this->fillComplete_);
    std::swap(crs_matrix.nonlocals_,     this->nonlocals_);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Teuchos::Comm<int> >
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getComm () const {
    return getCrsGraphRef ().getComm ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  isFillComplete () const {
    return fillComplete_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  isFillActive () const {
    return ! fillComplete_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  isStorageOptimized () const {
    return this->getCrsGraphRef ().isStorageOptimized ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  isLocallyIndexed () const {
    return getCrsGraphRef ().isLocallyIndexed ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  isGloballyIndexed () const {
    return getCrsGraphRef ().isGloballyIndexed ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  hasColMap () const {
    return getCrsGraphRef ().hasColMap ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getGlobalNumEntries () const {
    return getCrsGraphRef ().getGlobalNumEntries ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalNumEntries () const {
    return getCrsGraphRef ().getLocalNumEntries ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getGlobalNumRows () const {
    return getCrsGraphRef ().getGlobalNumRows ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getGlobalNumCols () const {
    return getCrsGraphRef ().getGlobalNumCols ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalNumRows () const {
    return getCrsGraphRef ().getLocalNumRows ();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalNumCols () const {
    return getCrsGraphRef ().getLocalNumCols ();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getNumEntriesInGlobalRow (GlobalOrdinal globalRow) const {
    return getCrsGraphRef ().getNumEntriesInGlobalRow (globalRow);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getNumEntriesInLocalRow (LocalOrdinal localRow) const {
    return getCrsGraphRef ().getNumEntriesInLocalRow (localRow);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getGlobalMaxNumRowEntries () const {
    return getCrsGraphRef ().getGlobalMaxNumRowEntries ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalMaxNumRowEntries () const {
    return getCrsGraphRef ().getLocalMaxNumRowEntries ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  GlobalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getIndexBase () const {
    return getRowMap ()->getIndexBase ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getRowMap () const {
    return getCrsGraphRef ().getRowMap ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getColMap () const {
    return getCrsGraphRef ().getColMap ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getDomainMap () const {
    return getCrsGraphRef ().getDomainMap ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getRangeMap () const {
    return getCrsGraphRef ().getRangeMap ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const RowGraph<LocalOrdinal, GlobalOrdinal, Node> >
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getGraph () const {
    if (staticGraph_ != Teuchos::null) {
      return staticGraph_;
    }
    return myGraph_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getCrsGraph () const {
    if (staticGraph_ != Teuchos::null) {
      return staticGraph_;
    }
    return myGraph_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const CrsGraph<LocalOrdinal, GlobalOrdinal, Node>&
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getCrsGraphRef () const
  {
#ifdef HAVE_TPETRA_DEBUG
    constexpr bool debug = true;
#else
    constexpr bool debug = false;
#endif // HAVE_TPETRA_DEBUG

    if (! this->staticGraph_.is_null ()) {
      return * (this->staticGraph_);
    }
    else {
      if (debug) {
        const char tfecfFuncName[] = "getCrsGraphRef: ";
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (this->myGraph_.is_null (), std::logic_error,
           "Both staticGraph_ and myGraph_ are null.  "
           "Please report this bug to the Tpetra developers.");
      }
      return * (this->myGraph_);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_device_type
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalMatrixDevice () const
  {
    auto numCols = staticGraph_->getColMap()->getLocalNumElements();
    return local_matrix_device_type("Tpetra::CrsMatrix::lclMatrixDevice",
                              numCols,
                              valuesPacked_wdv.getDeviceView(Access::ReadWrite),
                              staticGraph_->getLocalGraphDevice());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_host_type
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalMatrixHost () const
  {
    auto numCols = staticGraph_->getColMap()->getLocalNumElements();
    return local_matrix_host_type("Tpetra::CrsMatrix::lclMatrixHost", numCols,
                                valuesPacked_wdv.getHostView(Access::ReadWrite),
                                staticGraph_->getLocalGraphHost());
  }

#if KOKKOSKERNELS_VERSION < 40299
// KDDKDD NOT SURE WHY THIS MUST RETURN A SHARED_PTR
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::shared_ptr<typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_multiply_op_type>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalMultiplyOperator () const
  {
    auto localMatrix = getLocalMatrixDevice();
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) || defined(KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE) || defined(KOKKOSKERNELS_ENABLE_TPL_MKL)
    if(this->getLocalNumEntries() <= size_t(Teuchos::OrdinalTraits<LocalOrdinal>::max()))
    {
      if(this->ordinalRowptrs.data() == nullptr)
      {
        auto originalRowptrs = localMatrix.graph.row_map;
        //create LocalOrdinal-typed copy of the local graph's rowptrs.
        //This enables the LocalCrsMatrixOperator to use cuSPARSE SpMV.
        this->ordinalRowptrs = ordinal_rowptrs_type(
            Kokkos::ViewAllocateWithoutInitializing("CrsMatrix::ordinalRowptrs"), originalRowptrs.extent(0));
        auto ordinalRowptrs_ = this->ordinalRowptrs;  //don't want to capture 'this'
        Kokkos::parallel_for("CrsMatrix::getLocalMultiplyOperator::convertRowptrs",
            Kokkos::RangePolicy<execution_space>(0, originalRowptrs.extent(0)),
            KOKKOS_LAMBDA(LocalOrdinal i)
            {
              ordinalRowptrs_(i) = originalRowptrs(i);
            });
      }
      //return local operator using ordinalRowptrs
      return std::make_shared<local_multiply_op_type>(
          std::make_shared<local_matrix_device_type>(localMatrix), this->ordinalRowptrs);
    }
#endif
// KDDKDD NOT SURE WHY THIS MUST RETURN A SHARED_PTR
    return std::make_shared<local_multiply_op_type>(
                           std::make_shared<local_matrix_device_type>(localMatrix));
  }
#endif

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  isStaticGraph () const {
    return myGraph_.is_null ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  hasTransposeApply () const {
    return true;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  supportsRowViews () const {
    return true;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  allocateValues (ELocalGlobal lg, GraphAllocationStatus gas,
                  const bool verbose)
  {
    using Details::Behavior;
    using Details::ProfilingRegion;
    using std::endl;
    const char tfecfFuncName[] = "allocateValues: ";
    const char suffix[] =
      "  Please report this bug to the Tpetra developers.";
    ProfilingRegion region("Tpetra::CrsMatrix::allocateValues");

    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsMatrix", "allocateValues");
      std::ostringstream os;
      os << *prefix << "lg: "
         << (lg == LocalIndices ? "Local" : "Global") << "Indices"
         << ", gas: Graph"
         << (gas == GraphAlreadyAllocated ? "Already" : "NotYet")
         << "Allocated" << endl;
      std::cerr << os.str();
    }

    const bool debug = Behavior::debug("CrsMatrix");
    if (debug) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->staticGraph_.is_null (), std::logic_error,
         "staticGraph_ is null." << suffix);

      // If the graph indices are already allocated, then gas should be
      // GraphAlreadyAllocated.  Otherwise, gas should be
      // GraphNotYetAllocated.
      if ((gas == GraphAlreadyAllocated) !=
          staticGraph_->indicesAreAllocated ()) {
        const char err1[] = "The caller has asserted that the graph "
          "is ";
        const char err2[] = "already allocated, but the static graph "
          "says that its indices are ";
        const char err3[] = "already allocated.  ";
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (gas == GraphAlreadyAllocated &&
           ! staticGraph_->indicesAreAllocated (), std::logic_error,
           err1 << err2 << "not " << err3 << suffix);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (gas != GraphAlreadyAllocated &&
           staticGraph_->indicesAreAllocated (), std::logic_error,
           err1 << "not " << err2 << err3 << suffix);
      }

      // If the graph is unallocated, then it had better be a
      // matrix-owned graph.  ("Matrix-owned graph" means that the
      // matrix gets to define the graph structure.  If the CrsMatrix
      // constructor that takes an RCP<const CrsGraph> was used, then
      // the matrix does _not_ own the graph.)
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! this->staticGraph_->indicesAreAllocated () &&
         this->myGraph_.is_null (), std::logic_error,
         "The static graph says that its indices are not allocated, "
         "but the graph is not owned by the matrix." << suffix);
    }

    if (gas == GraphNotYetAllocated) {
      if (debug) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (this->myGraph_.is_null (), std::logic_error,
           "gas = GraphNotYetAllocated, but myGraph_ is null." << suffix);
      }
      try {
        this->myGraph_->allocateIndices (lg, verbose);
      }
      catch (std::exception& e) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (true, std::runtime_error, "CrsGraph::allocateIndices "
           "threw an exception: " << e.what ());
      }
      catch (...) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (true, std::runtime_error, "CrsGraph::allocateIndices "
           "threw an exception not a subclass of std::exception.");
      }
    }

    // Allocate matrix values.
    const size_t lclTotalNumEntries = this->staticGraph_->getLocalAllocationSize();
    if (debug) {
      const size_t lclNumRows = this->staticGraph_->getLocalNumRows ();
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->staticGraph_->getRowPtrsUnpackedHost()(lclNumRows) != lclTotalNumEntries, std::logic_error,
         "length of staticGraph's lclIndsUnpacked does not match final entry of rowPtrsUnapcked_host." << suffix);
    }

    // Allocate array of (packed???) matrix values.
    using values_type = typename local_matrix_device_type::values_type;
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Allocate values_wdv: Pre "
         << valuesUnpacked_wdv.extent(0) << ", post "
         << lclTotalNumEntries << endl;
      std::cerr << os.str();
    }
//    this->k_values1D_ =
    valuesUnpacked_wdv = values_wdv_type(
                                    values_type("Tpetra::CrsMatrix::values",
                                    lclTotalNumEntries));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  fillLocalGraphAndMatrix (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    using ::Tpetra::Details::computeOffsetsFromCounts;
    using ::Tpetra::Details::getEntryOnHost;
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    using row_map_type = typename local_graph_device_type::row_map_type;
    using lclinds_1d_type = typename Graph::local_graph_device_type::entries_type::non_const_type;
    using values_type = typename local_matrix_device_type::values_type;
    Details::ProfilingRegion regionFLGAM
      ("Tpetra::CrsMatrix::fillLocalGraphAndMatrix");

    const char tfecfFuncName[] = "fillLocalGraphAndMatrix (called from "
      "fillComplete or expertStaticFillComplete): ";
    const char suffix[] =
      "  Please report this bug to the Tpetra developers.";
    const bool debug = Details::Behavior::debug("CrsMatrix");
    const bool verbose = Details::Behavior::verbose("CrsMatrix");

    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsMatrix", "fillLocalGraphAndMatrix");
      std::ostringstream os;
      os << *prefix << endl;
      std::cerr << os.str ();
    }

    if (debug) {
      // fillComplete() only calls fillLocalGraphAndMatrix() if the
      // matrix owns the graph, which means myGraph_ is not null.
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (myGraph_.is_null (), std::logic_error, "The nonconst graph "
         "(myGraph_) is null.  This means that the matrix has a "
         "const (a.k.a. \"static\") graph.  fillComplete or "
         "expertStaticFillComplete should never call "
         "fillLocalGraphAndMatrix in that case." << suffix);
    }

    const size_t lclNumRows = this->getLocalNumRows ();

    // This method's goal is to fill in the three arrays (compressed
    // sparse row format) that define the sparse graph's and matrix's
    // structure, and the sparse matrix's values.
    //
    // Get references to the data in myGraph_, so we can modify them
    // as well.  Note that we only call fillLocalGraphAndMatrix() if
    // the matrix owns the graph, which means myGraph_ is not null.

    // NOTE: This does not work correctly w/ GCC 12.3 + CUDA due to a compiler bug.
    // See: https://github.com/trilinos/Trilinos/issues/12237
    //using row_entries_type = decltype (myGraph_->k_numRowEntries_); 
    using row_entries_type = typename crs_graph_type::num_row_entries_type;

    typename Graph::local_graph_device_type::row_map_type curRowOffsets = 
                                                   myGraph_->rowPtrsUnpacked_dev_;

    if (debug) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (curRowOffsets.extent (0) == 0, std::logic_error,
         "curRowOffsets.extent(0) == 0.");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (curRowOffsets.extent (0) != lclNumRows + 1, std::logic_error,
         "curRowOffsets.extent(0) = "
         << curRowOffsets.extent (0) << " != lclNumRows + 1 = "
         << (lclNumRows + 1) << ".");
      const size_t numOffsets = curRowOffsets.extent (0);
      const auto valToCheck = myGraph_->getRowPtrsUnpackedHost()(numOffsets - 1);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (numOffsets != 0 &&
         myGraph_->lclIndsUnpacked_wdv.extent (0) != valToCheck,
         std::logic_error, "numOffsets = " <<
         numOffsets << " != 0 and myGraph_->lclIndsUnpacked_wdv.extent(0) = "
         << myGraph_->lclIndsUnpacked_wdv.extent (0) << " != curRowOffsets("
         << numOffsets << ") = " << valToCheck << ".");
    }

    if (myGraph_->getLocalNumEntries() !=
        myGraph_->getLocalAllocationSize()) {

      // Use the nonconst version of row_map_type for k_ptrs,
      // because row_map_type is const and we need to modify k_ptrs here.
      typename row_map_type::non_const_type k_ptrs;
      row_map_type k_ptrs_const;
      lclinds_1d_type k_inds;
      values_type k_vals;

      if (verbose) {
        std::ostringstream os;
        const auto numEnt = myGraph_->getLocalNumEntries();
        const auto allocSize = myGraph_->getLocalAllocationSize();
        os << *prefix << "Unpacked 1-D storage: numEnt=" << numEnt
           << ", allocSize=" << allocSize << endl;
        std::cerr << os.str ();
      }
      // The matrix's current 1-D storage is "unpacked."  This means
      // the row offsets may differ from what the final row offsets
      // should be.  This could happen, for example, if the user
      // set an upper
      // bound on the number of entries per row, but didn't fill all
      // those entries.
      if (debug && curRowOffsets.extent (0) != 0) {
        const size_t numOffsets =
          static_cast<size_t> (curRowOffsets.extent (0));
        const auto valToCheck = myGraph_->getRowPtrsUnpackedHost()(numOffsets - 1);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (static_cast<size_t> (valToCheck) !=
           static_cast<size_t> (valuesUnpacked_wdv.extent (0)),
           std::logic_error, "(unpacked branch) Before "
           "allocating or packing, curRowOffsets(" << (numOffsets-1)
           << ") = " << valToCheck << " != valuesUnpacked_wdv.extent(0)"
           " = " << valuesUnpacked_wdv.extent (0) << ".");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (static_cast<size_t> (valToCheck) !=
           static_cast<size_t> (myGraph_->lclIndsUnpacked_wdv.extent (0)),
           std::logic_error, "(unpacked branch) Before "
           "allocating or packing, curRowOffsets(" << (numOffsets-1)
           << ") = " << valToCheck
           << " != myGraph_->lclIndsUnpacked_wdv.extent(0) = "
           << myGraph_->lclIndsUnpacked_wdv.extent (0) << ".");
      }
      // Pack the row offsets into k_ptrs, by doing a sum-scan of
      // the array of valid entry counts per row.

      // Total number of entries in the matrix on the calling
      // process.  We will compute this in the loop below.  It's
      // cheap to compute and useful as a sanity check.
      size_t lclTotalNumEntries = 0;
      {
        // Allocate the packed row offsets array.  We use a nonconst
        // temporary (packedRowOffsets) here, because k_ptrs is
        // const.  We will assign packedRowOffsets to k_ptrs below.
        if (verbose) {
          std::ostringstream os;
          os << *prefix << "Allocate packed row offsets: "
             << (lclNumRows+1) << endl;
          std::cerr << os.str ();
        }
        typename row_map_type::non_const_type
          packedRowOffsets ("Tpetra::CrsGraph::ptr", lclNumRows + 1);
        typename row_entries_type::const_type numRowEnt_h =
          myGraph_->k_numRowEntries_;
        // We're computing offsets on device.  This function can
        // handle numRowEnt_h being a host View.
        lclTotalNumEntries =
          computeOffsetsFromCounts (packedRowOffsets, numRowEnt_h);
        // packedRowOffsets is modifiable; k_ptrs isn't, so we have
        // to use packedRowOffsets in the loop above and assign here.
        k_ptrs = packedRowOffsets;
        k_ptrs_const = k_ptrs;
      }

      if (debug) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (static_cast<size_t> (k_ptrs.extent (0)) != lclNumRows + 1,
           std::logic_error,
           "(unpacked branch) After packing k_ptrs, "
           "k_ptrs.extent(0) = " << k_ptrs.extent (0) << " != "
           "lclNumRows+1 = " << (lclNumRows+1) << ".");
        const auto valToCheck = getEntryOnHost (k_ptrs, lclNumRows);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (valToCheck != lclTotalNumEntries, std::logic_error,
           "(unpacked branch) After filling k_ptrs, "
           "k_ptrs(lclNumRows=" << lclNumRows << ") = " << valToCheck
           << " != total number of entries on the calling process = "
           << lclTotalNumEntries << ".");
      }

      // Allocate the arrays of packed column indices and values.
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Allocate packed local column indices: "
           << lclTotalNumEntries << endl;
        std::cerr << os.str ();
      }
      k_inds = lclinds_1d_type ("Tpetra::CrsGraph::lclInds", lclTotalNumEntries);
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Allocate packed values: "
           << lclTotalNumEntries << endl;
        std::cerr << os.str ();
      }
      k_vals = values_type ("Tpetra::CrsMatrix::values", lclTotalNumEntries);

      // curRowOffsets (myGraph_->rowPtrsUnpacked_) (???), lclIndsUnpacked_wdv,
      // and valuesUnpacked_wdv are currently unpacked.  Pack them, using
      // the packed row offsets array k_ptrs that we created above.
      //
      // FIXME (mfh 06 Aug 2014) If "Optimize Storage" is false, we
      // need to keep around the unpacked row offsets, column
      // indices, and values arrays.

      // Pack the column indices from unpacked lclIndsUnpacked_wdv into
      // packed k_inds.  We will replace lclIndsUnpacked_wdv below.
      using inds_packer_type = pack_functor<
        typename Graph::local_graph_device_type::entries_type::non_const_type,
        typename Graph::local_inds_dualv_type::t_dev::const_type,
        typename Graph::local_graph_device_type::row_map_type::non_const_type,
        typename Graph::local_graph_device_type::row_map_type>;
      inds_packer_type indsPacker (
                  k_inds,
                  myGraph_->lclIndsUnpacked_wdv.getDeviceView(Access::ReadOnly),
                  k_ptrs, curRowOffsets);
      using exec_space = typename decltype (k_inds)::execution_space;
      using range_type = Kokkos::RangePolicy<exec_space, LocalOrdinal>;
      Kokkos::parallel_for
        ("Tpetra::CrsMatrix pack column indices",
         range_type (0, lclNumRows), indsPacker);

      // Pack the values from unpacked valuesUnpacked_wdv into packed
      // k_vals.  We will replace valuesPacked_wdv below.
      using vals_packer_type = pack_functor<
        typename values_type::non_const_type,
        typename values_type::const_type, 
        typename row_map_type::non_const_type, 
        typename row_map_type::const_type>;
      vals_packer_type valsPacker (
                       k_vals,
                       this->valuesUnpacked_wdv.getDeviceView(Access::ReadOnly),
                       k_ptrs, curRowOffsets);
      Kokkos::parallel_for ("Tpetra::CrsMatrix pack values",
                            range_type (0, lclNumRows), valsPacker);

      if (debug) {
        const char myPrefix[] = "(\"Optimize Storage\""
          "=true branch) After packing, ";
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (k_ptrs.extent (0) == 0, std::logic_error, myPrefix
           << "k_ptrs.extent(0) = 0.  This probably means that "
           "rowPtrsUnpacked_ was never allocated.");
        if (k_ptrs.extent (0) != 0) {
          const size_t numOffsets (k_ptrs.extent (0));
          const auto valToCheck =
            getEntryOnHost (k_ptrs, numOffsets - 1);
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (size_t (valToCheck) != k_vals.extent (0),
             std::logic_error, myPrefix <<
             "k_ptrs(" << (numOffsets-1) << ") = " << valToCheck <<
             " != k_vals.extent(0) = " << k_vals.extent (0) << ".");
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (size_t (valToCheck) != k_inds.extent (0),
             std::logic_error, myPrefix <<
             "k_ptrs(" << (numOffsets-1) << ") = " << valToCheck <<
             " != k_inds.extent(0) = " << k_inds.extent (0) << ".");
        }
      }
      // Build the local graph.
      myGraph_->setRowPtrsPacked(k_ptrs_const);
      myGraph_->lclIndsPacked_wdv = 
                typename crs_graph_type::local_inds_wdv_type(k_inds);
      valuesPacked_wdv = values_wdv_type(k_vals);
    }
    else { // We don't have to pack, so just set the pointers.
      // FIXME KDDKDD https://github.com/trilinos/Trilinos/issues/9657
      // FIXME? This is already done in the graph fill call - need to avoid the memcpy to host
      myGraph_->rowPtrsPacked_dev_ = myGraph_->rowPtrsUnpacked_dev_;
      myGraph_->rowPtrsPacked_host_ = myGraph_->rowPtrsUnpacked_host_;
      myGraph_->packedUnpackedRowPtrsMatch_ = true;
      myGraph_->lclIndsPacked_wdv = myGraph_->lclIndsUnpacked_wdv;
      valuesPacked_wdv = valuesUnpacked_wdv;

      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Storage already packed: rowPtrsUnpacked_: "
           << myGraph_->getRowPtrsUnpackedHost().extent(0) << ", lclIndsUnpacked_wdv: "
           << myGraph_->lclIndsUnpacked_wdv.extent(0) << ", valuesUnpacked_wdv: "
           << valuesUnpacked_wdv.extent(0) << endl;
        std::cerr << os.str();
      }

      if (debug) {
        const char myPrefix[] =
          "(\"Optimize Storage\"=false branch) ";
        auto rowPtrsUnpackedHost = myGraph_->getRowPtrsUnpackedHost();
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (myGraph_->rowPtrsUnpacked_dev_.extent (0) == 0, std::logic_error, myPrefix
           << "myGraph->rowPtrsUnpacked_dev_.extent(0) = 0.  This probably means "
           "that rowPtrsUnpacked_ was never allocated.");
        if (myGraph_->rowPtrsUnpacked_dev_.extent (0) != 0) {
          const size_t numOffsets = rowPtrsUnpackedHost.extent (0);
          const auto valToCheck = rowPtrsUnpackedHost(numOffsets - 1);
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (size_t (valToCheck) != valuesPacked_wdv.extent (0),
             std::logic_error, myPrefix <<
             "k_ptrs_const(" << (numOffsets-1) << ") = " << valToCheck
             << " != valuesPacked_wdv.extent(0) = " 
             << valuesPacked_wdv.extent (0) << ".");
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (size_t (valToCheck) != myGraph_->lclIndsPacked_wdv.extent (0),
             std::logic_error, myPrefix <<
             "k_ptrs_const(" << (numOffsets-1) << ") = " << valToCheck
             << " != myGraph_->lclIndsPacked.extent(0) = " 
             << myGraph_->lclIndsPacked_wdv.extent (0) << ".");
        }
      }
    }

    if (debug) {
      const char myPrefix[] = "After packing, ";
      auto rowPtrsPackedHost = myGraph_->getRowPtrsPackedHost();
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (size_t (rowPtrsPackedHost.extent (0)) != size_t (lclNumRows + 1),
         std::logic_error, myPrefix << "myGraph_->rowPtrsPacked_host_.extent(0) = "
         << rowPtrsPackedHost.extent (0) << " != lclNumRows+1 = " <<
         (lclNumRows+1) << ".");
      if (rowPtrsPackedHost.extent (0) != 0) {
        const size_t numOffsets (rowPtrsPackedHost.extent (0));
        const size_t valToCheck = rowPtrsPackedHost(numOffsets-1);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (valToCheck != size_t (valuesPacked_wdv.extent (0)),
           std::logic_error, myPrefix << "k_ptrs_const(" <<
           (numOffsets-1) << ") = " << valToCheck
           << " != valuesPacked_wdv.extent(0) = " 
           << valuesPacked_wdv.extent (0) << ".");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (valToCheck != size_t (myGraph_->lclIndsPacked_wdv.extent (0)),
           std::logic_error, myPrefix << "k_ptrs_const(" <<
           (numOffsets-1) << ") = " << valToCheck
           << " != myGraph_->lclIndsPacked_wdvk_inds.extent(0) = " 
           << myGraph_->lclIndsPacked_wdv.extent (0) << ".");
      }
    }

    // May we ditch the old allocations for the packed (and otherwise
    // "optimized") allocations, later in this routine?  Optimize
    // storage if the graph is not static, or if the graph already has
    // optimized storage.
    const bool defaultOptStorage =
      ! isStaticGraph () || staticGraph_->isStorageOptimized ();
    const bool requestOptimizedStorage =
      (! params.is_null () &&
       params->get ("Optimize Storage", defaultOptStorage)) ||
      (params.is_null () && defaultOptStorage);

    // The graph has optimized storage when indices are allocated,
    // myGraph_->k_numRowEntries_ is empty, and there are more than
    // zero rows on this process.  
    if (requestOptimizedStorage) {
      // Free the old, unpacked, unoptimized allocations.
      // Free graph data structures that are only needed for
      // unpacked 1-D storage.
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Optimizing storage: free k_numRowEntries_: "
           << myGraph_->k_numRowEntries_.extent(0) << endl;
        std::cerr << os.str();
      }

      myGraph_->k_numRowEntries_ = row_entries_type ();

      // Keep the new 1-D packed allocations.
      // FIXME KDDKDD https://github.com/trilinos/Trilinos/issues/9657
      // We directly set the memory spaces to avoid a memcpy from device to host
      myGraph_->rowPtrsUnpacked_dev_ = myGraph_->rowPtrsPacked_dev_;
      myGraph_->rowPtrsUnpacked_host_ = myGraph_->rowPtrsPacked_host_;
      myGraph_->packedUnpackedRowPtrsMatch_ = true;
      myGraph_->lclIndsUnpacked_wdv = myGraph_->lclIndsPacked_wdv;
      valuesUnpacked_wdv = valuesPacked_wdv;

      myGraph_->storageStatus_ = Details::STORAGE_1D_PACKED;
      this->storageStatus_ = Details::STORAGE_1D_PACKED;
    }
    else {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "User requested NOT to optimize storage"
           << endl;
        std::cerr << os.str();
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  fillLocalMatrix (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    using ::Tpetra::Details::ProfilingRegion;
    using Teuchos::ArrayRCP;
    using Teuchos::Array;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    using row_map_type = typename Graph::local_graph_device_type::row_map_type;
    using non_const_row_map_type = typename row_map_type::non_const_type;
    using values_type = typename local_matrix_device_type::values_type;
    ProfilingRegion regionFLM("Tpetra::CrsMatrix::fillLocalMatrix");
    const size_t lclNumRows = getLocalNumRows();

    const bool verbose = Details::Behavior::verbose("CrsMatrix");
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsMatrix", "fillLocalMatrix");
      std::ostringstream os;
      os << *prefix << "lclNumRows: " << lclNumRows << endl;
      std::cerr << os.str ();
    }

    // The goals of this routine are first, to allocate and fill
    // packed 1-D storage (see below for an explanation) in the vals
    // array, and second, to give vals to the local matrix and
    // finalize the local matrix.  We only need k_ptrs, the packed 1-D
    // row offsets, within the scope of this routine, since we're only
    // filling the local matrix here (use fillLocalGraphAndMatrix() to
    // fill both the graph and the matrix at the same time).

    // get data from staticGraph_
    size_t nodeNumEntries   = staticGraph_->getLocalNumEntries ();
    size_t nodeNumAllocated = staticGraph_->getLocalAllocationSize ();
    row_map_type k_rowPtrs = staticGraph_->rowPtrsPacked_dev_; 

    row_map_type k_ptrs; // "packed" row offsets array
    values_type k_vals; // "packed" values array

    // May we ditch the old allocations for the packed (and otherwise
    // "optimized") allocations, later in this routine?  Request
    // optimized storage by default.
    bool requestOptimizedStorage = true;
    const bool default_OptimizeStorage =
      ! isStaticGraph() || staticGraph_->isStorageOptimized();
    if (! params.is_null() &&
        ! params->get("Optimize Storage", default_OptimizeStorage)) {
      requestOptimizedStorage = false;
    }
    // If we're not allowed to change a static graph, then we can't
    // change the storage of the matrix, either.  This means that if
    // the graph's storage isn't already optimized, we can't optimize
    // the matrix's storage either.  Check and give warning, as
    // appropriate.
    if (! staticGraph_->isStorageOptimized () &&
        requestOptimizedStorage) {
      TPETRA_ABUSE_WARNING
        (true, std::runtime_error, "You requested optimized storage "
         "by setting the \"Optimize Storage\" flag to \"true\" in "
         "the ParameterList, or by virtue of default behavior.  "
         "However, the associated CrsGraph was filled separately and "
         "requested not to optimize storage. Therefore, the "
         "CrsMatrix cannot optimize storage.");
      requestOptimizedStorage = false;
    }

    // NOTE: This does not work correctly w/ GCC 12.3 + CUDA due to a compiler bug.
    // See: https://github.com/trilinos/Trilinos/issues/12237
    //using row_entries_type = decltype (staticGraph_->k_numRowEntries_);
    using row_entries_type = typename crs_graph_type::num_row_entries_type;

    // The matrix's values are currently
    // stored in a 1-D format.  However, this format is "unpacked";
    // it doesn't necessarily have the same row offsets as indicated
    // by the ptrs array returned by allocRowPtrs.  This could
    // happen, for example, if the user 
    // fixed the number of matrix entries in
    // each row, but didn't fill all those entries.
    //
    // As above, we don't need to keep the "packed" row offsets
    // array ptrs here, but we do need it here temporarily, so we
    // have to allocate it.  We'll free ptrs later in this method.
    //
    // Note that this routine checks whether storage has already
    // been packed.  This is a common case for solution of nonlinear
    // PDEs using the finite element method, as long as the
    // structure of the sparse matrix does not change between linear
    // solves.
    if (nodeNumEntries != nodeNumAllocated) {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Unpacked 1-D storage: numEnt="
           << nodeNumEntries << ", allocSize=" << nodeNumAllocated
           << endl;
        std::cerr << os.str();
      }
      // We have to pack the 1-D storage, since the user didn't fill
      // up all requested storage.
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Allocate packed row offsets: "
           << (lclNumRows+1) << endl;
        std::cerr << os.str();
      }
      non_const_row_map_type tmpk_ptrs ("Tpetra::CrsGraph::ptr",
                                        lclNumRows+1);
      // Total number of entries in the matrix on the calling
      // process.  We will compute this in the loop below.  It's
      // cheap to compute and useful as a sanity check.
      size_t lclTotalNumEntries = 0;
      k_ptrs = tmpk_ptrs;
      {
        typename row_entries_type::const_type numRowEnt_h =
          staticGraph_->k_numRowEntries_;
        // This function can handle the counts being a host View.
        lclTotalNumEntries =
          Details::computeOffsetsFromCounts (tmpk_ptrs, numRowEnt_h);
      }

      // Allocate the "packed" values array.
      // It has exactly the right number of entries.
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Allocate packed values: "
           << lclTotalNumEntries << endl;
        std::cerr << os.str ();
      }
      k_vals = values_type ("Tpetra::CrsMatrix::val", lclTotalNumEntries);

      // Pack values_wdv into k_vals.  We will replace values_wdv below.
      pack_functor<
        typename values_type::non_const_type,
        typename values_type::const_type, 
        typename row_map_type::non_const_type, 
        typename row_map_type::const_type> valsPacker
        (k_vals, valuesUnpacked_wdv.getDeviceView(Access::ReadOnly),
         tmpk_ptrs, k_rowPtrs);

      using exec_space = typename decltype (k_vals)::execution_space;
      using range_type = Kokkos::RangePolicy<exec_space, LocalOrdinal>;
      Kokkos::parallel_for ("Tpetra::CrsMatrix pack values",
                            range_type (0, lclNumRows), valsPacker);
      valuesPacked_wdv = values_wdv_type(k_vals);
    }
    else { // We don't have to pack, so just set the pointer.
      valuesPacked_wdv = valuesUnpacked_wdv;
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Storage already packed: "
           << "valuesUnpacked_wdv: " << valuesUnpacked_wdv.extent(0) << endl;
        std::cerr << os.str();
      }
    }

    // May we ditch the old allocations for the packed one?
    if (requestOptimizedStorage) {
      // The user requested optimized storage, so we can dump the
      // unpacked 1-D storage, and keep the packed storage.
      valuesUnpacked_wdv = valuesPacked_wdv;
//      k_values1D_ = valuesPacked_wdv.getDeviceView(Access::ReadWrite);
      this->storageStatus_ = Details::STORAGE_1D_PACKED;
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  insertIndicesAndValues (crs_graph_type& graph,
                          RowInfo& rowInfo,
                          const typename crs_graph_type::SLocalGlobalViews& newInds,
                          const Teuchos::ArrayView<impl_scalar_type>& oldRowVals,
                          const Teuchos::ArrayView<const impl_scalar_type>& newRowVals,
                          const ELocalGlobal lg,
                          const ELocalGlobal I)
  {
    const size_t oldNumEnt = rowInfo.numEntries;
    const size_t numInserted = graph.insertIndices (rowInfo, newInds, lg, I);

    // Use of memcpy here works around an issue with GCC >= 4.9.0,
    // that probably relates to scalar_type vs. impl_scalar_type
    // aliasing.  See history of Tpetra_CrsGraph_def.hpp for
    // details; look for GCC_WORKAROUND macro definition.
    if (numInserted > 0) {
      const size_t startOffset = oldNumEnt;
      memcpy ((void*) &oldRowVals[startOffset], &newRowVals[0],
              numInserted * sizeof (impl_scalar_type));
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  insertLocalValues (const LocalOrdinal lclRow,
                     const Teuchos::ArrayView<const LocalOrdinal>& indices,
                     const Teuchos::ArrayView<const Scalar>& values,
                     const CombineMode CM)
  {
    using std::endl;
    const char tfecfFuncName[] = "insertLocalValues: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->isFillActive (), std::runtime_error,
       "Fill is not active.  After calling fillComplete, you must call "
       "resumeFill before you may insert entries into the matrix again.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (this->isStaticGraph (), std::runtime_error,
       "Cannot insert indices with static graph; use replaceLocalValues() "
       "instead.");
    // At this point, we know that myGraph_ is nonnull.
    crs_graph_type& graph = * (this->myGraph_);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (graph.colMap_.is_null (), std::runtime_error,
       "Cannot insert local indices without a column map.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (graph.isGloballyIndexed (),
       std::runtime_error, "Graph indices are global; use "
       "insertGlobalValues().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (values.size () != indices.size (), std::runtime_error,
       "values.size() = " << values.size ()
       << " != indices.size() = " << indices.size () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! graph.rowMap_->isNodeLocalElement (lclRow), std::runtime_error,
      "Local row index " << lclRow << " does not belong to this process.");

    if (! graph.indicesAreAllocated ()) {
      // We only allocate values at most once per process, so it's OK
      // to check TPETRA_VERBOSE here.
      const bool verbose = Details::Behavior::verbose("CrsMatrix");
      this->allocateValues (LocalIndices, GraphNotYetAllocated, verbose);
    }

#ifdef HAVE_TPETRA_DEBUG
    const size_t numEntriesToAdd = static_cast<size_t> (indices.size ());
    // In a debug build, test whether any of the given column indices
    // are not in the column Map.  Keep track of the invalid column
    // indices so we can tell the user about them.
    {
      using Teuchos::toString;

      const map_type& colMap = * (graph.colMap_);
      Teuchos::Array<LocalOrdinal> badColInds;
      bool allInColMap = true;
      for (size_t k = 0; k < numEntriesToAdd; ++k) {
        if (! colMap.isNodeLocalElement (indices[k])) {
          allInColMap = false;
          badColInds.push_back (indices[k]);
        }
      }
      if (! allInColMap) {
        std::ostringstream os;
        os << "You attempted to insert entries in owned row " << lclRow
           << ", at the following column indices: " << toString (indices)
           << "." << endl;
        os << "Of those, the following indices are not in the column Map on "
          "this process: " << toString (badColInds) << "." << endl << "Since "
          "the matrix has a column Map already, it is invalid to insert "
          "entries at those locations.";
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (true, std::invalid_argument, os.str ());
      }
    }
#endif // HAVE_TPETRA_DEBUG

    RowInfo rowInfo = graph.getRowInfo (lclRow);

    auto valsView = this->getValuesViewHostNonConst(rowInfo);
    if (CM == ADD) {
      auto fun = [&](size_t const k, size_t const /*start*/, size_t const offset) {
        valsView[offset] += values[k]; };
      std::function<void(size_t const, size_t const, size_t const)> cb(std::ref(fun));
      graph.insertLocalIndicesImpl(lclRow, indices, cb);
    } else if (CM == INSERT) {
      auto fun = [&](size_t const k, size_t const /*start*/, size_t const offset) {
        valsView[offset] = values[k]; };
      std::function<void(size_t const, size_t const, size_t const)> cb(std::ref(fun));
      graph.insertLocalIndicesImpl(lclRow, indices, cb);
    } else {
      std::ostringstream os;
      os << "You attempted to use insertLocalValues with CombineMode " << combineModeToString(CM)
         << "but this has not been implemented." << endl;
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::invalid_argument, os.str ());
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  insertLocalValues (const LocalOrdinal localRow,
                     const LocalOrdinal numEnt,
                     const Scalar vals[],
                     const LocalOrdinal cols[],
                     const CombineMode CM)
  {
    Teuchos::ArrayView<const LocalOrdinal> colsT (cols, numEnt);
    Teuchos::ArrayView<const Scalar> valsT (vals, numEnt);
    this->insertLocalValues (localRow, colsT, valsT, CM);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  insertGlobalValuesImpl (crs_graph_type& graph,
                          RowInfo& rowInfo,
                          const GlobalOrdinal gblColInds[],
                          const impl_scalar_type vals[],
                          const size_t numInputEnt)
  {
#ifdef HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "insertGlobalValuesImpl: ";
    const size_t origNumEnt = graph.getNumEntriesInLocalRow (rowInfo.localRow);
    const size_t curNumEnt = rowInfo.numEntries;
#endif // HAVE_TPETRA_DEBUG

    if (! graph.indicesAreAllocated ()) {
      // We only allocate values at most once per process, so it's OK
      // to check TPETRA_VERBOSE here.
      using ::Tpetra::Details::Behavior;
      const bool verbose = Behavior::verbose("CrsMatrix");
      this->allocateValues (GlobalIndices, GraphNotYetAllocated, verbose);
      // mfh 23 Jul 2017: allocateValues invalidates existing
      // getRowInfo results.  Once we get rid of lazy graph
      // allocation, we'll be able to move the getRowInfo call outside
      // of this method.
      rowInfo = graph.getRowInfo (rowInfo.localRow);
    }

    auto valsView = this->getValuesViewHostNonConst(rowInfo);
    auto fun = [&](size_t const k, size_t const /*start*/, size_t const offset){
                 valsView[offset] += vals[k];
                 };
    std::function<void(size_t const, size_t const, size_t const)> cb(std::ref(fun));
#ifdef HAVE_TPETRA_DEBUG
    //numInserted is only used inside the debug code below.
    auto numInserted =
#endif
    graph.insertGlobalIndicesImpl(rowInfo, gblColInds, numInputEnt, cb);

#ifdef HAVE_TPETRA_DEBUG
    size_t newNumEnt = curNumEnt + numInserted;
    const size_t chkNewNumEnt =
      graph.getNumEntriesInLocalRow (rowInfo.localRow);
    if (chkNewNumEnt != newNumEnt) {
      std::ostringstream os;
      os << std::endl << "newNumEnt = " << newNumEnt
         << " != graph.getNumEntriesInLocalRow(" << rowInfo.localRow
         << ") = " << chkNewNumEnt << "." << std::endl
         << "\torigNumEnt: " << origNumEnt << std::endl
         << "\tnumInputEnt: " << numInputEnt << std::endl
         << "\tgblColInds: [";
      for (size_t k = 0; k < numInputEnt; ++k) {
        os << gblColInds[k];
        if (k + size_t (1) < numInputEnt) {
          os << ",";
        }
      }
      os << "]" << std::endl
         << "\tvals: [";
      for (size_t k = 0; k < numInputEnt; ++k) {
        os << vals[k];
        if (k + size_t (1) < numInputEnt) {
          os << ",";
        }
      }
      os << "]" << std::endl;

      if (this->supportsRowViews ()) {
        values_host_view_type vals2;
        if (this->isGloballyIndexed ()) {
          global_inds_host_view_type gblColInds2;
          const GlobalOrdinal gblRow =
            graph.rowMap_->getGlobalElement (rowInfo.localRow);
          if (gblRow == 
              Tpetra::Details::OrdinalTraits<GlobalOrdinal>::invalid ()) {
            os << "Local row index " << rowInfo.localRow << " is invalid!" 
               << std::endl;
          }
          else {
            bool getViewThrew = false;
            try {
              this->getGlobalRowView (gblRow, gblColInds2, vals2);
            }
            catch (std::exception& e) {
              getViewThrew = true;
              os << "getGlobalRowView threw exception:" << std::endl
                 << e.what () << std::endl;
            }
            if (! getViewThrew) {
              os << "\tNew global column indices: ";
              for (size_t jjj = 0; jjj < gblColInds2.extent(0); jjj++)
                 os << gblColInds2[jjj] << " ";
              os << std::endl;
              os << "\tNew values: "; 
              for (size_t jjj = 0; jjj < vals2.extent(0); jjj++)
                 os << vals2[jjj] << " ";
              os << std::endl;
            }
          }
        }
        else if (this->isLocallyIndexed ()) {
          local_inds_host_view_type lclColInds2;
          this->getLocalRowView (rowInfo.localRow, lclColInds2, vals2);
          os << "\tNew local column indices: ";
          for (size_t jjj = 0; jjj < lclColInds2.extent(0); jjj++)
             os << lclColInds2[jjj] << " ";
          os << std::endl;
          os << "\tNew values: "; 
          for (size_t jjj = 0; jjj < vals2.extent(0); jjj++)
             os << vals2[jjj] << " ";
          os << std::endl;
        }
      }

      os << "Please report this bug to the Tpetra developers.";
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::logic_error, os.str ());
    }
#endif // HAVE_TPETRA_DEBUG
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  insertGlobalValues (const GlobalOrdinal gblRow,
                      const Teuchos::ArrayView<const GlobalOrdinal>& indices,
                      const Teuchos::ArrayView<const Scalar>& values)
  {
    using Teuchos::toString;
    using std::endl;
    typedef impl_scalar_type IST;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef Tpetra::Details::OrdinalTraits<LO> OTLO;
    typedef typename Teuchos::ArrayView<const GO>::size_type size_type;
    const char tfecfFuncName[] = "insertGlobalValues: ";

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (values.size () != indices.size (), std::runtime_error,
      "values.size() = " << values.size () << " != indices.size() = "
      << indices.size () << ".");
#endif // HAVE_TPETRA_DEBUG

    // getRowMap() is not thread safe, because it increments RCP's
    // reference count.  getCrsGraphRef() is thread safe.
    const map_type& rowMap = * (this->getCrsGraphRef ().rowMap_);
    const LO lclRow = rowMap.getLocalElement (gblRow);

    if (lclRow == OTLO::invalid ()) {
      // Input row is _not_ owned by the calling process.
      //
      // See a note (now deleted) from mfh 14 Dec 2012: If input row
      // is not in the row Map, it doesn't matter whether or not the
      // graph is static; the data just get stashed for later use by
      // globalAssemble().
      this->insertNonownedGlobalValues (gblRow, indices, values);
    }
    else { // Input row _is_ owned by the calling process
      if (this->isStaticGraph ()) {
        // Uh oh!  Not allowed to insert into owned rows in that case.
        const int myRank = rowMap.getComm ()->getRank ();
        const int numProcs = rowMap.getComm ()->getSize ();
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (true, std::runtime_error,
           "The matrix was constructed with a constant (\"static\") graph, "
           "yet the given global row index " << gblRow << " is in the row "
           "Map on the calling process (with rank " << myRank << ", of " <<
           numProcs << " process(es)).  In this case, you may not insert "
           "new entries into rows owned by the calling process.");
      }

      crs_graph_type& graph = * (this->myGraph_);
      const IST* const inputVals =
        reinterpret_cast<const IST*> (values.getRawPtr ());
      const GO* const inputGblColInds = indices.getRawPtr ();
      const size_t numInputEnt = indices.size ();
      RowInfo rowInfo = graph.getRowInfo (lclRow);

      // If the matrix has a column Map, check at this point whether
      // the column indices belong to the column Map.
      //
      // FIXME (mfh 16 May 2013) We may want to consider deferring the
      // test to the CrsGraph method, since it may have to do this
      // anyway.
      if (! graph.colMap_.is_null ()) {
        const map_type& colMap = * (graph.colMap_);
        // In a debug build, keep track of the nonowned ("bad") column
        // indices, so that we can display them in the exception
        // message.  In a release build, just ditch the loop early if
        // we encounter a nonowned column index.
#ifdef HAVE_TPETRA_DEBUG
        Teuchos::Array<GO> badColInds;
#endif // HAVE_TPETRA_DEBUG
        const size_type numEntriesToInsert = indices.size ();
        bool allInColMap = true;
        for (size_type k = 0; k < numEntriesToInsert; ++k) {
          if (! colMap.isNodeGlobalElement (indices[k])) {
            allInColMap = false;
#ifdef HAVE_TPETRA_DEBUG
            badColInds.push_back (indices[k]);
#else
            break;
#endif // HAVE_TPETRA_DEBUG
          }
        }
        if (! allInColMap) {
          std::ostringstream os;
          os << "You attempted to insert entries in owned row " << gblRow
             << ", at the following column indices: " << toString (indices)
             << "." << endl;
#ifdef HAVE_TPETRA_DEBUG
          os << "Of those, the following indices are not in the column Map "
            "on this process: " << toString (badColInds) << "." << endl
             << "Since the matrix has a column Map already, it is invalid "
            "to insert entries at those locations.";
#else
          os << "At least one of those indices is not in the column Map "
            "on this process." << endl << "It is invalid to insert into "
            "columns not in the column Map on the process that owns the "
            "row.";
#endif // HAVE_TPETRA_DEBUG
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (true, std::invalid_argument, os.str ());
        }
      }

      this->insertGlobalValuesImpl (graph, rowInfo, inputGblColInds,
                                    inputVals, numInputEnt);
    }
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  insertGlobalValues (const GlobalOrdinal globalRow,
                      const LocalOrdinal numEnt,
                      const Scalar vals[],
                      const GlobalOrdinal inds[])
  {
    Teuchos::ArrayView<const GlobalOrdinal> indsT (inds, numEnt);
    Teuchos::ArrayView<const Scalar> valsT (vals, numEnt);
    this->insertGlobalValues (globalRow, indsT, valsT);
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  insertGlobalValuesFiltered(
    const GlobalOrdinal gblRow,
    const Teuchos::ArrayView<const GlobalOrdinal>& indices,
    const Teuchos::ArrayView<const Scalar>& values,
    const bool debug)
  {
    typedef impl_scalar_type IST;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef Tpetra::Details::OrdinalTraits<LO> OTLO;
    const char tfecfFuncName[] = "insertGlobalValuesFiltered: ";

    if (debug) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (values.size () != indices.size (), std::runtime_error,
         "values.size() = " << values.size () << " != indices.size() = "
         << indices.size () << ".");
    }

    // getRowMap() is not thread safe, because it increments RCP's
    // reference count.  getCrsGraphRef() is thread safe.
    const map_type& rowMap = * (this->getCrsGraphRef ().rowMap_);
    const LO lclRow = rowMap.getLocalElement (gblRow);
    if (lclRow == OTLO::invalid ()) {
      // Input row is _not_ owned by the calling process.
      //
      // See a note (now deleted) from mfh 14 Dec 2012: If input row
      // is not in the row Map, it doesn't matter whether or not the
      // graph is static; the data just get stashed for later use by
      // globalAssemble().
      this->insertNonownedGlobalValues (gblRow, indices, values);
    }
    else { // Input row _is_ owned by the calling process
      if (this->isStaticGraph ()) {
        // Uh oh!  Not allowed to insert into owned rows in that case.
        const int myRank = rowMap.getComm ()->getRank ();
        const int numProcs = rowMap.getComm ()->getSize ();
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (true, std::runtime_error,
           "The matrix was constructed with a constant (\"static\") graph, "
           "yet the given global row index " << gblRow << " is in the row "
           "Map on the calling process (with rank " << myRank << ", of " <<
           numProcs << " process(es)).  In this case, you may not insert "
           "new entries into rows owned by the calling process.");
      }

      crs_graph_type& graph = * (this->myGraph_);
      const IST* const inputVals =
        reinterpret_cast<const IST*> (values.getRawPtr ());
      const GO* const inputGblColInds = indices.getRawPtr ();
      const size_t numInputEnt = indices.size ();
      RowInfo rowInfo = graph.getRowInfo (lclRow);

      if (!graph.colMap_.is_null() && graph.isLocallyIndexed()) {
        // This branch is similar in function to the following branch, but for
        // the special case that the target graph is locally indexed.
        // In this case, we cannot simply filter
        // out global indices that don't exist on the receiving process and
        // insert the remaining (global) indices, but we must convert them (the
        // remaining global indices) to local and call `insertLocalValues`.
        const map_type& colMap = * (graph.colMap_);
        size_t curOffset = 0;
        while (curOffset < numInputEnt) {
          // Find a sequence of input indices that are in the column Map on the
          // calling process. Doing a sequence at a time, instead of one at a
          // time, amortizes some overhead.
          Teuchos::Array<LO> lclIndices;
          size_t endOffset = curOffset;
          for ( ; endOffset < numInputEnt; ++endOffset) {
            auto lclIndex = colMap.getLocalElement(inputGblColInds[endOffset]);
            if (lclIndex != OTLO::invalid())
              lclIndices.push_back(lclIndex);
            else
              break;
          }
          // curOffset, endOffset: half-exclusive range of indices in the column
          // Map on the calling process. If endOffset == curOffset, the range is
          // empty.
          const LO numIndInSeq = (endOffset - curOffset);
          if (numIndInSeq != 0) {
            this->insertLocalValues(lclRow, lclIndices(), values(curOffset, numIndInSeq));
          }
          // Invariant before the increment line: Either endOffset ==
          // numInputEnt, or inputGblColInds[endOffset] is not in the column Map
          // on the calling process.
          if (debug) {
            const bool invariant = endOffset == numInputEnt ||
              colMap.getLocalElement (inputGblColInds[endOffset]) == OTLO::invalid ();
            TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
              (! invariant, std::logic_error, std::endl << "Invariant failed!");
          }
          curOffset = endOffset + 1;
        }
      }
      else if (! graph.colMap_.is_null ()) { // We have a column Map.
        const map_type& colMap = * (graph.colMap_);
        size_t curOffset = 0;
        while (curOffset < numInputEnt) {
          // Find a sequence of input indices that are in the column
          // Map on the calling process.  Doing a sequence at a time,
          // instead of one at a time, amortizes some overhead.
          size_t endOffset = curOffset;
          for ( ; endOffset < numInputEnt &&
                  colMap.getLocalElement (inputGblColInds[endOffset]) != OTLO::invalid ();
                ++endOffset)
            {}
          // curOffset, endOffset: half-exclusive range of indices in
          // the column Map on the calling process.  If endOffset ==
          // curOffset, the range is empty.
          const LO numIndInSeq = (endOffset - curOffset);
          if (numIndInSeq != 0) {
            rowInfo = graph.getRowInfo(lclRow);  // KDD 5/19 Need fresh RowInfo in each loop iteration
            this->insertGlobalValuesImpl (graph, rowInfo,
                                          inputGblColInds + curOffset,
                                          inputVals + curOffset,
                                          numIndInSeq);
          }
          // Invariant before the increment line: Either endOffset ==
          // numInputEnt, or inputGblColInds[endOffset] is not in the
          // column Map on the calling process.
          if (debug) {
            const bool invariant = endOffset == numInputEnt ||
              colMap.getLocalElement (inputGblColInds[endOffset]) == OTLO::invalid ();
            TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
              (! invariant, std::logic_error, std::endl << "Invariant failed!");
          }
          curOffset = endOffset + 1;
        }
      }
      else { // we don't have a column Map.
        this->insertGlobalValuesImpl (graph, rowInfo, inputGblColInds,
                                      inputVals, numInputEnt);
      }
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  insertGlobalValuesFilteredChecked(
    const GlobalOrdinal gblRow,
    const Teuchos::ArrayView<const GlobalOrdinal>& indices,
    const Teuchos::ArrayView<const Scalar>& values,
    const char* const prefix,
    const bool debug,
    const bool verbose)
  {
    using Details::verbosePrintArray;
    using std::endl;

    try {
      insertGlobalValuesFiltered(gblRow, indices, values, debug);
    }
    catch(std::exception& e) {
      std::ostringstream os;
      if (verbose) {
        const size_t maxNumToPrint =
          Details::Behavior::verbosePrintCountThreshold();
        os << *prefix << ": insertGlobalValuesFiltered threw an "
          "exception: " << e.what() << endl
           << "Global row index: " << gblRow << endl;
        verbosePrintArray(os, indices, "Global column indices",
                          maxNumToPrint);
        os << endl;
        verbosePrintArray(os, values, "Values", maxNumToPrint);
        os << endl;
      }
      else {
        os << ": insertGlobalValuesFiltered threw an exception: "
           << e.what();
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, os.str());
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceLocalValuesImpl (impl_scalar_type rowVals[],
                          const crs_graph_type& graph,
                          const RowInfo& rowInfo,
                          const LocalOrdinal inds[],
                          const impl_scalar_type newVals[],
                          const LocalOrdinal numElts) 
  {
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const bool sorted = graph.isSorted ();

    size_t hint = 0; // Guess for the current index k into rowVals
    LO numValid = 0; // number of valid local column indices

    if (graph.isLocallyIndexed ()) {
      // Get a view of the column indices in the row.  This amortizes
      // the cost of getting the view over all the entries of inds.
      auto colInds = graph.getLocalIndsViewHost (rowInfo);

      for (LO j = 0; j < numElts; ++j) {
        const LO lclColInd = inds[j];
        const size_t offset =
          KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                       lclColInd, hint, sorted);
        if (offset != rowInfo.numEntries) {
          rowVals[offset] = newVals[j];
          hint = offset + 1;
          ++numValid;
        }
      }
    }
    else if (graph.isGloballyIndexed ()) {
      if (graph.colMap_.is_null ()) {
        return Teuchos::OrdinalTraits<LO>::invalid ();
      }
      const map_type colMap = * (graph.colMap_);

      // Get a view of the column indices in the row.  This amortizes
      // the cost of getting the view over all the entries of inds.
      auto colInds = graph.getGlobalIndsViewHost (rowInfo);

      for (LO j = 0; j < numElts; ++j) {
        const GO gblColInd = colMap.getGlobalElement (inds[j]);
        if (gblColInd != Teuchos::OrdinalTraits<GO>::invalid ()) {
          const size_t offset =
            KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                         gblColInd, hint, sorted);
          if (offset != rowInfo.numEntries) {
            rowVals[offset] = newVals[j];
            hint = offset + 1;
            ++numValid;
          }
        }
      }
    }
    // NOTE (mfh 26 Jun 2014, 26 Nov 2015) In the current version of
    // CrsGraph and CrsMatrix, it's possible for a matrix (or graph)
    // to be neither locally nor globally indexed on a process.
    // This means that the graph or matrix has no entries on that
    // process.  Epetra also works like this.  It's related to lazy
    // allocation (on first insertion, not at graph / matrix
    // construction).  Lazy allocation will go away because it is
    // not thread scalable.

    return numValid;
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceLocalValues (const LocalOrdinal localRow,
                      const Teuchos::ArrayView<const LocalOrdinal>& lclCols,
                      const Teuchos::ArrayView<const Scalar>& vals)
  {
    typedef LocalOrdinal LO;

    const LO numInputEnt = static_cast<LO> (lclCols.size ());
    if (static_cast<LO> (vals.size ()) != numInputEnt) {
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    const LO* const inputInds = lclCols.getRawPtr ();
    const Scalar* const inputVals = vals.getRawPtr ();
    return this->replaceLocalValues (localRow, numInputEnt,
                                     inputVals, inputInds);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    local_ordinal_type
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceLocalValues(
    const local_ordinal_type localRow,
    const Kokkos::View<const local_ordinal_type*, Kokkos::AnonymousSpace>& inputInds,
    const Kokkos::View<const impl_scalar_type*, Kokkos::AnonymousSpace>& inputVals)
  {
    using LO = local_ordinal_type;
    const LO numInputEnt = inputInds.extent(0);
    if (numInputEnt != static_cast<LO>(inputVals.extent(0))) {
      return Teuchos::OrdinalTraits<LO>::invalid();
    }
    const Scalar* const inVals =
      reinterpret_cast<const Scalar*>(inputVals.data());
    return this->replaceLocalValues(localRow, numInputEnt,
                                    inVals, inputInds.data());
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceLocalValues (const LocalOrdinal localRow,
                      const LocalOrdinal numEnt,
                      const Scalar inputVals[],
                      const LocalOrdinal inputCols[])
  {
    typedef impl_scalar_type IST;
    typedef LocalOrdinal LO;

    if (! this->isFillActive () || this->staticGraph_.is_null ()) {
      // Fill must be active and the "nonconst" graph must exist.
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    const crs_graph_type& graph = * (this->staticGraph_);
    const RowInfo rowInfo = graph.getRowInfo (localRow);

    if (rowInfo.localRow == Teuchos::OrdinalTraits<size_t>::invalid ()) {
      // The calling process does not own this row, so it is not
      // allowed to modify its values.
      return static_cast<LO> (0);
    }
    auto curRowVals = this->getValuesViewHostNonConst (rowInfo);
    const IST* const inVals = reinterpret_cast<const IST*> (inputVals);
    return this->replaceLocalValuesImpl (curRowVals.data (), graph, rowInfo,
                                         inputCols, inVals, numEnt);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceGlobalValuesImpl (impl_scalar_type rowVals[],
                           const crs_graph_type& graph,
                           const RowInfo& rowInfo,
                           const GlobalOrdinal inds[],
                           const impl_scalar_type newVals[],
                           const LocalOrdinal numElts)
  {
    Teuchos::ArrayView<const GlobalOrdinal> indsT(inds, numElts);
    auto fun =
      [&](size_t const k, size_t const /*start*/, size_t const offset) {
        rowVals[offset] = newVals[k];
      };
    std::function<void(size_t const, size_t const, size_t const)> cb(std::ref(fun));
    return graph.findGlobalIndices(rowInfo, indsT, cb);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceGlobalValues (const GlobalOrdinal globalRow,
                       const Teuchos::ArrayView<const GlobalOrdinal>& inputGblColInds,
                       const Teuchos::ArrayView<const Scalar>& inputVals)
  {
    typedef LocalOrdinal LO;

    const LO numInputEnt = static_cast<LO> (inputGblColInds.size ());
    if (static_cast<LO> (inputVals.size ()) != numInputEnt) {
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    return this->replaceGlobalValues (globalRow, numInputEnt,
                                      inputVals.getRawPtr (),
                                      inputGblColInds.getRawPtr ());
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceGlobalValues (const GlobalOrdinal globalRow,
                       const LocalOrdinal numEnt,
                       const Scalar inputVals[],
                       const GlobalOrdinal inputGblColInds[])
  {
    typedef impl_scalar_type IST;
    typedef LocalOrdinal LO;

    if (! this->isFillActive () || this->staticGraph_.is_null ()) {
      // Fill must be active and the "nonconst" graph must exist.
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    const crs_graph_type& graph = * (this->staticGraph_);

    const RowInfo rowInfo = graph.getRowInfoFromGlobalRowIndex (globalRow);
    if (rowInfo.localRow == Teuchos::OrdinalTraits<size_t>::invalid ()) {
      // The input local row is invalid on the calling process,
      // which means that the calling process summed 0 entries.
      return static_cast<LO> (0);
    }

    auto curRowVals = this->getValuesViewHostNonConst (rowInfo);
    const IST* const inVals = reinterpret_cast<const IST*> (inputVals);
    return this->replaceGlobalValuesImpl (curRowVals.data (), graph, rowInfo,
                                          inputGblColInds, inVals, numEnt);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    local_ordinal_type
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceGlobalValues(
    const global_ordinal_type globalRow,
    const Kokkos::View<const global_ordinal_type*, Kokkos::AnonymousSpace>& inputInds,
    const Kokkos::View<const impl_scalar_type*, Kokkos::AnonymousSpace>& inputVals)
  {
    // We use static_assert here to check the template parameters,
    // rather than std::enable_if (e.g., on the return value, to
    // enable compilation only if the template parameters match the
    // desired attributes).  This turns obscure link errors into
    // clear compilation errors.  It also makes the return value a
    // lot easier to see.
    using LO = local_ordinal_type;
    const LO numInputEnt = static_cast<LO>(inputInds.extent(0));
    if (static_cast<LO>(inputVals.extent(0)) != numInputEnt) {
      return Teuchos::OrdinalTraits<LO>::invalid();
    }
    const Scalar* const inVals =
      reinterpret_cast<const Scalar*>(inputVals.data());
    return this->replaceGlobalValues(globalRow, numInputEnt, inVals,
                                     inputInds.data());
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  sumIntoGlobalValuesImpl (impl_scalar_type rowVals[],
                           const crs_graph_type& graph,
                           const RowInfo& rowInfo,
                           const GlobalOrdinal inds[],
                           const impl_scalar_type newVals[],
                           const LocalOrdinal numElts,
                           const bool atomic)
  {
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;

    const bool sorted = graph.isSorted ();

    size_t hint = 0; // guess at the index's relative offset in the row
    LO numValid = 0; // number of valid input column indices

    if (graph.isLocallyIndexed ()) {
      // NOTE (mfh 04 Nov 2015) Dereferencing an RCP or reading its
      // pointer does NOT change its reference count.  Thus, this
      // code is still thread safe.
      if (graph.colMap_.is_null ()) {
        // NO input column indices are valid in this case, since if
        // the column Map is null on the calling process, then the
        // calling process owns no graph entries.
        return numValid;
      }
      const map_type& colMap = * (graph.colMap_);

      // Get a view of the column indices in the row.  This amortizes
      // the cost of getting the view over all the entries of inds.
      auto colInds = graph.getLocalIndsViewHost (rowInfo);
      const LO LINV = Teuchos::OrdinalTraits<LO>::invalid ();

      for (LO j = 0; j < numElts; ++j) {
        const LO lclColInd = colMap.getLocalElement (inds[j]);
        if (lclColInd != LINV) {
          const size_t offset =
            KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                         lclColInd, hint, sorted);
          if (offset != rowInfo.numEntries) {
            if (atomic) {
              Kokkos::atomic_add (&rowVals[offset], newVals[j]);
            }
            else {
              rowVals[offset] += newVals[j];
            }
            hint = offset + 1;
            numValid++;
          }
        }
      }
    }
    else if (graph.isGloballyIndexed ()) {
      // Get a view of the column indices in the row.  This amortizes
      // the cost of getting the view over all the entries of inds.
      auto colInds = graph.getGlobalIndsViewHost (rowInfo);

      for (LO j = 0; j < numElts; ++j) {
        const GO gblColInd = inds[j];
        const size_t offset =
          KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                       gblColInd, hint, sorted);
        if (offset != rowInfo.numEntries) {
          if (atomic) {
            Kokkos::atomic_add (&rowVals[offset], newVals[j]);
          }
          else {
            rowVals[offset] += newVals[j];
          }
          hint = offset + 1;
          numValid++;
        }
      }
    }
    // If the graph is neither locally nor globally indexed on the
    // calling process, that means the calling process has no graph
    // entries.  Thus, none of the input column indices are valid.

    return numValid;
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  sumIntoGlobalValues (const GlobalOrdinal gblRow,
                       const Teuchos::ArrayView<const GlobalOrdinal>& inputGblColInds,
                       const Teuchos::ArrayView<const Scalar>& inputVals,
                       const bool atomic)
  {
    typedef LocalOrdinal LO;

    const LO numInputEnt = static_cast<LO> (inputGblColInds.size ());
    if (static_cast<LO> (inputVals.size ()) != numInputEnt) {
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    return this->sumIntoGlobalValues (gblRow, numInputEnt,
                                      inputVals.getRawPtr (),
                                      inputGblColInds.getRawPtr (),
                                      atomic);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  sumIntoGlobalValues (const GlobalOrdinal gblRow,
                       const LocalOrdinal numInputEnt,
                       const Scalar inputVals[],
                       const GlobalOrdinal inputGblColInds[],
                       const bool atomic)
  {
    typedef impl_scalar_type IST;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;

    if (! this->isFillActive () || this->staticGraph_.is_null ()) {
      // Fill must be active and the "nonconst" graph must exist.
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    const crs_graph_type& graph = * (this->staticGraph_);

    const RowInfo rowInfo = graph.getRowInfoFromGlobalRowIndex (gblRow);
    if (rowInfo.localRow == Teuchos::OrdinalTraits<size_t>::invalid ()) {
      // mfh 23 Mar 2017, 26 Jul 2017: This branch may not be not
      // thread safe in a debug build, in part because it uses
      // Teuchos::ArrayView, and in part because of the data structure
      // used to stash outgoing entries.
      using Teuchos::ArrayView;
      ArrayView<const GO> inputGblColInds_av(
        numInputEnt == 0 ? nullptr : inputGblColInds,
        numInputEnt);
      ArrayView<const Scalar> inputVals_av(
        numInputEnt == 0 ? nullptr :
        inputVals, numInputEnt);
      // gblRow is not in the row Map on the calling process, so stash
      // the given entries away in a separate data structure.
      // globalAssemble() (called during fillComplete()) will exchange
      // that data and sum it in using sumIntoGlobalValues().
      this->insertNonownedGlobalValues (gblRow, inputGblColInds_av,
                                        inputVals_av);
      // FIXME (mfh 08 Jul 2014) It's not clear what to return here,
      // since we won't know whether the given indices were valid
      // until globalAssemble (called in fillComplete) is called.
      // That's why insertNonownedGlobalValues doesn't return
      // anything.  Just for consistency, I'll return the number of
      // entries that the user gave us.
      return numInputEnt;
    }
    else { // input row is in the row Map on the calling process
      auto curRowVals = this->getValuesViewHostNonConst (rowInfo);
      const IST* const inVals = reinterpret_cast<const IST*> (inputVals);
      return this->sumIntoGlobalValuesImpl (curRowVals.data (), graph, rowInfo,
                                            inputGblColInds, inVals,
                                            numInputEnt, atomic);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  transformLocalValues (const LocalOrdinal lclRow,
                        const LocalOrdinal numInputEnt,
                        const impl_scalar_type inputVals[],
                        const LocalOrdinal inputCols[],
                        std::function<impl_scalar_type (const impl_scalar_type&, const impl_scalar_type&) > f,
                        const bool atomic)
  {
    using Tpetra::Details::OrdinalTraits;
    typedef LocalOrdinal LO;

    if (! this->isFillActive () || this->staticGraph_.is_null ()) {
      // Fill must be active and the "nonconst" graph must exist.
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    const crs_graph_type& graph = * (this->staticGraph_);
    const RowInfo rowInfo = graph.getRowInfo (lclRow);

    if (rowInfo.localRow == OrdinalTraits<size_t>::invalid ()) {
      // The calling process does not own this row, so it is not
      // allowed to modify its values.
      return static_cast<LO> (0);
    }
    auto curRowVals = this->getValuesViewHostNonConst (rowInfo);
    return this->transformLocalValues (curRowVals.data (), graph,
                                       rowInfo, inputCols, inputVals,
                                       numInputEnt, f, atomic);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  transformGlobalValues (const GlobalOrdinal gblRow,
                         const LocalOrdinal numInputEnt,
                         const impl_scalar_type inputVals[],
                         const GlobalOrdinal inputCols[],
                         std::function<impl_scalar_type (const impl_scalar_type&, const impl_scalar_type&) > f,
                         const bool atomic)
  {
    using Tpetra::Details::OrdinalTraits;
    typedef LocalOrdinal LO;

    if (! this->isFillActive () || this->staticGraph_.is_null ()) {
      // Fill must be active and the "nonconst" graph must exist.
      return OrdinalTraits<LO>::invalid ();
    }
    const crs_graph_type& graph = * (this->staticGraph_);
    const RowInfo rowInfo = graph.getRowInfoFromGlobalRowIndex (gblRow);

    if (rowInfo.localRow == OrdinalTraits<size_t>::invalid ()) {
      // The calling process does not own this row, so it is not
      // allowed to modify its values.
      return static_cast<LO> (0);
    }
    auto curRowVals = this->getValuesViewHostNonConst (rowInfo);
    return this->transformGlobalValues (curRowVals.data (), graph,
                                        rowInfo, inputCols, inputVals,
                                        numInputEnt, f, atomic);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  transformLocalValues (impl_scalar_type rowVals[],
                        const crs_graph_type& graph,
                        const RowInfo& rowInfo,
                        const LocalOrdinal inds[],
                        const impl_scalar_type newVals[],
                        const LocalOrdinal numElts,
                        std::function<impl_scalar_type (const impl_scalar_type&, const impl_scalar_type&) > f,
                        const bool atomic)
  {
    typedef impl_scalar_type ST;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;

    //if (newVals.extent (0) != inds.extent (0)) {
    // The sizes of the input arrays must match.
    //return Tpetra::Details::OrdinalTraits<LO>::invalid ();
    //}
    //const LO numElts = static_cast<LO> (inds.extent (0));
    const bool sorted = graph.isSorted ();

    LO numValid = 0; // number of valid input column indices
    size_t hint = 0; // Guess for the current index k into rowVals

    if (graph.isLocallyIndexed ()) {
      // Get a view of the column indices in the row.  This amortizes
      // the cost of getting the view over all the entries of inds.
      auto colInds = graph.getLocalIndsViewHost (rowInfo);

      for (LO j = 0; j < numElts; ++j) {
        const LO lclColInd = inds[j];
        const size_t offset =
          KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                       lclColInd, hint, sorted);
        if (offset != rowInfo.numEntries) {
          if (atomic) {
            // NOTE (mfh 30 Nov 2015) The commented-out code is
            // wrong because another thread may have changed
            // rowVals[offset] between those two lines of code.
            //
            //const ST newVal = f (rowVals[offset], newVals[j]);
            //Kokkos::atomic_assign (&rowVals[offset], newVal);

            volatile ST* const dest = &rowVals[offset];
            (void) atomic_binary_function_update (dest, newVals[j], f);
          }
          else {
            // use binary function f
            rowVals[offset] = f (rowVals[offset], newVals[j]);
          }
          hint = offset + 1;
          ++numValid;
        }
      }
    }
    else if (graph.isGloballyIndexed ()) {
      // NOTE (mfh 26 Nov 2015) Dereferencing an RCP or reading its
      // pointer does NOT change its reference count.  Thus, this
      // code is still thread safe.
      if (graph.colMap_.is_null ()) {
        // NO input column indices are valid in this case.  Either
        // the column Map hasn't been set yet (so local indices
        // don't exist yet), or the calling process owns no graph
        // entries.
        return numValid;
      }
      const map_type& colMap = * (graph.colMap_);
      // Get a view of the column indices in the row.  This amortizes
      // the cost of getting the view over all the entries of inds.
      auto colInds = graph.getGlobalIndsViewHost (rowInfo);

      const GO GINV = Teuchos::OrdinalTraits<GO>::invalid ();
      for (LO j = 0; j < numElts; ++j) {
        const GO gblColInd = colMap.getGlobalElement (inds[j]);
        if (gblColInd != GINV) {
          const size_t offset =
            KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                         gblColInd, hint, sorted);
          if (offset != rowInfo.numEntries) {
            if (atomic) {
              // NOTE (mfh 30 Nov 2015) The commented-out code is
              // wrong because another thread may have changed
              // rowVals[offset] between those two lines of code.
              //
              //const ST newVal = f (rowVals[offset], newVals[j]);
              //Kokkos::atomic_assign (&rowVals[offset], newVal);

              volatile ST* const dest = &rowVals[offset];
              (void) atomic_binary_function_update (dest, newVals[j], f);
            }
            else {
              // use binary function f
              rowVals[offset] = f (rowVals[offset], newVals[j]);
            }
            hint = offset + 1;
            numValid++;
          }
        }
      }
    }
    // If the graph is neither locally nor globally indexed on the
    // calling process, that means the calling process has no graph
    // entries.  Thus, none of the input column indices are valid.

    return numValid;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  transformGlobalValues (impl_scalar_type rowVals[],
                         const crs_graph_type& graph,
                         const RowInfo& rowInfo,
                         const GlobalOrdinal inds[],
                         const impl_scalar_type newVals[],
                         const LocalOrdinal numElts,
                         std::function<impl_scalar_type (const impl_scalar_type&, const impl_scalar_type&) > f,
                         const bool atomic)
  {
    typedef impl_scalar_type ST;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;

    //if (newVals.extent (0) != inds.extent (0)) {
    // The sizes of the input arrays must match.
    //return Tpetra::Details::OrdinalTraits<LO>::invalid ();
    //}
    //const LO numElts = static_cast<LO> (inds.extent (0));
    const bool sorted = graph.isSorted ();

    LO numValid = 0; // number of valid input column indices
    size_t hint = 0; // Guess for the current index k into rowVals

    if (graph.isGloballyIndexed ()) {
      // Get a view of the column indices in the row.  This amortizes
      // the cost of getting the view over all the entries of inds.
      auto colInds = graph.getGlobalIndsViewHost (rowInfo);

      for (LO j = 0; j < numElts; ++j) {
        const GO gblColInd = inds[j];
        const size_t offset =
          KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                       gblColInd, hint, sorted);
        if (offset != rowInfo.numEntries) {
          if (atomic) {
            // NOTE (mfh 30 Nov 2015) The commented-out code is
            // wrong because another thread may have changed
            // rowVals[offset] between those two lines of code.
            //
            //const ST newVal = f (rowVals[offset], newVals[j]);
            //Kokkos::atomic_assign (&rowVals[offset], newVal);

            volatile ST* const dest = &rowVals[offset];
            (void) atomic_binary_function_update (dest, newVals[j], f);
          }
          else {
            // use binary function f
            rowVals[offset] = f (rowVals[offset], newVals[j]);
          }
          hint = offset + 1;
          ++numValid;
        }
      }
    }
    else if (graph.isLocallyIndexed ()) {
      // NOTE (mfh 26 Nov 2015) Dereferencing an RCP or reading its
      // pointer does NOT change its reference count.  Thus, this
      // code is still thread safe.
      if (graph.colMap_.is_null ()) {
        // NO input column indices are valid in this case.  Either the
        // column Map hasn't been set yet (so local indices don't
        // exist yet), or the calling process owns no graph entries.
        return numValid;
      }
      const map_type& colMap = * (graph.colMap_);
      // Get a view of the column indices in the row.  This amortizes
      // the cost of getting the view over all the entries of inds.
      auto colInds = graph.getLocalIndsViewHost (rowInfo);

      const LO LINV = Teuchos::OrdinalTraits<LO>::invalid ();
      for (LO j = 0; j < numElts; ++j) {
        const LO lclColInd = colMap.getLocalElement (inds[j]);
        if (lclColInd != LINV) {
          const size_t offset =
            KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                         lclColInd, hint, sorted);
          if (offset != rowInfo.numEntries) {
            if (atomic) {
              // NOTE (mfh 30 Nov 2015) The commented-out code is
              // wrong because another thread may have changed
              // rowVals[offset] between those two lines of code.
              //
              //const ST newVal = f (rowVals[offset], newVals[j]);
              //Kokkos::atomic_assign (&rowVals[offset], newVal);

              volatile ST* const dest = &rowVals[offset];
              (void) atomic_binary_function_update (dest, newVals[j], f);
            }
            else {
              // use binary function f
              rowVals[offset] = f (rowVals[offset], newVals[j]);
            }
            hint = offset + 1;
            numValid++;
          }
        }
      }
    }
    // If the graph is neither locally nor globally indexed on the
    // calling process, that means the calling process has no graph
    // entries.  Thus, none of the input column indices are valid.

    return numValid;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  sumIntoLocalValuesImpl (impl_scalar_type rowVals[],
                          const crs_graph_type& graph,
                          const RowInfo& rowInfo,
                          const LocalOrdinal inds[],
                          const impl_scalar_type newVals[],
                          const LocalOrdinal numElts,
                          const bool atomic)
  {
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;

    const bool sorted = graph.isSorted ();

    size_t hint = 0; // Guess for the current index k into rowVals
    LO numValid = 0; // number of valid local column indices

    if (graph.isLocallyIndexed ()) {
      // Get a view of the column indices in the row.  This amortizes
      // the cost of getting the view over all the entries of inds.
      auto colInds = graph.getLocalIndsViewHost (rowInfo);

      for (LO j = 0; j < numElts; ++j) {
        const LO lclColInd = inds[j];
        const size_t offset =
          KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                       lclColInd, hint, sorted);
        if (offset != rowInfo.numEntries) {
          if (atomic) {
            Kokkos::atomic_add (&rowVals[offset], newVals[j]);
          }
          else {
            rowVals[offset] += newVals[j];
          }
          hint = offset + 1;
          ++numValid;
        }
      }
    }
    else if (graph.isGloballyIndexed ()) {
      if (graph.colMap_.is_null ()) {
        return Teuchos::OrdinalTraits<LO>::invalid ();
      }
      const map_type colMap = * (graph.colMap_);

      // Get a view of the column indices in the row.  This amortizes
      // the cost of getting the view over all the entries of inds.
      auto colInds = graph.getGlobalIndsViewHost (rowInfo);

      for (LO j = 0; j < numElts; ++j) {
        const GO gblColInd = colMap.getGlobalElement (inds[j]);
        if (gblColInd != Teuchos::OrdinalTraits<GO>::invalid ()) {
          const size_t offset =
            KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                         gblColInd, hint, sorted);
          if (offset != rowInfo.numEntries) {
            if (atomic) {
              Kokkos::atomic_add (&rowVals[offset], newVals[j]);
            }
            else {
              rowVals[offset] += newVals[j];
            }
            hint = offset + 1;
            ++numValid;
          }
        }
      }
    }
    // NOTE (mfh 26 Jun 2014, 26 Nov 2015) In the current version of
    // CrsGraph and CrsMatrix, it's possible for a matrix (or graph)
    // to be neither locally nor globally indexed on a process.
    // This means that the graph or matrix has no entries on that
    // process.  Epetra also works like this.  It's related to lazy
    // allocation (on first insertion, not at graph / matrix
    // construction).  Lazy allocation will go away because it is
    // not thread scalable.

    return numValid;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  sumIntoLocalValues (const LocalOrdinal localRow,
                      const Teuchos::ArrayView<const LocalOrdinal>& indices,
                      const Teuchos::ArrayView<const Scalar>& values,
                      const bool atomic)
  {
    using LO = local_ordinal_type;
    const LO numInputEnt = static_cast<LO>(indices.size());
    if (static_cast<LO>(values.size()) != numInputEnt) {
      return Teuchos::OrdinalTraits<LO>::invalid();
    }
    const LO* const inputInds = indices.getRawPtr();
    const scalar_type* const inputVals = values.getRawPtr();
    return this->sumIntoLocalValues(localRow, numInputEnt,
                                    inputVals, inputInds, atomic);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    local_ordinal_type
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  sumIntoLocalValues(
    const local_ordinal_type localRow,
    const Kokkos::View<const local_ordinal_type*, Kokkos::AnonymousSpace>& inputInds,
    const Kokkos::View<const impl_scalar_type*, Kokkos::AnonymousSpace>& inputVals,
    const bool atomic)
  {
    using LO = local_ordinal_type;
    const LO numInputEnt = static_cast<LO>(inputInds.extent(0));
    if (static_cast<LO>(inputVals.extent(0)) != numInputEnt) {
      return Teuchos::OrdinalTraits<LO>::invalid();
    }
    const scalar_type* inVals =
      reinterpret_cast<const scalar_type*>(inputVals.data());
    return this->sumIntoLocalValues(localRow, numInputEnt, inVals,
                                    inputInds.data(), atomic);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  sumIntoLocalValues (const LocalOrdinal localRow,
                      const LocalOrdinal numEnt,
                      const Scalar vals[],
                      const LocalOrdinal cols[],
                      const bool atomic)
  {
    typedef impl_scalar_type IST;
    typedef LocalOrdinal LO;

    if (! this->isFillActive () || this->staticGraph_.is_null ()) {
      // Fill must be active and the "nonconst" graph must exist.
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    const crs_graph_type& graph = * (this->staticGraph_);
    const RowInfo rowInfo = graph.getRowInfo (localRow);

    if (rowInfo.localRow == Teuchos::OrdinalTraits<size_t>::invalid ()) {
      // The calling process does not own this row, so it is not
      // allowed to modify its values.
      return static_cast<LO> (0);
    }
    auto curRowVals = this->getValuesViewHostNonConst (rowInfo);
    const IST* const inputVals = reinterpret_cast<const IST*> (vals);
    return this->sumIntoLocalValuesImpl (curRowVals.data (), graph, rowInfo,
                                         cols, inputVals, numEnt, atomic);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
                    values_dualv_type::t_host::const_type
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getValuesViewHost (const RowInfo& rowinfo) const
  {
    if (rowinfo.allocSize == 0 || valuesUnpacked_wdv.extent(0) == 0)
      return typename values_dualv_type::t_host::const_type ();
    else
      return valuesUnpacked_wdv.getHostSubview(rowinfo.offset1D,
                                               rowinfo.allocSize,
                                               Access::ReadOnly);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
                    values_dualv_type::t_host
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getValuesViewHostNonConst (const RowInfo& rowinfo)
  {
    if (rowinfo.allocSize == 0 || valuesUnpacked_wdv.extent(0) == 0)
      return typename values_dualv_type::t_host ();
    else
      return valuesUnpacked_wdv.getHostSubview(rowinfo.offset1D,
                                               rowinfo.allocSize,
                                               Access::ReadWrite);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
                    values_dualv_type::t_dev::const_type
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getValuesViewDevice (const RowInfo& rowinfo) const
  {
    if (rowinfo.allocSize == 0 || valuesUnpacked_wdv.extent(0) == 0)
      return typename values_dualv_type::t_dev::const_type ();
    else
      return valuesUnpacked_wdv.getDeviceSubview(rowinfo.offset1D,
                                                 rowinfo.allocSize,
                                                 Access::ReadOnly);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
                    values_dualv_type::t_dev
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getValuesViewDeviceNonConst (const RowInfo& rowinfo)
  {
    if (rowinfo.allocSize == 0 || valuesUnpacked_wdv.extent(0) == 0)
      return typename values_dualv_type::t_dev ();
    else
      return valuesUnpacked_wdv.getDeviceSubview(rowinfo.offset1D,
                                                 rowinfo.allocSize,
                                                 Access::ReadWrite);
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalRowCopy (local_ordinal_type localRow,
                     nonconst_local_inds_host_view_type &indices,
                     nonconst_values_host_view_type &values,
                     size_t& numEntries) const 
 {
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    const char tfecfFuncName[] = "getLocalRowCopy: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->hasColMap (), std::runtime_error,
       "The matrix does not have a column Map yet.  This means we don't have "
       "local indices for columns yet, so it doesn't make sense to call this "
       "method.  If the matrix doesn't have a column Map yet, you should call "
       "fillComplete on it first.");

    const RowInfo rowinfo = staticGraph_->getRowInfo (localRow);
    const size_t theNumEntries = rowinfo.numEntries;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (static_cast<size_t> (indices.size ()) < theNumEntries ||
       static_cast<size_t> (values.size ()) < theNumEntries,
       std::runtime_error, "Row with local index " << localRow << " has " <<
       theNumEntries << " entry/ies, but indices.size() = " <<
       indices.size () << " and values.size() = " << values.size () << ".");
    numEntries = theNumEntries; // first side effect

    if (rowinfo.localRow != Teuchos::OrdinalTraits<size_t>::invalid ()) {
      if (staticGraph_->isLocallyIndexed ()) {
        auto curLclInds = staticGraph_->getLocalIndsViewHost(rowinfo);
        auto curVals = getValuesViewHost(rowinfo);

        for (size_t j = 0; j < theNumEntries; ++j) {
          values[j] = curVals[j];
          indices[j] = curLclInds(j);
        }
      }
      else if (staticGraph_->isGloballyIndexed ()) {
        // Don't call getColMap(), because it touches RCP's reference count.
        const map_type& colMap = * (staticGraph_->colMap_);
        auto curGblInds = staticGraph_->getGlobalIndsViewHost(rowinfo);
        auto curVals = getValuesViewHost(rowinfo);

        for (size_t j = 0; j < theNumEntries; ++j) {
          values[j] = curVals[j];
          indices[j] = colMap.getLocalElement (curGblInds(j));
        }
      }
    }
  }


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getGlobalRowCopy (global_ordinal_type globalRow,
                      nonconst_global_inds_host_view_type &indices,
                      nonconst_values_host_view_type &values,
                      size_t& numEntries) const
  {
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    const char tfecfFuncName[] = "getGlobalRowCopy: ";

    const RowInfo rowinfo =
      staticGraph_->getRowInfoFromGlobalRowIndex (globalRow);
    const size_t theNumEntries = rowinfo.numEntries;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (indices.size ()) < theNumEntries ||
      static_cast<size_t> (values.size ()) < theNumEntries,
      std::runtime_error, "Row with global index " << globalRow << " has "
      << theNumEntries << " entry/ies, but indices.size() = " <<
      indices.size () << " and values.size() = " << values.size () << ".");
    numEntries = theNumEntries; // first side effect

    if (rowinfo.localRow != Teuchos::OrdinalTraits<size_t>::invalid ()) {
      if (staticGraph_->isLocallyIndexed ()) {
        const map_type& colMap = * (staticGraph_->colMap_);
        auto curLclInds = staticGraph_->getLocalIndsViewHost(rowinfo);
        auto curVals = getValuesViewHost(rowinfo);

        for (size_t j = 0; j < theNumEntries; ++j) {
          values[j] = curVals[j];
          indices[j] = colMap.getGlobalElement (curLclInds(j));
        }
      }
      else if (staticGraph_->isGloballyIndexed ()) {
        auto curGblInds = staticGraph_->getGlobalIndsViewHost(rowinfo);
        auto curVals = getValuesViewHost(rowinfo);

        for (size_t j = 0; j < theNumEntries; ++j) {
          values[j] = curVals[j];
          indices[j] = curGblInds(j);
        }
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalRowView(LocalOrdinal localRow,
                  local_inds_host_view_type &indices,
                  values_host_view_type &values) const 
  {
    const char tfecfFuncName[] = "getLocalRowView: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isGloballyIndexed (), std::runtime_error, "The matrix currently stores "
      "its indices as global indices, so you cannot get a view with local "
      "column indices.  If the matrix has a column Map, you may call "
      "getLocalRowCopy() to get local column indices; otherwise, you may get "
      "a view with global column indices by calling getGlobalRowCopy().");

    const RowInfo rowInfo = staticGraph_->getRowInfo (localRow);
    if (rowInfo.localRow != Teuchos::OrdinalTraits<size_t>::invalid () &&
        rowInfo.numEntries > 0) {
      indices = staticGraph_->lclIndsUnpacked_wdv.getHostSubview(
                                                         rowInfo.offset1D,
                                                         rowInfo.numEntries,
                                                         Access::ReadOnly);
      values = valuesUnpacked_wdv.getHostSubview(rowInfo.offset1D,
                                                 rowInfo.numEntries,
                                                 Access::ReadOnly);
    }
    else {
      // This does the right thing (reports an empty row) if the input
      // row is invalid.
      indices = local_inds_host_view_type();
      values = values_host_view_type();
    }

#ifdef HAVE_TPETRA_DEBUG
    const char suffix[] = ".  This should never happen.  Please report this "
      "bug to the Tpetra developers.";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (static_cast<size_t> (indices.size ()) !=
       static_cast<size_t> (values.size ()), std::logic_error,
       "At the end of this method, for local row " << localRow << ", "
       "indices.size() = " << indices.size () << " != values.size () = "
       << values.size () << suffix);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (static_cast<size_t> (indices.size ()) !=
       static_cast<size_t> (rowInfo.numEntries), std::logic_error,
       "At the end of this method, for local row " << localRow << ", "
       "indices.size() = " << indices.size () << " != rowInfo.numEntries = "
       << rowInfo.numEntries << suffix);
    const size_t expectedNumEntries = getNumEntriesInLocalRow (localRow);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (rowInfo.numEntries != expectedNumEntries, std::logic_error, "At the end "
       "of this method, for local row " << localRow << ", rowInfo.numEntries = "
       << rowInfo.numEntries << " != getNumEntriesInLocalRow(localRow) = " <<
       expectedNumEntries << suffix);
#endif // HAVE_TPETRA_DEBUG
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getGlobalRowView (GlobalOrdinal globalRow,
                    global_inds_host_view_type &indices,
                    values_host_view_type &values) const
  {
    const char tfecfFuncName[] = "getGlobalRowView: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isLocallyIndexed (), std::runtime_error,
      "The matrix is locally indexed, so we cannot return a view of the row "
      "with global column indices.  Use getGlobalRowCopy() instead.");

    // This does the right thing (reports an empty row) if the input
    // row is invalid.
    const RowInfo rowInfo = 
          staticGraph_->getRowInfoFromGlobalRowIndex (globalRow);
    if (rowInfo.localRow != Teuchos::OrdinalTraits<size_t>::invalid () &&
        rowInfo.numEntries > 0) {
      indices = staticGraph_->gblInds_wdv.getHostSubview(rowInfo.offset1D,
                                                         rowInfo.numEntries,
                                                         Access::ReadOnly);
      values = valuesUnpacked_wdv.getHostSubview(rowInfo.offset1D,
                                                 rowInfo.numEntries,
                                                 Access::ReadOnly);
    }
    else {
      indices = global_inds_host_view_type();
      values = values_host_view_type();
    }

#ifdef HAVE_TPETRA_DEBUG
    const char suffix[] = ".  This should never happen.  Please report this "
      "bug to the Tpetra developers.";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (static_cast<size_t> (indices.size ()) !=
       static_cast<size_t> (values.size ()), std::logic_error,
       "At the end of this method, for global row " << globalRow << ", "
       "indices.size() = " << indices.size () << " != values.size () = "
       << values.size () << suffix);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (static_cast<size_t> (indices.size ()) !=
       static_cast<size_t> (rowInfo.numEntries), std::logic_error,
       "At the end of this method, for global row " << globalRow << ", "
       "indices.size() = " << indices.size () << " != rowInfo.numEntries = "
       << rowInfo.numEntries << suffix);
    const size_t expectedNumEntries = getNumEntriesInGlobalRow (globalRow);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (rowInfo.numEntries != expectedNumEntries, std::logic_error, "At the end "
       "of this method, for global row " << globalRow << ", rowInfo.numEntries "
       "= " << rowInfo.numEntries << " != getNumEntriesInGlobalRow(globalRow) ="
       " " << expectedNumEntries << suffix);
#endif // HAVE_TPETRA_DEBUG
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  scale (const Scalar& alpha)
  {
    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);

    const size_t nlrs = staticGraph_->getLocalNumRows ();
    const size_t numEntries = staticGraph_->getLocalNumEntries ();
    if (! staticGraph_->indicesAreAllocated () ||
        nlrs == 0 || numEntries == 0) {
      // do nothing
    }
    else {

      auto vals = valuesPacked_wdv.getDeviceView(Access::ReadWrite);
      KokkosBlas::scal(vals, theAlpha, vals);
   
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  setAllToScalar (const Scalar& alpha)
  {
    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);

    // replace all values in the matrix
    // it is easiest to replace all allocated values, instead of replacing only the ones with valid entries
    // however, if there are no valid entries, we can short-circuit
    // furthermore, if the values aren't allocated, we can short-circuit (no entry have been inserted so far)
    const size_t numEntries = staticGraph_->getLocalNumEntries();
    if (! staticGraph_->indicesAreAllocated () || numEntries == 0) {
      // do nothing
    }
    else {
      // DEEP_COPY REVIEW - VALUE-TO-DEVICE
      Kokkos::deep_copy (execution_space(), valuesUnpacked_wdv.getDeviceView(Access::OverwriteAll),
                         theAlpha);
      // CAG: This fence was found to be required on Cuda with UVM=on.
      Kokkos::fence("CrsMatrix::setAllToScalar");
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  setAllValues (const typename local_graph_device_type::row_map_type& rowPointers,
                const typename local_graph_device_type::entries_type::non_const_type& columnIndices,
                const typename local_matrix_device_type::values_type& values)
  {
    using ProfilingRegion=Details::ProfilingRegion;
    ProfilingRegion region ("Tpetra::CrsMatrix::setAllValues");
    const char tfecfFuncName[] = "setAllValues: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (columnIndices.size () != values.size (), std::invalid_argument,
       "columnIndices.size() = " << columnIndices.size () << " != values.size()"
       " = " << values.size () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (myGraph_.is_null (), std::runtime_error, "myGraph_ must not be null.");

    try {
      myGraph_->setAllIndices (rowPointers, columnIndices);
    }
    catch (std::exception &e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "myGraph_->setAllIndices() threw an "
         "exception: " << e.what ());
    }

    // Make sure that myGraph_ now has a local graph.  It may not be
    // fillComplete yet, so it's important to check.  We don't care
    // whether setAllIndices() did a shallow copy or a deep copy, so a
    // good way to check is to compare dimensions.
    auto lclGraph = myGraph_->getLocalGraphDevice ();
    const size_t numEnt = lclGraph.entries.extent (0);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (lclGraph.row_map.extent (0) != rowPointers.extent (0) ||
       numEnt != static_cast<size_t> (columnIndices.extent (0)),
       std::logic_error, "myGraph_->setAllIndices() did not correctly create "
       "local graph.  Please report this bug to the Tpetra developers.");

    valuesPacked_wdv = values_wdv_type(values);
    valuesUnpacked_wdv = valuesPacked_wdv;

    // Storage MUST be packed, since the interface doesn't give any
    // way to indicate any extra space at the end of each row.
    this->storageStatus_ = Details::STORAGE_1D_PACKED;

    checkInternalState ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  setAllValues ( const local_matrix_device_type& localDeviceMatrix)
  {
    using ProfilingRegion=Details::ProfilingRegion;
    ProfilingRegion region ("Tpetra::CrsMatrix::setAllValues from KokkosSparse::CrsMatrix");

    auto graph = localDeviceMatrix.graph;
    //FIXME how to check whether graph is allocated

    auto rows = graph.row_map;
    auto columns = graph.entries;
    auto values = localDeviceMatrix.values;

    setAllValues(rows,columns,values);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  setAllValues (const Teuchos::ArrayRCP<size_t>& ptr,
                const Teuchos::ArrayRCP<LocalOrdinal>& ind,
                const Teuchos::ArrayRCP<Scalar>& val)
  {
    using Kokkos::Compat::getKokkosViewDeepCopy;
    using Teuchos::ArrayRCP;
    using Teuchos::av_reinterpret_cast;
    typedef device_type DT;
    typedef impl_scalar_type IST;
    typedef typename local_graph_device_type::row_map_type row_map_type;
    //typedef typename row_map_type::non_const_value_type row_offset_type;
    const char tfecfFuncName[] = "setAllValues(ArrayRCP<size_t>, ArrayRCP<LO>, ArrayRCP<Scalar>): ";

    // The row offset type may depend on the execution space.  It may
    // not necessarily be size_t.  If it's not, we need to make a deep
    // copy.  We need to make a deep copy anyway so that Kokkos can
    // own the memory.  Regardless, ptrIn gets the copy.
    typename row_map_type::non_const_type ptrNative ("ptr", ptr.size ());
    Kokkos::View<const size_t*,
      typename row_map_type::array_layout,
      Kokkos::HostSpace,
      Kokkos::MemoryUnmanaged> ptrSizeT (ptr.getRawPtr (), ptr.size ());
    ::Tpetra::Details::copyOffsets (ptrNative, ptrSizeT);

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (ptrNative.extent (0) != ptrSizeT.extent (0),
       std::logic_error, "ptrNative.extent(0) = " <<
       ptrNative.extent (0) << " != ptrSizeT.extent(0) = "
       << ptrSizeT.extent (0) << ".  Please report this bug to the "
       "Tpetra developers.");

    auto indIn = getKokkosViewDeepCopy<DT> (ind ());
    auto valIn = getKokkosViewDeepCopy<DT> (av_reinterpret_cast<IST> (val ()));
    this->setAllValues (ptrNative, indIn, valIn);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalDiagOffsets (Teuchos::ArrayRCP<size_t>& offsets) const
  {
    const char tfecfFuncName[] = "getLocalDiagOffsets: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (staticGraph_.is_null (), std::runtime_error, "The matrix has no graph.");

    // mfh 11 May 2016: We plan to deprecate the ArrayRCP version of
    // this method in CrsGraph too, so don't call it (otherwise build
    // warnings will show up and annoy users).  Instead, copy results
    // in and out, if the memory space requires it.

    const size_t lclNumRows = staticGraph_->getLocalNumRows ();
    if (static_cast<size_t> (offsets.size ()) < lclNumRows) {
      offsets.resize (lclNumRows);
    }

    // The input ArrayRCP must always be a host pointer.  Thus, if
    // device_type::memory_space is Kokkos::HostSpace, it's OK for us
    // to write to that allocation directly as a Kokkos::View.
    if (std::is_same<memory_space, Kokkos::HostSpace>::value) {
      // It is always syntactically correct to assign a raw host
      // pointer to a device View, so this code will compile correctly
      // even if this branch never runs.
      typedef Kokkos::View<size_t*, device_type,
                           Kokkos::MemoryUnmanaged> output_type;
      output_type offsetsOut (offsets.getRawPtr (), lclNumRows);
      staticGraph_->getLocalDiagOffsets (offsetsOut);
    }
    else {
      Kokkos::View<size_t*, device_type> offsetsTmp ("diagOffsets", lclNumRows);
      staticGraph_->getLocalDiagOffsets (offsetsTmp);
      typedef Kokkos::View<size_t*, Kokkos::HostSpace,
                           Kokkos::MemoryUnmanaged> output_type;
      output_type offsetsOut (offsets.getRawPtr (), lclNumRows);
      // DEEP_COPY REVIEW - DEVICE-TO-HOST
      Kokkos::deep_copy (execution_space(), offsetsOut, offsetsTmp);
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalDiagCopy (Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& diag) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    const char tfecfFuncName[] = "getLocalDiagCopy (1-arg): ";
    typedef local_ordinal_type LO;


    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      staticGraph_.is_null (), std::runtime_error,
      "This method requires that the matrix have a graph.");
    auto rowMapPtr = this->getRowMap ();
    if (rowMapPtr.is_null () || rowMapPtr->getComm ().is_null ()) {
      // Processes on which the row Map or its communicator is null
      // don't participate.  Users shouldn't even call this method on
      // those processes.
      return;
    }
    auto colMapPtr = this->getColMap ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->hasColMap () || colMapPtr.is_null (), std::runtime_error,
       "This method requires that the matrix have a column Map.");
    const map_type& rowMap = * rowMapPtr;
    const map_type& colMap = * colMapPtr;
    const LO myNumRows = static_cast<LO> (this->getLocalNumRows ());

#ifdef HAVE_TPETRA_DEBUG
    // isCompatible() requires an all-reduce, and thus this check
    // should only be done in debug mode.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! diag.getMap ()->isCompatible (rowMap), std::runtime_error,
      "The input Vector's Map must be compatible with the CrsMatrix's row "
      "Map.  You may check this by using Map's isCompatible method: "
      "diag.getMap ()->isCompatible (A.getRowMap ());");
#endif // HAVE_TPETRA_DEBUG

    if (this->isFillComplete ()) {
      const auto D_lcl = diag.getLocalViewDevice(Access::OverwriteAll);
      // 1-D subview of the first (and only) column of D_lcl.
      const auto D_lcl_1d =
        Kokkos::subview (D_lcl, Kokkos::make_pair (LO (0), myNumRows), 0);

      const auto lclRowMap = rowMap.getLocalMap ();
      const auto lclColMap = colMap.getLocalMap ();
      using ::Tpetra::Details::getDiagCopyWithoutOffsets;
      (void) getDiagCopyWithoutOffsets (D_lcl_1d, lclRowMap,
                                        lclColMap,
                                        getLocalMatrixDevice ());
    }
    else {
      using ::Tpetra::Details::getLocalDiagCopyWithoutOffsetsNotFillComplete;
      (void) getLocalDiagCopyWithoutOffsetsNotFillComplete (diag, *this);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalDiagCopy (Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& diag,
                    const Kokkos::View<const size_t*, device_type,
                      Kokkos::MemoryUnmanaged>& offsets) const
  {
    typedef LocalOrdinal LO;

#ifdef HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "getLocalDiagCopy: ";
    const map_type& rowMap = * (this->getRowMap ());
    // isCompatible() requires an all-reduce, and thus this check
    // should only be done in debug mode.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! diag.getMap ()->isCompatible (rowMap), std::runtime_error,
      "The input Vector's Map must be compatible with (in the sense of Map::"
      "isCompatible) the CrsMatrix's row Map.");
#endif // HAVE_TPETRA_DEBUG

    // For now, we fill the Vector on the host and sync to device.
    // Later, we may write a parallel kernel that works entirely on
    // device.
    //
    // NOTE (mfh 21 Jan 2016): The host kernel here assumes UVM.  Once
    // we write a device kernel, it will not need to assume UVM.

    auto D_lcl = diag.getLocalViewDevice (Access::OverwriteAll);
    const LO myNumRows = static_cast<LO> (this->getLocalNumRows ());
    // Get 1-D subview of the first (and only) column of D_lcl.
    auto D_lcl_1d =
      Kokkos::subview (D_lcl, Kokkos::make_pair (LO (0), myNumRows), 0);

    KokkosSparse::getDiagCopy (D_lcl_1d, offsets,
                               getLocalMatrixDevice ());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalDiagCopy (Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& diag,
                    const Teuchos::ArrayView<const size_t>& offsets) const
  {
    using LO = LocalOrdinal;
    using host_execution_space = Kokkos::DefaultHostExecutionSpace;
    using IST = impl_scalar_type;

#ifdef HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "getLocalDiagCopy: ";
    const map_type& rowMap = * (this->getRowMap ());
    // isCompatible() requires an all-reduce, and thus this check
    // should only be done in debug mode.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! diag.getMap ()->isCompatible (rowMap), std::runtime_error,
      "The input Vector's Map must be compatible with (in the sense of Map::"
      "isCompatible) the CrsMatrix's row Map.");
#endif // HAVE_TPETRA_DEBUG

    // See #1510.  In case diag has already been marked modified on
    // device, we need to clear that flag, since the code below works
    // on host.
    //diag.clear_sync_state ();

    // For now, we fill the Vector on the host and sync to device.
    // Later, we may write a parallel kernel that works entirely on
    // device.
    auto lclVecHost = diag.getLocalViewHost(Access::OverwriteAll);
    // 1-D subview of the first (and only) column of lclVecHost.
    auto lclVecHost1d = Kokkos::subview (lclVecHost, Kokkos::ALL (), 0);

    using host_offsets_view_type =
      Kokkos::View<const size_t*, Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
    host_offsets_view_type h_offsets (offsets.getRawPtr (), offsets.size ());
    // Find the diagonal entries and put them in lclVecHost1d.
    using range_type = Kokkos::RangePolicy<host_execution_space, LO>;
    const LO myNumRows = static_cast<LO> (this->getLocalNumRows ());
    const size_t INV = Tpetra::Details::OrdinalTraits<size_t>::invalid ();

    auto rowPtrsPackedHost = staticGraph_->getRowPtrsPackedHost();
    auto valuesPackedHost = valuesPacked_wdv.getHostView(Access::ReadOnly);
    Kokkos::parallel_for
      ("Tpetra::CrsMatrix::getLocalDiagCopy",
       range_type (0, myNumRows),
       [&, INV, h_offsets] (const LO lclRow) { // Value capture is a workaround for cuda + gcc-7.2 compiler bug w/c++14
        lclVecHost1d(lclRow) = STS::zero (); // default value if no diag entry
        if (h_offsets[lclRow] != INV) {
          auto curRowOffset = rowPtrsPackedHost (lclRow);
          lclVecHost1d(lclRow) = 
            static_cast<IST> (valuesPackedHost(curRowOffset+h_offsets[lclRow]));
        }
      });
    //diag.sync_device ();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  leftScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x)
  {
    using ::Tpetra::Details::ProfilingRegion;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcpFromRef;
    using vec_type = Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    const char tfecfFuncName[] = "leftScale: ";

    ProfilingRegion region ("Tpetra::CrsMatrix::leftScale");

    RCP<const vec_type> xp;
    if (this->getRangeMap ()->isSameAs (* (x.getMap ()))) {
      // Take from Epetra: If we have a non-trivial exporter, we must
      // import elements that are permuted or are on other processors.
      auto exporter = this->getCrsGraphRef ().getExporter ();
      if (exporter.get () != nullptr) {
        RCP<vec_type> tempVec (new vec_type (this->getRowMap ()));
        tempVec->doImport (x, *exporter, REPLACE); // reverse mode
        xp = tempVec;
      }
      else {
        xp = rcpFromRef (x);
      }
    }
    else if (this->getRowMap ()->isSameAs (* (x.getMap ()))) {
      xp = rcpFromRef (x);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::invalid_argument, "x's Map must be the same as "
         "either the row Map or the range Map of the CrsMatrix.");
    }

    if (this->isFillComplete()) {
      auto x_lcl = xp->getLocalViewDevice (Access::ReadOnly);
      auto x_lcl_1d = Kokkos::subview (x_lcl, Kokkos::ALL (), 0);
      using ::Tpetra::Details::leftScaleLocalCrsMatrix;
      leftScaleLocalCrsMatrix (getLocalMatrixDevice (),
                               x_lcl_1d, false, false);
    }
    else {
      // 6/2020  Disallow leftScale of non-fillComplete matrices #7446
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "CrsMatrix::leftScale requires matrix to be"
         " fillComplete");
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  rightScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x)
  {
    using ::Tpetra::Details::ProfilingRegion;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcpFromRef;
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> vec_type;
    const char tfecfFuncName[] = "rightScale: ";

    ProfilingRegion region ("Tpetra::CrsMatrix::rightScale");

    RCP<const vec_type> xp;
    if (this->getDomainMap ()->isSameAs (* (x.getMap ()))) {
      // Take from Epetra: If we have a non-trivial exporter, we must
      // import elements that are permuted or are on other processors.
      auto importer = this->getCrsGraphRef ().getImporter ();
      if (importer.get () != nullptr) {
        RCP<vec_type> tempVec (new vec_type (this->getColMap ()));
        tempVec->doImport (x, *importer, REPLACE);
        xp = tempVec;
      }
      else {
        xp = rcpFromRef (x);
      }
    }
    else if (this->getColMap ()->isSameAs (* (x.getMap ()))) {
      xp = rcpFromRef (x);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "x's Map must be the same as "
         "either the domain Map or the column Map of the CrsMatrix.");
    }

    if (this->isFillComplete()) {
      auto x_lcl = xp->getLocalViewDevice (Access::ReadOnly);
      auto x_lcl_1d = Kokkos::subview (x_lcl, Kokkos::ALL (), 0);
      using ::Tpetra::Details::rightScaleLocalCrsMatrix;
      rightScaleLocalCrsMatrix (getLocalMatrixDevice (),
                                x_lcl_1d, false, false);
    }
    else {
      // 6/2020  Disallow rightScale of non-fillComplete matrices #7446
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "CrsMatrix::rightScale requires matrix to be"
         " fillComplete");
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mag_type
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getFrobeniusNorm () const
  {
    using Teuchos::ArrayView;
    using Teuchos::outArg;
    using Teuchos::REDUCE_SUM;
    using Teuchos::reduceAll;

    // FIXME (mfh 05 Aug 2014) Write a thread-parallel kernel for the
    // local part of this computation.  It could make sense to put
    // this operation in the Kokkos::CrsMatrix.

    // check the cache first
    mag_type mySum = STM::zero ();
    if (getLocalNumEntries() > 0) {
      if (isStorageOptimized ()) {
        // "Optimized" storage is packed storage.  That means we can
        // iterate in one pass through the 1-D values array.
        const size_t numEntries = getLocalNumEntries ();
        auto values = valuesPacked_wdv.getHostView(Access::ReadOnly);
        for (size_t k = 0; k < numEntries; ++k) {
          auto val = values[k];
          // Note (etp 06 Jan 2015) We need abs() here for composite types
          // (in general, if mag_type is on the left-hand-side, we need
          // abs() on the right-hand-side)
          const mag_type val_abs = STS::abs (val);
          mySum += val_abs * val_abs;
        }
      }
      else {
        const LocalOrdinal numRows =
          static_cast<LocalOrdinal> (this->getLocalNumRows ());
        for (LocalOrdinal r = 0; r < numRows; ++r) {
          const RowInfo rowInfo = myGraph_->getRowInfo (r);
          const size_t numEntries = rowInfo.numEntries;
          auto A_r = this->getValuesViewHost(rowInfo);
          for (size_t k = 0; k < numEntries; ++k) {
            const impl_scalar_type val = A_r[k];
            const mag_type val_abs = STS::abs (val);
            mySum += val_abs * val_abs;
          }
        }
      }
    }
    mag_type totalSum = STM::zero ();
    reduceAll<int, mag_type> (* (getComm ()), REDUCE_SUM,
                              mySum, outArg (totalSum));
    return STM::sqrt (totalSum);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceColMap (const Teuchos::RCP<const map_type>& newColMap)
  {
    const char tfecfFuncName[] = "replaceColMap: ";
    // FIXME (mfh 06 Aug 2014) What if the graph is locally indexed?
    // Then replacing the column Map might mean that we need to
    // reindex the column indices.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      myGraph_.is_null (), std::runtime_error,
      "This method does not work if the matrix has a const graph.  The whole "
      "idea of a const graph is that you are not allowed to change it, but "
      "this method necessarily must modify the graph, since the graph owns "
      "the matrix's column Map.");
    myGraph_->replaceColMap (newColMap);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  reindexColumns (crs_graph_type* const graph,
                  const Teuchos::RCP<const map_type>& newColMap,
                  const Teuchos::RCP<const import_type>& newImport,
                  const bool sortEachRow)
  {
    const char tfecfFuncName[] = "reindexColumns: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      graph == nullptr && myGraph_.is_null (), std::invalid_argument,
      "The input graph is null, but the matrix does not own its graph.");

    crs_graph_type& theGraph = (graph == nullptr) ? *myGraph_ : *graph;
    const bool sortGraph = false; // we'll sort graph & matrix together below

    theGraph.reindexColumns (newColMap, newImport, sortGraph);

    if (sortEachRow && theGraph.isLocallyIndexed () && ! theGraph.isSorted ()) {
      const LocalOrdinal lclNumRows =
        static_cast<LocalOrdinal> (theGraph.getLocalNumRows ());

      for (LocalOrdinal row = 0; row < lclNumRows; ++row) {

        const RowInfo rowInfo = theGraph.getRowInfo (row);
        auto lclColInds = theGraph.getLocalIndsViewHostNonConst (rowInfo);
        auto vals = this->getValuesViewHostNonConst (rowInfo);

        sort2 (lclColInds.data (),
               lclColInds.data () + rowInfo.numEntries,
               vals.data ());
      }
      theGraph.indicesAreSorted_ = true;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceDomainMap (const Teuchos::RCP<const map_type>& newDomainMap)
  {
    const char tfecfFuncName[] = "replaceDomainMap: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      myGraph_.is_null (), std::runtime_error,
      "This method does not work if the matrix has a const graph.  The whole "
      "idea of a const graph is that you are not allowed to change it, but this"
      " method necessarily must modify the graph, since the graph owns the "
      "matrix's domain Map and Import objects.");
    myGraph_->replaceDomainMap (newDomainMap);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceDomainMapAndImporter (const Teuchos::RCP<const map_type>& newDomainMap,
                               Teuchos::RCP<const import_type>& newImporter)
  {
    const char tfecfFuncName[] = "replaceDomainMapAndImporter: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      myGraph_.is_null (), std::runtime_error,
      "This method does not work if the matrix has a const graph.  The whole "
      "idea of a const graph is that you are not allowed to change it, but this"
      " method necessarily must modify the graph, since the graph owns the "
      "matrix's domain Map and Import objects.");
    myGraph_->replaceDomainMapAndImporter (newDomainMap, newImporter);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceRangeMap (const Teuchos::RCP<const map_type>& newRangeMap)
  {
    const char tfecfFuncName[] = "replaceRangeMap: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      myGraph_.is_null (), std::runtime_error,
      "This method does not work if the matrix has a const graph.  The whole "
      "idea of a const graph is that you are not allowed to change it, but this"
      " method necessarily must modify the graph, since the graph owns the "
      "matrix's domain Map and Import objects.");
    myGraph_->replaceRangeMap (newRangeMap);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceRangeMapAndExporter (const Teuchos::RCP<const map_type>& newRangeMap,
                              Teuchos::RCP<const export_type>& newExporter)
  {
    const char tfecfFuncName[] = "replaceRangeMapAndExporter: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      myGraph_.is_null (), std::runtime_error,
      "This method does not work if the matrix has a const graph.  The whole "
      "idea of a const graph is that you are not allowed to change it, but this"
      " method necessarily must modify the graph, since the graph owns the "
      "matrix's domain Map and Import objects.");
    myGraph_->replaceRangeMapAndExporter (newRangeMap, newExporter);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  insertNonownedGlobalValues (const GlobalOrdinal globalRow,
                              const Teuchos::ArrayView<const GlobalOrdinal>& indices,
                              const Teuchos::ArrayView<const Scalar>& values)
  {
    using Teuchos::Array;
    typedef GlobalOrdinal GO;
    typedef typename Array<GO>::size_type size_type;

    const size_type numToInsert = indices.size ();
    // Add the new data to the list of nonlocals.
    // This creates the arrays if they don't exist yet.
    std::pair<Array<GO>, Array<Scalar> >& curRow = nonlocals_[globalRow];
    Array<GO>& curRowInds = curRow.first;
    Array<Scalar>& curRowVals = curRow.second;
    const size_type newCapacity = curRowInds.size () + numToInsert;
    curRowInds.reserve (newCapacity);
    curRowVals.reserve (newCapacity);
    for (size_type k = 0; k < numToInsert; ++k) {
      curRowInds.push_back (indices[k]);
      curRowVals.push_back (values[k]);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  globalAssemble ()
  {
    using Details::Behavior;
    using Details::ProfilingRegion;
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::REDUCE_MAX;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using std::endl;
    typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
    //typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename Teuchos::Array<GO>::size_type size_type;
    const char tfecfFuncName[] = "globalAssemble: "; // for exception macro
    ProfilingRegion regionGlobalAssemble ("Tpetra::CrsMatrix::globalAssemble");

    const bool verbose = Behavior::verbose("CrsMatrix");
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsMatrix", "globalAssemble");
      std::ostringstream os;
      os << *prefix << "nonlocals_.size()=" << nonlocals_.size()
         << endl;
      std::cerr << os.str();
    }
    RCP<const Comm<int> > comm = getComm ();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! isFillActive (), std::runtime_error, "Fill must be active before "
       "you may call this method.");

    const size_t myNumNonlocalRows = nonlocals_.size ();

    // If no processes have nonlocal rows, then we don't have to do
    // anything.  Checking this is probably cheaper than constructing
    // the Map of nonlocal rows (see below) and noticing that it has
    // zero global entries.
    {
      const int iHaveNonlocalRows = (myNumNonlocalRows == 0) ? 0 : 1;
      int someoneHasNonlocalRows = 0;
      reduceAll<int, int> (*comm, REDUCE_MAX, iHaveNonlocalRows,
                           outArg (someoneHasNonlocalRows));
      if (someoneHasNonlocalRows == 0) {
        return; // no process has nonlocal rows, so nothing to do
      }
    }

    // 1. Create a list of the "nonlocal" rows on each process.  this
    //    requires iterating over nonlocals_, so while we do this,
    //    deduplicate the entries and get a count for each nonlocal
    //    row on this process.
    // 2. Construct a new row Map corresponding to those rows.  This
    //    Map is likely overlapping.  We know that the Map is not
    //    empty on all processes, because the above all-reduce and
    //    return exclude that case.

    RCP<const map_type> nonlocalRowMap;
    Teuchos::Array<size_t> numEntPerNonlocalRow (myNumNonlocalRows);
    {
      Teuchos::Array<GO> myNonlocalGblRows (myNumNonlocalRows);
      size_type curPos = 0;
      for (auto mapIter = nonlocals_.begin (); mapIter != nonlocals_.end ();
           ++mapIter, ++curPos) {
        myNonlocalGblRows[curPos] = mapIter->first;
        // Get the values and column indices by reference, since we
        // intend to change them in place (that's what "erase" does).
        Teuchos::Array<GO>& gblCols = (mapIter->second).first;
        Teuchos::Array<Scalar>& vals = (mapIter->second).second;

        // Sort both arrays jointly, using the column indices as keys,
        // then merge them jointly.  "Merge" here adds values
        // corresponding to the same column indices.  The first 2 args
        // of merge2 are output arguments that work just like the
        // return value of std::unique.
        sort2 (gblCols.begin (), gblCols.end (), vals.begin ());
        typename Teuchos::Array<GO>::iterator gblCols_newEnd;
        typename Teuchos::Array<Scalar>::iterator vals_newEnd;
        merge2 (gblCols_newEnd, vals_newEnd,
                gblCols.begin (), gblCols.end (),
                vals.begin (), vals.end ());
        gblCols.erase (gblCols_newEnd, gblCols.end ());
        vals.erase (vals_newEnd, vals.end ());
        numEntPerNonlocalRow[curPos] = gblCols.size ();
      }

      // Currently, Map requires that its indexBase be the global min
      // of all its global indices.  Map won't compute this for us, so
      // we must do it.  If our process has no nonlocal rows, set the
      // "min" to the max possible GO value.  This ensures that if
      // some process has at least one nonlocal row, then it will pick
      // that up as the min.  We know that at least one process has a
      // nonlocal row, since the all-reduce and return at the top of
      // this method excluded that case.
      GO myMinNonlocalGblRow = std::numeric_limits<GO>::max ();
      {
        auto iter = std::min_element (myNonlocalGblRows.begin (),
                                      myNonlocalGblRows.end ());
        if (iter != myNonlocalGblRows.end ()) {
          myMinNonlocalGblRow = *iter;
        }
      }
      GO gblMinNonlocalGblRow = 0;
      reduceAll<int, GO> (*comm, REDUCE_MIN, myMinNonlocalGblRow,
                          outArg (gblMinNonlocalGblRow));
      const GO indexBase = gblMinNonlocalGblRow;
      const global_size_t INV = Teuchos::OrdinalTraits<global_size_t>::invalid ();
      nonlocalRowMap = rcp (new map_type (INV, myNonlocalGblRows (), indexBase, comm));
    }

    // 3. Use the values and column indices for each nonlocal row, as
    //    stored in nonlocals_, to construct a CrsMatrix corresponding
    //    to nonlocal rows.  We have
    //    exact counts of the number of entries in each nonlocal row.

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Create nonlocal matrix" << endl;
      std::cerr << os.str();
    }
    RCP<crs_matrix_type> nonlocalMatrix =
      rcp (new crs_matrix_type (nonlocalRowMap, numEntPerNonlocalRow ()));
    {
      size_type curPos = 0;
      for (auto mapIter = nonlocals_.begin (); mapIter != nonlocals_.end ();
           ++mapIter, ++curPos) {
        const GO gblRow = mapIter->first;
        // Get values & column indices by ref, just to avoid copy.
        Teuchos::Array<GO>& gblCols = (mapIter->second).first;
        Teuchos::Array<Scalar>& vals = (mapIter->second).second;
        //const LO numEnt = static_cast<LO> (numEntPerNonlocalRow[curPos]);
        nonlocalMatrix->insertGlobalValues (gblRow, gblCols (), vals ());
      }
    }
    // There's no need to fill-complete the nonlocals matrix.
    // We just use it as a temporary container for the Export.

    // 4. If the original row Map is one to one, then we can Export
    //    directly from nonlocalMatrix into this.  Otherwise, we have
    //    to create a temporary matrix with a one-to-one row Map,
    //    Export into that, then Import from the temporary matrix into
    //    *this.

    auto origRowMap = this->getRowMap ();
    const bool origRowMapIsOneToOne = origRowMap->isOneToOne ();

    int isLocallyComplete = 1; // true by default

    if (origRowMapIsOneToOne) {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Original row Map is 1-to-1" << endl;
        std::cerr << os.str();
      }
      export_type exportToOrig (nonlocalRowMap, origRowMap);
      if (! exportToOrig.isLocallyComplete ()) {
        isLocallyComplete = 0;
      }
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "doExport from nonlocalMatrix" << endl;
        std::cerr << os.str();
      }
      this->doExport (*nonlocalMatrix, exportToOrig, Tpetra::ADD);
      // We're done at this point!
    }
    else {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Original row Map is NOT 1-to-1" << endl;
        std::cerr << os.str();
      }
      // If you ask a Map whether it is one to one, it does some
      // communication and stashes intermediate results for later use
      // by createOneToOne.  Thus, calling createOneToOne doesn't cost
      // much more then the original cost of calling isOneToOne.
      auto oneToOneRowMap = Tpetra::createOneToOne (origRowMap);
      export_type exportToOneToOne (nonlocalRowMap, oneToOneRowMap);
      if (! exportToOneToOne.isLocallyComplete ()) {
        isLocallyComplete = 0;
      }

      // Create a temporary matrix with the one-to-one row Map.
      //
      // TODO (mfh 09 Sep 2016, 12 Sep 2016) Estimate # entries in
      // each row, to avoid reallocation during the Export operation.
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Create & doExport into 1-to-1 matrix"
           << endl;
        std::cerr << os.str();
      }
      crs_matrix_type oneToOneMatrix (oneToOneRowMap, 0);
      // Export from matrix of nonlocals into the temp one-to-one matrix.
      oneToOneMatrix.doExport(*nonlocalMatrix, exportToOneToOne,
                              Tpetra::ADD);

      // We don't need the matrix of nonlocals anymore, so get rid of
      // it, to keep the memory high-water mark down.
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Free nonlocalMatrix" << endl;
        std::cerr << os.str();
      }
      nonlocalMatrix = Teuchos::null;

      // Import from the one-to-one matrix to the original matrix.
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "doImport from 1-to-1 matrix" << endl;
        std::cerr << os.str();
      }
      import_type importToOrig (oneToOneRowMap, origRowMap);
      this->doImport (oneToOneMatrix, importToOrig, Tpetra::ADD);
    }

    // It's safe now to clear out nonlocals_, since we've already
    // committed side effects to *this.  The standard idiom for
    // clearing a Container like std::map, is to swap it with an empty
    // Container and let the swapped Container fall out of scope.
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Free nonlocals_ (std::map)" << endl;
      std::cerr << os.str();
    }
    decltype (nonlocals_) newNonlocals;
    std::swap (nonlocals_, newNonlocals);

    // FIXME (mfh 12 Sep 2016) I don't like this all-reduce, and I
    // don't like throwing an exception here.  A local return value
    // would likely be more useful to users.  However, if users find
    // themselves exercising nonlocal inserts often, then they are
    // probably novice users who need the help.  See Gibhub Issues
    // #603 and #601 (esp. the latter) for discussion.

    int isGloballyComplete = 0; // output argument of reduceAll
    reduceAll<int, int> (*comm, REDUCE_MIN, isLocallyComplete,
                         outArg (isGloballyComplete));
    TEUCHOS_TEST_FOR_EXCEPTION
      (isGloballyComplete != 1, std::runtime_error, "On at least one process, "
       "you called insertGlobalValues with a global row index which is not in "
       "the matrix's row Map on any process in its communicator.");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  resumeFill (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    if (! isStaticGraph ()) { // Don't resume fill of a nonowned graph.
      myGraph_->resumeFill (params);
    }
#if KOKKOSKERNELS_VERSION >= 40299
    // Delete the apply helper (if it exists)
    applyHelper.reset();
#endif
    fillComplete_ = false;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  haveGlobalConstants() const {
    return getCrsGraphRef ().haveGlobalConstants ();
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  fillComplete (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    const char tfecfFuncName[] = "fillComplete(params): ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (this->getCrsGraph ().is_null (), std::logic_error,
      "getCrsGraph() returns null.  This should not happen at this point.  "
      "Please report this bug to the Tpetra developers.");

    const crs_graph_type& graph = this->getCrsGraphRef ();
    if (this->isStaticGraph () && graph.isFillComplete ()) {
      // If this matrix's graph is fill complete and the user did not
      // supply a domain or range Map, use the graph's domain and
      // range Maps.
      this->fillComplete (graph.getDomainMap (), graph.getRangeMap (), params);
    }
    else { // assume that user's row Map is the domain and range Map
      Teuchos::RCP<const map_type> rangeMap = graph.getRowMap ();
      Teuchos::RCP<const map_type> domainMap = rangeMap;
      this->fillComplete (domainMap, rangeMap, params);
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  fillComplete (const Teuchos::RCP<const map_type>& domainMap,
                const Teuchos::RCP<const map_type>& rangeMap,
                const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    using Details::Behavior;
    using Details::ProfilingRegion;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    const char tfecfFuncName[] = "fillComplete: ";
    ProfilingRegion regionFillComplete
      ("Tpetra::CrsMatrix::fillComplete");
    const bool verbose = Behavior::verbose("CrsMatrix");
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsMatrix", "fillComplete(dom,ran,p)");
      std::ostringstream os;
      os << *prefix << endl;
      std::cerr << os.str ();
    }
    Details::ProfilingRegion region(
      "Tpetra::CrsMatrix::fillCompete",
      "fillCompete");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->isFillActive () || this->isFillComplete (), std::runtime_error,
       "Matrix fill state must be active (isFillActive() "
       "must be true) before you may call fillComplete().");
    const int numProcs = this->getComm ()->getSize ();

    //
    // Read parameters from the input ParameterList.
    //
    {
      Details::ProfilingRegion region_fc("Tpetra::CrsMatrix::fillCompete", "ParameterList");

      // If true, the caller promises that no process did nonlocal
      // changes since the last call to fillComplete.
      bool assertNoNonlocalInserts = false;
      // If true, makeColMap sorts remote GIDs (within each remote
      // process' group).
      bool sortGhosts = true;

      if (! params.is_null ()) {
	assertNoNonlocalInserts = params->get ("No Nonlocal Changes",
					       assertNoNonlocalInserts);
	if (params->isParameter ("sort column map ghost gids")) {
	  sortGhosts = params->get ("sort column map ghost gids", sortGhosts);
	}
	else if (params->isParameter ("Sort column Map ghost GIDs")) {
	  sortGhosts = params->get ("Sort column Map ghost GIDs", sortGhosts);
	}
      }
      // We also don't need to do global assembly if there is only one
      // process in the communicator.
      const bool needGlobalAssemble = ! assertNoNonlocalInserts && numProcs > 1;
      // This parameter only matters if this matrix owns its graph.
      if (! this->myGraph_.is_null ()) {
	this->myGraph_->sortGhostsAssociatedWithEachProcessor_ = sortGhosts;
      }

      if (! this->getCrsGraphRef ().indicesAreAllocated ()) {
	if (this->hasColMap ()) { // use local indices
	  allocateValues(LocalIndices, GraphNotYetAllocated, verbose);
	}
	else { // no column Map, so use global indices
	  allocateValues(GlobalIndices, GraphNotYetAllocated, verbose);
	}
      }
      // Global assemble, if we need to.  This call only costs a single
      // all-reduce if we didn't need global assembly after all.
      if (needGlobalAssemble) {
	this->globalAssemble ();
      }
      else {
	TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
	  (numProcs == 1 && nonlocals_.size() > 0,
	   std::runtime_error, "Cannot have nonlocal entries on a serial run.  "
	   "An invalid entry (i.e., with row index not in the row Map) must have "
	   "been submitted to the CrsMatrix.");
      }
    }
    if (this->isStaticGraph ()) {
      Details::ProfilingRegion region_isg("Tpetra::CrsMatrix::fillCompete", "isStaticGraph");
      // FIXME (mfh 14 Nov 2016) In order to fix #843, I enable the
      // checks below only in debug mode.  It would be nicer to do a
      // local check, then propagate the error state in a deferred
      // way, whenever communication happens.  That would reduce the
      // cost of checking, to the point where it may make sense to
      // enable it even in release mode.
#ifdef HAVE_TPETRA_DEBUG
      // FIXME (mfh 18 Jun 2014) This check for correctness of the
      // input Maps incurs a penalty of two all-reduces for the
      // otherwise optimal const graph case.
      //
      // We could turn these (max) 2 all-reduces into (max) 1, by
      // fusing them.  We could do this by adding a "locallySameAs"
      // method to Map, which would return one of four states:
      //
      //   a. Certainly globally the same
      //   b. Certainly globally not the same
      //   c. Locally the same
      //   d. Locally not the same
      //
      // The first two states don't require further communication.
      // The latter two states require an all-reduce to communicate
      // globally, but we only need one all-reduce, since we only need
      // to check whether at least one of the Maps is wrong.
      const bool domainMapsMatch =
        this->staticGraph_->getDomainMap ()->isSameAs (*domainMap);
      const bool rangeMapsMatch =
        this->staticGraph_->getRangeMap ()->isSameAs (*rangeMap);

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! domainMapsMatch, std::runtime_error,
         "The CrsMatrix's domain Map does not match the graph's domain Map.  "
         "The graph cannot be changed because it was given to the CrsMatrix "
         "constructor as const.  You can fix this by passing in the graph's "
         "domain Map and range Map to the matrix's fillComplete call.");

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! rangeMapsMatch, std::runtime_error,
         "The CrsMatrix's range Map does not match the graph's range Map.  "
         "The graph cannot be changed because it was given to the CrsMatrix "
         "constructor as const.  You can fix this by passing in the graph's "
         "domain Map and range Map to the matrix's fillComplete call.");
#endif // HAVE_TPETRA_DEBUG

      // The matrix does _not_ own the graph, and the graph's
      // structure is already fixed, so just fill the local matrix.
      this->fillLocalMatrix (params);
    }
    else {
      Details::ProfilingRegion region_insg("Tpetra::CrsMatrix::fillCompete", "isNotStaticGraph");
      // Set the graph's domain and range Maps.  This will clear the
      // Import if the domain Map has changed (is a different
      // pointer), and the Export if the range Map has changed (is a
      // different pointer).
      this->myGraph_->setDomainRangeMaps (domainMap, rangeMap);

      // Make the graph's column Map, if necessary.
      Teuchos::Array<int> remotePIDs (0);
      const bool mustBuildColMap = ! this->hasColMap ();
      if (mustBuildColMap) {
        this->myGraph_->makeColMap (remotePIDs);
      }

      // Make indices local, if necessary.  The method won't do
      // anything if the graph is already locally indexed.
      const std::pair<size_t, std::string> makeIndicesLocalResult =
        this->myGraph_->makeIndicesLocal(verbose);
      // TODO (mfh 20 Jul 2017) Instead of throwing here, pass along
      // the error state to makeImportExport
      // which may do all-reduces and thus may
      // have the opportunity to communicate that error state.
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (makeIndicesLocalResult.first != 0, std::runtime_error,
         makeIndicesLocalResult.second);

      const bool sorted = this->myGraph_->isSorted ();
      const bool merged = this->myGraph_->isMerged ();
      this->sortAndMergeIndicesAndValues (sorted, merged);

      // Make Import and Export objects, if they haven't been made
      // already.  If we made a column Map above, reuse information
      // from that process to avoid communiation in the Import setup.
      this->myGraph_->makeImportExport (remotePIDs, mustBuildColMap);

      // The matrix _does_ own the graph, so fill the local graph at
      // the same time as the local matrix.
      this->fillLocalGraphAndMatrix (params);

      const bool callGraphComputeGlobalConstants = params.get () == nullptr ||
        params->get ("compute global constants", true);
      if (callGraphComputeGlobalConstants) {
        this->myGraph_->computeGlobalConstants ();
      }
      else {
        this->myGraph_->computeLocalConstants ();
      }
      this->myGraph_->fillComplete_ = true;
      this->myGraph_->checkInternalState ();
    }

    // FIXME (mfh 28 Aug 2014) "Preserve Local Graph" bool parameter no longer used.

    this->fillComplete_ = true; // Now we're fill complete!
    {
      Details::ProfilingRegion region_cis(
        "Tpetra::CrsMatrix::fillCompete", "checkInternalState"
      );
      this->checkInternalState ();
    }
  } //fillComplete(domainMap, rangeMap, params)

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  expertStaticFillComplete (const Teuchos::RCP<const map_type> & domainMap,
                            const Teuchos::RCP<const map_type> & rangeMap,
                            const Teuchos::RCP<const import_type>& importer,
                            const Teuchos::RCP<const export_type>& exporter,
                            const Teuchos::RCP<Teuchos::ParameterList> &params)
  {
#ifdef HAVE_TPETRA_MMM_TIMINGS
    std::string label;
    if(!params.is_null())
      label = params->get("Timer Label",label);
    std::string prefix = std::string("Tpetra ")+ label + std::string(": ");
    using Teuchos::TimeMonitor;

    Teuchos::TimeMonitor all(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-all")));
#endif

    const char tfecfFuncName[] = "expertStaticFillComplete: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( ! isFillActive() || isFillComplete(),
      std::runtime_error, "Matrix fill state must be active (isFillActive() "
      "must be true) before calling fillComplete().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      myGraph_.is_null (), std::logic_error, "myGraph_ is null.  This is not allowed.");

    {
#ifdef HAVE_TPETRA_MMM_TIMINGS
        Teuchos::TimeMonitor  graph(*TimeMonitor::getNewTimer(prefix + std::string("eSFC-M-Graph")));
#endif
        // We will presume globalAssemble is not needed, so we do the ESFC on the graph
        myGraph_->expertStaticFillComplete (domainMap, rangeMap, importer, exporter,params);
    }

    {
#ifdef HAVE_TPETRA_MMM_TIMINGS
        TimeMonitor  fLGAM(*TimeMonitor::getNewTimer(prefix + std::string("eSFC-M-fLGAM")));
#endif
        // Fill the local graph and matrix
        fillLocalGraphAndMatrix (params);
    }
    // FIXME (mfh 28 Aug 2014) "Preserve Local Graph" bool parameter no longer used.

    // Now we're fill complete!
    fillComplete_ = true;

    // Sanity checks at the end.
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isFillActive(), std::logic_error,
      ": We're at the end of fillComplete(), but isFillActive() is true.  "
      "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! isFillComplete(), std::logic_error,
      ": We're at the end of fillComplete(), but isFillActive() is true.  "
      "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
    {
#ifdef HAVE_TPETRA_MMM_TIMINGS
    Teuchos::TimeMonitor cIS(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-M-cIS")));
#endif

    checkInternalState();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  mergeRowIndicesAndValues (size_t rowLen, LocalOrdinal* cols, impl_scalar_type* vals)
  {
    impl_scalar_type* rowValueIter = vals;
    // beg,end define a half-exclusive interval over which to iterate.
    LocalOrdinal* beg = cols;
    LocalOrdinal* end = cols + rowLen;
    LocalOrdinal* newend = beg;
    if (beg != end) {
      LocalOrdinal* cur = beg + 1;
      impl_scalar_type* vcur = rowValueIter + 1;
      impl_scalar_type* vend = rowValueIter;
      cur = beg+1;
      while (cur != end) {
        if (*cur != *newend) {
          // new entry; save it
          ++newend;
          ++vend;
          (*newend) = (*cur);
          (*vend) = (*vcur);
        }
        else {
          // old entry; merge it
          //(*vend) = f (*vend, *vcur);
          (*vend) += *vcur;
        }
        ++cur;
        ++vcur;
      }
      ++newend; // one past the last entry, per typical [beg,end) semantics
    }
    return newend - beg;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  sortAndMergeIndicesAndValues (const bool sorted, const bool merged)
  {
    using ::Tpetra::Details::ProfilingRegion;
    typedef LocalOrdinal LO;
    typedef typename Kokkos::View<LO*, device_type>::HostMirror::execution_space
      host_execution_space;
    typedef Kokkos::RangePolicy<host_execution_space, LO> range_type;
    const char tfecfFuncName[] = "sortAndMergeIndicesAndValues: ";
    ProfilingRegion regionSAM ("Tpetra::CrsMatrix::sortAndMergeIndicesAndValues");

    if (! sorted || ! merged) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->isStaticGraph (), std::runtime_error, "Cannot sort or merge with "
         "\"static\" (const) graph, since the matrix does not own the graph.");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->myGraph_.is_null (), std::logic_error, "myGraph_ is null, but "
         "this matrix claims ! isStaticGraph().  "
         "Please report this bug to the Tpetra developers.");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->isStorageOptimized (), std::logic_error, "It is invalid to call "
         "this method if the graph's storage has already been optimized.  "
         "Please report this bug to the Tpetra developers.");

      crs_graph_type& graph = * (this->myGraph_);
      const LO lclNumRows = static_cast<LO> (this->getLocalNumRows ());
      size_t totalNumDups = 0;
      {
        //Accessing host unpacked (4-array CRS) local matrix.
        auto rowBegins_ = graph.getRowPtrsUnpackedHost();
        auto rowLengths_ = graph.k_numRowEntries_;
        auto vals_ = this->valuesUnpacked_wdv.getHostView(Access::ReadWrite);
        auto cols_ = graph.lclIndsUnpacked_wdv.getHostView(Access::ReadWrite);
        Kokkos::parallel_reduce ("sortAndMergeIndicesAndValues", range_type (0, lclNumRows),
          [=] (const LO lclRow, size_t& numDups) {
            size_t rowBegin = rowBegins_(lclRow);
            size_t rowLen = rowLengths_(lclRow);
            LO* cols = cols_.data() + rowBegin;
            impl_scalar_type* vals = vals_.data() + rowBegin;
            if (! sorted) {
              sort2 (cols, cols + rowLen, vals);
            }
            if (! merged) {
              size_t newRowLength = mergeRowIndicesAndValues (rowLen, cols, vals);
              rowLengths_(lclRow) = newRowLength;
              numDups += rowLen - newRowLength;
            }
          }, totalNumDups);
      }
      if (! sorted) {
        graph.indicesAreSorted_ = true; // we just sorted every row
      }
      if (! merged) {
        graph.noRedundancies_ = true; // we just merged every row
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  applyNonTranspose (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & X_in,
                     MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Y_in,
                     Scalar alpha,
                     Scalar beta) const
  {
    using Tpetra::Details::ProfilingRegion;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    using Teuchos::rcpFromRef;
    const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero ();
    const Scalar ONE = Teuchos::ScalarTraits<Scalar>::one ();

    // mfh 05 Jun 2014: Special case for alpha == 0.  I added this to
    // fix an Ifpack2 test (RILUKSingleProcessUnitTests), which was
    // failing only for the Kokkos refactor version of Tpetra.  It's a
    // good idea regardless to have the bypass.
    if (alpha == ZERO) {
      if (beta == ZERO) {
        Y_in.putScalar (ZERO);
      } else if (beta != ONE) {
        Y_in.scale (beta);
      }
      return;
    }

    // It's possible that X is a view of Y or vice versa.  We don't
    // allow this (apply() requires that X and Y not alias one
    // another), but it's helpful to detect and work around this case.
    // We don't try to to detect the more subtle cases (e.g., one is a
    // subview of the other, but their initial pointers differ).  We
    // only need to do this if this matrix's Import is trivial;
    // otherwise, we don't actually apply the operator from X into Y.

    RCP<const import_type> importer = this->getGraph ()->getImporter ();
    RCP<const export_type> exporter = this->getGraph ()->getExporter ();

    // If beta == 0, then the output MV will be overwritten; none of
    // its entries should be read.  (Sparse BLAS semantics say that we
    // must ignore any Inf or NaN entries in Y_in, if beta is zero.)
    // This matters if we need to do an Export operation; see below.
    const bool Y_is_overwritten = (beta == ZERO);

    // We treat the case of a replicated MV output specially.
    const bool Y_is_replicated =
      (! Y_in.isDistributed () && this->getComm ()->getSize () != 1);

    // This is part of the special case for replicated MV output.
    // We'll let each process do its thing, but do an all-reduce at
    // the end to sum up the results.  Setting beta=0 on all processes
    // but Proc 0 makes the math work out for the all-reduce.  (This
    // assumes that the replicated data is correctly replicated, so
    // that the data are the same on all processes.)
    if (Y_is_replicated && this->getComm ()->getRank () > 0) {
      beta = ZERO;
    }

    // Temporary MV for Import operation.  After the block of code
    // below, this will be an (Imported if necessary) column Map MV
    // ready to give to localApply(...).
    RCP<const MV> X_colMap;
    if (importer.is_null ()) {
      if (! X_in.isConstantStride ()) {
        // Not all sparse mat-vec kernels can handle an input MV with
        // nonconstant stride correctly, so we have to copy it in that
        // case into a constant stride MV.  To make a constant stride
        // copy of X_in, we force creation of the column (== domain)
        // Map MV (if it hasn't already been created, else fetch the
        // cached copy).  This avoids creating a new MV each time.
        RCP<MV> X_colMapNonConst = getColumnMapMultiVector (X_in, true);
        Tpetra::deep_copy (*X_colMapNonConst, X_in);
        X_colMap = rcp_const_cast<const MV> (X_colMapNonConst);
      }
      else {
        // The domain and column Maps are the same, so do the local
        // multiply using the domain Map input MV X_in.
        X_colMap = rcpFromRef (X_in);
      }
    }
    else { // need to Import source (multi)vector
      ProfilingRegion regionImport ("Tpetra::CrsMatrix::apply: Import");

      // We're doing an Import anyway, which will copy the relevant
      // elements of the domain Map MV X_in into a separate column Map
      // MV.  Thus, we don't have to worry whether X_in is constant
      // stride.
      RCP<MV> X_colMapNonConst = getColumnMapMultiVector (X_in);

      // Import from the domain Map MV to the column Map MV.
      X_colMapNonConst->doImport (X_in, *importer, INSERT);
      X_colMap = rcp_const_cast<const MV> (X_colMapNonConst);
    }

    // Temporary MV for doExport (if needed), or for copying a
    // nonconstant stride output MV into a constant stride MV.  This
    // is null if we don't need the temporary MV, that is, if the
    // Export is trivial (null).
    RCP<MV> Y_rowMap = getRowMapMultiVector (Y_in);

    // If we have a nontrivial Export object, we must perform an
    // Export.  In that case, the local multiply result will go into
    // the row Map multivector.  We don't have to make a
    // constant-stride version of Y_in in this case, because we had to
    // make a constant stride Y_rowMap MV and do an Export anyway.
    if (! exporter.is_null ()) {
      this->localApply (*X_colMap, *Y_rowMap, Teuchos::NO_TRANS, alpha, ZERO);
      {
        ProfilingRegion regionExport ("Tpetra::CrsMatrix::apply: Export");

        // If we're overwriting the output MV Y_in completely (beta ==
        // 0), then make sure that it is filled with zeros before we
        // do the Export.  Otherwise, the ADD combine mode will use
        // data in Y_in, which is supposed to be zero.
        if (Y_is_overwritten) {
          Y_in.putScalar (ZERO);
        }
        else {
          // Scale output MV by beta, so that doExport sums in the
          // mat-vec contribution: Y_in = beta*Y_in + alpha*A*X_in.
          Y_in.scale (beta);
        }
        // Do the Export operation.
        Y_in.doExport (*Y_rowMap, *exporter, ADD_ASSIGN);
      }
    }
    else { // Don't do an Export: row Map and range Map are the same.
      //
      // If Y_in does not have constant stride, or if the column Map
      // MV aliases Y_in, then we can't let the kernel write directly
      // to Y_in.  Instead, we have to use the cached row (== range)
      // Map MV as temporary storage.
      //
      // FIXME (mfh 05 Jun 2014) This test for aliasing only tests if
      // the user passed in the same MultiVector for both X and Y.  It
      // won't detect whether one MultiVector views the other.  We
      // should also check the MultiVectors' raw data pointers.
      if (! Y_in.isConstantStride () || X_colMap.getRawPtr () == &Y_in) {
        // Force creating the MV if it hasn't been created already.
        // This will reuse a previously created cached MV.
        Y_rowMap = getRowMapMultiVector (Y_in, true);

        // If beta == 0, we don't need to copy Y_in into Y_rowMap,
        // since we're overwriting it anyway.
        if (beta != ZERO) {
          Tpetra::deep_copy (*Y_rowMap, Y_in);
        }
        this->localApply (*X_colMap, *Y_rowMap, Teuchos::NO_TRANS, alpha, beta);
        Tpetra::deep_copy (Y_in, *Y_rowMap);
      }
      else {
        this->localApply (*X_colMap, Y_in, Teuchos::NO_TRANS, alpha, beta);
      }
    }

    // If the range Map is a locally replicated Map, sum up
    // contributions from each process.  We set beta = 0 on all
    // processes but Proc 0 initially, so this will handle the scaling
    // factor beta correctly.
    if (Y_is_replicated) {
      ProfilingRegion regionReduce ("Tpetra::CrsMatrix::apply: Reduce Y");
      Y_in.reduce ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  applyTranspose (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X_in,
                  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y_in,
                  const Teuchos::ETransp mode,
                  Scalar alpha,
                  Scalar beta) const
  {
    using Tpetra::Details::ProfilingRegion;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    using Teuchos::rcpFromRef;
    const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero ();

    // Take shortcuts for alpha == 0.
    if (alpha == ZERO) {
      // Follow the Sparse BLAS convention by ignoring both the matrix
      // and X_in, in this case.
      if (beta == ZERO) {
        // Follow the Sparse BLAS convention by overwriting any Inf or
        // NaN values in Y_in, in this case.
        Y_in.putScalar (ZERO);
      }
      else {
        Y_in.scale (beta);
      }
      return;
    }
    else if (beta == ZERO) {
      //Thyra was implicitly assuming that Y gets set to zero / or is overwritten
      //when bets==0. This was not the case with transpose in a multithreaded
      //environment where a multiplication with subsequent atomic_adds is used
      //since 0 is effectively not special cased. Doing the explicit set to zero here
      //This catches cases where Y is nan or inf.
      Y_in.putScalar (ZERO);
    }

    const size_t numVectors = X_in.getNumVectors ();

    // We don't allow X_in and Y_in to alias one another.  It's hard
    // to check this, because advanced users could create views from
    // raw pointers.  However, if X_in and Y_in reference the same
    // object, we will do the user a favor by copying X into new
    // storage (with a warning).  We only need to do this if we have
    // trivial importers; otherwise, we don't actually apply the
    // operator from X into Y.
    RCP<const import_type> importer = this->getGraph ()->getImporter ();
    RCP<const export_type> exporter = this->getGraph ()->getExporter ();
    // access X indirectly, in case we need to create temporary storage
    RCP<const MV> X;

    // some parameters for below
    const bool Y_is_replicated = (! Y_in.isDistributed () && this->getComm ()->getSize () != 1);
    const bool Y_is_overwritten = (beta == ZERO);
    if (Y_is_replicated && this->getComm ()->getRank () > 0) {
      beta = ZERO;
    }

    // The kernels do not allow input or output with nonconstant stride.
    if (! X_in.isConstantStride () && importer.is_null ()) {
      X = rcp (new MV (X_in, Teuchos::Copy)); // Constant-stride copy of X_in
    } else {
      X = rcpFromRef (X_in); // Reference to X_in
    }

    // Set up temporary multivectors for Import and/or Export.
    if (importer != Teuchos::null) {
      if (importMV_ != Teuchos::null && importMV_->getNumVectors() != numVectors) {
        importMV_ = null;
      }
      if (importMV_ == null) {
        importMV_ = rcp (new MV (this->getColMap (), numVectors));
      }
    }
    if (exporter != Teuchos::null) {
      if (exportMV_ != Teuchos::null && exportMV_->getNumVectors() != numVectors) {
        exportMV_ = null;
      }
      if (exportMV_ == null) {
        exportMV_ = rcp (new MV (this->getRowMap (), numVectors));
      }
    }

    // If we have a non-trivial exporter, we must import elements that
    // are permuted or are on other processors.
    if (! exporter.is_null ()) {
      ProfilingRegion regionImport ("Tpetra::CrsMatrix::apply (transpose): Import");
      exportMV_->doImport (X_in, *exporter, INSERT);
      X = exportMV_; // multiply out of exportMV_
    }

    // If we have a non-trivial importer, we must export elements that
    // are permuted or belong to other processors.  We will compute
    // solution into the to-be-exported MV; get a view.
    if (importer != Teuchos::null) {
      ProfilingRegion regionExport ("Tpetra::CrsMatrix::apply (transpose): Export");

      // FIXME (mfh 18 Apr 2015) Temporary fix suggested by Clark
      // Dohrmann on Fri 17 Apr 2015.  At some point, we need to go
      // back and figure out why this helps.  importMV_ SHOULD be
      // completely overwritten in the localApply(...) call
      // below, because beta == ZERO there.
      importMV_->putScalar (ZERO);
      // Do the local computation.
      this->localApply (*X, *importMV_, mode, alpha, ZERO);

      if (Y_is_overwritten) {
        Y_in.putScalar (ZERO);
      } else {
        Y_in.scale (beta);
      }
      Y_in.doExport (*importMV_, *importer, ADD_ASSIGN);
    }
    // otherwise, multiply into Y
    else {
      // can't multiply in-situ; can't multiply into non-strided multivector
      //
      // FIXME (mfh 05 Jun 2014) This test for aliasing only tests if
      // the user passed in the same MultiVector for both X and Y.  It
      // won't detect whether one MultiVector views the other.  We
      // should also check the MultiVectors' raw data pointers.
      if (! Y_in.isConstantStride () || X.getRawPtr () == &Y_in) {
        // Make a deep copy of Y_in, into which to write the multiply result.
        MV Y (Y_in, Teuchos::Copy);
        this->localApply (*X, Y, mode, alpha, beta);
        Tpetra::deep_copy (Y_in, Y);
      } else {
        this->localApply (*X, Y_in, mode, alpha, beta);
      }
    }

    // If the range Map is a locally replicated map, sum the
    // contributions from each process.  (That's why we set beta=0
    // above for all processes but Proc 0.)
    if (Y_is_replicated) {
      ProfilingRegion regionReduce ("Tpetra::CrsMatrix::apply (transpose): Reduce Y");
      Y_in.reduce ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  localApply (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
              MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
              const Teuchos::ETransp mode,
              const Scalar& alpha,
              const Scalar& beta) const
  {
    using Tpetra::Details::ProfilingRegion;
    using Teuchos::NO_TRANS;
    ProfilingRegion regionLocalApply ("Tpetra::CrsMatrix::localApply");

    auto X_lcl = X.getLocalViewDevice(Access::ReadOnly);
    auto Y_lcl = Y.getLocalViewDevice(Access::ReadWrite);

    const bool debug = ::Tpetra::Details::Behavior::debug ();
    if (debug) {
      const char tfecfFuncName[] = "localApply: ";
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (X.getNumVectors () != Y.getNumVectors (), std::runtime_error,
         "X.getNumVectors() = " << X.getNumVectors () << " != "
         "Y.getNumVectors() = " << Y.getNumVectors () << ".");
      const bool transpose = (mode != Teuchos::NO_TRANS);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! transpose && X.getLocalLength () !=
         getColMap ()->getLocalNumElements (), std::runtime_error,
         "NO_TRANS case: X has the wrong number of local rows.  "
         "X.getLocalLength() = " << X.getLocalLength () << " != "
         "getColMap()->getLocalNumElements() = " <<
         getColMap ()->getLocalNumElements () << ".");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! transpose && Y.getLocalLength () !=
         getRowMap ()->getLocalNumElements (), std::runtime_error,
         "NO_TRANS case: Y has the wrong number of local rows.  "
         "Y.getLocalLength() = " << Y.getLocalLength () << " != "
         "getRowMap()->getLocalNumElements() = " <<
         getRowMap ()->getLocalNumElements () << ".");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (transpose && X.getLocalLength () !=
         getRowMap ()->getLocalNumElements (), std::runtime_error,
         "TRANS or CONJ_TRANS case: X has the wrong number of local "
         "rows.  X.getLocalLength() = " << X.getLocalLength ()
         << " != getRowMap()->getLocalNumElements() = "
         << getRowMap ()->getLocalNumElements () << ".");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (transpose && Y.getLocalLength () !=
         getColMap ()->getLocalNumElements (), std::runtime_error,
         "TRANS or CONJ_TRANS case: X has the wrong number of local "
         "rows.  Y.getLocalLength() = " << Y.getLocalLength ()
         << " != getColMap()->getLocalNumElements() = "
         << getColMap ()->getLocalNumElements () << ".");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! isFillComplete (), std::runtime_error, "The matrix is not "
         "fill complete.  You must call fillComplete() (possibly with "
         "domain and range Map arguments) without an intervening "
         "resumeFill() call before you may call this method.");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! X.isConstantStride () || ! Y.isConstantStride (),
         std::runtime_error, "X and Y must be constant stride.");
      // If the two pointers are null, then they don't alias one
      // another, even though they are equal.
      // Kokkos does not guarantee that zero row-extent vectors 
      // point to different places, so we have to check that too.
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (X_lcl.data () == Y_lcl.data () && X_lcl.data () != nullptr
         && X_lcl.extent(0) != 0,
         std::runtime_error, "X and Y may not alias one another.");
    }

#if KOKKOSKERNELS_VERSION >= 40299
    auto A_lcl = getLocalMatrixDevice();

    if(!applyHelper.get()) {
      // The apply helper does not exist, so create it.
      // Decide now whether to use the imbalanced row path, or the default.
      bool useMergePath = false;
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
      //TODO: when https://github.com/kokkos/kokkos-kernels/issues/2166 is fixed and,
      //we can use SPMV_MERGE_PATH for the native spmv as well.
      //Take out this ifdef to enable that.
      //
      //Until then, only use SPMV_MERGE_PATH when calling cuSPARSE.
      if constexpr(std::is_same_v<execution_space, Kokkos::Cuda>) {
        LocalOrdinal nrows = getLocalNumRows();
        LocalOrdinal maxRowImbalance = 0;
        if(nrows != 0)
          maxRowImbalance = getLocalMaxNumRowEntries() - (getLocalNumEntries() / nrows);

        if(size_t(maxRowImbalance) >= Tpetra::Details::Behavior::rowImbalanceThreshold())
          useMergePath = true;
      }
#endif
      applyHelper = std::make_shared<ApplyHelper>(A_lcl.nnz(), A_lcl.graph.row_map,
          useMergePath ? KokkosSparse::SPMV_MERGE_PATH : KokkosSparse::SPMV_DEFAULT);
    }

    // Translate mode (Teuchos enum) to KokkosKernels (1-character string)
    const char* modeKK = nullptr;
    switch(mode)
    {
      case Teuchos::NO_TRANS:
        modeKK = KokkosSparse::NoTranspose;        break;
      case Teuchos::TRANS:
        modeKK = KokkosSparse::Transpose;          break;
      case Teuchos::CONJ_TRANS:
        modeKK = KokkosSparse::ConjugateTranspose; break;
      default:
        throw std::invalid_argument("Tpetra::CrsMatrix::localApply: invalid mode");
    }

    if(applyHelper->shouldUseIntRowptrs())
    {
      auto A_lcl_int_rowptrs = applyHelper->getIntRowptrMatrix(A_lcl);
      KokkosSparse::spmv(
          &applyHelper->handle_int, modeKK,
          impl_scalar_type(alpha), A_lcl_int_rowptrs, X_lcl, impl_scalar_type(beta), Y_lcl);
    }
    else
    {
      KokkosSparse::spmv(
          &applyHelper->handle, modeKK,
          impl_scalar_type(alpha), A_lcl, X_lcl, impl_scalar_type(beta), Y_lcl);
    }
#else
    LocalOrdinal nrows = getLocalNumRows();
    LocalOrdinal maxRowImbalance = 0;
    if(nrows != 0)
      maxRowImbalance = getLocalMaxNumRowEntries() - (getLocalNumEntries() / nrows);

    auto matrix_lcl = getLocalMultiplyOperator();
    if(size_t(maxRowImbalance) >= Tpetra::Details::Behavior::rowImbalanceThreshold())
      matrix_lcl->applyImbalancedRows (X_lcl, Y_lcl, mode, alpha, beta);
    else
      matrix_lcl->apply (X_lcl, Y_lcl, mode, alpha, beta);
#endif
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  apply (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &X,
         MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &Y,
         Teuchos::ETransp mode,
         Scalar alpha,
         Scalar beta) const
  {
    using Tpetra::Details::ProfilingRegion;
    const char fnName[] = "Tpetra::CrsMatrix::apply";

    TEUCHOS_TEST_FOR_EXCEPTION
      (! isFillComplete (), std::runtime_error,
       fnName << ": Cannot call apply() until fillComplete() "
       "has been called.");

    if (mode == Teuchos::NO_TRANS) {
      ProfilingRegion regionNonTranspose (fnName);
      this->applyNonTranspose (X, Y, alpha, beta);
    }
    else {
      ProfilingRegion regionTranspose ("Tpetra::CrsMatrix::apply (transpose)");
      this->applyTranspose (X, Y, mode, alpha, beta);
    }
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  template<class T>
  Teuchos::RCP<CrsMatrix<T, LocalOrdinal, GlobalOrdinal, Node> >
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  convert () const
  {
    using Teuchos::RCP;
    typedef CrsMatrix<T, LocalOrdinal, GlobalOrdinal, Node> output_matrix_type;
    const char tfecfFuncName[] = "convert: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->isFillComplete (), std::runtime_error, "This matrix (the source "
       "of the conversion) is not fill complete.  You must first call "
       "fillComplete() (possibly with the domain and range Map) without an "
       "intervening call to resumeFill(), before you may call this method.");

    RCP<output_matrix_type> newMatrix
      (new output_matrix_type (this->getCrsGraph ()));
    // Copy old values into new values.  impl_scalar_type and T may
    // differ, so we can't use Kokkos::deep_copy.
    using ::Tpetra::Details::copyConvert;
    copyConvert (newMatrix->getLocalMatrixDevice ().values,
                 this->getLocalMatrixDevice ().values);
    // Since newmat has a static (const) graph, the graph already has
    // a column Map, and Import and Export objects already exist (if
    // applicable).  Thus, calling fillComplete is cheap.
    newMatrix->fillComplete (this->getDomainMap (), this->getRangeMap ());

    return newMatrix;
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  checkInternalState () const
  {
    const bool debug = ::Tpetra::Details::Behavior::debug ("CrsGraph");
    if (debug) {
      const char tfecfFuncName[] = "checkInternalState: ";
      const char err[] = "Internal state is not consistent.  "
        "Please report this bug to the Tpetra developers.";

      // This version of the graph (RCP<const crs_graph_type>) must
      // always be nonnull.
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (staticGraph_.is_null (), std::logic_error, err);
      // myGraph == null means that the matrix has a const ("static")
      // graph.  Otherwise, the matrix has a dynamic graph (it owns its
      // graph).
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! myGraph_.is_null () && myGraph_ != staticGraph_,
         std::logic_error, err);
      // if matrix is fill complete, then graph must be fill complete
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (isFillComplete () && ! staticGraph_->isFillComplete (),
         std::logic_error, err << "  Specifically, the matrix is fill complete, "
         "but its graph is NOT fill complete.");
      // if values are allocated and they are non-zero in number, then
      // one of the allocations should be present
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (staticGraph_->indicesAreAllocated () &&
         staticGraph_->getLocalAllocationSize() > 0 &&
         staticGraph_->getLocalNumRows() > 0 &&
         valuesUnpacked_wdv.extent (0) == 0,
         std::logic_error, err);
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  description () const
  {
    std::ostringstream os;

    os << "Tpetra::CrsMatrix (Kokkos refactor): {";
    if (this->getObjectLabel () != "") {
      os << "Label: \"" << this->getObjectLabel () << "\", ";
    }
    if (isFillComplete ()) {
      os << "isFillComplete: true"
         << ", global dimensions: [" << getGlobalNumRows () << ", "
         << getGlobalNumCols () << "]"
         << ", global number of entries: " << getGlobalNumEntries ()
         << "}";
    }
    else {
      os << "isFillComplete: false"
         << ", global dimensions: [" << getGlobalNumRows () << ", "
         << getGlobalNumCols () << "]}";
    }
    return os.str ();
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using std::endl;
    using std::setw;
    using Teuchos::ArrayView;
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::TypeNameTraits;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;

    const Teuchos::EVerbosityLevel vl = (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

    if (vl == VERB_NONE) {
      return; // Don't print anything at all
    }

    // By convention, describe() always begins with a tab.
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = this->getComm();
    const int myRank = comm->getRank();
    const int numProcs = comm->getSize();
    size_t width = 1;
    for (size_t dec=10; dec<getGlobalNumRows(); dec *= 10) {
      ++width;
    }
    width = std::max<size_t> (width, static_cast<size_t> (11)) + 2;

    //    none: print nothing
    //     low: print O(1) info from node 0
    //  medium: print O(P) info, num entries per process
    //    high: print O(N) info, num entries per row
    // extreme: print O(NNZ) info: print indices and values
    //
    // for medium and higher, print constituent objects at specified verbLevel
    if (myRank == 0) {
      out << "Tpetra::CrsMatrix (Kokkos refactor):" << endl;
    }
    Teuchos::OSTab tab1 (out);

    if (myRank == 0) {
      if (this->getObjectLabel () != "") {
        out << "Label: \"" << this->getObjectLabel () << "\", ";
      }
      {
        out << "Template parameters:" << endl;
        Teuchos::OSTab tab2 (out);
        out << "Scalar: " << TypeNameTraits<Scalar>::name () << endl
            << "LocalOrdinal: " << TypeNameTraits<LocalOrdinal>::name () << endl
            << "GlobalOrdinal: " << TypeNameTraits<GlobalOrdinal>::name () << endl
            << "Node: " << TypeNameTraits<Node>::name () << endl;
      }
      if (isFillComplete()) {
        out << "isFillComplete: true" << endl
            << "Global dimensions: [" << getGlobalNumRows () << ", "
            << getGlobalNumCols () << "]" << endl
            << "Global number of entries: " << getGlobalNumEntries () << endl
            << endl << "Global max number of entries in a row: "
            << getGlobalMaxNumRowEntries () << endl;
      }
      else {
        out << "isFillComplete: false" << endl
            << "Global dimensions: [" << getGlobalNumRows () << ", "
            << getGlobalNumCols () << "]" << endl;
      }
    }

    if (vl < VERB_MEDIUM) {
      return; // all done!
    }

    // Describe the row Map.
    if (myRank == 0) {
      out << endl << "Row Map:" << endl;
    }
    if (getRowMap ().is_null ()) {
      if (myRank == 0) {
        out << "null" << endl;
      }
    }
    else {
      if (myRank == 0) {
        out << endl;
      }
      getRowMap ()->describe (out, vl);
    }

    // Describe the column Map.
    if (myRank == 0) {
      out << "Column Map: ";
    }
    if (getColMap ().is_null ()) {
      if (myRank == 0) {
        out << "null" << endl;
      }
    } else if (getColMap () == getRowMap ()) {
      if (myRank == 0) {
        out << "same as row Map" << endl;
      }
    } else {
      if (myRank == 0) {
        out << endl;
      }
      getColMap ()->describe (out, vl);
    }

    // Describe the domain Map.
    if (myRank == 0) {
      out << "Domain Map: ";
    }
    if (getDomainMap ().is_null ()) {
      if (myRank == 0) {
        out << "null" << endl;
      }
    } else if (getDomainMap () == getRowMap ()) {
      if (myRank == 0) {
        out << "same as row Map" << endl;
      }
    } else if (getDomainMap () == getColMap ()) {
      if (myRank == 0) {
        out << "same as column Map" << endl;
      }
    } else {
      if (myRank == 0) {
        out << endl;
      }
      getDomainMap ()->describe (out, vl);
    }

    // Describe the range Map.
    if (myRank == 0) {
      out << "Range Map: ";
    }
    if (getRangeMap ().is_null ()) {
      if (myRank == 0) {
        out << "null" << endl;
      }
    } else if (getRangeMap () == getDomainMap ()) {
      if (myRank == 0) {
        out << "same as domain Map" << endl;
      }
    } else if (getRangeMap () == getRowMap ()) {
      if (myRank == 0) {
        out << "same as row Map" << endl;
      }
    } else {
      if (myRank == 0) {
        out << endl;
      }
      getRangeMap ()->describe (out, vl);
    }

    // O(P) data
    for (int curRank = 0; curRank < numProcs; ++curRank) {
      if (myRank == curRank) {
        out << "Process rank: " << curRank << endl;
        Teuchos::OSTab tab2 (out);
        if (! staticGraph_->indicesAreAllocated ()) {
          out << "Graph indices not allocated" << endl;
        }
        else {
          out << "Number of allocated entries: "
              << staticGraph_->getLocalAllocationSize () << endl;
        }
        out << "Number of entries: " << getLocalNumEntries () << endl
            << "Max number of entries per row: " << getLocalMaxNumRowEntries ()
            << endl;
      }
      // Give output time to complete by executing some barriers.
      comm->barrier ();
      comm->barrier ();
      comm->barrier ();
    }

    if (vl < VERB_HIGH) {
      return; // all done!
    }

    // O(N) and O(NNZ) data
    for (int curRank = 0; curRank < numProcs; ++curRank) {
      if (myRank == curRank) {
        out << std::setw(width) << "Proc Rank"
            << std::setw(width) << "Global Row"
            << std::setw(width) << "Num Entries";
        if (vl == VERB_EXTREME) {
          out << std::setw(width) << "(Index,Value)";
        }
        out << endl;
        for (size_t r = 0; r < getLocalNumRows (); ++r) {
          const size_t nE = getNumEntriesInLocalRow(r);
          GlobalOrdinal gid = getRowMap()->getGlobalElement(r);
          out << std::setw(width) << myRank
              << std::setw(width) << gid
              << std::setw(width) << nE;
          if (vl == VERB_EXTREME) {
            if (isGloballyIndexed()) {
              global_inds_host_view_type rowinds;
              values_host_view_type rowvals;
              getGlobalRowView (gid, rowinds, rowvals);
              for (size_t j = 0; j < nE; ++j) {
                out << " (" << rowinds[j]
                    << ", " << rowvals[j]
                    << ") ";
              }
            }
            else if (isLocallyIndexed()) {
              local_inds_host_view_type rowinds;
              values_host_view_type rowvals;
              getLocalRowView (r, rowinds, rowvals);
              for (size_t j=0; j < nE; ++j) {
                out << " (" << getColMap()->getGlobalElement(rowinds[j])
                    << ", " << rowvals[j]
                    << ") ";
              }
            } // globally or locally indexed
          } // vl == VERB_EXTREME
          out << endl;
        } // for each row r on this process
      } // if (myRank == curRank)

      // Give output time to complete
      comm->barrier ();
      comm->barrier ();
      comm->barrier ();
    } // for each process p
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  checkSizes (const SrcDistObject& source)
  {
    // It's not clear what kind of compatibility checks on sizes can
    // be performed here.  Epetra_CrsGraph doesn't check any sizes for
    // compatibility.

    // Currently, the source object must be a RowMatrix with the same
    // four template parameters as the target CrsMatrix.  We might
    // relax this requirement later.
    const row_matrix_type* srcRowMat =
      dynamic_cast<const row_matrix_type*> (&source);
    return (srcRowMat != nullptr);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  applyCrsPadding(
    const typename crs_graph_type::padding_type& padding,
    const bool verbose)
  {
    using Details::ProfilingRegion;
    using Details::padCrsArrays;
    using std::endl;
    using LO = local_ordinal_type;
    using row_ptrs_type =
      typename local_graph_device_type::row_map_type::non_const_type;
    using range_policy =
      Kokkos::RangePolicy<execution_space, Kokkos::IndexType<LO>>;
    const char tfecfFuncName[] = "applyCrsPadding";
    const char suffix[] =
      ".  Please report this bug to the Tpetra developers.";
    ProfilingRegion regionCAP("Tpetra::CrsMatrix::applyCrsPadding");

    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsMatrix", tfecfFuncName);
      std::ostringstream os;
      os << *prefix << "padding: ";
      padding.print(os);
      os << endl;
      std::cerr << os.str();
    }
    const int myRank = ! verbose ? -1 : [&] () {
      auto map = this->getMap();
      if (map.is_null()) {
        return -1;
      }
      auto comm = map->getComm();
      if (comm.is_null()) {
        return -1;
      }
      return comm->getRank();
    } ();

    // NOTE (mfh 29 Jan 2020) This allocates the values array.
    if (! myGraph_->indicesAreAllocated()) {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Call allocateIndices" << endl;
        std::cerr << os.str();
      }
      allocateValues(GlobalIndices, GraphNotYetAllocated, verbose);
    }

    // FIXME (mfh 10 Feb 2020) We shouldn't actually reallocate
    // row_ptrs_beg or allocate row_ptrs_end unless the allocation
    // size needs to increase.  That should be the job of
    // padCrsArrays.

    // Making copies here because rowPtrsUnpacked_ has a const type. Otherwise, we
    // would use it directly.

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Allocate row_ptrs_beg: "
         << myGraph_->getRowPtrsUnpackedHost().extent(0) << endl;
      std::cerr << os.str();
    }
    using Kokkos::view_alloc;
    using Kokkos::WithoutInitializing;
    row_ptrs_type row_ptr_beg(view_alloc("row_ptr_beg", WithoutInitializing),
                              myGraph_->rowPtrsUnpacked_dev_.extent(0));
    // DEEP_COPY REVIEW - DEVICE-TO-DEVICE
    Kokkos::deep_copy(execution_space(),row_ptr_beg, myGraph_->rowPtrsUnpacked_dev_);
    
    const size_t N = row_ptr_beg.extent(0) == 0 ? size_t(0) :
      size_t(row_ptr_beg.extent(0) - 1);
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Allocate row_ptrs_end: " << N << endl;
      std::cerr << os.str();
    }
    row_ptrs_type row_ptr_end(
      view_alloc("row_ptr_end", WithoutInitializing), N);

    row_ptrs_type num_row_entries_d;

    const bool refill_num_row_entries =
      myGraph_->k_numRowEntries_.extent(0) != 0;
    
    if (refill_num_row_entries) { // unpacked storage
      // We can't assume correct *this capture until C++17, and it's
      // likely more efficient just to capture what we need anyway.
      num_row_entries_d = create_mirror_view_and_copy(memory_space(),
                                                 myGraph_->k_numRowEntries_);
      Kokkos::parallel_for
        ("Fill end row pointers", range_policy(0, N),
         KOKKOS_LAMBDA (const size_t i) {
          row_ptr_end(i) = row_ptr_beg(i) + num_row_entries_d(i);
        });
    }
    else {
      // FIXME (mfh 04 Feb 2020) Fix padCrsArrays so that if packed
      // storage, we don't need row_ptr_end to be separate allocation;
      // could just have it alias row_ptr_beg+1.
      Kokkos::parallel_for
        ("Fill end row pointers", range_policy(0, N),
         KOKKOS_LAMBDA (const size_t i) {
          row_ptr_end(i) = row_ptr_beg(i+1);
        });
    }

    if (myGraph_->isGloballyIndexed()) {
      padCrsArrays(row_ptr_beg, row_ptr_end,
                   myGraph_->gblInds_wdv,
                   valuesUnpacked_wdv, padding, myRank, verbose);
      const auto newValuesLen = valuesUnpacked_wdv.extent(0);
      const auto newColIndsLen = myGraph_->gblInds_wdv.extent(0);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (newValuesLen != newColIndsLen, std::logic_error,
         ": After padding, valuesUnpacked_wdv.extent(0)=" << newValuesLen
         << " != myGraph_->gblInds_wdv.extent(0)=" << newColIndsLen
         << suffix);
    }
    else {
      padCrsArrays(row_ptr_beg, row_ptr_end,
                   myGraph_->lclIndsUnpacked_wdv,
                   valuesUnpacked_wdv, padding, myRank, verbose);
      const auto newValuesLen = valuesUnpacked_wdv.extent(0);
      const auto newColIndsLen = myGraph_->lclIndsUnpacked_wdv.extent(0);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (newValuesLen != newColIndsLen, std::logic_error,
         ": After padding, valuesUnpacked_wdv.extent(0)=" << newValuesLen
         << " != myGraph_->lclIndsUnpacked_wdv.extent(0)=" << newColIndsLen
         << suffix);
    }

    if (refill_num_row_entries) {
      Kokkos::parallel_for
        ("Fill num entries", range_policy(0, N),
         KOKKOS_LAMBDA (const size_t i) {
          num_row_entries_d(i) = row_ptr_end(i) - row_ptr_beg(i);
        });
      Kokkos::deep_copy(myGraph_->k_numRowEntries_, num_row_entries_d);
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Assign myGraph_->rowPtrsUnpacked_; "
         << "old size: " << myGraph_->rowPtrsUnpacked_host_.extent(0)
         << ", new size: " << row_ptr_beg.extent(0) << endl;
      std::cerr << os.str();
      TEUCHOS_ASSERT( myGraph_->getRowPtrsUnpackedHost().extent(0) ==
                      row_ptr_beg.extent(0) );
    }
    myGraph_->setRowPtrsUnpacked(row_ptr_beg);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  copyAndPermuteStaticGraph(
    const RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& srcMat,
    const size_t numSameIDs,
    const LocalOrdinal permuteToLIDs[],
    const LocalOrdinal permuteFromLIDs[],
    const size_t numPermutes)
  {
    using Details::ProfilingRegion;
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using std::endl;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    const char tfecfFuncName[] = "copyAndPermuteStaticGraph";
    const char suffix[] =
      "  Please report this bug to the Tpetra developers.";
    ProfilingRegion regionCAP
      ("Tpetra::CrsMatrix::copyAndPermuteStaticGraph");

    const bool debug = Details::Behavior::debug("CrsGraph");
    const bool verbose = Details::Behavior::verbose("CrsGraph");
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsGraph", tfecfFuncName);
      std::ostringstream os;
      os << *prefix << "Start" << endl;
    }
    const char* const prefix_raw =
      verbose ? prefix.get()->c_str() : nullptr;

    const bool sourceIsLocallyIndexed = srcMat.isLocallyIndexed ();
    //
    // Copy the first numSame row from source to target (this matrix).
    // This involves copying rows corresponding to LIDs [0, numSame-1].
    //
    const map_type& srcRowMap = * (srcMat.getRowMap ());
    nonconst_global_inds_host_view_type rowInds;
    nonconst_values_host_view_type rowVals;
    const LO numSameIDs_as_LID = static_cast<LO> (numSameIDs);
    for (LO sourceLID = 0; sourceLID < numSameIDs_as_LID; ++sourceLID) {
      // Global ID for the current row index in the source matrix.
      // The first numSameIDs GIDs in the two input lists are the
      // same, so sourceGID == targetGID in this case.
      const GO sourceGID = srcRowMap.getGlobalElement (sourceLID);
      const GO targetGID = sourceGID;

      ArrayView<const GO>rowIndsConstView;
      ArrayView<const Scalar> rowValsConstView;

      if (sourceIsLocallyIndexed) {
        const size_t rowLength = srcMat.getNumEntriesInGlobalRow (sourceGID);
        if (rowLength > static_cast<size_t> (rowInds.size())) {
          Kokkos::resize(rowInds,rowLength);
          Kokkos::resize(rowVals,rowLength);
        }
        // Resizing invalidates an Array's views, so we must make new
        // ones, even if rowLength hasn't changed.
        nonconst_global_inds_host_view_type rowIndsView = Kokkos::subview(rowInds,std::make_pair((size_t)0, rowLength));
        nonconst_values_host_view_type rowValsView = Kokkos::subview(rowVals,std::make_pair((size_t)0, rowLength));

        // The source matrix is locally indexed, so we have to get a
        // copy.  Really it's the GIDs that have to be copied (because
        // they have to be converted from LIDs).
        size_t checkRowLength = 0;
        srcMat.getGlobalRowCopy (sourceGID, rowIndsView,
                                 rowValsView, checkRowLength);
        if (debug) {
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (rowLength != checkRowLength, std::logic_error, "For "
             "global row index " << sourceGID << ", the source "
             "matrix's getNumEntriesInGlobalRow returns a row length "
             "of " << rowLength << ", but getGlobalRowCopy reports "
             "a row length of " << checkRowLength << "." << suffix);
        }

        // KDDKDD UVM TEMPORARY:  refactor combineGlobalValues to take
        // KDDKDD UVM TEMPORARY:  Kokkos::View instead of ArrayView
        // KDDKDD UVM TEMPORARY:  For now, wrap the view in ArrayViews
        // KDDKDD UVM TEMPORARY:  Should be safe because we hold the KokkosViews
        rowIndsConstView = Teuchos::ArrayView<const GO> (  // BAD BAD BAD
                           rowIndsView.data(), rowIndsView.extent(0),
                           Teuchos::RCP_DISABLE_NODE_LOOKUP);
        rowValsConstView = Teuchos::ArrayView<const Scalar> (  // BAD BAD BAD
                           reinterpret_cast<const Scalar*>(rowValsView.data()), rowValsView.extent(0),
                           Teuchos::RCP_DISABLE_NODE_LOOKUP);
        // KDDKDD UVM TEMPORARY:  Add replace, sum, transform methods with
        // KDDKDD UVM TEMPORARY:  KokkosView interface
      }
      else { // source matrix is globally indexed.
        global_inds_host_view_type rowIndsView;
        values_host_view_type rowValsView;
        srcMat.getGlobalRowView(sourceGID, rowIndsView, rowValsView);
        // KDDKDD UVM TEMPORARY:  refactor combineGlobalValues to take
        // KDDKDD UVM TEMPORARY:  Kokkos::View instead of ArrayView
        // KDDKDD UVM TEMPORARY:  For now, wrap the view in ArrayViews
        // KDDKDD UVM TEMPORARY:  Should be safe because we hold the KokkosViews
        rowIndsConstView = Teuchos::ArrayView<const GO> (  // BAD BAD BAD
                           rowIndsView.data(), rowIndsView.extent(0),
                           Teuchos::RCP_DISABLE_NODE_LOOKUP);
        rowValsConstView = Teuchos::ArrayView<const Scalar> (  // BAD BAD BAD
                           reinterpret_cast<const Scalar*>(rowValsView.data()), rowValsView.extent(0),
                           Teuchos::RCP_DISABLE_NODE_LOOKUP);
        // KDDKDD UVM TEMPORARY:  Add replace, sum, transform methods with
        // KDDKDD UVM TEMPORARY:  KokkosView interface

      }

      // Applying a permutation to a matrix with a static graph
      // means REPLACE-ing entries.
      combineGlobalValues(targetGID, rowIndsConstView,
                          rowValsConstView, REPLACE,
                          prefix_raw, debug, verbose);
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Do permutes" << endl;
    }

    const map_type& tgtRowMap = * (this->getRowMap ());
    for (size_t p = 0; p < numPermutes; ++p) {
      const GO sourceGID = srcRowMap.getGlobalElement (permuteFromLIDs[p]);
      const GO targetGID = tgtRowMap.getGlobalElement (permuteToLIDs[p]);

      ArrayView<const GO> rowIndsConstView;
      ArrayView<const Scalar> rowValsConstView;

      if (sourceIsLocallyIndexed) {
        const size_t rowLength = srcMat.getNumEntriesInGlobalRow (sourceGID);
        if (rowLength > static_cast<size_t> (rowInds.size ())) {
          Kokkos::resize(rowInds,rowLength);
          Kokkos::resize(rowVals,rowLength);
        }
        // Resizing invalidates an Array's views, so we must make new
        // ones, even if rowLength hasn't changed.
        nonconst_global_inds_host_view_type rowIndsView = Kokkos::subview(rowInds,std::make_pair((size_t)0, rowLength));
        nonconst_values_host_view_type rowValsView = Kokkos::subview(rowVals,std::make_pair((size_t)0, rowLength));

        // The source matrix is locally indexed, so we have to get a
        // copy.  Really it's the GIDs that have to be copied (because
        // they have to be converted from LIDs).
        size_t checkRowLength = 0;
        srcMat.getGlobalRowCopy(sourceGID, rowIndsView,
                                rowValsView, checkRowLength);
        if (debug) {
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (rowLength != checkRowLength, std::logic_error, "For "
             "source matrix global row index " << sourceGID << ", "
             "getNumEntriesInGlobalRow returns a row length of " <<
             rowLength << ", but getGlobalRowCopy a row length of "
             << checkRowLength << "." << suffix);
        }

        // KDDKDD UVM TEMPORARY:  refactor combineGlobalValues to take
        // KDDKDD UVM TEMPORARY:  Kokkos::View instead of ArrayView
        // KDDKDD UVM TEMPORARY:  For now, wrap the view in ArrayViews
        // KDDKDD UVM TEMPORARY:  Should be safe because we hold the KokkosViews
        rowIndsConstView = Teuchos::ArrayView<const GO> (  // BAD BAD BAD
                           rowIndsView.data(), rowIndsView.extent(0),
                           Teuchos::RCP_DISABLE_NODE_LOOKUP);
        rowValsConstView = Teuchos::ArrayView<const Scalar> (  // BAD BAD BAD
                           reinterpret_cast<const Scalar*>(rowValsView.data()), rowValsView.extent(0),
                           Teuchos::RCP_DISABLE_NODE_LOOKUP);
        // KDDKDD UVM TEMPORARY:  Add replace, sum, transform methods with
        // KDDKDD UVM TEMPORARY:  KokkosView interface
      }
      else {
        global_inds_host_view_type rowIndsView;
        values_host_view_type rowValsView;
        srcMat.getGlobalRowView(sourceGID, rowIndsView, rowValsView);
        // KDDKDD UVM TEMPORARY:  refactor combineGlobalValues to take
        // KDDKDD UVM TEMPORARY:  Kokkos::View instead of ArrayView
        // KDDKDD UVM TEMPORARY:  For now, wrap the view in ArrayViews
        // KDDKDD UVM TEMPORARY:  Should be safe because we hold the KokkosViews
        rowIndsConstView = Teuchos::ArrayView<const GO> (  // BAD BAD BAD
                           rowIndsView.data(), rowIndsView.extent(0),
                           Teuchos::RCP_DISABLE_NODE_LOOKUP);
        rowValsConstView = Teuchos::ArrayView<const Scalar> (  // BAD BAD BAD
                           reinterpret_cast<const Scalar*>(rowValsView.data()), rowValsView.extent(0),
                           Teuchos::RCP_DISABLE_NODE_LOOKUP);
        // KDDKDD UVM TEMPORARY:  Add replace, sum, transform methods with
        // KDDKDD UVM TEMPORARY:  KokkosView interface
      }

      combineGlobalValues(targetGID, rowIndsConstView,
                          rowValsConstView, REPLACE,
                          prefix_raw, debug, verbose);
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  copyAndPermuteNonStaticGraph(
    const RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& srcMat,
    const size_t numSameIDs,
    const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& permuteToLIDs_dv,
    const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& permuteFromLIDs_dv,
    const size_t numPermutes)
  {
    using Details::ProfilingRegion;
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using std::endl;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    const char tfecfFuncName[] = "copyAndPermuteNonStaticGraph";
    const char suffix[] =
      "  Please report this bug to the Tpetra developers.";
    ProfilingRegion regionCAP
      ("Tpetra::CrsMatrix::copyAndPermuteNonStaticGraph");

    const bool debug = Details::Behavior::debug("CrsGraph");
    const bool verbose = Details::Behavior::verbose("CrsGraph");
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsGraph", tfecfFuncName);
      std::ostringstream os;
      os << *prefix << "Start" << endl;
    }
    const char* const prefix_raw =
      verbose ? prefix.get()->c_str() : nullptr;

    {
      using row_graph_type = RowGraph<LO, GO, Node>;
      const row_graph_type& srcGraph = *(srcMat.getGraph());
      auto padding =
        myGraph_->computeCrsPadding(srcGraph, numSameIDs,
          permuteToLIDs_dv, permuteFromLIDs_dv, verbose);
      applyCrsPadding(*padding, verbose);
    }
    const bool sourceIsLocallyIndexed = srcMat.isLocallyIndexed ();
    //
    // Copy the first numSame row from source to target (this matrix).
    // This involves copying rows corresponding to LIDs [0, numSame-1].
    //
    const map_type& srcRowMap = * (srcMat.getRowMap ());
    const LO numSameIDs_as_LID = static_cast<LO> (numSameIDs);
    using gids_type = nonconst_global_inds_host_view_type;
    using vals_type = nonconst_values_host_view_type;
    gids_type rowInds;
    vals_type rowVals;
    for (LO sourceLID = 0; sourceLID < numSameIDs_as_LID; ++sourceLID) {
      // Global ID for the current row index in the source matrix.
      // The first numSameIDs GIDs in the two input lists are the
      // same, so sourceGID == targetGID in this case.
      const GO sourceGID = srcRowMap.getGlobalElement (sourceLID);
      const GO targetGID = sourceGID;

      ArrayView<const GO> rowIndsConstView;
      ArrayView<const Scalar> rowValsConstView;

      if (sourceIsLocallyIndexed) {

        const size_t rowLength = srcMat.getNumEntriesInGlobalRow (sourceGID);
        if (rowLength > static_cast<size_t> (rowInds.extent(0))) {
          Kokkos::resize(rowInds,rowLength);
          Kokkos::resize(rowVals,rowLength);
        }
        // Resizing invalidates an Array's views, so we must make new
        // ones, even if rowLength hasn't changed.
        gids_type rowIndsView = Kokkos::subview(rowInds,std::make_pair((size_t)0, rowLength));
        vals_type rowValsView = Kokkos::subview(rowVals,std::make_pair((size_t)0, rowLength));

        // The source matrix is locally indexed, so we have to get a
        // copy.  Really it's the GIDs that have to be copied (because
        // they have to be converted from LIDs).
        size_t checkRowLength = 0;
        srcMat.getGlobalRowCopy (sourceGID, rowIndsView, rowValsView,
                                 checkRowLength);
        if (debug) {
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (rowLength != checkRowLength, std::logic_error, ": For "
             "global row index " << sourceGID << ", the source "
             "matrix's getNumEntriesInGlobalRow returns a row length "
             "of " << rowLength << ", but getGlobalRowCopy reports "
             "a row length of " << checkRowLength << "." << suffix);
        }
        rowIndsConstView = Teuchos::ArrayView<const GO>(rowIndsView.data(), rowLength);
        rowValsConstView = Teuchos::ArrayView<const Scalar>(reinterpret_cast<Scalar *>(rowValsView.data()), rowLength);
      }
      else { // source matrix is globally indexed.
        global_inds_host_view_type rowIndsView;
        values_host_view_type rowValsView;
        srcMat.getGlobalRowView(sourceGID, rowIndsView, rowValsView);

        // KDDKDD UVM TEMPORARY:  refactor combineGlobalValues to take
        // KDDKDD UVM TEMPORARY:  Kokkos::View instead of ArrayView
        // KDDKDD UVM TEMPORARY:  For now, wrap the view in ArrayViews
        // KDDKDD UVM TEMPORARY:  Should be safe because we hold the KokkosViews
        rowIndsConstView = Teuchos::ArrayView<const GO> (  // BAD BAD BAD
                           rowIndsView.data(), rowIndsView.extent(0),
                           Teuchos::RCP_DISABLE_NODE_LOOKUP);
        rowValsConstView = Teuchos::ArrayView<const Scalar> (  // BAD BAD BAD
                           reinterpret_cast<const Scalar*>(rowValsView.data()), rowValsView.extent(0),
                           Teuchos::RCP_DISABLE_NODE_LOOKUP);
        // KDDKDD UVM TEMPORARY:  Add replace, sum, transform methods with
        // KDDKDD UVM TEMPORARY:  KokkosView interface
      }

      // Combine the data into the target matrix.
      insertGlobalValuesFilteredChecked(targetGID, rowIndsConstView,
        rowValsConstView, prefix_raw, debug, verbose);
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Do permutes" << endl;
    }
    const LO* const permuteFromLIDs = permuteFromLIDs_dv.view_host().data();
    const LO* const permuteToLIDs = permuteToLIDs_dv.view_host().data();

    const map_type& tgtRowMap = * (this->getRowMap ());
    for (size_t p = 0; p < numPermutes; ++p) {
      const GO sourceGID = srcRowMap.getGlobalElement (permuteFromLIDs[p]);
      const GO targetGID = tgtRowMap.getGlobalElement (permuteToLIDs[p]);

      ArrayView<const GO> rowIndsConstView;
      ArrayView<const Scalar> rowValsConstView;

      if (sourceIsLocallyIndexed) {
        const size_t rowLength = srcMat.getNumEntriesInGlobalRow (sourceGID);
        if (rowLength > static_cast<size_t> (rowInds.extent(0))) {
          Kokkos::resize(rowInds,rowLength);
          Kokkos::resize(rowVals,rowLength);
        }
        // Resizing invalidates an Array's views, so we must make new
        // ones, even if rowLength hasn't changed.
        gids_type rowIndsView = Kokkos::subview(rowInds,std::make_pair((size_t)0, rowLength));
        vals_type rowValsView = Kokkos::subview(rowVals,std::make_pair((size_t)0, rowLength));

        // The source matrix is locally indexed, so we have to get a
        // copy.  Really it's the GIDs that have to be copied (because
        // they have to be converted from LIDs).
        size_t checkRowLength = 0;
        srcMat.getGlobalRowCopy(sourceGID, rowIndsView,
                                rowValsView, checkRowLength);
        if (debug) {
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (rowLength != checkRowLength, std::logic_error, "For "
             "source matrix global row index " << sourceGID << ", "
             "getNumEntriesInGlobalRow returns a row length of " <<
             rowLength << ", but getGlobalRowCopy a row length of "
             << checkRowLength << "." << suffix);
        }
        rowIndsConstView = Teuchos::ArrayView<const GO>(rowIndsView.data(), rowLength);
        rowValsConstView = Teuchos::ArrayView<const Scalar>(reinterpret_cast<Scalar *>(rowValsView.data()), rowLength);
      }
      else {
        global_inds_host_view_type rowIndsView;
        values_host_view_type rowValsView;
        srcMat.getGlobalRowView(sourceGID, rowIndsView, rowValsView);

        // KDDKDD UVM TEMPORARY:  refactor combineGlobalValues to take
        // KDDKDD UVM TEMPORARY:  Kokkos::View instead of ArrayView
        // KDDKDD UVM TEMPORARY:  For now, wrap the view in ArrayViews
        // KDDKDD UVM TEMPORARY:  Should be safe because we hold the KokkosViews
        rowIndsConstView = Teuchos::ArrayView<const GO> (  // BAD BAD BAD
                           rowIndsView.data(), rowIndsView.extent(0),
                           Teuchos::RCP_DISABLE_NODE_LOOKUP);
        rowValsConstView = Teuchos::ArrayView<const Scalar> (  // BAD BAD BAD
                           reinterpret_cast<const Scalar*>(rowValsView.data()), rowValsView.extent(0),
                           Teuchos::RCP_DISABLE_NODE_LOOKUP);
        // KDDKDD UVM TEMPORARY:  Add replace, sum, transform methods with
        // KDDKDD UVM TEMPORARY:  KokkosView interface
      }

      // Combine the data into the target matrix.
      insertGlobalValuesFilteredChecked(targetGID, rowIndsConstView,
        rowValsConstView, prefix_raw, debug, verbose);
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  copyAndPermute(
    const SrcDistObject& srcObj,
    const size_t numSameIDs,
    const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& permuteToLIDs,
    const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& permuteFromLIDs,
    const CombineMode /*CM*/)
  {
    using Details::Behavior;
    using Details::dualViewStatusToString;
    using Details::ProfilingRegion;
    using std::endl;

    // Method name string for TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC.
    const char tfecfFuncName[] = "copyAndPermute: ";
    ProfilingRegion regionCAP("Tpetra::CrsMatrix::copyAndPermute");

    const bool verbose = Behavior::verbose("CrsMatrix");
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsMatrix", "copyAndPermute");
      std::ostringstream os;
      os << *prefix << endl
         << *prefix << "  numSameIDs: " << numSameIDs << endl
         << *prefix << "  numPermute: " << permuteToLIDs.extent(0)
         << endl
         << *prefix << "  "
         << dualViewStatusToString (permuteToLIDs, "permuteToLIDs")
         << endl
         << *prefix << "  "
         << dualViewStatusToString (permuteFromLIDs, "permuteFromLIDs")
         << endl
         << *prefix << "  "
         << "isStaticGraph: " << (isStaticGraph() ? "true" : "false")
         << endl;
      std::cerr << os.str ();
    }

    const auto numPermute = permuteToLIDs.extent (0);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (numPermute != permuteFromLIDs.extent (0),
       std::invalid_argument, "permuteToLIDs.extent(0) = "
       << numPermute << "!= permuteFromLIDs.extent(0) = "
       << permuteFromLIDs.extent (0) << ".");

    // This dynamic cast should succeed, because we've already tested
    // it in checkSizes().
    using RMT = RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    const RMT& srcMat = dynamic_cast<const RMT&> (srcObj);
    if (isStaticGraph ()) {
      TEUCHOS_ASSERT( ! permuteToLIDs.need_sync_host () );
      auto permuteToLIDs_h = permuteToLIDs.view_host ();
      TEUCHOS_ASSERT( ! permuteFromLIDs.need_sync_host () );
      auto permuteFromLIDs_h = permuteFromLIDs.view_host ();

      copyAndPermuteStaticGraph(srcMat, numSameIDs,
                                permuteToLIDs_h.data(),
                                permuteFromLIDs_h.data(),
                                numPermute);
    }
    else {
      copyAndPermuteNonStaticGraph(srcMat, numSameIDs, permuteToLIDs,
                                   permuteFromLIDs, numPermute);
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str();
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  packAndPrepare
  (const SrcDistObject& source,
   const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& exportLIDs,
   Kokkos::DualView<char*, buffer_device_type>& exports,
   Kokkos::DualView<size_t*, buffer_device_type> numPacketsPerLID,
   size_t& constantNumPackets)
  {
    using Details::Behavior;
    using Details::dualViewStatusToString;
    using Details::ProfilingRegion;
    using Teuchos::outArg;
    using Teuchos::REDUCE_MAX;
    using Teuchos::reduceAll;
    using std::endl;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "packAndPrepare: ";
    ProfilingRegion regionPAP ("Tpetra::CrsMatrix::packAndPrepare");

    const bool debug = Behavior::debug("CrsMatrix");
    const bool verbose = Behavior::verbose("CrsMatrix");

    // Processes on which the communicator is null should not participate.
    Teuchos::RCP<const Teuchos::Comm<int> > pComm = this->getComm ();
    if (pComm.is_null ()) {
      return;
    }
    const Teuchos::Comm<int>& comm = *pComm;
    const int myRank = comm.getSize ();

    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsMatrix", "packAndPrepare");
      std::ostringstream os;
      os << *prefix << "Start" << endl
         << *prefix << "  "
         << dualViewStatusToString (exportLIDs, "exportLIDs")
         << endl
         << *prefix << "  "
         << dualViewStatusToString (exports, "exports")
         << endl
         << *prefix << "  "
         << dualViewStatusToString (numPacketsPerLID, "numPacketsPerLID")
         << endl;
      std::cerr << os.str ();
    }

    // Attempt to cast the source object to CrsMatrix.  If successful,
    // use the source object's packNew() method to pack its data for
    // communication.  Otherwise, attempt to cast to RowMatrix; if
    // successful, use the source object's pack() method.  Otherwise,
    // the source object doesn't have the right type.
    //
    // FIXME (mfh 30 Jun 2013, 11 Sep 2017) We don't even need the
    // RowMatrix to have the same Node type.  Unfortunately, we don't
    // have a way to ask if the RowMatrix is "a RowMatrix with any
    // Node type," since RowMatrix doesn't have a base class.  A
    // hypothetical RowMatrixBase<Scalar, LO, GO> class, which does
    // not currently exist, would satisfy this requirement.
    //
    // Why RowMatrixBase<Scalar, LO, GO>?  The source object's Scalar
    // type doesn't technically need to match the target object's
    // Scalar type, so we could just have RowMatrixBase<LO, GO>.  LO
    // and GO need not be the same, as long as there is no overflow of
    // the indices.  However, checking for index overflow is global
    // and therefore undesirable.

    std::ostringstream msg; // for collecting error messages
    int lclBad = 0; // to be set below

    using crs_matrix_type = CrsMatrix<Scalar, LO, GO, Node>;
    const crs_matrix_type* srcCrsMat =
      dynamic_cast<const crs_matrix_type*> (&source);
    if (srcCrsMat != nullptr) {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Source matrix same (CrsMatrix) type as target; "
          "calling packNew" << endl;
        std::cerr << os.str ();
      }
      try {
        srcCrsMat->packNew (exportLIDs, exports, numPacketsPerLID,
                            constantNumPackets);
      }
      catch (std::exception& e) {
        lclBad = 1;
        msg << "Proc " << myRank << ": " << e.what () << std::endl;
      }
    }
    else {
      using Kokkos::HostSpace;
      using Kokkos::subview;
      using exports_type = Kokkos::DualView<char*, buffer_device_type>;
      using range_type = Kokkos::pair<size_t, size_t>;

      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Source matrix NOT same (CrsMatrix) type as target"
           << endl;
        std::cerr << os.str ();
      }

      const row_matrix_type* srcRowMat =
        dynamic_cast<const row_matrix_type*> (&source);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (srcRowMat == nullptr, std::invalid_argument,
         "The source object of the Import or Export operation is neither a "
         "CrsMatrix (with the same template parameters as the target object), "
         "nor a RowMatrix (with the same first four template parameters as the "
         "target object).");

      // For the RowMatrix case, we need to convert from
      // Kokkos::DualView to Teuchos::Array*.  This doesn't need to be
      // so terribly efficient, since packing a non-CrsMatrix
      // RowMatrix for Import/Export into a CrsMatrix is not a
      // critical case.  Thus, we may allocate Teuchos::Array objects
      // here and copy to and from Kokkos::*View.

      // View exportLIDs's host data as a Teuchos::ArrayView.
      TEUCHOS_ASSERT( ! exportLIDs.need_sync_host () );
      auto exportLIDs_h = exportLIDs.view_host ();
      Teuchos::ArrayView<const LO> exportLIDs_av (exportLIDs_h.data (),
                                                  exportLIDs_h.size ());

      // pack() will allocate exports_a as needed.  We'll copy back
      // into exports (after (re)allocating exports if needed) below.
      Teuchos::Array<char> exports_a;

      // View exportLIDs' host data as a Teuchos::ArrayView.  We don't
      // need to sync, since we're doing write-only access, but we do
      // need to mark the DualView as modified on host.

      numPacketsPerLID.clear_sync_state (); // write-only access
      numPacketsPerLID.modify_host ();
      auto numPacketsPerLID_h = numPacketsPerLID.view_host ();
      Teuchos::ArrayView<size_t> numPacketsPerLID_av (numPacketsPerLID_h.data (),
                                                      numPacketsPerLID_h.size ());

      // Invoke RowMatrix's legacy pack() interface, using above
      // Teuchos::Array* objects.
      try {
        srcRowMat->pack (exportLIDs_av, exports_a, numPacketsPerLID_av,
                         constantNumPackets);
      }
      catch (std::exception& e) {
        lclBad = 1;
        msg << "Proc " << myRank << ": " << e.what () << std::endl;
      }

      // Allocate 'exports', and copy exports_a back into it.
      const size_t newAllocSize = static_cast<size_t> (exports_a.size ());
      if (static_cast<size_t> (exports.extent (0)) < newAllocSize) {
        const std::string oldLabel = exports.d_view.label ();
        const std::string newLabel = (oldLabel == "") ? "exports" : oldLabel;
        exports = exports_type (newLabel, newAllocSize);
      }
      // It's safe to assume that we're working on host anyway, so
      // just keep exports sync'd to host.
      // ignore current device contents
      exports.modify_host();

      auto exports_h = exports.view_host ();
      auto exports_h_sub = subview (exports_h, range_type (0, newAllocSize));

      // Kokkos::deep_copy needs a Kokkos::View input, so turn
      // exports_a into a nonowning Kokkos::View first before copying.
      typedef typename exports_type::t_host::execution_space HES;
      typedef Kokkos::Device<HES, HostSpace> host_device_type;
      Kokkos::View<const char*, host_device_type>
        exports_a_kv (exports_a.getRawPtr (), newAllocSize);
      // DEEP_COPY REVIEW - NOT TESTED
      Kokkos::deep_copy (exports_h_sub, exports_a_kv);
    }

    if (debug) {
      int gblBad = 0; // output argument; to be set below
      reduceAll<int, int> (comm, REDUCE_MAX, lclBad, outArg (gblBad));
      if (gblBad != 0) {
        Tpetra::Details::gathervPrint (std::cerr, msg.str (), comm);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (true, std::logic_error, "packNew() or pack() threw an exception on "
           "one or more participating processes.");
      }
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (lclBad != 0, std::logic_error, "packNew threw an exception on one "
         "or more participating processes.  Here is this process' error "
         "message: " << msg.str ());
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "packAndPrepare: Done!" << endl
         << *prefix << "  "
         << dualViewStatusToString (exportLIDs, "exportLIDs")
         << endl
         << *prefix << "  "
         << dualViewStatusToString (exports, "exports")
         << endl
         << *prefix << "  "
         << dualViewStatusToString (numPacketsPerLID, "numPacketsPerLID")
         << endl;
      std::cerr << os.str ();
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  packRow (char exports[],
           const size_t offset,
           const size_t numEnt,
           const GlobalOrdinal gidsIn[],
           const impl_scalar_type valsIn[],
           const size_t numBytesPerValue) const
  {
    using Kokkos::View;
    using Kokkos::subview;
    using Tpetra::Details::PackTraits;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef impl_scalar_type ST;

    if (numEnt == 0) {
      // Empty rows always take zero bytes, to ensure sparsity.
      return 0;
    }

    const GO gid = 0; // packValueCount wants this
    const LO numEntLO = static_cast<size_t> (numEnt);

    const size_t numEntBeg = offset;
    const size_t numEntLen = PackTraits<LO>::packValueCount (numEntLO);
    const size_t gidsBeg = numEntBeg + numEntLen;
    const size_t gidsLen = numEnt * PackTraits<GO>::packValueCount (gid);
    const size_t valsBeg = gidsBeg + gidsLen;
    const size_t valsLen = numEnt * numBytesPerValue;

    char* const numEntOut = exports + numEntBeg;
    char* const gidsOut = exports + gidsBeg;
    char* const valsOut = exports + valsBeg;

    size_t numBytesOut = 0;
    int errorCode = 0;
    numBytesOut += PackTraits<LO>::packValue (numEntOut, numEntLO);

    {
      Kokkos::pair<int, size_t> p;
      p = PackTraits<GO>::packArray (gidsOut, gidsIn, numEnt);
      errorCode += p.first;
      numBytesOut += p.second;

      p = PackTraits<ST>::packArray (valsOut, valsIn, numEnt);
      errorCode += p.first;
      numBytesOut += p.second;
    }

    const size_t expectedNumBytes = numEntLen + gidsLen + valsLen;
    TEUCHOS_TEST_FOR_EXCEPTION
      (numBytesOut != expectedNumBytes, std::logic_error, "packRow: "
       "numBytesOut = " << numBytesOut << " != expectedNumBytes = "
       << expectedNumBytes << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (errorCode != 0, std::runtime_error, "packRow: "
       "PackTraits::packArray returned a nonzero error code");

    return numBytesOut;
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  unpackRow (GlobalOrdinal gidsOut[],
             impl_scalar_type valsOut[],
             const char imports[],
             const size_t offset,
             const size_t numBytes,
             const size_t numEnt,
             const size_t numBytesPerValue)
  {
    using Kokkos::View;
    using Kokkos::subview;
    using Tpetra::Details::PackTraits;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef impl_scalar_type ST;

    Details::ProfilingRegion region_upack_row(
      "Tpetra::CrsMatrix::unpackRow",
      "Import/Export"
    );

    if (numBytes == 0) {
      // Rows with zero bytes should always have zero entries.
      if (numEnt != 0) {
        const int myRank = this->getMap ()->getComm ()->getRank ();
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::logic_error, "(Proc " << myRank << ") CrsMatrix::"
           "unpackRow: The number of bytes to unpack numBytes=0, but the "
           "number of entries to unpack (as reported by numPacketsPerLID) "
           "for this row numEnt=" << numEnt << " != 0.");
      }
      return 0;
    }

    if (numEnt == 0 && numBytes != 0) {
      const int myRank = this->getMap ()->getComm ()->getRank ();
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "(Proc " << myRank << ") CrsMatrix::"
         "unpackRow: The number of entries to unpack (as reported by "
         "numPacketsPerLID) numEnt=0, but the number of bytes to unpack "
         "numBytes=" << numBytes << " != 0.");
    }

    const GO gid = 0; // packValueCount wants this
    const LO lid = 0; // packValueCount wants this

    const size_t numEntBeg = offset;
    const size_t numEntLen = PackTraits<LO>::packValueCount (lid);
    const size_t gidsBeg = numEntBeg + numEntLen;
    const size_t gidsLen = numEnt * PackTraits<GO>::packValueCount (gid);
    const size_t valsBeg = gidsBeg + gidsLen;
    const size_t valsLen = numEnt * numBytesPerValue;

    const char* const numEntIn = imports + numEntBeg;
    const char* const gidsIn = imports + gidsBeg;
    const char* const valsIn = imports + valsBeg;

    size_t numBytesOut = 0;
    int errorCode = 0;
    LO numEntOut;
    numBytesOut += PackTraits<LO>::unpackValue (numEntOut, numEntIn);
    if (static_cast<size_t> (numEntOut) != numEnt ||
        numEntOut == static_cast<LO> (0)) {
      const int myRank = this->getMap ()->getComm ()->getRank ();
      std::ostringstream os;
      os << "(Proc " << myRank << ") CrsMatrix::unpackRow: ";
      bool firstErrorCondition = false;
      if (static_cast<size_t> (numEntOut) != numEnt) {
        os << "Number of entries from numPacketsPerLID numEnt=" << numEnt
           << " does not equal number of entries unpacked from imports "
          "buffer numEntOut=" << numEntOut << ".";
        firstErrorCondition = true;
      }
      if (numEntOut == static_cast<LO> (0)) {
        if (firstErrorCondition) {
          os << "  Also, ";
        }
        os << "Number of entries unpacked from imports buffer numEntOut=0, "
          "but number of bytes to unpack for this row numBytes=" << numBytes
           << " != 0.  This should never happen, since packRow should only "
          "ever pack rows with a nonzero number of entries.  In this case, "
          "the number of entries from numPacketsPerLID is numEnt=" << numEnt
           << ".";
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, os.str ());
    }

    {
      Kokkos::pair<int, size_t> p;
      p = PackTraits<GO>::unpackArray (gidsOut, gidsIn, numEnt);
      errorCode += p.first;
      numBytesOut += p.second;

      p = PackTraits<ST>::unpackArray (valsOut, valsIn, numEnt);
      errorCode += p.first;
      numBytesOut += p.second;
    }

    TEUCHOS_TEST_FOR_EXCEPTION
      (numBytesOut != numBytes, std::logic_error, "unpackRow: numBytesOut = "
       << numBytesOut << " != numBytes = " << numBytes << ".");

    const size_t expectedNumBytes = numEntLen + gidsLen + valsLen;
    TEUCHOS_TEST_FOR_EXCEPTION
      (numBytesOut != expectedNumBytes, std::logic_error, "unpackRow: "
       "numBytesOut = " << numBytesOut << " != expectedNumBytes = "
       << expectedNumBytes << ".");

    TEUCHOS_TEST_FOR_EXCEPTION
      (errorCode != 0, std::runtime_error, "unpackRow: "
       "PackTraits::unpackArray returned a nonzero error code");

    return numBytesOut;
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  allocatePackSpaceNew (Kokkos::DualView<char*, buffer_device_type>& exports,
                        size_t& totalNumEntries,
                        const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& exportLIDs) const
  {
    using Details::Behavior;
    using Details::dualViewStatusToString;
    using std::endl;
    typedef impl_scalar_type IST;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    //const char tfecfFuncName[] = "allocatePackSpaceNew: ";

    // mfh 18 Oct 2017: Set TPETRA_VERBOSE to true for copious debug
    // output to std::cerr on every MPI process.  This is unwise for
    // runs with large numbers of MPI processes.
    const bool verbose = Behavior::verbose("CrsMatrix");
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsMatrix", "allocatePackSpaceNew");
      std::ostringstream os;
      os << *prefix << "Before:"
         << endl
         << *prefix << "  "
         << dualViewStatusToString (exports, "exports")
         << endl
         << *prefix << "  "
         << dualViewStatusToString (exportLIDs, "exportLIDs")
         << endl;
      std::cerr << os.str ();
    }

    // The number of export LIDs must fit in LocalOrdinal, assuming
    // that the LIDs are distinct and valid on the calling process.
    const LO numExportLIDs = static_cast<LO> (exportLIDs.extent (0));

    TEUCHOS_ASSERT( ! exportLIDs.need_sync_host () );
    auto exportLIDs_h = exportLIDs.view_host ();

    // Count the total number of matrix entries to send.
    totalNumEntries = 0;
    for (LO i = 0; i < numExportLIDs; ++i) {
      const LO lclRow = exportLIDs_h[i];
      size_t curNumEntries = this->getNumEntriesInLocalRow (lclRow);
      // FIXME (mfh 25 Jan 2015) We should actually report invalid row
      // indices as an error.  Just consider them nonowned for now.
      if (curNumEntries == Teuchos::OrdinalTraits<size_t>::invalid ()) {
        curNumEntries = 0;
      }
      totalNumEntries += curNumEntries;
    }

    // FIXME (mfh 24 Feb 2013, 24 Mar 2017) This code is only correct
    // if sizeof(IST) is a meaningful representation of the amount of
    // data in a Scalar instance.  (LO and GO are always built-in
    // integer types.)
    //
    // Allocate the exports array.  It does NOT need padding for
    // alignment, since we use memcpy to write to / read from send /
    // receive buffers.
    const size_t allocSize =
      static_cast<size_t> (numExportLIDs) * sizeof (LO) +
      totalNumEntries * (sizeof (IST) + sizeof (GO));
    if (static_cast<size_t> (exports.extent (0)) < allocSize) {
      using exports_type = Kokkos::DualView<char*, buffer_device_type>;

      const std::string oldLabel = exports.d_view.label ();
      const std::string newLabel = (oldLabel == "") ? "exports" : oldLabel;
      exports = exports_type (newLabel, allocSize);
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "After:"
         << endl
         << *prefix << "  "
         << dualViewStatusToString (exports, "exports")
         << endl
         << *prefix << "  "
         << dualViewStatusToString (exportLIDs, "exportLIDs")
         << endl;
      std::cerr << os.str ();
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  packNew (const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& exportLIDs,
           Kokkos::DualView<char*, buffer_device_type>& exports,
           const Kokkos::DualView<size_t*, buffer_device_type>& numPacketsPerLID,
           size_t& constantNumPackets) const
  {
    // The call to packNew in packAndPrepare catches and handles any exceptions.
    Details::ProfilingRegion region_pack_new("Tpetra::CrsMatrix::packNew", "Import/Export");
    if (this->isStaticGraph ()) {
      using ::Tpetra::Details::packCrsMatrixNew;
      packCrsMatrixNew (*this, exports, numPacketsPerLID, exportLIDs,
                        constantNumPackets);
    }
    else {
      this->packNonStaticNew (exportLIDs, exports, numPacketsPerLID,
                              constantNumPackets);
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  packNonStaticNew (const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& exportLIDs,
                    Kokkos::DualView<char*, buffer_device_type>& exports,
                    const Kokkos::DualView<size_t*, buffer_device_type>& numPacketsPerLID,
                    size_t& constantNumPackets) const
  {
    using Details::Behavior;
    using Details::dualViewStatusToString;
    using Details::PackTraits;
    using Details::create_mirror_view_from_raw_host_array;
    using Kokkos::View;
    using std::endl;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    using ST = impl_scalar_type;
    const char tfecfFuncName[] = "packNonStaticNew: ";

    const bool verbose = Behavior::verbose("CrsMatrix");
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsMatrix", "packNonStaticNew");
      std::ostringstream os;
      os << *prefix << "Start" << endl;
      std::cerr << os.str ();
    }

    const size_t numExportLIDs = static_cast<size_t> (exportLIDs.extent (0));
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (numExportLIDs != static_cast<size_t> (numPacketsPerLID.extent (0)),
       std::invalid_argument, "exportLIDs.size() = " << numExportLIDs
       << " != numPacketsPerLID.size() = " << numPacketsPerLID.extent (0)
       << ".");

    // Setting this to zero tells the caller to expect a possibly
    // different ("nonconstant") number of packets per local index
    // (i.e., a possibly different number of entries per row).
    constantNumPackets = 0;

    // The pack buffer 'exports' enters this method possibly
    // unallocated.  Do the first two parts of "Count, allocate, fill,
    // compute."
    size_t totalNumEntries = 0;
    this->allocatePackSpaceNew (exports, totalNumEntries, exportLIDs);
    const size_t bufSize = static_cast<size_t> (exports.extent (0));

    // Write-only host access
    exports.clear_sync_state();
    exports.modify_host();
    auto exports_h = exports.view_host ();
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "After marking exports as modified on host, "
         << dualViewStatusToString (exports, "exports") << endl;
      std::cerr << os.str ();
    }

    // Read-only host access
    auto exportLIDs_h = exportLIDs.view_host ();

    // Write-only host access
    const_cast<Kokkos::DualView<size_t*, buffer_device_type>*>(&numPacketsPerLID)->clear_sync_state();
    const_cast<Kokkos::DualView<size_t*, buffer_device_type>*>(&numPacketsPerLID)->modify_host();
    auto numPacketsPerLID_h = numPacketsPerLID.view_host ();

    // Compute the number of "packets" (in this case, bytes) per
    // export LID (in this case, local index of the row to send), and
    // actually pack the data.
    auto maxRowNumEnt = this->getLocalMaxNumRowEntries();


    // Temporary buffer for global column indices.
    typename global_inds_host_view_type::non_const_type gidsIn_k;
    if (this->isLocallyIndexed()) { // Need storage for Global IDs
      gidsIn_k = 
        typename global_inds_host_view_type::non_const_type("packGids",
                                                            maxRowNumEnt);
    }

    size_t offset = 0; // current index into 'exports' array.
    for (size_t i = 0; i < numExportLIDs; ++i) {
      const LO lclRow = exportLIDs_h[i];

      size_t numBytes = 0;
      size_t numEnt = this->getNumEntriesInLocalRow (lclRow);

      // Only pack this row's data if it has a nonzero number of
      // entries.  We can do this because receiving processes get the
      // number of packets, and will know that zero packets means zero
      // entries.
      if (numEnt == 0) {
        numPacketsPerLID_h[i] = 0;
        continue;
      }

      if (this->isLocallyIndexed ()) {
        typename global_inds_host_view_type::non_const_type gidsIn; 
        values_host_view_type valsIn;
        // If the matrix is locally indexed on the calling process, we
        // have to use its column Map (which it _must_ have in this
        // case) to convert to global indices.
        local_inds_host_view_type lidsIn;
        this->getLocalRowView (lclRow, lidsIn, valsIn);
        const map_type& colMap = * (this->getColMap ());
        for (size_t k = 0; k < numEnt; ++k) {
          gidsIn_k[k] = colMap.getGlobalElement (lidsIn[k]);
        }
        gidsIn = Kokkos::subview(gidsIn_k, Kokkos::make_pair(GO(0),GO(numEnt)));

        const size_t numBytesPerValue =
          PackTraits<ST>::packValueCount (valsIn[0]);
        numBytes = this->packRow (exports_h.data (), offset, numEnt,
                                  gidsIn.data (), valsIn.data (), 
                                  numBytesPerValue);
      }
      else if (this->isGloballyIndexed ()) {
        global_inds_host_view_type gidsIn; 
        values_host_view_type valsIn;
        // If the matrix is globally indexed on the calling process,
        // then we can use the column indices directly.  However, we
        // have to get the global row index.  The calling process must
        // have a row Map, since otherwise it shouldn't be participating
        // in packing operations.
        const map_type& rowMap = * (this->getRowMap ());
        const GO gblRow = rowMap.getGlobalElement (lclRow);
        this->getGlobalRowView (gblRow, gidsIn, valsIn);

        const size_t numBytesPerValue =
          PackTraits<ST>::packValueCount (valsIn[0]);
        numBytes = this->packRow (exports_h.data (), offset, numEnt, 
                                  gidsIn.data (), valsIn.data (),
                                  numBytesPerValue);
      }
      // mfh 11 Sep 2017: Currently, if the matrix is neither globally
      // nor locally indexed, then it has no entries.  Therefore,
      // there is nothing to pack.  No worries!

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (offset > bufSize || offset + numBytes > bufSize, std::logic_error,
         "First invalid offset into 'exports' pack buffer at index i = " << i
         << ".  exportLIDs_h[i]: " << exportLIDs_h[i] << ", bufSize: " <<
         bufSize << ", offset: " << offset << ", numBytes: " << numBytes <<
         ".");
      // numPacketsPerLID_h[i] is the number of "packets" in the
      // current local row i.  Packet=char (really "byte") so use the
      // number of bytes of the packed data for that row.
      numPacketsPerLID_h[i] = numBytes;
      offset += numBytes;
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Tpetra::CrsMatrix::packNonStaticNew: After:" << endl
         << *prefix << "  "
         << dualViewStatusToString (exports, "exports")
         << endl
         << *prefix << "  "
         << dualViewStatusToString (exportLIDs, "exportLIDs")
         << endl;
      std::cerr << os.str ();
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  combineGlobalValuesRaw(const LocalOrdinal lclRow,
                         const LocalOrdinal numEnt,
                         const impl_scalar_type vals[],
                         const GlobalOrdinal cols[],
                         const Tpetra::CombineMode combMode,
                         const char* const prefix,
                         const bool debug,
                         const bool verbose)
  {
    using GO = GlobalOrdinal;

    // mfh 23 Mar 2017: This branch is not thread safe in a debug
    // build, due to use of Teuchos::ArrayView; see #229.
    const GO gblRow = myGraph_->rowMap_->getGlobalElement(lclRow);
    Teuchos::ArrayView<const GO> cols_av
      (numEnt == 0 ? nullptr : cols, numEnt);
    Teuchos::ArrayView<const Scalar> vals_av
      (numEnt == 0 ? nullptr : reinterpret_cast<const Scalar*> (vals), numEnt);

    // FIXME (mfh 23 Mar 2017) This is a work-around for less common
    // combine modes.  combineGlobalValues throws on error; it does
    // not return an error code.  Thus, if it returns, it succeeded.
    combineGlobalValues(gblRow, cols_av, vals_av, combMode,
                        prefix, debug, verbose);
    return numEnt;
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  combineGlobalValues(
    const GlobalOrdinal globalRowIndex,
    const Teuchos::ArrayView<const GlobalOrdinal>& columnIndices,
    const Teuchos::ArrayView<const Scalar>& values,
    const Tpetra::CombineMode combineMode,
    const char* const prefix,
    const bool debug,
    const bool verbose)
  {
    const char tfecfFuncName[] = "combineGlobalValues: ";

    if (isStaticGraph ()) {
      // INSERT doesn't make sense for a static graph, since you
      // aren't allowed to change the structure of the graph.
      // However, all the other combine modes work.
      if (combineMode == ADD) {
        sumIntoGlobalValues (globalRowIndex, columnIndices, values);
      }
      else if (combineMode == REPLACE) {
        replaceGlobalValues (globalRowIndex, columnIndices, values);
      }
      else if (combineMode == ABSMAX) {
        using ::Tpetra::Details::AbsMax;
        AbsMax<Scalar> f;
        this->template transformGlobalValues<AbsMax<Scalar> > (globalRowIndex,
                                                               columnIndices,
                                                               values, f);
      }
      else if (combineMode == INSERT) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (isStaticGraph() && combineMode == INSERT,
           std::invalid_argument, "INSERT combine mode is forbidden "
           "if the matrix has a static (const) graph (i.e., was "
           "constructed with the CrsMatrix constructor that takes a "
           "const CrsGraph pointer).");
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (true, std::logic_error, "Invalid combine mode; should "
           "never get here!  "
           "Please report this bug to the Tpetra developers.");
      }
    }
    else { // The matrix has a dynamic graph.
      if (combineMode == ADD || combineMode == INSERT) {
        // For a dynamic graph, all incoming column indices are
        // inserted into the target graph.  Duplicate indices will
        // have their values summed.  In this context, ADD and INSERT
        // are equivalent.  We need to call insertGlobalValues()
        // anyway if the column indices don't yet exist in this row,
        // so we just call insertGlobalValues() for both cases.
        insertGlobalValuesFilteredChecked(globalRowIndex,
          columnIndices, values, prefix, debug, verbose);
      }
      // FIXME (mfh 14 Mar 2012):
      //
      // Implementing ABSMAX or REPLACE for a dynamic graph would
      // require modifying assembly to attach a possibly different
      // combine mode to each inserted (i, j, A_ij) entry.  For
      // example, consider two different Export operations to the same
      // target CrsMatrix, the first with ABSMAX combine mode and the
      // second with REPLACE.  This isn't a common use case, so we
      // won't mess with it for now.
      else if (combineMode == ABSMAX) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          ! isStaticGraph () && combineMode == ABSMAX, std::logic_error,
          "ABSMAX combine mode when the matrix has a dynamic graph is not yet "
          "implemented.");
      }
      else if (combineMode == REPLACE) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          ! isStaticGraph () && combineMode == REPLACE, std::logic_error,
          "REPLACE combine mode when the matrix has a dynamic graph is not yet "
          "implemented.");
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          true, std::logic_error, "Should never get here!  Please report this "
          "bug to the Tpetra developers.");
      }
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  unpackAndCombine
  (const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& importLIDs,
   Kokkos::DualView<char*, buffer_device_type> imports,
   Kokkos::DualView<size_t*, buffer_device_type> numPacketsPerLID,
   const size_t constantNumPackets,
   const CombineMode combineMode)
  {
    using Details::Behavior;
    using Details::dualViewStatusToString;
    using Details::ProfilingRegion;
    using std::endl;
    const char tfecfFuncName[] = "unpackAndCombine: ";
    ProfilingRegion regionUAC ("Tpetra::CrsMatrix::unpackAndCombine");

    const bool debug = Behavior::debug("CrsMatrix");
    const bool verbose = Behavior::verbose("CrsMatrix");
    constexpr int numValidModes = 5;
    const CombineMode validModes[numValidModes] =
      {ADD, REPLACE, ABSMAX, INSERT, ZERO};
    const char* validModeNames[numValidModes] =
      {"ADD", "REPLACE", "ABSMAX", "INSERT", "ZERO"};

    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsMatrix", "unpackAndCombine");
      std::ostringstream os;
      os << *prefix << "Start:" << endl
         << *prefix << "  "
         << dualViewStatusToString (importLIDs, "importLIDs")
         << endl
         << *prefix << "  "
         << dualViewStatusToString (imports, "imports")
         << endl
         << *prefix << "  "
         << dualViewStatusToString (numPacketsPerLID, "numPacketsPerLID")
         << endl
         << *prefix << "  constantNumPackets: " << constantNumPackets
         << endl
         << *prefix << "  combineMode: " << combineModeToString (combineMode)
         << endl;
      std::cerr << os.str ();
    }

    if (debug) {
      if (std::find (validModes, validModes+numValidModes, combineMode) ==
          validModes+numValidModes) {
        std::ostringstream os;
        os << "Invalid combine mode.  Valid modes are {";
        for (int k = 0; k < numValidModes; ++k) {
          os << validModeNames[k];
          if (k < numValidModes - 1) {
            os << ", ";
          }
        }
        os << "}.";
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (true, std::invalid_argument, os.str ());
      }
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (importLIDs.extent(0) != numPacketsPerLID.extent(0),
         std::invalid_argument, "importLIDs.extent(0)="
         << importLIDs.extent(0)
         << " != numPacketsPerLID.extent(0)="
         << numPacketsPerLID.extent(0) << ".");
    }

    if (combineMode == ZERO) {
      return; // nothing to do
    }

    if (debug) {
      using Teuchos::reduceAll;
      std::unique_ptr<std::ostringstream> msg (new std::ostringstream ());
      int lclBad = 0;
      try {
        unpackAndCombineImpl(importLIDs, imports, numPacketsPerLID,
                             constantNumPackets, combineMode,
                             verbose);
      } catch (std::exception& e) {
        lclBad = 1;
        *msg << e.what ();
      }
      int gblBad = 0;
      const Teuchos::Comm<int>& comm = * (this->getComm ());
      reduceAll<int, int> (comm, Teuchos::REDUCE_MAX,
                           lclBad, Teuchos::outArg (gblBad));
      if (gblBad != 0) {
        // mfh 22 Oct 2017: 'prefix' might be null, since it is only
        // initialized in a debug build.  Thus, we get the process
        // rank again here.  This is an error message, so the small
        // run-time cost doesn't matter.  See #1887.
        std::ostringstream os;
        os << "Proc " << comm.getRank () << ": " << msg->str () << endl;
        msg = std::unique_ptr<std::ostringstream> (new std::ostringstream ());
        ::Tpetra::Details::gathervPrint (*msg, os.str (), comm);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (true, std::logic_error, std::endl << "unpackAndCombineImpl "
           "threw an exception on one or more participating processes: "
           << endl << msg->str ());
      }
    }
    else {
      unpackAndCombineImpl(importLIDs, imports, numPacketsPerLID,
                           constantNumPackets, combineMode,
                           verbose);
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done!" << endl
         << *prefix << "  "
         << dualViewStatusToString (importLIDs, "importLIDs")
         << endl
         << *prefix << "  "
         << dualViewStatusToString (imports, "imports")
         << endl
         << *prefix << "  "
         << dualViewStatusToString (numPacketsPerLID, "numPacketsPerLID")
         << endl;
      std::cerr << os.str ();
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  unpackAndCombineImpl(
    const Kokkos::DualView<const local_ordinal_type*,
      buffer_device_type>& importLIDs,
    Kokkos::DualView<char*, buffer_device_type> imports,
    Kokkos::DualView<size_t*, buffer_device_type> numPacketsPerLID,
    const size_t constantNumPackets,
    const CombineMode combineMode,
    const bool verbose)
  {
    Details::ProfilingRegion region_unpack_and_combine_impl(
      "Tpetra::CrsMatrix::unpackAndCombineImpl",
      "Import/Export"
    );
    using std::endl;
    const char tfecfFuncName[] = "unpackAndCombineImpl";
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsMatrix", tfecfFuncName);
      std::ostringstream os;
      os << *prefix << "isStaticGraph(): "
         << (isStaticGraph() ? "true" : "false")
         << ", importLIDs.extent(0): "
         << importLIDs.extent(0)
         << ", imports.extent(0): "
         << imports.extent(0)
         << ", numPacketsPerLID.extent(0): "
         << numPacketsPerLID.extent(0)
         << endl;
      std::cerr << os.str();
    }

    if (isStaticGraph ()) {
      using Details::unpackCrsMatrixAndCombineNew;
      unpackCrsMatrixAndCombineNew(*this, imports, numPacketsPerLID,
                                   importLIDs, constantNumPackets,
                                   combineMode);
    }
    else {
      {
        using padding_type = typename crs_graph_type::padding_type;
        std::unique_ptr<padding_type> padding;
        try {
          padding = myGraph_->computePaddingForCrsMatrixUnpack(
            importLIDs, imports, numPacketsPerLID, verbose);
        }
        catch (std::exception& e) {
          const auto rowMap = getRowMap();
          const auto comm = rowMap.is_null() ? Teuchos::null :
            rowMap->getComm();
          const int myRank = comm.is_null() ? -1 : comm->getRank();
          TEUCHOS_TEST_FOR_EXCEPTION
            (true, std::runtime_error, "Proc " << myRank << ": "
             "Tpetra::CrsGraph::computePaddingForCrsMatrixUnpack "
             "threw an exception: " << e.what());
        }
        if (verbose) {
          std::ostringstream os;
          os << *prefix << "Call applyCrsPadding" << endl;
          std::cerr << os.str();
        }
        applyCrsPadding(*padding, verbose);
      }
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Call unpackAndCombineImplNonStatic" << endl;
        std::cerr << os.str();
      }
      unpackAndCombineImplNonStatic(importLIDs, imports,
                                    numPacketsPerLID,
                                    constantNumPackets,
                                    combineMode);
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str();
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  unpackAndCombineImplNonStatic(
    const Kokkos::DualView<const local_ordinal_type*,
      buffer_device_type>& importLIDs,
    Kokkos::DualView<char*, buffer_device_type> imports,
    Kokkos::DualView<size_t*, buffer_device_type> numPacketsPerLID,
    const size_t constantNumPackets,
    const CombineMode combineMode)
  {
    using Kokkos::View;
    using Kokkos::subview;
    using Kokkos::MemoryUnmanaged;
    using Details::Behavior;
    using Details::castAwayConstDualView;
    using Details::create_mirror_view_from_raw_host_array;
    using Details::PackTraits;
    using Details::ScalarViewTraits;
    using std::endl;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    using ST = impl_scalar_type;
    using size_type = typename Teuchos::ArrayView<LO>::size_type;
    using HES =
      typename View<int*, device_type>::HostMirror::execution_space;
    using pair_type = std::pair<typename View<int*, HES>::size_type,
                                typename View<int*, HES>::size_type>;
    using gids_out_type = View<GO*, HES, MemoryUnmanaged>;
    using vals_out_type = View<ST*, HES, MemoryUnmanaged>;
    const char tfecfFuncName[] = "unpackAndCombineImplNonStatic";

    const bool debug = Behavior::debug("CrsMatrix");
    const bool verbose = Behavior::verbose("CrsMatrix");
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsMatrix", tfecfFuncName);
      std::ostringstream os;
      os << *prefix << endl; // we've already printed DualViews' statuses
      std::cerr << os.str ();
    }
    const char* const prefix_raw =
      verbose ? prefix.get()->c_str() : nullptr;

    const size_type numImportLIDs = importLIDs.extent (0);
    if (combineMode == ZERO || numImportLIDs == 0) {
      return; // nothing to do; no need to combine entries
    }

    Details::ProfilingRegion region_unpack_and_combine_impl_non_static(
      "Tpetra::CrsMatrix::unpackAndCombineImplNonStatic",
      "Import/Export"
    );

    // We're unpacking on host.  This is read-only host access.
    if (imports.need_sync_host()) {
      imports.sync_host ();
    }
    auto imports_h = imports.view_host();

    // Read-only host access.
    if (numPacketsPerLID.need_sync_host()) {
      numPacketsPerLID.sync_host ();
    }
    auto numPacketsPerLID_h = numPacketsPerLID.view_host();

    TEUCHOS_ASSERT( ! importLIDs.need_sync_host() );
    auto importLIDs_h = importLIDs.view_host();

    size_t numBytesPerValue;
    {
      // FIXME (mfh 17 Feb 2015, tjf 2 Aug 2017) What do I do about Scalar types
      // with run-time size?  We already assume that all entries in both the
      // source and target matrices have the same size.  If the calling process
      // owns at least one entry in either matrix, we can use that entry to set
      // the size.  However, it is possible that the calling process owns no
      // entries.  In that case, we're in trouble.  One way to fix this would be
      // for each row's data to contain the run-time size.  This is only
      // necessary if the size is not a compile-time constant.
      Scalar val;
      numBytesPerValue = PackTraits<ST>::packValueCount (val);
    }

    // Determine the maximum number of entries in any one row
    size_t offset = 0;
    size_t maxRowNumEnt = 0;
    for (size_type i = 0; i < numImportLIDs; ++i) {
      const size_t numBytes = numPacketsPerLID_h[i];
      if (numBytes == 0) {
        continue; // empty buffer for that row means that the row is empty
      }
      // We need to unpack a nonzero number of entries for this row.
      if (debug) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (offset + numBytes > size_t(imports_h.extent (0)),
           std::logic_error, ": At local row index importLIDs_h[i="
           << i << "]=" << importLIDs_h[i] << ", offset (=" << offset
           << ") + numBytes (=" << numBytes << ") > "
           "imports_h.extent(0)=" << imports_h.extent (0) << ".");
      }
      LO numEntLO = 0;

      if (debug) {
        const size_t theNumBytes =
          PackTraits<LO>::packValueCount (numEntLO);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (theNumBytes > numBytes, std::logic_error, ": theNumBytes="
           << theNumBytes << " > numBytes = " << numBytes << ".");
      }
      const char* const inBuf = imports_h.data () + offset;
      const size_t actualNumBytes =
        PackTraits<LO>::unpackValue (numEntLO, inBuf);

      if (debug) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (actualNumBytes > numBytes, std::logic_error, ": At i=" << i
           << ", actualNumBytes=" << actualNumBytes
           << " > numBytes=" << numBytes << ".");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (numEntLO == 0, std::logic_error, ": At local row index "
           "importLIDs_h[i=" << i << "]=" << importLIDs_h[i] << ", "
           "the number of entries read from the packed data is "
           "numEntLO=" << numEntLO << ", but numBytes=" << numBytes
           << " != 0.");
      }

      maxRowNumEnt = std::max(size_t(numEntLO), maxRowNumEnt);
      offset += numBytes;
    }

    // Temporary space to cache incoming global column indices and
    // values.  Column indices come in as global indices, in case the
    // source object's column Map differs from the target object's
    // (this's) column Map.
    View<GO*, HES> gblColInds;
    View<LO*, HES> lclColInds;
    View<ST*, HES> vals;
    {
      GO gid = 0;
      LO lid = 0;
      // FIXME (mfh 17 Feb 2015, tjf 2 Aug 2017) What do I do about Scalar types
      // with run-time size?  We already assume that all entries in both the
      // source and target matrices have the same size.  If the calling process
      // owns at least one entry in either matrix, we can use that entry to set
      // the size.  However, it is possible that the calling process owns no
      // entries.  In that case, we're in trouble.  One way to fix this would be
      // for each row's data to contain the run-time size.  This is only
      // necessary if the size is not a compile-time constant.
      Scalar val;
      gblColInds = ScalarViewTraits<GO, HES>::allocateArray(
        gid, maxRowNumEnt, "gids");
      lclColInds = ScalarViewTraits<LO, HES>::allocateArray(
        lid, maxRowNumEnt, "lids");
      vals = ScalarViewTraits<ST, HES>::allocateArray(
        val, maxRowNumEnt, "vals");
    }

    offset = 0;
    for (size_type i = 0; i < numImportLIDs; ++i) {
      const size_t numBytes = numPacketsPerLID_h[i];
      if (numBytes == 0) {
        continue; // empty buffer for that row means that the row is empty
      }
      LO numEntLO = 0;
      const char* const inBuf = imports_h.data () + offset;
      (void) PackTraits<LO>::unpackValue (numEntLO, inBuf);

      const size_t numEnt = static_cast<size_t>(numEntLO);;
      const LO lclRow = importLIDs_h[i];

      gids_out_type gidsOut = subview (gblColInds, pair_type (0, numEnt));
      vals_out_type valsOut = subview (vals, pair_type (0, numEnt));

      const size_t numBytesOut =
        unpackRow (gidsOut.data (), valsOut.data (), imports_h.data (),
                   offset, numBytes, numEnt, numBytesPerValue);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (numBytes != numBytesOut, std::logic_error, ": At i=" << i
         << ", numBytes=" << numBytes << " != numBytesOut="
         << numBytesOut << ".");

      const ST* const valsRaw = const_cast<const ST*> (valsOut.data ());
      const GO* const gidsRaw = const_cast<const GO*> (gidsOut.data ());
      combineGlobalValuesRaw(lclRow, numEnt, valsRaw, gidsRaw,
                             combineMode, prefix_raw, debug, verbose);
      // Don't update offset until current LID has succeeded.
      offset += numBytes;
    } // for each import LID i

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str();
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getColumnMapMultiVector (const MV& X_domainMap,
                           const bool force) const
  {
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! this->hasColMap (), std::runtime_error, "Tpetra::CrsMatrix::getColumn"
      "MapMultiVector: You may only call this method if the matrix has a "
      "column Map.  If the matrix does not yet have a column Map, you should "
      "first call fillComplete (with domain and range Map if necessary).");

    // If the graph is not fill complete, then the Import object (if
    // one should exist) hasn't been constructed yet.
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! this->getGraph ()->isFillComplete (), std::runtime_error, "Tpetra::"
      "CrsMatrix::getColumnMapMultiVector: You may only call this method if "
      "this matrix's graph is fill complete.");

    const size_t numVecs = X_domainMap.getNumVectors ();
    RCP<const import_type> importer = this->getGraph ()->getImporter ();
    RCP<const map_type> colMap = this->getColMap ();

    RCP<MV> X_colMap; // null by default

    // If the Import object is trivial (null), then we don't need a
    // separate column Map multivector.  Just return null in that
    // case.  The caller is responsible for knowing not to use the
    // returned null pointer.
    //
    // If the Import is nontrivial, then we do need a separate
    // column Map multivector for the Import operation.  Check in
    // that case if we have to (re)create the column Map
    // multivector.
    if (! importer.is_null () || force) {
      if (importMV_.is_null () || importMV_->getNumVectors () != numVecs) {
        X_colMap = rcp (new MV (colMap, numVecs));

        // Cache the newly created multivector for later reuse.
        importMV_ = X_colMap;
      }
      else { // Yay, we can reuse the cached multivector!
        X_colMap = importMV_;
        // mfh 09 Jan 2013: We don't have to fill with zeros first,
        // because the Import uses INSERT combine mode, which overwrites
        // existing entries.
        //
        //X_colMap->putScalar (ZERO);
      }
    }
    return X_colMap;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getRowMapMultiVector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y_rangeMap,
                        const bool force) const
  {
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;

    // If the graph is not fill complete, then the Export object (if
    // one should exist) hasn't been constructed yet.
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! this->getGraph ()->isFillComplete (), std::runtime_error, "Tpetra::"
      "CrsMatrix::getRowMapMultiVector: You may only call this method if this "
      "matrix's graph is fill complete.");

    const size_t numVecs = Y_rangeMap.getNumVectors ();
    RCP<const export_type> exporter = this->getGraph ()->getExporter ();
    // Every version of the constructor takes either a row Map, or a
    // graph (all of whose constructors take a row Map).  Thus, the
    // matrix always has a row Map.
    RCP<const map_type> rowMap = this->getRowMap ();

    RCP<MV> Y_rowMap; // null by default

    // If the Export object is trivial (null), then we don't need a
    // separate row Map multivector.  Just return null in that case.
    // The caller is responsible for knowing not to use the returned
    // null pointer.
    //
    // If the Export is nontrivial, then we do need a separate row
    // Map multivector for the Export operation.  Check in that case
    // if we have to (re)create the row Map multivector.
    if (! exporter.is_null () || force) {
      if (exportMV_.is_null () || exportMV_->getNumVectors () != numVecs) {
        Y_rowMap = rcp (new MV (rowMap, numVecs));
        exportMV_ = Y_rowMap; // Cache the newly created MV for later reuse.
      }
      else { // Yay, we can reuse the cached multivector!
        Y_rowMap = exportMV_;
      }
    }
    return Y_rowMap;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  removeEmptyProcessesInPlace (const Teuchos::RCP<const map_type>& newMap)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      myGraph_.is_null (), std::logic_error, "Tpetra::CrsMatrix::"
      "removeEmptyProcessesInPlace: This method does not work when the matrix "
      "was created with a constant graph (that is, when it was created using "
      "the version of its constructor that takes an RCP<const CrsGraph>).  "
      "This is because the matrix is not allowed to modify the graph in that "
      "case, but removing empty processes requires modifying the graph.");
    myGraph_->removeEmptyProcessesInPlace (newMap);
    // Even though CrsMatrix's row Map (as returned by getRowMap())
    // comes from its CrsGraph, CrsMatrix still implements DistObject,
    // so we also have to change the DistObject's Map.
    this->map_ = this->getRowMap ();
    // In the nonconst graph case, staticGraph_ is just a const
    // pointer to myGraph_.  This assignment is probably redundant,
    // but it doesn't hurt.
    staticGraph_ = Teuchos::rcp_const_cast<const Graph> (myGraph_);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  add (const Scalar& alpha,
       const RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
       const Scalar& beta,
       const Teuchos::RCP<const map_type>& domainMap,
       const Teuchos::RCP<const map_type>& rangeMap,
       const Teuchos::RCP<Teuchos::ParameterList>& params) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_implicit_cast;
    using Teuchos::sublist;
    using std::endl;
    using LO = local_ordinal_type;
    using GO = global_ordinal_type;
    using crs_matrix_type =
      CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    const char errPfx[] = "Tpetra::CrsMatrix::add: ";

    const bool debug = Details::Behavior::debug("CrsMatrix");
    const bool verbose = Details::Behavior::verbose("CrsMatrix");
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsMatrix", "add");
      std::ostringstream os;
      os << *prefix << "Start" << endl;
      std::cerr << os.str ();
    }

    const crs_matrix_type& B = *this; // a convenient abbreviation
    const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero();
    const Scalar ONE = Teuchos::ScalarTraits<Scalar>::one();

    // If the user didn't supply a domain or range Map, then try to
    // get one from B first (if it has them), then from A (if it has
    // them).  If we don't have any domain or range Maps, scold the
    // user.
    RCP<const map_type> A_domainMap = A.getDomainMap ();
    RCP<const map_type> A_rangeMap = A.getRangeMap ();
    RCP<const map_type> B_domainMap = B.getDomainMap ();
    RCP<const map_type> B_rangeMap = B.getRangeMap ();

    RCP<const map_type> theDomainMap = domainMap;
    RCP<const map_type> theRangeMap = rangeMap;

    if (domainMap.is_null ()) {
      if (B_domainMap.is_null ()) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          A_domainMap.is_null (), std::invalid_argument,
          "Tpetra::CrsMatrix::add: If neither A nor B have a domain Map, "
          "then you must supply a nonnull domain Map to this method.");
        theDomainMap = A_domainMap;
      } else {
        theDomainMap = B_domainMap;
      }
    }
    if (rangeMap.is_null ()) {
      if (B_rangeMap.is_null ()) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          A_rangeMap.is_null (), std::invalid_argument,
          "Tpetra::CrsMatrix::add: If neither A nor B have a range Map, "
          "then you must supply a nonnull range Map to this method.");
        theRangeMap = A_rangeMap;
      } else {
        theRangeMap = B_rangeMap;
      }
    }

    if (debug) {
      // In debug mode, check that A and B have matching domain and
      // range Maps, if they have domain and range Maps at all.  (If
      // they aren't fill complete, then they may not yet have them.)
      if (! A_domainMap.is_null() && ! A_rangeMap.is_null()) {
        if (! B_domainMap.is_null() && ! B_rangeMap.is_null()) {
          TEUCHOS_TEST_FOR_EXCEPTION
            (! B_domainMap->isSameAs(*A_domainMap),
             std::invalid_argument,
             errPfx << "The input RowMatrix A must have a domain Map "
             "which is the same as (isSameAs) this RowMatrix's "
             "domain Map.");
          TEUCHOS_TEST_FOR_EXCEPTION
            (! B_rangeMap->isSameAs(*A_rangeMap), std::invalid_argument,
             errPfx << "The input RowMatrix A must have a range Map "
             "which is the same as (isSameAs) this RowMatrix's range "
             "Map.");
          TEUCHOS_TEST_FOR_EXCEPTION
            (! domainMap.is_null() &&
             ! domainMap->isSameAs(*B_domainMap),
             std::invalid_argument,
             errPfx << "The input domain Map must be the same as "
             "(isSameAs) this RowMatrix's domain Map.");
          TEUCHOS_TEST_FOR_EXCEPTION
            (! rangeMap.is_null() &&
             ! rangeMap->isSameAs(*B_rangeMap),
             std::invalid_argument,
             errPfx << "The input range Map must be the same as "
             "(isSameAs) this RowMatrix's range Map.");
        }
      }
      else if (! B_domainMap.is_null() && ! B_rangeMap.is_null()) {
        TEUCHOS_TEST_FOR_EXCEPTION
          (! domainMap.is_null() &&
           ! domainMap->isSameAs(*B_domainMap),
           std::invalid_argument,
           errPfx << "The input domain Map must be the same as "
           "(isSameAs) this RowMatrix's domain Map.");
        TEUCHOS_TEST_FOR_EXCEPTION
          (! rangeMap.is_null() && ! rangeMap->isSameAs(*B_rangeMap),
           std::invalid_argument,
           errPfx << "The input range Map must be the same as "
           "(isSameAs) this RowMatrix's range Map.");
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION
          (domainMap.is_null() || rangeMap.is_null(),
           std::invalid_argument, errPfx << "If neither A nor B "
           "have a domain and range Map, then you must supply a "
           "nonnull domain and range Map to this method.");
      }
    }

    // What parameters do we pass to C's constructor?  Do we call
    // fillComplete on C after filling it?  And if so, what parameters
    // do we pass to C's fillComplete call?
    bool callFillComplete = true;
    RCP<ParameterList> constructorSublist;
    RCP<ParameterList> fillCompleteSublist;
    if (! params.is_null()) {
      callFillComplete =
        params->get("Call fillComplete", callFillComplete);
      constructorSublist = sublist(params, "Constructor parameters");
      fillCompleteSublist = sublist(params, "fillComplete parameters");
    }

    RCP<const map_type> A_rowMap = A.getRowMap ();
    RCP<const map_type> B_rowMap = B.getRowMap ();
    RCP<const map_type> C_rowMap = B_rowMap; // see discussion in documentation
    RCP<crs_matrix_type> C; // The result matrix.

    // If A and B's row Maps are the same, we can compute an upper
    // bound on the number of entries in each row of C, before
    // actually computing the sum.  A reasonable upper bound is the
    // sum of the two entry counts in each row.  
    if (A_rowMap->isSameAs (*B_rowMap)) {
      const LO localNumRows = static_cast<LO> (A_rowMap->getLocalNumElements ());
      Array<size_t> C_maxNumEntriesPerRow (localNumRows, 0);

      // Get the number of entries in each row of A.
      if (alpha != ZERO) {
        for (LO localRow = 0; localRow < localNumRows; ++localRow) {
          const size_t A_numEntries = A.getNumEntriesInLocalRow (localRow);
          C_maxNumEntriesPerRow[localRow] += A_numEntries;
        }
      }
      // Get the number of entries in each row of B.
      if (beta != ZERO) {
        for (LO localRow = 0; localRow < localNumRows; ++localRow) {
          const size_t B_numEntries = B.getNumEntriesInLocalRow (localRow);
          C_maxNumEntriesPerRow[localRow] += B_numEntries;
        }
      }
      // Construct the result matrix C.
      if (constructorSublist.is_null ()) {
        C = rcp (new crs_matrix_type (C_rowMap, C_maxNumEntriesPerRow ()));
      } else {
        C = rcp (new crs_matrix_type (C_rowMap, C_maxNumEntriesPerRow (),
                                      constructorSublist));
      }
      // Since A and B have the same row Maps, we could add them
      // together all at once and merge values before we call
      // insertGlobalValues.  However, we don't really need to, since
      // we've already allocated enough space in each row of C for C
      // to do the merge itself.
    }
    else { // the row Maps of A and B are not the same
      // Construct the result matrix C.
      // true: !A_rowMap->isSameAs (*B_rowMap)
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::invalid_argument, errPfx << "The row maps must "
         "be the same for statically allocated matrices, to ensure "
         "that there is sufficient space to do the addition.");
    }

    TEUCHOS_TEST_FOR_EXCEPTION
      (C.is_null (), std::logic_error,
       errPfx << "C should not be null at this point.  "
       "Please report this bug to the Tpetra developers.");

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Compute C = alpha*A + beta*B" << endl;
      std::cerr << os.str ();
    }
    using gids_type = nonconst_global_inds_host_view_type;
    using vals_type = nonconst_values_host_view_type;
    gids_type ind;
    vals_type val;

    if (alpha != ZERO) {
      const LO A_localNumRows = static_cast<LO> (A_rowMap->getLocalNumElements ());
      for (LO localRow = 0; localRow < A_localNumRows; ++localRow) {
        size_t A_numEntries = A.getNumEntriesInLocalRow (localRow);
        const GO globalRow = A_rowMap->getGlobalElement (localRow);
        if (A_numEntries > static_cast<size_t> (ind.size ())) {
          Kokkos::resize(ind,A_numEntries);
          Kokkos::resize(val,A_numEntries);
        }
        gids_type indView = Kokkos::subview(ind,std::make_pair((size_t)0, A_numEntries));
        vals_type valView = Kokkos::subview(val,std::make_pair((size_t)0, A_numEntries));
        A.getGlobalRowCopy (globalRow, indView, valView, A_numEntries);

        if (alpha != ONE) {
          for (size_t k = 0; k < A_numEntries; ++k) {
            valView[k] *= alpha;
          }
        }
        C->insertGlobalValues (globalRow, A_numEntries,
                               reinterpret_cast<Scalar *>(valView.data()),
                               indView.data());
      }
    }

    if (beta != ZERO) {
      const LO B_localNumRows = static_cast<LO> (B_rowMap->getLocalNumElements ());
      for (LO localRow = 0; localRow < B_localNumRows; ++localRow) {
        size_t B_numEntries = B.getNumEntriesInLocalRow (localRow);
        const GO globalRow = B_rowMap->getGlobalElement (localRow);
        if (B_numEntries > static_cast<size_t> (ind.size ())) {
          Kokkos::resize(ind,B_numEntries);
          Kokkos::resize(val,B_numEntries);
        }
        gids_type indView = Kokkos::subview(ind,std::make_pair((size_t)0, B_numEntries));
        vals_type valView = Kokkos::subview(val,std::make_pair((size_t)0, B_numEntries));
        B.getGlobalRowCopy (globalRow, indView, valView, B_numEntries);

        if (beta != ONE) {
          for (size_t k = 0; k < B_numEntries; ++k) {
            valView[k] *= beta;
          }
        }
        C->insertGlobalValues (globalRow, B_numEntries,
                               reinterpret_cast<Scalar *>(valView.data()),
                               indView.data());
      }
    }

    if (callFillComplete) {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Call fillComplete on C" << endl;
        std::cerr << os.str ();
      }
      if (fillCompleteSublist.is_null ()) {
        C->fillComplete (theDomainMap, theRangeMap);
      } else {
        C->fillComplete (theDomainMap, theRangeMap, fillCompleteSublist);
      }
    }
    else if (verbose) {
      std::ostringstream os;
      os << *prefix << "Do NOT call fillComplete on C" << endl;
      std::cerr << os.str ();
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str ();
    }
    return rcp_implicit_cast<row_matrix_type> (C);
  }



  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  transferAndFillComplete (Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > & destMat,
                           const ::Tpetra::Details::Transfer<LocalOrdinal, GlobalOrdinal, Node>& rowTransfer,
                           const Teuchos::RCP<const ::Tpetra::Details::Transfer<LocalOrdinal, GlobalOrdinal, Node> > & domainTransfer,
                           const Teuchos::RCP<const map_type>& domainMap,
                           const Teuchos::RCP<const map_type>& rangeMap,
                           const Teuchos::RCP<Teuchos::ParameterList>& params) const
  {
    using Details::Behavior;
    using Details::getArrayViewFromDualView;
    using Details::packCrsMatrixWithOwningPIDs;
    using Details::unpackAndCombineWithOwningPIDsCount;
    using Details::unpackAndCombineIntoCrsArrays;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::Comm;
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using std::endl;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef node_type NT;
    typedef CrsMatrix<Scalar, LO, GO, NT> this_CRS_type;
    typedef Vector<int, LO, GO, NT> IntVectorType;
    using Teuchos::as;

    const bool debug = Behavior::debug("CrsMatrix");
    const bool verbose = Behavior::verbose("CrsMatrix");
    int MyPID = getComm ()->getRank ();

    std::unique_ptr<std::string> verbosePrefix;
    if (verbose) {
      verbosePrefix =
        this->createPrefix("CrsMatrix", "transferAndFillComplete");
      std::ostringstream os;
      os << "Start" << endl;
      std::cerr << os.str();
    }

    //
    // Get the caller's parameters
    //
    bool isMM = false; // optimize for matrix-matrix ops.
    bool reverseMode = false; // Are we in reverse mode?
    bool restrictComm = false; // Do we need to restrict the communicator?

    int mm_optimization_core_count =
      Behavior::TAFC_OptimizationCoreCount();
    RCP<ParameterList> matrixparams; // parameters for the destination matrix
    bool overrideAllreduce = false;
    bool useKokkosPath = false;
    if (! params.is_null ()) {
      matrixparams = sublist (params, "CrsMatrix");
      reverseMode = params->get ("Reverse Mode", reverseMode);
      useKokkosPath = params->get ("TAFC: use kokkos path", useKokkosPath);
      restrictComm = params->get ("Restrict Communicator", restrictComm);
      auto & slist = params->sublist("matrixmatrix: kernel params",false);
      isMM = slist.get("isMatrixMatrix_TransferAndFillComplete",false);
      mm_optimization_core_count = slist.get("MM_TAFC_OptimizationCoreCount",mm_optimization_core_count);

      overrideAllreduce = slist.get("MM_TAFC_OverrideAllreduceCheck",false);
      if(getComm()->getSize() < mm_optimization_core_count && isMM)   isMM = false;
      if(reverseMode) isMM = false;
    }

   // Only used in the sparse matrix-matrix multiply (isMM) case.
   std::shared_ptr< ::Tpetra::Details::CommRequest> iallreduceRequest;
   int mismatch = 0;
   int reduced_mismatch = 0;
   if (isMM && !overrideAllreduce) {

     // Test for pathological matrix transfer
     const bool source_vals = ! getGraph ()->getImporter ().is_null();
     const bool target_vals = ! (rowTransfer.getExportLIDs ().size() == 0 ||
                                 rowTransfer.getRemoteLIDs ().size() == 0);
     mismatch = (source_vals != target_vals) ? 1 : 0;
     iallreduceRequest =
       ::Tpetra::Details::iallreduce (mismatch, reduced_mismatch,
                                      Teuchos::REDUCE_MAX, * (getComm ()));
   }

#ifdef HAVE_TPETRA_MMM_TIMINGS
    using Teuchos::TimeMonitor;
    std::string label;
    if(!params.is_null())
        label = params->get("Timer Label",label);
    std::string prefix = std::string("Tpetra ")+ label + std::string(": ");
    std::string tlstr;
    {
        std::ostringstream os;
        if(isMM) os<<":MMOpt";
        else os<<":MMLegacy";
        tlstr = os.str();
    }

    Teuchos::TimeMonitor MMall(*TimeMonitor::getNewTimer(prefix + std::string("TAFC All") +tlstr ));
#endif

    // Make sure that the input argument rowTransfer is either an
    // Import or an Export.  Import and Export are the only two
    // subclasses of Transfer that we defined, but users might
    // (unwisely, for now at least) decide to implement their own
    // subclasses.  Exclude this possibility.
    const import_type* xferAsImport = dynamic_cast<const import_type*> (&rowTransfer);
    const export_type* xferAsExport = dynamic_cast<const export_type*> (&rowTransfer);
    TEUCHOS_TEST_FOR_EXCEPTION(
      xferAsImport == nullptr && xferAsExport == nullptr, std::invalid_argument,
      "Tpetra::CrsMatrix::transferAndFillComplete: The 'rowTransfer' input "
      "argument must be either an Import or an Export, and its template "
      "parameters must match the corresponding template parameters of the "
      "CrsMatrix.");

    // Make sure that the input argument domainTransfer is either an
    // Import or an Export.  Import and Export are the only two
    // subclasses of Transfer that we defined, but users might
    // (unwisely, for now at least) decide to implement their own
    // subclasses.  Exclude this possibility.
    Teuchos::RCP<const import_type> xferDomainAsImport = Teuchos::rcp_dynamic_cast<const import_type> (domainTransfer);
    Teuchos::RCP<const export_type> xferDomainAsExport = Teuchos::rcp_dynamic_cast<const export_type> (domainTransfer);

    if(! domainTransfer.is_null()) {
      TEUCHOS_TEST_FOR_EXCEPTION(
         (xferDomainAsImport.is_null() && xferDomainAsExport.is_null()), std::invalid_argument,
        "Tpetra::CrsMatrix::transferAndFillComplete: The 'domainTransfer' input "
        "argument must be either an Import or an Export, and its template "
        "parameters must match the corresponding template parameters of the "
        "CrsMatrix.");

      TEUCHOS_TEST_FOR_EXCEPTION(
         ( xferAsImport != nullptr || ! xferDomainAsImport.is_null() ) &&
         (( xferAsImport != nullptr &&   xferDomainAsImport.is_null() ) ||
          ( xferAsImport == nullptr && ! xferDomainAsImport.is_null() )), std::invalid_argument,
         "Tpetra::CrsMatrix::transferAndFillComplete: The 'rowTransfer' and 'domainTransfer' input "
         "arguments must be of the same type (either Import or Export).");

      TEUCHOS_TEST_FOR_EXCEPTION(
         ( xferAsExport != nullptr || ! xferDomainAsExport.is_null() ) &&
         (( xferAsExport != nullptr &&   xferDomainAsExport.is_null() ) ||
          ( xferAsExport == nullptr && ! xferDomainAsExport.is_null() )), std::invalid_argument,
         "Tpetra::CrsMatrix::transferAndFillComplete: The 'rowTransfer' and 'domainTransfer' input "
         "arguments must be of the same type (either Import or Export).");
    } // domainTransfer != null


    // FIXME (mfh 15 May 2014) Wouldn't communication still be needed,
    // if the source Map is not distributed but the target Map is?
    const bool communication_needed = rowTransfer.getSourceMap ()->isDistributed ();

    // Get the new domain and range Maps.  We need some of them for
    // error checking, now that we have the reverseMode parameter.
    RCP<const map_type> MyRowMap = reverseMode ?
      rowTransfer.getSourceMap () : rowTransfer.getTargetMap ();
    RCP<const map_type> MyColMap; // create this below
    RCP<const map_type> MyDomainMap = ! domainMap.is_null () ?
      domainMap : getDomainMap ();
    RCP<const map_type> MyRangeMap = ! rangeMap.is_null () ?
      rangeMap : getRangeMap ();
    RCP<const map_type> BaseRowMap = MyRowMap;
    RCP<const map_type> BaseDomainMap = MyDomainMap;

    // If the user gave us a nonnull destMat, then check whether it's
    // "pristine."  That means that it has no entries.
    //
    // FIXME (mfh 15 May 2014) If this is not true on all processes,
    // then this exception test may hang.  It would be better to
    // forward an error flag to the next communication phase.
    if (! destMat.is_null ()) {
      // FIXME (mfh 15 May 2014): The Epetra idiom for checking
      // whether a graph or matrix has no entries on the calling
      // process, is that it is neither locally nor globally indexed.
      // This may change eventually with the Kokkos refactor version
      // of Tpetra, so it would be better just to check the quantity
      // of interest directly.  Note that with the Kokkos refactor
      // version of Tpetra, asking for the total number of entries in
      // a graph or matrix that is not fill complete might require
      // computation (kernel launch), since it is not thread scalable
      // to update a count every time an entry is inserted.
      const bool NewFlag = ! destMat->getGraph ()->isLocallyIndexed () &&
        ! destMat->getGraph ()->isGloballyIndexed ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! NewFlag, std::invalid_argument, "Tpetra::CrsMatrix::"
        "transferAndFillComplete: The input argument 'destMat' is only allowed "
        "to be nonnull, if its graph is empty (neither locally nor globally "
        "indexed).");
      // FIXME (mfh 15 May 2014) At some point, we want to change
      // graphs and matrices so that their DistObject Map
      // (this->getMap()) may differ from their row Map.  This will
      // make redistribution for 2-D distributions more efficient.  I
      // hesitate to change this check, because I'm not sure how much
      // the code here depends on getMap() and getRowMap() being the
      // same.
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! destMat->getRowMap ()->isSameAs (*MyRowMap), std::invalid_argument,
        "Tpetra::CrsMatrix::transferAndFillComplete: The (row) Map of the "
        "input argument 'destMat' is not the same as the (row) Map specified "
        "by the input argument 'rowTransfer'.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! destMat->checkSizes (*this), std::invalid_argument,
        "Tpetra::CrsMatrix::transferAndFillComplete: You provided a nonnull "
        "destination matrix, but checkSizes() indicates that it is not a legal "
        "legal target for redistribution from the source matrix (*this).  This "
        "may mean that they do not have the same dimensions.");
    }

    // If forward mode (the default), then *this's (row) Map must be
    // the same as the source Map of the Transfer.  If reverse mode,
    // then *this's (row) Map must be the same as the target Map of
    // the Transfer.
    //
    // FIXME (mfh 15 May 2014) At some point, we want to change graphs
    // and matrices so that their DistObject Map (this->getMap()) may
    // differ from their row Map.  This will make redistribution for
    // 2-D distributions more efficient.  I hesitate to change this
    // check, because I'm not sure how much the code here depends on
    // getMap() and getRowMap() being the same.
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! (reverseMode || getRowMap ()->isSameAs (*rowTransfer.getSourceMap ())),
      std::invalid_argument, "Tpetra::CrsMatrix::transferAndFillComplete: "
      "rowTransfer->getSourceMap() must match this->getRowMap() in forward mode.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! (! reverseMode || getRowMap ()->isSameAs (*rowTransfer.getTargetMap ())),
      std::invalid_argument, "Tpetra::CrsMatrix::transferAndFillComplete: "
      "rowTransfer->getTargetMap() must match this->getRowMap() in reverse mode.");

    // checks for domainTransfer
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! xferDomainAsImport.is_null() && ! xferDomainAsImport->getTargetMap()->isSameAs(*domainMap),
      std::invalid_argument,
      "Tpetra::CrsMatrix::transferAndFillComplete: The target map of the 'domainTransfer' input "
      "argument must be the same as the rebalanced domain map 'domainMap'");

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! xferDomainAsExport.is_null() && ! xferDomainAsExport->getSourceMap()->isSameAs(*domainMap),
      std::invalid_argument,
      "Tpetra::CrsMatrix::transferAndFillComplete: The source map of the 'domainTransfer' input "
      "argument must be the same as the rebalanced domain map 'domainMap'");

    // The basic algorithm here is:
    //
    // 1. Call the moral equivalent of "Distor.do" to handle the import.
    // 2. Copy all the Imported and Copy/Permuted data into the raw
    //    CrsMatrix / CrsGraphData pointers, still using GIDs.
    // 3. Call an optimized version of MakeColMap that avoids the
    //    Directory lookups (since the importer knows who owns all the
    //    GIDs) AND reindexes to LIDs.
    // 4. Call expertStaticFillComplete()

    // Get information from the Importer
    const size_t NumSameIDs = rowTransfer.getNumSameIDs();
    ArrayView<const LO> ExportLIDs = reverseMode ?
      rowTransfer.getRemoteLIDs () : rowTransfer.getExportLIDs ();
    auto RemoteLIDs = reverseMode ?
      rowTransfer.getExportLIDs_dv() : rowTransfer.getRemoteLIDs_dv();
    auto PermuteToLIDs = reverseMode ?
      rowTransfer.getPermuteFromLIDs_dv() : rowTransfer.getPermuteToLIDs_dv();
    auto PermuteFromLIDs = reverseMode ?
      rowTransfer.getPermuteToLIDs_dv() : rowTransfer.getPermuteFromLIDs_dv();
    Distributor& Distor = rowTransfer.getDistributor ();

    // Owning PIDs
    Teuchos::Array<int> SourcePids;

    // Temp variables for sub-communicators
    RCP<const map_type> ReducedRowMap, ReducedColMap,
      ReducedDomainMap, ReducedRangeMap;
    RCP<const Comm<int> > ReducedComm;

    // If the user gave us a null destMat, then construct the new
    // destination matrix.  We will replace its column Map later.
    if (destMat.is_null ()) {
      destMat = rcp (new this_CRS_type (MyRowMap, 0, matrixparams));
    }

    /***************************************************/
    /***** 1) First communicator restriction phase ****/
    /***************************************************/
    if (restrictComm) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC restrictComm")));
#endif
      ReducedRowMap = MyRowMap->removeEmptyProcesses ();
      ReducedComm = ReducedRowMap.is_null () ?
        Teuchos::null :
        ReducedRowMap->getComm ();
      destMat->removeEmptyProcessesInPlace (ReducedRowMap);

      ReducedDomainMap = MyRowMap.getRawPtr () == MyDomainMap.getRawPtr () ?
        ReducedRowMap :
        MyDomainMap->replaceCommWithSubset (ReducedComm);
      ReducedRangeMap = MyRowMap.getRawPtr () == MyRangeMap.getRawPtr () ?
        ReducedRowMap :
        MyRangeMap->replaceCommWithSubset (ReducedComm);

      // Reset the "my" maps
      MyRowMap    = ReducedRowMap;
      MyDomainMap = ReducedDomainMap;
      MyRangeMap  = ReducedRangeMap;

      // Update my PID, if we've restricted the communicator
      if (! ReducedComm.is_null ()) {
        MyPID = ReducedComm->getRank ();
      }
      else {
        MyPID = -2; // For debugging
      }
    }
    else {
      ReducedComm = MyRowMap->getComm ();
    }



    /***************************************************/
    /***** 2) From Tpetra::DistObject::doTransfer() ****/
    /***************************************************/
    // Get the owning PIDs
    RCP<const import_type> MyImporter = getGraph ()->getImporter ();

    // check whether domain maps of source matrix and base domain map is the same
    bool bSameDomainMap = BaseDomainMap->isSameAs (*getDomainMap ());

    if (! restrictComm && ! MyImporter.is_null () && bSameDomainMap ) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC getOwningPIDs same map")));
#endif
      // Same domain map as source matrix
      //
      // NOTE: This won't work for restrictComm (because the Import
      // doesn't know the restricted PIDs), though writing an
      // optimized version for that case would be easy (Import an
      // IntVector of the new PIDs).  Might want to add this later.
      Import_Util::getPids (*MyImporter, SourcePids, false);
    }
    else if (restrictComm && ! MyImporter.is_null () && bSameDomainMap) {
      // Same domain map as source matrix (restricted communicator)
      // We need one import from the domain to the column map
#ifdef HAVE_TPETRA_MMM_TIMINGS
      Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC getOwningPIDs restricted comm")));
#endif
      IntVectorType SourceDomain_pids(getDomainMap (),true);
      IntVectorType SourceCol_pids(getColMap());
      // SourceDomain_pids contains the restricted pids
      SourceDomain_pids.putScalar(MyPID);

      SourceCol_pids.doImport (SourceDomain_pids, *MyImporter, INSERT);
      SourcePids.resize (getColMap ()->getLocalNumElements ());
      SourceCol_pids.get1dCopy (SourcePids ());
    }
    else if (MyImporter.is_null ()) {
      // Matrix has no off-process entries
#ifdef HAVE_TPETRA_MMM_TIMINGS
      Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC getOwningPIDs all local entries")));
#endif
      SourcePids.resize (getColMap ()->getLocalNumElements ());
      SourcePids.assign (getColMap ()->getLocalNumElements (), MyPID);
    }
    else if ( ! MyImporter.is_null () &&
              ! domainTransfer.is_null () ) {
      // general implementation for rectangular matrices with
      // domain map different than SourceMatrix domain map.
      // User has to provide a DomainTransfer object. We need
      // to communications (import/export)
#ifdef HAVE_TPETRA_MMM_TIMINGS
      Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC getOwningPIDs rectangular case")));
#endif

      // TargetDomain_pids lives on the rebalanced new domain map
      IntVectorType TargetDomain_pids (domainMap);
      TargetDomain_pids.putScalar (MyPID);

      // SourceDomain_pids lives on the non-rebalanced old domain map
      IntVectorType SourceDomain_pids (getDomainMap ());

      // SourceCol_pids lives on the non-rebalanced old column map
      IntVectorType SourceCol_pids (getColMap ());

      if (! reverseMode && ! xferDomainAsImport.is_null() ) {
        SourceDomain_pids.doExport (TargetDomain_pids, *xferDomainAsImport, INSERT);
      }
      else if (reverseMode && ! xferDomainAsExport.is_null() ) {
        SourceDomain_pids.doExport (TargetDomain_pids, *xferDomainAsExport, INSERT);
      }
      else if (! reverseMode && ! xferDomainAsExport.is_null() ) {
        SourceDomain_pids.doImport (TargetDomain_pids, *xferDomainAsExport, INSERT);
      }
      else if (reverseMode && ! xferDomainAsImport.is_null() ) {
        SourceDomain_pids.doImport (TargetDomain_pids, *xferDomainAsImport, INSERT);
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error, "Tpetra::CrsMatrix::"
          "transferAndFillComplete: Should never get here!  "
          "Please report this bug to a Tpetra developer.");
      }
      SourceCol_pids.doImport (SourceDomain_pids, *MyImporter, INSERT);
      SourcePids.resize (getColMap ()->getLocalNumElements ());
      SourceCol_pids.get1dCopy (SourcePids ());
    }
    else if ( ! MyImporter.is_null () &&
             BaseDomainMap->isSameAs (*BaseRowMap) &&
             getDomainMap ()->isSameAs (*getRowMap ())) {
      // We can use the rowTransfer + SourceMatrix's Import to find out who owns what.
#ifdef HAVE_TPETRA_MMM_TIMINGS
      Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC getOwningPIDs query import")));
#endif

      IntVectorType TargetRow_pids (domainMap);
      IntVectorType SourceRow_pids (getRowMap ());
      IntVectorType SourceCol_pids (getColMap ());

      TargetRow_pids.putScalar (MyPID);
      if (! reverseMode && xferAsImport != nullptr) {
        SourceRow_pids.doExport (TargetRow_pids, *xferAsImport, INSERT);
      }
      else if (reverseMode && xferAsExport != nullptr) {
        SourceRow_pids.doExport (TargetRow_pids, *xferAsExport, INSERT);
      }
      else if (! reverseMode && xferAsExport != nullptr) {
        SourceRow_pids.doImport (TargetRow_pids, *xferAsExport, INSERT);
      }
      else if (reverseMode && xferAsImport != nullptr) {
        SourceRow_pids.doImport (TargetRow_pids, *xferAsImport, INSERT);
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error, "Tpetra::CrsMatrix::"
          "transferAndFillComplete: Should never get here!  "
          "Please report this bug to a Tpetra developer.");
      }

      SourceCol_pids.doImport (SourceRow_pids, *MyImporter, INSERT);
      SourcePids.resize (getColMap ()->getLocalNumElements ());
      SourceCol_pids.get1dCopy (SourcePids ());
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::CrsMatrix::"
        "transferAndFillComplete: This method only allows either domainMap == "
        "getDomainMap (), or (domainMap == rowTransfer.getTargetMap () and "
        "getDomainMap () == getRowMap ()).");
    }

    // Tpetra-specific stuff
    size_t constantNumPackets = destMat->constantNumberOfPackets ();
    {
#ifdef HAVE_TPETRA_MMM_TIMINGS
    Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC reallocate buffers")));
#endif
    if (constantNumPackets == 0) {
      destMat->reallocArraysForNumPacketsPerLid (ExportLIDs.size (),
                                                 RemoteLIDs.view_host().size ());
    }
    else {
      // There are a constant number of packets per element.  We
      // already know (from the number of "remote" (incoming)
      // elements) how many incoming elements we expect, so we can
      // resize the buffer accordingly.
      const size_t rbufLen = RemoteLIDs.view_host().size() * constantNumPackets;
      destMat->reallocImportsIfNeeded (rbufLen, false, nullptr);
    }
    }

    // Pack & Prepare w/ owning PIDs
    {
#ifdef HAVE_TPETRA_MMM_TIMINGS
    Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC pack and prepare")));
#endif
    if (debug) {
      using Teuchos::outArg;
      using Teuchos::REDUCE_MAX;
      using Teuchos::reduceAll;
      using std::cerr;
      using std::endl;
      RCP<const Teuchos::Comm<int> > comm = this->getComm ();
      const int myRank = comm->getRank ();

      std::ostringstream errStrm;
      int lclErr = 0;
      int gblErr = 0;

      Teuchos::ArrayView<size_t> numExportPacketsPerLID;
      try {
        // packAndPrepare* methods modify numExportPacketsPerLID_.
        destMat->numExportPacketsPerLID_.modify_host ();
        numExportPacketsPerLID =
          getArrayViewFromDualView (destMat->numExportPacketsPerLID_);
      }
      catch (std::exception& e) {
        errStrm << "Proc " << myRank << ": getArrayViewFromDualView threw: "
                << e.what () << std::endl;
        lclErr = 1;
      }
      catch (...) {
        errStrm << "Proc " << myRank << ": getArrayViewFromDualView threw "
          "an exception not a subclass of std::exception" << std::endl;
        lclErr = 1;
      }

      if (! comm.is_null ()) {
        reduceAll<int, int> (*comm, REDUCE_MAX, lclErr, outArg (gblErr));
      }
      if (gblErr != 0) {
        ::Tpetra::Details::gathervPrint (cerr, errStrm.str (), *comm);
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error, "getArrayViewFromDualView threw an "
          "exception on at least one process.");
      }

      if (verbose) {
        std::ostringstream os;
        os << *verbosePrefix << "Calling packCrsMatrixWithOwningPIDs"
           << std::endl;
        std::cerr << os.str ();
      }
      try {
        packCrsMatrixWithOwningPIDs (*this,
                                     destMat->exports_,
                                     numExportPacketsPerLID,
                                     ExportLIDs,
                                     SourcePids,
                                     constantNumPackets);
      }
      catch (std::exception& e) {
        errStrm << "Proc " << myRank << ": packCrsMatrixWithOwningPIDs threw: "
           << e.what () << std::endl;
        lclErr = 1;
      }
      catch (...) {
        errStrm << "Proc " << myRank << ": packCrsMatrixWithOwningPIDs threw "
          "an exception not a subclass of std::exception" << std::endl;
        lclErr = 1;
      }

      if (verbose) {
        std::ostringstream os;
        os << *verbosePrefix << "Done with packCrsMatrixWithOwningPIDs"
           << std::endl;
        std::cerr << os.str ();
      }

      if (! comm.is_null ()) {
        reduceAll<int, int> (*comm, REDUCE_MAX, lclErr, outArg (gblErr));
      }
      if (gblErr != 0) {
        ::Tpetra::Details::gathervPrint (cerr, errStrm.str (), *comm);
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error, "packCrsMatrixWithOwningPIDs threw an "
          "exception on at least one process.");
      }
    }
    else {
      // packAndPrepare* methods modify numExportPacketsPerLID_.
      destMat->numExportPacketsPerLID_.modify_host ();
      Teuchos::ArrayView<size_t> numExportPacketsPerLID =
        getArrayViewFromDualView (destMat->numExportPacketsPerLID_);
      if (verbose) {
        std::ostringstream os;
        os << *verbosePrefix << "Calling packCrsMatrixWithOwningPIDs"
           << std::endl;
        std::cerr << os.str ();
      }
      packCrsMatrixWithOwningPIDs (*this,
                                   destMat->exports_,
                                   numExportPacketsPerLID,
                                   ExportLIDs,
                                   SourcePids,
                                   constantNumPackets);
      if (verbose) {
        std::ostringstream os;
        os << *verbosePrefix << "Done with packCrsMatrixWithOwningPIDs"
           << std::endl;
        std::cerr << os.str ();
      }
    }
    }

    // Do the exchange of remote data.
    {
#ifdef HAVE_TPETRA_MMM_TIMINGS
    Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC getOwningPIDs exchange remote data")));
#endif
    if (! communication_needed) {
      if (verbose) {
        std::ostringstream os;
        os << *verbosePrefix << "Communication not needed" << std::endl;
        std::cerr << os.str ();
      }
    }
    else {
      if (reverseMode) {
        if (constantNumPackets == 0) { // variable number of packets per LID
          if (verbose) {
            std::ostringstream os;
            os << *verbosePrefix << "Reverse mode, variable # packets / LID"
               << std::endl;
            std::cerr << os.str ();
          }
          // Make sure that host has the latest version, since we're
          // using the version on host.  If host has the latest
          // version, syncing to host does nothing.
          destMat->numExportPacketsPerLID_.sync_host ();
          Teuchos::ArrayView<const size_t> numExportPacketsPerLID =
            getArrayViewFromDualView (destMat->numExportPacketsPerLID_);
          destMat->numImportPacketsPerLID_.sync_host ();
          Teuchos::ArrayView<size_t> numImportPacketsPerLID =
            getArrayViewFromDualView (destMat->numImportPacketsPerLID_);

          if (verbose) {
            std::ostringstream os;
            os << *verbosePrefix << "Calling 3-arg doReversePostsAndWaits"
               << std::endl;
            std::cerr << os.str ();
          }
          Distor.doReversePostsAndWaits(destMat->numExportPacketsPerLID_.view_host(), 1,
                                            destMat->numImportPacketsPerLID_.view_host());
          if (verbose) {
            std::ostringstream os;
            os << *verbosePrefix << "Finished 3-arg doReversePostsAndWaits"
               << std::endl;
            std::cerr << os.str ();
          }

          size_t totalImportPackets = 0;
          for (Array_size_type i = 0; i < numImportPacketsPerLID.size (); ++i) {
            totalImportPackets += numImportPacketsPerLID[i];
          }

          // Reallocation MUST go before setting the modified flag,
          // because it may clear out the flags.
          destMat->reallocImportsIfNeeded (totalImportPackets, verbose,
                                           verbosePrefix.get ());
          destMat->imports_.modify_host ();
          auto hostImports = destMat->imports_.view_host();
          // This is a legacy host pack/unpack path, so use the host
          // version of exports_.
          destMat->exports_.sync_host ();
          auto hostExports = destMat->exports_.view_host();
          if (verbose) {
            std::ostringstream os;
            os << *verbosePrefix << "Calling 4-arg doReversePostsAndWaits"
               << std::endl;
            std::cerr << os.str ();
          }
          Distor.doReversePostsAndWaits (hostExports,
                                         numExportPacketsPerLID,
                                         hostImports,
                                         numImportPacketsPerLID);
          if (verbose) {
            std::ostringstream os;
            os << *verbosePrefix << "Finished 4-arg doReversePostsAndWaits"
               << std::endl;
            std::cerr << os.str ();
          }
        }
        else { // constant number of packets per LID
          if (verbose) {
            std::ostringstream os;
            os << *verbosePrefix << "Reverse mode, constant # packets / LID"
               << std::endl;
            std::cerr << os.str ();
          }
          destMat->imports_.modify_host ();
          auto hostImports = destMat->imports_.view_host();
          // This is a legacy host pack/unpack path, so use the host
          // version of exports_.
          destMat->exports_.sync_host ();
          auto hostExports = destMat->exports_.view_host();
          if (verbose) {
            std::ostringstream os;
            os << *verbosePrefix << "Calling 3-arg doReversePostsAndWaits"
               << std::endl;
            std::cerr << os.str ();
          }
          Distor.doReversePostsAndWaits (hostExports,
                                         constantNumPackets,
                                         hostImports);
          if (verbose) {
            std::ostringstream os;
            os << *verbosePrefix << "Finished 3-arg doReversePostsAndWaits"
               << std::endl;
            std::cerr << os.str ();
          }
        }
      }
      else { // forward mode (the default)
        if (constantNumPackets == 0) { // variable number of packets per LID
          if (verbose) {
            std::ostringstream os;
            os << *verbosePrefix << "Forward mode, variable # packets / LID"
               << std::endl;
            std::cerr << os.str ();
          }
          // Make sure that host has the latest version, since we're
          // using the version on host.  If host has the latest
          // version, syncing to host does nothing.
          destMat->numExportPacketsPerLID_.sync_host ();
          Teuchos::ArrayView<const size_t> numExportPacketsPerLID =
            getArrayViewFromDualView (destMat->numExportPacketsPerLID_);
          destMat->numImportPacketsPerLID_.sync_host ();
          Teuchos::ArrayView<size_t> numImportPacketsPerLID =
            getArrayViewFromDualView (destMat->numImportPacketsPerLID_);
          if (verbose) {
            std::ostringstream os;
            os << *verbosePrefix << "Calling 3-arg doPostsAndWaits"
               << std::endl;
            std::cerr << os.str ();
          }
          Distor.doPostsAndWaits(destMat->numExportPacketsPerLID_.view_host(), 1,
                                      destMat->numImportPacketsPerLID_.view_host());
          if (verbose) {
            std::ostringstream os;
            os << *verbosePrefix << "Finished 3-arg doPostsAndWaits"
               << std::endl;
            std::cerr << os.str ();
          }

          size_t totalImportPackets = 0;
          for (Array_size_type i = 0; i < numImportPacketsPerLID.size (); ++i) {
            totalImportPackets += numImportPacketsPerLID[i];
          }

          // Reallocation MUST go before setting the modified flag,
          // because it may clear out the flags.
          destMat->reallocImportsIfNeeded (totalImportPackets, verbose,
                                           verbosePrefix.get ());
          destMat->imports_.modify_host ();
          auto hostImports = destMat->imports_.view_host();
          // This is a legacy host pack/unpack path, so use the host
          // version of exports_.
          destMat->exports_.sync_host ();
          auto hostExports = destMat->exports_.view_host();
          if (verbose) {
            std::ostringstream os;
            os << *verbosePrefix << "Calling 4-arg doPostsAndWaits"
               << std::endl;
            std::cerr << os.str ();
          }
          Distor.doPostsAndWaits (hostExports,
                                  numExportPacketsPerLID,
                                  hostImports,
                                  numImportPacketsPerLID);
          if (verbose) {
            std::ostringstream os;
            os << *verbosePrefix << "Finished 4-arg doPostsAndWaits"
               << std::endl;
            std::cerr << os.str ();
          }
        }
        else { // constant number of packets per LID
          if (verbose) {
            std::ostringstream os;
            os << *verbosePrefix << "Forward mode, constant # packets / LID"
               << std::endl;
            std::cerr << os.str ();
          }
          destMat->imports_.modify_host ();
          auto hostImports = destMat->imports_.view_host();
          // This is a legacy host pack/unpack path, so use the host
          // version of exports_.
          destMat->exports_.sync_host ();
          auto hostExports = destMat->exports_.view_host();
          if (verbose) {
            std::ostringstream os;
            os << *verbosePrefix << "Calling 3-arg doPostsAndWaits"
               << std::endl;
            std::cerr << os.str ();
          }
          Distor.doPostsAndWaits (hostExports,
                                  constantNumPackets,
                                  hostImports);
          if (verbose) {
            std::ostringstream os;
            os << *verbosePrefix << "Finished 3-arg doPostsAndWaits"
               << std::endl;
            std::cerr << os.str ();
          }
        }
      }
    }
    }

    /*********************************************************************/
    /**** 3) Copy all of the Same/Permute/Remote data into CSR_arrays ****/
    /*********************************************************************/

    bool runOnHost = std::is_same_v<typename device_type::memory_space, Kokkos::HostSpace> && !useKokkosPath;

    Teuchos::Array<int> RemotePids;
    if (runOnHost) {
      Teuchos::Array<int> TargetPids;
      // Backwards compatibility measure.  We'll use this again below.
  
      // TODO JHU Need to track down why numImportPacketsPerLID_ has not been corrently marked as modified on host (which it has been)
      // TODO JHU somewhere above, e.g., call to Distor.doPostsAndWaits().
      // TODO JHU This only becomes apparent as we begin to convert TAFC to run on device.
      destMat->numImportPacketsPerLID_.modify_host(); //FIXME
  
#  ifdef HAVE_TPETRA_MMM_TIMINGS
      RCP<TimeMonitor> tmCopySPRdata = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("TAFC unpack-count-resize + copy same-perm-remote data"))));
#  endif
      ArrayRCP<size_t> CSR_rowptr;
      ArrayRCP<GO> CSR_colind_GID;
      ArrayRCP<LO> CSR_colind_LID;
      ArrayRCP<Scalar> CSR_vals;
  
      destMat->imports_.sync_device ();
      destMat->numImportPacketsPerLID_.sync_device ();
  
      size_t N = BaseRowMap->getLocalNumElements ();

      auto RemoteLIDs_d = RemoteLIDs.view_device();
      auto PermuteToLIDs_d = PermuteToLIDs.view_device();
      auto PermuteFromLIDs_d = PermuteFromLIDs.view_device();

      Details::unpackAndCombineIntoCrsArrays(
                                     *this, 
                                     RemoteLIDs_d,
                                     destMat->imports_.view_device(),                //hostImports
                                     destMat->numImportPacketsPerLID_.view_device(), //numImportPacketsPerLID
                                     NumSameIDs,
                                     PermuteToLIDs_d,
                                     PermuteFromLIDs_d,
                                     N,
                                     MyPID,
                                     CSR_rowptr,
                                     CSR_colind_GID,
                                     CSR_vals,
                                     SourcePids(),
                                     TargetPids);
  
      // If LO and GO are the same, we can reuse memory when
      // converting the column indices from global to local indices.
      if (typeid (LO) == typeid (GO)) {
        CSR_colind_LID = Teuchos::arcp_reinterpret_cast<LO> (CSR_colind_GID);
      }
      else {
        CSR_colind_LID.resize (CSR_colind_GID.size());
      }
      CSR_colind_LID.resize (CSR_colind_GID.size());
  
      // On return from unpackAndCombineIntoCrsArrays TargetPids[i] == -1 for locally
      // owned entries.  Convert them to the actual PID.
      // JHU FIXME This can be done within unpackAndCombineIntoCrsArrays with a parallel_for.
      for(size_t i=0; i<static_cast<size_t>(TargetPids.size()); i++)
      {
        if(TargetPids[i] == -1) TargetPids[i] = MyPID;
      }
#ifdef HAVE_TPETRA_MMM_TIMINGS
      tmCopySPRdata = Teuchos::null;
#endif
      /**************************************************************/
      /**** 4) Call Optimized MakeColMap w/ no Directory Lookups ****/
      /**************************************************************/
      // Call an optimized version of makeColMap that avoids the
      // Directory lookups (since the Import object knows who owns all
      // the GIDs).
      if (verbose) {
        std::ostringstream os;
        os << *verbosePrefix << "Calling lowCommunicationMakeColMapAndReindex"
           << std::endl;
        std::cerr << os.str ();
      }
      {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC makeColMap")));
#endif
      Import_Util::lowCommunicationMakeColMapAndReindexSerial(CSR_rowptr (),
                                                        CSR_colind_LID (),
                                                        CSR_colind_GID (),
                                                        BaseDomainMap,
                                                        TargetPids,
                                                        RemotePids,
                                                        MyColMap);
      }

      if (verbose) {
        std::ostringstream os;
        os << *verbosePrefix << "restrictComm="
           << (restrictComm ? "true" : "false") << std::endl;
        std::cerr << os.str ();
      }
  
      /*******************************************************/
      /**** 4) Second communicator restriction phase      ****/
      /*******************************************************/
      {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC restrict colmap")));
#endif
      if (restrictComm) {
        ReducedColMap = (MyRowMap.getRawPtr () == MyColMap.getRawPtr ()) ?
          ReducedRowMap :
          MyColMap->replaceCommWithSubset (ReducedComm);
        MyColMap = ReducedColMap; // Reset the "my" maps
      }
  
      // Replace the col map
      if (verbose) {
        std::ostringstream os;
        os << *verbosePrefix << "Calling replaceColMap" << std::endl;
        std::cerr << os.str ();
      }
      destMat->replaceColMap (MyColMap);
  
      // Short circuit if the processor is no longer in the communicator
      //
      // NOTE: Epetra replaces modifies all "removed" processes so they
      // have a dummy (serial) Map that doesn't touch the original
      // communicator.  Duplicating that here might be a good idea.
      if (ReducedComm.is_null ()) {
        if (verbose) {
          std::ostringstream os;
          os << *verbosePrefix << "I am no longer in the communicator; "
            "returning" << std::endl;
          std::cerr << os.str ();
        }
        return;
      }
      }
  
      /***************************************************/
      /**** 5) Sort                                   ****/
      /***************************************************/
      if ((! reverseMode && xferAsImport != nullptr) ||
          (reverseMode && xferAsExport != nullptr)) {
        if (verbose) {
          std::ostringstream os;
          os << *verbosePrefix << "Calling sortCrsEntries" << endl;
          std::cerr << os.str ();
        }
#ifdef HAVE_TPETRA_MMM_TIMINGS
        Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC sortCrsEntries")));
#endif
        Import_Util::sortCrsEntries (CSR_rowptr(),
                                     CSR_colind_LID(),
                                     CSR_vals());
      }
      else if ((! reverseMode && xferAsExport != nullptr) ||
               (reverseMode && xferAsImport != nullptr)) {
        if (verbose) {
          std::ostringstream os;
          os << *verbosePrefix << "Calling sortAndMergeCrsEntries"
             << endl;
          std::cerr << os.str();
        }
#ifdef HAVE_TPETRA_MMM_TIMINGS
        Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC sortAndMergeCrsEntries")));
#endif
        Import_Util::sortAndMergeCrsEntries (CSR_rowptr(),
                                             CSR_colind_LID(),
                                             CSR_vals());
        if (CSR_rowptr[N] != static_cast<size_t>(CSR_vals.size())) {
          CSR_colind_LID.resize (CSR_rowptr[N]);
          CSR_vals.resize (CSR_rowptr[N]);
        }
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error, "Tpetra::CrsMatrix::"
          "transferAndFillComplete: Should never get here!  "
          "Please report this bug to a Tpetra developer.");
      }
      /***************************************************/
      /**** 6) Reset the colmap and the arrays        ****/
      /***************************************************/
  
      if (verbose) {
        std::ostringstream os;
        os << *verbosePrefix << "Calling destMat->setAllValues" << endl;
        std::cerr << os.str ();
      }
  
      // Call constructor for the new matrix (restricted as needed)
      //
      // NOTE (mfh 15 May 2014) This should work fine for the Kokkos
      // refactor version of CrsMatrix, though it reserves the right to
      // make a deep copy of the arrays.
      {
#ifdef HAVE_TPETRA_MMM_TIMINGS
        Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC setAllValues")));
#endif
        destMat->setAllValues (CSR_rowptr, CSR_colind_LID, CSR_vals);
      }

    } else {
      // run on device
  
  
      // Backwards compatibility measure.  We'll use this again below.
  
      // TODO JHU Need to track down why numImportPacketsPerLID_ has not been corrently marked as modified on host (which it has been)
      // TODO JHU somewhere above, e.g., call to Distor.doPostsAndWaits().
      // TODO JHU This only becomes apparent as we begin to convert TAFC to run on device.
      destMat->numImportPacketsPerLID_.modify_host(); //FIXME
  
#  ifdef HAVE_TPETRA_MMM_TIMINGS
      RCP<TimeMonitor> tmCopySPRdata = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("TAFC unpack-count-resize + copy same-perm-remote data"))));
#  endif
      ArrayRCP<size_t> CSR_rowptr;
      ArrayRCP<GO> CSR_colind_GID;
      ArrayRCP<LO> CSR_colind_LID;
      ArrayRCP<Scalar> CSR_vals;
  
      destMat->imports_.sync_device ();
      destMat->numImportPacketsPerLID_.sync_device ();
  
      size_t N = BaseRowMap->getLocalNumElements ();
  
      auto RemoteLIDs_d = RemoteLIDs.view_device();
      auto PermuteToLIDs_d = PermuteToLIDs.view_device();
      auto PermuteFromLIDs_d = PermuteFromLIDs.view_device();
  
      Kokkos::View<size_t*,device_type> CSR_rowptr_d;
      Kokkos::View<GO*,device_type>     CSR_colind_GID_d;
      Kokkos::View<LO*,device_type>     CSR_colind_LID_d;
      Kokkos::View<impl_scalar_type*,device_type> CSR_vals_d;
      Kokkos::View<int*,device_type>    TargetPids_d;
  
      Details::unpackAndCombineIntoCrsArrays(
                                     *this, 
                                     RemoteLIDs_d,
                                     destMat->imports_.view_device(),                //hostImports
                                     destMat->numImportPacketsPerLID_.view_device(), //numImportPacketsPerLID
                                     NumSameIDs,
                                     PermuteToLIDs_d,
                                     PermuteFromLIDs_d,
                                     N,
                                     MyPID,
                                     CSR_rowptr_d,
                                     CSR_colind_GID_d,
                                     CSR_vals_d,
                                     SourcePids(),
                                     TargetPids_d);
  
      Kokkos::resize (CSR_colind_LID_d, CSR_colind_GID_d.size());
  
#ifdef HAVE_TPETRA_MMM_TIMINGS
      tmCopySPRdata = Teuchos::null;
#endif
      /**************************************************************/
      /**** 4) Call Optimized MakeColMap w/ no Directory Lookups ****/
      /**************************************************************/
      // Call an optimized version of makeColMap that avoids the
      // Directory lookups (since the Import object knows who owns all
      // the GIDs).
      if (verbose) {
        std::ostringstream os;
        os << *verbosePrefix << "Calling lowCommunicationMakeColMapAndReindex"
           << std::endl;
        std::cerr << os.str ();
      }
      {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC makeColMap")));
#endif
      Import_Util::lowCommunicationMakeColMapAndReindex(CSR_rowptr_d,
                                                        CSR_colind_LID_d,
                                                        CSR_colind_GID_d,
                                                        BaseDomainMap,
                                                        TargetPids_d,
                                                        RemotePids,
                                                        MyColMap);
      }
  
      if (verbose) {
        std::ostringstream os;
        os << *verbosePrefix << "restrictComm="
           << (restrictComm ? "true" : "false") << std::endl;
        std::cerr << os.str ();
      }
  
      /*******************************************************/
      /**** 4) Second communicator restriction phase      ****/
      /*******************************************************/
      {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC restrict colmap")));
#endif
      if (restrictComm) {
        ReducedColMap = (MyRowMap.getRawPtr () == MyColMap.getRawPtr ()) ?
          ReducedRowMap :
          MyColMap->replaceCommWithSubset (ReducedComm);
        MyColMap = ReducedColMap; // Reset the "my" maps
      }
  
      // Replace the col map
      if (verbose) {
        std::ostringstream os;
        os << *verbosePrefix << "Calling replaceColMap" << std::endl;
        std::cerr << os.str ();
      }
      destMat->replaceColMap (MyColMap);
  
      // Short circuit if the processor is no longer in the communicator
      //
      // NOTE: Epetra replaces modifies all "removed" processes so they
      // have a dummy (serial) Map that doesn't touch the original
      // communicator.  Duplicating that here might be a good idea.
      if (ReducedComm.is_null ()) {
        if (verbose) {
          std::ostringstream os;
          os << *verbosePrefix << "I am no longer in the communicator; "
            "returning" << std::endl;
          std::cerr << os.str ();
        }
        return;
      }
      }
  
      /***************************************************/
      /**** 5) Sort                                   ****/
      /***************************************************/

      if ((! reverseMode && xferAsImport != nullptr) ||
          (reverseMode && xferAsExport != nullptr)) {
        if (verbose) {
          std::ostringstream os;
          os << *verbosePrefix << "Calling sortCrsEntries" << endl;
          std::cerr << os.str ();
        }
#ifdef HAVE_TPETRA_MMM_TIMINGS
        Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC sortCrsEntries")));
#endif
        Import_Util::sortCrsEntries (CSR_rowptr_d,
                                     CSR_colind_LID_d,
                                     CSR_vals_d);
      }
      else if ((! reverseMode && xferAsExport != nullptr) ||
               (reverseMode && xferAsImport != nullptr)) {
        if (verbose) {
          std::ostringstream os;
          os << *verbosePrefix << "Calling sortAndMergeCrsEntries"
             << endl;
          std::cerr << os.str();
        }
#ifdef HAVE_TPETRA_MMM_TIMINGS
        Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC sortAndMergeCrsEntries")));
#endif
        Import_Util::sortAndMergeCrsEntries (CSR_rowptr_d,
                                             CSR_colind_LID_d,
                                             CSR_vals_d);
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error, "Tpetra::CrsMatrix::"
          "transferAndFillComplete: Should never get here!  "
          "Please report this bug to a Tpetra developer.");
      }

      /***************************************************/
      /**** 6) Reset the colmap and the arrays        ****/
      /***************************************************/
  
      if (verbose) {
        std::ostringstream os;
        os << *verbosePrefix << "Calling destMat->setAllValues" << endl;
        std::cerr << os.str ();
      }
  
      {
#ifdef HAVE_TPETRA_MMM_TIMINGS
        Teuchos::TimeMonitor MMrc(*TimeMonitor::getNewTimer(prefix + std::string("TAFC setAllValues")));
#endif
        destMat->setAllValues (CSR_rowptr_d, CSR_colind_LID_d, CSR_vals_d);
      }
  
    } //if (runOnHost) .. else ..

    /***************************************************/
    /**** 7) Build Importer & Call ESFC             ****/
    /***************************************************/
#ifdef HAVE_TPETRA_MMM_TIMINGS
    RCP<TimeMonitor> tmIESFC = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("TAFC build importer and esfc"))));
#endif
    // Pre-build the importer using the existing PIDs
    Teuchos::ParameterList esfc_params;

    RCP<import_type> MyImport;

    // Fulfull the non-blocking allreduce on reduced_mismatch.
    if (iallreduceRequest.get () != nullptr) {
      if (verbose) {
        std::ostringstream os;
        os << *verbosePrefix << "Calling iallreduceRequest->wait()"
           << endl;
        std::cerr << os.str ();
      }
      iallreduceRequest->wait ();
      if (reduced_mismatch != 0) {
        isMM = false;
      }
    }

    if( isMM ) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
        Teuchos::TimeMonitor MMisMM (*TimeMonitor::getNewTimer(prefix + std::string("isMM Block")));
#endif
        // Combine all type1/2/3 lists, [filter them], then call the expert import constructor.

        if (verbose) {
          std::ostringstream os;
          os << *verbosePrefix << "Getting CRS pointers" << endl;
          std::cerr << os.str ();
        }

        Teuchos::ArrayRCP<LocalOrdinal> type3LIDs;
        Teuchos::ArrayRCP<int>          type3PIDs;
        auto rowptr = getCrsGraph()->getLocalRowPtrsHost();
        auto colind = getCrsGraph()->getLocalIndicesHost();

        if (verbose) {
          std::ostringstream os;
          os << *verbosePrefix << "Calling reverseNeighborDiscovery" << std::endl;
          std::cerr << os.str ();
        }

        {
#ifdef HAVE_TPETRA_MMM_TIMINGS
            TimeMonitor tm_rnd (*TimeMonitor::getNewTimer(prefix + std::string("isMMrevNeighDis")));
#endif
            Import_Util::reverseNeighborDiscovery(*this,
                                                  rowptr,
                                                  colind,
                                                  rowTransfer,
                                                  MyImporter,
                                                  MyDomainMap,
                                                  type3PIDs,
                                                  type3LIDs,
                                                  ReducedComm);
        }

        if (verbose) {
          std::ostringstream os;
          os << *verbosePrefix << "Done with reverseNeighborDiscovery" << std::endl;
          std::cerr << os.str ();
        }

        Teuchos::ArrayView<const int>  EPID1 = MyImporter.is_null() ? Teuchos::ArrayView<const int>() : MyImporter->getExportPIDs();
        Teuchos::ArrayView<const LO>   ELID1 = MyImporter.is_null() ? Teuchos::ArrayView<const LO>() : MyImporter->getExportLIDs();

        Teuchos::ArrayView<const int>  TEPID2  =  rowTransfer.getExportPIDs(); // row matrix
        Teuchos::ArrayView<const LO>   TELID2  =  rowTransfer.getExportLIDs();

        const int numCols = getGraph()->getColMap()->getLocalNumElements(); // may be dup
        // from EpetraExt_MMHelpers.cpp: build_type2_exports
        std::vector<bool> IsOwned(numCols,true);
        std::vector<int>  SentTo(numCols,-1);
        if (! MyImporter.is_null ()) {
          for (auto && rlid : MyImporter->getRemoteLIDs()) { // the remoteLIDs must be from sourcematrix
            IsOwned[rlid]=false;
          }
        }

        std::vector<std::pair<int,GO> > usrtg;
        usrtg.reserve(TEPID2.size());

        {
          const auto& colMap = * (this->getColMap ()); // *this is sourcematrix
          for (Array_size_type i = 0; i < TEPID2.size (); ++i) {
            const LO  row = TELID2[i];
            const int pid = TEPID2[i];
            for (auto j = rowptr[row]; j < rowptr[row+1]; ++j) {
              const int col = colind[j];
              if (IsOwned[col] && SentTo[col] != pid) {
                SentTo[col]    = pid;
                GO gid = colMap.getGlobalElement (col);
                usrtg.push_back (std::pair<int,GO> (pid, gid));
              }
            }
          }
        }

// This sort can _not_ be omitted.[
        std::sort(usrtg.begin(),usrtg.end()); // default comparator does the right thing, now sorted in gid order
        auto eopg = std ::unique(usrtg.begin(),usrtg.end());
        // 25 Jul 2018: Could just ignore the entries at and after eopg.
        usrtg.erase(eopg,usrtg.end());

        const Array_size_type type2_us_size = usrtg.size();
        Teuchos::ArrayRCP<int>  EPID2=Teuchos::arcp(new int[type2_us_size],0,type2_us_size,true);
        Teuchos::ArrayRCP< LO>  ELID2=Teuchos::arcp(new  LO[type2_us_size],0,type2_us_size,true);

        int pos=0;
        for(auto && p : usrtg) {
            EPID2[pos]= p.first;
            ELID2[pos]= this->getDomainMap()->getLocalElement(p.second);
            pos++;
        }

        Teuchos::ArrayView<int>  EPID3  = type3PIDs();
        Teuchos::ArrayView< LO>  ELID3  = type3LIDs();
        GO InfGID = std::numeric_limits<GO>::max();
        int InfPID = INT_MAX;
#ifdef TPETRA_MIN3
#  undef TPETRA_MIN3
#endif // TPETRA_MIN3
#define TPETRA_MIN3(x,y,z) ((x)<(y)?(std::min(x,z)):(std::min(y,z)))
        int i1=0, i2=0, i3=0;
        int Len1 = EPID1.size();
        int Len2 = EPID2.size();
        int Len3 = EPID3.size();

        int MyLen=Len1+Len2+Len3;
        Teuchos::ArrayRCP<LO>  userExportLIDs = Teuchos::arcp(new LO[MyLen],0,MyLen,true);
        Teuchos::ArrayRCP<int> userExportPIDs = Teuchos::arcp(new int[MyLen],0,MyLen,true);
        int iloc = 0; // will be the size of the userExportLID/PIDs

        while(i1 < Len1 || i2 < Len2 || i3 < Len3){
            int PID1 = (i1<Len1)?(EPID1[i1]):InfPID;
            int PID2 = (i2<Len2)?(EPID2[i2]):InfPID;
            int PID3 = (i3<Len3)?(EPID3[i3]):InfPID;

            GO GID1 = (i1<Len1)?getDomainMap()->getGlobalElement(ELID1[i1]):InfGID;
            GO GID2 = (i2<Len2)?getDomainMap()->getGlobalElement(ELID2[i2]):InfGID;
            GO GID3 = (i3<Len3)?getDomainMap()->getGlobalElement(ELID3[i3]):InfGID;

            int MIN_PID = TPETRA_MIN3(PID1,PID2,PID3);
            GO  MIN_GID = TPETRA_MIN3( ((PID1==MIN_PID)?GID1:InfGID), ((PID2==MIN_PID)?GID2:InfGID), ((PID3==MIN_PID)?GID3:InfGID));
#ifdef TPETRA_MIN3
#  undef TPETRA_MIN3
#endif // TPETRA_MIN3
            bool added_entry=false;

            if(PID1 == MIN_PID && GID1 == MIN_GID){
                userExportLIDs[iloc]=ELID1[i1];
                userExportPIDs[iloc]=EPID1[i1];
                i1++;
                added_entry=true;
                iloc++;
            }
            if(PID2 == MIN_PID && GID2 == MIN_GID){
                if(!added_entry) {
                    userExportLIDs[iloc]=ELID2[i2];
                    userExportPIDs[iloc]=EPID2[i2];
                    added_entry=true;
                    iloc++;
                }
                i2++;
            }
            if(PID3 == MIN_PID && GID3 == MIN_GID){
                if(!added_entry) {
                    userExportLIDs[iloc]=ELID3[i3];
                    userExportPIDs[iloc]=EPID3[i3];
                    iloc++;
                }
                i3++;
            }
        }

        if (verbose) {
          std::ostringstream os;
          os << *verbosePrefix << "Create Import" << std::endl;
          std::cerr << os.str ();
        }

#ifdef HAVE_TPETRA_MMM_TIMINGS
        auto ismmIctor(*TimeMonitor::getNewTimer(prefix + std::string("isMMIportCtor")));
#endif
        Teuchos::RCP<Teuchos::ParameterList> plist = rcp(new Teuchos::ParameterList());
        // 25 Jul 2018: Test for equality with the non-isMM path's Import object.
        if ((MyDomainMap != MyColMap) && (!MyDomainMap->isSameAs(*MyColMap)))
          MyImport = rcp ( new import_type (MyDomainMap,
                                            MyColMap,
                                            RemotePids,
                                            userExportLIDs.view(0,iloc).getConst(),
                                            userExportPIDs.view(0,iloc).getConst(),
                                            plist)
            );

        if (verbose) {
          std::ostringstream os;
          os << *verbosePrefix << "Call expertStaticFillComplete" << std::endl;
          std::cerr << os.str ();
        }

        {
#ifdef HAVE_TPETRA_MMM_TIMINGS
            TimeMonitor esfc (*TimeMonitor::getNewTimer(prefix + std::string("isMM::destMat->eSFC")));
            esfc_params.set("Timer Label",label+std::string("isMM eSFC"));
#endif
            if(!params.is_null())
                esfc_params.set("compute global constants",params->get("compute global constants",true));
            destMat->expertStaticFillComplete (MyDomainMap, MyRangeMap, MyImport,Teuchos::null,rcp(new Teuchos::ParameterList(esfc_params)));

        }

    }  // if(isMM)
    else {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      TimeMonitor MMnotMMblock (*TimeMonitor::getNewTimer(prefix + std::string("TAFC notMMblock")));
#endif
      if (verbose) {
        std::ostringstream os;
        os << *verbosePrefix << "Create Import" << std::endl;
        std::cerr << os.str ();
      }

#ifdef HAVE_TPETRA_MMM_TIMINGS
      TimeMonitor  notMMIcTor(*TimeMonitor::getNewTimer(prefix + std::string("TAFC notMMCreateImporter")));
#endif
      Teuchos::RCP<Teuchos::ParameterList> mypars = rcp(new Teuchos::ParameterList);
      mypars->set("Timer Label","notMMFrom_tAFC");
      if ((MyDomainMap != MyColMap) && (!MyDomainMap->isSameAs(*MyColMap)))
        MyImport = rcp (new import_type (MyDomainMap, MyColMap, RemotePids, mypars));

      if (verbose) {
        std::ostringstream os;
        os << *verbosePrefix << "Call expertStaticFillComplete" << endl;
        std::cerr << os.str ();
      }

#ifdef HAVE_TPETRA_MMM_TIMINGS
      TimeMonitor  esfcnotmm(*TimeMonitor::getNewTimer(prefix + std::string("notMMdestMat->expertStaticFillComplete")));
      esfc_params.set("Timer Label",prefix+std::string("notMM eSFC"));
#else
      esfc_params.set("Timer Label",std::string("notMM eSFC"));
#endif

      if (!params.is_null ()) {
        esfc_params.set ("compute global constants",
                         params->get ("compute global constants", true));
      }
      destMat->expertStaticFillComplete (MyDomainMap, MyRangeMap,
                                         MyImport, Teuchos::null,
                                         rcp (new Teuchos::ParameterList (esfc_params)));
    }

#ifdef HAVE_TPETRA_MMM_TIMINGS
    tmIESFC = Teuchos::null;
#endif

    if (verbose) {
      std::ostringstream os;
      os << *verbosePrefix << "Done" << endl;
      std::cerr << os.str ();
    }
  } //transferAndFillComplete


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  importAndFillComplete (Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& destMatrix,
                         const import_type& importer,
                         const Teuchos::RCP<const map_type>& domainMap,
                         const Teuchos::RCP<const map_type>& rangeMap,
                         const Teuchos::RCP<Teuchos::ParameterList>& params) const
  {
    transferAndFillComplete (destMatrix, importer, Teuchos::null, domainMap, rangeMap, params);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  importAndFillComplete (Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& destMatrix,
                         const import_type& rowImporter,
                         const import_type& domainImporter,
                         const Teuchos::RCP<const map_type>& domainMap,
                         const Teuchos::RCP<const map_type>& rangeMap,
                         const Teuchos::RCP<Teuchos::ParameterList>& params) const
  {
    transferAndFillComplete (destMatrix, rowImporter, Teuchos::rcpFromRef(domainImporter), domainMap, rangeMap, params);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  exportAndFillComplete (Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& destMatrix,
                         const export_type& exporter,
                         const Teuchos::RCP<const map_type>& domainMap,
                         const Teuchos::RCP<const map_type>& rangeMap,
                         const Teuchos::RCP<Teuchos::ParameterList>& params) const
  {
    transferAndFillComplete (destMatrix, exporter, Teuchos::null, domainMap, rangeMap, params);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  exportAndFillComplete (Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& destMatrix,
                         const export_type& rowExporter,
                         const export_type& domainExporter,
                         const Teuchos::RCP<const map_type>& domainMap,
                         const Teuchos::RCP<const map_type>& rangeMap,
                         const Teuchos::RCP<Teuchos::ParameterList>& params) const
  {
    transferAndFillComplete (destMatrix, rowExporter, Teuchos::rcpFromRef(domainExporter), domainMap, rangeMap, params);
  }

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_CRSMATRIX_MATRIX_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class CrsMatrix< SCALAR , LO , GO , NODE >;

#define TPETRA_CRSMATRIX_CONVERT_INSTANT(SO,SI,LO,GO,NODE) \
  \
  template Teuchos::RCP< CrsMatrix< SO , LO , GO , NODE > >   \
                CrsMatrix< SI , LO , GO , NODE >::convert< SO > () const;

#define TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(SCALAR, LO, GO, NODE) \
  template<>                                                                        \
  Teuchos::RCP<CrsMatrix<SCALAR, LO, GO, NODE> >                        \
  importAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrix<SCALAR, LO, GO, NODE> >& sourceMatrix, \
                                  const Import<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,  \
                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,  \
                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type>& importer, \
                                  const Teuchos::RCP<const Map<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,      \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,     \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type> >& domainMap, \
                                  const Teuchos::RCP<const Map<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,      \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,     \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type> >& rangeMap,  \
                                                               const Teuchos::RCP<Teuchos::ParameterList>& params);

#define TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT_TWO(SCALAR, LO, GO, NODE) \
  template<>                                                                        \
  Teuchos::RCP<CrsMatrix<SCALAR, LO, GO, NODE> >                        \
  importAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrix<SCALAR, LO, GO, NODE> >& sourceMatrix, \
                                  const Import<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,  \
                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,  \
                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type>& rowImporter, \
                                  const Import<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,  \
                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,  \
                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type>& domainImporter, \
                                  const Teuchos::RCP<const Map<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,      \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,     \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type> >& domainMap, \
                                  const Teuchos::RCP<const Map<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,      \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,     \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type> >& rangeMap,  \
                                                               const Teuchos::RCP<Teuchos::ParameterList>& params);


#define TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(SCALAR, LO, GO, NODE) \
  template<>                                                                        \
  Teuchos::RCP<CrsMatrix<SCALAR, LO, GO, NODE> >                        \
  exportAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrix<SCALAR, LO, GO, NODE> >& sourceMatrix, \
                                  const Export<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,  \
                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,  \
                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type>& exporter, \
                                  const Teuchos::RCP<const Map<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,      \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,     \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type> >& domainMap, \
                                  const Teuchos::RCP<const Map<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,      \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,     \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type> >& rangeMap,  \
                                                               const Teuchos::RCP<Teuchos::ParameterList>& params);

#define TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT_TWO(SCALAR, LO, GO, NODE) \
  template<>                                                                        \
  Teuchos::RCP<CrsMatrix<SCALAR, LO, GO, NODE> >                        \
  exportAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrix<SCALAR, LO, GO, NODE> >& sourceMatrix, \
                                  const Export<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,  \
                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,  \
                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type>& rowExporter, \
                                  const Export<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,  \
                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,  \
                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type>& domainExporter, \
                                  const Teuchos::RCP<const Map<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,      \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,     \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type> >& domainMap, \
                                  const Teuchos::RCP<const Map<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,      \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,     \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type> >& rangeMap,  \
                                                               const Teuchos::RCP<Teuchos::ParameterList>& params);


#define TPETRA_CRSMATRIX_INSTANT(SCALAR, LO, GO ,NODE)                    \
  TPETRA_CRSMATRIX_MATRIX_INSTANT(SCALAR, LO, GO, NODE)                   \
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(SCALAR, LO, GO, NODE) \
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(SCALAR, LO, GO, NODE) \
  TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT_TWO(SCALAR, LO, GO, NODE) \
  TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT_TWO(SCALAR, LO, GO, NODE)

#endif // TPETRA_CRSMATRIX_DEF_HPP

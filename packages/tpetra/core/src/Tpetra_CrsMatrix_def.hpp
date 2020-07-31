// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
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
#include "KokkosSparse_getDiagCopy.hpp"
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

#include <memory>
#include <sstream>
#include <typeinfo>
#include <utility>
#include <vector>

using Teuchos::rcpFromRef;

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
             const ProfileType pftype,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
  {
    const char tfecfFuncName[] = "CrsMatrix(RCP<const Map>, size_t, "
      "ProfileType[, RCP<ParameterList>]): ";
    Teuchos::RCP<crs_graph_type> graph;
    try {
      graph = Teuchos::rcp (new crs_graph_type (rowMap, maxNumEntriesPerRow,
                                                pftype, params));
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "CrsGraph constructor (RCP<const Map>, "
         "size_t, ProfileType[, RCP<ParameterList>]) threw an exception: "
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
             const ProfileType pftype,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
  {
    const char tfecfFuncName[] = "CrsMatrix(RCP<const Map>, "
      "ArrayView<const size_t>, ProfileType[, RCP<ParameterList>]): ";
    Teuchos::RCP<crs_graph_type> graph;
    try {
      using Teuchos::rcp;
      graph = rcp(new crs_graph_type(rowMap, numEntPerRowToAlloc,
                                     pftype, params));
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "CrsGraph constructor "
         "(RCP<const Map>, ArrayView<const size_t>, "
         "ProfileType[, RCP<ParameterList>]) threw an exception: "
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
             const ProfileType pftype,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
  {
    const char tfecfFuncName[] = "CrsMatrix(RCP<const Map>, "
      "RCP<const Map>, size_t, ProfileType[, RCP<ParameterList>]): ";
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
                                                pftype, params));
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "CrsGraph constructor (RCP<const Map>, "
         "RCP<const Map>, size_t, ProfileType[, RCP<ParameterList>]) threw an "
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
             const ProfileType pftype,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
  {
    const char tfecfFuncName[] =
      "CrsMatrix(RCP<const Map>, RCP<const Map>, "
      "ArrayView<const size_t>, ProfileType[, RCP<ParameterList>]): ";
    Teuchos::RCP<crs_graph_type> graph;
    try {
      graph = Teuchos::rcp (new crs_graph_type (rowMap, colMap,
                                                numEntPerRowToAlloc,
                                                pftype, params));
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "CrsGraph constructor (RCP<const Map>, "
         "RCP<const Map>, ArrayView<const size_t>, ProfileType[, "
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
    typedef typename local_matrix_type::values_type values_type;
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

    const size_t numCols = graph->getColMap ()->getNodeNumElements ();
    auto lclGraph = graph->getLocalGraph ();
    const size_t numEnt = lclGraph.entries.extent (0);
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Allocate values: " << numEnt << endl;
      std::cerr << os.str ();
    }
    values_type val ("Tpetra::CrsMatrix::val", numEnt);

    auto lclMat = std::make_shared<local_matrix_type>
      ("Tpetra::CrsMatrix::lclMatrix_", numCols, val, lclGraph);
    lclMatrix_ = std::make_shared<local_multiply_op_type> (lclMat);

    // FIXME (22 Jun 2016) I would very much like to get rid of
    // k_values1D_ at some point.  I find it confusing to have all
    // these extra references lying around.
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Assign k_values1D_: old="
         << k_values1D_.extent(0) << ", new="
         << lclMat->values.extent(0) << endl;
      std::cerr << os.str ();
    }
    k_values1D_ = lclMat->values;

    checkInternalState ();

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str ();
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  CrsMatrix (const Teuchos::RCP<const crs_graph_type>& graph,
             const typename local_matrix_type::values_type& values,
             const Teuchos::RCP<Teuchos::ParameterList>& /* params */) :
    dist_object_type (graph->getRowMap ()),
    staticGraph_ (graph),
    storageStatus_ (Details::STORAGE_1D_PACKED)
  {
    const char tfecfFuncName[] = "CrsMatrix(RCP<const CrsGraph>, "
      "local_matrix_type::values_type, "
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

    const size_t numCols = graph->getColMap ()->getNodeNumElements ();
    auto lclGraph = graph->getLocalGraph ();

    auto lclMat = std::make_shared<local_matrix_type>
      ("Tpetra::CrsMatrix::lclMatrix_", numCols, values, lclGraph);
    lclMatrix_ = std::make_shared<local_multiply_op_type> (lclMat);

    // FIXME (22 Jun 2016) I would very much like to get rid of
    // k_values1D_ at some point.  I find it confusing to have all
    // these extra references lying around.
    k_values1D_ = lclMat->values;

    checkInternalState ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
             const Teuchos::RCP<const map_type>& colMap,
             const typename local_matrix_type::row_map_type& rowPointers,
             const typename local_graph_type::entries_type::non_const_type& columnIndices,
             const typename local_matrix_type::values_type& values,
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
    auto lclGraph = graph->getLocalGraph ();
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

    const size_t numCols = graph->getColMap ()->getNodeNumElements ();

    auto lclMat = std::make_shared<local_matrix_type>
      ("Tpetra::CrsMatrix::lclMatrix_", numCols, values, lclGraph);
    lclMatrix_ = std::make_shared<local_multiply_op_type> (lclMat);

    auto newValues = lclMat->values;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (newValues.extent (0) != values.extent (0),
       std::logic_error, "Local matrix's constructor did not set the "
       "values correctly.  newValues.extent(0) = " <<
       newValues.extent (0) << " != values.extent(0) = " <<
       values.extent (0) << suffix);

    // FIXME (22 Jun 2016) I would very much like to get rid of
    // k_values1D_ at some point.  I find it confusing to have all
    // these extra references lying around.
    this->k_values1D_ = newValues;

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
    using values_type = typename local_matrix_type::values_type;
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
    auto lclGraph = staticGraph_->getLocalGraph ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (size_t (lclGraph.row_map.extent (0)) != size_t (ptr.size ()) ||
       size_t (lclGraph.entries.extent (0)) != size_t (ind.size ()),
       std::logic_error, "CrsGraph's constructor (rowMap, colMap, "
       "ptr, ind[, params]) did not set the local graph correctly.  "
       "Please report this bug to the Tpetra developers.");

    const size_t numCols =
      staticGraph_->getColMap ()->getNodeNumElements ();
    values_type valIn =
      getKokkosViewDeepCopy<device_type> (av_reinterpret_cast<IST> (val ()));

    auto lclMat = std::make_shared<local_matrix_type>
      ("Tpetra::CrsMatrix::lclMatrix_", numCols, valIn, lclGraph);
    lclMatrix_ = std::make_shared<local_multiply_op_type> (lclMat);

    // FIXME (22 Jun 2016) I would very much like to get rid of
    // k_values1D_ at some point.  I find it confusing to have all
    // these extra references lying around.
    this->k_values1D_ = lclMat->values;

    checkInternalState ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
             const Teuchos::RCP<const map_type>& colMap,
             const local_matrix_type& lclMatrix,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap),
    lclMatrix_ (std::make_shared<local_multiply_op_type>
                (std::make_shared<local_matrix_type> (lclMatrix))),
    k_values1D_ (lclMatrix.values),
    storageStatus_ (Details::STORAGE_1D_PACKED),
    fillComplete_ (true)
  {
    const char tfecfFuncName[] = "Tpetra::CrsMatrix(RCP<const Map>, "
      "RCP<const Map>, local_matrix_type[, RCP<ParameterList>]): ";
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
         "RCP<const Map>, local_graph_type[, RCP<ParameterList>]) threw an "
         "exception: " << e.what ());
    }
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (!graph->isFillComplete (), std::logic_error, "CrsGraph constructor (RCP"
       "<const Map>, RCP<const Map>, local_graph_type[, RCP<ParameterList>]) "
       "did not produce a fill-complete graph.  Please report this bug to the "
       "Tpetra developers.");
    // myGraph_ not null means that the matrix owns the graph.  This
    // is true because the column indices come in as nonconst through
    // the matrix, implying shared ownership.
    myGraph_ = graph;
    staticGraph_ = graph;

    const bool callComputeGlobalConstants = params.get () == nullptr ||
      params->get ("compute global constants", true);
    if (callComputeGlobalConstants) {
      this->computeGlobalConstants ();
    }

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
  CrsMatrix (const local_matrix_type& lclMatrix,
             const Teuchos::RCP<const map_type>& rowMap,
             const Teuchos::RCP<const map_type>& colMap,
             const Teuchos::RCP<const map_type>& domainMap,
             const Teuchos::RCP<const map_type>& rangeMap,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap),
    lclMatrix_ (std::make_shared<local_multiply_op_type>
                (std::make_shared<local_matrix_type> (lclMatrix))),
    k_values1D_ (lclMatrix.values),
    storageStatus_ (Details::STORAGE_1D_PACKED),
    fillComplete_ (true)
  {
    const char tfecfFuncName[] = "Tpetra::CrsMatrix(RCP<const Map>, "
      "RCP<const Map>, RCP<const Map>, RCP<const Map>, "
      "local_matrix_type[, RCP<ParameterList>]): ";
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
         "RCP<const Map>, RCP<const Map>, RCP<const Map>, local_graph_type[, "
         "RCP<ParameterList>]) threw an exception: " << e.what ());
    }
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! graph->isFillComplete (), std::logic_error, "CrsGraph "
       "constructor (RCP<const Map>, RCP<const Map>, RCP<const Map>, "
       "RCP<const Map>, local_graph_type[, RCP<ParameterList>]) did "
       "not produce a fillComplete graph." << suffix);
    // myGraph_ not null means that the matrix owns the graph.  This
    // is true because the column indices come in as nonconst through
    // the matrix, implying shared ownership.
    myGraph_ = graph;
    staticGraph_ = graph;

    const bool callComputeGlobalConstants = params.get () == nullptr ||
      params->get ("compute global constants", true);
    if (callComputeGlobalConstants) {
      this->computeGlobalConstants ();
    }

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
  CrsMatrix (const local_matrix_type& lclMatrix,
             const Teuchos::RCP<const map_type>& rowMap,
             const Teuchos::RCP<const map_type>& colMap,
             const Teuchos::RCP<const map_type>& domainMap,
             const Teuchos::RCP<const map_type>& rangeMap,
             const Teuchos::RCP<const import_type>& importer,
             const Teuchos::RCP<const export_type>& exporter,
             const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap),
    lclMatrix_ (std::make_shared<local_multiply_op_type>
                (std::make_shared<local_matrix_type> (lclMatrix))),
    k_values1D_ (lclMatrix.values),
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
         "(local_graph_type, Map, Map, Map, Map, Import, Export, "
         "params) threw: " << e.what ());
    }
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (!graph->isFillComplete (), std::logic_error, "CrsGraph "
       "constructor (local_graph_type, Map, Map, Map, Map, Import, "
       "Export, params) did not produce a fill-complete graph.  "
       "Please report this bug to the Tpetra developers.");
    // myGraph_ not null means that the matrix owns the graph.  This
    // is true because the column indices come in as nonconst through
    // the matrix, implying shared ownership.
    myGraph_ = graph;
    staticGraph_ = graph;

    const bool callComputeGlobalConstants = params.get () == nullptr ||
      params->get ("compute global constants", true);
    if (callComputeGlobalConstants) {
      this->computeGlobalConstants ();
    }

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
             const Teuchos::DataAccess copyOrView)
    : CrsMatrix (source.getCrsGraph (), source.getLocalValuesView ())
  {
    const char tfecfFuncName[] = "Tpetra::CrsMatrix("
      "const CrsMatrix&, const Teuchos::DataAccess): ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! source.isFillComplete (), std::invalid_argument,
       "Source graph must be fillComplete().");

    if (copyOrView == Teuchos::Copy) {
      using values_type = typename local_matrix_type::values_type;
      values_type vals = source.getLocalValuesView ();
      using Kokkos::view_alloc;
      using Kokkos::WithoutInitializing;
      values_type newvals (view_alloc ("val", WithoutInitializing),
                           vals.extent (0));
      Kokkos::deep_copy (newvals, vals);
      k_values1D_ = newvals;
      if (source.isFillComplete ()) {
        fillComplete (source.getDomainMap (), source.getRangeMap ());
      }
    }
    else if (copyOrView == Teuchos::View) {
      return;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::invalid_argument, "Second argument 'copyOrView' "
         "has an invalid value " << copyOrView << ".  Valid values "
         "include Teuchos::Copy = " << Teuchos::Copy << " and "
         "Teuchos::View = " << Teuchos::View << ".");
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  swap(CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & crs_matrix)
  {
    std::swap(crs_matrix.importMV_,      this->importMV_);        // mutable Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
    std::swap(crs_matrix.exportMV_,      this->exportMV_);        // mutable Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
    std::swap(crs_matrix.staticGraph_,   this->staticGraph_);     // Teuchos::RCP<const CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>
    std::swap(crs_matrix.myGraph_,       this->myGraph_);         // Teuchos::RCP<      CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>
    std::swap(crs_matrix.lclMatrix_,     this->lclMatrix_);       // KokkosSparse::CrsMatrix<impl_scalar_type, LocalOrdinal, execution_space, void, typename local_graph_type::size_type>
    std::swap(crs_matrix.k_values1D_,    this->k_values1D_);      // KokkosSparse::CrsMatrix<impl_scalar_type, LocalOrdinal, execution_space, void, typename local_graph_type::size_type>::values_type
    std::swap(crs_matrix.storageStatus_, this->storageStatus_);   // ::Tpetra::Details::EStorageStatus (enum f/m Tpetra_CrsGraph_decl.hpp)
    std::swap(crs_matrix.fillComplete_,  this->fillComplete_);    // bool
    std::swap(crs_matrix.nonlocals_,     this->nonlocals_);       // std::map<GO, pair<Teuchos::Array<GO>,Teuchos::Array<Scalar>>
    std::swap(crs_matrix.frobNorm_,      this->frobNorm_);        // mutable Kokkos::Details::ArithTraits<impl_scalar_type>::mag_type
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Teuchos::Comm<int> >
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getComm () const {
    return getCrsGraphRef ().getComm ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  ProfileType
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getProfileType () const {
    return this->getCrsGraphRef ().getProfileType ();
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
  getNodeNumEntries () const {
    return getCrsGraphRef ().getNodeNumEntries ();
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
  getNodeNumRows () const {
    return getCrsGraphRef ().getNodeNumRows ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getNodeNumCols () const {
    return getCrsGraphRef ().getNodeNumCols ();
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
  getNodeMaxNumRowEntries () const {
    return getCrsGraphRef ().getNodeMaxNumRowEntries ();
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
  typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalMatrix () const
  {
    return lclMatrix_.get () == nullptr ?
      local_matrix_type () :
      lclMatrix_->getLocalMatrix ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::shared_ptr<typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_multiply_op_type>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalMultiplyOperator () const
  {
    return lclMatrix_;
  }

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
    // "Static profile" means that the number of matrix entries in
    // each row was fixed at the time the CrsMatrix constructor was
    // called.  This lets us use 1-D storage for the matrix's
    // values.  ("1-D storage" means the same as that used by the
    // three arrays in the compressed sparse row storage format.)

    if (debug) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->staticGraph_.is_null (), std::logic_error,
         "this->getProfileType() == StaticProfile, but staticGraph_ "
         "is null." << suffix);
    }

    const size_t lclNumRows = this->staticGraph_->getNodeNumRows ();
    typename Graph::local_graph_type::row_map_type k_ptrs =
      this->staticGraph_->k_rowPtrs_;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (k_ptrs.extent (0) != lclNumRows+1, std::logic_error,
      "With StaticProfile, row offsets array has length "
      << k_ptrs.extent (0) << " != (lclNumRows+1) = "
      << (lclNumRows+1) << ".");

    const size_t lclTotalNumEntries =
      ::Tpetra::Details::getEntryOnHost (k_ptrs, lclNumRows);

    // Allocate array of (packed???) matrix values.
    using values_type = typename local_matrix_type::values_type;
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Allocate k_values1D_: Pre "
         << k_values1D_.extent(0) << ", post "
         << lclTotalNumEntries << endl;
      std::cerr << os.str();
    }
    this->k_values1D_ =
      values_type ("Tpetra::CrsMatrix::val", lclTotalNumEntries);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getAllValues (Teuchos::ArrayRCP<const size_t>& rowPointers,
                Teuchos::ArrayRCP<const LocalOrdinal>& columnIndices,
                Teuchos::ArrayRCP<const Scalar>& values) const
  {
    using Teuchos::RCP;
    const char tfecfFuncName[] = "getAllValues: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      columnIndices.size () != values.size (), std::runtime_error,
      "Requires that columnIndices and values are the same size.");

    RCP<const crs_graph_type> relevantGraph = getCrsGraph ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      relevantGraph.is_null (), std::runtime_error,
      "Requires that getCrsGraph() is not null.");
    try {
      rowPointers = relevantGraph->getNodeRowPtrs ();
    }
    catch (std::exception &e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        true, std::runtime_error,
        "Caught exception while calling graph->getNodeRowPtrs(): "
        << e.what ());
    }
    try {
      columnIndices = relevantGraph->getNodePackedIndices ();
    }
    catch (std::exception &e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        true, std::runtime_error,
        "Caught exception while calling graph->getNodePackedIndices(): "
        << e.what ());
    }
    Teuchos::ArrayRCP<const impl_scalar_type> vals =
      Kokkos::Compat::persistingView (k_values1D_);
    values = Teuchos::arcp_reinterpret_cast<const Scalar> (vals);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  fillLocalGraphAndMatrix (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    using ::Tpetra::Details::computeOffsetsFromCounts;
    using ::Tpetra::Details::getEntryOnHost;
    using Kokkos::create_mirror_view;
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    using row_map_type = typename local_matrix_type::row_map_type;
    using lclinds_1d_type = typename Graph::local_graph_type::entries_type::non_const_type;
    using values_type = typename local_matrix_type::values_type;
    Details::ProfilingRegion regionFLGAM
      ("Tpetra::CrsGraph::fillLocalGraphAndMatrix");

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

    const size_t lclNumRows = this->getNodeNumRows ();

    // This method's goal is to fill in the three arrays (compressed
    // sparse row format) that define the sparse graph's and matrix's
    // structure, and the sparse matrix's values.
    //
    // Use the nonconst version of row_map_type for k_ptrs,
    // because row_map_type is const and we need to modify k_ptrs here.
    typename row_map_type::non_const_type k_ptrs;
    row_map_type k_ptrs_const;
    lclinds_1d_type k_inds;
    values_type k_vals;

    // Get references to the data in myGraph_, so we can modify them
    // as well.  Note that we only call fillLocalGraphAndMatrix() if
    // the matrix owns the graph, which means myGraph_ is not null.
    lclinds_1d_type k_lclInds1D_ = myGraph_->k_lclInds1D_;

    typedef decltype (myGraph_->k_numRowEntries_) row_entries_type;

    // StaticProfile means that the matrix's column indices and
    // values are currently stored in a 1-D format, with row offsets
    // in k_rowPtrs_ and local column indices in k_lclInds1D_.

    // StaticProfile also means that the graph's array of row
    // offsets must already be allocated.
    typename Graph::local_graph_type::row_map_type curRowOffsets =
      myGraph_->k_rowPtrs_;

    if (debug) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (curRowOffsets.extent (0) == 0, std::logic_error,
         "(StaticProfile branch) curRowOffsets.extent(0) == 0.");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (curRowOffsets.extent (0) != lclNumRows + 1, std::logic_error,
         "(StaticProfile branch) curRowOffsets.extent(0) = "
         << curRowOffsets.extent (0) << " != lclNumRows + 1 = "
         << (lclNumRows + 1) << ".");
      const size_t numOffsets = curRowOffsets.extent (0);
      const auto valToCheck =
        getEntryOnHost (curRowOffsets, numOffsets - 1);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (numOffsets != 0 &&
         myGraph_->k_lclInds1D_.extent (0) != valToCheck,
         std::logic_error, "(StaticProfile branch) numOffsets = " <<
         numOffsets << " != 0 and myGraph_->k_lclInds1D_.extent(0) = "
         << myGraph_->k_lclInds1D_.extent (0) << " != curRowOffsets("
         << numOffsets << ") = " << valToCheck << ".");
    }

    if (myGraph_->getNodeNumEntries() !=
        myGraph_->getNodeAllocationSize()) {
      if (verbose) {
        std::ostringstream os;
        const auto numEnt = myGraph_->getNodeNumEntries();
        const auto allocSize = myGraph_->getNodeAllocationSize();
        os << *prefix << "Unpacked 1-D storage: numEnt=" << numEnt
           << ", allocSize=" << allocSize << endl;
        std::cerr << os.str ();
      }
      // The matrix's current 1-D storage is "unpacked."  This means
      // the row offsets may differ from what the final row offsets
      // should be.  This could happen, for example, if the user
      // specified StaticProfile in the constructor and set an upper
      // bound on the number of entries per row, but didn't fill all
      // those entries.
      if (debug && curRowOffsets.extent (0) != 0) {
        const size_t numOffsets =
          static_cast<size_t> (curRowOffsets.extent (0));
        const auto valToCheck =
          getEntryOnHost (curRowOffsets, numOffsets - 1);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (static_cast<size_t> (valToCheck) !=
           static_cast<size_t> (k_values1D_.extent (0)),
           std::logic_error, "(StaticProfile unpacked branch) Before "
           "allocating or packing, curRowOffsets(" << (numOffsets-1)
           << ") = " << valToCheck << " != k_values1D_.extent(0)"
           " = " << k_values1D_.extent (0) << ".");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (static_cast<size_t> (valToCheck) !=
           static_cast<size_t> (myGraph_->k_lclInds1D_.extent (0)),
           std::logic_error, "(StaticProfile unpacked branch) Before "
           "allocating or packing, curRowOffsets(" << (numOffsets-1)
           << ") = " << valToCheck
           << " != myGraph_->k_lclInds1D_.extent(0) = "
           << myGraph_->k_lclInds1D_.extent (0) << ".");
      }
      // Pack the row offsets into k_ptrs, by doing a sum-scan of
      // the array of valid entry counts per row.

      // Total number of entries in the matrix on the calling
      // process.  We will compute this in the loop below.  It's
      // cheap to compute and useful as a sanity check.
      size_t lclTotalNumEntries = 0;
      // This will be a host view of packed row offsets.
      typename row_map_type::non_const_type::HostMirror h_ptrs;
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
           "(StaticProfile unpacked branch) After packing k_ptrs, "
           "k_ptrs.extent(0) = " << k_ptrs.extent (0) << " != "
           "lclNumRows+1 = " << (lclNumRows+1) << ".");
        const auto valToCheck = getEntryOnHost (k_ptrs, lclNumRows);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (valToCheck != lclTotalNumEntries, std::logic_error,
           "(StaticProfile unpacked branch) After filling k_ptrs, "
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
      k_inds = lclinds_1d_type ("Tpetra::CrsGraph::ind", lclTotalNumEntries);
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Allocate packed values: "
           << lclTotalNumEntries << endl;
        std::cerr << os.str ();
      }
      k_vals = values_type ("Tpetra::CrsMatrix::val", lclTotalNumEntries);

      // curRowOffsets (myGraph_->k_rowPtrs_) (???), k_lclInds1D_,
      // and k_values1D_ are currently unpacked.  Pack them, using
      // the packed row offsets array k_ptrs that we created above.
      //
      // FIXME (mfh 06 Aug 2014) If "Optimize Storage" is false, we
      // need to keep around the unpacked row offsets, column
      // indices, and values arrays.

      // Pack the column indices from unpacked k_lclInds1D_ into
      // packed k_inds.  We will replace k_lclInds1D_ below.
      using inds_packer_type = pack_functor<
        typename Graph::local_graph_type::entries_type::non_const_type,
        typename Graph::local_graph_type::row_map_type>;
      inds_packer_type indsPacker (k_inds, myGraph_->k_lclInds1D_,
                                   k_ptrs, curRowOffsets);
      using exec_space = typename decltype (k_inds)::execution_space;
      using range_type = Kokkos::RangePolicy<exec_space, LocalOrdinal>;
      Kokkos::parallel_for
        ("Tpetra::CrsMatrix pack column indices",
         range_type (0, lclNumRows), indsPacker);

      // Pack the values from unpacked k_values1D_ into packed
      // k_vals.  We will replace k_values1D_ below.
      using vals_packer_type = pack_functor<values_type, row_map_type>;
      vals_packer_type valsPacker (k_vals, this->k_values1D_,
                                   k_ptrs, curRowOffsets);
      Kokkos::parallel_for ("Tpetra::CrsMatrix pack values",
                            range_type (0, lclNumRows), valsPacker);

      if (debug) {
        const char myPrefix[] = "(StaticProfile \"Optimize Storage\""
          "=true branch) After packing, ";
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (k_ptrs.extent (0) == 0, std::logic_error, myPrefix
           << "k_ptrs.extent(0) = 0.  This probably means that "
           "k_rowPtrs_ was never allocated.");
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
    }
    else { // We don't have to pack, so just set the pointers.
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Storage already packed: k_rowPtrs_: "
           << myGraph_->k_rowPtrs_.extent(0) << ", k_lclInds1D_: "
           << myGraph_->k_lclInds1D_.extent(0) << ", k_values1D_: "
           << k_values1D_.extent(0) << endl;
        std::cerr << os.str();
      }
      k_ptrs_const = myGraph_->k_rowPtrs_;
      k_inds = myGraph_->k_lclInds1D_;
      k_vals = this->k_values1D_;

      if (debug) {
        const char myPrefix[] =
          "(StaticProfile \"Optimize Storage\"=false branch) ";
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (k_ptrs_const.extent (0) == 0, std::logic_error, myPrefix
           << "k_ptrs_const.extent(0) = 0.  This probably means "
           "that k_rowPtrs_ was never allocated.");
        if (k_ptrs_const.extent (0) != 0) {
          const size_t numOffsets (k_ptrs_const.extent (0));
          const auto valToCheck =
            getEntryOnHost (k_ptrs_const, numOffsets - 1);
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (size_t (valToCheck) != k_vals.extent (0),
             std::logic_error, myPrefix <<
             "k_ptrs_const(" << (numOffsets-1) << ") = " << valToCheck
             << " != k_vals.extent(0) = " << k_vals.extent (0) << ".");
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (size_t (valToCheck) != k_inds.extent (0),
             std::logic_error, myPrefix <<
             "k_ptrs_const(" << (numOffsets-1) << ") = " << valToCheck
             << " != k_inds.extent(0) = " << k_inds.extent (0) << ".");
        }
      }
    }

    if (debug) {
      const char myPrefix[] = "After packing, ";
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (size_t (k_ptrs_const.extent (0)) != size_t (lclNumRows + 1),
         std::logic_error, myPrefix << "k_ptrs_const.extent(0) = "
         << k_ptrs_const.extent (0) << " != lclNumRows+1 = " <<
         (lclNumRows+1) << ".");
      if (k_ptrs_const.extent (0) != 0) {
        const size_t numOffsets (k_ptrs_const.extent (0));
        const size_t k_ptrs_const_numOffsetsMinus1 =
          getEntryOnHost (k_ptrs_const, numOffsets - 1);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (k_ptrs_const_numOffsetsMinus1 != size_t (k_vals.extent (0)),
           std::logic_error, myPrefix << "k_ptrs_const(" <<
           (numOffsets-1) << ") = " << k_ptrs_const_numOffsetsMinus1
           << " != k_vals.extent(0) = " << k_vals.extent (0) << ".");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (k_ptrs_const_numOffsetsMinus1 != size_t (k_inds.extent (0)),
           std::logic_error, myPrefix << "k_ptrs_const(" <<
           (numOffsets-1) << ") = " << k_ptrs_const_numOffsetsMinus1
           << " != k_inds.extent(0) = " << k_inds.extent (0) << ".");
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
    // zero rows on this process.  It's impossible for the graph to
    // have dynamic profile (getProfileType() == DynamicProfile) and
    // be optimized (isStorageOptimized()).
    if (requestOptimizedStorage) {
      // Free the old, unpacked, unoptimized allocations.
      // Change the graph from dynamic to static allocation profile

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
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Assign k_rowPtrs_: old="
           << myGraph_->k_rowPtrs_.extent(0) << ", new="
           << k_ptrs_const.extent(0) << endl;
        std::cerr << os.str();
      }
      myGraph_->k_rowPtrs_ = k_ptrs_const;
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Assign k_lclInds1D_: old="
           << myGraph_->k_lclInds1D_.extent(0) << ", new="
           << k_inds.extent(0) << endl;
        std::cerr << os.str();
      }
      myGraph_->k_lclInds1D_ = k_inds;
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Assign k_values1D_: old="
           << k_values1D_.extent(0) << ", new="
           << k_vals.extent(0) << endl;
        std::cerr << os.str();
      }
      this->k_values1D_ = k_vals;

      myGraph_->storageStatus_ = Details::STORAGE_1D_PACKED;
      this->storageStatus_ = Details::STORAGE_1D_PACKED;
    }
    else {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "User requestetd NOT to optimize storage"
           << endl;
        std::cerr << os.str();
      }
    }

    // Make the local graph, using the arrays of row offsets and
    // column indices that we built above.  The local graph should be
    // null, but we delete it first so that any memory can be freed
    // before we allocate the new one.
    //
    // FIXME (mfh 06,28 Aug 2014) It would make more sense for
    // Tpetra::CrsGraph to have a protected method that accepts k_inds
    // and k_ptrs, and creates the local graph lclGraph_.
    myGraph_->lclGraph_ =
      typename Graph::local_graph_type (k_inds, k_ptrs_const);

    // Make the local matrix, using the local graph and vals array.
    auto lclMat = std::make_shared<local_matrix_type>
      ("Tpetra::CrsMatrix::lclMatrix_", getNodeNumCols (),
       k_vals, myGraph_->lclGraph_);
    lclMatrix_ = std::make_shared<local_multiply_op_type> (lclMat);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  fillLocalMatrix (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    using ::Tpetra::Details::ProfilingRegion;
    using Kokkos::create_mirror_view;
    using Teuchos::ArrayRCP;
    using Teuchos::Array;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    using row_map_type = typename Graph::local_graph_type::row_map_type;
    using non_const_row_map_type = typename row_map_type::non_const_type;
    using values_type = typename local_matrix_type::values_type;
    ProfilingRegion regionFLM("Tpetra::CrsMatrix::fillLocalMatrix");
    const size_t lclNumRows = getNodeNumRows();

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
    size_t nodeNumEntries   = staticGraph_->getNodeNumEntries ();
    size_t nodeNumAllocated = staticGraph_->getNodeAllocationSize ();
    row_map_type k_rowPtrs_ = staticGraph_->lclGraph_.row_map;

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

    using row_entries_type = decltype (staticGraph_->k_numRowEntries_);

    // StaticProfile means that the matrix's values are currently
    // stored in a 1-D format.  However, this format is "unpacked";
    // it doesn't necessarily have the same row offsets as indicated
    // by the ptrs array returned by allocRowPtrs.  This could
    // happen, for example, if the user specified StaticProfile in
    // the constructor and fixed the number of matrix entries in
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
        typename row_entries_type::const_type numRowEnt_d =
          staticGraph_->k_numRowEntries_;
        // This function can handle the counts being a host View.
        lclTotalNumEntries =
          Details::computeOffsetsFromCounts (tmpk_ptrs, numRowEnt_d);
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

      // Pack k_values1D_ into k_vals.  We will replace k_values1D_ below.
      pack_functor<values_type, row_map_type> valsPacker
        (k_vals, k_values1D_, tmpk_ptrs, k_rowPtrs_);

      using exec_space = typename decltype (k_vals)::execution_space;
      using range_type = Kokkos::RangePolicy<exec_space, LocalOrdinal>;
      Kokkos::parallel_for ("Tpetra::CrsMatrix pack values",
                            range_type (0, lclNumRows), valsPacker);
    }
    else { // We don't have to pack, so just set the pointer.
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Storage already packed: "
           << "k_values1D_: " << k_values1D_.extent(0) << endl;
        std::cerr << os.str();
      }
      k_vals = k_values1D_;
    }

    // May we ditch the old allocations for the packed one?
    if (requestOptimizedStorage) {
      // The user requested optimized storage, so we can dump the
      // unpacked 1-D storage, and keep the packed storage.
      k_values1D_ = k_vals;
      this->storageStatus_ = Details::STORAGE_1D_PACKED;
    }

    // Build the local sparse matrix object.  At this point, the local
    // matrix certainly has a column Map.  Remember that the local
    // matrix's number of columns comes from the column Map, not the
    // domain Map.
    auto lclMat = std::make_shared<local_matrix_type>
      ("Tpetra::CrsMatrix::lclMatrix_",
       getColMap ()->getNodeNumElements (),
       k_vals, staticGraph_->getLocalGraph ());
    lclMatrix_ = std::make_shared<local_multiply_op_type> (lclMat);
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
      memcpy (&oldRowVals[startOffset], &newRowVals[0],
              numInserted * sizeof (impl_scalar_type));
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  insertLocalValues (const LocalOrdinal lclRow,
                     const Teuchos::ArrayView<const LocalOrdinal>& indices,
                     const Teuchos::ArrayView<const Scalar>& values)
  {
    using std::endl;
    typedef impl_scalar_type IST;
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

    Teuchos::ArrayView<IST> valsView = this->getViewNonConst(rowInfo);
    auto fun = [&](size_t const k, size_t const /*start*/, size_t const offset) {
                 valsView[offset] += values[k]; };
    std::function<void(size_t const, size_t const, size_t const)> cb(std::ref(fun));
    graph.insertLocalIndicesImpl(lclRow, indices, cb);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  insertLocalValues (const LocalOrdinal localRow,
                     const LocalOrdinal numEnt,
                     const Scalar vals[],
                     const LocalOrdinal cols[])
  {
    Teuchos::ArrayView<const LocalOrdinal> colsT (cols, numEnt);
    Teuchos::ArrayView<const Scalar> valsT (vals, numEnt);
    this->insertLocalValues (localRow, colsT, valsT);
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
    typedef impl_scalar_type IST;
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

    Teuchos::ArrayView<IST> valsView = this->getViewNonConst(rowInfo);
    auto fun = [&](size_t const k, size_t const /*start*/, size_t const offset) {
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
        Teuchos::ArrayView<const Scalar> vals2;
        if (this->isGloballyIndexed ()) {
          Teuchos::ArrayView<const GlobalOrdinal> gblColInds2;
          const GlobalOrdinal gblRow =
            graph.rowMap_->getGlobalElement (rowInfo.localRow);
          if (gblRow == Tpetra::Details::OrdinalTraits<GlobalOrdinal>::invalid ()) {
            os << "Local row index " << rowInfo.localRow << " is invalid!" << std::endl;
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
              os << "\tNew global column indices: "
                 << Teuchos::toString (gblColInds2) << std::endl
                 << "\tNew values: " << Teuchos::toString (vals2) << std::endl;
            }
          }
        }
        else if (this->isLocallyIndexed ()) {
          Teuchos::ArrayView<const LocalOrdinal> lclColInds2;
          this->getLocalRowView (rowInfo.localRow, lclColInds2, vals2);
          os << "\tNew local column indices: " << Teuchos::toString (lclColInds2)
             << std::endl;
          os << "\tNew values: " << Teuchos::toString (vals2) << std::endl;
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
        // the special case that the target graph is locally indexed (and the
        // profile type is StaticProfile). In this case, we cannot simply filter
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
                          const LocalOrdinal numElts) const
  {
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const bool sorted = graph.isSorted ();

    size_t hint = 0; // Guess for the current index k into rowVals
    LO numValid = 0; // number of valid local column indices

    // NOTE (mfh 11 Oct 2015) This method assumes UVM.  More
    // accurately, it assumes that the host execution space can
    // access data in both InputMemorySpace and ValsMemorySpace.

    if (graph.isLocallyIndexed ()) {
      // Get a view of the column indices in the row.  This amortizes
      // the cost of getting the view over all the entries of inds.
      auto colInds = graph.getLocalKokkosRowView (rowInfo);

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
      auto colInds = graph.getGlobalKokkosRowView (rowInfo);

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
                      const Teuchos::ArrayView<const Scalar>& vals) const
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
    const Kokkos::View<const impl_scalar_type*, Kokkos::AnonymousSpace>& inputVals) const
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
                      const LocalOrdinal inputCols[]) const
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
    auto curRowVals = this->getRowViewNonConst (rowInfo);
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
                           const LocalOrdinal numElts) const
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
                       const Teuchos::ArrayView<const Scalar>& inputVals) const
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
                       const GlobalOrdinal inputGblColInds[]) const
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

    auto curRowVals = this->getRowViewNonConst (rowInfo);
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
    const Kokkos::View<const impl_scalar_type*, Kokkos::AnonymousSpace>& inputVals) const
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
                           const bool atomic) const
  {
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;

    const bool sorted = graph.isSorted ();

    size_t hint = 0; // guess at the index's relative offset in the row
    LO numValid = 0; // number of valid input column indices

    // NOTE (mfh 11 Oct 2015) This method assumes UVM.  More
    // accurately, it assumes that the host execution space can
    // access data in both InputMemorySpace and ValsMemorySpace.

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
      auto colInds = graph.getLocalKokkosRowView (rowInfo);
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
      auto colInds = graph.getGlobalKokkosRowView (rowInfo);

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
      auto curRowVals = this->getRowViewNonConst (rowInfo);
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
                        const bool atomic) const
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
    auto curRowVals = this->getRowViewNonConst (rowInfo);
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
                         const bool atomic) const
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
    auto curRowVals = this->getRowViewNonConst (rowInfo);
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
                        const bool atomic) const
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
      auto colInds = graph.getLocalKokkosRowView (rowInfo);

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
      auto colInds = graph.getGlobalKokkosRowView (rowInfo);

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
                         const bool atomic) const
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
      auto colInds = graph.getGlobalKokkosRowView (rowInfo);

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
      auto colInds = graph.getLocalKokkosRowView (rowInfo);

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
                          const bool atomic) const
  {
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;

    const bool sorted = graph.isSorted ();

    size_t hint = 0; // Guess for the current index k into rowVals
    LO numValid = 0; // number of valid local column indices

    // NOTE (mfh 11 Oct 2015) This method assumes UVM.  More
    // accurately, it assumes that the host execution space can
    // access data in both InputMemorySpace and ValsMemorySpace.

    if (graph.isLocallyIndexed ()) {
      // Get a view of the column indices in the row.  This amortizes
      // the cost of getting the view over all the entries of inds.
      auto colInds = graph.getLocalKokkosRowView (rowInfo);

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
      auto colInds = graph.getGlobalKokkosRowView (rowInfo);

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
                      const bool atomic) const
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
    const bool atomic) const
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
                      const bool atomic) const
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
    auto curRowVals = this->getRowViewNonConst (rowInfo);
    const IST* const inputVals = reinterpret_cast<const IST*> (vals);
    return this->sumIntoLocalValuesImpl (curRowVals.data (), graph, rowInfo,
                                         cols, inputVals, numEnt, atomic);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayView<const typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::impl_scalar_type>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getView (RowInfo rowinfo) const
  {
    using Kokkos::MemoryUnmanaged;
    using Kokkos::View;
    using Teuchos::ArrayView;
    using ST = impl_scalar_type;
    using range_type = std::pair<size_t, size_t>;

    if (k_values1D_.extent (0) != 0 && rowinfo.allocSize > 0) {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        rowinfo.offset1D + rowinfo.allocSize > k_values1D_.extent (0),
        std::range_error, "Tpetra::CrsMatrix::getView: Invalid access "
        "to 1-D storage of values." << std::endl << "rowinfo.offset1D (" <<
        rowinfo.offset1D << ") + rowinfo.allocSize (" << rowinfo.allocSize <<
        ") > k_values1D_.extent(0) (" << k_values1D_.extent (0) << ").");
#endif // HAVE_TPETRA_DEBUG
      range_type range (rowinfo.offset1D, rowinfo.offset1D + rowinfo.allocSize);
      typedef View<const ST*, execution_space, MemoryUnmanaged> subview_type;
      // mfh 23 Nov 2015: Don't just create a subview of k_values1D_
      // directly, because that first creates a _managed_ subview,
      // then returns an unmanaged version of that.  That touches the
      // reference count, which costs performance in a measurable way.
      // Instead, we create a temporary unmanaged view, then create
      // the subview from that.
      subview_type sv = Kokkos::subview (subview_type (k_values1D_), range);
      const ST* const sv_raw = (rowinfo.allocSize == 0) ? nullptr : sv.data ();
      return ArrayView<const ST> (sv_raw, rowinfo.allocSize);
    }
    else {
      return ArrayView<impl_scalar_type> ();
    }
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getViewRawConst (const impl_scalar_type*& vals,
                   LocalOrdinal& numEnt,
                   const RowInfo& rowinfo) const
  {
#ifdef HAVE_TPETRA_DEBUG
    constexpr bool debug = true;
#else
    constexpr bool debug = false;
#endif // HAVE_TPETRA_DEBUG

    if (k_values1D_.extent (0) != 0 && rowinfo.allocSize > 0) {
      if (debug) {
        if (rowinfo.offset1D + rowinfo.allocSize > k_values1D_.extent (0)) {
          vals = nullptr;
          numEnt = 0;
          return Teuchos::OrdinalTraits<LocalOrdinal>::invalid ();
        }
      }
      vals = k_values1D_.data () + rowinfo.offset1D;
      numEnt = rowinfo.allocSize;
    }
    else {
      vals = nullptr;
      numEnt = 0;
    }

    return static_cast<LocalOrdinal> (0);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getViewRaw (impl_scalar_type*& vals,
              LocalOrdinal& numEnt,
              const RowInfo& rowinfo) const
  {
    const impl_scalar_type* valsConst;
    const LocalOrdinal err = this->getViewRawConst (valsConst, numEnt, rowinfo);
    vals = const_cast<impl_scalar_type*> (valsConst);
    return err;
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Kokkos::View<const typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::impl_scalar_type*,
               typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::execution_space,
               Kokkos::MemoryUnmanaged>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getRowView (const RowInfo& rowInfo) const
  {
    using Kokkos::MemoryUnmanaged;
    using Kokkos::View;
    typedef impl_scalar_type ST;
    typedef View<const ST*, execution_space, MemoryUnmanaged> subview_type;
    typedef std::pair<size_t, size_t> range_type;

    if (k_values1D_.extent (0) != 0 && rowInfo.allocSize > 0) {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION
        (rowInfo.offset1D + rowInfo.allocSize > this->k_values1D_.extent (0),
         std::range_error, "Tpetra::CrsMatrix::getRowView: Invalid access "
         "to 1-D storage of values.  rowInfo.offset1D ("
         << rowInfo.offset1D << ") + rowInfo.allocSize (" << rowInfo.allocSize
         << ") > this->k_values1D_.extent(0) ("
         << this->k_values1D_.extent (0) << ").");
#endif // HAVE_TPETRA_DEBUG
      range_type range (rowInfo.offset1D, rowInfo.offset1D + rowInfo.allocSize);
      // mfh 23 Nov 2015: Don't just create a subview of k_values1D_
      // directly, because that first creates a _managed_ subview,
      // then returns an unmanaged version of that.  That touches the
      // reference count, which costs performance in a measurable way.
      // Instead, we create a temporary unmanaged view, then create
      // the subview from that.
      return Kokkos::subview (subview_type (this->k_values1D_), range);
    }
    else {
      return subview_type ();
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Kokkos::View<typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::impl_scalar_type*,
               typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::execution_space,
               Kokkos::MemoryUnmanaged>
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getRowViewNonConst (const RowInfo& rowInfo) const
  {
    using Kokkos::MemoryUnmanaged;
    using Kokkos::View;
    typedef impl_scalar_type ST;
    typedef View<ST*, execution_space, MemoryUnmanaged> subview_type;
    typedef std::pair<size_t, size_t> range_type;

    if (k_values1D_.extent (0) != 0 && rowInfo.allocSize > 0) {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION
        (rowInfo.offset1D + rowInfo.allocSize > this->k_values1D_.extent (0),
         std::range_error, "Tpetra::CrsMatrix::getRowViewNonConst: Invalid "
         "access to 1-D storage of values.  rowInfo.offset1D ("
         << rowInfo.offset1D << ") + rowInfo.allocSize (" << rowInfo.allocSize
         << ") > this->k_values1D_.extent(0) ("
         << this->k_values1D_.extent (0) << ").");
#endif // HAVE_TPETRA_DEBUG
      range_type range (rowInfo.offset1D, rowInfo.offset1D + rowInfo.allocSize);
      // mfh 23 Nov 2015: Don't just create a subview of k_values1D_
      // directly, because that first creates a _managed_ subview,
      // then returns an unmanaged version of that.  That touches the
      // reference count, which costs performance in a measurable way.
      // Instead, we create a temporary unmanaged view, then create
      // the subview from that.
      return Kokkos::subview (subview_type (this->k_values1D_), range);
    }
    else {
      return subview_type ();
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayView<typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::impl_scalar_type>
  CrsMatrix<Scalar, LocalOrdinal,GlobalOrdinal, Node>::
  getViewNonConst (const RowInfo& rowinfo) const
  {
    return Teuchos::av_const_cast<impl_scalar_type> (this->getView (rowinfo));
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalRowCopy (LocalOrdinal localRow,
                   const Teuchos::ArrayView<LocalOrdinal>& indices,
                   const Teuchos::ArrayView<Scalar>& values,
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
        const LocalOrdinal* curLclInds;
        const impl_scalar_type* curVals;
        LocalOrdinal numSpots; // includes both current entries and extra space

        // If we got this far, rowinfo should be correct and should
        // refer to a valid local row.  Thus, these error checks are
        // superfluous, but we retain them in a debug build.
#ifdef HAVE_TPETRA_DEBUG
        int err =
          staticGraph_->getLocalViewRawConst (curLclInds, numSpots, rowinfo);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (err != static_cast<LocalOrdinal> (0), std::logic_error,
           "staticGraph_->getLocalViewRawConst returned nonzero error code "
           << err << ".");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (static_cast<size_t> (numSpots) < theNumEntries, std::logic_error,
           "numSpots = " << numSpots << " < theNumEntries = " << theNumEntries
           << ".");
        const LocalOrdinal numSpotsBefore = numSpots;
        err = getViewRawConst (curVals, numSpots, rowinfo);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (err != static_cast<LocalOrdinal> (0), std::logic_error,
           "getViewRaw returned nonzero error code " << err << ".");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (numSpotsBefore != numSpots, std::logic_error,
           "numSpotsBefore = " << numSpotsBefore << " != numSpots = "
           << numSpots << ".");
#else
        (void) staticGraph_->getLocalViewRawConst (curLclInds, numSpots, rowinfo);
        (void) getViewRawConst (curVals, numSpots, rowinfo);
#endif // HAVE_TPETRA_DEBUG

        for (size_t j = 0; j < theNumEntries; ++j) {
          values[j] = curVals[j];
          indices[j] = curLclInds[j];
        }
      }
      else if (staticGraph_->isGloballyIndexed ()) {
        // Don't call getColMap(), because it touches RCP's reference count.
        const map_type& colMap = * (staticGraph_->colMap_);
        const GlobalOrdinal* curGblInds;
        const impl_scalar_type* curVals;
        LocalOrdinal numSpots; // includes both current entries and extra space

        // If we got this far, rowinfo should be correct and should
        // refer to a valid local row.  Thus, these error checks are
        // superfluous, but we retain them in a debug build.
#ifdef HAVE_TPETRA_DEBUG
        int err =
          staticGraph_->getGlobalViewRawConst (curGblInds, numSpots, rowinfo);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (err != static_cast<LocalOrdinal> (0), std::logic_error,
           "staticGraph_->getGlobalViewRawConst returned nonzero error code "
           << err << ".");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (static_cast<size_t> (numSpots) < theNumEntries, std::logic_error,
           "numSpots = " << numSpots << " < theNumEntries = " << theNumEntries
           << ".");
        const LocalOrdinal numSpotsBefore = numSpots;
        err = getViewRawConst (curVals, numSpots, rowinfo);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (err != static_cast<LocalOrdinal> (0), std::logic_error,
           "getViewRawConst returned nonzero error code " << err << ".");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (numSpotsBefore != numSpots, std::logic_error,
           "numSpotsBefore = " << numSpotsBefore << " != numSpots = "
           << numSpots << ".");
#else
        (void) staticGraph_->getGlobalViewRawConst (curGblInds, numSpots, rowinfo);
        (void) getViewRawConst (curVals, numSpots, rowinfo);
#endif //HAVE_TPETRA_DEBUG

        for (size_t j = 0; j < theNumEntries; ++j) {
          values[j] = curVals[j];
          indices[j] = colMap.getLocalElement (curGblInds[j]);
        }
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getGlobalRowCopy (GlobalOrdinal globalRow,
                    const Teuchos::ArrayView<GlobalOrdinal>& indices,
                    const Teuchos::ArrayView<Scalar>& values,
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
        const LocalOrdinal* curLclInds;
        const impl_scalar_type* curVals;
        LocalOrdinal numSpots; // includes both current entries and extra space

        // If we got this far, rowinfo should be correct and should
        // refer to a valid local row.  Thus, these error checks are
        // superfluous, but we retain them in a debug build.
#ifdef HAVE_TPETRA_DEBUG
        int err =
          staticGraph_->getLocalViewRawConst (curLclInds, numSpots, rowinfo);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (err != static_cast<LocalOrdinal> (0), std::logic_error,
           "staticGraph_->getLocalViewRawConst returned nonzero error code "
           << err << ".");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (static_cast<size_t> (numSpots) < theNumEntries, std::logic_error,
           "numSpots = " << numSpots << " < theNumEntries = " << theNumEntries
           << ".");
        const LocalOrdinal numSpotsBefore = numSpots;
        err = getViewRawConst (curVals, numSpots, rowinfo);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (err != static_cast<LocalOrdinal> (0), std::logic_error,
           "getViewRaw returned nonzero error code " << err << ".");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (numSpotsBefore != numSpots, std::logic_error,
           "numSpotsBefore = " << numSpotsBefore << " != numSpots = "
           << numSpots << ".");
#else
        (void) staticGraph_->getLocalViewRawConst (curLclInds, numSpots, rowinfo);
        (void) getViewRawConst (curVals, numSpots, rowinfo);
#endif //HAVE_TPETRA_DEBUG

        for (size_t j = 0; j < theNumEntries; ++j) {
          values[j] = curVals[j];
          indices[j] = colMap.getGlobalElement (curLclInds[j]);
        }
      }
      else if (staticGraph_->isGloballyIndexed ()) {
        const GlobalOrdinal* curGblInds;
        const impl_scalar_type* curVals;
        LocalOrdinal numSpots; // includes both current entries and extra space

        // If we got this far, rowinfo should be correct and should
        // refer to a valid local row.  Thus, these error checks are
        // superfluous, but we retain them in a debug build.
#ifdef HAVE_TPETRA_DEBUG
        int err =
          staticGraph_->getGlobalViewRawConst (curGblInds, numSpots, rowinfo);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (err != static_cast<LocalOrdinal> (0), std::logic_error,
           "staticGraph_->getGlobalViewRawConst returned nonzero error code "
           << err << ".");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (static_cast<size_t> (numSpots) < theNumEntries, std::logic_error,
           "numSpots = " << numSpots << " < theNumEntries = " << theNumEntries
           << ".");
        const LocalOrdinal numSpotsBefore = numSpots;
        err = getViewRawConst (curVals, numSpots, rowinfo);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (err != static_cast<LocalOrdinal> (0), std::logic_error,
           "getViewRawConst returned nonzero error code " << err << ".");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (numSpotsBefore != numSpots, std::logic_error,
           "numSpotsBefore = " << numSpotsBefore << " != numSpots = "
           << numSpots << ".");
#else
        (void) staticGraph_->getGlobalViewRawConst (curGblInds, numSpots, rowinfo);
        (void) getViewRawConst (curVals, numSpots, rowinfo);
#endif //HAVE_TPETRA_DEBUG

        for (size_t j = 0; j < theNumEntries; ++j) {
          values[j] = curVals[j];
          indices[j] = curGblInds[j];
        }
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalRowView (LocalOrdinal localRow,
                   Teuchos::ArrayView<const LocalOrdinal>& indices,
                   Teuchos::ArrayView<const Scalar>& values) const
  {
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    typedef LocalOrdinal LO;
    const char tfecfFuncName[] = "getLocalRowView: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isGloballyIndexed (), std::runtime_error, "The matrix currently stores "
      "its indices as global indices, so you cannot get a view with local "
      "column indices.  If the matrix has a column Map, you may call "
      "getLocalRowCopy() to get local column indices; otherwise, you may get "
      "a view with global column indices by calling getGlobalRowCopy().");
    indices = Teuchos::null;
    values = Teuchos::null;
    const RowInfo rowinfo = staticGraph_->getRowInfo (localRow);
    if (rowinfo.localRow != Teuchos::OrdinalTraits<size_t>::invalid () &&
        rowinfo.numEntries > 0) {
      ArrayView<const LO> indTmp = staticGraph_->getLocalView (rowinfo);
      ArrayView<const Scalar> valTmp =
        av_reinterpret_cast<const Scalar> (this->getView (rowinfo));
      indices = indTmp (0, rowinfo.numEntries);
      values = valTmp (0, rowinfo.numEntries);
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
       static_cast<size_t> (rowinfo.numEntries), std::logic_error,
       "At the end of this method, for local row " << localRow << ", "
       "indices.size() = " << indices.size () << " != rowinfo.numEntries = "
       << rowinfo.numEntries << suffix);
    const size_t expectedNumEntries = getNumEntriesInLocalRow (localRow);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (rowinfo.numEntries != expectedNumEntries, std::logic_error, "At the end "
       "of this method, for local row " << localRow << ", rowinfo.numEntries = "
       << rowinfo.numEntries << " != getNumEntriesInLocalRow(localRow) = " <<
       expectedNumEntries << suffix);
#endif // HAVE_TPETRA_DEBUG
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalRowView (const LocalOrdinal lclRow,
                   LocalOrdinal& numEnt,
                   const impl_scalar_type*& val,
                   const LocalOrdinal*& ind) const
  {
    typedef LocalOrdinal LO;

    // Don't call getCrsGraph(), because that modfies an RCP reference
    // count, which is not thread safe.  Checking whether an RCP is
    // null does NOT modify its reference count, and is therefore
    // thread safe.  Note that isGloballyIndexed() calls
    // getCrsGraph(), so we have to go to the graph directly.
    if (staticGraph_.is_null () || staticGraph_->isGloballyIndexed ()) {
      return Tpetra::Details::OrdinalTraits<LO>::invalid ();
    }
    else {
      const RowInfo rowInfo = staticGraph_->getRowInfo (lclRow);
      if (rowInfo.localRow == Tpetra::Details::OrdinalTraits<size_t>::invalid ()) {
        numEnt = 0; // no valid entries in this row on the calling process
        val = nullptr;
        ind = nullptr;
        // First argument (lclRow) invalid, so make 1 the error code.
        return static_cast<LO> (1);
      }
      else {
        numEnt = static_cast<LO> (rowInfo.numEntries);
        auto lclColInds = staticGraph_->getLocalKokkosRowView (rowInfo);
        ind = lclColInds.data (); // FIXME (mfh 18 Jul 2016) UVM
        const LO err = this->getViewRawConst (val, numEnt, rowInfo);
        return err;
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalRowViewRaw (const LocalOrdinal lclRow,
                      LocalOrdinal& numEnt,
                      const LocalOrdinal*& lclColInds,
                      const Scalar*& vals) const
  {
    const impl_scalar_type* vals_ist = nullptr;
    const LocalOrdinal errCode =
      this->getLocalRowView (lclRow, numEnt, vals_ist, lclColInds);
    vals = reinterpret_cast<const Scalar*> (vals_ist);
    return errCode;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getGlobalRowView (GlobalOrdinal globalRow,
                    Teuchos::ArrayView<const GlobalOrdinal>& indices,
                    Teuchos::ArrayView<const Scalar>& values) const
  {
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "getGlobalRowView: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isLocallyIndexed (), std::runtime_error,
      "The matrix is locally indexed, so we cannot return a view of the row "
      "with global column indices.  Use getGlobalRowCopy() instead.");
    indices = Teuchos::null;
    values  = Teuchos::null;
    const RowInfo rowinfo =
      staticGraph_->getRowInfoFromGlobalRowIndex (globalRow);
    if (rowinfo.localRow != Teuchos::OrdinalTraits<size_t>::invalid () &&
        rowinfo.numEntries > 0) {
      ArrayView<const GO> indTmp = staticGraph_->getGlobalView (rowinfo);
      ArrayView<const Scalar> valTmp =
        av_reinterpret_cast<const Scalar> (this->getView (rowinfo));
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (static_cast<size_t> (indTmp.size ()) < rowinfo.numEntries ||
         static_cast<size_t> (valTmp.size ()) < rowinfo.numEntries,
         std::logic_error, std::endl << "rowinfo.numEntries not accurate.  "
         << std::endl << "indTmp.size() = " << indTmp.size ()
         << ", valTmp.size() = " << valTmp.size ()
         << ", rowinfo.numEntries = " << rowinfo.numEntries << ".");
#endif // HAVE_TPETRA_DEBUG
      indices = indTmp (0, rowinfo.numEntries);
      values = valTmp (0, rowinfo.numEntries);
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
       static_cast<size_t> (rowinfo.numEntries), std::logic_error,
       "At the end of this method, for global row " << globalRow << ", "
       "indices.size() = " << indices.size () << " != rowinfo.numEntries = "
       << rowinfo.numEntries << suffix);
    const size_t expectedNumEntries = getNumEntriesInGlobalRow (globalRow);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (rowinfo.numEntries != expectedNumEntries, std::logic_error, "At the end "
       "of this method, for global row " << globalRow << ", rowinfo.numEntries "
       "= " << rowinfo.numEntries << " != getNumEntriesInGlobalRow(globalRow) ="
       " " << expectedNumEntries << suffix);
#endif // HAVE_TPETRA_DEBUG
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  scale (const Scalar& alpha)
  {
    typedef LocalOrdinal LO;
    const char tfecfFuncName[] = "scale: ";
    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillActive (), std::runtime_error,
      "Fill must be active before you may call this method.  "
      "Please call resumeFill() to make fill active.");

    const size_t nlrs = staticGraph_->getNodeNumRows ();
    const size_t numEntries = staticGraph_->getNodeNumEntries ();
    if (! staticGraph_->indicesAreAllocated () ||
        nlrs == 0 || numEntries == 0) {
      // do nothing
    }
    else {
      auto lclMat = this->getLocalMatrix ();

      const LO lclNumRows = lclMat.numRows ();
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        auto row_i = lclMat.row (lclRow);
        for (LO k = 0; k < row_i.length; ++k) {
          // FIXME (mfh 02 Jan 2015) This assumes CUDA UVM.
          row_i.value (k) *= theAlpha;
        }
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  setAllToScalar (const Scalar& alpha)
  {
    const char tfecfFuncName[] = "setAllToScalar: ";
    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillActive (), std::runtime_error,
      "Fill must be active before you may call this method.  "
      "Please call resumeFill() to make fill active.");

    // replace all values in the matrix
    // it is easiest to replace all allocated values, instead of replacing only the ones with valid entries
    // however, if there are no valid entries, we can short-circuit
    // furthermore, if the values aren't allocated, we can short-circuit (no entry have been inserted so far)
    const size_t numEntries = staticGraph_->getNodeNumEntries();
    if (! staticGraph_->indicesAreAllocated () || numEntries == 0) {
      // do nothing
    }
    else {
      // FIXME (mfh 24 Dec 2014) Once CrsMatrix implements DualView
      // semantics, this would be the place to mark memory as
      // modified.
      Kokkos::deep_copy (k_values1D_, theAlpha);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  setAllValues (const typename local_matrix_type::row_map_type& rowPointers,
                const typename local_graph_type::entries_type::non_const_type& columnIndices,
                const typename local_matrix_type::values_type& values)
  {
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
    auto lclGraph = myGraph_->getLocalGraph ();
    const size_t numEnt = lclGraph.entries.extent (0);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (lclGraph.row_map.extent (0) != rowPointers.extent (0) ||
       numEnt != static_cast<size_t> (columnIndices.extent (0)),
       std::logic_error, "myGraph_->setAllIndices() did not correctly create "
       "local graph.  Please report this bug to the Tpetra developers.");

    const size_t numCols = myGraph_->getColMap ()->getNodeNumElements ();

    auto lclMat = std::make_shared<local_matrix_type>
      ("Tpetra::CrsMatrix::lclMatrix_", numCols, values, lclGraph);
    lclMatrix_ = std::make_shared<local_multiply_op_type> (lclMat);

    // FIXME (22 Jun 2016) I would very much like to get rid of
    // k_values1D_ at some point.  I find it confusing to have all
    // these extra references lying around.
    k_values1D_ = lclMat->values;

    // Storage MUST be packed, since the interface doesn't give any
    // way to indicate any extra space at the end of each row.
    this->storageStatus_ = Details::STORAGE_1D_PACKED;

    checkInternalState ();
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
    typedef typename local_matrix_type::row_map_type row_map_type;
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

    const size_t lclNumRows = staticGraph_->getNodeNumRows ();
    if (static_cast<size_t> (offsets.size ()) < lclNumRows) {
      offsets.resize (lclNumRows);
    }

    // The input ArrayRCP must always be a host pointer.  Thus, if
    // device_type::memory_space is Kokkos::HostSpace, it's OK for us
    // to write to that allocation directly as a Kokkos::View.
    typedef typename device_type::memory_space memory_space;
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
      Kokkos::deep_copy (offsetsOut, offsetsTmp);
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
    const LO myNumRows = static_cast<LO> (this->getNodeNumRows ());

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
      diag.template modify<device_type> ();
      const auto D_lcl = diag.template getLocalView<device_type> ();
      // 1-D subview of the first (and only) column of D_lcl.
      const auto D_lcl_1d =
        Kokkos::subview (D_lcl, Kokkos::make_pair (LO (0), myNumRows), 0);

      const auto lclRowMap = rowMap.getLocalMap ();
      const auto lclColMap = colMap.getLocalMap ();
      using ::Tpetra::Details::getDiagCopyWithoutOffsets;
      (void) getDiagCopyWithoutOffsets (D_lcl_1d, lclRowMap,
                                        lclColMap,
                                        lclMatrix_->getLocalMatrix ());
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

    diag.template modify<device_type> ();
    auto D_lcl = diag.template getLocalView<device_type> ();
    const LO myNumRows = static_cast<LO> (this->getNodeNumRows ());
    // Get 1-D subview of the first (and only) column of D_lcl.
    auto D_lcl_1d =
      Kokkos::subview (D_lcl, Kokkos::make_pair (LO (0), myNumRows), 0);

    KokkosSparse::getDiagCopy (D_lcl_1d, offsets,
                               lclMatrix_->getLocalMatrix ());
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
    diag.clear_sync_state ();

    // For now, we fill the Vector on the host and sync to device.
    // Later, we may write a parallel kernel that works entirely on
    // device.
    diag.modify_host ();
    auto lclVecHost = diag.getLocalViewHost ();
    // 1-D subview of the first (and only) column of lclVecHost.
    auto lclVecHost1d = Kokkos::subview (lclVecHost, Kokkos::ALL (), 0);

    using host_offsets_view_type =
      Kokkos::View<const size_t*, Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
    host_offsets_view_type h_offsets (offsets.getRawPtr (), offsets.size ());
    // Find the diagonal entries and put them in lclVecHost1d.
    using range_type = Kokkos::RangePolicy<host_execution_space, LO>;
    const LO myNumRows = static_cast<LO> (this->getNodeNumRows ());
    const size_t INV = Tpetra::Details::OrdinalTraits<size_t>::invalid ();

    local_matrix_type lclMat = lclMatrix_->getLocalMatrix ();
    Kokkos::parallel_for
      ("Tpetra::CrsMatrix::getLocalDiagCopy",
       range_type (0, myNumRows),
       [&] (const LO lclRow) {
        lclVecHost1d(lclRow) = STS::zero (); // default value if no diag entry
        if (h_offsets[lclRow] != INV) {
          auto curRow = lclMat.rowConst (lclRow);
          lclVecHost1d(lclRow) = static_cast<IST> (curRow.value(h_offsets[lclRow]));
        }
      });
    diag.sync_device ();
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
      using dev_memory_space = typename device_type::memory_space;
      if (xp->template need_sync<dev_memory_space> ()) {
        using Teuchos::rcp_const_cast;
        rcp_const_cast<vec_type> (xp)->template sync<dev_memory_space> ();
      }
      auto x_lcl = xp->template getLocalView<dev_memory_space> ();
      auto x_lcl_1d = Kokkos::subview (x_lcl, Kokkos::ALL (), 0);
      using ::Tpetra::Details::leftScaleLocalCrsMatrix;
      leftScaleLocalCrsMatrix (lclMatrix_->getLocalMatrix (),
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
      using dev_memory_space = typename device_type::memory_space;
      if (xp->template need_sync<dev_memory_space> ()) {
        using Teuchos::rcp_const_cast;
        rcp_const_cast<vec_type> (xp)->template sync<dev_memory_space> ();
      }
      auto x_lcl = xp->template getLocalView<dev_memory_space> ();
      auto x_lcl_1d = Kokkos::subview (x_lcl, Kokkos::ALL (), 0);
      using ::Tpetra::Details::rightScaleLocalCrsMatrix;
      rightScaleLocalCrsMatrix (lclMatrix_->getLocalMatrix (),
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
    typedef typename Teuchos::ArrayRCP<const impl_scalar_type>::size_type size_type;

    // FIXME (mfh 05 Aug 2014) Write a thread-parallel kernel for the
    // local part of this computation.  It could make sense to put
    // this operation in the Kokkos::CrsMatrix.

    // check the cache first
    mag_type frobNorm = frobNorm_;
    if (frobNorm == -STM::one ()) {
      mag_type mySum = STM::zero ();
      if (getNodeNumEntries() > 0) {
        if (isStorageOptimized ()) {
          // "Optimized" storage is packed storage.  That means we can
          // iterate in one pass through the 1-D values array.
          const size_type numEntries =
            static_cast<size_type> (getNodeNumEntries ());
          for (size_type k = 0; k < numEntries; ++k) {
            // FIXME (mfh 05 Aug 2014) This assumes UVM.
            const impl_scalar_type val = k_values1D_(k);
            // Note (etp 06 Jan 2015) We need abs() here for composite types
            // (in general, if mag_type is on the left-hand-side, we need
            // abs() on the right-hand-side)
            const mag_type val_abs = STS::abs (val);
            mySum += val_abs * val_abs;
          }
        }
        else {
          const LocalOrdinal numRows =
            static_cast<LocalOrdinal> (this->getNodeNumRows ());
          for (LocalOrdinal r = 0; r < numRows; ++r) {
            const RowInfo rowInfo = myGraph_->getRowInfo (r);
            const size_type numEntries =
              static_cast<size_type> (rowInfo.numEntries);
            ArrayView<const impl_scalar_type> A_r =
              this->getView (rowInfo).view (0, numEntries);
            for (size_type k = 0; k < numEntries; ++k) {
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
      frobNorm = STM::sqrt (totalSum);
    }
    if (isFillComplete ()) {
      // Only cache the result if the matrix is fill complete.
      // Otherwise, the values might still change.  resumeFill clears
      // the cache.
      frobNorm_ = frobNorm;
    }
    return frobNorm;
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
        static_cast<LocalOrdinal> (theGraph.getNodeNumRows ());
      for (LocalOrdinal row = 0; row < lclNumRows; ++row) {
        const RowInfo rowInfo = theGraph.getRowInfo (row);
        auto lclColInds = theGraph.getLocalKokkosRowViewNonConst (rowInfo);
        auto vals = this->getRowViewNonConst (rowInfo);
        // FIXME (mfh 09 May 2017) This assumes CUDA UVM, at least for
        // lclColInds, if not also for values.
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
    // Keep this for CrsGraph's constructor, so we can use StaticProfile.
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
    //    to nonlocal rows.  We may use StaticProfile, since we have
    //    exact counts of the number of entries in each nonlocal row.

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Create nonlocal matrix" << endl;
      std::cerr << os.str();
    }
    RCP<crs_matrix_type> nonlocalMatrix =
      rcp (new crs_matrix_type (nonlocalRowMap, numEntPerNonlocalRow (),
                                StaticProfile));
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
    clearGlobalConstants ();
    fillComplete_ = false;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  computeGlobalConstants ()
  {
    // This method doesn't do anything.  The analogous method in
    // CrsGraph does actually compute something.
    //
    // Oddly enough, clearGlobalConstants() clears frobNorm_ (by
    // setting it to -1), but computeGlobalConstants() does _not_
    // compute the Frobenius norm; this is done on demand in
    // getFrobeniusNorm(), and the result is cached there.
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
  clearGlobalConstants () {
    // We use -1 to indicate that the Frobenius norm needs to be
    // recomputed, since the values might change between now and the
    // next fillComplete call.
    //
    // Oddly enough, clearGlobalConstants() clears frobNorm_, but
    // computeGlobalConstants() does _not_ compute the Frobenius norm;
    // this is done on demand in getFrobeniusNorm(), and the result is
    // cached there.
    frobNorm_ = -STM::one ();
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
      // the error state to makeImportExport or
      // computeGlobalConstants, which may do all-reduces and thus may
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

    {
      Details::ProfilingRegion region_ccgc(
        "Tpetra::CrsMatrix::fillCompete", "callComputeGlobalConstamnts"
      );
      const bool callComputeGlobalConstants = params.get () == nullptr ||
	params->get ("compute global constants", true);
      if (callComputeGlobalConstants) {
	this->computeGlobalConstants ();
      }
    }

    // FIXME (mfh 28 Aug 2014) "Preserve Local Graph" bool parameter no longer used.

    this->fillComplete_ = true; // Now we're fill complete!
    {
      Details::ProfilingRegion region_cis(
        "Tpetra::CrsMatrix::fillCompete", "checkInternalState"
      );
      this->checkInternalState ();
    }
  }

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

    const bool callComputeGlobalConstants = params.get () == nullptr ||
      params->get ("compute global constants", true);
    if (callComputeGlobalConstants) {
        this->computeGlobalConstants ();
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
  size_t
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  mergeRowIndicesAndValues (crs_graph_type& graph,
                            const RowInfo& rowInfo)
  {
#ifdef HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "mergeRowIndicesAndValues: ";
#endif // HAVE_TPETRA_DEBUG

    auto rowValues = this->getRowViewNonConst (rowInfo);
    typedef typename std::decay<decltype (rowValues[0]) >::type value_type;
    value_type* rowValueIter = rowValues.data ();
    auto inds_view = graph.getLocalKokkosRowViewNonConst (rowInfo);

    // beg,end define a half-exclusive interval over which to iterate.
    LocalOrdinal* beg = inds_view.data ();
    LocalOrdinal* end = inds_view.data () + rowInfo.numEntries;

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (rowInfo.allocSize != static_cast<size_t> (inds_view.extent (0)) ||
       rowInfo.allocSize != static_cast<size_t> (rowValues.extent (0)),
       std::runtime_error, "rowInfo.allocSize = " << rowInfo.allocSize
       << " != inds_view.extent(0) = " << inds_view.extent (0)
       << " || rowInfo.allocSize = " << rowInfo.allocSize
       << " != rowValues.extent(0) = " << rowValues.extent (0) << ".");
#endif // HAVE_TPETRA_DEBUG

    LocalOrdinal* newend = beg;
    if (beg != end) {
      LocalOrdinal* cur = beg + 1;
      value_type* vcur = rowValueIter + 1;
      value_type* vend = rowValueIter;
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
    const size_t mergedEntries = newend - beg;
    graph.k_numRowEntries_(rowInfo.localRow) = mergedEntries;
    const size_t numDups = rowInfo.numEntries - mergedEntries;
    return numDups;
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
    //typedef Kokkos::RangePolicy<Kokkos::Serial, LO> range_type;
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
      const LO lclNumRows = static_cast<LO> (this->getNodeNumRows ());
      size_t totalNumDups = 0;
      // FIXME (mfh 10 May 2017) This may assume CUDA UVM.
      Kokkos::parallel_reduce (range_type (0, lclNumRows),
        [this, &graph, sorted, merged] (const LO& lclRow, size_t& numDups) {
          const RowInfo rowInfo = graph.getRowInfo (lclRow);
          if (! sorted) {
            auto lclColInds = graph.getLocalKokkosRowViewNonConst (rowInfo);
            auto vals = this->getRowViewNonConst (rowInfo);
            // FIXME (mfh 09 May 2017) This assumes CUDA UVM, at least
            // for lclColInds, if not also for values.
            sort2 (lclColInds.data (),
                   lclColInds.data () + rowInfo.numEntries,
                   vals.data ());
          }
          if (! merged) {
            numDups += this->mergeRowIndicesAndValues (graph, rowInfo);
          }
        }, totalNumDups);
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
        Y_in.doExport (*Y_rowMap, *exporter, ADD);
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
    const bool Y_is_replicated = ! Y_in.isDistributed ();
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
      Y_in.doExport (*importMV_, *importer, ADD);
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

    auto X_lcl = X.getLocalViewDevice ();
    auto Y_lcl = Y.getLocalViewDevice ();
    // TODO (24 Jul 2019) uncomment later; this line of code wasn't
    // here before, so we need to test it separately before pushing.
    //
    // Y.modify_device ();

    const bool debug = ::Tpetra::Details::Behavior::debug ();
    if (debug) {
      const char tfecfFuncName[] = "localApply: ";
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (lclMatrix_.get () == nullptr, std::logic_error,
         "lclMatrix_ not created yet.");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (X.getNumVectors () != Y.getNumVectors (), std::runtime_error,
         "X.getNumVectors() = " << X.getNumVectors () << " != "
         "Y.getNumVectors() = " << Y.getNumVectors () << ".");
      const bool transpose = (mode != Teuchos::NO_TRANS);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! transpose && X.getLocalLength () !=
         getColMap ()->getNodeNumElements (), std::runtime_error,
         "NO_TRANS case: X has the wrong number of local rows.  "
         "X.getLocalLength() = " << X.getLocalLength () << " != "
         "getColMap()->getNodeNumElements() = " <<
         getColMap ()->getNodeNumElements () << ".");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! transpose && Y.getLocalLength () !=
         getRowMap ()->getNodeNumElements (), std::runtime_error,
         "NO_TRANS case: Y has the wrong number of local rows.  "
         "Y.getLocalLength() = " << Y.getLocalLength () << " != "
         "getRowMap()->getNodeNumElements() = " <<
         getRowMap ()->getNodeNumElements () << ".");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (transpose && X.getLocalLength () !=
         getRowMap ()->getNodeNumElements (), std::runtime_error,
         "TRANS or CONJ_TRANS case: X has the wrong number of local "
         "rows.  X.getLocalLength() = " << X.getLocalLength ()
         << " != getRowMap()->getNodeNumElements() = "
         << getRowMap ()->getNodeNumElements () << ".");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (transpose && Y.getLocalLength () !=
         getColMap ()->getNodeNumElements (), std::runtime_error,
         "TRANS or CONJ_TRANS case: X has the wrong number of local "
         "rows.  Y.getLocalLength() = " << Y.getLocalLength ()
         << " != getColMap()->getNodeNumElements() = "
         << getColMap ()->getNodeNumElements () << ".");
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
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (X_lcl.data () == Y_lcl.data () && X_lcl.data () != nullptr,
         std::runtime_error, "X and Y may not alias one another.");
    }

    LocalOrdinal nrows = getNodeNumRows();
    LocalOrdinal maxRowImbalance = 0;
    if(nrows != 0)
      maxRowImbalance = getNodeMaxNumRowEntries() - (getNodeNumEntries() / nrows);
    if(size_t(maxRowImbalance) >= Tpetra::Details::Behavior::rowImbalanceThreshold())
      lclMatrix_->applyImbalancedRows (X_lcl, Y_lcl, mode, alpha, beta);
    else
      lclMatrix_->apply (X_lcl, Y_lcl, mode, alpha, beta);
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

      //Thyra was implicitly assuming that Y gets set to zero / or is overwritten
      //when bets==0. This was not the case with transpose in a multithreaded
      //environment where a multiplication with subsequent atomic_adds is used
      //since 0 is effectively not special cased. Doing the explicit set to zero here
      //This catches cases where Y is nan or inf.
      const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero ();
      if (beta == ZERO) {
        Y.putScalar (ZERO);
      }
      this->applyTranspose (X, Y, mode, alpha, beta);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  gaussSeidel (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
               MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
               const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& D,
               const Scalar& dampingFactor,
               const ESweepDirection direction,
               const int numSweeps) const
  {
    reorderedGaussSeidel (B, X, D, Teuchos::null, dampingFactor, direction, numSweeps);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  reorderedGaussSeidel (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
                        MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                        const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& D,
                        const Teuchos::ArrayView<LocalOrdinal>& rowIndices,
                        const Scalar& dampingFactor,
                        const ESweepDirection direction,
                        const int numSweeps) const
  {
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    using Teuchos::rcpFromRef;
    typedef Scalar ST;

    TEUCHOS_TEST_FOR_EXCEPTION(
      isFillComplete() == false, std::runtime_error,
      "Tpetra::CrsMatrix::gaussSeidel: cannot call this method until "
      "fillComplete() has been called.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      numSweeps < 0,
      std::invalid_argument,
      "Tpetra::CrsMatrix::gaussSeidel: The number of sweeps must be , "
      "nonnegative but you provided numSweeps = " << numSweeps << " < 0.");

    // Translate from global to local sweep direction.
    // While doing this, validate the input.
    ESweepDirection localDirection;
    if (direction == Forward) {
      localDirection = Forward;
    }
    else if (direction == Backward) {
      localDirection = Backward;
    }
    else if (direction == Symmetric) {
      // We'll control local sweep direction manually.
      localDirection = Forward;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        "Tpetra::CrsMatrix::gaussSeidel: The 'direction' enum does not have "
        "any of its valid values: Forward, Backward, or Symmetric.");
    }

    if (numSweeps == 0) {
      return; // Nothing to do.
    }

    // We don't need the Export object because this method assumes
    // that the row, domain, and range Maps are the same.  We do need
    // the Import object, if there is one, though.
    RCP<const import_type> importer = this->getGraph()->getImporter();
    RCP<const export_type> exporter = this->getGraph()->getExporter();
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! exporter.is_null (), std::runtime_error,
      "Tpetra's gaussSeidel implementation requires that the row, domain, "
      "and range Maps be the same.  This cannot be the case, because the "
      "matrix has a nontrivial Export object.");

    RCP<const map_type> domainMap = this->getDomainMap ();
    RCP<const map_type> rangeMap = this->getRangeMap ();
    RCP<const map_type> rowMap = this->getGraph ()->getRowMap ();
    RCP<const map_type> colMap = this->getGraph ()->getColMap ();

#ifdef HAVE_TEUCHOS_DEBUG
    {
      // The relation 'isSameAs' is transitive.  It's also a
      // collective, so we don't have to do a "shared" test for
      // exception (i.e., a global reduction on the test value).
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! X.getMap ()->isSameAs (*domainMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the input "
        "multivector X be in the domain Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! B.getMap ()->isSameAs (*rangeMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the input "
        "B be in the range Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! D.getMap ()->isSameAs (*rowMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the input "
        "D be in the row Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! rowMap->isSameAs (*rangeMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the row Map and the "
        "range Map be the same (in the sense of Tpetra::Map::isSameAs).");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! domainMap->isSameAs (*rangeMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the domain Map and "
        "the range Map of the matrix be the same.");
    }
#else
    // Forestall any compiler warnings for unused variables.
    (void) rangeMap;
    (void) rowMap;
#endif // HAVE_TEUCHOS_DEBUG

    // If B is not constant stride, copy it into a constant stride
    // multivector.  We'l handle the right-hand side B first and deal
    // with X right before the sweeps, to improve locality of the
    // first sweep.  (If the problem is small enough, then that will
    // hopefully keep more of the entries of X in cache.  This
    // optimizes for the typical case of a small number of sweeps.)
    RCP<const MV> B_in;
    if (B.isConstantStride()) {
      B_in = rcpFromRef (B);
    }
    else {
      // The range Map and row Map are the same in this case, so we
      // can use the (possibly cached) row Map multivector to store a
      // constant stride copy of B.  We don't have to copy back, since
      // Gauss-Seidel won't modify B.
      RCP<MV> B_in_nonconst = getRowMapMultiVector (B, true);
      deep_copy (*B_in_nonconst, B); // Copy from B into B_in(_nonconst).
      B_in = rcp_const_cast<const MV> (B_in_nonconst);

      TPETRA_EFFICIENCY_WARNING(
        ! B.isConstantStride (),
        std::runtime_error,
        "gaussSeidel: The current implementation of the Gauss-Seidel kernel "
        "requires that X and B both have constant stride.  Since B does not "
        "have constant stride, we had to make a copy.  This is a limitation of "
        "the current implementation and not your fault, but we still report it "
        "as an efficiency warning for your information.");
    }

    // If X is not constant stride, copy it into a constant stride
    // multivector.  Also, make the column Map multivector X_colMap,
    // and its domain Map view X_domainMap.  (X actually must be a
    // domain Map view of a column Map multivector; exploit this, if X
    // has constant stride.)

    RCP<MV> X_domainMap;
    RCP<MV> X_colMap;
    bool copiedInput = false;

    if (importer.is_null ()) { // Domain and column Maps are the same.
      if (X.isConstantStride ()) {
        X_domainMap = rcpFromRef (X);
        X_colMap = X_domainMap;
        copiedInput = false;
      }
      else {
        // Get a temporary column Map multivector, make a domain Map
        // view of it, and copy X into the domain Map view.  We have
        // to copy here because we won't be doing Import operations.
        X_colMap = getColumnMapMultiVector (X, true);
        X_domainMap = X_colMap; // Domain and column Maps are the same.
        deep_copy (*X_domainMap, X); // Copy X into the domain Map view.
        copiedInput = true;
        TPETRA_EFFICIENCY_WARNING(
          ! X.isConstantStride (), std::runtime_error,
          "Tpetra::CrsMatrix::gaussSeidel: The current implementation of the "
          "Gauss-Seidel kernel requires that X and B both have constant "
          "stride.  Since X does not have constant stride, we had to make a "
          "copy.  This is a limitation of the current implementation and not "
          "your fault, but we still report it as an efficiency warning for "
          "your information.");
      }
    }
    else { // We will be doing Import operations in the sweeps.
      if (X.isConstantStride ()) {
        X_domainMap = rcpFromRef (X);
        // This kernel assumes that X is a domain Map view of a column
        // Map multivector.  We will only check if this is valid if
        // the CMake configure Teuchos_ENABLE_DEBUG is ON.
        X_colMap = X_domainMap->offsetViewNonConst (colMap, 0);

        // FIXME (mfh 19 Mar 2013) Do we need to fill the remote
        // entries of X_colMap with zeros?  Do we need to fill all of
        // X_domainMap initially with zeros?  Ifpack
        // (Ifpack_PointRelaxation.cpp, line 906) creates an entirely
        // new MultiVector each time.

        // Do the first Import for the first sweep.  This simplifies
        // the logic in the sweeps.
        X_colMap->doImport (X, *importer, INSERT);
        copiedInput = false;
      }
      else {
        // Get a temporary column Map multivector X_colMap, and make a
        // domain Map view X_domainMap of it.  Instead of copying, we
        // do an Import from X into X_domainMap.  This saves us a
        // copy, since the Import has to copy the data anyway.
        X_colMap = getColumnMapMultiVector (X, true);
        X_domainMap = X_colMap->offsetViewNonConst (domainMap, 0);
        X_colMap->doImport (X, *importer, INSERT);
        copiedInput = true;
        TPETRA_EFFICIENCY_WARNING(
          ! X.isConstantStride (), std::runtime_error,
          "Tpetra::CrsMatrix::gaussSeidel: The current implementation of the "
          "Gauss-Seidel kernel requires that X and B both have constant stride.  "
          "Since X does not have constant stride, we had to make a copy.  "
          "This is a limitation of the current implementation and not your fault, "
          "but we still report it as an efficiency warning for your information.");
      }
    }

    for (int sweep = 0; sweep < numSweeps; ++sweep) {
      if (! importer.is_null () && sweep > 0) {
        // We already did the first Import for the zeroth sweep.
        X_colMap->doImport (*X_domainMap, *importer, INSERT);
      }

      // Do local Gauss-Seidel.
      if (direction != Symmetric) {
        if (rowIndices.is_null ()) {
          this->template localGaussSeidel<ST, ST> (*B_in, *X_colMap, D,
                                                   dampingFactor,
                                                   localDirection);
        }
        else {
          this->template reorderedLocalGaussSeidel<ST, ST> (*B_in, *X_colMap,
                                                            D, rowIndices,
                                                            dampingFactor,
                                                            localDirection);
        }
      }
      else { // direction == Symmetric
        const bool doImportBetweenDirections = false;
        if (rowIndices.is_null ()) {
          this->template localGaussSeidel<ST, ST> (*B_in, *X_colMap, D,
                                                   dampingFactor,
                                                   Forward);
          // mfh 18 Mar 2013: Aztec's implementation of "symmetric
          // Gauss-Seidel" does _not_ do an Import between the forward
          // and backward sweeps.  This makes sense, because Aztec
          // considers "symmetric Gauss-Seidel" a subdomain solver.
          if (doImportBetweenDirections) {
            // Communicate again before the Backward sweep.
            if (! importer.is_null ()) {
              X_colMap->doImport (*X_domainMap, *importer, INSERT);
            }
          }
          this->template localGaussSeidel<ST, ST> (*B_in, *X_colMap, D,
                                                   dampingFactor,
                                                   Backward);
        }
        else {
          this->template reorderedLocalGaussSeidel<ST, ST> (*B_in, *X_colMap,
                                                            D, rowIndices,
                                                            dampingFactor,
                                                            Forward);
          if (doImportBetweenDirections) {
            // Communicate again before the Backward sweep.
            if (! importer.is_null ()) {
              X_colMap->doImport (*X_domainMap, *importer, INSERT);
            }
          }
          this->template reorderedLocalGaussSeidel<ST, ST> (*B_in, *X_colMap,
                                                            D, rowIndices,
                                                            dampingFactor,
                                                            Backward);
        }
      }
    }

    if (copiedInput) {
      deep_copy (X, *X_domainMap); // Copy back from X_domainMap to X.
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  gaussSeidelCopy (MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                   const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
                   const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& D,
                   const Scalar& dampingFactor,
                   const ESweepDirection direction,
                   const int numSweeps,
                   const bool zeroInitialGuess) const
  {
    reorderedGaussSeidelCopy (X, B, D, Teuchos::null, dampingFactor, direction,
                              numSweeps, zeroInitialGuess);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  reorderedGaussSeidelCopy (MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                            const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
                            const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& D,
                            const Teuchos::ArrayView<LocalOrdinal>& rowIndices,
                            const Scalar& dampingFactor,
                            const ESweepDirection direction,
                            const int numSweeps,
                            const bool zeroInitialGuess) const
  {
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcpFromRef;
    using Teuchos::rcp_const_cast;
    typedef Scalar ST;
    const char prefix[] = "Tpetra::CrsMatrix::(reordered)gaussSeidelCopy: ";
    const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! isFillComplete (), std::runtime_error,
      prefix << "The matrix is not fill complete.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      numSweeps < 0, std::invalid_argument,
      prefix << "The number of sweeps must be nonnegative, "
      "but you provided numSweeps = " << numSweeps << " < 0.");

    // Translate from global to local sweep direction.
    // While doing this, validate the input.
    ESweepDirection localDirection;
    if (direction == Forward) {
      localDirection = Forward;
    }
    else if (direction == Backward) {
      localDirection = Backward;
    }
    else if (direction == Symmetric) {
      // We'll control local sweep direction manually.
      localDirection = Forward;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument,
        prefix << "The 'direction' enum does not have any of its valid "
        "values: Forward, Backward, or Symmetric.");
    }

    if (numSweeps == 0) {
      return;
    }

    RCP<const import_type> importer = this->getGraph ()->getImporter ();
    RCP<const export_type> exporter = this->getGraph ()->getExporter ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! exporter.is_null (), std::runtime_error,
      "This method's implementation currently requires that the matrix's row, "
      "domain, and range Maps be the same.  This cannot be the case, because "
      "the matrix has a nontrivial Export object.");

    RCP<const map_type> domainMap = this->getDomainMap ();
    RCP<const map_type> rangeMap = this->getRangeMap ();
    RCP<const map_type> rowMap = this->getGraph ()->getRowMap ();
    RCP<const map_type> colMap = this->getGraph ()->getColMap ();

#ifdef HAVE_TEUCHOS_DEBUG
    {
      // The relation 'isSameAs' is transitive.  It's also a
      // collective, so we don't have to do a "shared" test for
      // exception (i.e., a global reduction on the test value).
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! X.getMap ()->isSameAs (*domainMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the input "
        "multivector X be in the domain Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! B.getMap ()->isSameAs (*rangeMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the input "
        "B be in the range Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! D.getMap ()->isSameAs (*rowMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the input "
        "D be in the row Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! rowMap->isSameAs (*rangeMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the row Map and the "
        "range Map be the same (in the sense of Tpetra::Map::isSameAs).");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! domainMap->isSameAs (*rangeMap), std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidelCopy requires that the domain Map and "
        "the range Map of the matrix be the same.");
    }
#else
    // Forestall any compiler warnings for unused variables.
    (void) rangeMap;
    (void) rowMap;
#endif // HAVE_TEUCHOS_DEBUG

    // Fetch a (possibly cached) temporary column Map multivector
    // X_colMap, and a domain Map view X_domainMap of it.  Both have
    // constant stride by construction.  We know that the domain Map
    // must include the column Map, because our Gauss-Seidel kernel
    // requires that the row Map, domain Map, and range Map are all
    // the same, and that each process owns all of its own diagonal
    // entries of the matrix.

    RCP<MV> X_colMap;
    RCP<MV> X_domainMap;
    bool copyBackOutput = false;
    if (importer.is_null ()) {
      if (X.isConstantStride ()) {
        X_colMap = rcpFromRef (X);
        X_domainMap = rcpFromRef (X);
        // Column Map and domain Map are the same, so there are no
        // remote entries.  Thus, if we are not setting the initial
        // guess to zero, we don't have to worry about setting remote
        // entries to zero, even though we are not doing an Import in
        // this case.
        if (zeroInitialGuess) {
          X_colMap->putScalar (ZERO);
        }
        // No need to copy back to X at end.
      }
      else { // We must copy X into a constant stride multivector.
        // Just use the cached column Map multivector for that.
        // force=true means fill with zeros, so no need to fill
        // remote entries (not in domain Map) with zeros.
        X_colMap = getColumnMapMultiVector (X, true);
        // X_domainMap is always a domain Map view of the column Map
        // multivector.  In this case, the domain and column Maps are
        // the same, so X_domainMap _is_ X_colMap.
        X_domainMap = X_colMap;
        if (! zeroInitialGuess) { // Don't copy if zero initial guess
          try {
            deep_copy (*X_domainMap , X); // Copy X into constant stride MV
          } catch (std::exception& e) {
            std::ostringstream os;
            os << "Tpetra::CrsMatrix::reorderedGaussSeidelCopy: "
              "deep_copy(*X_domainMap, X) threw an exception: "
               << e.what () << ".";
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, e.what ());
          }
        }
        copyBackOutput = true; // Don't forget to copy back at end.
        TPETRA_EFFICIENCY_WARNING(
          ! X.isConstantStride (),
          std::runtime_error,
          "gaussSeidelCopy: The current implementation of the Gauss-Seidel "
          "kernel requires that X and B both have constant stride.  Since X "
          "does not have constant stride, we had to make a copy.  This is a "
          "limitation of the current implementation and not your fault, but we "
          "still report it as an efficiency warning for your information.");
      }
    }
    else { // Column Map and domain Map are _not_ the same.
      X_colMap = getColumnMapMultiVector (X);
      X_domainMap = X_colMap->offsetViewNonConst (domainMap, 0);

#ifdef HAVE_TPETRA_DEBUG
      auto X_colMap_host_view = X_colMap->getLocalViewHost ();
      auto X_domainMap_host_view = X_domainMap->getLocalViewHost ();

      if (X_colMap->getLocalLength () != 0 && X_domainMap->getLocalLength ()) {
        TEUCHOS_TEST_FOR_EXCEPTION
          (X_colMap_host_view.data () != X_domainMap_host_view.data (),
           std::logic_error, "Tpetra::CrsMatrix::gaussSeidelCopy: Pointer to "
           "start of column Map view of X is not equal to pointer to start of "
           "(domain Map view of) X.  This may mean that Tpetra::MultiVector::"
           "offsetViewNonConst is broken.  "
           "Please report this bug to the Tpetra developers.");
      }

      TEUCHOS_TEST_FOR_EXCEPTION(
        X_colMap_host_view.extent (0) < X_domainMap_host_view.extent (0) ||
        X_colMap->getLocalLength () < X_domainMap->getLocalLength (),
        std::logic_error, "Tpetra::CrsMatrix::gaussSeidelCopy: "
        "X_colMap has fewer local rows than X_domainMap.  "
        "X_colMap_host_view.extent(0) = " << X_colMap_host_view.extent (0)
        << ", X_domainMap_host_view.extent(0) = "
        << X_domainMap_host_view.extent (0)
        << ", X_colMap->getLocalLength() = " << X_colMap->getLocalLength ()
        << ", and X_domainMap->getLocalLength() = "
        << X_domainMap->getLocalLength ()
        << ".  This means that Tpetra::MultiVector::offsetViewNonConst "
        "is broken.  Please report this bug to the Tpetra developers.");

      TEUCHOS_TEST_FOR_EXCEPTION(
        X_colMap->getNumVectors () != X_domainMap->getNumVectors (),
        std::logic_error, "Tpetra::CrsMatrix::gaussSeidelCopy: "
        "X_colMap has a different number of columns than X_domainMap.  "
        "X_colMap->getNumVectors() = " << X_colMap->getNumVectors ()
        << " != X_domainMap->getNumVectors() = "
        << X_domainMap->getNumVectors ()
        << ".  This means that Tpetra::MultiVector::offsetViewNonConst "
        "is broken.  Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

      if (zeroInitialGuess) {
        // No need for an Import, since we're filling with zeros.
        X_colMap->putScalar (ZERO);
      } else {
        // We could just copy X into X_domainMap.  However, that
        // wastes a copy, because the Import also does a copy (plus
        // communication).  Since the typical use case for
        // Gauss-Seidel is a small number of sweeps (2 is typical), we
        // don't want to waste that copy.  Thus, we do the Import
        // here, and skip the first Import in the first sweep.
        // Importing directly from X effects the copy into X_domainMap
        // (which is a view of X_colMap).
        X_colMap->doImport (X, *importer, INSERT);
      }
      copyBackOutput = true; // Don't forget to copy back at end.
    } // if column and domain Maps are (not) the same

    // The Gauss-Seidel / SOR kernel expects multivectors of constant
    // stride.  X_colMap is by construction, but B might not be.  If
    // it's not, we have to make a copy.
    RCP<const MV> B_in;
    if (B.isConstantStride ()) {
      B_in = rcpFromRef (B);
    }
    else {
      // Range Map and row Map are the same in this case, so we can
      // use the cached row Map multivector to store a constant stride
      // copy of B.
      RCP<MV> B_in_nonconst = getRowMapMultiVector (B, true);
      try {
        deep_copy (*B_in_nonconst, B);
      } catch (std::exception& e) {
        std::ostringstream os;
        os << "Tpetra::CrsMatrix::reorderedGaussSeidelCopy: "
          "deep_copy(*B_in_nonconst, B) threw an exception: "
           << e.what () << ".";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, e.what ());
      }
      B_in = rcp_const_cast<const MV> (B_in_nonconst);

      TPETRA_EFFICIENCY_WARNING(
        ! B.isConstantStride (),
        std::runtime_error,
        "gaussSeidelCopy: The current implementation requires that B have "
        "constant stride.  Since B does not have constant stride, we had to "
        "copy it into a separate constant-stride multivector.  This is a "
        "limitation of the current implementation and not your fault, but we "
        "still report it as an efficiency warning for your information.");
    }

    for (int sweep = 0; sweep < numSweeps; ++sweep) {
      if (! importer.is_null () && sweep > 0) {
        // We already did the first Import for the zeroth sweep above,
        // if it was necessary.
        X_colMap->doImport (*X_domainMap, *importer, INSERT);
      }

      // Do local Gauss-Seidel.
      if (direction != Symmetric) {
        if (rowIndices.is_null ()) {
          this->template localGaussSeidel<ST, ST> (*B_in, *X_colMap, D,
                                                   dampingFactor,
                                                   localDirection);
        }
        else {
          this->template reorderedLocalGaussSeidel<ST, ST> (*B_in, *X_colMap,
                                                            D, rowIndices,
                                                            dampingFactor,
                                                            localDirection);
        }
      }
      else { // direction == Symmetric
        if (rowIndices.is_null ()) {
          this->template localGaussSeidel<ST, ST> (*B_in, *X_colMap, D,
                                                   dampingFactor,
                                                   Forward);
          // mfh 18 Mar 2013: Aztec's implementation of "symmetric
          // Gauss-Seidel" does _not_ do an Import between the forward
          // and backward sweeps.  This makes symmetric Gauss-Seidel a
          // symmetric preconditioner if the matrix A is symmetric.  We
          // imitate Aztec's behavior here.
          this->template localGaussSeidel<ST, ST> (*B_in, *X_colMap, D,
                                                   dampingFactor,
                                                   Backward);
        }
        else {
          this->template reorderedLocalGaussSeidel<ST, ST> (*B_in, *X_colMap,
                                                            D, rowIndices,
                                                            dampingFactor,
                                                            Forward);
          this->template reorderedLocalGaussSeidel<ST, ST> (*B_in, *X_colMap,
                                                            D, rowIndices,
                                                            dampingFactor,
                                                            Backward);

        }
      }
    }

    if (copyBackOutput) {
      try {
        deep_copy (X , *X_domainMap); // Copy result back into X.
      } catch (std::exception& e) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error, prefix << "deep_copy(X, *X_domainMap) "
          "threw an exception: " << e.what ());
      }
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
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->isStaticGraph (), std::logic_error, "This matrix (the source "
       "of the conversion) claims to be fill complete, but does not have a "
       "static (i.e., constant) graph.  Please report this bug to the Tpetra "
       "developers.");

    RCP<output_matrix_type> newMatrix
      (new output_matrix_type (this->getCrsGraph ()));
    // Copy old values into new values.  impl_scalar_type and T may
    // differ, so we can't use Kokkos::deep_copy.
    using ::Tpetra::Details::copyConvert;
    copyConvert (newMatrix->lclMatrix_->getLocalMatrix ().values,
                 this->lclMatrix_->getLocalMatrix ().values);
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
         staticGraph_->getNodeAllocationSize() > 0 &&
         staticGraph_->getNodeNumRows() > 0 &&
         k_values1D_.extent (0) == 0,
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
              << staticGraph_->getNodeAllocationSize () << endl;
        }
        out << "Number of entries: " << getNodeNumEntries () << endl
            << "Max number of entries per row: " << getNodeMaxNumRowEntries ()
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
        for (size_t r = 0; r < getNodeNumRows (); ++r) {
          const size_t nE = getNumEntriesInLocalRow(r);
          GlobalOrdinal gid = getRowMap()->getGlobalElement(r);
          out << std::setw(width) << myRank
              << std::setw(width) << gid
              << std::setw(width) << nE;
          if (vl == VERB_EXTREME) {
            if (isGloballyIndexed()) {
              ArrayView<const GlobalOrdinal> rowinds;
              ArrayView<const Scalar> rowvals;
              getGlobalRowView (gid, rowinds, rowvals);
              for (size_t j = 0; j < nE; ++j) {
                out << " (" << rowinds[j]
                    << ", " << rowvals[j]
                    << ") ";
              }
            }
            else if (isLocallyIndexed()) {
              ArrayView<const LocalOrdinal> rowinds;
              ArrayView<const Scalar> rowvals;
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
    typedef RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> row_matrix_type;
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
    using execution_space = typename device_type::execution_space;
    using row_ptrs_type =
      typename local_graph_type::row_map_type::non_const_type;
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

    // Making copies here because k_rowPtrs_ has a const type. Otherwise, we
    // would use it directly.

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Allocate row_ptrs_beg: "
         << myGraph_->k_rowPtrs_.extent(0) << endl;
      std::cerr << os.str();
    }
    using Kokkos::view_alloc;
    using Kokkos::WithoutInitializing;
    row_ptrs_type row_ptr_beg(
      view_alloc("row_ptr_beg", WithoutInitializing),
      myGraph_->k_rowPtrs_.extent(0));
    Kokkos::deep_copy(row_ptr_beg, myGraph_->k_rowPtrs_);

    const size_t N = row_ptr_beg.extent(0) == 0 ? size_t(0) :
      size_t(row_ptr_beg.extent(0) - 1);
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Allocate row_ptrs_end: " << N << endl;
      std::cerr << os.str();
    }
    row_ptrs_type row_ptr_end(
      view_alloc("row_ptr_end", WithoutInitializing), N);

    const bool refill_num_row_entries =
      myGraph_->k_numRowEntries_.extent(0) != 0;

    if (refill_num_row_entries) { // unpacked storage
      // We can't assume correct *this capture until C++17, and it's
      // likely more efficient just to capture what we need anyway.
      auto num_row_entries = myGraph_->k_numRowEntries_;
      Kokkos::parallel_for
        ("Fill end row pointers", range_policy(0, N),
         KOKKOS_LAMBDA (const size_t i) {
          row_ptr_end(i) = row_ptr_beg(i) + num_row_entries(i);
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
      padCrsArrays(row_ptr_beg, row_ptr_end, myGraph_->k_gblInds1D_,
                   k_values1D_, padding, myRank, verbose);
      const auto newValuesLen = k_values1D_.extent(0);
      const auto newColIndsLen = myGraph_->k_gblInds1D_.extent(0);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (newValuesLen != newColIndsLen, std::logic_error,
         ": After padding, k_values1D_.extent(0)=" << newValuesLen
         << " != myGraph_->k_gblInds1D_.extent(0)=" << newColIndsLen
         << suffix);
    }
    else {
      padCrsArrays(row_ptr_beg, row_ptr_end, myGraph_->k_lclInds1D_,
                   k_values1D_, padding, myRank, verbose);
      const auto newValuesLen = k_values1D_.extent(0);
      const auto newColIndsLen = myGraph_->k_lclInds1D_.extent(0);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (newValuesLen != newColIndsLen, std::logic_error,
         ": After padding, k_values1D_.extent(0)=" << newValuesLen
         << " != myGraph_->k_lclInds1D_.extent(0)=" << newColIndsLen
         << suffix);
    }

    if (refill_num_row_entries) {
      auto num_row_entries = myGraph_->k_numRowEntries_;
      Kokkos::parallel_for
        ("Fill num entries", range_policy(0, N),
         KOKKOS_LAMBDA (const size_t i) {
          num_row_entries(i) = row_ptr_end(i) - row_ptr_beg(i);
        });
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Assign myGraph_->k_rowPtrs_; "
         << "old size: " << myGraph_->k_rowPtrs_.extent(0)
         << ", new size: " << row_ptr_beg.extent(0) << endl;
      std::cerr << os.str();
      TEUCHOS_ASSERT( myGraph_->k_rowPtrs_.extent(0) ==
                      row_ptr_beg.extent(0) );
    }
    myGraph_->k_rowPtrs_ = row_ptr_beg;
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
    Array<GO> rowInds;
    Array<Scalar> rowVals;
    const LO numSameIDs_as_LID = static_cast<LO> (numSameIDs);
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
        if (rowLength > static_cast<size_t> (rowInds.size())) {
          rowInds.resize (rowLength);
          rowVals.resize (rowLength);
        }
        // Resizing invalidates an Array's views, so we must make new
        // ones, even if rowLength hasn't changed.
        ArrayView<GO> rowIndsView = rowInds.view (0, rowLength);
        ArrayView<Scalar> rowValsView = rowVals.view (0, rowLength);

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
        rowIndsConstView = rowIndsView.view (0, rowLength);
        rowValsConstView = rowValsView.view (0, rowLength);
      }
      else { // source matrix is globally indexed.
        srcMat.getGlobalRowView(sourceGID, rowIndsConstView,
                                rowValsConstView);
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
          rowInds.resize (rowLength);
          rowVals.resize (rowLength);
        }
        // Resizing invalidates an Array's views, so we must make new
        // ones, even if rowLength hasn't changed.
        ArrayView<GO> rowIndsView = rowInds.view (0, rowLength);
        ArrayView<Scalar> rowValsView = rowVals.view (0, rowLength);

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
        rowIndsConstView = rowIndsView.view (0, rowLength);
        rowValsConstView = rowValsView.view (0, rowLength);
      }
      else {
        srcMat.getGlobalRowView(sourceGID, rowIndsConstView,
                                rowValsConstView);
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
    Array<GO> rowInds;
    Array<Scalar> rowVals;
    const LO numSameIDs_as_LID = static_cast<LO> (numSameIDs);
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
        if (rowLength > static_cast<size_t> (rowInds.size())) {
          rowInds.resize (rowLength);
          rowVals.resize (rowLength);
        }
        // Resizing invalidates an Array's views, so we must make new
        // ones, even if rowLength hasn't changed.
        ArrayView<GO> rowIndsView = rowInds.view (0, rowLength);
        ArrayView<Scalar> rowValsView = rowVals.view (0, rowLength);

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
        rowIndsConstView = rowIndsView.view (0, rowLength);
        rowValsConstView = rowValsView.view (0, rowLength);
      }
      else { // source matrix is globally indexed.
        srcMat.getGlobalRowView(sourceGID, rowIndsConstView,
                                rowValsConstView);
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
        if (rowLength > static_cast<size_t> (rowInds.size ())) {
          rowInds.resize (rowLength);
          rowVals.resize (rowLength);
        }
        // Resizing invalidates an Array's views, so we must make new
        // ones, even if rowLength hasn't changed.
        ArrayView<GO> rowIndsView = rowInds.view (0, rowLength);
        ArrayView<Scalar> rowValsView = rowVals.view (0, rowLength);

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
        rowIndsConstView = rowIndsView.view (0, rowLength);
        rowValsConstView = rowValsView.view (0, rowLength);
      }
      else {
        srcMat.getGlobalRowView(sourceGID, rowIndsConstView,
                                rowValsConstView);
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
    const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& permuteFromLIDs)
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
   size_t& constantNumPackets,
   Distributor& distor)
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
                            constantNumPackets, distor);
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

      using row_matrix_type = RowMatrix<Scalar, LO, GO, Node>;
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
                         constantNumPackets, distor);
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
           size_t& constantNumPackets,
           Distributor& dist) const
  {
    // The call to packNew in packAndPrepare catches and handles any exceptions.
    Details::ProfilingRegion region_pack_new("Tpetra::CrsMatrix::packNew", "Import/Export");
    if (this->isStaticGraph ()) {
      using ::Tpetra::Details::packCrsMatrixNew;
      packCrsMatrixNew (*this, exports, numPacketsPerLID, exportLIDs,
                        constantNumPackets, dist);
    }
    else {
      this->packNonStaticNew (exportLIDs, exports, numPacketsPerLID,
                              constantNumPackets, dist);
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  packNonStaticNew (const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& exportLIDs,
                    Kokkos::DualView<char*, buffer_device_type>& exports,
                    const Kokkos::DualView<size_t*, buffer_device_type>& numPacketsPerLID,
                    size_t& constantNumPackets,
                    Distributor& /* distor */) const
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
    using HES =
      typename View<int*, device_type>::HostMirror::execution_space;
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
    size_t offset = 0; // current index into 'exports' array.
    for (size_t i = 0; i < numExportLIDs; ++i) {
      const LO lclRow = exportLIDs_h[i];

      size_t numEnt;
      numEnt = this->getNumEntriesInLocalRow (lclRow);

      // Only pack this row's data if it has a nonzero number of
      // entries.  We can do this because receiving processes get the
      // number of packets, and will know that zero packets means zero
      // entries.
      if (numEnt == 0) {
        numPacketsPerLID_h[i] = 0;
        continue;
      }

      // Temporary buffer for global column indices.
      using Details::ScalarViewTraits;
      View<GO*, HES> gidsIn_k =
        ScalarViewTraits<GO, HES>::allocateArray (GO (0), numEnt, "gids");

      Teuchos::ArrayView<const Scalar> valsIn;
      if (this->isLocallyIndexed ()) {
        // If the matrix is locally indexed on the calling process, we
        // have to use its column Map (which it _must_ have in this
        // case) to convert to global indices.
        Teuchos::ArrayView<const LO> lidsIn;
        this->getLocalRowView (lclRow, lidsIn, valsIn);
        const map_type& colMap = * (this->getColMap ());
        for (size_t k = 0; k < numEnt; ++k) {
          gidsIn_k[k] = colMap.getGlobalElement (lidsIn[k]);
        }
      }
      else if (this->isGloballyIndexed ()) {
        // If the matrix is globally indexed on the calling process,
        // then we can use the column indices directly.  However, we
        // have to get the global row index.  The calling process must
        // have a row Map, since otherwise it shouldn't be participating
        // in packing operations.
        Teuchos::ArrayView<const GO> gblIndView;;
        const map_type& rowMap = * (this->getRowMap ());
        const GO gblRow = rowMap.getGlobalElement (lclRow);
        this->getGlobalRowView (gblRow, gblIndView, valsIn);
        for (size_t k = 0; k < numEnt; ++k) {
          gidsIn_k[k] = gblIndView[k];
        }
      }
      // mfh 11 Sep 2017: Currently, if the matrix is neither globally
      // nor locally indexed, then it has no entries.  Therefore,
      // there is nothing to pack.  No worries!

      typename HES::device_type outputDevice;
      auto valsIn_k =
        create_mirror_view_from_raw_host_array (outputDevice,
                                                reinterpret_cast<const ST*> (valsIn.getRawPtr ()),
                                                valsIn.size (),
                                                true, "valsIn");
      const size_t numBytesPerValue =
        PackTraits<ST>::packValueCount (valsIn[0]);
      const size_t numBytes =
        this->packRow (exports_h.data (), offset, numEnt, gidsIn_k.data (),
                       valsIn_k.data (), numBytesPerValue);
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
   Distributor& distor,
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
                             constantNumPackets, distor, combineMode,
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
                           constantNumPackets, distor, combineMode,
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
    Distributor & distor,
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
                                   distor, combineMode);
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
                                    distor, combineMode);
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
    Distributor& distor,
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
    using row_matrix_type =
      RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
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
    // sum of the two entry counts in each row.  If we choose this as
    // the actual per-row upper bound, we can use static profile.
    if (A_rowMap->isSameAs (*B_rowMap)) {
      const LO localNumRows = static_cast<LO> (A_rowMap->getNodeNumElements ());
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
        C = rcp (new crs_matrix_type (C_rowMap, C_maxNumEntriesPerRow (),
                                      StaticProfile));
      } else {
        C = rcp (new crs_matrix_type (C_rowMap, C_maxNumEntriesPerRow (),
                                      StaticProfile, constructorSublist));
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

    Array<GO> ind;
    Array<Scalar> val;

    if (alpha != ZERO) {
      const LO A_localNumRows = static_cast<LO> (A_rowMap->getNodeNumElements ());
      for (LO localRow = 0; localRow < A_localNumRows; ++localRow) {
        size_t A_numEntries = A.getNumEntriesInLocalRow (localRow);
        const GO globalRow = A_rowMap->getGlobalElement (localRow);
        if (A_numEntries > static_cast<size_t> (ind.size ())) {
          ind.resize (A_numEntries);
          val.resize (A_numEntries);
        }
        ArrayView<GO> indView = ind (0, A_numEntries);
        ArrayView<Scalar> valView = val (0, A_numEntries);
        A.getGlobalRowCopy (globalRow, indView, valView, A_numEntries);

        if (alpha != ONE) {
          for (size_t k = 0; k < A_numEntries; ++k) {
            valView[k] *= alpha;
          }
        }
        C->insertGlobalValues (globalRow, indView, valView);
      }
    }

    if (beta != ZERO) {
      const LO B_localNumRows = static_cast<LO> (B_rowMap->getNodeNumElements ());
      for (LO localRow = 0; localRow < B_localNumRows; ++localRow) {
        size_t B_numEntries = B.getNumEntriesInLocalRow (localRow);
        const GO globalRow = B_rowMap->getGlobalElement (localRow);
        if (B_numEntries > static_cast<size_t> (ind.size ())) {
          ind.resize (B_numEntries);
          val.resize (B_numEntries);
        }
        ArrayView<GO> indView = ind (0, B_numEntries);
        ArrayView<Scalar> valView = val (0, B_numEntries);
        B.getGlobalRowCopy (globalRow, indView, valView, B_numEntries);

        if (beta != ONE) {
          for (size_t k = 0; k < B_numEntries; ++k) {
            valView[k] *= beta;
          }
        }
        C->insertGlobalValues (globalRow, indView, valView);
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
    typedef CrsMatrix<Scalar, LO, GO, NT> this_type;
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
    if (! params.is_null ()) {
      matrixparams = sublist (params, "CrsMatrix");
      reverseMode = params->get ("Reverse Mode", reverseMode);
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
    // 1. Call the moral equivalent of "distor.do" to handle the import.
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
    ArrayView<const LO> RemoteLIDs = reverseMode ?
      rowTransfer.getExportLIDs () : rowTransfer.getRemoteLIDs ();
    ArrayView<const LO> PermuteToLIDs = reverseMode ?
      rowTransfer.getPermuteFromLIDs () : rowTransfer.getPermuteToLIDs ();
    ArrayView<const LO> PermuteFromLIDs = reverseMode ?
      rowTransfer.getPermuteToLIDs () : rowTransfer.getPermuteFromLIDs ();
    Distributor& Distor = rowTransfer.getDistributor ();

    // Owning PIDs
    Teuchos::Array<int> SourcePids;
    Teuchos::Array<int> TargetPids;

    // Temp variables for sub-communicators
    RCP<const map_type> ReducedRowMap, ReducedColMap,
      ReducedDomainMap, ReducedRangeMap;
    RCP<const Comm<int> > ReducedComm;

    // If the user gave us a null destMat, then construct the new
    // destination matrix.  We will replace its column Map later.
    if (destMat.is_null ()) {
      destMat = rcp (new this_type (MyRowMap, 0, StaticProfile, matrixparams));
    }

    /***************************************************/
    /***** 1) First communicator restriction phase ****/
    /***************************************************/
    if (restrictComm) {
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
    /***** 2) From Tpera::DistObject::doTransfer() ****/
    /***************************************************/
    // Get the owning PIDs
    RCP<const import_type> MyImporter = getGraph ()->getImporter ();

    // check whether domain maps of source matrix and base domain map is the same
    bool bSameDomainMap = BaseDomainMap->isSameAs (*getDomainMap ());

    if (! restrictComm && ! MyImporter.is_null () && bSameDomainMap ) {
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
      IntVectorType SourceDomain_pids(getDomainMap (),true);
      IntVectorType SourceCol_pids(getColMap());
      // SourceDomain_pids contains the restricted pids
      SourceDomain_pids.putScalar(MyPID);

      SourceCol_pids.doImport (SourceDomain_pids, *MyImporter, INSERT);
      SourcePids.resize (getColMap ()->getNodeNumElements ());
      SourceCol_pids.get1dCopy (SourcePids ());
    }
    else if (MyImporter.is_null () && bSameDomainMap) {
      // Matrix has no off-process entries
      SourcePids.resize (getColMap ()->getNodeNumElements ());
      SourcePids.assign (getColMap ()->getNodeNumElements (), MyPID);
    }
    else if ( ! MyImporter.is_null () &&
              ! domainTransfer.is_null () ) {
      // general implementation for rectangular matrices with
      // domain map different than SourceMatrix domain map.
      // User has to provide a DomainTransfer object. We need
      // to communications (import/export)

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
      SourcePids.resize (getColMap ()->getNodeNumElements ());
      SourceCol_pids.get1dCopy (SourcePids ());
    }
    else if ( ! MyImporter.is_null () &&
             BaseDomainMap->isSameAs (*BaseRowMap) &&
             getDomainMap ()->isSameAs (*getRowMap ())) {
      // We can use the rowTransfer + SourceMatrix's Import to find out who owns what.

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
      SourcePids.resize (getColMap ()->getNodeNumElements ());
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
    if (constantNumPackets == 0) {
      destMat->reallocArraysForNumPacketsPerLid (ExportLIDs.size (),
                                                 RemoteLIDs.size ());
    }
    else {
      // There are a constant number of packets per element.  We
      // already know (from the number of "remote" (incoming)
      // elements) how many incoming elements we expect, so we can
      // resize the buffer accordingly.
      const size_t rbufLen = RemoteLIDs.size() * constantNumPackets;
      destMat->reallocImportsIfNeeded (rbufLen, false, nullptr);
    }

    // Pack & Prepare w/ owning PIDs
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
                                     constantNumPackets,
                                     Distor);
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
                                   constantNumPackets,
                                   Distor);
      if (verbose) {
        std::ostringstream os;
        os << *verbosePrefix << "Done with packCrsMatrixWithOwningPIDs"
           << std::endl;
        std::cerr << os.str ();
      }
    }

    // Do the exchange of remote data.
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
          Distor.doReversePostsAndWaits (numExportPacketsPerLID, 1,
                                         numImportPacketsPerLID);
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
          Teuchos::ArrayView<char> hostImports =
            getArrayViewFromDualView (destMat->imports_);
          // This is a legacy host pack/unpack path, so use the host
          // version of exports_.
          destMat->exports_.sync_host ();
          Teuchos::ArrayView<const char> hostExports =
            getArrayViewFromDualView (destMat->exports_);
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
          Teuchos::ArrayView<char> hostImports =
            getArrayViewFromDualView (destMat->imports_);
          // This is a legacy host pack/unpack path, so use the host
          // version of exports_.
          destMat->exports_.sync_host ();
          Teuchos::ArrayView<const char> hostExports =
            getArrayViewFromDualView (destMat->exports_);
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
          Distor.doPostsAndWaits (numExportPacketsPerLID, 1,
                                  numImportPacketsPerLID);
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
          Teuchos::ArrayView<char> hostImports =
            getArrayViewFromDualView (destMat->imports_);
          // This is a legacy host pack/unpack path, so use the host
          // version of exports_.
          destMat->exports_.sync_host ();
          Teuchos::ArrayView<const char> hostExports =
            getArrayViewFromDualView (destMat->exports_);
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
          Teuchos::ArrayView<char> hostImports =
            getArrayViewFromDualView (destMat->imports_);
          // This is a legacy host pack/unpack path, so use the host
          // version of exports_.
          destMat->exports_.sync_host ();
          Teuchos::ArrayView<const char> hostExports =
            getArrayViewFromDualView (destMat->exports_);
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

    /*********************************************************************/
    /**** 3) Copy all of the Same/Permute/Remote data into CSR_arrays ****/
    /*********************************************************************/

    // Backwards compatibility measure.  We'll use this again below.
    destMat->numImportPacketsPerLID_.sync_host ();
    Teuchos::ArrayView<const size_t> numImportPacketsPerLID =
      getArrayViewFromDualView (destMat->numImportPacketsPerLID_);
    destMat->imports_.sync_host ();
    Teuchos::ArrayView<const char> hostImports =
      getArrayViewFromDualView (destMat->imports_);

    if (verbose) {
      std::ostringstream os;
      os << *verbosePrefix << "Calling unpackAndCombineWithOwningPIDsCount"
         << std::endl;
      std::cerr << os.str ();
    }
    size_t mynnz =
      unpackAndCombineWithOwningPIDsCount (*this,
                                           RemoteLIDs,
                                           hostImports,
                                           numImportPacketsPerLID,
                                           constantNumPackets,
                                           Distor,
                                           INSERT,
                                           NumSameIDs,
                                           PermuteToLIDs,
                                           PermuteFromLIDs);
    if (verbose) {
      std::ostringstream os;
      os << *verbosePrefix << "unpackAndCombineWithOwningPIDsCount returned "
         << mynnz << std::endl;
      std::cerr << os.str ();
    }
    size_t N = BaseRowMap->getNodeNumElements ();

    // Allocations
    ArrayRCP<size_t> CSR_rowptr(N+1);
    ArrayRCP<GO> CSR_colind_GID;
    ArrayRCP<LO> CSR_colind_LID;
    ArrayRCP<Scalar> CSR_vals;
    CSR_colind_GID.resize (mynnz);
    CSR_vals.resize (mynnz);

    // If LO and GO are the same, we can reuse memory when
    // converting the column indices from global to local indices.
    if (typeid (LO) == typeid (GO)) {
      CSR_colind_LID = Teuchos::arcp_reinterpret_cast<LO> (CSR_colind_GID);
    }
    else {
      CSR_colind_LID.resize (mynnz);
    }

    if (verbose) {
      std::ostringstream os;
      os << *verbosePrefix << "Calling unpackAndCombineIntoCrsArrays"
         << std::endl;
      std::cerr << os.str ();
    }
    // FIXME (mfh 15 May 2014) Why can't we abstract this out as an
    // unpackAndCombine method on a "CrsArrays" object?  This passing
    // in a huge list of arrays is icky.  Can't we have a bit of an
    // abstraction?  Implementing a concrete DistObject subclass only
    // takes five methods.
    unpackAndCombineIntoCrsArrays (*this,
                                   RemoteLIDs,
                                   hostImports,
                                   numImportPacketsPerLID,
                                   constantNumPackets,
                                   Distor,
                                   INSERT,
                                   NumSameIDs,
                                   PermuteToLIDs,
                                   PermuteFromLIDs,
                                   N,
                                   mynnz,
                                   MyPID,
                                   CSR_rowptr (),
                                   CSR_colind_GID (),
                                   Teuchos::av_reinterpret_cast<impl_scalar_type> (CSR_vals ()),
                                   SourcePids (),
                                   TargetPids);

    /**************************************************************/
    /**** 4) Call Optimized MakeColMap w/ no Directory Lookups ****/
    /**************************************************************/
    // Call an optimized version of makeColMap that avoids the
    // Directory lookups (since the Import object knows who owns all
    // the GIDs).
    Teuchos::Array<int> RemotePids;
    if (verbose) {
      std::ostringstream os;
      os << *verbosePrefix << "Calling lowCommunicationMakeColMapAndReindex"
         << std::endl;
      std::cerr << os.str ();
    }
    Import_Util::lowCommunicationMakeColMapAndReindex (CSR_rowptr (),
                                                       CSR_colind_LID (),
                                                       CSR_colind_GID (),
                                                       BaseDomainMap,
                                                       TargetPids,
                                                       RemotePids,
                                                       MyColMap);

    if (verbose) {
      std::ostringstream os;
      os << *verbosePrefix << "restrictComm="
         << (restrictComm ? "true" : "false") << std::endl;
      std::cerr << os.str ();
    }

    /*******************************************************/
    /**** 4) Second communicator restriction phase      ****/
    /*******************************************************/
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
      Import_Util::sortCrsEntries (CSR_rowptr (),
                                   CSR_colind_LID (),
                                   CSR_vals ());
    }
    else if ((! reverseMode && xferAsExport != nullptr) ||
             (reverseMode && xferAsImport != nullptr)) {
      if (verbose) {
        std::ostringstream os;
        os << *verbosePrefix << "Calling sortAndMergeCrsEntries"
           << endl;
        std::cerr << os.str();
      }
      Import_Util::sortAndMergeCrsEntries (CSR_rowptr (),
                                           CSR_colind_LID (),
                                           CSR_vals ());
      if (CSR_rowptr[N] != mynnz) {
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
    destMat->setAllValues (CSR_rowptr, CSR_colind_LID, CSR_vals);

    /***************************************************/
    /**** 7) Build Importer & Call ESFC             ****/
    /***************************************************/
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
          os << *verbosePrefix << "Calling getAllValues" << endl;
          std::cerr << os.str ();
        }

        Teuchos::ArrayRCP<LocalOrdinal> type3LIDs;
        Teuchos::ArrayRCP<int>          type3PIDs;
        Teuchos::ArrayRCP<const size_t> rowptr;
        Teuchos::ArrayRCP<const LO> colind;
        Teuchos::ArrayRCP<const Scalar> vals;
        {
#ifdef HAVE_TPETRA_MMM_TIMINGS
            TimeMonitor tm_getAllValues (*TimeMonitor::getNewTimer(prefix + std::string("isMMgetAllValues")));
#endif
            getAllValues(rowptr,colind,vals);
        }

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

        const int numCols = getGraph()->getColMap()->getNodeNumElements(); // may be dup
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

    if (verbose) {
      std::ostringstream os;
      os << *verbosePrefix << "Done" << endl;
      std::cerr << os.str ();
    }
  }


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
  template class CrsMatrix< SCALAR , LO , GO , NODE >; \
  template Teuchos::RCP< CrsMatrix< SCALAR , LO , GO , NODE > >   \
                CrsMatrix< SCALAR , LO , GO , NODE >::convert< SCALAR > () const;

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

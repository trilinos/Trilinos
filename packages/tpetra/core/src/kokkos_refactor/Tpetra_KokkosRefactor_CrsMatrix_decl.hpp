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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_KOKKOSREFACTOR_CRSMATRIX_DECL_HPP
#define TPETRA_KOKKOSREFACTOR_CRSMATRIX_DECL_HPP

// This file gets included by Tpetra_CrsMatrix_decl.hpp, so it
// inherits all of that file's includes.  At some point, it will be
// subsumed into Tpetra_CrsMatrix_decl.hpp.

#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#include <Kokkos_CrsMatrix.hpp>

namespace Tpetra {

  /// \brief Partial specialization of CrsMatrix for the new Kokkos
  ///   wrapper Nodes.
  ///
  /// This implements the "Kokkos refactor" version of CrsMatrix.
  /// For full documentation, see the "classic" version of CrsMatrix.
  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  class CrsMatrix<Scalar,
                  LocalOrdinal,
                  GlobalOrdinal,
                  Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>,
                  false> :
    public RowMatrix<Scalar,
                     LocalOrdinal,
                     GlobalOrdinal,
                     Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >,
    public DistObject<char,
                      LocalOrdinal,
                      GlobalOrdinal,
                      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >
  {
  public:
    //! @name Typedefs
    //@{

    /// \brief This class' first template parameter; the type of each
    ///   entry in the matrix.
    typedef Scalar scalar_type;
    /// \brief The type used internally in place of \c Scalar.
    ///
    /// Some \c Scalar types might not work with Kokkos on all
    /// execution spaces, due to missing CUDA device macros or
    /// volatile overloads.  The C++ standard type std::complex<T> has
    /// this problem.  To fix this, we replace std::complex<T> values
    /// internally with the (usually) bitwise identical type
    /// Kokkos::complex<T>.  The latter is the \c impl_scalar_type
    /// corresponding to \c Scalar = std::complex.
    typedef typename Kokkos::Details::ArithTraits<Scalar>::val_type impl_scalar_type;
    //! This class' second template parameter; the type of local indices.
    typedef LocalOrdinal local_ordinal_type;
    //! This class' third template parameter; the type of global indices.
    typedef GlobalOrdinal global_ordinal_type;
    //! This class' fourth template parameter; the Kokkos device type.
    typedef DeviceType execution_space;

    /// \brief Type of a norm result.
    ///
    /// This is usually the same as the type of the magnitude
    /// (absolute value) of <tt>Scalar</tt>, but may differ for
    /// certain <tt>Scalar</tt> types.
    typedef typename Kokkos::Details::ArithTraits<impl_scalar_type>::mag_type mag_type;

    /// \brief The Kokkos Node type used by this class.
    ///
    /// This type depends on the DeviceType template parameter.  In
    /// this, the Kokkos refactor version of Tpetra, it exists only
    /// for backwards compatibility.
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> node_type;

    //! The Map specialization suitable for this CrsMatrix specialization.
    typedef Map<LocalOrdinal, GlobalOrdinal, node_type> map_type;

    //! The Import specialization suitable for this CrsMatrix specialization.
    typedef Import<LocalOrdinal, GlobalOrdinal, node_type> import_type;

    //! The Export specialization suitable for this CrsMatrix specialization.
    typedef Export<LocalOrdinal, GlobalOrdinal, node_type> export_type;

    //! The CrsGraph specialization suitable for this CrsMatrix specialization.
    typedef CrsGraph<LocalOrdinal, GlobalOrdinal, node_type> crs_graph_type;

    /// \brief The specialization of Kokkos::CrsMatrix that represents
    ///   the part of the sparse matrix on each MPI process.
    typedef Kokkos::CrsMatrix<impl_scalar_type, LocalOrdinal, execution_space, void, size_t> local_matrix_type;

    //! DEPRECATED; use <tt>local_matrix_type::row_map_type</tt> instead.
    typedef typename local_matrix_type::row_map_type t_RowPtrs TPETRA_DEPRECATED;
    //! DEPRECATED; use <tt>local_matrix_type::row_map_type::non_const_type</tt> instead.
    typedef typename local_matrix_type::row_map_type::non_const_type t_RowPtrsNC TPETRA_DEPRECATED;

    typedef typename crs_graph_type::t_LocalOrdinal_1D t_LocalOrdinal_1D;
    typedef typename local_matrix_type::values_type t_ValuesType;

    //! DEPRECATED; use local_matrix_type instead.
    typedef local_matrix_type k_local_matrix_type TPETRA_DEPRECATED;

    //@}
    //! @name Constructors and destructor
    //@{

    /// \brief Constructor specifying fixed number of entries for each row.
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param maxNumEntriesPerRow [in] Maximum number of matrix
    ///   entries per row.  If pftype==DynamicProfile, this is only a
    ///   hint, and you can set this to zero without affecting
    ///   correctness.  If pftype==StaticProfile, this sets the amount
    ///   of storage allocated, and you cannot exceed this number of
    ///   entries in any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               size_t maxNumEntriesPerRow,
               ProfileType pftype = DynamicProfile,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying (possibly different) number of entries in each row.
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param NumEntriesPerRowToAlloc [in] Maximum number of matrix
    ///   entries to allocate for each row.  If
    ///   pftype==DynamicProfile, this is only a hint.  If
    ///   pftype==StaticProfile, this sets the amount of storage
    ///   allocated, and you cannot exceed the allocated number of
    ///   entries for any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::ArrayRCP<const size_t>& NumEntriesPerRowToAlloc,
               ProfileType pftype = DynamicProfile,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying column Map and fixed number of entries for each row.
    ///
    /// The column Map will be used to filter any matrix entries
    /// inserted using insertLocalValues() or insertGlobalValues().
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param colMap [in] Distribution of columns of the matrix.
    ///
    /// \param maxNumEntriesPerRow [in] Maximum number of matrix
    ///   entries per row.  If pftype==DynamicProfile, this is only a
    ///   hint, and you can set this to zero without affecting
    ///   correctness.  If pftype==StaticProfile, this sets the amount
    ///   of storage allocated, and you cannot exceed this number of
    ///   entries in any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap,
               size_t maxNumEntriesPerRow,
               ProfileType pftype = DynamicProfile,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying column Map and number of entries in each row.
    ///
    /// The column Map will be used to filter any matrix indices
    /// inserted using insertLocalValues() or insertGlobalValues().
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param colMap [in] Distribution of columns of the matrix.
    ///
    /// \param NumEntriesPerRowToAlloc [in] Maximum number of matrix
    ///   entries to allocate for each row.  If
    ///   pftype==DynamicProfile, this is only a hint.  If
    ///   pftype==StaticProfile, this sets the amount of storage
    ///   allocated, and you cannot exceed the allocated number of
    ///   entries for any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap,
               const Teuchos::ArrayRCP<const size_t>& NumEntriesPerRowToAlloc,
               ProfileType pftype = DynamicProfile,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying a previously constructed graph.
    ///
    /// Calling this constructor fixes the graph structure of the
    /// sparse matrix.  We say in this case that the matrix has a
    /// "static graph."  If you create a CrsMatrix with this
    /// constructor, you are not allowed to insert new entries into
    /// the matrix, but you are allowed to change values in the
    /// matrix.
    ///
    /// The given graph must be fill complete.  Note that calling
    /// resumeFill() on the graph makes it not fill complete, even if
    /// you had previously called fillComplete() on the graph.  In
    /// that case, you must call fillComplete() on the graph again
    /// before invoking this CrsMatrix constructor.
    ///
    /// This constructor is marked \c explicit so that you can't
    /// create a CrsMatrix by accident when passing a CrsGraph into a
    /// function that takes a CrsMatrix.
    ///
    /// \param graph [in] The graph structure of the sparse matrix.
    ///   The graph <i>must</i> be fill complete.
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    explicit CrsMatrix (const Teuchos::RCP<const crs_graph_type>& graph,
                        const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying column Map and arrays containing
    ///   the matrix in sorted local indices.
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param colMap [in] Distribution of columns of the matrix.
    ///
    /// \param rowPointers [in] The beginning of each row in the matrix,
    ///   as in a CSR "rowptr" array.  The length of this vector should be
    ///   equal to the number of rows in the graph, plus one.  This last
    ///   entry should store the nunber of nonzeros in the matrix.
    ///
    /// \param columnIndices [in] The local indices of the columns,
    ///   as in a CSR "colind" array.  The length of this vector
    ///   should be equal to the number of unknowns in the matrix.
    ///
    /// \param values [in] The local entries in the matrix,
    ///   as in a CSR "vals" array.  The length of this vector
    ///   should be equal to the number of unknowns in the matrix.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap,
               const typename local_matrix_type::row_map_type& rowPointers,
               const t_LocalOrdinal_1D & columnIndices,
               const t_ValuesType & values,
               const Teuchos::RCP<Teuchos::ParameterList>& params = null);

    /// \brief Constructor specifying column Map and arrays containing
    ///   the matrix in sorted, local ids.
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param colMap [in] Distribution of columns of the matrix.
    ///
    /// \param rowPointers [in] The beginning of each row in the matrix,
    ///   as in a CSR "rowptr" array.  The length of this vector should be
    ///   equal to the number of rows in the graph, plus one.  This last
    ///   entry should store the nunber of nonzeros in the matrix.
    ///
    /// \param columnIndices [in] The local indices of the columns,
    ///   as in a CSR "colind" array.  The length of this vector
    ///   should be equal to the number of unknowns in the matrix.
    ///
    /// \param values [in] The local entries in the matrix,
    ///   as in a CSR "vals" array.  The length of this vector
    ///   should be equal to the number of unknowns in the matrix.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap,
               const Teuchos::ArrayRCP<size_t>& rowPointers,
               const Teuchos::ArrayRCP<LocalOrdinal>& columnIndices,
               const Teuchos::ArrayRCP<Scalar>& values,
               const Teuchos::RCP<ParameterList>& params = null);

    /// \brief Constructor specifying column Map and a local matrix,
    ///   which the resulting CrsMatrix views.
    ///
    /// Unlike most other CrsMatrix constructors, successful
    /// completion of this constructor will result in a fill-complete
    /// matrix.
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param colMap [in] Distribution of columns of the matrix.
    ///
    /// \param lclMatrix [in] A local CrsMatrix containing all local
    ///    matrix values as well as a local graph.  The graph's local
    ///    row indices must come from the specified row Map, and its
    ///    local column indices must come from the specified column
    ///    Map.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap,
               const local_matrix_type& lclMatrix,
               const Teuchos::RCP<Teuchos::ParameterList>& params = null);

    // This friend declaration makes the clone() method work.
    template <class S2, class LO2, class GO2, class N2, const bool isClassic>
    friend class CrsMatrix;

    /// \brief Create a deep copy of this CrsMatrix, where the copy
    ///   may have a different Node type.
    ///
    /// \param node2 [in] Kokkos Node instance for the returned copy.
    /// \param params [in/out] Optional list of parameters. If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    ///
    /// Parameters to \c params:
    /// - "Static profile clone" [boolean, default: true] If \c true,
    ///   create the copy with a static allocation profile. If false,
    ///   use a dynamic allocation profile.
    /// - "Locally indexed clone" [boolean] If \c true, fill clone
    ///   using this matrix's column Map and local indices.  This
    ///   matrix must have a column Map in order for this to work.  If
    ///   false, fill clone using global indices.  By default, this
    ///   will use local indices only if this matrix is using local
    ///   indices.
    /// - "fillComplete clone" [boolean, default: true] If \c true,
    ///   call fillComplete() on the cloned CrsMatrix object, with
    ///   parameters from the input parameters' "CrsMatrix" sublist
    ///   The domain Map and range Map passed to fillComplete() are
    ///   those of the map being cloned, if they exist. Otherwise, the
    ///   row Map is used.
    template <class Node2>
    Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node2> >
    clone (const Teuchos::RCP<Node2>& node2,
           const Teuchos::RCP<Teuchos::ParameterList>& params = null) const
    {
      using Teuchos::ArrayRCP;
      using Teuchos::null;
      using Teuchos::ParameterList;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::sublist;
      typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node2> CrsMatrix2;
      typedef Map<LocalOrdinal, GlobalOrdinal, Node2> Map2;
      const char tfecfFuncName[] = "clone";

      // Get parameter values.  Set them initially to their default values.
      bool fillCompleteClone = true;
      bool useLocalIndices = this->hasColMap ();
      ProfileType pftype = StaticProfile;
      if (! params.is_null ()) {
        fillCompleteClone = params->get ("fillComplete clone", fillCompleteClone);
        useLocalIndices = params->get ("Locally indexed clone", useLocalIndices);

        bool staticProfileClone = true;
        staticProfileClone = params->get ("Static profile clone", staticProfileClone);
        pftype = staticProfileClone ? StaticProfile : DynamicProfile;
      }

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        ! this->hasColMap () && useLocalIndices, std::runtime_error,
        ": You requested that the returned clone have local indices, but the "
        "the source matrix does not have a column Map yet.");

      RCP<const Map2> clonedRowMap = this->getRowMap ()->template clone<Node2> (node2);

      // Get an upper bound on the number of entries per row.
      RCP<CrsMatrix2> clonedMatrix;
      ArrayRCP<const size_t> numEntriesPerRow;
      size_t numEntriesForAll = 0;
      bool boundSameForAllLocalRows = false;
      staticGraph_->getNumEntriesPerLocalRowUpperBound (numEntriesPerRow,
                                                        numEntriesForAll,
                                                        boundSameForAllLocalRows);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        numEntriesForAll != 0 &&
        static_cast<size_t> (numEntriesPerRow.size ()) != 0,
        std::logic_error, ": getNumEntriesPerLocalRowUpperBound returned a "
        "nonzero numEntriesForAll = " << numEntriesForAll << " , as well as a "
        "numEntriesPerRow array of nonzero length " << numEntriesPerRow.size ()
        << ".  This should never happen.  Please report this bug to the Tpetra "
        "developers.");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        numEntriesForAll != 0 && ! boundSameForAllLocalRows,
        std::logic_error, ": getNumEntriesPerLocalRowUpperBound returned a "
        "nonzero numEntriesForAll = " << numEntriesForAll << " , but claims "
        "(via its third output value) that the upper bound is not the same for "
        "all rows.  This should never happen.  Please report this bug to the "
        "Tpetra developers.");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        numEntriesPerRow.size () != 0 && boundSameForAllLocalRows,
        std::logic_error, ": getNumEntriesPerLocalRowUpperBound returned a "
        "numEntriesPerRow array of nonzero length " << numEntriesPerRow.size ()
        << ", but claims (via its third output value) that the upper bound is "
        "not the same for all rows.  This should never happen.  Please report "
        "this bug to the Tpetra developers.");

      RCP<ParameterList> matParams =
        params.is_null () ? null : sublist (params,"CrsMatrix");
      if (useLocalIndices) {
        RCP<const Map2> clonedColMap =
          this->getColMap ()->template clone<Node2> (node2);
        if (numEntriesPerRow.is_null ()) {
          clonedMatrix = rcp (new CrsMatrix2 (clonedRowMap, clonedColMap,
                                              numEntriesForAll, pftype,
                                              matParams));
        }
        else {
          clonedMatrix = rcp (new CrsMatrix2 (clonedRowMap, clonedColMap,
                                              numEntriesPerRow, pftype,
                                              matParams));
        }
      }
      else {
        if (numEntriesPerRow.is_null ()) {
          clonedMatrix = rcp (new CrsMatrix2 (clonedRowMap, numEntriesForAll,
                                              pftype, matParams));
        }
        else {
          clonedMatrix = rcp (new CrsMatrix2 (clonedRowMap, numEntriesPerRow,
                                              pftype, matParams));
        }
      }
      // done with these
      numEntriesPerRow = Teuchos::null;
      numEntriesForAll = 0;

      if (useLocalIndices) {
        clonedMatrix->allocateValues (LocalIndices,
                                      CrsMatrix2::GraphNotYetAllocated);
        if (this->isLocallyIndexed ()) {
          ArrayView<const LocalOrdinal> linds;
          ArrayView<const Scalar> vals;
          for (LocalOrdinal lrow = clonedRowMap->getMinLocalIndex ();
               lrow <= clonedRowMap->getMaxLocalIndex ();
               ++lrow) {
            this->getLocalRowView (lrow, linds, vals);
            if (linds.size ()) {
              clonedMatrix->insertLocalValues (lrow, linds, vals);
            }
          }
        }
        else { // this->isGloballyIndexed()
          Array<LocalOrdinal> linds;
          Array<Scalar> vals;
          for (LocalOrdinal lrow = clonedRowMap->getMinLocalIndex ();
               lrow <= clonedRowMap->getMaxLocalIndex ();
               ++lrow) {
            size_t theNumEntries = this->getNumEntriesInLocalRow (lrow);
            if (theNumEntries > static_cast<size_t> (linds.size ())) {
              linds.resize (theNumEntries);
            }
            if (theNumEntries > static_cast<size_t> (vals.size ())) {
              vals.resize (theNumEntries);
            }
            this->getLocalRowCopy (clonedRowMap->getGlobalElement (lrow),
                                   linds (), vals (), theNumEntries);
            if (theNumEntries != 0) {
              clonedMatrix->insertLocalValues (lrow, linds (0, theNumEntries),
                                               vals (0, theNumEntries));
            }
          }
        }
      }
      else { // useGlobalIndices
        clonedMatrix->allocateValues (GlobalIndices,
                                      CrsMatrix2::GraphNotYetAllocated);
        if (this->isGloballyIndexed ()) {
          ArrayView<const GlobalOrdinal> ginds;
          ArrayView<const Scalar> vals;
          for (GlobalOrdinal grow = clonedRowMap->getMinGlobalIndex ();
               grow <= clonedRowMap->getMaxGlobalIndex ();
               ++grow) {
            this->getGlobalRowView (grow, ginds, vals);
            if (ginds.size () > 0) {
              clonedMatrix->insertGlobalValues (grow, ginds, vals);
            }
          }
        }
        else { // this->isLocallyIndexed()
          Array<GlobalOrdinal> ginds;
          Array<Scalar> vals;
          for (GlobalOrdinal grow = clonedRowMap->getMinGlobalIndex ();
               grow <= clonedRowMap->getMaxGlobalIndex ();
               ++grow) {
            size_t theNumEntries = this->getNumEntriesInGlobalRow (grow);
            if (theNumEntries > static_cast<size_t> (ginds.size ())) {
              ginds.resize (theNumEntries);
            }
            if (theNumEntries > static_cast<size_t> (vals.size ())) {
              vals.resize (theNumEntries);
            }
            this->getGlobalRowCopy (grow, ginds (), vals (), theNumEntries);
            if (theNumEntries != 0) {
              clonedMatrix->insertGlobalValues (grow, ginds (0, theNumEntries),
                                                vals (0, theNumEntries));
            }
          }
        }
      }

      if (fillCompleteClone) {
        RCP<ParameterList> fillparams =
          params.is_null () ? Teuchos::null : sublist (params, "fillComplete");
        try {
          RCP<const Map2> clonedRangeMap;
          RCP<const Map2> clonedDomainMap;
          if (! this->getRangeMap ().is_null () &&
              this->getRangeMap () != clonedRowMap) {
            clonedRangeMap  = this->getRangeMap ()->template clone<Node2> (node2);
          }
          else {
            clonedRangeMap = clonedRowMap;
          }
          if (! this->getDomainMap ().is_null () &&
              this->getDomainMap () != clonedRowMap) {
            clonedDomainMap = this->getDomainMap ()->template clone<Node2> (node2);
          }
          else {
            clonedDomainMap = clonedRowMap;
          }
          clonedMatrix->fillComplete (clonedDomainMap, clonedRangeMap,
                                      fillparams);
        }
        catch (std::exception &e) {
          const bool caughtExceptionOnClone = true;
          TEUCHOS_TEST_FOR_EXCEPTION(
            caughtExceptionOnClone, std::runtime_error,
            Teuchos::typeName (*this) << std::endl << "clone: " << std::endl <<
            "Caught the following exception while calling fillComplete() on a "
            "clone of type" << std::endl << Teuchos::typeName (*clonedMatrix)
            << ": " << std::endl << e.what () << std::endl);
        }
      }
      return clonedMatrix;
    }

    //! Destructor.
    virtual ~CrsMatrix ();

    //@}
    //! @name Methods for inserting, modifying, or removing entries
    //@{

    /// \brief Insert one or more entries into the matrix, using global indices.
    ///
    /// \param globalRow [in] Global index of the row into which to
    ///   insert the entries.
    /// \param cols [in] Global indices of the columns into which
    ///   to insert the entries.
    /// \param values [in] Values to insert into the above columns.
    ///
    /// For all k in 0, ..., <tt>col.size()-1</tt>, insert the value
    /// <tt>values[k]</tt> into entry <tt>(globalRow, cols[k])</tt> of
    /// the matrix.  If that entry already exists, add the new value
    /// to the old value.
    ///
    /// This is a local operation.  It does not communicate (using
    /// MPI).  If row \c globalRow is owned by the calling process,
    /// the entries will be inserted immediately.  Otherwise, if that
    /// row is <i>not</i> owned by the calling process, then the
    /// entries will be stored locally for now, and only communicated
    /// to the process that owns the row when either fillComplete() or
    /// globalAssemble() is called.  If that process already has an
    /// entry, the incoming value will be added to the old value, just
    /// as if it were inserted on the owning process.
    //
    /// If the matrix has a column Map (<tt>hasColMap() == true</tt>),
    /// and if globalRow is owned by process p, then it is forbidden
    /// to insert column indices that are not in the column Map on
    /// process p.  Tpetra will test the input column indices to
    /// ensure this is the case, but if \c globalRow is not owned by
    /// the calling process, the test will be deferred until the next
    /// call to globalAssemble() or fillComplete().
    ///
    /// \warning The behavior described in the above paragraph differs
    ///   from that of Epetra.  If the matrix has a column Map,
    ///   Epetra_CrsMatrix "filters" column indices not in the column
    ///   Map.  Many users found this confusing, so we changed it so
    ///   that nonowned column indices are forbidden.
    ///
    /// It is legal to call this method whether the matrix's column
    /// indices are globally or locally indexed.  If the matrix's
    /// column indices are locally indexed (<tt>isLocallyIndexed() ==
    /// true</tt>), then this method will convert the input global
    /// column indices to local column indices.
    ///
    /// For better performance when filling entries into a sparse
    /// matrix, consider the following tips:
    /// <ol>
    /// <li>Use local indices (e.g., insertLocalValues()) if you know
    ///   the column Map in advance.  Converting global indices to
    ///   local indices is expensive.  Of course, if you don't know
    ///   the column Map in advance, you must use global indices.</li>
    /// <li>When invoking the CrsMatrix constructor, give the best
    ///   possible upper bounds on the number of entries in each row
    ///   of the matrix.  This will avoid expensive reallocation if
    ///   your bound was not large enough.</li>
    /// <li>If your upper bound on the number of entries in each row
    ///   will always be correct, create the matrix with
    ///   StaticProfile.  This uses a faster and more compact data
    ///   structure to store the matrix.</li>
    /// <li>If you plan to reuse a matrix's graph structure, but
    ///   change its values, in repeated fillComplete() / resumeFill()
    ///   cycles, you can get the best performance by creating the
    ///   matrix with a const CrsGraph.  Do this by using the
    ///   CrsMatrix constructor that accepts an RCP of a const
    ///   CrsGraph.  If you do this, you must use the "replace" or
    ///   "sumInto" methods to change the values of the matrix; you
    ///   may not use insertGlobalValues() or
    ///   insertLocalValues().</li>
    /// </ol>
    void
    insertGlobalValues (const GlobalOrdinal globalRow,
                        const Teuchos::ArrayView<const GlobalOrdinal>& cols,
                        const Teuchos::ArrayView<const Scalar>& vals);

    /// \brief Insert one or more entries into the matrix, using local indices.
    ///
    /// \param LocalRow [in] Local index of the row into which to
    ///   insert the entries.  It must be owned by the row Map on the
    ///   calling process.
    /// \param cols [in] Local indices of the columns into which to
    ///   insert the entries.  All of the column indices must be owned
    ///   by the column Map on the calling process.
    /// \param values [in] Values to insert into the above columns.
    ///
    /// For all k in 0, ..., <tt>cols.size()-1</tt>, insert the value
    /// <tt>values[k]</tt> into entry <tt>(globalRow, cols[k])</tt> of
    /// the matrix.  If that entry already exists, add the new value
    /// to the old value.
    ///
    /// In order to call this method, the matrix must be locally
    /// indexed, and it must have a column Map.
    ///
    /// For better performance when filling entries into a sparse
    /// matrix, consider the following tips:
    /// <ol>
    /// <li>When invoking the CrsMatrix constructor, give the best
    ///   possible upper bounds on the number of entries in each row
    ///   of the matrix.  This will avoid expensive reallocation if
    ///   your bound was not large enough.</li>
    /// <li>If your upper bound on the number of entries in each row
    ///   will always be correct, create the matrix with
    ///   StaticProfile.  This uses a faster and more compact data
    ///   structure to store the matrix.</li>
    /// <li>If you plan to reuse a matrix's graph structure, but
    ///   change its values, in repeated fillComplete() / resumeFill()
    ///   cycles, you can get the best performance by creating the
    ///   matrix with a const CrsGraph.  Do this by using the
    ///   CrsMatrix constructor that accepts an RCP of a const
    ///   CrsGraph.  If you do this, you must use the "replace" or
    ///   "sumInto" methods to change the values of the matrix; you
    ///   may not use insertGlobalValues() or
    ///   insertLocalValues().</li>
    /// </ol>
    void
    insertLocalValues (const LocalOrdinal localRow,
                       const ArrayView<const LocalOrdinal> &cols,
                       const ArrayView<const Scalar> &vals);

    /// \brief Replace one or more entries' values, using global indices.
    ///
    /// \param globalRow [in] Global index of the row in which to
    ///   replace the entries.  This row <i>must</i> be owned by the
    ///   calling process.
    /// \param cols [in] Global indices of the columns in which to
    ///   replace the entries.
    /// \param vals [in] Values to use for replacing the entries.
    ///
    /// For all k in 0, ..., <tt>cols.size()-1</tt>, replace the value
    /// at entry <tt>(globalRow, cols[k])</tt> of the matrix with
    /// <tt>vals[k]</tt>.  That entry must exist in the matrix
    /// already.
    ///
    /// If <tt>(globalRow, cols[k])</tt> corresponds to an entry that
    /// is duplicated in this matrix row (likely because it was
    /// inserted more than once and fillComplete() has not been called
    /// in the interim), the behavior of this method is not defined.
    ///
    /// \return The number of indices for which values were actually
    ///   replaced; the number of "correct" indices.
    ///
    /// If the returned value N satisfies
    ///
    /// <tt>0 <= N < cols.size()</tt>,
    ///
    /// then <tt>cols.size() - N</tt> of the entries of <tt>cols</tt>
    /// are not valid global column indices.  If the returned value is
    /// Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), then at least
    /// one of the following is true:
    ///   <ul>
    ///   <li> <tt>! isFillActive ()</tt> </li>
    ///   <li> <tt> cols.size () != vals.size ()</tt> </li>
    ///   </ul>
    LocalOrdinal
    replaceGlobalValues (GlobalOrdinal globalRow,
                         const ArrayView<const GlobalOrdinal>& cols,
                         const ArrayView<const Scalar>& vals);

    /// \brief Replace one or more entries' values, using local indices.
    ///
    /// \param localRow [in] local index of the row in which to
    ///   replace the entries.  This row <i>must</i> be owned by the
    ///   calling process.
    /// \param cols [in] Local indices of the columns in which to
    ///   replace the entries.
    /// \param vals [in] Values to use for replacing the entries.
    ///
    /// For all k in 0, ..., <tt>cols.size()-1</tt>, replace the value
    /// at entry <tt>(localRow, cols[k])</tt> of the matrix with
    /// <tt>vals[k]</tt>.  That entry must exist in the matrix
    /// already.
    ///
    /// \return The number of indices for which values were actually
    ///   replaced; the number of "correct" indices.
    ///
    /// If the returned value N satisfies
    ///
    /// <tt>0 <= N < cols.size()</tt>,
    ///
    /// then <tt>cols.size() - N</tt> of the entries of <tt>cols</tt>
    /// are not valid local column indices.  If the returned value is
    /// Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), then at least
    /// one of the following is true:
    ///   <ul>
    ///   <li> <tt>! isFillActive ()</tt> </li>
    ///   <li> <tt>! hasColMap ()</tt> </li>
    ///   <li> <tt> cols.size () != vals.size ()</tt> </li>
    ///   </ul>
    LocalOrdinal
    replaceLocalValues (const LocalOrdinal localRow,
                        const ArrayView<const LocalOrdinal>& cols,
                        const ArrayView<const Scalar>& vals);

    /// \brief Sum into one or more sparse matrix entries, using global indices.
    ///
    /// This is a local operation; it does not involve communication.
    /// However, if you sum into rows not owned by the calling
    /// process, it may result in future communication in
    /// globalAssemble() (which is called by fillComplete()).
    ///
    /// If \c globalRow is owned by the calling process, then this
    /// method performs the sum-into operation right away.  Otherwise,
    /// if the row is <i>not</i> owned by the calling process, this
    /// method defers the sum-into operation until globalAssemble().
    /// That method communicates data for nonowned rows to the
    /// processes that own those rows.  Then, globalAssemble() does
    /// one of the following:
    /// - It calls insertGlobalValues() for that data if the matrix
    ///   has a dynamic graph.
    /// - It calls sumIntoGlobalValues() for that data if the matrix
    ///   has a static graph.  The matrix silently ignores
    ///   (row,column) pairs that do not exist in the graph.
    ///
    /// \param globalRow [in] The global index of the row in which to
    ///   sum into the matrix entries.
    /// \param cols [in] One or more column indices.
    /// \param vals [in] One or more values corresponding to those
    ///   column indices.  <tt>vals[k]</tt> corresponds to
    ///   <tt>cols[k]</tt>.
    ///
    /// \return The number of indices for which values were actually
    ///   modified; the number of "correct" indices.
    ///
    /// This method has the same preconditions and return value
    /// meaning as replaceGlobalValues() (which see).
    LocalOrdinal
    sumIntoGlobalValues (const GlobalOrdinal globalRow,
                         const ArrayView<const GlobalOrdinal> &cols,
                         const ArrayView<const Scalar>        &vals);

    /// \brief Sum into one or more sparse matrix entries, using local indices.
    ///
    /// \param localRow [in] Local index of a row.  This row
    ///   <i>must</i> be owned by the calling process.
    /// \param cols [in] Local indices of the columns whose entries we
    ///   want to modify.
    /// \param vals [in] Values corresponding to the above column
    ///   indices.  <tt>vals[k]</tt> corresponds to <tt>cols[k]</tt>.
    ///
    /// \return The number of indices for which values were actually
    ///   modified; the number of "correct" indices.
    ///
    /// This method has the same preconditions and return value
    /// meaning as replaceLocalValues() (which see).
    LocalOrdinal
    sumIntoLocalValues (const LocalOrdinal localRow,
                        const ArrayView<const LocalOrdinal>& cols,
                        const ArrayView<const Scalar>& vals);

    //! Set all matrix entries equal to \c alpha.
    void setAllToScalar (const Scalar &alpha);

    //! Scale the matrix's values: <tt>this := alpha*this</tt>.
    void scale (const Scalar &alpha);

    //! Sets the 1D pointer arrays of the graph.
    /**
       \pre <tt>hasColMap() == true</tt>
       \pre <tt>getGraph() != Teuchos::null</tt>
       \pre No insert/sum routines have been called

       \warning This method is intended for expert developer use only, and should never be called by user code.
    */
    void
    setAllValues (const typename local_matrix_type::row_map_type& rowPointers,
                  const t_LocalOrdinal_1D& columnIndices,
                  const t_ValuesType& values);

    //! Sets the 1D pointer arrays of the graph.
    /**
       \pre <tt>hasColMap() == true</tt>
       \pre <tt>getGraph() != Teuchos::null</tt>
       \pre No insert/sum routines have been called

       FIXME (mfh 24 Feb 2014) Why is the third prerequisites above
       different than the third prerequisite from the original class?
       The original is that fillComplete() must have been called.

       \warning This method is intended for expert developer use only, and should never be called by user code.
    */
    void
    setAllValues (const Teuchos::ArrayRCP<size_t>& rowPointers,
                  const Teuchos::ArrayRCP<LocalOrdinal>& columnIndices,
                  const Teuchos::ArrayRCP<Scalar>& values);

    void
    getAllValues (Teuchos::ArrayRCP<const size_t>& rowPointers,
                  Teuchos::ArrayRCP<const LocalOrdinal>& columnIndices,
                  Teuchos::ArrayRCP<const Scalar>& values) const;

    //@}
    //! @name Transformational methods
    //@{

    /// \brief Communicate nonlocal contributions to other processes.
    ///
    /// Users do not normally need to call this method.  fillComplete
    /// always calls this method, unless you specifically tell
    /// fillComplete to do otherwise by setting its "No Nonlocal
    /// Changes" parameter to \c true.  Thus, it suffices to call
    /// fillComplete.
    ///
    /// Methods like insertGlobalValues and sumIntoGlobalValues let
    /// you add or modify entries in rows that are not owned by the
    /// calling process.  These entries are called "nonlocal
    /// contributions."  The methods that allow nonlocal contributions
    /// store the entries on the calling process, until globalAssemble
    /// is called.  globalAssemble sends these nonlocal contributions
    /// to the process(es) that own them, where they then become part
    /// of the matrix.
    ///
    /// This method only does global assembly if there are nonlocal
    /// entries on at least one process.  It does an all-reduce to
    /// find that out.  If not, it returns early, without doing any
    /// more communication or work.
    ///
    /// If you previously inserted into a row which is not owned by
    /// <i>any</i> process in the row Map, the behavior of this method
    /// is undefined.  It may detect the invalid row indices and throw
    /// an exception, or it may silently drop the entries inserted
    /// into invalid rows.  Behavior may vary, depending on whether
    /// Tpetra was built with debug checking enabled.
    void globalAssemble();

    /// \brief Resume operations that may change the values or
    ///   structure of the matrix.
    ///
    /// This method must be called as a collective operation.
    ///
    /// Calling fillComplete "freezes" both the values and the
    /// structure of the matrix.  If you want to modify the matrix
    /// again, you must first call resumeFill.  You then may not call
    /// resumeFill again on that matrix until you first call
    /// fillComplete.  You may make sequences of fillComplete,
    /// resumeFill calls as many times as you wish.
    ///
    /// \post <tt>isFillActive() && ! isFillComplete()</tt>
    void resumeFill (const RCP<ParameterList>& params = null);

    /*! \brief Signal that data entry is complete, specifying domain and range maps.

      Off-node indices are distributed (via globalAssemble()), indices are sorted, redundant indices are eliminated, and global indices are transformed to local indices.

      \pre  <tt>isFillActive() == true<tt>
      \pre <tt>isFillComplete()() == false<tt>

      \post <tt>isFillActive() == false<tt>
      \post <tt>isFillComplete() == true<tt>

      Parameters:

      - "No Nonlocal Changes" (\c bool): Default is false.  If true,
        the caller promises that no modifications to nonowned rows
        have happened on any process since the last call to
        fillComplete.  This saves a global all-reduce to check whether
        any process did a nonlocal insert.  Nonlocal changes include
        any sumIntoGlobalValues or insertGlobalValues call with a row
        index that is not in the row Map of the calling process.

      - "Sort column Map ghost GIDs" (\c bool): Default is true.
        makeColMap() (which fillComplete may call) always groups
        remote GIDs by process rank, so that all remote GIDs with the
        same owning rank occur contiguously.  By default, it always
        sorts remote GIDs in increasing order within those groups.
        This behavior differs from Epetra, which does not sort remote
        GIDs with the same owning process.  If you don't want to sort
        (for compatibility with Epetra), set this parameter to \c
        false.  This parameter only takes effect if the matrix owns
        the graph.  This is an expert mode parameter ONLY.  We make no
        promises about backwards compatibility of this parameter.  It
        may change or disappear at any time.
    */
    void
    fillComplete (const RCP<const map_type>& domainMap,
                  const RCP<const map_type>& rangeMap,
                  const RCP<ParameterList>& params = null);

    /*! \brief Signal that data entry is complete.

      Off-node entries are distributed (via globalAssemble()), repeated entries are summed, and global indices are transformed to local indices.

      \note This method calls fillComplete( getRowMap(), getRowMap(), os ).

      \pre  <tt>isFillActive() == true<tt>
      \pre <tt>isFillComplete()() == false<tt>

      \post <tt>isFillActive() == false<tt>
      \post <tt>isFillComplete() == true<tt>
    */
    void fillComplete (const RCP<ParameterList>& params = null);

    /// \brief Perform a fillComplete on a matrix that already has data.
    ///
    /// The matrix must already have filled local 1-D storage
    /// (k_clInds1D_ and k_rowPtrs_ for the graph, and k_values1D_ in
    /// the matrix).  If the matrix has been constructed in any other
    /// way, this method will throw an exception.  This routine is
    /// needed to support other Trilinos packages and should not be
    /// called by ordinary users.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    void
    expertStaticFillComplete (const RCP<const map_type>& domainMap,
                              const RCP<const map_type>& rangeMap,
                              const RCP<const import_type>& importer = Teuchos::null,
                              const RCP<const export_type>& exporter = Teuchos::null,
                              const RCP<ParameterList>& params = Teuchos::null);

    /// \brief Replace the matrix's column Map with the given Map.
    ///
    /// \param newColMap [in] New column Map.  Must be nonnull.
    ///
    /// \pre The matrix must have no entries inserted yet, on any
    ///   process in the row Map's communicator.
    ///
    /// \pre The matrix must not have been created with a constant
    ///   (a.k.a. "static") CrsGraph.
    void
    replaceColMap (const Teuchos::RCP<const map_type>& newColMap);

    /// \brief Reindex the column indices in place, and replace the
    ///   column Map.  Optionally, replace the Import object as well.
    ///
    /// \pre The matrix is <i>not</i> fill complete:
    ///   <tt>! this->isFillComplete() </tt>.
    /// \pre Either the input graph is \c NULL, or it is <i>not</i>
    ///   fill complete:
    ///   <tt>graph == NULL || ! graph->isFillComplete()</tt>.
    /// \pre On every calling process, every index owned by the
    ///   current column Map must also be owned by the new column Map.
    /// \pre If the new Import object is provided, the new Import
    ///   object's source Map must be the same as the current domain
    ///   Map, and the new Import's target Map must be the same as the
    ///   new column Map.
    ///
    /// \param graph [in] The matrix's graph.  If you don't provide
    ///   this (i.e., if <tt>graph == NULL</tt>), then the matrix must
    ///   own its graph, which will be modified in place.  (That is,
    ///   you must <i>not</i> have created the matrix with a constant
    ///   graph.)  If you <i>do</i> provide this, then the method will
    ///   assume that it is the same graph as the matrix's graph, and
    ///   the provided graph will be modified in place.
    /// \param newColMap [in] New column Map.  Must be nonnull.
    /// \param newImport [in] New Import object.  Optional; computed
    ///   if not provided or if null.  Computing an Import is
    ///   expensive, so it is worth providing this if you can.
    /// \param sortEachRow [in] If true, sort the indices (and their
    ///   corresponding values) in each row after reindexing.
    ///
    /// Why would you want to use this method?  Well, for example, you
    /// might need to use an Ifpack2 preconditioner that only accepts
    /// a matrix with a certain kind of column Map.  Your matrix has
    /// the wrong kind of column Map, but you know how to compute the
    /// right kind of column Map.  You might also know an efficient
    /// way to compute an Import object from the current domain Map to
    /// the new column Map.  (For an instance of the latter, see the
    /// Details::makeOptimizedColMapAndImport function in
    /// Tpetra_Details_makeOptimizedColMap.hpp.)
    ///
    /// Suppose that you created this CrsMatrix with a constant graph;
    /// that is, that you called the CrsMatrix constructor that takes
    /// a CrsGraph as input:
    ///
    /// \code
    /// RCP<CrsGraph<> > G (new CrsGraph<> (rowMap, origColMap, ...));
    /// // ... fill G ...
    /// G->fillComplete (domMap, ranMap);
    /// CrsMatrix<> A (G);
    /// // ... fill A ...
    /// \endcode
    ///
    /// Now suppose that you want to give A to a preconditioner that
    /// can't handle a matrix with an arbitrary column Map (in the
    /// example above, <tt>origColMap</tt>).  You first must create a
    /// new suitable column Map <tt>newColMap</tt>, and optionally a
    /// new Import object <tt>newImport</tt> from the matrix's current
    /// domain Map to the new column Map.  Then, call this method,
    /// passing in G (which must <i>not</i> be fill complete) while
    /// the matrix is <i>not</i> fill complete.  Be sure to save the
    /// graph's <i>original</i> Import object; you'll need that later.
    ///
    /// \code
    /// RCP<const CrsGraph<>::import_type> origImport = G->getImporter ();
    /// G->resumeFill ();
    /// A.reindexColumns (G.getRawPtr (), newColMap, newImport);
    /// G.fillComplete (domMap, ranMap);
    /// A.fillComplete (domMap, ranMap);
    /// \endcode
    ///
    /// Now you may give the matrix A to the preconditioner in
    /// question.  After doing so, and after you solve the linear
    /// system using the preconditioner, you might want to put the
    /// matrix back like it originally was.  You can do that, too!
    ///
    /// \code
    /// A.resumeFill ();
    /// G->resumeFill ();
    /// A.reindexColumns (G.getRawPtr (), origColMap, origImport);
    /// G->fillComplete (domMap, ranMap);
    /// A->fillComplete (domMap, ranMap);
    /// \endcode
    void
    reindexColumns (crs_graph_type* const graph,
                    const Teuchos::RCP<const map_type>& newColMap,
                    const Teuchos::RCP<const import_type>& newImport = Teuchos::null,
                    const bool sortEachRow = true);

    /// \brief Replace the current domain Map and Import with the given objects.
    ///
    /// \param newDomainMap [in] New domain Map.  Must be nonnull.
    /// \param newImporter [in] Optional Import object.  If null, we
    ///   will compute it.
    ///
    /// \pre The matrix must be fill complete:
    ///   <tt>isFillComplete() == true</tt>.
    /// \pre If the Import is provided, its target Map must be the
    ///   same as the column Map of the matrix.
    /// \pre If the Import is provided, its source Map must be the
    ///   same as the provided new domain Map.
    void
    replaceDomainMapAndImporter (const Teuchos::RCP<const map_type>& newDomainMap,
                                 Teuchos::RCP<const import_type>& newImporter);

    /// \brief Remove processes owning zero rows from the Maps and their communicator.
    ///
    /// \warning This method is ONLY for use by experts.  We highly
    ///   recommend using the nonmember function of the same name
    ///   defined in Tpetra_DistObject_decl.hpp.
    ///
    /// \warning We make NO promises of backwards compatibility.
    ///   This method may change or disappear at any time.
    ///
    /// \param newMap [in] This <i>must</i> be the result of calling
    ///   the removeEmptyProcesses() method on the row Map.  If it
    ///   is not, this method's behavior is undefined.  This pointer
    ///   will be null on excluded processes.
    virtual void
    removeEmptyProcessesInPlace (const Teuchos::RCP<const map_type>& newMap);

    //@}
    //! @name Methods implementing RowMatrix
    //@{

    //! The communicator over which the matrix is distributed.
    Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

    //! The Kokkos Node instance.
    Teuchos::RCP<node_type> getNode () const;

    //! The Map that describes the row distribution in this matrix.
    Teuchos::RCP<const map_type> getRowMap () const;

    //! The Map that describes the column distribution in this matrix.
    Teuchos::RCP<const map_type> getColMap () const;

    //! This matrix's graph, as a RowGraph.
    Teuchos::RCP<const RowGraph<LocalOrdinal,GlobalOrdinal,node_type> > getGraph () const;

    //! This matrix's graph, as a CrsGraph.
    Teuchos::RCP<const crs_graph_type> getCrsGraph () const;

    //! The local sparse matrix.
    local_matrix_type getLocalMatrix () const {return lclMatrix_; }

    /// \brief Number of global elements in the row map of this matrix.
    ///
    /// This is <it>not</it> the number of rows in the matrix as a
    /// mathematical object.  This method returns the global sum of
    /// the number of local elements in the row map on each processor,
    /// which is the row map's getGlobalNumElements().  Since the row
    /// map is not one-to-one in general, that global sum could be
    /// different than the number of rows in the matrix.  If you want
    /// the number of rows in the matrix, ask the range map for its
    /// global number of elements, using the following code:
    /// <code>
    /// global_size_t globalNumRows = getRangeMap()->getGlobalNumElements();
    /// </code>
    /// This method retains the behavior of Epetra, which also asks
    /// the row map for the global number of rows, rather than asking
    /// the range map.
    ///
    /// \warning Undefined if isFillActive().
    ///
    global_size_t getGlobalNumRows() const;

    /// \brief The number of global columns in the matrix.
    ///
    /// This equals the number of entries in the matrix's domain Map.
    ///
    /// \warning Undefined if isFillActive().
    global_size_t getGlobalNumCols() const;

    /// \brief The number of matrix rows owned by the calling process.
    ///
    /// Note that the sum of all the return values over all processes
    /// in the row Map's communicator does not necessarily equal the
    /// global number of rows in the matrix, if the row Map is
    /// overlapping.
    size_t getNodeNumRows() const;

    /// \brief The number of columns connected to the locally owned rows of this matrix.
    ///
    /// Throws std::runtime_error if <tt>! hasColMap ()</tt>.
    size_t getNodeNumCols() const;

    //! The index base for global indices for this matrix.
    GlobalOrdinal getIndexBase() const;

    //! The global number of entries in this matrix.
    global_size_t getGlobalNumEntries() const;

    //! The local number of entries in this matrix.
    size_t getNodeNumEntries() const;

    //! \brief Returns the current number of entries on this node in the specified global row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this matrix. */
    size_t getNumEntriesInGlobalRow (GlobalOrdinal globalRow) const;

    //! Returns the current number of entries on this node in the specified local row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
    size_t getNumEntriesInLocalRow (LocalOrdinal localRow) const;

    //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons.
    /** Undefined if isFillActive().
     */
    global_size_t getGlobalNumDiags() const;

    //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons.
    /** Undefined if isFillActive().
     */
    size_t getNodeNumDiags() const;

    //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
    /** Undefined if isFillActive().
     */
    size_t getGlobalMaxNumRowEntries() const;

    //! \brief Returns the maximum number of entries across all rows/columns on this node.
    /** Undefined if isFillActive().
     */
    size_t getNodeMaxNumRowEntries() const;

    //! \brief Indicates whether the matrix has a well-defined column map.
    bool hasColMap() const;

    //! \brief Indicates whether the matrix is lower triangular.
    /** Undefined if isFillActive().
     */
    bool isLowerTriangular() const;

    //! \brief Indicates whether the matrix is upper triangular.
    /** Undefined if isFillActive().
     */
    bool isUpperTriangular() const;

    /// \brief Whether the matrix is locally indexed on the calling process.
    ///
    /// The matrix is locally indexed on the calling process if and
    /// only if all of the following hold:
    /// <ol>
    /// <li> The matrix is not empty on the calling process </li>
    /// <li> The matrix has a column Map </li>
    /// </ol>
    ///
    /// The following is always true:
    /// \code
    /// (! locallyIndexed() && ! globallyIndexed()) || (locallyIndexed() || globallyIndexed());
    /// \endcode
    /// That is, a matrix may be neither locally nor globally indexed,
    /// but it can never be both.  Furthermore a matrix that is not
    /// fill complete, might have some processes that are neither
    /// locally nor globally indexed, and some processes that are
    /// globally indexed.  The processes that are neither do not have
    /// any entries.
    bool isLocallyIndexed() const;

    /// \brief Whether the matrix is globally indexed on the calling process.
    ///
    /// The matrix is globally indexed on the calling process if and
    /// only if all of the following hold:
    /// <ol>
    /// <li> The matrix is not empty on the calling process </li>
    /// <li> The matrix does not yet have a column Map </li>
    /// </ol>
    ///
    /// The following is always true:
    /// \code
    /// (! locallyIndexed() && ! globallyIndexed()) || (locallyIndexed() || globallyIndexed());
    /// \endcode
    /// That is, a matrix may be neither locally nor globally indexed,
    /// but it can never be both.  Furthermore a matrix that is not
    /// fill complete, might have some processes that are neither
    /// locally nor globally indexed, and some processes that are
    /// globally indexed.  The processes that are neither do not have
    /// any entries.
    bool isGloballyIndexed() const;

    /// \brief Whether the matrix is fill complete.
    ///
    /// A matrix is <i>fill complete</i> (or "in compute mode") when
    /// fillComplete() has been called without an intervening call to
    /// resumeFill().  A matrix must be fill complete in order to call
    /// computational kernels like sparse matrix-vector multiply and
    /// sparse triangular solve.  A matrix must be <i>not</i> fill
    /// complete ("in edit mode") in order to call methods that
    /// insert, modify, or remove entries.
    ///
    /// The following are always true:
    /// <ul>
    /// <li> <tt> isFillActive() == ! isFillComplete() </tt>
    /// <li> <tt> isFillActive() || isFillComplete() </tt>
    /// </ul>
    ///
    /// A matrix starts out (after its constructor returns) as not
    /// fill complete.  It becomes fill complete after fillComplete()
    /// returns, and becomes not fill complete again if resumeFill()
    /// is called.  Some methods like clone() and some of the
    /// "nonmember constructors" (like importAndFillComplete() and
    /// exportAndFillComplete()) may return a fill-complete matrix.
    bool isFillComplete() const;

    /// \brief Whether the matrix is not fill complete.
    ///
    /// A matrix is <i>fill complete</i> (or "in compute mode") when
    /// fillComplete() has been called without an intervening call to
    /// resumeFill().  A matrix must be fill complete in order to call
    /// computational kernels like sparse matrix-vector multiply and
    /// sparse triangular solve.  A matrix must be <i>not</i> fill
    /// complete ("in edit mode") in order to call methods that
    /// insert, modify, or remove entries.
    ///
    /// The following are always true:
    /// <ul>
    /// <li> <tt> isFillActive() == ! isFillComplete() </tt>
    /// <li> <tt> isFillActive() || isFillComplete() </tt>
    /// </ul>
    ///
    /// A matrix starts out (after its constructor returns) as not
    /// fill complete.  It becomes fill complete after fillComplete()
    /// returns, and becomes not fill complete again if resumeFill()
    /// is called.  Some methods like clone() and some of the
    /// "nonmember constructors" (like importAndFillComplete() and
    /// exportAndFillComplete()) may return a fill-complete matrix.
    bool isFillActive() const;

    //! \brief Returns \c true if storage has been optimized.
    /**
       Optimized storage means that the allocation of each row is equal to the
       number of entries. The effect is that a pass through the matrix, i.e.,
       during a mat-vec, requires minimal memory traffic. One limitation of
       optimized storage is that no new indices can be added to the matrix.
    */
    bool isStorageOptimized () const;

    //! Returns \c true if the matrix was allocated with static data structures.
    ProfileType getProfileType () const;

    //! Indicates that the graph is static, so that new entries cannot be added to this matrix.
    bool isStaticGraph () const;

    /// \brief Compute and return the Frobenius norm of the matrix.
    ///
    /// The Frobenius norm of the matrix is defined as
    /// \f\[
    ///   \|A\|_F = \sqrt{\sum_{i,j} \|\a_{ij}\|^2}.
    /// \f\].
    ///
    /// If the matrix is fill complete, then the computed value is
    /// cached; the cache is cleared whenever resumeFill() is called.
    /// Otherwise, the value is computed every time the method is
    /// called.
    mag_type getFrobeniusNorm () const;

    /// \brief Return \c true if getLocalRowView() and
    ///   getGlobalRowView() are valid for this object.
    virtual bool supportsRowViews () const;

    /// \brief Fill given arrays with a deep copy of the locally owned
    ///   entries of the matrix in a given row, using global column
    ///   indices.
    ///
    /// \param GlobalRow [in] Global index of the row for which to
    ///   return entries.
    /// \param Indices [out] Global column indices corresponding to
    ///   values.
    /// \param Values [out] Matrix values.
    /// \param NumEntries [out] Number of entries.
    ///
    /// \note To Tpetra developers: Discussion of whether to use
    ///   <tt>Scalar</tt> or <tt>impl_scalar_type</tt> for output
    ///   array of matrix values.
    ///
    /// If \c Scalar differs from <tt>impl_scalar_type</tt>, as for
    /// example with std::complex<T> and Kokkos::complex<T>, we must
    /// choose which type to use.  We must make the same choice as
    /// RowMatrix does, else CrsMatrix won't compile, because it won't
    /// implement a pure virtual method.  We choose <tt>Scalar</tt>,
    /// for the following reasons.  First, <tt>Scalar</tt> is the
    /// user's preferred type, and <tt>impl_scalar_type</tt> an
    /// implementation detail that makes Tpetra work with Kokkos.
    /// Second, Tpetra's public interface provides a host-only
    /// interface, which eliminates some reasons for requiring
    /// implementation-specific types like Kokkos::complex.
    ///
    /// We do eventually want to put Tpetra methods in Kokkos kernels,
    /// but we only <i>need</i> to put them in host kernels, since
    /// Tpetra is a host-only interface.  Users can still manually
    /// handle conversion from <tt>Scalar</tt> to
    /// <tt>impl_scalar_type</tt> for reductions.
    ///
    /// The right thing to do would be to rewrite RowMatrix so that
    /// getGlobalRowCopy is NOT inherited, but is implemented by a
    /// pure virtual "hook" getGlobalRowCopyImpl.  The latter takes
    /// raw pointers.  That would give us the freedom to overload
    /// getGlobalRowCopy, which one normally can't do with virtual
    /// methods.  It would make sense for one getGlobalRowCopyImpl
    /// method to implement both Teuchos::ArrayView and Kokos::View
    /// versions of getGlobalRowCopy.
    ///
    /// Note: A std::runtime_error exception is thrown if either
    /// <tt>Indices</tt> or <tt>Values</tt> is not large enough to
    /// hold the data associated with row \c GlobalRow. If row
    /// <tt>GlobalRow</tt> is not owned by the calling process, then
    /// \c Indices and \c Values are unchanged and \c NumIndices is
    /// returned as Teuchos::OrdinalTraits<size_t>::invalid().
    void
    getGlobalRowCopy (GlobalOrdinal GlobalRow,
                      const Teuchos::ArrayView<GlobalOrdinal>& Indices,
                      const Teuchos::ArrayView<Scalar>& Values,
                      size_t& NumEntries) const;

    /// \brief Fill given arrays with a deep copy of the locally owned
    ///   entries of the matrix in a given row, using local column
    ///   indices.
    ///
    /// \param localRow [in] Local index of the row for which to
    ///   return entries.
    /// \param colInds [out] Local column indices corresponding to
    ///   values.
    /// \param vals [out] Matrix values.
    /// \param numEntries [out] Number of entries returned.
    ///
    /// Note: A std::runtime_error exception is thrown if either
    /// <tt>colInds</tt> or \c vals is not large enough to hold the
    /// data associated with row \c localRow. If row \c localRow is
    /// not owned by the calling process, then <tt>colInds</tt> and
    /// <tt>vals</tt> are unchanged and <tt>numEntries</tt> is
    /// returned as Teuchos::OrdinalTraits<size_t>::invalid().
    void
    getLocalRowCopy (LocalOrdinal localRow,
                     const Teuchos::ArrayView<LocalOrdinal>& colInds,
                     const Teuchos::ArrayView<Scalar>& vals,
                     size_t& numEntries) const;

    /// \brief Get a constant, nonpersisting view of a row of this
    ///   matrix, using global row and column indices.
    ///
    /// \param GlobalRow [in] Global index of the row to view.
    /// \param indices [out] On output: view of the global column
    ///   indices in the row.
    /// \param values [out] On output: view of the values in the row.
    ///
    /// \pre <tt>isLocallyIndexed () == false</tt>
    /// \post <tt>indices.size () == this->getNumEntriesInGlobalRow (GlobalRow)</tt>
    ///
    /// If \c GlobalRow is not a valid global row index on the calling
    /// process, then \c indices is set to null.
    void
    getGlobalRowView (GlobalOrdinal GlobalRow,
                      Teuchos::ArrayView<const GlobalOrdinal>& indices,
                      Teuchos::ArrayView<const Scalar>& values) const;

    /// \brief Get a constant, nonpersisting view of a row of this
    ///   matrix, using local row and column indices.
    ///
    /// \param LocalRow [in] Local index of the row to view.
    /// \param indices [out] On output: view of the local column
    ///   indices in the row.
    /// \param values [out] On output: view of the values in the row.
    ///
    /// \pre <tt>isGloballyIndexed () == false</tt>
    /// \post <tt>indices.size () == this->getNumEntriesInLocalRow (LocalRow)</tt>
    ///
    /// If \c LocalRow is not a valid local row index on the calling
    /// process, then \c indices is set to null.
    void
    getLocalRowView (LocalOrdinal LocalRow,
                     Teuchos::ArrayView<const LocalOrdinal>& indices,
                     Teuchos::ArrayView<const Scalar>& values) const;

    /// \brief Get a copy of the diagonal entries of the matrix.
    ///
    /// This method returns a Vector with the same Map as this
    /// matrix's row Map.  On each process, it contains the diagonal
    /// entries owned by the calling process.
    void getLocalDiagCopy (Vector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& diag) const;

    /// \brief Get offsets of the diagonal entries in the matrix.
    ///
    /// \warning This method is only for expert users.
    /// \warning We make no promises about backwards compatibility
    ///   for this method.  It may disappear or change at any time.
    /// \warning This method must be called collectively.  We reserve
    ///   the right to do extra checking in a debug build that will
    ///   require collectives.
    ///
    /// \pre The matrix must be locally indexed (which means that it
    ///   has a column Map).
    /// \pre All diagonal entries of the matrix's graph must be
    ///   populated on this process.  Results are undefined otherwise.
    /// \post <tt>offsets.size() == getNodeNumRows()</tt>
    ///
    /// This method creates an array of offsets of the local diagonal
    /// entries in the matrix.  This array is suitable for use in the
    /// two-argument version of getLocalDiagCopy().  However, its
    /// contents are not defined in any other context.  For example,
    /// you should not rely on offsets[i] being the index of the
    /// diagonal entry in the views returned by getLocalRowView().
    /// This may be the case, but it need not be.  (For example, we
    /// may choose to optimize the lookups down to the optimized
    /// storage level, in which case the offsets will be computed with
    /// respect to the underlying storage format, rather than with
    /// respect to the views.)
    ///
    /// Calling any of the following invalidates the output array:
    /// - insertGlobalValues()
    /// - insertLocalValues()
    /// - fillComplete() (with a dynamic graph)
    /// - resumeFill() (with a dynamic graph)
    ///
    /// If the matrix has a const ("static") graph, and if that graph
    /// is fill complete, then the offsets array remains valid through
    /// calls to fillComplete() and resumeFill().  "Invalidates" means
    /// that you must call this method again to recompute the offsets.
    void getLocalDiagOffsets (Teuchos::ArrayRCP<size_t>& offsets) const;

    /// \brief Variant of getLocalDiagCopy() that uses precomputed offsets.
    ///
    /// \warning This method is only for expert users.
    /// \warning We make no promises about backwards compatibility
    ///   for this method.  It may disappear or change at any time.
    ///
    /// This method uses the offsets of the diagonal entries, as
    /// precomputed by getLocalDiagOffsets(), to speed up copying the
    /// diagonal of the matrix.  The offsets must be recomputed if any
    /// of the following methods are called:
    /// - insertGlobalValues()
    /// - insertLocalValues()
    /// - fillComplete() (with a dynamic graph)
    /// - resumeFill() (with a dynamic graph)
    ///
    /// If the matrix has a const ("static") graph, and if that graph
    /// is fill complete, then the offsets array remains valid through
    /// calls to fillComplete() and resumeFill().
    void
    getLocalDiagCopy (Vector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& diag,
                      const Teuchos::ArrayView<const size_t>& offsets) const;

    /** \brief . */
    void
    leftScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& x);

    /** \brief . */
    void
    rightScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& x);

    //@}
    //! @name Advanced templated methods
    //@{

    /// \brief Compute a sparse matrix-MultiVector product local to each process.
    ///
    /// This method computes the <i>local</i> part of <tt>Y := beta*Y
    /// + alpha*Op(A)*X</tt>, where <tt>Op(A)</tt> is either \f$A\f$,
    /// \f$A^T\f$ (the transpose), or \f$A^H\f$ (the conjugate
    /// transpose).  "The local part" means that this method does no
    /// communication between processes, even if this is necessary for
    /// correctness of the matrix-vector multiply.  Use the apply()
    /// method if you want to compute the mathematical sparse
    /// matrix-vector multiply.
    ///
    /// This method is mainly of use to Tpetra developers, though some
    /// users may find it helpful if they plan to reuse the result of
    /// doing an Import on the input MultiVector for several sparse
    /// matrix-vector multiplies with matrices that have the same
    /// column Map.
    ///
    /// When <tt>Op(A)</tt> is \f$A\f$ (<tt>trans ==
    /// Teuchos::NO_TRANS</tt>), then X's Map must be the same as the
    /// column Map of this matrix, and Y's Map must be the same as the
    /// row Map of this matrix.  We say in this case that X is
    /// "post-Imported," and Y is "pre-Exported."  When <tt>Op(A)</tt>
    /// is \f$A^T\f$ or \f$A^H\f$ (\c trans is <tt>Teuchos::TRANS</tt>
    /// or <tt>Teuchos::CONJ_TRANS</tt>, then X's Map must be the same
    /// as the row Map of this matrix, and Y's Map must be the same as
    /// the column Map of this matrix.
    ///
    /// Both X and Y must have constant stride, and they may not alias
    /// one another (that is, occupy overlapping space in memory).  We
    /// may not necessarily check for aliasing, and if we do, we will
    /// only do this in a debug build.  Aliasing X and Y may cause
    /// nondeterministically incorrect results.
    ///
    /// This method is templated on the type of entries in both the
    /// input MultiVector (\c DomainScalar) and the output MultiVector
    /// (\c RangeScalar).  Thus, this method works for MultiVector
    /// objects of arbitrary type.  However, this method only performs
    /// computation local to each MPI process.  Use
    /// CrsMatrixMultiplyOp to handle global communication (the Import
    /// and Export operations for the input resp. output MultiVector),
    /// if you have a matrix with entries of a different type than the
    /// input and output MultiVector objects.
    ///
    /// If <tt>beta == 0</tt>, this operation will enjoy overwrite
    /// semantics: Y will be overwritten with the result of the
    /// multiplication, even if it contains <tt>NaN</tt>
    /// (not-a-number) floating-point entries.  Otherwise, the
    /// multiply result will be accumulated into \c Y.
    template <class DomainScalar, class RangeScalar>
    void
    localMultiply (const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,node_type>& X,
                   MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,node_type>& Y,
                   Teuchos::ETransp trans,
                   RangeScalar alpha,
                   RangeScalar beta) const;

    /// \brief Gauss-Seidel or SOR on \f$B = A X\f$.
    ///
    /// Apply a forward or backward sweep of Gauss-Seidel or
    /// Successive Over-Relaxation (SOR) to the linear system(s) \f$B
    /// = A X\f$.  For Gauss-Seidel, set the damping factor \c omega
    /// to 1.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in B.
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector B.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param B [in] Right-hand side(s).
    /// \param X [in/out] On input: initial guess(es).  On output:
    ///   result multivector(s).
    /// \param D [in] Inverse of diagonal entries of the matrix A.
    /// \param omega [in] SOR damping factor.  omega = 1 results in
    ///   Gauss-Seidel.
    /// \param direction [in] Sweep direction: KokkosClassic::Forward or
    ///   KokkosClassic::Backward.  ("Symmetric" requires interprocess
    ///   communication (before each sweep), which is not part of the
    ///   local kernel.)
    template <class DomainScalar, class RangeScalar>
    void
    localGaussSeidel (const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,node_type> &B,
                      MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,node_type> &X,
                      const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type> &D,
                      const RangeScalar& dampingFactor,
                      const KokkosClassic::ESweepDirection direction) const;

    /// \brief Reordered Gauss-Seidel or SOR on \f$B = A X\f$.
    ///
    /// Apply a forward or backward sweep of reordered Gauss-Seidel or
    /// Successive Over-Relaxation (SOR) to the linear system(s) \f$B
    /// = A X\f$.  For Gauss-Seidel, set the damping factor \c omega
    /// to 1.  The ordering can be a partial one, in which case the Gauss-Seidel is only
    /// executed on a local subset of unknowns.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in B.
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector B.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param B [in] Right-hand side(s).
    /// \param X [in/out] On input: initial guess(es).  On output:
    ///   result multivector(s).
    /// \param D [in] Inverse of diagonal entries of the matrix A.
    /// \param rowIndices [in] Ordered list of indices on which to execute GS.
    /// \param omega [in] SOR damping factor.  omega = 1 results in
    ///   Gauss-Seidel.
    /// \param direction [in] Sweep direction: KokkosClassic::Forward or
    ///   KokkosClassic::Backward.  ("Symmetric" requires interprocess
    ///   communication (before each sweep), which is not part of the
    ///   local kernel.)
    template <class DomainScalar, class RangeScalar>
    void
    reorderedLocalGaussSeidel (const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,node_type>& B,
                               MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,node_type>& X,
                               const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& D,
                               const Teuchos::ArrayView<LocalOrdinal>& rowIndices,
                               const RangeScalar& dampingFactor,
                               const KokkosClassic::ESweepDirection direction) const;

    /// \brief Solves a linear system when the underlying matrix is triangular.
    ///
    /// X is required to be post-imported, i.e., described by the
    /// column map of the matrix. Y is required to be pre-exported,
    /// i.e., described by the row map of the matrix.
    ///
    /// This method is templated on the scalar type of MultiVector
    /// objects, allowing this method to be applied to MultiVector
    /// objects of arbitrary type. However, it is recommended that
    /// solve() not be called directly; instead, use the
    /// CrsMatrixSolveOp, as it will handle the Import/Export operations
    /// required to apply a matrix with non-trivial communication needs.
    ///
    /// Both X and Y are required to have constant stride. However,
    /// unlike multiply(), it is permissible for <tt>&X == &Y</tt>. No
    /// run-time checking will be performed in a non-debug build.
    template <class DomainScalar, class RangeScalar>
    void
    localSolve (const MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,node_type>& Y,
                MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,node_type>& X,
                Teuchos::ETransp trans) const;

    /// \return Return another CrsMatrix with the same entries, but
    ///   converted into a different Scalar type \c T.
    template <class T>
    Teuchos::RCP<CrsMatrix<T, LocalOrdinal, GlobalOrdinal, node_type> >
    convert () const;

    //@}
    //! @name Methods implementing Operator
    //@{

    /// \brief Compute a sparse matrix-MultiVector multiply.
    ///
    /// This method computes <tt>Y := beta*Y + alpha*Op(A)*X</tt>,
    /// where <tt>Op(A)</tt> is either \f$A\f$, \f$A^T\f$ (the
    /// transpose), or \f$A^H\f$ (the conjugate transpose).
    ///
    /// If <tt>beta == 0</tt>, this operation will enjoy overwrite
    /// semantics: Y's entries will be ignored, and Y will be
    /// overwritten with the result of the multiplication, even if it
    /// contains <tt>NaN</tt> (not-a-number) floating-point entries.
    void
    apply (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& X,
           MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>&Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           Scalar alpha = ScalarTraits<Scalar>::one(),
           Scalar beta = ScalarTraits<Scalar>::zero()) const;

    //! Whether apply() allows applying the transpose or conjugate transpose.
    bool hasTransposeApply () const;

    /// \brief The domain Map of this matrix.
    ///
    /// This method implements Tpetra::Operator.  If fillComplete()
    /// has not yet been called at least once on this matrix, or if
    /// the matrix was not constructed with a domain Map, then this
    /// method returns Teuchos::null.
    Teuchos::RCP<const map_type> getDomainMap () const;

    /// \brief The range Map of this matrix.
    ///
    /// This method implements Tpetra::Operator.  If fillComplete()
    /// has not yet been called at least once on this matrix, or if
    /// the matrix was not constructed with a domain Map, then this
    /// method returns Teuchos::null.
    Teuchos::RCP<const map_type> getRangeMap () const;

    //@}
    //! @name Other "apply"-like methods
    //@{

    /// \brief "Hybrid" Jacobi + (Gauss-Seidel or SOR) on \f$B = A X\f$.
    ///
    /// "Hybrid" means Successive Over-Relaxation (SOR) or
    /// Gauss-Seidel within an (MPI) process, but Jacobi between
    /// processes.  Gauss-Seidel is a special case of SOR, where the
    /// damping factor is one.
    ///
    /// The Forward or Backward sweep directions have their usual SOR
    /// meaning within the process.  Interprocess communication occurs
    /// once before the sweep, as it normally would in Jacobi.
    ///
    /// The Symmetric sweep option means two sweeps: first Forward,
    /// then Backward.  Interprocess communication occurs before each
    /// sweep, as in Jacobi.  Thus, Symmetric results in two
    /// interprocess communication steps.
    ///
    /// \param B [in] Right-hand side(s).
    /// \param X [in/out] On input: initial guess(es).  On output:
    ///   result multivector(s).
    /// \param D [in] Inverse of diagonal entries of the matrix A.
    /// \param dampingFactor [in] SOR damping factor.  A damping
    ///   factor of one results in Gauss-Seidel.
    /// \param direction [in] Sweep direction: Forward, Backward, or
    ///   Symmetric.
    /// \param numSweeps [in] Number of sweeps.  We count each
    ///   Symmetric sweep (including both its Forward and its Backward
    ///   sweep) as one.
    ///
    /// \section Tpetra_KR_CrsMatrix_gaussSeidel_req Requirements
    ///
    /// This method has the following requirements:
    ///
    /// 1. X is in the domain Map of the matrix.
    /// 2. The domain and row Maps of the matrix are the same.
    /// 3. The column Map contains the domain Map, and both start at the same place.
    /// 4. The row Map is uniquely owned.
    /// 5. D is in the row Map of the matrix.
    /// 6. X is actually a view of a column Map multivector.
    /// 7. Neither B nor D alias X.
    ///
    /// #1 is just the usual requirement for operators: the input
    /// multivector must always be in the domain Map.  The
    /// Gauss-Seidel kernel imposes additional requirements, since it
    ///
    /// - overwrites the input multivector with the output (which
    ///   implies #2), and
    /// - uses the same local indices for the input and output
    ///   multivector (which implies #2 and #3).
    ///
    /// #3 is reasonable if the matrix constructed the column Map,
    /// because the method that does this (CrsGraph::makeColMap) puts
    /// the local GIDs (those in the domain Map) in front and the
    /// remote GIDs (not in the domain Map) at the end of the column
    /// Map.  However, if you constructed the column Map yourself, you
    /// are responsible for maintaining this invariant.  #6 lets us do
    /// the Import from the domain Map to the column Map in place.
    ///
    /// The Gauss-Seidel kernel also assumes that each process has the
    /// entire value (not a partial value to sum) of all the diagonal
    /// elements in the rows in its row Map.  (We guarantee this anyway
    /// though the separate D vector.)  This is because each element of
    /// the output multivector depends nonlinearly on the diagonal
    /// elements.  Shared ownership of off-diagonal elements would
    /// produce different results.
    void
    gaussSeidel (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type> &B,
                 MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type> &X,
                 const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type> &D,
                 const Scalar& dampingFactor,
                 const ESweepDirection direction,
                 const int numSweeps) const;

    /// \brief Reordered "Hybrid" Jacobi + (Gauss-Seidel or SOR) on \f$B = A X\f$.
    ///
    /// "Hybrid" means Successive Over-Relaxation (SOR) or
    /// Gauss-Seidel within an (MPI) process, but Jacobi between
    /// processes.  Gauss-Seidel is a special case of SOR, where the
    /// damping factor is one.  The ordering can be a partial one, in which case the Gauss-Seidel is only
    /// executed on a local subset of unknowns.
    ///
    /// The Forward or Backward sweep directions have their usual SOR
    /// meaning within the process.  Interprocess communication occurs
    /// once before the sweep, as it normally would in Jacobi.
    ///
    /// The Symmetric sweep option means two sweeps: first Forward,
    /// then Backward.  Interprocess communication occurs before each
    /// sweep, as in Jacobi.  Thus, Symmetric results in two
    /// interprocess communication steps.
    ///
    /// \param B [in] Right-hand side(s).
    /// \param X [in/out] On input: initial guess(es).  On output:
    ///   result multivector(s).
    /// \param D [in] Inverse of diagonal entries of the matrix A.
    /// \param rowIndices [in] Ordered list of indices on which to execute GS.
    /// \param dampingFactor [in] SOR damping factor.  A damping
    ///   factor of one results in Gauss-Seidel.
    /// \param direction [in] Sweep direction: Forward, Backward, or
    ///   Symmetric.
    /// \param numSweeps [in] Number of sweeps.  We count each
    ///   Symmetric sweep (including both its Forward and its Backward
    ///   sweep) as one.
    ///
    /// \section Tpetra_KR_CrsMatrix_reorderedGaussSeidel_req Requirements
    ///
    /// This method has the following requirements:
    ///
    /// 1. X is in the domain Map of the matrix.
    /// 2. The domain and row Maps of the matrix are the same.
    /// 3. The column Map contains the domain Map, and both start at the same place.
    /// 4. The row Map is uniquely owned.
    /// 5. D is in the row Map of the matrix.
    /// 6. X is actually a view of a column Map multivector.
    /// 7. Neither B nor D alias X.
    ///
    /// #1 is just the usual requirement for operators: the input
    /// multivector must always be in the domain Map.  The
    /// Gauss-Seidel kernel imposes additional requirements, since it
    ///
    /// - overwrites the input multivector with the output (which
    ///   implies #2), and
    /// - uses the same local indices for the input and output
    ///   multivector (which implies #2 and #3).
    ///
    /// #3 is reasonable if the matrix constructed the column Map,
    /// because the method that does this (CrsGraph::makeColMap) puts
    /// the local GIDs (those in the domain Map) in front and the
    /// remote GIDs (not in the domain Map) at the end of the column
    /// Map.  However, if you constructed the column Map yourself, you
    /// are responsible for maintaining this invariant.  #6 lets us do
    /// the Import from the domain Map to the column Map in place.
    ///
    /// The Gauss-Seidel kernel also assumes that each process has the
    /// entire value (not a partial value to sum) of all the diagonal
    /// elements in the rows in its row Map.  (We guarantee this anyway
    /// though the separate D vector.)  This is because each element of
    /// the output multivector depends nonlinearly on the diagonal
    /// elements.  Shared ownership of off-diagonal elements would
    /// produce different results.
    void
    reorderedGaussSeidel (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& B,
                          MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& X,
                          const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& D,
                          const Teuchos::ArrayView<LocalOrdinal>& rowIndices,
                          const Scalar& dampingFactor,
                          const ESweepDirection direction,
                          const int numSweeps) const;

    /// \brief Version of gaussSeidel(), with fewer requirements on X.
    ///
    /// This method is just like gaussSeidel(), except that X need
    /// only be in the domain Map.  This method does not require that
    /// X be a domain Map view of a column Map multivector.  As a
    /// result, this method must copy X into a domain Map multivector
    /// before operating on it.
    ///
    /// \param X [in/out] On input: initial guess(es).  On output:
    ///   result multivector(s).
    /// \param B [in] Right-hand side(s), in the range Map.
    /// \param D [in] Inverse of diagonal entries of the matrix,
    ///   in the row Map.
    /// \param dampingFactor [in] SOR damping factor.  A damping
    ///   factor of one results in Gauss-Seidel.
    /// \param direction [in] Sweep direction: Forward, Backward, or
    ///   Symmetric.
    /// \param numSweeps [in] Number of sweeps.  We count each
    ///   Symmetric sweep (including both its Forward and its
    ///   Backward sweep) as one.
    /// \param zeroInitialGuess [in] If true, this method will fill X
    ///   with zeros initially.  If false, this method will assume
    ///   that X contains a possibly nonzero initial guess on input.
    ///   Note that a nonzero initial guess may impose an additional
    ///   nontrivial communication cost (an additional Import).
    ///
    /// \pre Domain, range, and row Maps of the sparse matrix are all the same.
    /// \pre No other argument aliases X.
    void
    gaussSeidelCopy (MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type> &X,
                     const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type> &B,
                     const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type> &D,
                     const Scalar& dampingFactor,
                     const ESweepDirection direction,
                     const int numSweeps,
                     const bool zeroInitialGuess) const;

    /// \brief Version of reorderedGaussSeidel(), with fewer requirements on X.
    ///
    /// This method is just like reorderedGaussSeidel(), except that X need
    /// only be in the domain Map.  This method does not require that
    /// X be a domain Map view of a column Map multivector.  As a
    /// result, this method must copy X into a domain Map multivector
    /// before operating on it.
    ///
    /// \param X [in/out] On input: initial guess(es).  On output:
    ///   result multivector(s).
    /// \param B [in] Right-hand side(s), in the range Map.
    /// \param D [in] Inverse of diagonal entries of the matrix,
    ///   in the row Map.
    /// \param rowIndices [in] Ordered list of indices on which to execute GS.
    /// \param dampingFactor [in] SOR damping factor.  A damping
    ///   factor of one results in Gauss-Seidel.
    /// \param direction [in] Sweep direction: Forward, Backward, or
    ///   Symmetric.
    /// \param numSweeps [in] Number of sweeps.  We count each
    ///   Symmetric sweep (including both its Forward and its
    ///   Backward sweep) as one.
    /// \param zeroInitialGuess [in] If true, this method will fill X
    ///   with zeros initially.  If false, this method will assume
    ///   that X contains a possibly nonzero initial guess on input.
    ///   Note that a nonzero initial guess may impose an additional
    ///   nontrivial communication cost (an additional Import).
    ///
    /// \pre Domain, range, and row Maps of the sparse matrix are all the same.
    /// \pre No other argument aliases X.
    void
    reorderedGaussSeidelCopy (MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& X,
                              const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& B,
                              const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>& D,
                              const Teuchos::ArrayView<LocalOrdinal>& rowIndices,
                              const Scalar& dampingFactor,
                              const ESweepDirection direction,
                              const int numSweeps,
                              const bool zeroInitialGuess) const;

    /// \brief Implementation of RowMatrix::add: return <tt>alpha*A + beta*this</tt>.
    ///
    /// This override of the default implementation ensures that, when
    /// called on a CrsMatrix, this method always returns a CrsMatrix
    /// of exactly the same type as <tt>*this</tt>.  "Exactly the same
    /// type" means that all the template parameters match, including
    /// the fifth template parameter.  The input matrix A need not
    /// necessarily be a CrsMatrix or a CrsMatrix of the same type as
    /// <tt>*this</tt>, though this method may be able to optimize
    /// further in that case.
    virtual Teuchos::RCP<RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, node_type> >
    add (const Scalar& alpha,
         const RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& A,
         const Scalar& beta,
         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, node_type> >& domainMap,
         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, node_type> >& rangeMap,
         const Teuchos::RCP<Teuchos::ParameterList>& params) const;

    //@}
    //! @name Implementation of Teuchos::Describable interface
    //@{

    //! A one-line description of this object.
    std::string description () const;

    //! Print the object with some verbosity level to an FancyOStream object.
    void
    describe (Teuchos::FancyOStream &out,
              const Teuchos::EVerbosityLevel verbLevel =
              Teuchos::Describable::verbLevel_default) const;

    //@}
    //! @name Implementation of DistObject interface
    //@{

    virtual bool
    checkSizes (const SrcDistObject& source);

    virtual void
    copyAndPermute (const SrcDistObject& source,
                    size_t numSameIDs,
                    const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs,
                    const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs);

    virtual void
    packAndPrepare (const SrcDistObject& source,
                    const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
                    Teuchos::Array<char>& exports,
                    const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                    size_t& constantNumPackets,
                    Distributor& distor);

  private:
    /// \brief Unpack the imported column indices and values, and
    ///   combine into matrix.
    void
    unpackAndCombineImpl (const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                          const Teuchos::ArrayView<const char> &imports,
                          const Teuchos::ArrayView<size_t> &numPacketsPerLID,
                          size_t constantNumPackets,
                          Distributor& distor,
                          CombineMode combineMode);

  public:
    /// \brief Unpack the imported column indices and values, and combine into matrix.
    ///
    /// \warning The allowed \c combineMode depends on whether the
    ///   matrix's graph is static or dynamic.  ADD, REPLACE, and
    ///   ABSMAX are valid for a static graph, but INSERT is not.
    ///   ADD and INSERT are valid for a dynamic graph; ABSMAX and
    ///   REPLACE have not yet been implemented (and would require
    ///   serious changes to matrix assembly in order to implement
    ///   sensibly).
    void
    unpackAndCombine (const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                      const Teuchos::ArrayView<const char> &imports,
                      const Teuchos::ArrayView<size_t> &numPacketsPerLID,
                      size_t constantNumPackets,
                      Distributor& distor,
                      CombineMode combineMode);
    //@}
    //! @name Implementation of Packable interface
    //@{

    /// \brief Pack this object's data for an Import or Export.
    ///
    /// \warning To be called only by the packAndPrepare method of
    ///   appropriate classes of DistObject.
    ///
    /// \param exportLIDs [in] Local indices of the rows to pack.
    /// \param exports [out] On output: array of packed matrix
    ///   entries; allocated by method.
    /// \param numPacketsPerLID [out] On output: numPacketsPerLID[i]
    ///   is the number of bytes of the \c exports array used for
    ///   storing packed local row \c exportLIDs[i].
    /// \param constantNumPackets [out] If zero on output, the packed
    ///   rows may have different numbers of entries.  If nonzero on
    ///   output, then that number gives the constant number of
    ///   entries for all packed rows <i>on all processes in the
    ///   matrix's communicator</i>.
    /// \param distor [in/out] The Distributor object which implements
    ///   the Import or Export operation that is calling this method.
    ///
    /// \subsection Tpetra_KR_CrsMatrix_pack_summary Packing scheme
    ///
    /// The number of "packets" per row is the number of bytes per
    /// row.  Each row has the following storage format:
    ///
    /// <tt>[numEnt, vals, inds]</tt>,
    ///
    /// where:
    /// <ul>
    /// <li> \c numEnt (\c LocalOrdinal): number of entries in the
    ///      row. </li>
    /// <li> \c vals: array of \c Scalar.  For the k-th entry in the
    ///      row, \c vals[k] is its value and \c inds[k] its global
    ///      column index. </li>
    /// <li> \c inds: array of \c GlobalOrdinal.  For the k-th entry
    ///      in the row, \c vals[k] is its value and \c inds[k] its
    ///      global column index. </li>
    /// </ul>
    ///
    /// We reserve the right to pad for alignment in the future.  In
    /// that case, the number of bytes reported by \c numPacketsPerLID
    /// will reflect padding to align each datum to its size, and the
    /// row will have final padding as well to ensure that the
    /// <i>next</i> row is aligned.  Rows with zero entries will still
    /// take zero bytes, however.
    ///
    /// RowMatrix::pack will always use the same packing scheme as
    /// this method.  This ensures correct Import / Export from a
    /// RowMatrix to a CrsMatrix.
    ///
    /// We do <i>not</i> recommend relying on the details of this
    /// packing scheme.  We describe it here more for Tpetra
    /// developers and less for users.
    ///
    /// \subsection Tpetra_KR_CrsMatrix_pack_disc Discussion
    ///
    /// DistObject requires packing an object's entries as type
    /// <tt>Packet</tt>, which is the first template parameter of
    /// DistObject.  Since sparse matrices have both values and
    /// indices, we use <tt>Packet=char</tt> and pack them into
    /// buffers of <tt>char</tt> (really "byte").  Indices are stored
    /// as global indices, in case the source and target matrices have
    /// different column Maps (or don't have a column Map yet).
    ///
    /// Currently, we only pack values and column indices.  Row
    /// indices are stored implicitly as the local indices (LIDs) to
    /// pack (see \c exportLIDs).  This is because a DistObject
    /// instance only has one Map, and currently we use the row Map
    /// for CrsMatrix (and RowMatrix).  This makes redistribution of
    /// matrices with 2-D distributions less efficient, but it works
    /// for now.  This may change in the future.
    ///
    /// On output, \c numPacketsPerLID[i] gives the number of bytes
    /// used to pack local row \c exportLIDs[i] of \c this object (the
    /// source object of an Import or Export).  If \c offset is the
    /// exclusive prefix sum-scan of \c numPacketsPerLID, then on
    /// output, <tt>exports[offset[i] .. offset[i+1]]</tt>
    /// (half-exclusive range) contains the packed entries for local
    /// row \c exportLIDs[i].
    ///
    /// Entries for each row use a "struct of arrays" pattern to match
    /// how sparse matrices actually store their data.  The number of
    /// entries in the row goes first, all values go next, and all
    /// column indices (stored as global indices) go last.  Values and
    /// column indices occur in the same order.  Rows with zero
    /// entries always take zero bytes (we do not store their number
    /// of entries explicitly).  This ensures sparsity of storage and
    /// communication in case most rows are empty.
    ///
    /// \subsection Tpetra_KR_CrsMatrix_pack_why Justification
    ///
    /// GCC >= 4.9 and recent-future versions of the Intel compiler
    /// implement stricter aliasing rules that forbid unaligned type
    /// punning.  If we were to pack as an "array of structs" -- in
    /// this case, an array of <tt>(Scalar, GlobalOrdinal)</tt> pairs
    /// -- then we would either have to pad each matrix entry for
    /// alignment, or call memcpy twice per matrix entry to pack and
    /// unpack.  The "struct of arrays" storage scheme reduces the
    /// padding requirement to a constant per row, or reduces the
    /// number of memcpy calls to two per row.
    ///
    /// We include the number of entries in each row in that row's
    /// packed data, to make unpacking easier.  This saves us from an
    /// error-prone computation to find the number of entries from the
    /// number of bytes.  That computation gets even more difficult if
    /// we have to introduce padding for alignment in the future.
    /// Knowing the number of entries for each row also makes
    /// parallelizing packing and unpacking easier.
    ///
    /// \subsection Tpetra_KR_CrsMatrix_pack_assum Technical assumptions
    ///
    /// <ul>
    /// <li> \c sizeof(Scalar) says how much data were used to
    ///      represent a \c Scalar in its packed form. </li>
    /// <li> \c sizeof returns the same value on all processes for
    ///      <tt>Scalar</tt>, \c LocalOrdinal, and \c GlobalOrdinal.
    ///      </li>
    /// </ul>
    virtual void
    pack (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
          Teuchos::Array<char>& exports,
          const Teuchos::ArrayView<size_t>& numPacketsPerLID,
          size_t& constantNumPackets,
          Distributor& distor) const;

  private:
    /// \brief Pack data for the current row to send.
    ///
    /// \param numEntOut [out] Where to write the number of entries in
    ///   the row.
    /// \param valOut [out] Output (packed) array of matrix values.
    /// \param indOut [out] Output (packed) array of matrix column
    ///   indices (as global indices).
    /// \param numEnt [in] Number of entries in the row.
    /// \param lclRow [in] Local index of the row.
    ///
    /// This method does not allocate temporary storage.  We intend
    /// for this to be safe to call in a thread-parallel way at some
    /// point, though it is currently not, due to thread safety issues
    /// with Teuchos::RCP (always) and Teuchos::ArrayView (in a debug
    /// build).
    ///
    /// \return \c true if the method succeeded, else \c false.
    bool
    packRow (char* const numEntOut,
             char* const valOut,
             char* const indOut,
             const size_t numEnt,
             const LocalOrdinal lclRow) const;

    /// \brief Unpack and combine received data for the current row.
    ///
    /// \pre <tt>tmpSize >= numEnt</tt>
    ///
    /// \param valInTmp [out] Temporary storage for values.  Has
    ///   tmpSize entries.
    /// \param indInTmp [out] Temporary storage for indices.  Has
    ///   tmpSize entries.
    /// \param tmpNumEnt [in] Number of entries (not bytes!) in each
    ///   of valInTmp and indInTmp.
    /// \param valIn [in] Pointer to where values live in receive
    ///   buffer.  Not necessarily aligned to sizeof(Scalar) (so must
    ///   memcpy into temporary storage).
    /// \param indIn [out] Pointer to where indices live in receive
    ///   buffer.  Not necessarily aligned to sizeof(GlobalOrdinal)
    ///   (so must memcpy into temporary storage).
    /// \param numEnt [in] Number of entries in the row.
    /// \param lclRow [in] Local index of the row.
    /// \param combineMode [in] Combine mode (how to merge entries in
    ///   the same row with the same column index).
    ///
    /// \return \c true if the method succeeded, else \c false.
    bool
    unpackRow (Scalar* const valInTmp,
               GlobalOrdinal* const indInTmp,
               const size_t tmpNumEnt,
               const char* const valIn,
               const char* const indIn,
               const size_t numEnt,
               const LocalOrdinal lclRow,
               const Tpetra::CombineMode combineMode);

    /// \brief Allocate space for pack() to pack entries to send.
    ///
    /// \param exports [in/out] Pack buffer to (re)allocate.
    /// \param totalNumEntries [out] Total number of entries to send.
    /// \param exportLIDs [in] Local indices of the rows to send.
    void
    allocatePackSpace (Teuchos::Array<char>& exports,
                       size_t& totalNumEntries,
                       const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs) const;
    //@}

  public:
    //! Get the Kokkos local values
    t_ValuesType getLocalValuesView () const { return k_values1D_; }

  private:
    // Friend declaration for nonmember function.
    template<class CrsMatrixType>
    friend Teuchos::RCP<CrsMatrixType>
    importAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrixType>& sourceMatrix,
                                    const Import<typename CrsMatrixType::local_ordinal_type,
                                                 typename CrsMatrixType::global_ordinal_type,
                                                 typename CrsMatrixType::node_type>& importer,
                                    const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                                 typename CrsMatrixType::global_ordinal_type,
                                                                 typename CrsMatrixType::node_type> >& domainMap,
                                    const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                                 typename CrsMatrixType::global_ordinal_type,
                                                                 typename CrsMatrixType::node_type> >& rangeMap,
                                    const Teuchos::RCP<Teuchos::ParameterList>& params);

    // Friend declaration for nonmember function.
    template<class CrsMatrixType>
    friend Teuchos::RCP<CrsMatrixType>
    exportAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrixType>& sourceMatrix,
                                    const Export<typename CrsMatrixType::local_ordinal_type,
                                                 typename CrsMatrixType::global_ordinal_type,
                                                 typename CrsMatrixType::node_type>& exporter,
                                    const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                                 typename CrsMatrixType::global_ordinal_type,
                                                                 typename CrsMatrixType::node_type> >& domainMap,
                                    const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                                 typename CrsMatrixType::global_ordinal_type,
                                                                 typename CrsMatrixType::node_type> >& rangeMap,
                                    const Teuchos::RCP<Teuchos::ParameterList>& params);

  public:
    /// \brief Import from <tt>this</tt> to the given destination
    ///   matrix, and make the result fill complete.
    ///
    /// If destMatrix.is_null(), this creates a new matrix as the
    /// destination.  (This is why destMatrix is passed in by nonconst
    /// reference to RCP.)  Otherwise it checks for "pristine" status
    /// and throws if that is not the case.  "Pristine" means that the
    /// matrix has no entries and is not fill complete.
    ///
    /// Use of the "non-member constructor" version of this method,
    /// exportAndFillCompleteCrsMatrix, is preferred for user
    /// applications.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    void
    importAndFillComplete (Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, node_type> >& destMatrix,
                           const import_type& importer,
                           const Teuchos::RCP<const map_type>& domainMap,
                           const Teuchos::RCP<const map_type>& rangeMap,
                           const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) const;

    /// \brief Export from <tt>this</tt> to the given destination
    ///   matrix, and make the result fill complete.
    ///
    /// If destMatrix.is_null(), this creates a new matrix as the
    /// destination.  (This is why destMatrix is passed in by nonconst
    /// reference to RCP.)  Otherwise it checks for "pristine" status
    /// and throws if that is not the case.  "Pristine" means that the
    /// matrix has no entries and is not fill complete.
    ///
    /// Use of the "non-member constructor" version of this method,
    /// exportAndFillCompleteCrsMatrix, is preferred for user
    /// applications.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    void
    exportAndFillComplete (Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, node_type> >& destMatrix,
                           const export_type& exporter,
                           const Teuchos::RCP<const map_type>& domainMap = Teuchos::null,
                           const Teuchos::RCP<const map_type>& rangeMap = Teuchos::null,
                           const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) const;

  private:
    /// \brief Transfer (e.g. Import/Export) from <tt>this</tt> to the
    ///   given destination matrix, and make the result fill complete.
    ///
    /// If destMat.is_null(), this creates a new matrix, otherwise it
    /// checks for "pristine" status and throws if that is not the
    /// case.  This method implements importAndFillComplete and
    /// exportAndFillComplete, which in turn implemment the nonmember
    /// "constructors" importAndFillCompleteCrsMatrix and
    /// exportAndFillCompleteCrsMatrix.  It's convenient to put those
    /// nonmember constructors' implementations inside the CrsMatrix
    /// class, so that we don't have to put much code in the _decl
    /// header file.
    ///
    /// The point of this method is to fuse three tasks:
    ///
    ///   1. Create a destination matrix (CrsMatrix constructor)
    ///   2. Import or Export this matrix to the destination matrix
    ///   3. Call fillComplete on the destination matrix
    ///
    /// Fusing these tasks can avoid some communication and work.
    void
    transferAndFillComplete (Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, node_type> >& destMatrix,
                             const ::Tpetra::Details::Transfer<LocalOrdinal, GlobalOrdinal, node_type>& rowTransfer,
                             const Teuchos::RCP<const map_type>& domainMap = Teuchos::null,
                             const Teuchos::RCP<const map_type>& rangeMap = Teuchos::null,
                             const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) const;

    // We forbid copy construction by declaring this method private
    // and not implementing it.
    CrsMatrix (const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& rhs);

    // We forbid assignment (operator=) by declaring this method
    // private and not implementing it.
    CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, node_type>&
    operator= (const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& rhs);

    /// \brief Like insertGlobalValues(), but with column filtering.
    ///
    /// "Column filtering" means that if the matrix has a column Map,
    /// then this method ignores entries in columns that are not in
    /// the column Map.
    ///
    /// See discussion in the documentation of getGlobalRowCopy()
    /// about why we use \c Scalar and not \c impl_scalar_type here
    /// for the input array type.
    void
    insertGlobalValuesFiltered (const GlobalOrdinal globalRow,
                                const Teuchos::ArrayView<const GlobalOrdinal>& indices,
                                const Teuchos::ArrayView<const Scalar>& values);

    /// \brief Like insertLocalValues(), but with column filtering.
    ///
    /// "Column filtering" means that if the matrix has a column Map,
    /// then this method ignores entries in columns that are not in
    /// the column Map.
    ///
    /// See discussion in the documentation of getGlobalRowCopy()
    /// about why we use \c Scalar and not \c impl_scalar_type here
    /// for the input array type.
    void
    insertLocalValuesFiltered (const LocalOrdinal localRow,
                               const Teuchos::ArrayView<const LocalOrdinal>& indices,
                               const Teuchos::ArrayView<const Scalar>& values);

    /// \brief Combine in the data using the given combine mode.
    ///
    /// The copyAndPermute() and unpackAndCombine() methods use this
    /// function to combine incoming entries from the source matrix
    /// with the target matrix's current data.  This method's behavior
    /// depends on whether the target matrix (that is, this matrix)
    /// has a static graph.
    ///
    /// See discussion in the documentation of getGlobalRowCopy()
    /// about why we use \c Scalar and not \c impl_scalar_type here
    /// for the input array type.
    void
    combineGlobalValues (const GlobalOrdinal globalRowIndex,
                         const Teuchos::ArrayView<const GlobalOrdinal>& columnIndices,
                         const Teuchos::ArrayView<const Scalar>& values,
                         const Tpetra::CombineMode combineMode);

    /// \brief Transform CrsMatrix entries, using global indices.
    ///
    /// For every entry \f$A(i,j)\f$ to transform, if \f$v_{ij}\f$ is
    /// the corresponding entry of the \c values array, then we apply
    /// the binary function f to \f$A(i,j)\f$ as follows:
    /// \f[
    ///   A(i,j) := f(A(i,j), v_{ij}).
    /// \f]
    /// For example, BinaryFunction = std::plus<Scalar> does the same
    /// thing as sumIntoLocalValues(), and BinaryFunction =
    /// project2nd<Scalar,Scalar> does the same thing as
    /// replaceLocalValues().
    ///
    /// \tparam BinaryFunction The type of binary function to apply.
    ///   std::binary_function is a model for this.
    ///
    /// \param globalRow [in] (Global) index of the row to modify.
    ///   This row <i>must</t> be owned by the calling process.
    ///
    /// \param indices [in] (Global) indices in the row to modify.
    ///   Indices not in the column Map (if the matrix already has a
    ///   column Map) and their corresponding values will be ignored.
    ///
    /// \param values [in] Values to use for modification.
    ///
    /// This method works whether indices are local or global.
    /// However, it will cost more if indices are local, since it will
    /// have to convert the local indices to global indices in that
    /// case.
    ///
    /// See discussion in the documentation of getGlobalRowCopy()
    /// about why we use \c Scalar and not \c impl_scalar_type here
    /// for the input array type.
    template<class BinaryFunction>
    LocalOrdinal
    transformGlobalValues (GlobalOrdinal globalRow,
                           const Teuchos::ArrayView<const GlobalOrdinal>& indices,
                           const Teuchos::ArrayView<const Scalar>& values,
                           BinaryFunction f)
    {
      using Teuchos::Array;
      using Teuchos::ArrayView;
      typedef LocalOrdinal LO;
      typedef GlobalOrdinal GO;
      typedef impl_scalar_type ST;
      typedef typename ArrayView<const GO>::size_type size_type;

      if (! isFillActive ()) {
        // Fill must be active in order to call this method.
        return Teuchos::OrdinalTraits<LO>::invalid ();
      }
      else if (values.size () != indices.size ()) {
        // The sizes of values and indices must match.
        return Teuchos::OrdinalTraits<LO>::invalid ();
      }

      const LO lrow = this->getRowMap ()->getLocalElement (globalRow);
      if (lrow == Teuchos::OrdinalTraits<LO>::invalid ()) {
        // We don't own the row, so we're not allowed to modify its values.
        return Teuchos::OrdinalTraits<LO>::invalid ();
      }

      if (staticGraph_.is_null ()) {
        return Teuchos::OrdinalTraits<LO>::invalid ();
      }
      const crs_graph_type& graph = *staticGraph_;
      RowInfo rowInfo = graph.getRowInfo (lrow);
      if (indices.size () == 0) {
        return static_cast<LO> (0);
      }
      else {
        ArrayView<const ST> valsIn =
          Teuchos::av_reinterpret_cast<const ST> (values);
        ArrayView<ST> curVals = this->getViewNonConst (rowInfo);
        if (isLocallyIndexed ()) {
          // Convert the given global indices to local indices.
          //
          // FIXME (mfh 08 Jul 2014) Why can't we ask the graph to do
          // that?  It could do the conversions in place, so that we
          // wouldn't need temporary storage.
          const map_type& colMap = * (this->getColMap ());
          const size_type numInds = indices.size ();
          Array<LO> lclInds (numInds);
          for (size_type k = 0; k < numInds; ++k) {
            // There is no need to filter out indices not in the
            // column Map.  Those that aren't will be mapped to
            // invalid(), which the graph's transformGlobalValues()
            // will filter out (but not count in its return value).
            lclInds[k] = colMap.getLocalElement (indices[k]);
          }
          return graph.template transformLocalValues<ST, BinaryFunction> (rowInfo,
                                                                          curVals,
                                                                          lclInds (),
                                                                          valsIn, f);
        }
        else if (isGloballyIndexed ()) {
          return graph.template transformGlobalValues<ST, BinaryFunction> (rowInfo,
                                                                           curVals,
                                                                           indices,
                                                                           valsIn, f);
        }
        else {
          // If the graph is neither locally nor globally indexed on
          // the calling process, that means that the calling process
          // can't possibly have any entries in the owned row.  Thus,
          // there are no entries to transform, so we return zero.
          return static_cast<LO> (0);
        }
      }
    }

  private:
    /// \brief Special case of insertGlobalValues for when globalRow
    ///   is <i>not<i> owned by the calling process.
    ///
    /// See discussion in the documentation of getGlobalRowCopy()
    /// about why we use \c Scalar and not \c impl_scalar_type here
    /// for the input array type.
    void
    insertNonownedGlobalValues (const GlobalOrdinal globalRow,
                                const Teuchos::ArrayView<const GlobalOrdinal>& indices,
                                const Teuchos::ArrayView<const Scalar>& values);

    //! Type of the DistObject specialization from which this class inherits.
    typedef DistObject<char, LocalOrdinal, GlobalOrdinal, node_type> dist_object_type;

  protected:
    // useful typedefs
    typedef Teuchos::OrdinalTraits<LocalOrdinal> OTL;
    typedef Kokkos::Details::ArithTraits<impl_scalar_type> STS;
    typedef Kokkos::Details::ArithTraits<mag_type> STM;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type> MV;
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>      V;
    typedef crs_graph_type Graph;

    // Enums
    enum GraphAllocationStatus {
      GraphAlreadyAllocated,
      GraphNotYetAllocated
    };

    /// \brief Allocate values (and optionally indices) using the Node.
    ///
    /// \param gas [in] If GraphNotYetAllocated, allocate the
    ///   indices of \c myGraph_ via \c allocateIndices(lg) before
    ///   allocating values.
    ///
    /// \param lg [in] Argument passed into \c
    ///   myGraph_->allocateIndices(), if applicable.
    ///
    /// \pre If the graph (that is, staticGraph_) indices are
    ///   already allocated, then gas must be GraphAlreadyAllocated.
    ///   Otherwise, gas must be GraphNotYetAllocated.  We only
    ///   check for this precondition in debug mode.
    ///
    /// \pre If the graph indices are not already allocated, then
    ///   the graph must be owned by the matrix.
    void allocateValues (ELocalGlobal lg, GraphAllocationStatus gas);

    /// \brief Sort the entries of each row by their column indices.
    ///
    /// This only does anything if the graph isn't already sorted
    /// (i.e., ! myGraph_->isSorted ()).  This method is called in
    /// fillComplete().
    void sortEntries();

    /// \brief Merge entries in each row with the same column indices.
    ///
    /// This only does anything if the graph isn't already merged
    /// (i.e., ! myGraph_->isMerged ()).  This method is called in
    /// fillComplete().
    void mergeRedundantEntries();

    /// \brief Clear matrix properties that require collectives.
    ///
    /// This clears whatever computeGlobalConstants() (which see)
    /// computed, in preparation for changes to the matrix.  The
    /// current implementation of this method does nothing.
    ///
    /// This method is called in resumeFill().
    void clearGlobalConstants();

    /// \brief Compute matrix properties that require collectives.
    ///
    /// The corresponding Epetra_CrsGraph method computes things
    /// like the global number of nonzero entries, that require
    /// collectives over the matrix's communicator.  The current
    /// Tpetra implementation of this method does nothing.
    ///
    /// This method is called in fillComplete().
    void computeGlobalConstants();

    /// \brief Column Map MultiVector used in apply() and gaussSeidel().
    ///
    /// This is a column Map MultiVector.  It is used as the target of
    /// the forward mode Import operation (if necessary) in apply()
    /// and gaussSeidel(), and the source of the reverse mode Export
    /// operation (if necessary) in these methods.  Both of these
    /// methods create this MultiVector on demand if needed, and reuse
    /// it (if possible) for subsequent calls.
    ///
    /// This is declared <tt>mutable</tt> because the methods in
    /// question are const, yet want to cache the MultiVector for
    /// later use.
    mutable Teuchos::RCP<MV> importMV_;

    /// \brief Row Map MultiVector used in apply().
    ///
    /// This is a row Map MultiVector.  It is uses as the source of
    /// the forward mode Export operation (if necessary) in apply()
    /// and gaussSeidel(), and the target of the reverse mode Import
    /// operation (if necessary) in these methods.  Both of these
    /// methods create this MultiVector on demand if needed, and reuse
    /// it (if possible) for subsequent calls.
    ///
    /// This is declared <tt>mutable</tt> because the methods in
    /// question are const, yet want to cache the MultiVector for
    /// later use.
    mutable Teuchos::RCP<MV> exportMV_;

    /// \brief Create a (or fetch a cached) column Map MultiVector.
    ///
    /// \param X_domainMap [in] A domain Map Multivector.  The
    ///   returned MultiVector, if nonnull, will have the same number
    ///   of columns as Y_domainMap.
    ///
    /// \param force [in] Force creating the MultiVector if it hasn't
    ///   been created already.
    ///
    /// The \c force parameter is helpful when the domain Map and the
    /// column Map are the same (so that normally we wouldn't need the
    /// column Map MultiVector), but the following (for example)
    /// holds:
    ///
    /// 1. The kernel needs a constant stride input MultiVector, but
    ///    the given input MultiVector is not constant stride.
    ///
    /// We don't test for the above in this method, because it depends
    /// on the specific kernel.
    Teuchos::RCP<MV>
    getColumnMapMultiVector (const MV& X_domainMap,
                             const bool force = false) const;

    /// \brief Create a (or fetch a cached) row Map MultiVector.
    ///
    /// \param Y_rangeMap [in] A range Map Multivector.  The returned
    ///   MultiVector, if nonnull, will have the same number of
    ///   columns as Y_rangeMap.
    ///
    /// \param force [in] Force creating the MultiVector if it hasn't
    ///   been created already.
    ///
    /// The \c force parameter is helpful when the range Map and the
    /// row Map are the same (so that normally we wouldn't need the
    /// row Map MultiVector), but one of the following holds:
    ///
    /// 1. The kernel needs a constant stride output MultiVector,
    ///    but the given output MultiVector is not constant stride.
    ///
    /// 2. The kernel does not permit aliasing of its input and output
    ///    MultiVector arguments, but they do alias each other.
    ///
    /// We don't test for the above in this method, because it depends
    /// on the specific kernel.
    Teuchos::RCP<MV>
    getRowMapMultiVector (const MV& Y_rangeMap,
                          const bool force = false) const;

    //! Special case of apply() for <tt>mode == Teuchos::NO_TRANS</tt>.
    void
    applyNonTranspose (const MV& X_in,
                       MV& Y_in,
                       Scalar alpha,
                       Scalar beta) const;

    //! Special case of apply() for <tt>mode != Teuchos::NO_TRANS</tt>.
    void
    applyTranspose (const MV& X_in,
                    MV& Y_in,
                    const Teuchos::ETransp mode,
                    Scalar alpha,
                    Scalar beta) const;

    // matrix data accessors

    /// \brief Constant view of all entries (including extra space) in
    ///   the given row.
    ///
    /// Unlike getGlobalRowView(), this method returns
    /// <tt>impl_scalar_type</tt>, not \c Scalar.  This is because
    /// this method is <i>not</i> part of the public interface of
    /// CrsMatrix.
    Teuchos::ArrayView<const impl_scalar_type> getView (RowInfo rowinfo) const;

    /// \brief Nonconst view of all entries (including extra space) in
    ///   the given row.
    ///
    /// Unlike getGlobalRowView(), this method returns
    /// <tt>impl_scalar_type</tt>, not \c Scalar.  This is because
    /// this method is <i>not</i> part of the public interface of
    /// CrsMatrix.
    Teuchos::ArrayView<impl_scalar_type> getViewNonConst (RowInfo rowinfo);

    /// \brief Fill data into the local matrix.
    ///
    /// This method is only called in fillComplete(), and it is only
    /// called if the graph's structure is already fixed (that is, if
    /// the matrix does not own the graph).
    void fillLocalMatrix (const Teuchos::RCP<Teuchos::ParameterList>& params);

    /// \brief Fill data into the local graph and matrix.
    ///
    /// This method is only called in fillComplete(), and it is only
    /// called if the graph's structure is <i>not</i> already fixed
    /// (that is, if the matrix <i>does</i> own the graph).
    void fillLocalGraphAndMatrix (const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Check that this object's state is sane; throw if it's not.
    void checkInternalState () const;

    /// \name (Global) graph pointers
    ///
    /// We keep two graph pointers in order to maintain const
    /// correctness.  myGraph_ is a graph which we create internally.
    /// Operations that change the sparsity structure also modify
    /// myGraph_.  If myGraph_ != null, then staticGraph_ == myGraph_
    /// pointerwise (we set the pointers equal to each other when we
    /// create myGraph_).  myGraph_ is only null if this CrsMatrix was
    /// created using the constructor with a const CrsGraph input
    /// argument.  In this case, staticGraph_ is set to the input
    /// CrsGraph.
    //@{
    Teuchos::RCP<const Graph> staticGraph_;
    Teuchos::RCP<      Graph>     myGraph_;
    //@}

    //! The local sparse matrix.
    local_matrix_type lclMatrix_;

    /// \name Sparse matrix values.
    ///
    /// k_values1D_ represents the values assuming "1-D" compressed
    /// sparse row storage.  values2D_ represents the values as an
    /// array of arrays, one (inner) array per row of the sparse
    /// matrix.
    ///
    /// Before allocation, both arrays are null.  After allocation,
    /// one is null.  If static allocation, then values2D_ is null.
    /// If dynamic allocation, then k_values1D_ is null.  The
    /// allocation always matches that of graph_, as the graph does
    /// the allocation for the matrix.
    //@{
    t_ValuesType k_values1D_;
    Teuchos::ArrayRCP<Teuchos::Array<impl_scalar_type> > values2D_;
    //@}

    /// \brief Status of the matrix's storage, when not in a
    ///   fill-complete state.
    ///
    /// The phrase "When not in a fill-complete state" is important.
    /// When the matrix is fill complete, it <i>always</i> uses 1-D
    /// "packed" storage.  However, if the "Optimize Storage"
    /// parameter to fillComplete was false, the matrix may keep
    /// unpacked 1-D or 2-D storage around and resume it on the next
    /// resumeFill call.
    Details::EStorageStatus storageStatus_;

    //! Whether the matrix is fill complete.
    bool fillComplete_;

    /// \brief Nonlocal data added using insertGlobalValues().
    ///
    /// These data are cleared by globalAssemble(), once it finishes
    /// redistributing them to their owning processes.
    ///
    /// For a given nonowned global row gRow which was given to
    /// insertGlobalValues() or sumIntoGlobalValues(),
    /// <tt>nonlocals_[gRow].first[k]</tt> is the column index of an
    /// inserted entry, and <tt>nonlocals_[gRow].second[k]</tt> is its
    /// value.  Duplicate column indices for the same row index are
    /// allowed and will be summed during globalAssemble().
    ///
    /// This used to be a map from GlobalOrdinal to (GlobalOrdinal,
    /// Scalar) pairs.  This makes gcc issue a "note" about the ABI of
    /// structs containing std::complex members changing.  CDash
    /// reports this as a warning, even though it's a "note," not a
    /// warning.  However, I don't want it to show up, so I rearranged
    /// the map's value type to a pair of arrays, rather than an array
    /// of pairs.
    ///
    /// \note For Epetra developers: Tpetra::CrsMatrix corresponds
    ///   more to Epetra_FECrsMatrix than to Epetra_CrsMatrix.  The
    ///   insertGlobalValues() method in Tpetra::CrsMatrix, unlike
    ///   its corresponding method in Epetra_CrsMatrix, allows
    ///   insertion into rows which are not owned by the calling
    ///   process.  The globalAssemble() method redistributes these
    ///   to their owning processes.
    std::map<GlobalOrdinal, std::pair<Teuchos::Array<GlobalOrdinal>,
                                      Teuchos::Array<Scalar> > > nonlocals_;

    /// \brief Cached Frobenius norm of the (global) matrix.
    ///
    /// The value -1 means that the norm has not yet been computed, or
    /// that the values in the matrix may have changed and the norm
    /// must be recomputed.
    mutable mag_type frobNorm_;

  public:
    // FIXME (mfh 24 Feb 2014) Is it _really_ necessary to make this a
    // public inner class of CrsMatrix?  It looks like it doesn't
    // depend on any implementation details of CrsMatrix at all.  It
    // should really be declared and defined outside of CrsMatrix.
    template<class ViewType, class OffsetViewType>
    struct pack_functor {
      typedef typename ViewType::execution_space execution_space;
      ViewType src_;
      ViewType dst_;
      OffsetViewType src_offset_;
      OffsetViewType dst_offset_;
      typedef typename OffsetViewType::non_const_value_type scalar_index_type;

      pack_functor (ViewType dst, ViewType src,
                    OffsetViewType dst_offset, OffsetViewType src_offset) :
        src_ (src),
        dst_ (dst),
        src_offset_ (src_offset),
        dst_offset_ (dst_offset)
      {}

      KOKKOS_INLINE_FUNCTION
      void operator () (const LocalOrdinal row) const {
        scalar_index_type srcPos = src_offset_(row);
        const scalar_index_type dstEnd = dst_offset_(row+1);
        scalar_index_type dstPos = dst_offset_(row);
        for ( ; dstPos < dstEnd; ++dstPos, ++srcPos) {
          dst_(dstPos) = src_(srcPos);
        }
      }
    };
  }; // class CrsMatrix

} // namespace Tpetra

#endif // TPETRA_KOKKOSREFACTOR_CRSMATRIX_DECL_HPP

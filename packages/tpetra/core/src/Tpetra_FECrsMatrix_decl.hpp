// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_FECRSMATRIX_DECL_HPP
#define TPETRA_FECRSMATRIX_DECL_HPP


/// \file Tpetra_FECrsMatrix_decl.hpp
/// \brief Declaration of the Tpetra::FECrsMatrix class

#include "Tpetra_CrsMatrix_decl.hpp"
#include "Tpetra_FECrsGraph.hpp"

namespace Tpetra {



/// \class FECrsMatrix
/// \brief Sparse matrix that presents a row-oriented interface that lets
///        users read or modify entries.
template<class Scalar        = ::Tpetra::Details::DefaultTypes::scalar_type,
         class LocalOrdinal  = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
         class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
         class Node          = ::Tpetra::Details::DefaultTypes::node_type>
class FECrsMatrix :
    public CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>
{
  friend class CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
public:
    //! @name Typedefs
    //@{

    /// \brief This class' first template parameter; the type of each
    ///        entry in the matrix.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::scalar_type scalar_type;

    /// \brief The type used internally in place of \c Scalar.
    ///
    /// Some \c Scalar types might not work with Kokkos on all
    /// execution spaces, due to missing CUDA device macros or
    /// volatile overloads.  The C++ standard type std::complex<T> has
    /// this problem.  To fix this, we replace std::complex<T> values
    /// internally with the (usually) bitwise identical type
    /// Kokkos::complex<T>.  The latter is the \c impl_scalar_type
    /// corresponding to \c Scalar = std::complex.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::impl_scalar_type impl_scalar_type;

    //! This class' second template parameter; the type of local indices.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_ordinal_type local_ordinal_type;

    //! This class' third template parameter; the type of global indices.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::global_ordinal_type global_ordinal_type;

    //! This class' fourth template parameter; the Kokkos device type.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::node_type node_type;

    //! The Kokkos device type.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::device_type device_type;

    //! The Kokkos execution space.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::execution_space execution_space;

    /// \brief Type of a norm result.
    ///
    /// This is usually the same as the type of the magnitude
    /// (absolute value) of <tt>Scalar</tt>, but may differ for
    /// certain <tt>Scalar</tt> types.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mag_type mag_type;

    //! The Map specialization suitable for this CrsMatrix specialization.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::map_type map_type;

    //! The Import specialization suitable for this CrsMatrix specialization
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::import_type import_type;

    //! The Export specialization suitable for this CrsMatrix specialization.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::export_type export_type;

    //! The CrsGraph specialization suitable for this CrsMatrix specialization.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::crs_graph_type crs_graph_type;

    //! The part of the sparse matrix's graph on each MPI process.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_graph_device_type local_graph_device_type;

    /// \brief The specialization of Kokkos::CrsMatrix that represents
    ///        the part of the sparse matrix for each MPI process on device.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_device_type local_matrix_device_type;

    /// \brief The specialization of Kokkos::CrsMatrix that represents
    ///        the part of the sparse matrix for each MPI process on host.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_host_type local_matrix_host_type;

    /// \brief Parent CrsMatrix type using the same scalars
    typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;

    //! The CrsGraph specialization suitable for this CrsMatrix specialization.
    typedef FECrsGraph<LocalOrdinal, GlobalOrdinal, Node> fe_crs_graph_type;

    //@}
    //! @name Constructors and destructor
    //@{


    /// \brief Constructor specifying one or two previously constructed graphs.
    ///
    /// Calling this constructor fixes the graph structure of the
    /// sparse matrix.  We say in this case that the matrix has a
    /// "static graph."  If you create a FECrsMatrix with this
    /// constructor, you are not allowed to insert new entries into
    /// the matrix, but you are allowed to change values in the
    /// matrix.
    ///
    /// The given graphs must be fill complete.  Note that calling
    /// resumeFill() on the graph makes it not fill complete, even if
    /// you had previously called fillComplete() on the graph.  In
    /// that case, you must call fillComplete() on the graph again
    /// before invoking this CrsMatrix constructor.
    ///
    /// The matrix is constructed in "closed" fill state and, if the the
    /// optional "start owned" parameter is set to true, "owned" state. You must
    /// first call beginAssembly or beginModify if you wish to change its
    /// entries (the difference between these two can be found elsewhere in the
    /// documentation).
    ///
    /// This constructor is marked \c explicit so that you can't
    /// create a FECrsMatrix by accident when passing a FECrsGraph into a
    /// function that takes a FECrsMatrix.
    ///
    /// \param graph [in] The FECrsGraph structure of the sparse matrix.
    ///   The graph <i>must</i> have endFill() called
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    explicit FECrsMatrix (const Teuchos::RCP<const fe_crs_graph_type>& graph,
                          const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    //! Copy constructor (forbidden).
    FECrsMatrix (const FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&) = delete;

    //! Move constructor (forbidden).
    FECrsMatrix (FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&&) = delete;

    //! Copy assignment (forbidden).
    FECrsMatrix&
    operator= (const FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&) = delete;
    //! Move assignment (forbidden).
    FECrsMatrix&
    operator= (FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&&) = delete;

    /// \brief Destructor (virtual for memory safety of derived classes).
    ///
    /// \note To Tpetra developers: See the C++ Core Guidelines C.21
    ///   ("If you define or <tt>=delete</tt> any default operation,
    ///   define or <tt>=delete</tt> them all"), in particular the
    ///   AbstractBase example, for why this destructor declaration
    ///   implies that we need the above four <tt>=delete</tt>
    ///   declarations for copy construction, move construction, copy
    ///   assignment, and move assignment.
    virtual ~FECrsMatrix () = default;

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
    void globalAssemble () {endFill();}

    //! Migrates data to the owned mode
    void endAssembly();

    /// \brief Activates the owned+shared mode for assembly
    ///
    /// This routine begins the finite element assembly in owned+shared mode.
    void beginAssembly();

    //! Closes modification phase
    void endModify();

    /// \brief Activates the owned mode for modifying local values
    ///
    /// This is for use in modifying local values *only*, not information
    /// which needs to be migrated across ranks.  Applying Dirichlet boundary
    /// conditions to owned degrees-of-freedom is the architypical example
    /// of a use-case for beginModify()
    void beginModify();

  private:

    //! Called by beginAssembly, endAssembly.
    void endFill();
    void beginFill();

    /// \brief Whether sumIntoLocalValues and sumIntoGlobalValues
    ///   should use atomic updates by default.
    ///
    /// \warning This is an implementation detail.
    static const bool useAtomicUpdatesByDefault =
#ifdef KOKKOS_ENABLE_SERIAL
      ! std::is_same<execution_space, Kokkos::Serial>::value;
#else
      true;
#endif // KOKKOS_ENABLE_SERIAL

  protected:

    //! Overloads of modification methods
    LocalOrdinal
    replaceGlobalValuesImpl (impl_scalar_type rowVals[],
                             const crs_graph_type& graph,
                             const RowInfo& rowInfo,
                             const GlobalOrdinal inds[],
                             const impl_scalar_type newVals[],
                             const LocalOrdinal numElts);

    LocalOrdinal
    replaceLocalValuesImpl (impl_scalar_type rowVals[],
                            const crs_graph_type& graph,
                            const RowInfo& rowInfo,
                            const LocalOrdinal inds[],
                            const impl_scalar_type newVals[],
                            const LocalOrdinal numElts);

    LocalOrdinal
    sumIntoGlobalValuesImpl (impl_scalar_type rowVals[],
                             const crs_graph_type& graph,
                             const RowInfo& rowInfo,
                             const GlobalOrdinal inds[],
                             const impl_scalar_type newVals[],
                             const LocalOrdinal numElts,
                             const bool atomic = useAtomicUpdatesByDefault);

    LocalOrdinal
    sumIntoLocalValuesImpl (impl_scalar_type rowVals[],
                            const crs_graph_type& graph,
                            const RowInfo& rowInfo,
                            const LocalOrdinal inds[],
                            const impl_scalar_type newVals[],
                            const LocalOrdinal numElts,
                            const bool atomic = useAtomicUpdatesByDefault);

    void
    insertGlobalValuesImpl (crs_graph_type& graph,
                            RowInfo& rowInfo,
                            const GlobalOrdinal gblColInds[],
                            const impl_scalar_type vals[],
                            const size_t numInputEnt);

  protected:
    /// \brief Migrate data from the owned+shared to the owned matrix
    /// Since this is non-unique -> unique, we need a combine mode.
    /// Precondition: Must be FE_ACTIVE_OWNED_PLUS_SHARED mode
    void doOwnedPlusSharedToOwned(const CombineMode CM=Tpetra::ADD);

    /// \brief Migrate data from the owned to the owned+shared matrix
    /// Precondition: Must be FE_ACTIVE_OWNED mode
    void doOwnedToOwnedPlusShared(const CombineMode CM=Tpetra::ADD);

    //! Switches which CrsGraph is active (without migrating data)
    void switchActiveCrsMatrix();
    //@}

  private:
    // Enum for activity
    enum FEWhichActive
    {
      FE_ACTIVE_OWNED,
      FE_ACTIVE_OWNED_PLUS_SHARED
    };

    // The FECrsGraph from construction time
    Teuchos::RCP<const FECrsGraph<LocalOrdinal, GlobalOrdinal, Node> > feGraph_;

    // This is whichever multivector isn't currently active
    Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > inactiveCrsMatrix_;
    // This is in RCP to make shallow copies of the FECrsMatrix work correctly
    Teuchos::RCP<FEWhichActive> activeCrsMatrix_;

    enum class FillState
    {
      open,  // matrix is "open".  Values can freely summed in to and replaced
      modify,  // matrix is open for modification.  *local* values can be replaced
      closed
    };
    Teuchos::RCP<FillState> fillState_;

};    // end class FECrsMatrix



}     // end namespace Tpetra



#endif // TPETRA_FECRSMATRIX_DECL_HPP

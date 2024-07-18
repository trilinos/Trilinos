// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_FEMULTIVECTOR_DECL_HPP
#define TPETRA_FEMULTIVECTOR_DECL_HPP

/// \file Tpetra_FEMultiVector_decl.hpp
/// \brief Declaration of the Tpetra::MultiVector class

#include "Tpetra_FEMultiVector_fwd.hpp"
#include "Tpetra_MultiVector_decl.hpp"

namespace Tpetra {

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class FEMultiVector :
    public MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>
  {
  private:
    using base_type = ::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    friend base_type;

    //! Declare the beginning of a phase of owned+shared modifications.
    void beginFill ();

    //! Declare the end of a phase of owned+shared modifications.
    void endFill ();

  public:
    //! @name Typedefs to facilitate template metaprogramming.
    //@{

    //! The type of each entry of the object.
    using scalar_type = Scalar;
    //! The type of the object's local indices.
    using local_ordinal_type = LocalOrdinal;
    //! The type of the object's global indices.
    using global_ordinal_type = GlobalOrdinal;

    //! The object's Kokkos::Device specialization.
    using device_type = typename base_type::device_type;

    //! The object's Kokkos execution space.
    using execution_space = typename base_type::execution_space;

    //! Legacy typedef that will eventually disappear.
    using node_type = Node;

    //! Specialization of dual_view_type that this object may use.
    using dual_view_type = typename base_type::dual_view_type;

    //! The type of the Map specialization used by this class.
    using map_type = typename base_type::map_type;

    //! The storage type of each entry of the object.
    using impl_scalar_type = typename base_type::impl_scalar_type;

    //! The type of this object's dot product result.
    using dot_type = typename base_type::dot_type;

    //! The type of this object's norm result.
    using mag_type = typename base_type::mag_type;

    //@}
    //! @name Constructors and destructor
    //@{

    //! Default constructor (forbidden).
    FEMultiVector () = delete;

    /// \brief Basic constructor.
    ///
    /// The vector is constructed in "owned plus shared" and "closed" fill state.
    /// You must first call beginAssembly or beginModify if you wish to change
    /// its entries (the difference between these two can be found elsewhere in
    /// the documentation).
    ///
    /// \param map [in] Map describing the distribution of rows of the
    ///   resulting MultiVector.  If the importer is not null, this
    ///   must be the same as importer->getSourceMap().
    ///
    /// \param importer [in] Import describing the distribution of
    ///   rows of the the multivector in two separate modes.  In the
    ///   case of finite-element assembly, importer->getSourceMap()
    ///   should correspond to the uniquely owned entries in the
    ///   domain Map of the finite element problem.
    ///   importer->getTargetMap() should correspond to the overlapped
    ///   entries needed for assembly.  The canonical way to get this
    ///   Import object is to take the Import object associated with
    ///   the CrsGraph for your finite element matrix.
    ///
    /// \param numVectors [in] Number of columns of the resulting
    ///   MultiVector.
    ///
    /// \note A Map AND an Import need to be arguments because in
    ///   serial, the importer will be null.  This will default to
    ///   importer->getTargetMap() being the active MultiVector.
    FEMultiVector (const Teuchos::RCP<const map_type>& map,
                   const Teuchos::RCP<const Import<local_ordinal_type, global_ordinal_type, node_type>>& importer,
                   const size_t numVecs,
                   const bool zeroOut = true);

    //! Copy constructor (forbidden).
    FEMultiVector (const FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&) = delete;

    //! Move constructor (forbidden).
    FEMultiVector (FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&&) = delete;

    //! Copy assigment (forbidden).
    FEMultiVector&
    operator= (const FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&) = delete;

    //! Move assigment (forbidden).
    FEMultiVector&
    operator= (FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&&) = delete;

    /// \brief Destructor (virtual for memory safety of derived classes).
    ///
    /// \note To Tpetra developers: See the C++ Core Guidelines C.21
    ///   ("If you define or <tt>=delete</tt> any default operation,
    ///   define or <tt>=delete</tt> them all"), in particular the
    ///   AbstractBase example, for why this destructor declaration
    ///   implies that we need the above four <tt>=delete</tt>
    ///   declarations for copy construction, move construction, copy
    ///   assignment, and move assignment.
    virtual ~FEMultiVector () = default;

    //@}
    //! @name Post-construction modification routines
    //@{

    //! Declare the beginning of a phase of owned+shared modifications.
    void beginAssembly ();

    //! Declare the end of a phase of owned+shared modifications.
    void endAssembly ();

    void beginModify ();
    void endModify ();

    /// \brief Declare the end of a phase of owned+shared
    ///   modifications; same as endFill().
    void globalAssemble ();

    /// \brief Migrate data from the owned+shared to the owned multivector
    /// Since this is non-unique -> unique, we need a combine mode.
    /// Precondition: Must be FE_ACTIVE_OWNED_PLUS_SHARED mode
    /// \warning This is for EXPERT USE ONLY.  We make NO PROMISES of
    ///   backwards compatibility.
    void doOwnedPlusSharedToOwned(const CombineMode CM=Tpetra::ADD);

    /// \brief Migrate data from the owned to the owned+shared multivector
    /// Precondition: Must be FE_ACTIVE_OWNED mode
    /// \warning This is for EXPERT USE ONLY.  We make NO PROMISES of
    ///   backwards compatibility.
    void doOwnedToOwnedPlusShared(const CombineMode CM=Tpetra::ADD);

    //! Switches which Multivector is active (without migrating data)
    void switchActiveMultiVector();

  protected:
    /// \brief Replace the underlying Map in place.
    ///
    /// \warning FEMultiVector does not allow this and will throw if
    ///   you call this method.
    void replaceMap (const Teuchos::RCP<const map_type>& map);

    //! Enum for activity
    enum FEWhichActive
    {
      FE_ACTIVE_OWNED_PLUS_SHARED,
      FE_ACTIVE_OWNED
    };

    enum class FillState
    {
      open,  // matrix is "open".  Values can freely summed in to and replaced
      modify,  // matrix is open for modification.  *local* values can be replaced
      closed
    };
    Teuchos::RCP<FillState> fillState_;

    //! Whichever MultiVector is <i>not</i> currently active.
    Teuchos::RCP<base_type> inactiveMultiVector_;

    /// \brief Whether the owned MultiVector or the owned plus shared
    ///   MultiVector is active.
    ///
    /// This is an RCP in order to make shallow copies of the
    /// FEMultiVector work correctly.
    Teuchos::RCP<FEWhichActive> activeMultiVector_;

    //! Import object used for communication between the two MultiVectors.
    Teuchos::RCP<const Import<local_ordinal_type, global_ordinal_type, node_type>> importer_;
  };

} // namespace Tpetra

#endif // TPETRA_FEMULTIVECTOR_DECL_HPP

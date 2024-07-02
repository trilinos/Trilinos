// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_BLOCKTRIDICONTAINER_DECL_HPP
#define IFPACK2_BLOCKTRIDICONTAINER_DECL_HPP

/// \file Ifpack2_BlockTriDiContainer_decl.hpp
/// \brief Ifpack2::BlockTriDiContainer class declaration

#include "Ifpack2_config.h"
#include "Ifpack2_Container.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_BlockCrsMatrix_decl.hpp"
#include <type_traits>
#include <string>

namespace Ifpack2 {

  /// \class BlockTriDiContainer
  /// \brief Store and solve local block tridiagonal linear problems.
  /// \tparam MatrixType A specialization of Tpetra::RowMatrix.
  ///
  /// This class can be used as a Container for BlockRelaxation, in which case
  /// please refer to the documentation of the Container interface, or
  /// standalone. In standalone use, there are special constructor,
  /// <tt>compute</tt>, and <tt>applyInverseJacobi</tt> functions that may be
  /// called, with extra non-ParameterList inputs.
  ///
  /// If the partitioner returns an empty set of partitions, then this class
  /// performs pure block-Jacobi, i.e., where the preconditioner is the diagonal
  /// of little blocks or, in other words, block tridiagonal matrices of size one
  /// block.
  ///
  /// BlockTriDiContainer requires the Tpetra::RowMatrix to be
  /// Tpetra::BlockCrsMatrix. 
  ///
  /// This class currently assumes the following about the column and
  /// row Maps of the input matrix:
  /// <ol>
  /// <li> The domain and range maps are the same.
  /// <li> On all processes, all off-process indices in the column Map of the
  ///      input matrix occur after that initial set.</li>
  /// </ol>
  /// These assumptions may be violated if the input matrix was constructed with a
  /// user-provided column Map.
  ///
  /// Currently, this class is expected to perform well on conventional CPU and
  /// Intel Xeon Phi (reaching the ceiling of the bandwidth utilization) and 
  /// perform reasonably well on GPU (comparable to Intel Xeon Phi performance).
  /// The main performance issue on GPU is block sparse matrix vector multiplication
  /// which does not use a SIMD format.
  ///
  /// Implementation specific comments:
  ///  - When ETI is not enabled, do not use something like "using namepsace KokkosBatched::Experimental".
  ///    Albany (or any applications) can use the same struct name (in this case Albany uses Side).
  ///  - Use an impl pointer to hide details. If you use an object inside of an Ifpack container,  
  ///    it requires a complete definition of the member object, which needs to expose an impl details 
  ///    header. However, a pointer does not require a complete definition of the member objects.
  ///  - Always test with complex even if this code is not used with complex. 
  ///  - Always check a non MPI build to check MPI is guarded by #ifdef HAVE_IFPACK2_MPI
  ///  - Do not trust CMake variables and macros because you see the variables in the file. If you use
  ///    CMake varialbes and macro definitions, check Ifpack2_config.h. 
  ///  - Always better remove warnings (shadows and signed/unsinged comparison). If the code is used by 
  ///    other customers, they may have a different software quality standard. It is better to follow 
  ///    a higher quality standard.

  ///
  /// Impl Tag
  ///
  namespace BlockTriDiContainerDetails {
    ///
    /// impl tag to distinguish built-in types and sacado types
    ///
    struct ImplNotAvailTag {};
    struct ImplSimdTag {};
    struct ImplSacadoTag {};

    template<typename T> struct ImplTag                        { typedef ImplNotAvailTag type; };
    template<>           struct ImplTag<float>                 { typedef ImplSimdTag type;     };
    template<>           struct ImplTag<double>                { typedef ImplSimdTag type;     };
    template<>           struct ImplTag<std::complex<float> >  { typedef ImplSimdTag type;     };
    template<>           struct ImplTag<std::complex<double> > { typedef ImplSimdTag type;     };

    /// forward declaration 
    template<typename MatrixType> struct ImplObject;
  }
  
  ///
  /// Primary declation
  ///
  template <typename MatrixType, 
            typename ImplTagType = typename BlockTriDiContainerDetails::ImplTag<typename MatrixType::scalar_type>::type>
  class BlockTriDiContainer;
  
  ///
  /// Partial specialization with SIMD<ScalarType> internally used.
  /// The code does not support Sacado UQ types. As UQ type also uses SIMD like structure,
  /// it needs an additional specialization.
  ///
  template <typename MatrixType>
  class BlockTriDiContainer<MatrixType,BlockTriDiContainerDetails::ImplSimdTag> 
    : public Container<MatrixType> {
    //! @name Internal typedefs (private)
    //@{
  private:
    /// \brief The first template parameter of this class.
    ///
    /// This must be a Tpetra::RowMatrix specialization.  It may have
    /// entirely different template parameters (e.g., \c scalar_type)
    /// than \c InverseType.
    typedef MatrixType matrix_type;
      
    //! The type of entries in the input (global) matrix.
    typedef typename MatrixType::scalar_type scalar_type;
    //! The magnitude of entries in the input (global) matrix.
    typedef typename Kokkos::ArithTraits<scalar_type>::magnitudeType magnitude_type;
    //! The type of local indices in the input (global) matrix.
    typedef typename Container<MatrixType>::local_ordinal_type local_ordinal_type;
    //! The type of global indices in the input (global) matrix.
    typedef typename Container<MatrixType>::global_ordinal_type global_ordinal_type;
    //! The Node type of the input (global) matrix.
    typedef typename Container<MatrixType>::node_type node_type;

    typedef typename Container<MatrixType>::mv_type mv_type;
    typedef typename Container<MatrixType>::map_type map_type;
    typedef typename Container<MatrixType>::vector_type vector_type;
    typedef typename Container<MatrixType>::import_type import_type;

    typedef typename Container<MatrixType>::HostView host_view_type;
    typedef typename Container<MatrixType>::ConstHostView const_host_view_type;
    typedef host_view_type HostView;
    typedef const_host_view_type ConstHostView;
    //typedef Tpetra::MultiVector<local_scalar_type, local_ordinal_type, global_ordinal_type, node_type> local_mv_type;
    //typedef typename Kokkos::View<local_scalar_type**, Kokkos::HostSpace> HostViewLocal;

    typedef Tpetra::CrsMatrix
    <scalar_type,local_ordinal_type,global_ordinal_type,node_type> crs_matrix_type;
    typedef Tpetra::BlockCrsMatrix
    <scalar_type,local_ordinal_type,global_ordinal_type,node_type> block_crs_matrix_type;

    const Teuchos::Array<Teuchos::Array<local_ordinal_type> > partitions_;

    /// \brief The (base class) type of the input matrix.
    ///
    /// The input matrix to the constructor must be a Tpetra::BlockCrsMatrix.
    /// However, we want to make the constructor as general as possible, so we
    /// always accept the matrix as a Tpetra::RowMatrix.  This typedef is the
    /// appropriate specialization of Tpetra::RowMatrix.
    typedef typename Container<MatrixType>::row_matrix_type row_matrix_type;

    static_assert (std::is_same<MatrixType, row_matrix_type>::value,
                   "Ifpack2::BlockTriDiContainer: MatrixType must be a Tpetra::RowMatrix specialization.");
    //@}
  public:
    //! \name Constructor and destructor
    //@{

    /// \brief Constructor.
    ///
    /// \brief matrix [in] The original input matrix.  This Container
    ///   will construct a local diagonal block from the rows given by
    ///   <tt>localRows</tt>.
    ///
    /// \param partitions [in] The set of (local) rows assigned to this
    ///   container.  <tt>partitions[i] == j</tt>, where i (from 0 to
    ///   <tt>getNumRows() - 1</tt>) indicates the SparseContainer's
    ///   row, and j indicates the local row in the calling process.
    ///   <tt>partitions.size()</tt> gives the number of rows in the
    ///   local matrix on each process.  This may be different on
    ///   different processes. If partitions is empty, then a local
    ///   block-Jacobi preconditioner is formed; this is equivalent to
    ///   every part of partitions having size 1.
    BlockTriDiContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                         const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                         const Teuchos::RCP<const import_type>& importer,
                         bool pointIndexed);

    /// \brief Constructor for bare construction.
    ///
    /// Constructor for users who want to use BlockTriDiContainer directly. See
    /// main constructor for documentation. This constructor removes the
    /// Container-general arguments that are not used in BlockTriDiContainer.
    ///
    /// \param overlapCommAndComp [in] Overlap communication and computation. This
    ///   is not always better; it depends on (at least) the MPI implementation
    ///   and the machine architecture. Defaults to false. It has to be specified
    ///   at construction because the core data structures depend on whether
    ///   overlapCommAndComp is true or false.
    ///
    /// \param useSequentialMethod [in] This is for development and testing
    ///   purposes only.
    BlockTriDiContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                         const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                         const int n_subparts_per_part = 1,
                         bool overlapCommAndComp = false, 
                         bool useSequentialMethod = false,
                         const int block_size = -1,
                         const bool explicitConversion = false);

    //! Destructor (declared virtual for memory safety of derived classes).
    ~BlockTriDiContainer () override;

    struct ComputeParameters {
      //! addRadiallyToDiagonal [in] Add a constant to each diagonal entry of the
      //! matrix. It is added according to the (possibly complex) sign of the
      //! digonal entry at the time that the entry is the pivot in the LU
      //! factorization, such that the entry is moved radially outward in the
      //! complex plane. N.B. that this constant modifies the matrix in the linear
      //! equation, not simply the diagonal preconditioner.
      magnitude_type addRadiallyToDiagonal = Kokkos::ArithTraits<magnitude_type>::zero();
    };

    //! Input arguments to <tt>applyInverseJacobi</tt>
    struct ApplyParameters {
      //! Set this to true if <tt>Y</tt> is meant to be treated as 0 on
      //! entry. Defaults to false.
      bool zeroStartingSolution = false;
      //! Damping factor. Defaults to 1.
      scalar_type dampingFactor = Kokkos::ArithTraits<scalar_type>::one();
      //! The maximum number of sweeps. If the norm-based criterion is not used,
      //! it's exactly the number of sweeps. Defaults to 1.
      int maxNumSweeps = 1;
      //! The 2-norm-based termination criterion. If it is larger than 0, then it
      //! is used. Termination is then based on <tt>maxNumSweeps</tt> and this
      //! criterion, whichever is reached first. In the case of multiple right
      //! hand sides, norm-based termination occurs only when each vector
      //! satisfies the criterion. The criterion for a vector is <tt>f(D^{-1} (x -
      //! A*y)) <= tolerance * f(D^{-1} (x - R*y_0))</tt>, where <tt>y_0</tt> is
      //! the input <tt>y</tt>, often 0, and <tt>f</tt> is the maximum of the
      //! 2-norms of each degree of freedom, i.e., index of a block. Defaults to
      //! 0.
      magnitude_type tolerance = Kokkos::ArithTraits<magnitude_type>::zero();
      //! Check the norm-based termination criterion every
      //! <tt>checkToleranceEvery</tt> iterations. Defaults to 1. A norm
      //! computation requires a global reduction, which is expensive. Hence it
      //! can make sense to check for termination only every 2 or 3 iterations if
      //! it is anticipated that the number of iterations required to terminate
      //! will be sufficiently large to make 1 or 2 extra iterations relatively
      //! cheap. Defaults to 1.
      int checkToleranceEvery = 1;
    };

    //@}
    //! \name Get and set methods
    //@{

    //! Set all necessary parameters.
    void setParameters(const Teuchos::ParameterList& List) override;

    void clearBlocks() override;

    //@}
    //! \name Mathematical functions
    //@{

    //! Do all set-up operations that only require matrix structure.
    void initialize () override;

    //! Extract the local tridiagonal block and prepare the solver.
    void compute () override;

    // \brief Compute <tt>Y := D^{-1} (X - R*Y)</tt>.
    void applyInverseJacobi (const mv_type& X, mv_type& Y,
                             scalar_type dampingFactor,
                             bool zeroStartingSolution = false,
                             int numSweeps = 1) const override;

    /// \brief Create a ComputeParameters struct with default values.
    ComputeParameters createDefaultComputeParameters () const;

    /// \brief Extract the local tridiagonal block and prepare the solver.
    ///
    /// This version of <tt>applyInverseJacobi</tt> is meant to be called by
    /// direct users of this class, rather than by <tt>BlockRelaxation</tt>.
    ///
    /// \param addRadiallyToDiagonal [in] Add a constant to each diagonal entry of
    ///   the matrix. It is added according to the (possibly complex) sign of the
    ///   digonal entry at the time that the entry is the pivot in the LU
    ///   factorization, such that the entry is moved radially outward in the
    ///   complex plane. N.B. that this constant modifies the matrix in the linear
    ///   equation, not simply the diagonal preconditioner.
    void compute (const ComputeParameters& input);

    /// \brief Create an ApplyParameters struct with default values.
    ApplyParameters createDefaultApplyParameters () const;

    /// \brief Compute <tt>Y := D^{-1} (X - R*Y)</tt>.
    ///
    /// This version of <tt>applyInverseJacobi</tt> is meant to be called by
    /// direct users of this class, rather than by <tt>BlockRelaxation</tt>. It
    /// supports 2-norm-based termination. It returns the number of sweeps
    /// performed.
    int applyInverseJacobi (const mv_type& X, mv_type& Y, 
                            const ApplyParameters& input) const;

    /// \brief If a norm-based method was used, return a L2 norm of all rhs 
    ///        at the first iteration; otherwise return a minus one indicating
    ///        norm is not requested.
    const magnitude_type getNorms0 () const;

    /// \brief If a norm-based method was used, return a L2 norm of all rhs;
    ///        otherwise return zero.
    const magnitude_type getNormsFinal () const;
  
    //! Compute <tt>Y := (1 - a) Y + a D^{-1} (X - R*Y)</tt>. Not supported. Call
    //! <tt>applyInverseJacobi</tt> instead.
    void
    apply (const_host_view_type X,
           host_view_type Y,
           int blockIndex,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
           scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const override;

    //! Compute <tt>Y := alpha * diag(D) * M^{-1} (diag(D) * X) + beta*Y</tt>. Not
    //! supported.
    void
    weightedApply (const_host_view_type X,
                   host_view_type Y,
                   const_host_view_type W,
                   int blockIndex,
                   Teuchos::ETransp mode = Teuchos::NO_TRANS,
                   scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
                   scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const override;

    //@}
    //! \name Miscellaneous methods
    //@{

    /// \brief Print information about this object to the given output stream.
    ///
    /// operator<< uses this method.
    std::ostream& print (std::ostream& os) const override;

    //@}
    //! @name Implementation of Teuchos::Describable
    //@{

    //! A one-line description of this object.
    std::string description () const override;

    //! Print the object with some verbosity level to the given FancyOStream.
    void
    describe (Teuchos::FancyOStream &out,
              const Teuchos::EVerbosityLevel verbLevel =
              Teuchos::Describable::verbLevel_default) const override;

    //@}

    /// \brief Get the name of this container type for Details::constructContainer()
    static std::string getName();

  private:
    //! Copy constructor: Declared but not implemented, to forbid copy construction.
    BlockTriDiContainer (const BlockTriDiContainer<MatrixType>& rhs);

    // hide details of impl using ImplObj; finally I understand why AMB did that way.
    Teuchos::RCP<BlockTriDiContainerDetails::ImplObject<MatrixType> > impl_;
    int n_subparts_per_part_;
    int block_size_ = -1;
    
    // initialize distributed and local objects
    void initInternal (const Teuchos::RCP<const row_matrix_type>& matrix,
                       const Teuchos::RCP<const import_type> &importer,
                       const bool overlapCommAndComp,
                       const bool useSeqMethod,
                       const int block_size = -1,
                       const bool explicitConversion = false);

    void clearInternal();
  };
  
  ///
  /// ImplNotAvailTag
  ///   This container does not support UQ types; however, the UQ types are required for 
  ///   Stokhos ETI. To prevent linking errors, we provide an empty implementation with
  ///   ImplNotAvailTag. Upon the request to support UQ types, we need to specialize
  ///   the impl function and interface with ImplSacadoTag.
  ///
  template <typename MatrixType>
  class BlockTriDiContainer<MatrixType,BlockTriDiContainerDetails::ImplNotAvailTag> 
    : public Container<MatrixType> {
  private:
    typedef typename MatrixType::scalar_type scalar_type;
    typedef typename Kokkos::ArithTraits<scalar_type>::magnitudeType magnitude_type;
    typedef typename Container<MatrixType>::local_ordinal_type local_ordinal_type;
    typedef typename Container<MatrixType>::global_ordinal_type global_ordinal_type;

    typedef typename Container<MatrixType>::mv_type mv_type;
    typedef typename Container<MatrixType>::import_type import_type;

    typedef typename Container<MatrixType>::HostView host_view_type;
    typedef typename Container<MatrixType>::ConstHostView const_host_view_type;
    typedef typename Container<MatrixType>::row_matrix_type row_matrix_type;

    static_assert (std::is_same<MatrixType, row_matrix_type>::value,
                   "Ifpack2::BlockTriDiContainer: MatrixType must be a Tpetra::RowMatrix specialization.");
  public:

    BlockTriDiContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                         const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                         const Teuchos::RCP<const import_type>& importer,
                         bool pointIndexed)
      : Container<MatrixType>(matrix, partitions, pointIndexed) {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Error: BlockTriDiContainer is not available for this scalar_type");
    }

    void setParameters(const Teuchos::ParameterList& List) override {}
    void clearBlocks() override {}

    void initialize () override {}
    void compute () override {}
    void applyInverseJacobi (const mv_type& X, mv_type& Y,
                             scalar_type dampingFactor,
                             bool zeroStartingSolution = false,
                             int numSweeps = 1) const override {}
    
    void
    apply (const_host_view_type X,
           host_view_type Y,
           int blockIndex,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
           scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const override {}

    void
    weightedApply (const_host_view_type X,
                   host_view_type Y,
                   const_host_view_type W,
                   int blockIndex,
                   Teuchos::ETransp mode = Teuchos::NO_TRANS,
                   scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
                   scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const override {}

    std::ostream& print (std::ostream& os) const override { 
      return os << "Ifpack2::BlockTriDiContainer::ImplNotAvailTag"; 
    }

    std::string description () const override {
      return "Ifpack2::BlockTriDiContainer::ImplNotAvailTag";
    }

    void
    describe (Teuchos::FancyOStream &out,
              const Teuchos::EVerbosityLevel verbLevel =
              Teuchos::Describable::verbLevel_default) const override {
      out << "Ifpack2::BlockTriDiContainer::ImplNotAvailTag";
    }
    
    static std::string getName() {
      return "Ifpack2::BlockTriDiContainer::ImplNotAvailTag";
    }
  };


} // namespace Ifpack2

#endif // IFPACK2_BLOCKTRIDICONTAINER_DECL_HPP

// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef AMESOS2_SOLVER_UQ_PCE_HPP
#define AMESOS2_SOLVER_UQ_PCE_HPP

#include "Amesos2_Solver.hpp"
#include "Amesos2_Factory.hpp"
#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"
#include "Stokhos_Tpetra_Utilities_UQ_PCE.hpp"

namespace Amesos2 {

  template <class S, class LO, class GO, class NO>
  typename Sacado::UQ::PCE<S>::cijk_type
  get_pce_cijk(
    const Teuchos::RCP<const Tpetra::CrsMatrix<Sacado::UQ::PCE<S>, LO, GO, NO > >& A = Teuchos::null,
    const Teuchos::RCP<Tpetra::MultiVector<Sacado::UQ::PCE<S>, LO, GO, NO > >& X = Teuchos::null,
    const Teuchos::RCP<const Tpetra::MultiVector<Sacado::UQ::PCE<S>, LO, GO, NO > >& B = Teuchos::null)
  {
    if (A != Teuchos::null) {
      return Kokkos::cijk(A->getLocalValuesDevice(Tpetra::Access::ReadOnly));
    }
    else if (X != Teuchos::null) {
      return Kokkos::cijk(X->getLocalViewDevice(Tpetra::Access::ReadOnly));
    }
    else if (B != Teuchos::null) {
      return Kokkos::cijk(B->getLocalViewDevice(Tpetra::Access::ReadOnly));
    }
    return typename Sacado::UQ::PCE<S>::cijk_type();
  }

  /// \brief Amesos2 solver adapter for UQ::PCE scalar type
  ///
  /// This adapter enables Amesos2 solvers to work with Tpetra matrices
  /// and vectors of the Sacado::UQ::PCE scalar type by "flattening"
  /// these matrices and vectors into ones with a standard (e.g., double)
  /// scalar type.
  template <class Storage, class LocalOrdinal, class GlobalOrdinal,
            class Node, template<class,class> class ConcreteSolver>
  class PCESolverAdapter :
    public Solver< Tpetra::CrsMatrix<Sacado::UQ::PCE<Storage>,
                                     LocalOrdinal,
                                     GlobalOrdinal,
                                     Node >,
                   Tpetra::MultiVector<Sacado::UQ::PCE<Storage>,
                                       LocalOrdinal,
                                       GlobalOrdinal,
                                       Node >
                   >
  {
  public:

    typedef Sacado::UQ::PCE<Storage> Scalar;
    typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Matrix;
    typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Vector;

    typedef typename Scalar::value_type BaseScalar;
    typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Map;
    typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> FlatGraph;
    typedef Tpetra::CrsMatrix<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> FlatMatrix;
    typedef Tpetra::MultiVector<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> FlatVector;
    typedef ConcreteSolver<FlatMatrix,FlatVector> FlatConcreteSolver;
    typedef Solver<FlatMatrix,FlatVector> FlatSolver;

    typedef Solver<Matrix,Vector> solver_type;
    typedef typename solver_type::type type;
    typedef typename Scalar::cijk_type cijk_type;

    /// Constructor
    PCESolverAdapter(
      const Teuchos::RCP<const Matrix>& A_,
      const Teuchos::RCP<Vector>& X_,
      const Teuchos::RCP<const Vector>& B_) :
      A(A_), X(X_), B(B_) {
      cijk = get_pce_cijk(A, X, B);
      const size_t pce_size =
        Kokkos::dimension_scalar(A->getLocalMatrixDevice().values);
      flat_graph =
        Stokhos::create_flat_pce_graph(*(A->getCrsGraph()),
                                       cijk,
                                       flat_X_map,
                                       flat_B_map,
                                       cijk_graph,
                                       pce_size);
      if (A != Teuchos::null)
        flat_A = Stokhos::create_flat_matrix(*A, flat_graph, cijk_graph, cijk);
      if (X != Teuchos::null)
        flat_X = Stokhos::create_flat_vector_view(*X, flat_X_map);
      if (B != Teuchos::null)
        flat_B = Stokhos::create_flat_vector_view(*B, flat_B_map);
      flat_solver =
        create_solver_with_supported_type<ConcreteSolver,FlatMatrix,FlatVector>::apply(flat_A, flat_X, flat_B);
    }

    /// \name Mathematical Functions
    //@{

    /** \brief Pre-orders the matrix.
     *
     * Uses the default solver option, unless a solver-specific
     * pre-ordering parameter is given.
     *
     * \sa setParameters
     */
    virtual type& preOrdering( void ) {
      flat_solver->preOrdering();
      return *this;
    }


    /** \brief Performs symbolic factorization on the matrix
     *
     * \pre
     *  - The matrix A must not be \c null
     */
    virtual type& symbolicFactorization( void ) {
      flat_solver->symbolicFactorization();
      return *this;
    }


    /** \brief Performs numeric factorization on the matrix
     *
     * numericFactorization checks first that symbolicFactorization
     * has successfully been called, and if not, calls it before
     * continuing.
     *
     * \pre
     *  - The matrix A must not be \c null
     *
     * \post
     *  - The factors L and U of A are computed
     */
    virtual type& numericFactorization( void ) {
      flat_solver->numericFactorization();
      return *this;
    }


    /** \brief Solves \f$ A X = B\f$ (or \f$ A^T X = B\f$ )
     *
     * solve checks first that numericFactorization has successfully
     * been called, and if not, calls it before continuing.
     *
     * \pre
     *  - The (multi)vectors \c X and \c B must not be \c null
     *
     * \post
     *  - The (multi)vector \c X (given at construction time) contains
     *    the solution to the system.
     */
    virtual void solve( void ) {
      flat_solver->solve();
    }


    /** \brief Solve \f$ A X = B\f$ using the given X and B vectors.
     *
     * This overload of solve uses the given X and B vectors when
     * solving.  This X and B are used in place of any X and B that
     * were given upon construction of the Amesos2 solver instance and
     * are used only for this solve.
     *
     * If a permanent change of X and B are required, see the setX()
     * and setB() methods.
     *
     * \post
     *  - The (multi)vector \c XX contains the solution to the system
     *  - The \c XX and \c BB given at construction time (if any) are unchanged.
     */
    virtual void solve(const Teuchos::Ptr<Vector>       XX,
                       const Teuchos::Ptr<const Vector> BB) const {
      flat_solver->solve(
        Stokhos::create_flat_vector_view(*XX, flat_X_map).get(),
        Stokhos::create_flat_vector_view(*BB, flat_B_map).get() );
    }


    /** \brief Solve \f$ A X = B\f$ using the given X and B vectors.
     *
     * This overload of solve uses the given X and B vectors when
     * solving.  This X and B are used in place of any X and B that
     * were given upon construction of the Amesos2 solver instance and
     * are used only for this solve.
     *
     * If a permanent change of X and B are required, see the setX()
     * and setB() methods.
     *
     * \post
     *  - The (multi)vector \c XX contains the solution to the system
     *  - The \c XX and \c BB given at construction time (if any) are unchanged.
     */
    virtual void solve(Vector* XX, const Vector* BB) const {
      flat_solver->solve(
        Stokhos::create_flat_vector_view(*XX, flat_X_map).get(),
        Stokhos::create_flat_vector_view(*BB, flat_B_map).get() );
    }

    //@} End Mathematical Functions


    /** \name Parameter Methods
     * @{
     */

    /** \brief Set/update internal variables and solver options.
     *
     * Expects that parameterList be named "Amesos2".  That list may
     * contain Amesos2-specific parameters.  In addition, it may
     * contain sublist for solver-specific parameters.  These sublists
     * should be named according to what is returned by the name()
     * function (i.e. The solver's name when enabling for Amesos2
     * during configuration).
     *
     * See each solver interface directly for a list of the supported
     * parameters for that solver.
     */
    virtual type& setParameters(
      const Teuchos::RCP<Teuchos::ParameterList> & parameterList ) {
      flat_solver->setParameters(parameterList);
      return *this;
    }


    /**
     * \brief Return a const parameter list of all of the valid parameters that
     * this->setParameterList(...)  will accept.
     */
    virtual Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters( void ) const {
      return flat_solver->getValidParameters();
    }

    /// @} End Parameter Methods


    /** \name Accessor Methods
     * @{
     */

    /** \brief Sets the matrix A of this solver
     *
     * \param [in] a          An RCP to a matrix will will be used for
     *                        future computation steps
     *
     * \param [in] keep_phase This parameter tells the solver what
     *                        state it should keep.  For example, you
     *                        may want to replace the matrix but keep
     *                        the symbolic factorization because you
     *                        know the structure of the new matrix is
     *                        the same as the structure of the old
     *                        matrix.  In this case you would pass
     *                        Amesos2::SYMBFACT as this parameter.
     *
     * The default value for the second parameter is Amesos2::CLEAN,
     * which means that the internal state of the solver will be
     * completely reset.  It will be as if no previous computational
     * steps were performed.
     */
    virtual void setA( const Teuchos::RCP<const Matrix> a,
                       EPhase keep_phase = CLEAN ) {
      A = a;

      // Rebuild flat matrix/graph
      cijk = get_pce_cijk(A);
      if (keep_phase <= CLEAN) {
        flat_X_map = Teuchos::null;
        flat_B_map = Teuchos::null;
        flat_graph = Teuchos::null;
        flat_graph =
          Stokhos::create_flat_pce_graph(*(A->getCrsGraph()),
                                         cijk,
                                         flat_X_map,
                                         flat_B_map,
                                         cijk_graph,
                                         Kokkos::dimension_scalar(A->getLocalMatrixDevice().values));
      }
      if (keep_phase <= SYMBFACT) // should this by NUMFACT???
        flat_A = Stokhos::create_flat_matrix(*a, flat_graph, cijk_graph, cijk);

      flat_solver->setA(flat_A, keep_phase);
    }

    /** \brief Sets the matrix A of this solver
     *
     * \param [in] a          An raw C pointer to a matrix will will
     *                        be used for future computation steps.
     *
     * \param [in] keep_phase This parameter tells the solver what
     *                        state it should keep.  For example, you
     *                        may want to replace the matrix but keep
     *                        the symbolic factorization because you
     *                        know the structure of the new matrix is
     *                        the same as the structure of the old
     *                        matrix.  In this case you would pass
     *                        Amesos2::SYMBFACT as this parameter.
     *
     * The default value for the second parameter is Amesos2::CLEAN,
     * which means that the internal state of the solver will be
     * completely reset.  It will be as if no previous computational
     * steps were performed.
     */
    virtual void setA( const Matrix* a, EPhase keep_phase = CLEAN ) {
      this->setA(Teuchos::rcp(a,false), keep_phase);
    }


    /// Returns \c true if the solver can handle the matrix shape
    virtual bool matrixShapeOK( void ) {
      return flat_solver->matrixShapeOK();
    }


    /// Sets the LHS vector X
    virtual void setX( const Teuchos::RCP<Vector> x ) {
      X = x;
      if (x != Teuchos::null)
        flat_X = Stokhos::create_flat_vector_view(*x, flat_X_map);
      else
        flat_X = Teuchos::null;
      flat_solver->setX(flat_X);
    }


    /// Sets the LHS vector X using a raw pointer
    virtual void setX( Vector* x ) {
      if (x != 0) {
        X = Teuchos::rcp(x, false);
        flat_X = Stokhos::create_flat_vector_view(*x, flat_X_map);
      }
      else {
        X = Teuchos::null;
        flat_X = Teuchos::null;
      }
      flat_solver->setX(flat_X);
    }


    /// Returns the vector that is the LHS of the linear system
    virtual const Teuchos::RCP<Vector> getX( void ) {
      return X;
    }


    /// Returns a raw pointer to the LHS of the linear system
    virtual Vector* getXRaw( void ) {
      return X.get();
    }


    /// Sets the RHS vector B
    virtual void setB( const Teuchos::RCP<const Vector> b ) {
      B = b;
      if (b != Teuchos::null)
        flat_B = Stokhos::create_flat_vector_view(*b, flat_B_map);
      else
        flat_B = Teuchos::null;
      flat_solver->setB(flat_B);
    }


    /// Sets the RHS vector B using a raw pointer
    virtual void setB( const Vector* b ) {
      if (b != 0) {
        B = Teuchos::rcp(b, false);
        flat_B = Stokhos::create_flat_vector_view(*b, flat_B_map);
      }
      else {
        B = Teuchos::null;
        flat_B = Teuchos::null;
      }
      flat_solver->setB(flat_B);
    }


    /// Returns the vector that is the RHS of the linear system
    virtual const Teuchos::RCP<const Vector> getB( void ) {
      return B;
    }


    /// Returns a raw pointer to the RHS of the linear system
    virtual const Vector* getBRaw( void ) {
      return B.get();
    }


    /// Returns a pointer to the Teuchos::Comm communicator with this matrix
    virtual Teuchos::RCP<const Teuchos::Comm<int> > getComm( void ) const {
      return flat_solver->getComm();
    }


    /// Returns a reference to this solver's internal status object
    virtual Status& getStatus() const {
      return flat_solver->getStatus();
    }


    /// Return the name of this solver.
    virtual std::string name( void ) const {
      return flat_solver->name();
    }

    /// @} End Accessor Methods


    /** \name Methods implementing Describable
     * @{
     */

    /// Returns a short description of this Solver
    virtual std::string description( void ) const {
      return flat_solver->description();
    }


    /// Prints the status information about the current solver with some level
    /// of verbosity.
    virtual void describe( Teuchos::FancyOStream &out,
                           const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default ) const {
      flat_solver->describe(out, verbLevel);
    }

    /// @} End Methods implementing Describable


    /** \name Performance and Timing
     * @{
     */

    /// Prints timing information about the current solver.
    virtual void printTiming( Teuchos::FancyOStream &out,
                              const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default ) const{
      flat_solver->printTiming(out, verbLevel);
    }


    /**
     * \brief Extracts timing information from the current solver.
     *
     * Results are placed in the parameter list \c timingParameterList
     *
     * \param timingParameterList Accepts timing information from the
     * current solver
     */
    virtual void getTiming( Teuchos::ParameterList& timingParameterList ) const{
      flat_solver->getTiming(timingParameterList);
    }

    /// @} End Performance and Timing

  protected:

    Teuchos::RCP<const Matrix> A;
    Teuchos::RCP<Vector> X;
    Teuchos::RCP<const Vector> B;
    Teuchos::RCP<const Map> flat_X_map, flat_B_map;
    Teuchos::RCP<const FlatGraph> flat_graph, cijk_graph;
    Teuchos::RCP<const FlatMatrix> flat_A;
    Teuchos::RCP<FlatVector> flat_X;
    Teuchos::RCP<const FlatVector> flat_B;
    Teuchos::RCP<FlatSolver> flat_solver;
    cijk_type cijk;

  };

  // Specialization of create_solver_with_supported_type for
  // Sacado::UQ::PCE where we create PCESolverAdapter wrapping
  // each solver
  template < template <class,class> class ConcreteSolver,
             class ST, class LO, class GO, class NO >
  struct create_solver_with_supported_type<
    ConcreteSolver,
    Tpetra::CrsMatrix<Sacado::UQ::PCE<ST>,LO,GO,NO >,
    Tpetra::MultiVector<Sacado::UQ::PCE<ST>,LO,GO,NO > > {
    typedef Sacado::UQ::PCE<ST> SC;
    typedef Tpetra::CrsMatrix<SC,LO,GO,NO> Matrix;
    typedef Tpetra::MultiVector<SC,LO,GO,NO> Vector;
    static Teuchos::RCP<Solver<Matrix,Vector> >
    apply(Teuchos::RCP<const Matrix> A,
          Teuchos::RCP<Vector>       X,
          Teuchos::RCP<const Vector> B )
    {
      ctassert<
        std::is_same_v<
          typename MatrixTraits<Matrix>::scalar_t,
          typename MultiVecAdapter<Vector>::scalar_t
        >
      > same_scalar_assertion;
      (void)same_scalar_assertion; // This stops the compiler from warning about unused declared variables

      // If our assertion did not fail, then create and return a new solver
      return Teuchos::rcp( new PCESolverAdapter<ST,LO,GO,NO,ConcreteSolver>(A, X, B) );
    }
  };

  // Specialization for solver_supports_scalar for Sacado::UQ::PCE<Storage>
  // value == true if and only if
  // solver_supprts_scalar<ConcreteSolver,Storage::value_type> == true
  template <template <class,class> class ConcreteSolver,
            typename Storage>
  struct solver_supports_scalar<ConcreteSolver, Sacado::UQ::PCE<Storage> > {
    typedef Sacado::UQ::PCE<Storage> Scalar;
    typedef typename Scalar::value_type BaseScalar;
    typedef typename solver_traits<ConcreteSolver>::supported_scalars supported_scalars;
    static const bool value =
      std::conditional_t<std::is_same_v<supported_scalars, Meta::nil_t>,
                         std::true_type,
                         Meta::type_list_contains<supported_scalars,
                                                  BaseScalar> >::value;
  };


} // namespace Amesos2

#endif // AMESOS2_SOLVER_UQ_PCE_HPP

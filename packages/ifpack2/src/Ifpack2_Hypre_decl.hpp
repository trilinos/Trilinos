// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_HYPRE_DECL_HPP
#define IFPACK2_HYPRE_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#if defined(HAVE_IFPACK2_HYPRE) && defined(HAVE_IFPACK2_MPI)

#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Exceptions.hpp"

#include "Ifpack2_Hypre_FunctionParameters.hpp"

#include <map>

// Hypre forward declarations (to avoid downstream header pollution)
struct hypre_IJMatrix_struct;
typedef struct hypre_IJMatrix_struct *HYPRE_IJMatrix;
struct hypre_IJVector_struct;
typedef struct hypre_IJVector_struct *HYPRE_IJVector;
struct hypre_ParCSRMatrix_struct;
typedef struct hypre_ParCSRMatrix_struct* HYPRE_ParCSRMatrix;
struct hypre_ParVector_struct;
typedef struct hypre_ParVector_struct * HYPRE_ParVector;
struct hypre_Solver_struct;
typedef struct hypre_Solver_struct *HYPRE_Solver;
struct hypre_ParVector_struct;
typedef struct hypre_ParVector_struct hypre_ParVector;
//struct hypre_Vector;

namespace Ifpack2 {


#ifndef HYPRE_ENUMS
#define HYPRE_ENUMS
//! This enumerated type defines the allowed solvers and preconditioners in Hypre. Some can be used as both solver and preconditioner.
  enum Hypre_Solver{
    BoomerAMG,
    ParaSails,
    Euclid,
    AMS,
    Hybrid,
    PCG,
    GMRES,
    FlexGMRES,
    LGMRES,
    BiCGSTAB
  };
  
  //! This enumerated type defines the two options for applying inverse, either solve or apply the preconditioner.
  enum Hypre_Chooser{
    Hypre_Is_Solver,
    Hypre_Is_Preconditioner
  };
#endif //HYPRE_ENUMS


class NotImplemented : public Teuchos::ExceptionBase {
public:
  NotImplemented(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
};

//! Ifpack2::Hypre: A class for constructing and using an AMG preconditioner for a given Tpetra::RowMatrix, using the Hypre library by Lawrence Livermore National Laboratory.

/*!
Class Ifpack2::Hypre: A class for using methods of Hypre with Tpetra objects.
*/

template<class MatrixType>
class Hypre:
    virtual public Ifpack2::Preconditioner<typename MatrixType::scalar_type,
                                           typename MatrixType::local_ordinal_type,
                                           typename MatrixType::global_ordinal_type,
                                           typename MatrixType::node_type>,
    virtual public Ifpack2::Details::CanChangeMatrix<Tpetra::RowMatrix<typename MatrixType::scalar_type,
                                                                       typename MatrixType::local_ordinal_type,
                                                                       typename MatrixType::global_ordinal_type,
                                                                       typename MatrixType::node_type> >
{
public:

  //! The type of the entries of the input MatrixType.
  typedef typename MatrixType::scalar_type scalar_type;

  //! The type of local indices in the input MatrixType.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;

  //! The type of global indices in the input MatrixType.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;

  //! The Kokkos::Device specialization used by the input MatrixType.
  typedef typename MatrixType::node_type::device_type device_type;

  //! The Node type used by the input MatrixType.
  typedef typename MatrixType::node_type node_type;

  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  /// \brief The Tpetra::RowMatrix specialization matching MatrixType.
  ///
  /// MatrixType must be a Tpetra::RowMatrix specialization.  This
  /// typedef will always be a Tpetra::RowMatrix specialization.
  typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type> row_matrix_type;

  static_assert (std::is_same<MatrixType, row_matrix_type>::value,
                 "Ifpack2::Hypre: MatrixType must be a Tpetra::RowMatrix "
                 "specialization.  Don't use Tpetra::CrsMatrix here.");


  //! Tpetra::CrsMatrix specialization matching MatrixType
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type> crs_matrix_type;

  //! The Tpetra::Map specialization matching MatrixType.
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

  /// \brief The Tpetra::Vector specialization matching MatrixType.
  ///
  /// If you wish to supply setParameters() a precomputed vector of
  /// diagonal entries of the matrix, use a pointer to an object of
  /// this type.
  typedef Tpetra::Vector<scalar_type, local_ordinal_type,
                         global_ordinal_type, node_type> vector_type;

  /// \brief The Tpetra::MultiVector specialization matching MatrixType.
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
                         global_ordinal_type, node_type> multivector_type;


  explicit Hypre(const Teuchos::RCP<const row_matrix_type>& A) { throw NotImplemented("Ifpack2::Hypre only works when instantiated on <HYPRE_REAL, LocalOrdinal, HYPRE_Int, Node>"); }

  // @}
  // @{ Construction methods
  //! Initialize the preconditioner, does not touch matrix values.
  void initialize() {}

  //! Returns \c true if the preconditioner has been successfully initialized.
  bool isInitialized() const{ return false;}

  //! Compute ILU factors L and U using the specified graph, diagonal perturbation thresholds and relaxation parameters.
  /*! This function computes the ILU(k) factors.
   */
  void compute() {}

  //! If factor is completed, this query returns true, otherwise it returns false.
  bool isComputed() const{ return false;}

  void setParameters(const Teuchos::ParameterList& parameterlist) {}


  //! \name Implementation of Ifpack2::Details::CanChangeMatrix
  //@{

  /// \brief Change the matrix to be preconditioned.
  ///
  /// \param[in] A The new matrix.
  ///
  /// \post <tt>! isInitialized ()</tt>
  /// \post <tt>! isComputed ()</tt>
  ///
  /// Calling this method resets the preconditioner's state.  After
  /// calling this method with a nonnull input, you must first call
  /// initialize() and compute() (in that order) before you may call
  /// apply().
  ///
  /// You may call this method with a null input.  If A is null, then
  /// you may not call initialize() or compute() until you first call
  /// this method again with a nonnull input.  This method invalidates
  /// any previous factorization whether or not A is null, so calling
  /// setMatrix() with a null input is one way to clear the
  /// preconditioner's state (and free any memory that it may be
  /// using).
  ///
  /// The new matrix A need not necessarily have the same Maps or even
  /// the same communicator as the original matrix.
  virtual void
  setMatrix (const Teuchos::RCP<const row_matrix_type>& A) {}
  //@}

  /// \brief Apply the preconditioner to X, returning the result in Y.
  ///
  /// This method actually computes Y = beta*Y + alpha*(M*X), where
  /// M*X represents the result of Chebyshev iteration on X, using the
  /// matrix Op(A).  Op(A) is either A itself, its transpose
  /// \f$A^T\f$, or its Hermitian transpose \f$A^H\f$, depending on
  /// the <tt>mode</tt> argument.  Since this class currently requires
  /// A to be real and symmetric positive definite, it should always
  /// be the case that \f$A = A^T = A^H\f$, but we will still respect
  /// the <tt>mode</tt> argument.
  ///
  /// \warning If you did not set the "chebyshev: zero starting
  ///   solution" parameter to true, then this method will use X as
  ///   the starting guess for Chebyshev iteration.  If you did not
  ///   initialize X before calling this method, then the resulting
  ///   solution will be undefined, since it will be computed using
  ///   uninitialized data.
  ///
  /// \param[in] X  A (multi)vector to which to apply the preconditioner.
  /// \param[in,out] Y A (multi)vector containing the result of
  ///   applying the preconditioner to X.
  /// \param[in] mode  If <tt>Teuchos::NO_TRANS</tt>, apply the matrix
  ///   A.  If <tt>mode</tt> is <tt>Teuchos::NO_TRANS</tt>, apply its
  ///   transpose \f$A^T\f$.  If <tt>Teuchos::CONJ_TRANS</tt>, apply
  ///   its Hermitian transpose \f$A^H\f$.
  /// \param[in] alpha  Scaling factor for the result of Chebyshev
  ///   iteration.  The default is 1.
  /// \param[in] beta  Scaling factor for Y.  The default is 0.
  void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const {}

  //! The Tpetra::Map representing the domain of this operator.
  Teuchos::RCP<const map_type> getDomainMap() const { return Teuchos::null; }

  //! The Tpetra::Map representing the range of this operator.
  Teuchos::RCP<const map_type> getRangeMap() const { return Teuchos::null; }

  //! Whether it's possible to apply the transpose of this operator.
  bool hasTransposeApply() const { return false; }

  /// \brief Compute Y = Op(A)*X, where Op(A) is either A, \f$A^T\f$, or \f$A^H\f$.
  ///
  /// \param[in] X  Input (multi)vector of sparse matrix-vector
  ///   multiply.  If mode == Teuchos::NO_TRANS, X must be in the
  ///   domain Map of the matrix A.  Otherwise, X must be in the range
  ///   Map of A.
  /// \param[out] Y  Output (multi)vector of sparse matrix-vector
  ///   multiply.  If mode == Teuchos::NO_TRANS, Y must be in the
  ///   range Map of the matrix A.  Otherwise, Y must be in the domain
  ///   Map of A.
  /// \param[in] mode  Whether to apply the matrix A, its transpose
  ///   \f$A^T\f$, or its conjugate transpose \f$A^H\f$.  This method
  ///   applies A if <tt>mode</tt> is <tt>Teuchos::NO_TRANS</tt>,
  ///   \f$A^T\f$ if <tt>mode</tt> is <tt>Teuchos::TRANS</tt>, and
  ///   \f$A^H\f$ (the Hermitian transpose) if <tt>mode</tt> is
  ///   <tt>Teuchos::CONJ_TRANS</tt>.
  ///
  /// Since this class currently requires A to be real and symmetric
  /// positive definite, setting <tt>mode</tt> should not affect the
  /// result.
  void
  applyMat (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
            Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
            Teuchos::ETransp mode = Teuchos::NO_TRANS) const {}

  Teuchos::RCP<const row_matrix_type> getMatrix() const { return Teuchos::null; }

  //@}

  //! The number of calls to initialize().
  int getNumInitialize() const { return 0; }

  //! The number of calls to compute().
  int getNumCompute() const { return 0; }

  //! The number of calls to apply().
  int getNumApply() const { return 0; }

  //! The time (in seconds) spent in initialize().
  double getInitializeTime() const { return 0.0; }

  //! The time (in seconds) spent in compute().
  double getComputeTime() const { return 0.0; }

  //! The time (in seconds) spent in apply().
  double getApplyTime() const { return 0.0; }


};

template<class LocalOrdinal, class Node>
class Hypre<Tpetra::RowMatrix<HYPRE_Real, LocalOrdinal, HYPRE_Int, Node> >:
    virtual public Ifpack2::Preconditioner<HYPRE_Real,
                                           LocalOrdinal,
                                           HYPRE_Int,
                                           Node>,
    virtual public Ifpack2::Details::CanChangeMatrix<Tpetra::RowMatrix<HYPRE_Real,
                                                                       LocalOrdinal,
                                                                       HYPRE_Int,
                                                                       Node> >
{
public:
  //! \name Typedefs
  //@{

  //! The template parameter of this class.
  typedef Tpetra::RowMatrix<HYPRE_Real, LocalOrdinal, HYPRE_Int, Node> MatrixType;

  //! The type of the entries of the input MatrixType.
  typedef typename MatrixType::scalar_type scalar_type;

  //! The type of local indices in the input MatrixType.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;

  //! The type of global indices in the input MatrixType.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;

  //! The Kokkos::Device specialization used by the input MatrixType.
  typedef typename MatrixType::node_type::device_type device_type;

  //! The Node type used by the input MatrixType.
  typedef typename MatrixType::node_type node_type;

  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  /// \brief The Tpetra::RowMatrix specialization matching MatrixType.
  ///
  /// MatrixType must be a Tpetra::RowMatrix specialization.  This
  /// typedef will always be a Tpetra::RowMatrix specialization.
  typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type> row_matrix_type;

  static_assert (std::is_same<MatrixType, row_matrix_type>::value,
                 "Ifpack2::Hypre: MatrixType must be a Tpetra::RowMatrix "
                 "specialization.  Don't use Tpetra::CrsMatrix here.");


  //! Tpetra::CrsMatrix specialization matching MatrixType
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type> crs_matrix_type;

  //! The Tpetra::Map specialization matching MatrixType.
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

  /// \brief The Tpetra::Vector specialization matching MatrixType.
  ///
  /// If you wish to supply setParameters() a precomputed vector of
  /// diagonal entries of the matrix, use a pointer to an object of
  /// this type.
  typedef Tpetra::Vector<scalar_type, local_ordinal_type,
                         global_ordinal_type, node_type> vector_type;

  /// \brief The Tpetra::MultiVector specialization matching MatrixType.
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
                         global_ordinal_type, node_type> multivector_type;
  
  // Hypre Specs
  // This will need to be either int or long long depending on how Hypre was built
  //    typedef global_ordinal_type global_ordinal_type;

  typedef global_ordinal_type (*HYPRE_PtrToParSolverFcn)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);

  //@}
  // \name Constructors and destructors
  //@{

  /// \brief Constructor.
  ///
  /// \param[in] A The sparse matrix to which to apply Hypre
  ///   iteration.  The matrix A must be square, and its domain Map
  ///   and range Map must be the same.  The latter means that the
  ///   vectors x and y in the sparse matrix-vector product y = A*x
  ///   must both have the same distribution over process(es).
  ///
  /// We do <i>not</i> require that the row Map and the range Map of A
  /// be the same.  However, set-up will take less time if they are
  /// identical (in terms of pointer equality).  This is because we
  /// have to extract the diagonal entries of A as a row Map vector:
  /// if the row and range Maps are not identical, we have to
  /// redistribute the vector from the row Map to the range Map.
  ///
  /// The constructor will only check the requirements on the various
  /// Maps of A if the CMake configuration option
  /// <tt>Teuchos_ENABLE_DEBUG</tt> was set to <tt>ON</tt> before
  /// building Trilinos.  The checks require \f$O(1)\f$ global
  /// reductions over all processes in A's communicator, so we prefer
  /// to avoid them if we can.
  explicit Hypre(const Teuchos::RCP<const row_matrix_type>& A);

  //! Destructor
  ~Hypre();

  // @}
  // @{ Construction methods
  //! Initialize the preconditioner, does not touch matrix values.
  void initialize();

  //! Returns \c true if the preconditioner has been successfully initialized.
  bool isInitialized() const{ return(IsInitialized_);}

  //! Compute ILU factors L and U using the specified graph, diagonal perturbation thresholds and relaxation parameters.
  /*! This function computes the ILU(k) factors.
   */
  void compute();

  //! If factor is completed, this query returns true, otherwise it returns false.
  bool isComputed() const{ return(IsComputed_);}
 
  //! Set parameters using a Teuchos::ParameterList object.
  /* This method is only available if the Teuchos package is enabled.
     This method recognizes six parameter names: Solver,
     Preconditioner, SolveOrPrecondition, SetPreconditioner, NumFunctions and Functions. These names are
     case sensitive. Solver requires an enumerated parameter of type Hypre_Solver. Preconditioner is similar
     except requires the type be a preconditioner. The options are listed below:
                       Solvers                            Preconditioners
                       BoomerAMG                          BoomerAMG
                       AMS                                ParaSails
                       Hybrid                             AMS
                       PCG (Default)                      Euclid (Default)
                       GMRES
                       FlexGMRES
                       LGMRES
                       BiCGSTAB
     SolveOrPrecondition takes enumerated type Hypre_Chooser, Solver will solve the system, Preconditioner will apply the preconditioner.
     SetPreconditioner takes a boolean, true means the solver will use the preconditioner.
     NumFunctions takes an int that describes how many parameters will be passed into Functions. (This needs to be correct.)
     Functions takes an array of Ref Counted Pointers to an object called FunctionParameter. This class is implemented in Ifpack_Hypre.h.
     The object takes whether it is Solver or Preconditioner that we are setting a parameter for.
     The function in Hypre that sets the parameter, and the parameters for that function. An example is below:

     RCP<FunctionParameter> functs[2];
     functs[0] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetMaxIter, 1000)); // max iterations
     functs[1] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetTol, 1e-7)); // conv. tolerance
     list.set("NumFunctions", 2);
     list.set<RCP<FunctionParameter>*>("Functions", functs);
     NOTE: SetParameters() must be called to use ApplyInverse(), the solvers will not be created otherwise. An empty list is acceptable to use defaults.
  */
  void setParameters(const Teuchos::ParameterList& parameterlist);

    //! Set a parameter that takes a single int.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set maximum iterations would be &HYPRE_BoomerAMGSetMaxIter
    \param parameter (In) -The integer parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre_Chooser chooser, HYPRE_Int (*pt2Func)(HYPRE_Solver, HYPRE_Int), HYPRE_Int parameter);

    //! Set a parameter that takes a single double.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set tolerance would be &HYPRE_BoomerAMGSetTol
    \param parameter (In) -The double parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre_Chooser chooser, HYPRE_Int (*pt2Func)(HYPRE_Solver, HYPRE_Real), HYPRE_Real parameter);

    //! Set a parameter that takes a double then an int.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set relaxation weight for a given level would be &HYPRE_BoomerAMGSetLevelRelaxWt
    \param parameter1 (In) -The double parameter being set.
    \param parameter2 (In) - The integer parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre_Chooser chooser, HYPRE_Int (*pt2Func)(HYPRE_Solver, HYPRE_Real, HYPRE_Int), HYPRE_Real parameter1, HYPRE_Int parameter2);

    //! Set a parameter that takes an int then a double.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set relaxation weight for a given level would be &HYPRE_BoomerAMGSetLevelRelaxWt
    \param parameter1 (In) -The int parameter being set.
    \param parameter2 (In) - The double parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre_Chooser chooser, HYPRE_Int (*pt2Func)(HYPRE_Solver, HYPRE_Int, HYPRE_Real), HYPRE_Int parameter1, HYPRE_Real parameter2);

    //! Set a parameter that takes two int parameters.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set relaxation type for a given level would be &HYPRE_BoomerAMGSetCycleRelaxType
    \param parameter1 (In) -The first integer parameter being set.
    \param parameter2 (In) - The second integer parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre_Chooser chooser, HYPRE_Int (*pt2Func)(HYPRE_Solver, HYPRE_Int, HYPRE_Int), HYPRE_Int parameter1, HYPRE_Int parameter2);

    //! Set a parameter that takes a double*.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set relaxation weight would be &HYPRE_BoomerAMGSetRelaxWeight
    \param parameter (In) -The double* parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre_Chooser chooser, HYPRE_Int (*pt2Func)(HYPRE_Solver, HYPRE_Real*), HYPRE_Real* parameter);

    //! Set a parameter that takes an int*.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set grid relax type would be &HYPRE_BoomerAMGSetGridRelaxType
    \param parameter (In) -The int* parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre_Chooser chooser, HYPRE_Int (*pt2Func)(HYPRE_Solver, HYPRE_Int*), HYPRE_Int* parameter);

    //! Set a parameter that takes an int**.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set the order in which the points are relaxed would be
      &HYPRE_BoomerAMGSetGridRelaxPoints used primarily for AIR AMG.
    \param parameter (In) -The int** parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre_Chooser chooser, HYPRE_Int (*pt2Func)(HYPRE_Solver, HYPRE_Int**), HYPRE_Int** parameter);

    //! Sets the solver that is used by the Solve() and ApplyInverse() methods. Until this is called, the default solver is PCG.
    /*!
    \param chooser (In) - A Hypre_Chooser enumerated type. If Solver, then we are selecting which solver, if Preconditioner, we are choosing which preconditioner to use.
    \param Solver (In) -A Hypre_Solver enumerated type to select the solver or preconditioner. Options for solver are:
    BoomerAMG, AMS, Hybrid, PCG, GMRES, FlexGMRES, LGMRES, and BiCGSTAB. See Hypre Ref Manual for more info on the solvers.
    Options for Preconditioner are: BoomerAMG, ParaSails, Euclid, and AMS.

    \return Integer error code, set to 0 if successful.
  */

    int SetParameter(Hypre_Chooser chooser, Hypre_Solver Solver);

    //! Sets the solver to use the selected preconditioner.
    /*!
    \param UsePreconditioner (In) -A boolean, true use preconditioner, false do not use the supplied preconditioner with the solver.
    The solver and preconditioner must have been selected and the solver must be one of the following solvers:
      Hybrid, PCG, GMRES, FlexGMRES, LGMRES, BiCGSTAB.

    \return Integer error code, set to 0 if successful.
  */

    int SetParameter(bool UsePreconditioner){ UsePreconditioner_ = UsePreconditioner; return 0;}

    //! Choose to solve the problem or apply the preconditioner.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type, either Solver or Preconditioner.
    The chosen type must have been selected before this method is called.

    \return Integer error code, set to 0 if successful.
  */
    int SetParameter(Hypre_Chooser chooser) { SolveOrPrec_ = chooser; return 0;}

  //! Set coordinates
  int SetCoordinates(Teuchos::RCP<multivector_type> coords);

  //! Set discrete gradient
  int SetDiscreteGradient(Teuchos::RCP<const crs_matrix_type> G);

  //! Call all the function pointers stored in this object.
    int CallFunctions() const;

//@}
  //! \name Implementation of Ifpack2::Details::CanChangeMatrix
  //@{

  /// \brief Change the matrix to be preconditioned.
  ///
  /// \param[in] A The new matrix.
  ///
  /// \post <tt>! isInitialized ()</tt>
  /// \post <tt>! isComputed ()</tt>
  ///
  /// Calling this method resets the preconditioner's state.  After
  /// calling this method with a nonnull input, you must first call
  /// initialize() and compute() (in that order) before you may call
  /// apply().
  ///
  /// You may call this method with a null input.  If A is null, then
  /// you may not call initialize() or compute() until you first call
  /// this method again with a nonnull input.  This method invalidates
  /// any previous factorization whether or not A is null, so calling
  /// setMatrix() with a null input is one way to clear the
  /// preconditioner's state (and free any memory that it may be
  /// using).
  ///
  /// The new matrix A need not necessarily have the same Maps or even
  /// the same communicator as the original matrix.
  virtual void
  setMatrix (const Teuchos::RCP<const row_matrix_type>& A);
  //@}

  /// \brief Apply the preconditioner to X, returning the result in Y.
  ///
  /// This method actually computes Y = beta*Y + alpha*(M*X), where
  /// M*X represents the result of Chebyshev iteration on X, using the
  /// matrix Op(A).  Op(A) is either A itself, its transpose
  /// \f$A^T\f$, or its Hermitian transpose \f$A^H\f$, depending on
  /// the <tt>mode</tt> argument.  Since this class currently requires
  /// A to be real and symmetric positive definite, it should always
  /// be the case that \f$A = A^T = A^H\f$, but we will still respect
  /// the <tt>mode</tt> argument.
  ///
  /// \warning If you did not set the "chebyshev: zero starting
  ///   solution" parameter to true, then this method will use X as
  ///   the starting guess for Chebyshev iteration.  If you did not
  ///   initialize X before calling this method, then the resulting
  ///   solution will be undefined, since it will be computed using
  ///   uninitialized data.
  ///
  /// \param[in] X  A (multi)vector to which to apply the preconditioner.
  /// \param[in,out] Y A (multi)vector containing the result of
  ///   applying the preconditioner to X.
  /// \param[in] mode  If <tt>Teuchos::NO_TRANS</tt>, apply the matrix
  ///   A.  If <tt>mode</tt> is <tt>Teuchos::NO_TRANS</tt>, apply its
  ///   transpose \f$A^T\f$.  If <tt>Teuchos::CONJ_TRANS</tt>, apply
  ///   its Hermitian transpose \f$A^H\f$.
  /// \param[in] alpha  Scaling factor for the result of Chebyshev
  ///   iteration.  The default is 1.
  /// \param[in] beta  Scaling factor for Y.  The default is 0.
  void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! The Tpetra::Map representing the domain of this operator.
  Teuchos::RCP<const map_type> getDomainMap() const;

  //! The Tpetra::Map representing the range of this operator.
  Teuchos::RCP<const map_type> getRangeMap() const;

  //! Whether it's possible to apply the transpose of this operator.
  bool hasTransposeApply() const;

  /// \brief Compute Y = Op(A)*X, where Op(A) is either A, \f$A^T\f$, or \f$A^H\f$.
  ///
  /// \param[in] X  Input (multi)vector of sparse matrix-vector
  ///   multiply.  If mode == Teuchos::NO_TRANS, X must be in the
  ///   domain Map of the matrix A.  Otherwise, X must be in the range
  ///   Map of A.
  /// \param[out] Y  Output (multi)vector of sparse matrix-vector
  ///   multiply.  If mode == Teuchos::NO_TRANS, Y must be in the
  ///   range Map of the matrix A.  Otherwise, Y must be in the domain
  ///   Map of A.
  /// \param[in] mode  Whether to apply the matrix A, its transpose
  ///   \f$A^T\f$, or its conjugate transpose \f$A^H\f$.  This method
  ///   applies A if <tt>mode</tt> is <tt>Teuchos::NO_TRANS</tt>,
  ///   \f$A^T\f$ if <tt>mode</tt> is <tt>Teuchos::TRANS</tt>, and
  ///   \f$A^H\f$ (the Hermitian transpose) if <tt>mode</tt> is
  ///   <tt>Teuchos::CONJ_TRANS</tt>.
  ///
  /// Since this class currently requires A to be real and symmetric
  /// positive definite, setting <tt>mode</tt> should not affect the
  /// result.
  void
  applyMat (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
            Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
            Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //@}
  //! \name Attribute accessor methods
  //@{

  //! The communicator over which the matrix is distributed.
  Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

  //! The matrix for which this is a preconditioner.
  Teuchos::RCP<const row_matrix_type> getMatrix() const;

  /// \brief Attempt to return the matrix A as a Tpetra::CrsMatrix.
  ///
  /// This class does not require that A be a Tpetra::CrsMatrix.
  /// If it is NOT, this method will return Teuchos::null.
  Teuchos::RCP<const Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> >
  getCrsMatrix() const;

  //! The number of calls to initialize().
  int getNumInitialize() const;

  //! The number of calls to compute().
  int getNumCompute() const;

  //! The number of calls to apply().
  int getNumApply() const;

  //! The time (in seconds) spent in initialize().
  double getInitializeTime() const;

  //! The time (in seconds) spent in compute().
  double getComputeTime() const;

  //! The time (in seconds) spent in apply().
  double getApplyTime() const;

  //@}
  //! @name Implementation of Teuchos::Describable
  //@{

  //! A simple one-line description of this object.
  std::string description() const;

  //! Print the object with some verbosity level to a Teuchos::FancyOStream.
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  //@}

private:
  // @{ Private methods

  //! Abbreviation for the Teuchos::ScalarTraits specialization for scalar_type.
  typedef Teuchos::ScalarTraits<typename MatrixType::scalar_type> STS;

  //! Abbreviation for the Tpetra::MultiVector specialization used in methods like apply().
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> MV;

  //! Copy constructor (use is syntactically forbidden)
  Hypre(const Hypre<MatrixType>&);

  //! Assignment operator (use is syntactically forbidded)
  Hypre<MatrixType>& operator= (const Hypre<MatrixType>&);

  //! Sets the solver type to be the passed in solver type.
  int SetSolverType(Hypre_Solver solver);

  //! Sets the preconditioner type to be the passed in type.
  int SetPrecondType(Hypre_Solver precond);

  //! Create the solver.
  int CreateSolver();

  //! Create the Preconditioner.
  int CreatePrecond();

  //! Copies matrix data from Tpetra matrix to Hypre matrix.
  int CopyTpetraToHypre();

  //! Add a function to be called in Compute()
  int AddFunToList(Teuchos::RCP<FunctionParameter> NewFun);

  //! Create a BoomerAMG solver.
  HYPRE_Int Hypre_BoomerAMGCreate(MPI_Comm comm, HYPRE_Solver *solver);

  //! Create a ParaSails solver.
  HYPRE_Int Hypre_ParaSailsCreate(MPI_Comm comm, HYPRE_Solver *solver);

  //! Create a Euclid solver.
  HYPRE_Int Hypre_EuclidCreate(MPI_Comm comm, HYPRE_Solver *solver);

  //! Create an AMS solver.
  HYPRE_Int Hypre_AMSCreate(MPI_Comm comm, HYPRE_Solver *solver);

  //! Create a Hybrid solver.
  HYPRE_Int Hypre_ParCSRHybridCreate(MPI_Comm comm, HYPRE_Solver *solver);

  //! Create a PCG solver.
  HYPRE_Int Hypre_ParCSRPCGCreate(MPI_Comm comm, HYPRE_Solver *solver);

  //! Create a GMRES solver.
  HYPRE_Int Hypre_ParCSRGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver);

  //! Create a FlexGMRES solver.
  HYPRE_Int Hypre_ParCSRFlexGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver);

  //! Create a LGMRES solver.
  HYPRE_Int Hypre_ParCSRLGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver);

  //! Create a BiCGSTAB solver.
  HYPRE_Int Hypre_ParCSRBiCGSTABCreate(MPI_Comm comm, HYPRE_Solver *solver);
  
  //! Map generation function
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, HYPRE_Int, Node> >
  MakeContiguousColumnMap(Teuchos::RCP<const crs_matrix_type> &Matrix) const;

  //! Destroy
  void Destroy();

  // @}
  // @{ Internal data
  //! Pointer to the Tpetra::RowMatrix;
  Teuchos::RCP<const MatrixType> A_;
  //! This objects copy of the ParameterList
  Teuchos::ParameterList List_;

  //! If \c true, initialize() has completed successfully.
  bool IsInitialized_;
  //! If \c true, compute() has completed successfully.
  bool IsComputed_;
  //! The total number of successful calls to initialize().
  int NumInitialize_;
  //! The total number of successful calls to compute().
  int NumCompute_;
  /// \brief The total number of successful calls to apply().
  ///
  /// This is "mutable" because apply() is a const method; apply() is
  /// const because it is declared this way in Tpetra::Operator.
  mutable int NumApply_;
  //! The total time in seconds over all calls to initialize().
  double InitializeTime_;
  //! The total time in seconds over all calls to compute().
  double ComputeTime_;
  /// \brief The total time in seconds over all calls to apply().
  ///
  /// This is "mutable" because apply() is a const method; apply() is
  /// const because it is declared this way in Tpetra::Operator.
  mutable double ApplyTime_;
  //! The total number of floating-point operations over all calls to compute().
  double ComputeFlops_;
  /// \brief The total number of floating-point operations over all calls to apply().
  ///
  /// This is "mutable" because apply() is a const method; apply() is
  /// const because it is declared this way in Tpetra::Operator.
  mutable double ApplyFlops_;  
 

  //! The Hypre matrix created in initialize()
  mutable HYPRE_IJMatrix HypreA_;
  //! Pointer to the CSR (same matrix)
  mutable HYPRE_ParCSRMatrix ParMatrix_;

  //! Tpetra copy of discrete gradient
  Teuchos::RCP<const crs_matrix_type> G_;
  //! The Hypre matrix created in SetDiscreteGradient)
  mutable HYPRE_IJMatrix HypreG_;
  //! Pointer to the CSR (same matrix)
  mutable HYPRE_ParCSRMatrix ParMatrixG_;

  //! The Hypre Vector for input
  mutable HYPRE_IJVector XHypre_;
  //! The Hypre Vector for output
  mutable HYPRE_IJVector YHypre_;
  mutable HYPRE_ParVector ParX_;
  mutable HYPRE_ParVector ParY_;
  mutable Teuchos::RCP<hypre_ParVector> XVec_;
  mutable Teuchos::RCP<hypre_ParVector> YVec_;

  Teuchos::RCP<multivector_type> Coords_;
  mutable HYPRE_IJVector xHypre_;
  mutable HYPRE_IJVector yHypre_;
  mutable HYPRE_IJVector zHypre_;
  mutable HYPRE_ParVector xPar_;
  mutable HYPRE_ParVector yPar_;
  mutable HYPRE_ParVector zPar_;

  //! The Hypre Solver if doing a solve
  mutable HYPRE_Solver Solver_;
  //! The Hypre Solver if applying preconditioner
  mutable HYPRE_Solver Preconditioner_;
  //  The following are pointers to functions to use the solver and preconditioner.
  HYPRE_Int (Hypre::*SolverCreatePtr_)(MPI_Comm, HYPRE_Solver*);
  HYPRE_Int (*SolverDestroyPtr_)(HYPRE_Solver);
  HYPRE_Int (*SolverSetupPtr_)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*SolverSolvePtr_)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*SolverPrecondPtr_)(HYPRE_Solver, HYPRE_PtrToParSolverFcn, HYPRE_PtrToParSolverFcn, HYPRE_Solver);
  HYPRE_Int (Hypre::*PrecondCreatePtr_)(MPI_Comm, HYPRE_Solver*);
  HYPRE_Int (*PrecondDestroyPtr_)(HYPRE_Solver);
  HYPRE_Int (*PrecondSetupPtr_)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*PrecondSolvePtr_)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);

  bool IsSolverCreated_;
  bool IsPrecondCreated_;
  //! Is the system to be solved or apply preconditioner
  Hypre_Chooser SolveOrPrec_;
  //! These are linear maps that meet the needs of Hypre
  Teuchos::RCP<const map_type> GloballyContiguousRowMap_;
  Teuchos::RCP<const map_type> GloballyContiguousColMap_;
  Teuchos::RCP<const map_type> GloballyContiguousNodeRowMap_;
  Teuchos::RCP<const map_type> GloballyContiguousNodeColMap_;
  //! Counter of the number of parameters set
  int NumFunsToCall_;
  //! Which solver was chosen
  Hypre_Solver SolverType_;
  //! Which preconditioner was chosen
  Hypre_Solver PrecondType_;
  //! Should the preconditioner be used in the solver
  bool UsePreconditioner_;
  //! This contains a list of function pointers that will be called in compute
  std::vector<Teuchos::RCP<FunctionParameter> > FunsToCall_;
  //! Should information be dumped to files
  bool Dump_;
  //! Dummy vector for caching
  mutable Teuchos::ArrayRCP<scalar_type> VectorCache_;

};


}//end Ifpack2 namespace

#endif // HAVE_IFPACK2_HYPRE && HAVE_IFPACK2_MPI
#endif /* IFPACK2_HYPRE_DECL_HPP */

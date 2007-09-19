
// File: index.xml

// File: classAnasazi_1_1AnasaziError.xml
%feature("docstring") Anasazi::AnasaziError "

An exception class parent to all Anasazi exceptions.

C++ includes: AnasaziTypes.hpp ";

%feature("docstring")  Anasazi::AnasaziError::AnasaziError "Anasazi::AnasaziError::AnasaziError(const std::string &what_arg) ";


// File: classAnasazi_1_1BasicEigenproblem.xml
%feature("docstring") Anasazi::BasicEigenproblem "

This provides a basic implementation for defining standard or
generalized eigenvalue problems.

C++ includes: AnasaziBasicEigenproblem.hpp ";


// File: classAnasazi_1_1BasicOrthoManager.xml
%feature("docstring") Anasazi::BasicOrthoManager "

An implementation of the Anasazi::MatOrthoManager that performs
orthogonalization using (potentially) multiple steps of classical
Gram-Schmidt.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziBasicOrthoManager.hpp ";


// File: classAnasazi_1_1BasicOutputManager.xml
%feature("docstring") Anasazi::BasicOutputManager "

Anasazi's basic output manager for sending information of select
verbosity levels to the appropriate output stream.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziBasicOutputManager.hpp ";


// File: classAnasazi_1_1BasicSort.xml
%feature("docstring") Anasazi::BasicSort "

An implementation of the Anasazi::SortManager that performs a
collection of common sorting techniques.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziBasicSort.hpp ";

%feature("docstring")  Anasazi::BasicSort::BasicSort "Anasazi::BasicSort< ScalarType, MV, OP >::BasicSort(const std::string
which=\"LM\")

Constructor.

Parameters:
-----------

which:  [in] The eigenvalues of interest for this eigenproblem. \"LM\"
- Largest Magnitude [ default ]

\"SM\" - Smallest Magnitude

\"LR\" - Largest Real

\"SR\" - Smallest Real

\"LI\" - Largest Imaginary

\"SI\" - Smallest Imaginary ";

%feature("docstring")  Anasazi::BasicSort::~BasicSort "virtual
Anasazi::BasicSort< ScalarType, MV, OP >::~BasicSort()

Destructor. ";

%feature("docstring")  Anasazi::BasicSort::setSortType "void
Anasazi::BasicSort< ScalarType, MV, OP >::setSortType(const
std::string which)

Set sort type.

Parameters:
-----------

which:  [in] The eigenvalues of interest for this eigenproblem. \"LM\"
- Largest Magnitude [ default ]

\"SM\" - Smallest Magnitude

\"LR\" - Largest Real

\"SR\" - Smallest Real

\"LI\" - Largest Imaginary

\"SI\" - Smallest Imaginary ";

%feature("docstring")  Anasazi::BasicSort::sort "void
Anasazi::BasicSort< ScalarType, MV, OP >::sort(Eigensolver<
ScalarType, MV, OP > *solver, const int n, std::vector< typename
Teuchos::ScalarTraits< ScalarType >::magnitudeType > &evals,
std::vector< int > *perm=0) const

Sort the vector of eigenvalues, optionally returning the permutation
vector.

Parameters:
-----------

solver:  [in] Eigensolver that is calling the sorting routine

n:  [in] Number of values in evals to be sorted.

evals:  [in/out] Vector of length n containing the eigenvalues to be
sorted

perm:  [out] Vector of length n to store the permutation index
(optional) ";

%feature("docstring")  Anasazi::BasicSort::sort "void
Anasazi::BasicSort< ScalarType, MV, OP >::sort(Eigensolver<
ScalarType, MV, OP > *solver, const int n, std::vector< typename
Teuchos::ScalarTraits< ScalarType >::magnitudeType > &r_evals,
std::vector< typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType > &i_evals, std::vector< int > *perm=0) const

Sort the vectors of eigenpairs, optionally returning the permutation
vector.

This routine takes two vectors, one for each part of a complex
eigenvalue. This is helpful for solving real, non-symmetric eigenvalue
problems.

Parameters:
-----------

solver:  [in] Eigensolver that is calling the sorting routine

n:  [in] Number of values in r_evals,i_evals to be sorted.

r_evals:  [in/out] Vector of length n containing the real part of the
eigenvalues to be sorted

i_evals:  [in/out] Vector of length n containing the imaginary part of
the eigenvalues to be sorted

perm:  [out] Vector of length n to store the permutation index
(optional) ";


// File: classAnasazi_1_1BlockDavidson.xml
%feature("docstring") Anasazi::BlockDavidson "

This class implements a Block Davidson iteration, a preconditioned
iteration for solving linear Hermitian eigenproblems.

This method is described in A Comparison of Eigensolvers for Large-
scale 3D Modal Analysis Using AMG-Preconditioned Iterative Methods, P.
Arbenz, U. L. Hetmaniuk, R. B. Lehoucq, R. S. Tuminaro, Internat. J.
for Numer. Methods Engrg., 64, pp. 204-236 (2005)

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist

C++ includes: AnasaziBlockDavidson.hpp ";


// File: structAnasazi_1_1BlockDavidson_1_1CheckList.xml


// File: classAnasazi_1_1BlockDavidsonInitFailure.xml
%feature("docstring") Anasazi::BlockDavidsonInitFailure "

BlockDavidsonInitFailure is thrown when the BlockDavidson solver is
unable to generate an initial iterate in the
BlockDavidson::initialize() routine.

This exception is thrown from the BlockDavidson::initialize() method,
which is called by the user or from the BlockDavidson::iterate()
method if isInitialized() == false.

In the case that this exception is thrown,
BlockDavidson::isInitialized() will be false and the user will need to
provide a new initial iterate to the solver.

C++ includes: AnasaziBlockDavidson.hpp ";

%feature("docstring")
Anasazi::BlockDavidsonInitFailure::BlockDavidsonInitFailure "Anasazi::BlockDavidsonInitFailure::BlockDavidsonInitFailure(const
std::string &what_arg) ";


// File: classAnasazi_1_1BlockDavidsonOrthoFailure.xml
%feature("docstring") Anasazi::BlockDavidsonOrthoFailure "

BlockDavidsonOrthoFailure is thrown when the orthogonalization manager
is unable to orthogonalize the preconditioned residual against (a.k.a.
H) the current basis (a.k.a. V).

This exception is thrown from the BlockDavidson::iterate() method.

C++ includes: AnasaziBlockDavidson.hpp ";

%feature("docstring")
Anasazi::BlockDavidsonOrthoFailure::BlockDavidsonOrthoFailure "Anasazi::BlockDavidsonOrthoFailure::BlockDavidsonOrthoFailure(const
std::string &what_arg) ";


// File: classAnasazi_1_1BlockDavidsonSolMgr.xml
%feature("docstring") Anasazi::BlockDavidsonSolMgr "

The BlockDavidsonSolMgr provides a powerful solver manager over the
BlockDavidson eigensolver.

This solver manager implements a hard-locking mechanism, whereby
eigenpairs designated to be locked are moved from the eigensolver and
placed in auxilliary storage. The eigensolver is then restarted and
continues to iterate, orthogonal to the locked eigenvectors.

The solver manager provides to the solver a StatusTestCombo object
constructed as follows: combo = globaltest OR lockingtest OR debugtest
where  globaltest terminates computation when global convergence has
been detected.  It is encapsulated in a StatusTestWithOrdering object,
to ensure that computation is terminated only after the most
significant eigenvalues/eigenvectors have met the convergence
criteria.  If not specified via setGlobalStatusTest(), globaltest is a
StatusTestResNorm object which tests the M-norms of the direct
residuals relative to the Ritz values.

lockingtest halts BlockDavidson::iterate() in order to deflate
converged eigenpairs for locking.  It will query the underlying
BlockDavidson eigensolver to determine when eigenvectors should be
locked.  If not specified via setLockingStatusTest(), lockingtest is a
StatusTestResNorm object.

debugtest allows a user to specify additional monitoring of the
iteration, encapsulated in a StatusTest object  If not specified via
setDebugStatusTest(), debugtest is ignored.  In most cases, it should
return Failed; if it returns Passed, solve() will throw an
AnasaziError exception.

Additionally, the solver manager will terminate solve() after a
specified number of restarts.

Much of this behavior is controlled via parameters and options passed
to the solver manager. For more information, see
BlockDavidsonSolMgr().

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist

C++ includes: AnasaziBlockDavidsonSolMgr.hpp ";


// File: structAnasazi_1_1BlockDavidsonState.xml
%feature("docstring") Anasazi::BlockDavidsonState "

Structure to contain pointers to BlockDavidson state variables.

This struct is utilized by BlockDavidson::initialize() and
BlockDavidson::getState().

C++ includes: AnasaziBlockDavidson.hpp ";

%feature("docstring")  Anasazi::BlockDavidsonState::BlockDavidsonState
"Anasazi::BlockDavidsonState< ScalarType, MV >::BlockDavidsonState()
";


// File: classAnasazi_1_1BlockKrylovSchur.xml
%feature("docstring") Anasazi::BlockKrylovSchur "

This class implements the block Krylov-Schur iteration, for solving
linear eigenvalue problems.

This method is a block version of the iteration presented by G.W.
Stewart in \"A Krylov-Schur Algorithm for Large Eigenproblems\", SIAM
J. Matrix Anal. Appl., Vol 23(2001), No. 3, pp. 601-614.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist

C++ includes: AnasaziBlockKrylovSchur.hpp ";


// File: structAnasazi_1_1BlockKrylovSchur_1_1CheckList.xml


// File: classAnasazi_1_1BlockKrylovSchurInitFailure.xml
%feature("docstring") Anasazi::BlockKrylovSchurInitFailure "

BlockKrylovSchurInitFailure is thrown when the BlockKrylovSchur solver
is unable to generate an initial iterate in the
BlockKrylovSchur::initialize() routine.

This exception is thrown from the BlockKrylovSchur::initialize()
method, which is called by the user or from the
BlockKrylovSchur::iterate() method if isInitialized() == false.

In the case that this exception is thrown,
BlockKrylovSchur::isInitialized() will be false and the user will need
to provide a new initial iterate to the solver.

C++ includes: AnasaziBlockKrylovSchur.hpp ";

%feature("docstring")
Anasazi::BlockKrylovSchurInitFailure::BlockKrylovSchurInitFailure "Anasazi::BlockKrylovSchurInitFailure::BlockKrylovSchurInitFailure(const
std::string &what_arg) ";


// File: classAnasazi_1_1BlockKrylovSchurOrthoFailure.xml
%feature("docstring") Anasazi::BlockKrylovSchurOrthoFailure "

BlockKrylovSchurOrthoFailure is thrown when the orthogonalization
manager is unable to generate orthonormal columns from the new basis
vectors.

This exception is thrown from the BlockKrylovSchur::iterate() method.

C++ includes: AnasaziBlockKrylovSchur.hpp ";

%feature("docstring")
Anasazi::BlockKrylovSchurOrthoFailure::BlockKrylovSchurOrthoFailure "Anasazi::BlockKrylovSchurOrthoFailure::BlockKrylovSchurOrthoFailure(const
std::string &what_arg) ";


// File: classAnasazi_1_1BlockKrylovSchurSolMgr.xml
%feature("docstring") Anasazi::BlockKrylovSchurSolMgr "

The Anasazi::BlockKrylovSchurSolMgr provides a flexible solver manager
over the BlockKrylovSchur eigensolver.

The solver manager provides to the solver a StatusTestCombo object
constructed as follows: combo = globaltest OR debugtest  where
globaltest terminates computation when global convergence has been
detected.  It is encapsulated in a StatusTestWithOrdering object, to
ensure that computation is terminated only after the most significant
eigenvalues/eigenvectors have met the convergence criteria.  If not
specified via setGlobalStatusTest(), this test is a StatusTestResNorm
object which tests the 2-norms of the Ritz residuals relative to the
Ritz values.

debugtest allows a user to specify additional monitoring of the
iteration, encapsulated in a StatusTest object  If not specified via
setDebugStatusTest(), debugtest is ignored.  In most cases, it should
return Failed; if it returns Passed, solve() will throw an
AnasaziError exception.

Additionally, the solver manager will terminate solve() after a
specified number of restarts.

Much of this behavior is controlled via parameters and options passed
to the solver manager. For more information, see
BlockKrylovSchurSolMgr().

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist

C++ includes: AnasaziBlockKrylovSchurSolMgr.hpp ";


// File: structAnasazi_1_1BlockKrylovSchurState.xml
%feature("docstring") Anasazi::BlockKrylovSchurState "

Structure to contain pointers to BlockKrylovSchur state variables.

This struct is utilized by BlockKrylovSchur::initialize() and
BlockKrylovSchur::getState().

C++ includes: AnasaziBlockKrylovSchur.hpp ";

%feature("docstring")
Anasazi::BlockKrylovSchurState::BlockKrylovSchurState "Anasazi::BlockKrylovSchurState< ScalarType, MulVec
>::BlockKrylovSchurState() ";


// File: classAnasazi_1_1DenseMatTraits.xml
%feature("docstring") Anasazi::DenseMatTraits "

Virtual base class which defines basic traits for the multi-vector
type.

An adapter for this traits class must exist for the DM type. If not,
this class will produce a compile-time error.

C++ includes: AnasaziDenseMatTraits.hpp ";


// File: classAnasazi_1_1DirectSolver.xml
%feature("docstring") Anasazi::DirectSolver "

Anasazi's templated abstract base class providing solver capabilities
for projected eigenproblems.

This class provides concrete, templated implementations of solvers for
projected eigenproblems.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziDirectSolver.hpp ";


// File: classAnasazi_1_1Eigenproblem.xml
%feature("docstring") Anasazi::Eigenproblem "

This class defines the interface required by an eigensolver and status
test class to compute solutions to an eigenproblem.

C++ includes: AnasaziEigenproblem.hpp ";


// File: structAnasazi_1_1Eigensolution.xml
%feature("docstring") Anasazi::Eigensolution "

Struct for storing an eigenproblem solution.

C++ includes: AnasaziTypes.hpp ";

%feature("docstring")  Anasazi::Eigensolution::Eigensolution "Anasazi::Eigensolution< ScalarType, MV >::Eigensolution() ";


// File: classAnasazi_1_1Eigensolver.xml
%feature("docstring") Anasazi::Eigensolver "

The Eigensolver is a templated virtual base class that defines the
basic interface that any eigensolver will support.

This interface is mainly concerned with providing a set of eigensolver
status method that can be requested from any eigensolver by an
StatusTest object.

C++ includes: AnasaziEigensolverDecl.hpp ";


// File: classAnasazi_1_1EpetraGenOp.xml
%feature("docstring") Anasazi::EpetraGenOp "

Adapter class for creating an operators often used in solving
generalized eigenproblems.

This class will apply the operation $A^{-1}M$ [default] or $AM$, for
the Apply method of the Epetra_Operator / Anasazi::Operator. The
Anasazi::EpetraGenOp operator is useful when spectral transformations
are used within eigensolvers. For instance, $A^{-1}M$ is a shift and
invert spectral transformation commonly used with
Anasazi::BlockKrylovSchur to compute the smallest-magnitude
eigenvalues for the eigenproblem $Ax = \\\\lambda Mx$.

The Epetra package performs double-precision arithmetic, so the use of
Epetra with Anasazi will only provide a double-precision eigensolver.

C++ includes: AnasaziEpetraAdapter.hpp ";

%feature("docstring")  Anasazi::EpetraGenOp::EpetraGenOp "Anasazi::EpetraGenOp::EpetraGenOp(const Teuchos::RCP< Epetra_Operator
> &AOp, const Teuchos::RCP< Epetra_Operator > &MOp, bool
isAInverse=true)

Basic constructor for applying operator $A^{-1}M$ [default] or $AM$.

If isAInverse is true this operator will apply $A^{-1}M$, else it will
apply $AM$. ";

%feature("docstring")  Anasazi::EpetraGenOp::~EpetraGenOp "Anasazi::EpetraGenOp::~EpetraGenOp()

Destructor. ";

%feature("docstring")  Anasazi::EpetraGenOp::Apply "void
Anasazi::EpetraGenOp::Apply(const MultiVec< double > &X, MultiVec<
double > &Y) const

Apply method [inherited from Anasazi::Operator class].

This method will apply $A^{-1}M$ or $AM$ to X, returning Y. ";

%feature("docstring")  Anasazi::EpetraGenOp::Apply "int
Anasazi::EpetraGenOp::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Apply method [inherited from Epetra_Operator class].

This method will apply $A^{-1}M$ or $AM$ to X, returning Y. ";

%feature("docstring")  Anasazi::EpetraGenOp::ApplyInverse "int
Anasazi::EpetraGenOp::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Apply inverse method [inherited from Epetra_Operator class].

This method will apply $(A^{-1}M)^{-1}$ or $(AM)^{-1}$ to X, returning
Y. ";

%feature("docstring")  Anasazi::EpetraGenOp::Label "const char*
Anasazi::EpetraGenOp::Label() const

Returns a character string describing the operator. ";

%feature("docstring")  Anasazi::EpetraGenOp::UseTranspose "bool
Anasazi::EpetraGenOp::UseTranspose() const

Returns the current UseTranspose setting [always false for this
operator]. ";

%feature("docstring")  Anasazi::EpetraGenOp::SetUseTranspose "int
Anasazi::EpetraGenOp::SetUseTranspose(bool UseTranspose)

If set true, the transpose of this operator will be applied [not
functional for this operator]. ";

%feature("docstring")  Anasazi::EpetraGenOp::HasNormInf "bool
Anasazi::EpetraGenOp::HasNormInf() const

Returns true if this object can provide an approximate inf-norm
[always false for this operator]. ";

%feature("docstring")  Anasazi::EpetraGenOp::NormInf "double
Anasazi::EpetraGenOp::NormInf() const

Returns the infinity norm of the global matrix [not functional for
this operator]. ";

%feature("docstring")  Anasazi::EpetraGenOp::Comm "const Epetra_Comm&
Anasazi::EpetraGenOp::Comm() const

Returns the Epetra_Comm communicator associated with this operator. ";

%feature("docstring")  Anasazi::EpetraGenOp::OperatorDomainMap "const
Epetra_Map& Anasazi::EpetraGenOp::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Anasazi::EpetraGenOp::OperatorRangeMap "const
Epetra_Map& Anasazi::EpetraGenOp::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this
operator. ";


// File: classAnasazi_1_1EpetraMultiVec.xml
%feature("docstring") Anasazi::EpetraMultiVec "

Basic adapter class for Anasazi::MultiVec that uses
Epetra_MultiVector.

The Epetra package performs double-precision arithmetic, so the use of
Epetra with Anasazi will only provide a double-precision eigensolver.

C++ includes: AnasaziEpetraAdapter.hpp ";


// File: classAnasazi_1_1EpetraMultiVecFailure.xml
%feature("docstring") Anasazi::EpetraMultiVecFailure "

EpetraMultiVecFailure is thrown when a return value from an Epetra
call on an Epetra_MultiVector is non-zero.

C++ includes: AnasaziEpetraAdapter.hpp ";

%feature("docstring")
Anasazi::EpetraMultiVecFailure::EpetraMultiVecFailure "Anasazi::EpetraMultiVecFailure::EpetraMultiVecFailure(const
std::string &what_arg) ";


// File: classAnasazi_1_1EpetraOp.xml
%feature("docstring") Anasazi::EpetraOp "

Basic adapter class for Anasazi::Operator that uses Epetra_Operator.

The Epetra package performs double-precision arithmetic, so the use of
Epetra with Anasazi will only provide a double-precision eigensolver.

C++ includes: AnasaziEpetraAdapter.hpp ";


// File: classAnasazi_1_1EpetraOpFailure.xml
%feature("docstring") Anasazi::EpetraOpFailure "

EpetraOpFailure is thrown when a return value from an Epetra call on
an Epetra_Operator is non-zero.

C++ includes: AnasaziEpetraAdapter.hpp ";

%feature("docstring")  Anasazi::EpetraOpFailure::EpetraOpFailure "Anasazi::EpetraOpFailure::EpetraOpFailure(const std::string &what_arg)
";


// File: classAnasazi_1_1EpetraSymMVOp.xml
%feature("docstring") Anasazi::EpetraSymMVOp "

Adapter class for creating a symmetric operator from an
Epetra_MultiVector.

This class will apply the operation $A^TA$ [default] or $AA^T$, for
the Apply method of the Epetra_Operator / Anasazi::Operator. The
Anasazi::EpetraSymMvOp operator is useful when trying to compute a few
singular values of the Epetra_MultiVector $A$. The singular values are
the square-root of the eigenvalues of $A^TA$ and $AA^T$.

The Epetra package performs double-precision arithmetic, so the use of
Epetra with Anasazi will only provide a double-precision eigensolver.

C++ includes: AnasaziEpetraAdapter.hpp ";

%feature("docstring")  Anasazi::EpetraSymMVOp::EpetraSymMVOp "Anasazi::EpetraSymMVOp::EpetraSymMVOp(const Teuchos::RCP< const
Epetra_MultiVector > &MV, bool isTrans=false)

Basic constructor for applying operator $A^TA$ [default] or $AA^T$.

If isTrans is false this operator will apply $A^TA$, else it will
apply $AA^T$. ";

%feature("docstring")  Anasazi::EpetraSymMVOp::~EpetraSymMVOp "Anasazi::EpetraSymMVOp::~EpetraSymMVOp()

Destructor. ";

%feature("docstring")  Anasazi::EpetraSymMVOp::Apply "void
Anasazi::EpetraSymMVOp::Apply(const MultiVec< double > &X, MultiVec<
double > &Y) const

Apply method.

This method will apply $A^TA$ or $AA^T$ to X, returning Y. ";


// File: classAnasazi_1_1EpetraSymOp.xml
%feature("docstring") Anasazi::EpetraSymOp "

Adapter class for creating a symmetric operator from an
Epetra_Operator.

This class will apply the operation $A^TA$ [default] or $AA^T$, for
the Apply method of the Epetra_Operator / Anasazi::Operator. The
Anasazi::EpetraSymOp operator is useful when trying to compute a few
singular values of the operator $A$. The singular values are the
square-root of the eigenvalues of $A^TA$ and $AA^T$.

The Epetra package performs double-precision arithmetic, so the use of
Epetra with Anasazi will only provide a double-precision eigensolver.

C++ includes: AnasaziEpetraAdapter.hpp ";

%feature("docstring")  Anasazi::EpetraSymOp::EpetraSymOp "Anasazi::EpetraSymOp::EpetraSymOp(const Teuchos::RCP< Epetra_Operator
> &Op, bool isTrans=false)

Basic constructor for applying operator $A^TA$ [default] or $AA^T$.

If isTrans is false this operator will apply $A^TA$, else it will
apply $AA^T$. ";

%feature("docstring")  Anasazi::EpetraSymOp::~EpetraSymOp "Anasazi::EpetraSymOp::~EpetraSymOp()

Destructor. ";

%feature("docstring")  Anasazi::EpetraSymOp::Apply "void
Anasazi::EpetraSymOp::Apply(const MultiVec< double > &X, MultiVec<
double > &Y) const

Apply method [inherited from Anasazi::Operator class].

This method will apply $A^TA$ or $AA^T$ to X, returning Y. ";

%feature("docstring")  Anasazi::EpetraSymOp::Apply "int
Anasazi::EpetraSymOp::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Apply method [inherited from Epetra_Operator class].

This method will apply $A^TA$ or $AA^T$ to X, returning Y. ";

%feature("docstring")  Anasazi::EpetraSymOp::ApplyInverse "int
Anasazi::EpetraSymOp::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Apply inverse method [inherited from Epetra_Operator class].

This method will apply $(A^TA)^{-1}$ or $(AA^T)^{-1}$ to X, returning
Y. This method is only defined if $A^{-1}$ is defined for the given
Epetra_Operator. ";

%feature("docstring")  Anasazi::EpetraSymOp::Label "const char*
Anasazi::EpetraSymOp::Label() const

Returns a character string describing the operator. ";

%feature("docstring")  Anasazi::EpetraSymOp::UseTranspose "bool
Anasazi::EpetraSymOp::UseTranspose() const

Returns the current UseTranspose setting [always false for this
operator]. ";

%feature("docstring")  Anasazi::EpetraSymOp::SetUseTranspose "int
Anasazi::EpetraSymOp::SetUseTranspose(bool UseTranspose)

If set true, the transpose of this operator will be applied [not
functional for this operator]. ";

%feature("docstring")  Anasazi::EpetraSymOp::HasNormInf "bool
Anasazi::EpetraSymOp::HasNormInf() const

Returns true if this object can provide an approximate inf-norm
[always false for this operator]. ";

%feature("docstring")  Anasazi::EpetraSymOp::NormInf "double
Anasazi::EpetraSymOp::NormInf() const

Returns the infinity norm of the global matrix [not functional for
this operator]. ";

%feature("docstring")  Anasazi::EpetraSymOp::Comm "const Epetra_Comm&
Anasazi::EpetraSymOp::Comm() const

Returns the Epetra_Comm communicator associated with this operator. ";

%feature("docstring")  Anasazi::EpetraSymOp::OperatorDomainMap "const
Epetra_Map& Anasazi::EpetraSymOp::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Anasazi::EpetraSymOp::OperatorRangeMap "const
Epetra_Map& Anasazi::EpetraSymOp::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this
operator. ";


// File: classAnasazi_1_1EpetraWSymMVOp.xml
%feature("docstring") Anasazi::EpetraWSymMVOp "

Adapter class for creating a weighted symmetric operator from an
Epetra_MultiVector and Epetra_Operator.

This class will apply the operation $(WA)^T*WA$ for the Apply method
of the Anasazi::Operator. The Anasazi::EpetraWSymMvOp operator is
useful when trying to compute a few singular values of the
Epetra_MultiVector $A$ under the weighting matrix $W$. The singular
values are the square-root of the eigenvalues of $(WA)^T*WA$.

The Epetra package performs double-precision arithmetic, so the use of
Epetra with Anasazi will only provide a double-precision eigensolver.

C++ includes: AnasaziEpetraAdapter.hpp ";

%feature("docstring")  Anasazi::EpetraWSymMVOp::EpetraWSymMVOp "Anasazi::EpetraWSymMVOp::EpetraWSymMVOp(const Teuchos::RCP< const
Epetra_MultiVector > &MV, const Teuchos::RCP< Epetra_Operator > &OP)

Basic constructor for applying operator $A^TA$ [default] or $AA^T$.

If isTrans is false this operator will apply $A^TA$, else it will
apply $AA^T$. ";

%feature("docstring")  Anasazi::EpetraWSymMVOp::~EpetraWSymMVOp "Anasazi::EpetraWSymMVOp::~EpetraWSymMVOp()

Destructor. ";

%feature("docstring")  Anasazi::EpetraWSymMVOp::Apply "void
Anasazi::EpetraWSymMVOp::Apply(const MultiVec< double > &X, MultiVec<
double > &Y) const

Apply method.

This method will apply $(WA)^T*WA$ to X, returning Y. ";


// File: classAnasazi_1_1HelperTraits.xml
%feature("docstring") Anasazi::HelperTraits "

Class which defines basic traits for working with different scalar
types.

An adapter for this traits class must exist for the ScalarType. If
not, this class will produce a compile-time error.

C++ includes: AnasaziHelperTraits.hpp ";


// File: classAnasazi_1_1LOBPCG.xml
%feature("docstring") Anasazi::LOBPCG "

This class provides the Locally Optimal Block Preconditioned Conjugate
Gradient (LOBPCG) iteration, a preconditioned iteration for solving
linear Hermitian eigenproblems.

This implementation is a modification of the one found in A. Knyazev,
\"Toward the optimal preconditioned eigensolver: Locally optimal block
preconditioner conjugate gradient method\", SIAM J. Sci. Comput., vol
23, n 2, pp. 517-541.

The modification consists of the orthogonalization steps recommended
in U. Hetmaniuk and R. Lehoucq, \"Basis Selection in LOBPCG\", Journal
of Computational Physics.

These modifcation are referred to as full orthogonalization, and
consist of also conducting the local optimization using an orthonormal
basis.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist

C++ includes: AnasaziLOBPCG.hpp ";


// File: structAnasazi_1_1LOBPCG_1_1CheckList.xml


// File: classAnasazi_1_1LOBPCGInitFailure.xml
%feature("docstring") Anasazi::LOBPCGInitFailure "

LOBPCGInitFailure is thrown when the LOBPCG solver is unable to
generate an initial iterate in the LOBPCG::initialize() routine.

This exception is thrown from the LOBPCG::initialize() method, which
is called by the user or from the LOBPCG::iterate() method when
isInitialized() == false.

In the case that this exception is thrown, LOBPCG::hasP() and
LOBPCG::isInitialized() will be false and the user will need to
provide a new initial iterate to the solver.

C++ includes: AnasaziLOBPCG.hpp ";

%feature("docstring")  Anasazi::LOBPCGInitFailure::LOBPCGInitFailure "Anasazi::LOBPCGInitFailure::LOBPCGInitFailure(const std::string
&what_arg) ";


// File: classAnasazi_1_1LOBPCGOrthoFailure.xml
%feature("docstring") Anasazi::LOBPCGOrthoFailure "

LOBPCGOrthoFailure is thrown when an orthogonalization attempt fails.

This is thrown in one of two scenarstd::ios. After preconditioning the
residual, the orthogonalization manager is asked to orthogonalize the
preconditioned residual (H) against the auxiliary vectors. If full
orthogonalization is enabled, H is also orthogonalized against X and P
and normalized.

The second scenario involves the generation of new X and P from the
basis [X H P]. When full orthogonalization is enabled, an attempt is
made to select coefficients for X and P so that they will be mutually
orthogonal and orthonormal.

If either of these attempts fail, the solver throws an
LOBPCGOrthoFailure exception.

C++ includes: AnasaziLOBPCG.hpp ";

%feature("docstring")  Anasazi::LOBPCGOrthoFailure::LOBPCGOrthoFailure
"Anasazi::LOBPCGOrthoFailure::LOBPCGOrthoFailure(const std::string
&what_arg) ";


// File: classAnasazi_1_1LOBPCGRitzFailure.xml
%feature("docstring") Anasazi::LOBPCGRitzFailure "

LOBPCGRitzFailure is thrown when the LOBPCG solver is unable to
continue a call to LOBPCG::iterate() due to a failure of the
algorithm.

This signals that the Rayleigh-Ritz analysis over the subspace
colsp([X H P]) detected ill-conditioning of the projected mass matrix
and the inability to generate a set of orthogonal eigenvectors for the
projected problem.

This exception is only thrown from the LOBPCG::iterate() routine.
After catching this exception, the user can recover the subspace via
LOBPCG::getState(). This information can be used to restart the
solver.

C++ includes: AnasaziLOBPCG.hpp ";

%feature("docstring")  Anasazi::LOBPCGRitzFailure::LOBPCGRitzFailure "Anasazi::LOBPCGRitzFailure::LOBPCGRitzFailure(const std::string
&what_arg) ";


// File: classAnasazi_1_1LOBPCGSolMgr.xml
%feature("docstring") Anasazi::LOBPCGSolMgr "

The LOBPCGSolMgr provides a powerful solver manager over the LOBPCG
eigensolver.

This solver manager exists to provide a flexible manager over the
LOBPCG eigensolver intended for general use. Features provided by this
solver manager include: locking of converged eigenpairs

global convergence on only the significant eigenpairs (instead of any
eigenpairs with low residual)

recovery from LOBPCGRitzFailure when full orthogonalization is
disabled

The solver manager provides to the solver a StatusTestCombo object
constructed as follows: combo = maxiterstest OR globaltest OR
lockingtest OR debugtest  where  maxiters terminates computation when
a maximum number of iterations have been performed maxiters is a
StatusTestMaxIters object

globaltest terminates computation when global convergence has been
detected.  It is encapsulated in a StatusTestWithOrdering object, to
ensure that computation is terminated only after the most significant
eigenvalues/eigenvectors have met the convergence criteria.  If not
specified via setGlobalStatusTest(), globaltest is a StatusTestResNorm
object which tests the M-norms of the direct residuals relative to the
Ritz values.

lockingtest halts LOBPCG::iterate() in order to deflate converged
eigenpairs for locking.  It will query the underlying LOBPCG
eigensolver to determine when eigenvectors should be locked.  If not
specified via setLockingStatusTest(), lockingtest is a
StatusTestResNorm object.

debugtest allows a user to specify additional monitoring of the
iteration, encapsulated in a StatusTest object  If not specified via
setDebugStatusTest(), debugtest is ignored.  In most cases, it should
return Failed; if it returns Passed, solve() will throw an
AnasaziError exception.

Much of this behavior is controlled via parameters and options passed
to the solver manager. For more information, see LOBPCGSolMgr().

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist

C++ includes: AnasaziLOBPCGSolMgr.hpp ";


// File: structAnasazi_1_1LOBPCGState.xml
%feature("docstring") Anasazi::LOBPCGState "

Structure to contain pointers to Anasazi state variables.

This struct is utilized by LOBPCG::initialize() and
LOBPCG::getState().

C++ includes: AnasaziLOBPCG.hpp ";

%feature("docstring")  Anasazi::LOBPCGState::LOBPCGState "Anasazi::LOBPCGState< ScalarType, MultiVector >::LOBPCGState() ";


// File: classAnasazi_1_1MatOrthoManager.xml
%feature("docstring") Anasazi::MatOrthoManager "

Anasazi's templated virtual class for providing routines for
orthogonalization and orthonormalization of multivectors using matrix-
based inner products.

This class extends Anasazi::OrthoManager by providing extra calling
arguments to orthogonalization routines, to reduce the cost of
applying the inner product in cases where the user already has the
image of the source multivector under the inner product matrix.

A concrete implementation of this class is necessary. The user can
create their own implementation if those supplied are not suitable for
their needs.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziMatOrthoManager.hpp ";


// File: classAnasazi_1_1MultiVec.xml
%feature("docstring") Anasazi::MultiVec "

Anasazi's templated virtual class for constructing a multi-vector that
can interface with the MultiVecTraits class used by the eigensolvers.

A concrete implementation of this class is necessary. The user can
create their own implementation if those supplied are not suitable for
their needs.

Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziMultiVec.hpp ";


// File: classAnasazi_1_1MultiVecTraits.xml
%feature("docstring") Anasazi::MultiVecTraits "

Virtual base class which defines basic traits for the multi-vector
type.

An adapter for this traits class must exist for the MV type. If not,
this class will produce a compile-time error.

C++ includes: AnasaziMultiVecTraits.hpp ";


// File: classAnasazi_1_1MultiVecTraits_3_01double_00_01Epetra__MultiVector_01_4.xml
%feature("docstring") Anasazi::MultiVecTraits< double,
Epetra_MultiVector > "

Template specialization of Anasazi::MultiVecTraits class using the
Epetra_MultiVector class.

This interface will ensure that any Epetra_MultiVector will be
accepted by the Anasazi templated solvers.

The Epetra package performs double-precision arithmetic, so the use of
Epetra with Anasazi will only provide a double-precision eigensolver.

C++ includes: AnasaziEpetraAdapter.hpp ";


// File: classAnasazi_1_1MultiVecTraits_3_01ScalarType_00_01MultiVec_3_01ScalarType_01_4_01_4.xml
%feature("docstring") Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > > "

Template specialization of Anasazi::MultiVecTraits class using the
Anasazi::MultiVec virtual base class.

Any class that inherits from Anasazi::MultiVec will be accepted by the
Anasazi templated solvers due to this interface to the
Anasazi::MultiVecTraits class.

C++ includes: AnasaziMultiVec.hpp ";


// File: classAnasazi_1_1Operator.xml
%feature("docstring") Anasazi::Operator "

Anasazi's templated virtual class for constructing an operator that
can interface with the OperatorTraits class used by the eigensolvers.

A concrete implementation of this class is necessary. The user can
create their own implementation if those supplied are not suitable for
their needs.

Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziOperator.hpp ";


// File: classAnasazi_1_1OperatorError.xml
%feature("docstring") Anasazi::OperatorError "

Exceptions thrown to signal error in operator application.

C++ includes: AnasaziOperatorTraits.hpp ";

%feature("docstring")  Anasazi::OperatorError::OperatorError "Anasazi::OperatorError::OperatorError(const std::string &what_arg) ";


// File: classAnasazi_1_1OperatorTraits.xml
%feature("docstring") Anasazi::OperatorTraits "

Virtual base class which defines basic traits for the operator type.

An adapter for this traits class must exist for the MV and OP types.
If not, this class will produce a compile-time error.

C++ includes: AnasaziOperatorTraits.hpp ";


// File: classAnasazi_1_1OperatorTraits_3_01double_00_01Epetra__MultiVector_00_01Epetra__Operator_01_4.xml
%feature("docstring") Anasazi::OperatorTraits< double,
Epetra_MultiVector, Epetra_Operator > "

Template specialization of Anasazi::OperatorTraits class using the
Epetra_Operator virtual base class and Epetra_MultiVector class.

This interface will ensure that any Epetra_Operator and
Epetra_MultiVector will be accepted by the Anasazi templated solvers.

The Epetra package performs double-precision arithmetic, so the use of
Epetra with Anasazi will only provide a double-precision eigensolver.

C++ includes: AnasaziEpetraAdapter.hpp ";


// File: classAnasazi_1_1OperatorTraits_3_01ScalarType_00_01MultiVec_3_01ScalarType_01_4_00_01Operator_3_01ScalarType_01_4_01_4.xml
%feature("docstring") Anasazi::OperatorTraits< ScalarType, MultiVec<
ScalarType >, Operator< ScalarType > > "

Template specialization of Anasazi::OperatorTraits class using
Anasazi::Operator and Anasazi::MultiVec virtual base classes.

Any class that inherits from Anasazi::Operator will be accepted by the
Anasazi templated solvers due to this interface to the
Anasazi::OperatorTraits class.

C++ includes: AnasaziOperator.hpp ";


// File: classAnasazi_1_1OrthoError.xml
%feature("docstring") Anasazi::OrthoError "

Exception thrown to signal error in an orthogonalization manager
method.

C++ includes: AnasaziOrthoManager.hpp ";

%feature("docstring")  Anasazi::OrthoError::OrthoError "Anasazi::OrthoError::OrthoError(const std::string &what_arg) ";


// File: classAnasazi_1_1OrthoManager.xml
%feature("docstring") Anasazi::OrthoManager "

Anasazi's templated virtual class for providing routines for
orthogonalization and orthonormalization of multivectors.

This class defines concepts of orthogonality through the definition of
an inner product. It also provides computational routines for
orthogonalization.

A concrete implementation of this class is necessary. The user can
create their own implementation if those supplied are not suitable for
their needs.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziOrthoManager.hpp ";


// File: classAnasazi_1_1OutputManager.xml
%feature("docstring") Anasazi::OutputManager "

Output managers remove the need for the eigensolver to know any
information about the required output. Calling isVerbosity( MsgType
type ) informs the solver if it is supposed to output the information
corresponding to the message type.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziOutputManager.hpp ";


// File: classAnasazi_1_1ResNormNaNError.xml
%feature("docstring") Anasazi::ResNormNaNError "

ResNormNaNError is thrown from StatusTestResNorm::checkStatus() when a
NaN (\"not a number\") is detected among the residual norms returned
by the eigensolver.

This behavior is optional and is controlled by flag to
StatusTestResNorm::StatusTestResNorm().

C++ includes: AnasaziStatusTestResNorm.hpp ";

%feature("docstring")  Anasazi::ResNormNaNError::ResNormNaNError "Anasazi::ResNormNaNError::ResNormNaNError(const std::string &what_arg)
";


// File: classAnasazi_1_1SimpleLOBPCGSolMgr.xml
%feature("docstring") Anasazi::SimpleLOBPCGSolMgr "

The Anasazi::SimpleLOBPCGSolMgr provides a simple solver manager over
the LOBPCG eigensolver.

Anasazi::SimpleLOBPCGSolMgr allows the user to specify convergence
tolerance, verbosity level and block size. When block size is less
than the number of requested eigenvalues specified in the
eigenproblem, checkpointing is activated.

The purpose of this solver manager was to provide an example of a
simple solver manager, useful for demonstration as well as a jumping-
off point for solvermanager development. Also, the solver manager is
useful for testing some of the features of the Anasazi::LOBPCG
eigensolver, principally the use of auxiliary vectors.

This solver manager does not verify before quitting that the nev
eigenvectors that have converged are also the smallest nev
eigenvectors that are known.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist

C++ includes: AnasaziSimpleLOBPCGSolMgr.hpp ";


// File: classAnasazi_1_1SolverManager.xml
%feature("docstring") Anasazi::SolverManager "

The Anasazi::SolverManager is a templated virtual base class that
defines the basic interface that any solver manager will support.

C++ includes: AnasaziSolverManager.hpp ";


// File: classAnasazi_1_1SolverUtils.xml
%feature("docstring") Anasazi::SolverUtils "

Anasazi's templated, static class providing utilities for the solvers.

This class provides concrete, templated implementations of utilities
necessary for the solvers. These utilities include sorting,
orthogonalization, projecting/solving local eigensystems, and sanity
checking. These are internal utilties, so the user should not alter
this class.

Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziSolverUtils.hpp ";


// File: classAnasazi_1_1SortManager.xml
%feature("docstring") Anasazi::SortManager "

Anasazi's templated pure virtual class for managing the sorting of
approximate eigenvalues computed by the eigensolver.

A concrete implementation of this class is necessary. The user can
create their own implementation if those supplied are not suitable for
their needs.

Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziSortManager.hpp ";

%feature("docstring")  Anasazi::SortManager::SortManager "Anasazi::SortManager< ScalarType, MV, OP >::SortManager()

Default Constructor. ";

%feature("docstring")  Anasazi::SortManager::~SortManager "virtual
Anasazi::SortManager< ScalarType, MV, OP >::~SortManager()

Destructor. ";

%feature("docstring")  Anasazi::SortManager::sort "virtual void
Anasazi::SortManager< ScalarType, MV, OP >::sort(Eigensolver<
ScalarType, MV, OP > *solver, const int n, std::vector< typename
Teuchos::ScalarTraits< ScalarType >::magnitudeType > &evals,
std::vector< int > *perm=0) const =0

Sort the vector of eigenvalues, optionally returning the permutation
vector.

Parameters:
-----------

solver:  [in] Eigensolver that is calling the sorting routine

n:  [in] Number of values in evals to be sorted.

evals:  [in/out] Vector of length n containing the eigenvalues to be
sorted

perm:  [out] Vector of length n to store the permutation index
(optional) ";

%feature("docstring")  Anasazi::SortManager::sort "virtual void
Anasazi::SortManager< ScalarType, MV, OP >::sort(Eigensolver<
ScalarType, MV, OP > *solver, const int n, std::vector< typename
Teuchos::ScalarTraits< ScalarType >::magnitudeType > &r_evals,
std::vector< typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType > &i_evals, std::vector< int > *perm=0) const =0

Sort the vectors of eigenpairs, optionally returning the permutation
vector.

This routine takes two vectors, one for each part of a complex
eigenvalue. This is helpful for solving real, non-symmetric eigenvalue
problems.

Parameters:
-----------

solver:  [in] Eigensolver that is calling the sorting routine

n:  [in] Number of values in r_evals,i_evals to be sorted.

r_evals:  [in/out] Vector of length n containing the real part of the
eigenvalues to be sorted

i_evals:  [in/out] Vector of length n containing the imaginary part of
the eigenvalues to be sorted

perm:  [out] Vector of length n to store the permutation index
(optional) ";


// File: classAnasazi_1_1SortManagerError.xml
%feature("docstring") Anasazi::SortManagerError "

SortManagerError is thrown when the Anasazi::SortManager is unable to
sort the numbers, due to some failure of the sort method or error in
calling it.

C++ includes: AnasaziSortManager.hpp ";

%feature("docstring")  Anasazi::SortManagerError::SortManagerError "Anasazi::SortManagerError::SortManagerError(const std::string
&what_arg) ";


// File: classAnasazi_1_1StatusTest.xml
%feature("docstring") Anasazi::StatusTest "

A pure virtual class for defining the status tests for the Anasazi
iterative solvers.

StatusTest is an interface that can be implemented to create
convergence tests for all Anasazi solvers. Almost any kind of test can
be expressed using this mechanism, including composite tests (see
StatusTestCombo).

C++ includes: AnasaziStatusTestDecl.hpp ";


// File: classAnasazi_1_1StatusTestCombo.xml
%feature("docstring") Anasazi::StatusTestCombo "

Status test for forming logical combinations of other status tests.

Test types include StatusTestCombo::OR, StatusTestCombo::AND,
StatusTestCombo::SEQOR and StatusTestCombo::SEQAND. The
StatusTestCombo::OR and StatusTestCombo::AND tests evaluate all of the
tests, in the order they were passed to the StatusTestCombo. The
StatusTestCombo::SEQOR and StatusTestCombo::SEQAND run only the tests
necessary to determine the final outcome, short- circuiting on the
first test that conclusively decides the outcome. More formally,
StatusTestCombo::SEQAND runs the tests in the order they were given to
the StatusTestCombo class and stops after the first test that
evaluates Failed. StatusTestCombo::SEQOR run the tests in the order
they were given to the StatusTestCombo class and stops after the first
test that evaluates Passed.

C++ includes: AnasaziStatusTestCombo.hpp ";


// File: classAnasazi_1_1StatusTestError.xml
%feature("docstring") Anasazi::StatusTestError "

Exception thrown to signal error in a status test during
Anasazi::StatusTest::checkStatus().

C++ includes: AnasaziStatusTest.hpp ";

%feature("docstring")  Anasazi::StatusTestError::StatusTestError "Anasazi::StatusTestError::StatusTestError(const std::string &what_arg)
";


// File: classAnasazi_1_1StatusTestMaxIters.xml
%feature("docstring") Anasazi::StatusTestMaxIters "

A status test for testing the number of iterations.

Anasazi::StatusTestMaxIters will test true when an eigensolver has
reached some number of iterations. Specifically,
{ Passed,  if solver->getNumIters() >= maxIter status(solver) = {
{ Failed,  if solver->getNumIters()  < maxIter where maxIter is the
parameter given to the status tester.

This status test also supports negation, so that it negates the need
for a StatusTestMinIters status tester. In this way, all tests on the
range of iterations can be constructed through the appropriate use of
StatusTestMaxIters and StatusTestCombo.

C++ includes: AnasaziStatusTestMaxIters.hpp ";


// File: classAnasazi_1_1StatusTestOutput.xml
%feature("docstring") Anasazi::StatusTestOutput "

A special StatusTest for printing other status tests.

StatusTestOutput is a wrapper around another StatusTest that calls
StatusTest::print() on the underlying object on calls to
StatusTestOutput::checkStatus(). The frequency and occasion of the
printing can be dictated according to some parameters passed to
StatusTestOutput::StatusTestOutput().

C++ includes: AnasaziStatusTestOutput.hpp ";


// File: classAnasazi_1_1StatusTestResNorm.xml
%feature("docstring") Anasazi::StatusTestResNorm "

A status test for testing the norm of the eigenvectors residuals.

StatusTestResNorm was designed to be used as a test for convergence.
The tester compares the norms of the residual vectors against a user
specified tolerance.

In addition to specifying the tolerance, the user may specify: the
norm to be used: 2-norm or OrthoManager::norm() or
Eigensolver::getRitzRes2Norms()

the scale: absolute or relative to magnitude of Ritz value

the quorum: the number of vectors required for the test to evaluate as
Passed.

C++ includes: AnasaziStatusTestResNorm.hpp ";


// File: classAnasazi_1_1StatusTestWithOrdering.xml
%feature("docstring") Anasazi::StatusTestWithOrdering "

A status test for testing the norm of the eigenvectors residuals along
with a set of auxiliary eigenvalues.

The test evaluates to Passed when then the most significant of the
eigenvalues all have a residual below a certain threshhold. The
purpose of the test is to not only test convergence for some number of
eigenvalues, but to test convergence for the correct ones.

In addition to specifying the tolerance, the user may specify: the
norm to be used: 2-norm or OrthoManager::norm() or getRitzRes2Norms()

the scale: absolute or relative to magnitude of Ritz value

the quorum: the number of vectors required for the test to evaluate as
Passed.

Finally, the user must specify the Anasazi::SortManager used for
deciding significance.

C++ includes: AnasaziStatusTestWithOrdering.hpp ";


// File: classAnasazi_1_1SVQBOrthoManager.xml
%feature("docstring") Anasazi::SVQBOrthoManager "

An implementation of the Anasazi::MatOrthoManager that performs
orthogonalization using the SVQB iterative orthogonalization technique
described by Stathapoulos and Wu. This orthogonalization routine,
while not returning the upper triangular factors of the popular Gram-
Schmidt method, has a communication cost (measured in number of
communication calls) that is independent of the number of columns in
the basis.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziSVQBOrthoManager.hpp ";


// File: structAnasazi_1_1UndefinedDenseMatTraits.xml
%feature("docstring") Anasazi::UndefinedDenseMatTraits "

This is the default struct used by DenseMatrixTraits<OrdinalType,
ScalarType> class to produce a compile time error when the
specialization does not exist for dense matrix type DM.

C++ includes: AnasaziDenseMatTraits.hpp ";


// File: structAnasazi_1_1UndefinedMultiVecTraits.xml
%feature("docstring") Anasazi::UndefinedMultiVecTraits "

This is the default struct used by MultiVecTraits<ScalarType, MV>
class to produce a compile time error when the specialization does not
exist for multivector type MV.

C++ includes: AnasaziMultiVecTraits.hpp ";


// File: structAnasazi_1_1UndefinedOperatorTraits.xml
%feature("docstring") Anasazi::UndefinedOperatorTraits "

This is the default struct used by OperatorTraits<ScalarType, MV, OP>
class to produce a compile time error when the specialization does not
exist for operator type OP.

C++ includes: AnasaziOperatorTraits.hpp ";


// File: structAnasazi_1_1Value.xml
%feature("docstring") Anasazi::Value "

This struct is used for storing eigenvalues and Ritz values, as a pair
of real values.

C++ includes: AnasaziTypes.hpp ";

%feature("docstring")  Anasazi::Value::set "void Anasazi::Value<
ScalarType >::set(const typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType &rp, const typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType &ip) ";


// File: namespaceAnasazi.xml
%feature("docstring")  Anasazi::Anasazi_Version "std::string
Anasazi::Anasazi_Version() ";

%feature("docstring")  Anasazi::TestMultiVecTraits "bool
Anasazi::TestMultiVecTraits(const Teuchos::RCP< OutputManager<
ScalarType > > &om, const Teuchos::RCP< const MV > &A)

This is a function to test the correctness of a MultiVecTraits
specialization and multivector implementation.

Status of the test: true is success, false is error ";

%feature("docstring")  Anasazi::TestOperatorTraits "bool
Anasazi::TestOperatorTraits(const Teuchos::RCP< OutputManager<
ScalarType > > &om, const Teuchos::RCP< const MV > &A, const
Teuchos::RCP< const OP > &M)

This function tests the correctness of an operator implementation with
respect to an OperatorTraits specialization.

Status of the test: true is successful, false otherwise. ";


// File: AnasaziBasicEigenproblem_8hpp.xml


// File: AnasaziBasicOrthoManager_8hpp.xml


// File: AnasaziBasicOutputManager_8hpp.xml


// File: AnasaziBasicSort_8hpp.xml


// File: AnasaziBlockDavidson_8hpp.xml


// File: AnasaziBlockDavidsonSolMgr_8hpp.xml


// File: AnasaziBlockKrylovSchur_8hpp.xml


// File: AnasaziBlockKrylovSchurSolMgr_8hpp.xml


// File: AnasaziConfigDefs_8hpp.xml


// File: AnasaziDenseMatTraits_8hpp.xml


// File: AnasaziDirectSolver_8hpp.xml


// File: AnasaziEigenproblem_8hpp.xml


// File: AnasaziEigensolver_8hpp.xml


// File: AnasaziEigensolverDecl_8hpp.xml


// File: AnasaziEpetraAdapter_8cpp.xml


// File: AnasaziEpetraAdapter_8hpp.xml


// File: AnasaziHelperTraits_8hpp.xml


// File: AnasaziLOBPCG_8hpp.xml


// File: AnasaziLOBPCGSolMgr_8hpp.xml


// File: AnasaziMatOrthoManager_8hpp.xml


// File: AnasaziMultiVec_8hpp.xml


// File: AnasaziMultiVecTraits_8hpp.xml


// File: AnasaziMVOPTester_8hpp.xml


// File: AnasaziOperator_8hpp.xml


// File: AnasaziOperatorTraits_8hpp.xml


// File: AnasaziOrthoManager_8hpp.xml


// File: AnasaziOutputManager_8hpp.xml


// File: AnasaziSimpleLOBPCGSolMgr_8hpp.xml


// File: AnasaziSolverManager_8hpp.xml


// File: AnasaziSolverUtils_8hpp.xml


// File: AnasaziSortManager_8hpp.xml


// File: AnasaziStatusTest_8hpp.xml


// File: AnasaziStatusTestCombo_8hpp.xml


// File: AnasaziStatusTestDecl_8hpp.xml


// File: AnasaziStatusTestMaxIters_8hpp.xml


// File: AnasaziStatusTestOutput_8hpp.xml


// File: AnasaziStatusTestResNorm_8hpp.xml


// File: AnasaziStatusTestWithOrdering_8hpp.xml


// File: AnasaziSVQBOrthoManager_8hpp.xml


// File: AnasaziTypes_8hpp.xml


// File: AnasaziVersion_8cpp.xml


// File: dir_dc5a21cdf9ae5c2aac185b6aad0b5f42.xml


// File: dir_7c2139b3455d9e99d84c8cdfb950729b.xml


// File: BlockDavidson_2BlockDavidsonEpetraEx_8cpp-example.xml


// File: BlockDavidson_2BlockDavidsonEpetraExGen_8cpp-example.xml


// File: BlockDavidson_2BlockDavidsonEpetraExGenPrecIfpack_8cpp-example.xml


// File: BlockKrylovSchur_2BlockKrylovSchurEpetraEx_8cpp-example.xml


// File: BlockKrylovSchur_2BlockKrylovSchurEpetraExGenAmesos_8cpp-example.xml


// File: BlockKrylovSchur_2BlockKrylovSchurEpetraExGenAztecOO_8cpp-example.xml


// File: BlockKrylovSchur_2BlockKrylovSchurEpetraExGenBelos_8cpp-example.xml


// File: BlockKrylovSchur_2BlockKrylovSchurEpetraExSVD_8cpp-example.xml


// File: LOBPCG_2LOBPCGEpetraEx_8cpp-example.xml


// File: LOBPCG_2LOBPCGEpetraExGen_8cpp-example.xml


// File: LOBPCG_2LOBPCGEpetraExGenPrecIfpack_8cpp-example.xml


// File: LOBPCG_2LOBPCGEpetraExSimple_8cpp-example.xml


// File: MVOPTester_2MVOPTesterEx_8cpp-example.xml


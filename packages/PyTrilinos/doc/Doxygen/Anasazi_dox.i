
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

/*  Constructors/Destructor  */

%feature("docstring")  Anasazi::BasicEigenproblem::BasicEigenproblem "Anasazi::BasicEigenproblem< ScalarType, MV, OP >::BasicEigenproblem()

Empty constructor - allows Anasazi::BasicEigenproblem to be described
at a later time through \"Set Methods\". ";

%feature("docstring")  Anasazi::BasicEigenproblem::BasicEigenproblem "Anasazi::BasicEigenproblem< ScalarType, MV, OP
>::BasicEigenproblem(const Teuchos::RCP< const OP > &Op, const
Teuchos::RCP< MV > &InitVec)

Standard Eigenvalue Problem Constructor. ";

%feature("docstring")  Anasazi::BasicEigenproblem::BasicEigenproblem "Anasazi::BasicEigenproblem< ScalarType, MV, OP
>::BasicEigenproblem(const Teuchos::RCP< const OP > &Op, const
Teuchos::RCP< const OP > &B, const Teuchos::RCP< MV > &InitVec)

Generalized Eigenvalue Problem Constructor. ";

%feature("docstring")  Anasazi::BasicEigenproblem::BasicEigenproblem "Anasazi::BasicEigenproblem< ScalarType, MV, OP
>::BasicEigenproblem(const BasicEigenproblem< ScalarType, MV, OP >
&Problem)

Copy Constructor. ";

%feature("docstring")  Anasazi::BasicEigenproblem::~BasicEigenproblem
"virtual Anasazi::BasicEigenproblem< ScalarType, MV, OP
>::~BasicEigenproblem()

Destructor. ";

/*  Set Methods  */

%feature("docstring")  Anasazi::BasicEigenproblem::setOperator "void
Anasazi::BasicEigenproblem< ScalarType, MV, OP >::setOperator(const
Teuchos::RCP< const OP > &Op)

Set the operator for which eigenvalues will be computed.

This may be different from the A if a spectral transformation is
employed. For example, this operator may apply the operation
$(A-\\\\sigma I)^{-1}$ if you are looking for eigenvalues of A around
$\\\\sigma$. ";

%feature("docstring")  Anasazi::BasicEigenproblem::setA "void
Anasazi::BasicEigenproblem< ScalarType, MV, OP >::setA(const
Teuchos::RCP< const OP > &A)

Set the operator A of the eigenvalue problem $Ax=Mx\\\\lambda$. ";

%feature("docstring")  Anasazi::BasicEigenproblem::setM "void
Anasazi::BasicEigenproblem< ScalarType, MV, OP >::setM(const
Teuchos::RCP< const OP > &M)

Set the operator M of the eigenvalue problem $Ax = Mx\\\\lambda$. ";

%feature("docstring")  Anasazi::BasicEigenproblem::setPrec "void
Anasazi::BasicEigenproblem< ScalarType, MV, OP >::setPrec(const
Teuchos::RCP< const OP > &Prec)

Set the preconditioner for this eigenvalue problem $Ax =
Mx\\\\lambda$. ";

%feature("docstring")  Anasazi::BasicEigenproblem::setInitVec "void
Anasazi::BasicEigenproblem< ScalarType, MV, OP >::setInitVec(const
Teuchos::RCP< MV > &InitVec)

Set the initial guess.

This vector is required to create all the space needed by Anasazi to
solve the eigenvalue problem.

Even if an initial guess is not known by the user, an initial vector
must be passed in. ";

%feature("docstring")  Anasazi::BasicEigenproblem::setAuxVecs "void
Anasazi::BasicEigenproblem< ScalarType, MV, OP >::setAuxVecs(const
Teuchos::RCP< const MV > &AuxVecs)

Set auxiliary vectors.

This multivector can have any number of columns, and most likely will
contain vectors that will be used by the eigensolver to orthogonalize
against. ";

%feature("docstring")  Anasazi::BasicEigenproblem::setNEV "void
Anasazi::BasicEigenproblem< ScalarType, MV, OP >::setNEV(int nev)

Specify the number of eigenvalues (NEV) that are requested. ";

%feature("docstring")  Anasazi::BasicEigenproblem::setHermitian "void
Anasazi::BasicEigenproblem< ScalarType, MV, OP >::setHermitian(bool
isSym)

Specify the symmetry of this eigenproblem.

This knowledge may allow the solver to take advantage of the
eigenproblems' symmetry. Some computational work can be avoided by
setting this properly. ";

%feature("docstring")  Anasazi::BasicEigenproblem::setProblem "bool
Anasazi::BasicEigenproblem< ScalarType, MV, OP >::setProblem()

Specify that this eigenproblem is fully defined.

This routine serves multiple purpose: sanity check that the
eigenproblem has been fully and consistently defined

opportunity for the eigenproblem to allocate internal storage for
eigenvalues and eigenvectors (to be used by eigensolvers and solver
managers)

This method reallocates internal storage, so that any previously
retrieved references to internal storage (eigenvectors or eigenvalues)
are invalidated.

The user MUST call this routine before they send the eigenproblem to
any solver or solver manager.

true signifies success, false signifies error. ";

%feature("docstring")  Anasazi::BasicEigenproblem::setSolution "void
Anasazi::BasicEigenproblem< ScalarType, MV, OP >::setSolution(const
Eigensolution< ScalarType, MV > &sol)

Set the solution to the eigenproblem.

This mechanism allows an Eigensolution struct to be associated with an
Eigenproblem object. setSolution() is usually called by a solver
manager at the end of its SolverManager::solve() routine. ";

/*  Accessor Methods  */

%feature("docstring")  Anasazi::BasicEigenproblem::getOperator "Teuchos::RCP<const OP> Anasazi::BasicEigenproblem< ScalarType, MV, OP
>::getOperator() const

Get a pointer to the operator for which eigenvalues will be computed.
";

%feature("docstring")  Anasazi::BasicEigenproblem::getA "Teuchos::RCP<const OP> Anasazi::BasicEigenproblem< ScalarType, MV, OP
>::getA() const

Get a pointer to the operator A of the eigenproblem $Ax=\\\\lambda
Mx$. ";

%feature("docstring")  Anasazi::BasicEigenproblem::getM "Teuchos::RCP<const OP> Anasazi::BasicEigenproblem< ScalarType, MV, OP
>::getM() const

Get a pointer to the operator M of the eigenproblem $Ax=\\\\lambda
Mx$. ";

%feature("docstring")  Anasazi::BasicEigenproblem::getPrec "Teuchos::RCP<const OP> Anasazi::BasicEigenproblem< ScalarType, MV, OP
>::getPrec() const

Get a pointer to the preconditioner of the eigenproblem $Ax=\\\\lambda
Mx$. ";

%feature("docstring")  Anasazi::BasicEigenproblem::getInitVec "Teuchos::RCP<const MV> Anasazi::BasicEigenproblem< ScalarType, MV, OP
>::getInitVec() const

Get a pointer to the initial vector. ";

%feature("docstring")  Anasazi::BasicEigenproblem::getAuxVecs "Teuchos::RCP<const MV> Anasazi::BasicEigenproblem< ScalarType, MV, OP
>::getAuxVecs() const

Get a pointer to the auxiliary vector. ";

%feature("docstring")  Anasazi::BasicEigenproblem::getNEV "int
Anasazi::BasicEigenproblem< ScalarType, MV, OP >::getNEV() const

Get the number of eigenvalues (NEV) that are required by this
eigenproblem. ";

%feature("docstring")  Anasazi::BasicEigenproblem::isHermitian "bool
Anasazi::BasicEigenproblem< ScalarType, MV, OP >::isHermitian() const

Get the symmetry information for this eigenproblem. ";

%feature("docstring")  Anasazi::BasicEigenproblem::isProblemSet "bool
Anasazi::BasicEigenproblem< ScalarType, MV, OP >::isProblemSet() const

If the problem has been set, this method will return true. ";

%feature("docstring")  Anasazi::BasicEigenproblem::getSolution "const
Eigensolution<ScalarType,MV>& Anasazi::BasicEigenproblem< ScalarType,
MV, OP >::getSolution() const

Get the solution to the eigenproblem.

There is no computation associated with this method. It only provides
a mechanism for associating an Eigensolution with a Eigenproblem. ";


// File: classAnasazi_1_1BasicOrthoManager.xml
%feature("docstring") Anasazi::BasicOrthoManager "

An implementation of the Anasazi::MatOrthoManager that performs
orthogonalization using (potentially) multiple steps of classical
Gram-Schmidt.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziBasicOrthoManager.hpp ";

/*  Constructor/Destructor  */

%feature("docstring")  Anasazi::BasicOrthoManager::BasicOrthoManager "Anasazi::BasicOrthoManager< ScalarType, MV, OP
>::BasicOrthoManager(Teuchos::RCP< const OP > Op=Teuchos::null,
typename Teuchos::ScalarTraits< ScalarType >::magnitudeType
kappa=1.41421356, typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType eps=0.0, typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType tol=0.20)

Constructor specifying re-orthogonalization tolerance. ";

%feature("docstring")  Anasazi::BasicOrthoManager::~BasicOrthoManager
"Anasazi::BasicOrthoManager< ScalarType, MV, OP
>::~BasicOrthoManager()

Destructor. ";

/*  Methods implementing Anasazi::MatOrthoManager  */

%feature("docstring")  Anasazi::BasicOrthoManager::projectMat "void
Anasazi::BasicOrthoManager< ScalarType, MV, OP >::projectMat(MV &X,
Teuchos::Array< Teuchos::RCP< const MV > > Q, Teuchos::Array<
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > > >
C=Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix< int,
ScalarType > >(Teuchos::null)), Teuchos::RCP< MV > MX=Teuchos::null,
Teuchos::Array< Teuchos::RCP< const MV > >
MQ=Teuchos::tuple(Teuchos::RCP< const MV >(Teuchos::null))) const

Given a list of mutually orthogonal and internally orthonormal bases
Q, this method projects a multivector X onto the space orthogonal to
the individual Q[i], optionally returning the coefficients of X for
the individual Q[i]. All of this is done with respect to the inner
product innerProd().

After calling this routine, X will be orthogonal to each of the Q[i].

Parameters:
-----------

X:  [in/out] The multivector to be modified.  On output, the columns
of X will be orthogonal to each Q[i], satisfying \\\\[ X_{out} =
X_{in} - \\\\sum_i Q[i] \\\\langle Q[i], X_{in} \\\\rangle \\\\]

MX:  [in/out] The image of X under the inner product operator Op. If $
MX != 0$: On input, this is expected to be consistent with Op  X. On
output, this is updated consistent with updates to X. If $ MX == 0$ or
$ Op == 0$: MX is not referenced.

C:  [out] The coefficients of X in the bases Q[i]. If C[i] is a non-
null pointer and C[i] matches the dimensions of X and Q[i], then the
coefficients computed during the orthogonalization routine will be
stored in the matrix C[i], similar to calling If C[i] points to a
Teuchos::SerialDenseMatrix with size inconsistent with X and  Q[i],
then a std::invalid_argument exception will be thrown. Otherwise, if
C.size() < i or C[i] is a null pointer, the caller will not have
access to the computed coefficients.

Q:  [in] A list of multivector bases specifying the subspaces to be
orthogonalized against, satisfying \\\\[ \\\\langle Q[i], Q[j]
\\\\rangle = I \\\\quad\\\\textrm{if}\\\\quad i=j \\\\] and \\\\[
\\\\langle Q[i], Q[j] \\\\rangle = 0 \\\\quad\\\\textrm{if}\\\\quad i
\\\\neq j\\\\ . \\\\] ";

%feature("docstring")  Anasazi::BasicOrthoManager::normalizeMat "int
Anasazi::BasicOrthoManager< ScalarType, MV, OP >::normalizeMat(MV &X,
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > >
B=Teuchos::null, Teuchos::RCP< MV > MX=Teuchos::null) const

This method takes a multivector X and attempts to compute an
orthonormal basis for $colspan(X)$, with respect to innerProd().

The method uses classical Gram-Schmidt with selective
reorthogonalization. As a result, the coefficient matrix B is upper
triangular.

This routine returns an integer rank stating the rank of the computed
basis. If X does not have full rank and the normalize() routine does
not attempt to augment the subspace, then rank may be smaller than the
number of columns in X. In this case, only the first rank columns of
output X and first rank rows of B will be valid.

The method attempts to find a basis with dimension equal to the number
of columns in X. It does this by augmenting linearly dependent vectors
in X with random directions. A finite number of these attempts will be
made; therefore, it is possible that the dimension of the computed
basis is less than the number of vectors in X.

Parameters:
-----------

X:  [in/out] The multivector to be modified.  On output, the first
rank columns of X satisfy \\\\[ \\\\langle X[i], X[j] \\\\rangle =
\\\\delta_{ij}\\\\ . \\\\] Also, \\\\[ X_{in}(1:m,1:n) =
X_{out}(1:m,1:rank) B(1:rank,1:n) \\\\] where m is the number of rows
in X and n is the number of columns in X.

MX:  [in/out] The image of X under the inner product operator Op. If $
MX != 0$: On input, this is expected to be consistent with Op  X. On
output, this is updated consistent with updates to X. If $ MX == 0$ or
$ Op == 0$: MX is not referenced.

B:  [out] The coefficients of the original X with respect to the
computed basis. If B is a non-null pointer and B matches the
dimensions of B, then the coefficients computed during the
orthogonalization routine will be stored in B, similar to calling If B
points to a Teuchos::SerialDenseMatrix with size inconsistent with X,
then a std::invalid_argument exception will be thrown. Otherwise, if B
is null, the caller will not have access to the computed coefficients.
This matrix is not necessarily triangular (as in a QR factorization);
see the documentation of specific orthogonalization managers.  The
first rows in B corresponding to the valid columns in X will be upper
triangular.

Rank of the basis computed by this method, less than or equal to the
number of columns in X. This specifies how many columns in the
returned X and rows in the returned B are valid. ";

%feature("docstring")
Anasazi::BasicOrthoManager::projectAndNormalizeMat "int
Anasazi::BasicOrthoManager< ScalarType, MV, OP
>::projectAndNormalizeMat(MV &X, Teuchos::Array< Teuchos::RCP< const
MV > > Q, Teuchos::Array< Teuchos::RCP< Teuchos::SerialDenseMatrix<
int, ScalarType > > > C=Teuchos::tuple(Teuchos::RCP<
Teuchos::SerialDenseMatrix< int, ScalarType > >(Teuchos::null)),
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > >
B=Teuchos::null, Teuchos::RCP< MV > MX=Teuchos::null, Teuchos::Array<
Teuchos::RCP< const MV > > MQ=Teuchos::tuple(Teuchos::RCP< const MV
>(Teuchos::null))) const

Given a set of bases Q[i] and a multivector X, this method computes an
orthonormal basis for $colspan(X) - \\\\sum_i colspan(Q[i])$.

This routine returns an integer rank stating the rank of the computed
basis. If the subspace $colspan(X) - \\\\sum_i colspan(Q[i])$ does not
have dimension as large as the number of columns of X and the
orthogonalization manager doe not attempt to augment the subspace,
then rank may be smaller than the number of columns of X. In this
case, only the first rank columns of output X and first rank rows of B
will be valid.

The method attempts to find a basis with dimension the same as the
number of columns in X. It does this by augmenting linearly dependent
vectors with random directions. A finite number of these attempts will
be made; therefore, it is possible that the dimension of the computed
basis is less than the number of vectors in X.

Parameters:
-----------

X:  [in/out] The multivector to be modified.  On output, the first
rank columns of X satisfy \\\\[ \\\\langle X[i], X[j] \\\\rangle =
\\\\delta_{ij} \\\\quad \\\\textrm{and} \\\\quad \\\\langle X, Q[i]
\\\\rangle = 0\\\\ . \\\\] Also, \\\\[ X_{in}(1:m,1:n) =
X_{out}(1:m,1:rank) B(1:rank,1:n) + \\\\sum_i Q[i] C[i] \\\\] where m
is the number of rows in X and n is the number of columns in X.

MX:  [in/out] The image of X under the inner product operator Op. If $
MX != 0$: On input, this is expected to be consistent with Op  X. On
output, this is updated consistent with updates to X. If $ MX == 0$ or
$ Op == 0$: MX is not referenced.

C:  [out] The coefficients of X in the Q[i]. If C[i] is a non-null
pointer and C[i] matches the dimensions of X and Q[i], then the
coefficients computed during the orthogonalization routine will be
stored in the matrix C[i], similar to calling If C[i] points to a
Teuchos::SerialDenseMatrix with size inconsistent with X and  Q[i],
then a std::invalid_argument exception will be thrown. Otherwise, if
C.size() < i or C[i] is a null pointer, the caller will not have
access to the computed coefficients.

B:  [out] The coefficients of the original X with respect to the
computed basis. If B is a non-null pointer and B matches the
dimensions of B, then the coefficients computed during the
orthogonalization routine will be stored in B, similar to calling If B
points to a Teuchos::SerialDenseMatrix with size inconsistent with X,
then a std::invalid_argument exception will be thrown. Otherwise, if B
is null, the caller will not have access to the computed coefficients.
This matrix is not necessarily triangular (as in a QR factorization);
see the documentation of specific orthogonalization managers.  The
first rows in B corresponding to the valid columns in X will be upper
triangular.

Q:  [in] A list of multivector bases specifying the subspaces to be
orthogonalized against, satisfying \\\\[ \\\\langle Q[i], Q[j]
\\\\rangle = I \\\\quad\\\\textrm{if}\\\\quad i=j \\\\] and \\\\[
\\\\langle Q[i], Q[j] \\\\rangle = 0 \\\\quad\\\\textrm{if}\\\\quad i
\\\\neq j\\\\ . \\\\]

Rank of the basis computed by this method, less than or equal to the
number of columns in X. This specifies how many columns in the
returned X and rows in the returned B are valid. ";

/*  Error methods  */

%feature("docstring")  Anasazi::BasicOrthoManager::orthonormErrorMat "Teuchos::ScalarTraits< ScalarType >::magnitudeType
Anasazi::BasicOrthoManager< ScalarType, MV, OP
>::orthonormErrorMat(const MV &X, Teuchos::RCP< const MV >
MX=Teuchos::null) const

This method computes the error in orthonormality of a multivector,
measured as the Frobenius norm of the difference innerProd(X,Y) - I.
The method has the option of exploiting a caller-provided MX. ";

%feature("docstring")  Anasazi::BasicOrthoManager::orthogErrorMat "Teuchos::ScalarTraits< ScalarType >::magnitudeType
Anasazi::BasicOrthoManager< ScalarType, MV, OP >::orthogErrorMat(const
MV &X1, const MV &X2, Teuchos::RCP< const MV > MX1, Teuchos::RCP<
const MV > MX2) const

This method computes the error in orthogonality of two multivectors,
measured as the Frobenius norm of innerProd(X,Y). The method has the
option of exploiting a caller-provided MX. ";

/*  Accessor routines  */

%feature("docstring")  Anasazi::BasicOrthoManager::setKappa "void
Anasazi::BasicOrthoManager< ScalarType, MV, OP >::setKappa(typename
Teuchos::ScalarTraits< ScalarType >::magnitudeType kappa)

Set parameter for re-orthogonalization threshold. ";

%feature("docstring")  Anasazi::BasicOrthoManager::getKappa "Teuchos::ScalarTraits<ScalarType>::magnitudeType
Anasazi::BasicOrthoManager< ScalarType, MV, OP >::getKappa() const

Return parameter for re-orthogonalization threshold. ";


// File: classAnasazi_1_1BasicOutputManager.xml
%feature("docstring") Anasazi::BasicOutputManager "

Anasazi's basic output manager for sending information of select
verbosity levels to the appropriate output stream.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziBasicOutputManager.hpp ";

/*  Constructors/Destructor  */

%feature("docstring")  Anasazi::BasicOutputManager::BasicOutputManager
"Anasazi::BasicOutputManager< ScalarType >::BasicOutputManager(int
vb=Anasazi::Errors, Teuchos::RCP< ostream >
os=Teuchos::rcpFromRef(std::cout), int printingRank=0)

Default constructor. ";

%feature("docstring")
Anasazi::BasicOutputManager::~BasicOutputManager "virtual
Anasazi::BasicOutputManager< ScalarType >::~BasicOutputManager()

Destructor. ";

/*  Set/Get methods  */

%feature("docstring")  Anasazi::BasicOutputManager::setOStream "void
Anasazi::BasicOutputManager< ScalarType >::setOStream(Teuchos::RCP<
ostream > os)

Set the output stream for this manager. ";

%feature("docstring")  Anasazi::BasicOutputManager::getOStream "Teuchos::RCP< ostream > Anasazi::BasicOutputManager< ScalarType
>::getOStream()

Get the output stream for this manager. ";

/*  Output methods  */

%feature("docstring")  Anasazi::BasicOutputManager::isVerbosity "bool
Anasazi::BasicOutputManager< ScalarType >::isVerbosity(MsgType type)
const

Find out whether we need to print out information for this message
type.

This method is used by the solver to determine whether computations
are necessary for this message type. ";

%feature("docstring")  Anasazi::BasicOutputManager::print "void
Anasazi::BasicOutputManager< ScalarType >::print(MsgType type, const
std::string output)

Send some output to this output stream. ";

%feature("docstring")  Anasazi::BasicOutputManager::stream "ostream &
Anasazi::BasicOutputManager< ScalarType >::stream(MsgType type)

Return a stream for outputting to. ";

/*  Undefined methods  */


// File: classAnasazi_1_1BasicSort.xml
%feature("docstring") Anasazi::BasicSort "

An implementation of the Anasazi::SortManager that performs a
collection of common sorting techniques.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziBasicSort.hpp ";

%feature("docstring")  Anasazi::BasicSort::BasicSort "Anasazi::BasicSort< MagnitudeType >::BasicSort(Teuchos::ParameterList
&pl)

Parameter list driven constructor.

This constructor accepts a paramter list with the following options:
\"Sort Strategy\" - a string specifying the desired sorting strategy.
See setSortType() for valid options. ";

%feature("docstring")  Anasazi::BasicSort::BasicSort "Anasazi::BasicSort< MagnitudeType >::BasicSort(const std::string
&which=\"LM\")

String driven constructor.

Directly pass the string specifying sort strategy. See setSortType()
for valid options. ";

%feature("docstring")  Anasazi::BasicSort::~BasicSort "Anasazi::BasicSort< MagnitudeType >::~BasicSort()

Destructor. ";

%feature("docstring")  Anasazi::BasicSort::setSortType "void
Anasazi::BasicSort< MagnitudeType >::setSortType(const std::string
&which)

Set sort type.

Parameters:
-----------

which:  [in] The eigenvalues of interest for this eigenproblem.
\"LM\" - Largest Magnitude [ default ]

\"SM\" - Smallest Magnitude

\"LR\" - Largest Real

\"SR\" - Smallest Real

\"LI\" - Largest Imaginary

\"SI\" - Smallest Imaginary ";

%feature("docstring")  Anasazi::BasicSort::sort "void
Anasazi::BasicSort< MagnitudeType >::sort(std::vector< MagnitudeType >
&evals, Teuchos::RCP< std::vector< int > > perm=Teuchos::null, int
n=-1) const

Sort real eigenvalues, optionally returning the permutation vector.

This method is not valid when the sort manager is configured for
\"LI\" or \"SI\" sorting (i.e., sorting by the imaginary components).
Calling this method in that scenario will result in a SortManagerError
exception.

Parameters:
-----------

evals:  [in/out] Vector of length at least n containing the
eigenvalues to be sorted.  On output, the first n eigenvalues will be
sorted. The rest will be unchanged.

perm:  [out] Vector of length at least n to store the permutation
index (optional).  If specified, on output the first n eigenvalues
will contain the permutation indices, in the range [0,n-1], such that
evals_out[i] = evals_in[perm[i]]

n:  [in] Number of values in evals to be sorted. If n == -1, all
values will be sorted. ";

%feature("docstring")  Anasazi::BasicSort::sort "void
Anasazi::BasicSort< MagnitudeType >::sort(std::vector< MagnitudeType >
&r_evals, std::vector< MagnitudeType > &i_evals, Teuchos::RCP<
std::vector< int > > perm=Teuchos::null, int n=-1) const

Sort complex eigenvalues, optionally returning the permutation vector.

This routine takes two vectors, one for each part of a complex
eigenvalue. This is helpful for solving real, non-symmetric eigenvalue
problems.

Parameters:
-----------

r_evals:  [in/out] Vector of length at least n containing the real
part of the eigenvalues to be sorted.  On output, the first n
eigenvalues will be sorted. The rest will be unchanged.

i_evals:  [in/out] Vector of length at least n containing the
imaginary part of the eigenvalues to be sorted.  On output, the first
n eigenvalues will be sorted. The rest will be unchanged.

perm:  [out] Vector of length at least n to store the permutation
index (optional).  If specified, on output the first n eigenvalues
will contain the permutation indices, in the range [0,n-1], such that
r_evals_out[i] = r_evals_in[perm[i]] and similarly for i_evals.

n:  [in] Number of values in r_evals, i_evals to be sorted. If n ==
-1, all values will be sorted, as decided by the minimum of the length
of r_evals and the length of i_evals. ";


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

/*  Constructor/Destructor  */

%feature("docstring")  Anasazi::BlockDavidson::BlockDavidson "Anasazi::BlockDavidson< ScalarType, MV, OP >::BlockDavidson(const
Teuchos::RCP< Eigenproblem< ScalarType, MV, OP > > &problem, const
Teuchos::RCP< SortManager< typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType > > &sorter, const Teuchos::RCP< OutputManager<
ScalarType > > &printer, const Teuchos::RCP< StatusTest< ScalarType,
MV, OP > > &tester, const Teuchos::RCP< MatOrthoManager< ScalarType,
MV, OP > > &ortho, Teuchos::ParameterList &params)

BlockDavidson constructor with eigenproblem, solver utilities, and
parameter list of solver options.

This constructor takes pointers required by the eigensolver, in
addition to a parameter list of options for the eigensolver. These
options include the following: \"Block Size\" - an int specifying the
block size used by the algorithm. This can also be specified using the
setBlockSize() method.

\"Num Blocks\" - an int specifying the maximum number of blocks
allocated for the solver basis. ";

%feature("docstring")  Anasazi::BlockDavidson::~BlockDavidson "Anasazi::BlockDavidson< ScalarType, MV, OP >::~BlockDavidson()

Anasazi::BlockDavidson destructor. ";

/*  Solver methods  */

%feature("docstring")  Anasazi::BlockDavidson::iterate "void
Anasazi::BlockDavidson< ScalarType, MV, OP >::iterate()

This method performs BlockDavidson iterations until the status test
indicates the need to stop or an error occurs (in which case, an
appropriate exception is thrown).

iterate() will first determine whether the solver is uninitialized; if
not, it will call initialize(). After initialization, the solver
performs block Davidson iterations until the status test evaluates as
::Passed, at which point the method returns to the caller.

The block Davidson iteration proceeds as follows: The current residual
(R) is preconditioned to form H

H is orthogonalized against the auxiliary vectors and the previous
basis vectors, and made orthonormal.

The current basis is expanded with H and used to project the problem
matrix.

The projected eigenproblem is solved, and the desired eigenvectors and
eigenvalues are selected.

These are used to form the new eigenvector estimates (X).

The new residual (R) is formed.

The status test is queried at the beginning of the iteration.

Possible exceptions thrown include std::invalid_argument or one of the
BlockDavidson-specific exceptions. ";

%feature("docstring")  Anasazi::BlockDavidson::initialize "void
Anasazi::BlockDavidson< ScalarType, MV, OP
>::initialize(BlockDavidsonState< ScalarType, MV > newstate)

Initialize the solver to an iterate, optionally providing the current
basis and projected problem matrix, the current Ritz vectors and
values, and the current residual.

The BlockDavidson eigensolver contains a certain amount of state,
including the current Krylov basis, the current eigenvectors, the
current residual, etc. (see getState())

initialize() gives the user the opportunity to manually set these,
although this must be done with caution, as the validity of the user
input will not be checked.

Only the first newstate.curDim columns of newstate.V and newstate.KK
and the first newstate.curDim rows of newstate.KK will be used.

If newstate.V == getState().V, then the data is not copied. The same
holds for newstate.KK, newstate.X, newstate.KX, newstate.MX, and
newstate.R Only the upper triangular half of newstate.KK is used to
initialize the state of the solver.

isInitialized() == true (see post-conditions of isInitialize())  The
user has the option of specifying any component of the state using
initialize(). However, these arguments are assumed to match the post-
conditions specified under isInitialized(). Any component of the state
(i.e., KX) not given to initialize() will be generated.

Note, for any pointer in newstate which directly points to the
multivectors in the solver, the data is not copied. ";

%feature("docstring")  Anasazi::BlockDavidson::initialize "void
Anasazi::BlockDavidson< ScalarType, MV, OP >::initialize()

Initialize the solver with the initial vectors from the eigenproblem
or random data. ";

%feature("docstring")  Anasazi::BlockDavidson::isInitialized "bool
Anasazi::BlockDavidson< ScalarType, MV, OP >::isInitialized() const

Indicates whether the solver has been initialized or not.

bool indicating the state of the solver.

If isInitialized() == true:  getCurSubspaceDim() > 0 and is a multiple
of getBlockSize()

the first getCurSubspaceDim() vectors of V are orthogonal to auxiliary
vectors and have orthonormal columns

the principal submatrix of order getCurSubspaceDim() of KK contains
the project eigenproblem matrix

X contains the Ritz vectors with respect to the current Krylov basis

T contains the Ritz values with respect to the current Krylov basis

KX == Op*X

MX == M*X if M != Teuchos::null  Otherwise, MX == Teuchos::null

R contains the residual vectors with respect to X ";

%feature("docstring")  Anasazi::BlockDavidson::getState "BlockDavidsonState< ScalarType, MV > Anasazi::BlockDavidson<
ScalarType, MV, OP >::getState() const

Get access to the current state of the eigensolver.

The data is only valid if isInitialized() == true.

The data for the preconditioned residual is only meaningful in the
scenario that the solver throws a ::BlockDavidsonRitzFailure exception
during iterate().

A BlockDavidsonState object containing const pointers to the current
solver state. Note, these are direct pointers to the multivectors;
they are not pointers to views of the multivectors. ";

/*  Status methods  */

%feature("docstring")  Anasazi::BlockDavidson::getNumIters "int
Anasazi::BlockDavidson< ScalarType, MV, OP >::getNumIters() const

Get the current iteration count. ";

%feature("docstring")  Anasazi::BlockDavidson::resetNumIters "void
Anasazi::BlockDavidson< ScalarType, MV, OP >::resetNumIters()

Reset the iteration count. ";

%feature("docstring")  Anasazi::BlockDavidson::getRitzVectors "Teuchos::RCP< const MV > Anasazi::BlockDavidson< ScalarType, MV, OP
>::getRitzVectors()

Get access to the current Ritz vectors.

A multivector with getBlockSize() vectors containing the sorted Ritz
vectors corresponding to the most significant Ritz values. The i-th
vector of the return corresponds to the i-th Ritz vector; there is no
need to use getRitzIndex(). ";

%feature("docstring")  Anasazi::BlockDavidson::getRitzValues "std::vector< Value< ScalarType > > Anasazi::BlockDavidson< ScalarType,
MV, OP >::getRitzValues()

Get the Ritz values for the previous iteration.

A vector of length getCurSubspaceDim() containing the Ritz values from
the previous projected eigensolve. ";

%feature("docstring")  Anasazi::BlockDavidson::getRitzIndex "std::vector< int > Anasazi::BlockDavidson< ScalarType, MV, OP
>::getRitzIndex()

Get the index used for extracting individual Ritz vectors from
getRitzVectors().

Because BlockDavidson is a Hermitian solver, all Ritz values are real
and all Ritz vectors can be represented in a single column of a
multivector. Therefore, getRitzIndex() is not needed when using the
output from getRitzVectors().

An int vector of size getCurSubspaceDim() composed of zeros. ";

%feature("docstring")  Anasazi::BlockDavidson::getResNorms "std::vector< typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType > Anasazi::BlockDavidson< ScalarType, MV, OP
>::getResNorms()

Get the current residual norms, computing the norms if they are not
up-to-date with the current residual vectors.

A vector of length getCurSubspaceDim() containing the norms of the
residuals, with respect to the orthogonalization manager's norm()
method. ";

%feature("docstring")  Anasazi::BlockDavidson::getRes2Norms "std::vector< typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType > Anasazi::BlockDavidson< ScalarType, MV, OP
>::getRes2Norms()

Get the current residual 2-norms, computing the norms if they are not
up-to-date with the current residual vectors.

A vector of length getCurSubspaceDim() containing the 2-norms of the
current residuals. ";

%feature("docstring")  Anasazi::BlockDavidson::getRitzRes2Norms "std::vector< typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType > Anasazi::BlockDavidson< ScalarType, MV, OP
>::getRitzRes2Norms()

Get the 2-norms of the residuals.

The Ritz residuals are not defined for the LOBPCG iteration. Hence,
this method returns the 2-norms of the direct residuals, and is
equivalent to calling getRes2Norms().

A vector of length getBlockSize() containing the 2-norms of the direct
residuals. ";

%feature("docstring")  Anasazi::BlockDavidson::getCurSubspaceDim "int
Anasazi::BlockDavidson< ScalarType, MV, OP >::getCurSubspaceDim()
const

Get the dimension of the search subspace used to generate the current
eigenvectors and eigenvalues.

An integer specifying the rank of the Krylov subspace currently in use
by the eigensolver. If isInitialized() == false, the return is 0.
Otherwise, it will be some strictly positive multiple of
getBlockSize(). ";

%feature("docstring")  Anasazi::BlockDavidson::getMaxSubspaceDim "int
Anasazi::BlockDavidson< ScalarType, MV, OP >::getMaxSubspaceDim()
const

Get the maximum dimension allocated for the search subspace. For
BlockDavidson, this always returns numBlocks*blockSize. ";

/*  Accessor routines from Eigensolver  */

%feature("docstring")  Anasazi::BlockDavidson::setStatusTest "void
Anasazi::BlockDavidson< ScalarType, MV, OP
>::setStatusTest(Teuchos::RCP< StatusTest< ScalarType, MV, OP > >
test)

Set a new StatusTest for the solver. ";

%feature("docstring")  Anasazi::BlockDavidson::getStatusTest "Teuchos::RCP< StatusTest< ScalarType, MV, OP > >
Anasazi::BlockDavidson< ScalarType, MV, OP >::getStatusTest() const

Get the current StatusTest used by the solver. ";

%feature("docstring")  Anasazi::BlockDavidson::getProblem "const
Eigenproblem< ScalarType, MV, OP > & Anasazi::BlockDavidson<
ScalarType, MV, OP >::getProblem() const

Get a constant reference to the eigenvalue problem. ";

%feature("docstring")  Anasazi::BlockDavidson::setBlockSize "void
Anasazi::BlockDavidson< ScalarType, MV, OP >::setBlockSize(int
blockSize)

Set the blocksize.

This method is required to support the interface provided by
Eigensolver. However, the preferred method of setting the allocated
size for the BlockDavidson eigensolver is setSize(). In fact,
setBlockSize() simply calls setSize(), maintaining the current number
of blocks.

The block size determines the number of Ritz vectors and values that
are computed on each iteration, thereby determining the increase in
the Krylov subspace at each iteration. ";

%feature("docstring")  Anasazi::BlockDavidson::getBlockSize "int
Anasazi::BlockDavidson< ScalarType, MV, OP >::getBlockSize() const

Get the blocksize used by the iterative solver. ";

%feature("docstring")  Anasazi::BlockDavidson::setAuxVecs "void
Anasazi::BlockDavidson< ScalarType, MV, OP >::setAuxVecs(const
Teuchos::Array< Teuchos::RCP< const MV > > &auxvecs)

Set the auxiliary vectors for the solver.

Because the current basis V cannot be assumed orthogonal to the new
auxiliary vectors, a call to setAuxVecs() will reset the solver to the
uninitialized state. This happens only in the case where the new
auxiliary vectors have a combined dimension of greater than zero.

In order to preserve the current state, the user will need to extract
it from the solver using getState(), orthogonalize it against the new
auxiliary vectors, and reinitialize using initialize(). ";

%feature("docstring")  Anasazi::BlockDavidson::getAuxVecs "Teuchos::Array< Teuchos::RCP< const MV > > Anasazi::BlockDavidson<
ScalarType, MV, OP >::getAuxVecs() const

Get the auxiliary vectors for the solver. ";

/*  BlockDavidson-specific accessor routines  */

%feature("docstring")  Anasazi::BlockDavidson::setSize "void
Anasazi::BlockDavidson< ScalarType, MV, OP >::setSize(int blockSize,
int numBlocks)

Set the blocksize and number of blocks to be used by the iterative
solver in solving this eigenproblem.

Changing either the block size or the number of blocks will reset the
solver to an uninitialized state.

The requested block size must be strictly positive; the number of
blocks must be greater than one. Invalid arguments will result in a
std::invalid_argument exception. ";

/*  Output methods  */

%feature("docstring")  Anasazi::BlockDavidson::currentStatus "void
Anasazi::BlockDavidson< ScalarType, MV, OP
>::currentStatus(std::ostream &os)

This method requests that the solver print out its current status to
the given output stream. ";


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
return ::Failed; if it returns ::Passed, solve() will throw an
AnasaziError exception.

Additionally, the solver manager will terminate solve() after a
specified number of restarts.

Much of this behavior is controlled via parameters and options passed
to the solver manager. For more information, see
BlockDavidsonSolMgr().

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist

C++ includes: AnasaziBlockDavidsonSolMgr.hpp ";

/*  Constructors/Destructor  */

%feature("docstring")
Anasazi::BlockDavidsonSolMgr::BlockDavidsonSolMgr "Anasazi::BlockDavidsonSolMgr< ScalarType, MV, OP
>::BlockDavidsonSolMgr(const Teuchos::RCP< Eigenproblem< ScalarType,
MV, OP > > &problem, Teuchos::ParameterList &pl)

Basic constructor for BlockDavidsonSolMgr.

This constructor accepts the Eigenproblem to be solved in addition to
a parameter list of options for the solver manager. These options
include the following: Solver parameters  \"Which\" - a string
specifying the desired eigenvalues: SM, LM, SR or LR. Default: \"SR\"

\"Block Size\" - a int specifying the block size to be used by the
underlying block Davidson solver. Default: problem->getNEV()

\"Num Blocks\" - a int specifying the number of blocks allocated for
the Krylov basis. Default: 2

\"Maximum Restarts\" - a int specifying the maximum number of restarts
the underlying solver is allowed to perform. Default: 20

\"Verbosity\" - a sum of MsgType specifying the verbosity. Default:
::Errors

Convergence parameters (if using default convergence test; see
setGlobalStatusTest())  \"Convergence Tolerance\" - a MagnitudeType
specifying the level that residual norms must reach to decide
convergence. Default: machine precision.

\"Relative Convergence Tolerance\" - a bool specifying whether
residuals norms should be scaled by their eigenvalues for the
purposing of deciding convergence. Default: true

\"Convergence Norm\" - a string specifying the norm for convergence
testing: \"2\" or \"M\"

Locking parameters (if using default locking test; see
setLockingStatusTest())  \"Use Locking\" - a bool specifying whether
the algorithm should employ locking of converged eigenpairs. Default:
false

\"Max Locked\" - a int specifying the maximum number of eigenpairs to
be locked. Default: problem->getNEV()

\"Locking Quorum\" - a int specifying the number of eigenpairs that
must meet the locking criteria before locking actually occurs.
Default: 1

\"Locking Tolerance\" - a MagnitudeType specifying the level that
residual norms must reach to decide locking. Default: 0.1*convergence
tolerance

\"Relative Locking Tolerance\" - a bool specifying whether residuals
norms should be scaled by their eigenvalues for the purposing of
deciding locking. Default: true

\"Locking Norm\" - a string specifying the norm for locking testing:
\"2\" or \"M\" ";

%feature("docstring")
Anasazi::BlockDavidsonSolMgr::~BlockDavidsonSolMgr "virtual
Anasazi::BlockDavidsonSolMgr< ScalarType, MV, OP
>::~BlockDavidsonSolMgr()

Destructor. ";

/*  Accessor methods  */

%feature("docstring")  Anasazi::BlockDavidsonSolMgr::getProblem "const Eigenproblem<ScalarType,MV,OP>& Anasazi::BlockDavidsonSolMgr<
ScalarType, MV, OP >::getProblem() const

Return the eigenvalue problem. ";

%feature("docstring")  Anasazi::BlockDavidsonSolMgr::getNumIters "int
Anasazi::BlockDavidsonSolMgr< ScalarType, MV, OP >::getNumIters()
const

Get the iteration count for the most recent call to  solve(). ";

%feature("docstring")  Anasazi::BlockDavidsonSolMgr::getTimers "Teuchos::Array<Teuchos::RCP<Teuchos::Time> >
Anasazi::BlockDavidsonSolMgr< ScalarType, MV, OP >::getTimers() const

Return the timers for this object.

The timers are ordered as follows: time spent in solve() routine

time spent restarting

time spent locking converged eigenvectors ";

/*  Solver application methods  */

%feature("docstring")  Anasazi::BlockDavidsonSolMgr::solve "ReturnType Anasazi::BlockDavidsonSolMgr< ScalarType, MV, OP >::solve()

This method performs possibly repeated calls to the underlying
eigensolver's iterate() routine until the problem has been solved (as
decided by the solver manager) or the solver manager decides to quit.

This method calls BlockDavidson::iterate(), which will return either
because a specially constructed status test evaluates to ::Passed or
an exception is thrown.

A return from BlockDavidson::iterate() signifies one of the following
scenarios: the maximum number of restarts has been exceeded. In this
scenario, the solver manager will place  all converged eigenpairs into
the eigenproblem and return ::Unconverged.

the locking conditions have been met. In this scenario, some of the
current eigenpairs will be removed  from the eigensolver and placed
into auxiliary storage. The eigensolver will be restarted with the
remaining part of the Krylov subspace  and some random information to
replace the removed subspace.

global convergence has been met. In this case, the most significant
NEV eigenpairs in the solver and locked storage  have met the
convergence criterion. (Here, NEV refers to the number of eigenpairs
requested by the Eigenproblem.)  In this scenario, the solver manager
will return ::Converged.

::ReturnType specifying: ::Converged: the eigenproblem was solved to
the specification required by the solver manager.

::Unconverged: the eigenproblem was not solved to the specification
desired by the solver manager. ";

%feature("docstring")
Anasazi::BlockDavidsonSolMgr::setGlobalStatusTest "void
Anasazi::BlockDavidsonSolMgr< ScalarType, MV, OP
>::setGlobalStatusTest(const Teuchos::RCP< StatusTest< ScalarType, MV,
OP > > &global)

Set the status test defining global convergence. ";

%feature("docstring")
Anasazi::BlockDavidsonSolMgr::getGlobalStatusTest "const
Teuchos::RCP< StatusTest< ScalarType, MV, OP > > &
Anasazi::BlockDavidsonSolMgr< ScalarType, MV, OP
>::getGlobalStatusTest() const

Get the status test defining global convergence. ";

%feature("docstring")
Anasazi::BlockDavidsonSolMgr::setLockingStatusTest "void
Anasazi::BlockDavidsonSolMgr< ScalarType, MV, OP
>::setLockingStatusTest(const Teuchos::RCP< StatusTest< ScalarType,
MV, OP > > &locking)

Set the status test defining locking. ";

%feature("docstring")
Anasazi::BlockDavidsonSolMgr::getLockingStatusTest "const
Teuchos::RCP< StatusTest< ScalarType, MV, OP > > &
Anasazi::BlockDavidsonSolMgr< ScalarType, MV, OP
>::getLockingStatusTest() const

Get the status test defining locking. ";

%feature("docstring")
Anasazi::BlockDavidsonSolMgr::setDebugStatusTest "void
Anasazi::BlockDavidsonSolMgr< ScalarType, MV, OP
>::setDebugStatusTest(const Teuchos::RCP< StatusTest< ScalarType, MV,
OP > > &debug)

Set the status test for debugging. ";

%feature("docstring")
Anasazi::BlockDavidsonSolMgr::getDebugStatusTest "const Teuchos::RCP<
StatusTest< ScalarType, MV, OP > > & Anasazi::BlockDavidsonSolMgr<
ScalarType, MV, OP >::getDebugStatusTest() const

Get the status test for debugging. ";


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

/*  Constructor/Destructor  */

%feature("docstring")  Anasazi::BlockKrylovSchur::BlockKrylovSchur "Anasazi::BlockKrylovSchur< ScalarType, MV, OP
>::BlockKrylovSchur(const Teuchos::RCP< Eigenproblem< ScalarType, MV,
OP > > &problem, const Teuchos::RCP< SortManager< typename
Teuchos::ScalarTraits< ScalarType >::magnitudeType > > &sorter, const
Teuchos::RCP< OutputManager< ScalarType > > &printer, const
Teuchos::RCP< StatusTest< ScalarType, MV, OP > > &tester, const
Teuchos::RCP< OrthoManager< ScalarType, MV > > &ortho,
Teuchos::ParameterList &params)

BlockKrylovSchur constructor with eigenproblem, solver utilities, and
parameter list of solver options.

This constructor takes pointers required by the eigensolver, in
addition to a parameter list of options for the eigensolver. These
options include the following: \"Block Size\" - an int specifying the
block size used by the algorithm. This can also be specified using the
setBlockSize() method. Default: 1

\"Num Blocks\" - an int specifying the maximum number of blocks
allocated for the solver basis. Default: 3*problem->getNEV()

\"Step Size\" - an int specifying how many iterations are performed
between computations of eigenvalues and eigenvectors.  Note: This
parameter is mandatory.

\"Number of Ritz Vectors\" - an int specifying how many Ritz vectors
are computed on calls to getRitzVectors(). Default: 0

\"Print Number of Ritz Values\" - an int specifying how many Ritz
values are printed on calls to currentStatus(). Default: \"Block
Size\" ";

%feature("docstring")  Anasazi::BlockKrylovSchur::~BlockKrylovSchur "virtual Anasazi::BlockKrylovSchur< ScalarType, MV, OP
>::~BlockKrylovSchur()

BlockKrylovSchur destructor. ";

/*  Solver methods  */

%feature("docstring")  Anasazi::BlockKrylovSchur::iterate "void
Anasazi::BlockKrylovSchur< ScalarType, MV, OP >::iterate()

This method performs Block Krylov-Schur iterations until the status
test indicates the need to stop or an error occurs (in which case, an
exception is thrown).

iterate() will first determine whether the solver is inintialized; if
not, it will call initialize() using default arguments. After
initialization, the solver performs Block Krylov-Schur iterations
until the status test evaluates as ::Passed, at which point the method
returns to the caller.

The Block Krylov-Schur iteration proceeds as follows: The operator
problem->getOperator() is applied to the newest blockSize vectors in
the Krylov basis.

The resulting vectors are orthogonalized against the auxiliary vectors
and the previous basis vectors, and made orthonormal.

The Hessenberg matrix is updated.

If we have performed stepSize iterations since the last update, update
the Ritz values and Ritz residuals.

The status test is queried at the beginning of the iteration.

Possible exceptions thrown include the BlockKrylovSchurOrthoFailure.
";

%feature("docstring")  Anasazi::BlockKrylovSchur::initialize "void
Anasazi::BlockKrylovSchur< ScalarType, MV, OP
>::initialize(BlockKrylovSchurState< ScalarType, MV > state)

Initialize the solver to an iterate, providing a Krylov basis and
Hessenberg matrix.

The BlockKrylovSchur eigensolver contains a certain amount of state,
consisting of the current Krylov basis and the associated Hessenberg
matrix.

initialize() gives the user the opportunity to manually set these,
although this must be done with caution, abiding by the rules given
below. All notions of orthogonality and orthonormality are derived
from the inner product specified by the orthogonalization manager.

isInitialized() == true (see post-conditions of isInitialize())  The
user has the option of specifying any component of the state using
initialize(). However, these arguments are assumed to match the post-
conditions specified under isInitialized(). Any necessary component of
the state not given to initialize() will be generated.

Note, for any pointer in newstate which directly points to the
multivectors in the solver, the data is not copied. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::initialize "void
Anasazi::BlockKrylovSchur< ScalarType, MV, OP >::initialize()

Initialize the solver with the initial vectors from the eigenproblem
or random data. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::isInitialized "bool
Anasazi::BlockKrylovSchur< ScalarType, MV, OP >::isInitialized() const

Indicates whether the solver has been initialized or not.

bool indicating the state of the solver.

If isInitialized() == true: the first getCurSubspaceDim() vectors of V
are orthogonal to auxiliary vectors and have orthonormal columns

the principal Hessenberg submatrix of of H contains the Hessenberg
matrix associated with V ";

%feature("docstring")  Anasazi::BlockKrylovSchur::getState "BlockKrylovSchurState<ScalarType,MV> Anasazi::BlockKrylovSchur<
ScalarType, MV, OP >::getState() const

Get the current state of the eigensolver.

The data is only valid if isInitialized() == true.

A BlockKrylovSchurState object containing const pointers to the
current solver state. ";

/*  Status methods  */

%feature("docstring")  Anasazi::BlockKrylovSchur::getNumIters "int
Anasazi::BlockKrylovSchur< ScalarType, MV, OP >::getNumIters() const

Get the current iteration count. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::resetNumIters "void
Anasazi::BlockKrylovSchur< ScalarType, MV, OP >::resetNumIters()

Reset the iteration count. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::getRitzVectors "Teuchos::RCP<const MV> Anasazi::BlockKrylovSchur< ScalarType, MV, OP
>::getRitzVectors()

Get the Ritz vectors.

A multivector of columns not exceeding the maximum dimension of the
subspace containing the Ritz vectors from the most recent call to
computeRitzVectors().

To see if the returned Ritz vectors are current, call
isRitzVecsCurrent(). ";

%feature("docstring")  Anasazi::BlockKrylovSchur::getRitzValues "std::vector<Value<ScalarType> > Anasazi::BlockKrylovSchur< ScalarType,
MV, OP >::getRitzValues()

Get the Ritz values.

A vector of length not exceeding the maximum dimension of the subspace
containing the Ritz values from the most recent Schur form update.

To see if the returned Ritz values are current, call
isRitzValsCurrent(). ";

%feature("docstring")  Anasazi::BlockKrylovSchur::getRitzIndex "std::vector<int> Anasazi::BlockKrylovSchur< ScalarType, MV, OP
>::getRitzIndex()

Get the Ritz index vector.

A vector of length not exceeding the maximum dimension of the subspace
containing the index vector for the Ritz values and Ritz vectors, if
they are computed. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::getResNorms "std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>
Anasazi::BlockKrylovSchur< ScalarType, MV, OP >::getResNorms()

Get the current residual norms.

Block Krylov-Schur cannot provide this so a zero length vector will be
returned. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::getRes2Norms "std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>
Anasazi::BlockKrylovSchur< ScalarType, MV, OP >::getRes2Norms()

Get the current residual 2-norms.

Block Krylov-Schur cannot provide this so a zero length vector will be
returned. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::getRitzRes2Norms "std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>
Anasazi::BlockKrylovSchur< ScalarType, MV, OP >::getRitzRes2Norms()

Get the current Ritz residual 2-norms.

A vector of length blockSize containing the 2-norms of the Ritz
residuals. ";

/*  Accessor routines  */

%feature("docstring")  Anasazi::BlockKrylovSchur::setStatusTest "void
Anasazi::BlockKrylovSchur< ScalarType, MV, OP
>::setStatusTest(Teuchos::RCP< StatusTest< ScalarType, MV, OP > >
test)

Set a new StatusTest for the solver. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::getStatusTest "Teuchos::RCP< StatusTest< ScalarType, MV, OP > >
Anasazi::BlockKrylovSchur< ScalarType, MV, OP >::getStatusTest() const

Get the current StatusTest used by the solver. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::getProblem "const
Eigenproblem<ScalarType,MV,OP>& Anasazi::BlockKrylovSchur< ScalarType,
MV, OP >::getProblem() const

Get a constant reference to the eigenvalue problem. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::setSize "void
Anasazi::BlockKrylovSchur< ScalarType, MV, OP >::setSize(int
blockSize, int numBlocks)

Set the blocksize and number of blocks to be used by the iterative
solver in solving this eigenproblem.

Changing either the block size or the number of blocks will reset the
solver to an uninitialized state. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::setBlockSize "void
Anasazi::BlockKrylovSchur< ScalarType, MV, OP >::setBlockSize(int
blockSize)

Set the blocksize. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::setStepSize "void
Anasazi::BlockKrylovSchur< ScalarType, MV, OP >::setStepSize(int
stepSize)

Set the step size. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::setNumRitzVectors "void Anasazi::BlockKrylovSchur< ScalarType, MV, OP
>::setNumRitzVectors(int numRitzVecs)

Set the number of Ritz vectors to compute. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::getStepSize "int
Anasazi::BlockKrylovSchur< ScalarType, MV, OP >::getStepSize() const

Get the step size. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::getBlockSize "int
Anasazi::BlockKrylovSchur< ScalarType, MV, OP >::getBlockSize() const

Get the blocksize to be used by the iterative solver in solving this
eigenproblem. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::getNumRitzVectors "int Anasazi::BlockKrylovSchur< ScalarType, MV, OP
>::getNumRitzVectors() const

Get the number of Ritz vectors to compute. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::getCurSubspaceDim "int Anasazi::BlockKrylovSchur< ScalarType, MV, OP
>::getCurSubspaceDim() const

Get the dimension of the search subspace used to generate the current
eigenvectors and eigenvalues.

An integer specifying the rank of the Krylov subspace currently in use
by the eigensolver. If isInitialized() == false, the return is 0. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::getMaxSubspaceDim "int Anasazi::BlockKrylovSchur< ScalarType, MV, OP
>::getMaxSubspaceDim() const

Get the maximum dimension allocated for the search subspace. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::setAuxVecs "void
Anasazi::BlockKrylovSchur< ScalarType, MV, OP >::setAuxVecs(const
Teuchos::Array< Teuchos::RCP< const MV > > &auxvecs)

Set the auxiliary vectors for the solver.

Because the current Krylov subspace cannot be assumed orthogonal to
the new auxiliary vectors, a call to setAuxVecs() will reset the
solver to the uninitialized state. This happens only in the case where
the new auxiliary vectors have a combined dimension of greater than
zero.

In order to preserve the current state, the user will need to extract
it from the solver using getState(), orthogonalize it against the new
auxiliary vectors, and reinitialize using initialize(). ";

%feature("docstring")  Anasazi::BlockKrylovSchur::getAuxVecs "Teuchos::Array<Teuchos::RCP<const MV> > Anasazi::BlockKrylovSchur<
ScalarType, MV, OP >::getAuxVecs() const

Get the auxiliary vectors for the solver. ";

/*  Output methods  */

%feature("docstring")  Anasazi::BlockKrylovSchur::currentStatus "void
Anasazi::BlockKrylovSchur< ScalarType, MV, OP
>::currentStatus(std::ostream &os)

This method requests that the solver print out its current status to
screen. ";

/*  Block-Krylov Schur status routines  */

%feature("docstring")  Anasazi::BlockKrylovSchur::isRitzVecsCurrent "bool Anasazi::BlockKrylovSchur< ScalarType, MV, OP
>::isRitzVecsCurrent() const

Get the status of the Ritz vectors currently stored in the
eigensolver. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::isRitzValsCurrent "bool Anasazi::BlockKrylovSchur< ScalarType, MV, OP
>::isRitzValsCurrent() const

Get the status of the Ritz values currently stored in the eigensolver.
";

%feature("docstring")  Anasazi::BlockKrylovSchur::isSchurCurrent "bool Anasazi::BlockKrylovSchur< ScalarType, MV, OP >::isSchurCurrent()
const

Get the status of the Schur form currently stored in the eigensolver.
";

/*  Block-Krylov Schur compute routines  */

%feature("docstring")  Anasazi::BlockKrylovSchur::computeRitzVectors "void Anasazi::BlockKrylovSchur< ScalarType, MV, OP
>::computeRitzVectors()

Compute the Ritz vectors using the current Krylov factorization. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::computeRitzValues "void Anasazi::BlockKrylovSchur< ScalarType, MV, OP
>::computeRitzValues()

Compute the Ritz values using the current Krylov factorization. ";

%feature("docstring")  Anasazi::BlockKrylovSchur::computeSchurForm "void Anasazi::BlockKrylovSchur< ScalarType, MV, OP
>::computeSchurForm(const bool sort=true)

Compute the Schur form of the projected eigenproblem from the current
Krylov factorization. ";


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
return ::Failed; if it returns ::Passed, solve() will throw an
AnasaziError exception.

Additionally, the solver manager will terminate solve() after a
specified number of restarts.

Much of this behavior is controlled via parameters and options passed
to the solver manager. For more information, see
BlockKrylovSchurSolMgr().

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist

C++ includes: AnasaziBlockKrylovSchurSolMgr.hpp ";

/*  Constructors/Destructor  */

%feature("docstring")
Anasazi::BlockKrylovSchurSolMgr::BlockKrylovSchurSolMgr "Anasazi::BlockKrylovSchurSolMgr< ScalarType, MV, OP
>::BlockKrylovSchurSolMgr(const Teuchos::RCP< Eigenproblem<
ScalarType, MV, OP > > &problem, Teuchos::ParameterList &pl)

Basic constructor for BlockKrylovSchurSolMgr.

This constructor accepts the Eigenproblem to be solved in addition to
a parameter list of options for the solver manager. These options
include the following: Solver parameters  \"Which\" - a string
specifying the desired eigenvalues: SM, LM, SR or LR. Default: \"LM\"

\"Block Size\" - a int specifying the block size to be used by the
underlying block Krylov- Schur solver. Default: 1

\"Num Blocks\" - a int specifying the number of blocks allocated for
the Krylov basis. Default: 3*nev

\"Extra NEV Blocks\" - a int specifying the number of extra blocks the
solver should keep in addition to those required to compute the number
of eigenvalues requested. Default: 0

\"Maximum Restarts\" - a int specifying the maximum number of restarts
the underlying solver is allowed to perform. Default: 20

\"Orthogonalization\" - a string specifying the desired
orthogonalization: DGKS and SVQB. Default: \"SVQB\"

\"Verbosity\" - a sum of MsgType specifying the verbosity. Default:
Anasazi::Errors

Convergence parameters  \"Convergence Tolerance\" - a MagnitudeType
specifying the level that residual norms must reach to decide
convergence. Default: machine precision.

\"Relative Convergence Tolerance\" - a bool specifying whether
residuals norms should be scaled by their eigenvalues for the
purposing of deciding convergence. Default: true ";

%feature("docstring")
Anasazi::BlockKrylovSchurSolMgr::~BlockKrylovSchurSolMgr "virtual
Anasazi::BlockKrylovSchurSolMgr< ScalarType, MV, OP
>::~BlockKrylovSchurSolMgr()

Destructor. ";

/*  Accessor methods  */

%feature("docstring")  Anasazi::BlockKrylovSchurSolMgr::getProblem "const Eigenproblem<ScalarType,MV,OP>& Anasazi::BlockKrylovSchurSolMgr<
ScalarType, MV, OP >::getProblem() const

Return the eigenvalue problem. ";

%feature("docstring")  Anasazi::BlockKrylovSchurSolMgr::getNumIters "int Anasazi::BlockKrylovSchurSolMgr< ScalarType, MV, OP
>::getNumIters() const

Get the iteration count for the most recent call to  solve(). ";

%feature("docstring")  Anasazi::BlockKrylovSchurSolMgr::getRitzValues
"std::vector<Value<ScalarType> > Anasazi::BlockKrylovSchurSolMgr<
ScalarType, MV, OP >::getRitzValues() const

Return the Ritz values from the most recent solve. ";

%feature("docstring")  Anasazi::BlockKrylovSchurSolMgr::getTimers "Teuchos::Array<Teuchos::RCP<Teuchos::Time> >
Anasazi::BlockKrylovSchurSolMgr< ScalarType, MV, OP >::getTimers()
const

Return the timers for this object.

The timers are ordered as follows: time spent in solve() routine

time spent restarting ";

/*  Solver application methods  */

%feature("docstring")  Anasazi::BlockKrylovSchurSolMgr::solve "ReturnType Anasazi::BlockKrylovSchurSolMgr< ScalarType, MV, OP
>::solve()

This method performs possibly repeated calls to the underlying
eigensolver's iterate() routine until the problem has been solved (as
decided by the solver manager) or the solver manager decides to quit.

This method calls BlockKrylovSchur::iterate(), which will return
either because a specially constructed status test evaluates to
::Passed or an exception is thrown.

A return from BlockKrylovSchur::iterate() signifies one of the
following scenarios: the maximum number of restarts has been exceeded.
In this scenario, the solver manager will place  all converged
eigenpairs into the eigenproblem and return ::Unconverged.

global convergence has been met. In this case, the most significant
NEV eigenpairs in the solver and locked storage  have met the
convergence criterion. (Here, NEV refers to the number of eigenpairs
requested by the Eigenproblem.)  In this scenario, the solver manager
will return ::Converged.

::ReturnType specifying: ::Converged: the eigenproblem was solved to
the specification required by the solver manager.

::Unconverged: the eigenproblem was not solved to the specification
desired by the solver manager. ";

%feature("docstring")
Anasazi::BlockKrylovSchurSolMgr::setGlobalStatusTest "void
Anasazi::BlockKrylovSchurSolMgr< ScalarType, MV, OP
>::setGlobalStatusTest(const Teuchos::RCP< StatusTest< ScalarType, MV,
OP > > &global)

Set the status test defining global convergence. ";

%feature("docstring")
Anasazi::BlockKrylovSchurSolMgr::getGlobalStatusTest "const
Teuchos::RCP< StatusTest< ScalarType, MV, OP > > &
Anasazi::BlockKrylovSchurSolMgr< ScalarType, MV, OP
>::getGlobalStatusTest() const

Get the status test defining global convergence. ";

%feature("docstring")
Anasazi::BlockKrylovSchurSolMgr::setDebugStatusTest "void
Anasazi::BlockKrylovSchurSolMgr< ScalarType, MV, OP
>::setDebugStatusTest(const Teuchos::RCP< StatusTest< ScalarType, MV,
OP > > &debug)

Set the status test for debugging. ";

%feature("docstring")
Anasazi::BlockKrylovSchurSolMgr::getDebugStatusTest "const
Teuchos::RCP< StatusTest< ScalarType, MV, OP > > &
Anasazi::BlockKrylovSchurSolMgr< ScalarType, MV, OP
>::getDebugStatusTest() const

Get the status test for debugging. ";


// File: structAnasazi_1_1BlockKrylovSchurState.xml
%feature("docstring") Anasazi::BlockKrylovSchurState "

Structure to contain pointers to BlockKrylovSchur state variables.

This struct is utilized by BlockKrylovSchur::initialize() and
BlockKrylovSchur::getState().

C++ includes: AnasaziBlockKrylovSchur.hpp ";

%feature("docstring")
Anasazi::BlockKrylovSchurState::BlockKrylovSchurState "Anasazi::BlockKrylovSchurState< ScalarType, MulVec
>::BlockKrylovSchurState() ";


// File: structAnasazi_1_1LOBPCG_1_1CheckList.xml


// File: structAnasazi_1_1BlockKrylovSchur_1_1CheckList.xml


// File: structAnasazi_1_1RTRBase_1_1CheckList.xml


// File: structAnasazi_1_1BlockDavidson_1_1CheckList.xml


// File: structAnasazi_1_1BasicSort_1_1compAlg.xml


// File: structAnasazi_1_1BasicSort_1_1compMag.xml


// File: structAnasazi_1_1BasicSort_1_1compMag2.xml


// File: classAnasazi_1_1DenseMatTraits.xml
%feature("docstring") Anasazi::DenseMatTraits "

Virtual base class which defines basic traits for the multi-vector
type.

An adapter for this traits class must exist for the DM type. If not,
this class will produce a compile-time error.

C++ includes: AnasaziDenseMatTraits.hpp ";

%feature("docstring")  Anasazi::DenseMatTraits::Clone "static
Teuchos::RCP<DM> Anasazi::DenseMatTraits< ScalarType, DM
>::Clone(const DM &dm, const int numrows, const int numcols)

Creates a new empty DM containing numvecs columns.

Reference-counted pointer to a new dense matrix of type DM. ";

%feature("docstring")  Anasazi::DenseMatTraits::CloneCopy "static
Teuchos::RCP<DM> Anasazi::DenseMatTraits< ScalarType, DM
>::CloneCopy(const DM &dm)

Creates a new DM and copies the contents of dm into the new matrix
(deep copy).

Reference-counted pointer to the new matrix of type DM. ";

%feature("docstring")  Anasazi::DenseMatTraits::CloneCopy "static
Teuchos::RCP<DM> Anasazi::DenseMatTraits< ScalarType, DM
>::CloneCopy(const DM &dm, const int numrows, const int numcols, const
int firstrow, const int firstcol)

Creates a new DM and copies the selected contents of dm into the new
matrix (deep copy).

Reference-counted pointer to the new matrix of type DM. ";

%feature("docstring")  Anasazi::DenseMatTraits::CloneView "static
Teuchos::RCP<DM> Anasazi::DenseMatTraits< ScalarType, DM
>::CloneView(DM &dm, const int numrows, const int numcols, const int
firstrow, const int firstcol)

Creates a new DM that shares the selected contents of dm (shallow
copy).

Reference-counted pointer to the new matrix of type DM. ";

%feature("docstring")  Anasazi::DenseMatTraits::CloneView "static
Teuchos::RCP<const DM> Anasazi::DenseMatTraits< ScalarType, DM
>::CloneView(const DM &dm, const int numrows, const int numcols, const
int firstrow, const int firstcol)

Creates a new DM that shares the selected contents of dm (shallow
copy).

Reference-counted pointer to the new matrix of type DM. ";

%feature("docstring")  Anasazi::DenseMatTraits::GetNumRows "static
int Anasazi::DenseMatTraits< ScalarType, DM >::GetNumRows(const DM
&dm)

Obtain the number of rows of dm. ";

%feature("docstring")  Anasazi::DenseMatTraits::GetNumCols "static
int Anasazi::DenseMatTraits< ScalarType, DM >::GetNumCols(const DM
&dm)

Obtain the number of columns of dm. ";

%feature("docstring")  Anasazi::DenseMatTraits::value "static
ScalarType& Anasazi::DenseMatTraits< ScalarType, DM >::value(DM &dm,
const int i, const int j)

Access a reference to the (i,j) entry of dm, e_i^T dm e_j. ";

%feature("docstring")  Anasazi::DenseMatTraits::value "static const
ScalarType& Anasazi::DenseMatTraits< ScalarType, DM >::value(const DM
&dm, const int i, const int j)

Access a const reference to the (i,j) entry of dm, e_i^T dm e_j. ";

%feature("docstring")  Anasazi::DenseMatTraits::values "static
ScalarType* Anasazi::DenseMatTraits< ScalarType, DM >::values(DM &dm,
int *stride, bool *cor)

Access the pointers to the data in dm, information enabling C-style
access to the data.

The return value is a pointer to the data, stored sequentially. *cor
denotes whether the data is stored column-oriented ( $*cor == true$)
or row-oriented ($ == false$). *stride denotes the stride between
columns/rows. ";

%feature("docstring")  Anasazi::DenseMatTraits::MvTimesMatAddMv "static void Anasazi::DenseMatTraits< ScalarType, DM
>::MvTimesMatAddMv(bool transA, bool transB, const ScalarType alpha,
const DM &A, const DM &B, const ScalarType beta, DM &dm)

Update dm with $ \\\\alpha op(A) op(B) + \\\\beta dm $, where op(A) :=
A^H if transA == true op(A) := A if transA == false.

and

op(B) := B^H if transB == true op(B) := B if transB == false ";

%feature("docstring")  Anasazi::DenseMatTraits::MvAddMv "static void
Anasazi::DenseMatTraits< ScalarType, DM >::MvAddMv(const ScalarType
alpha, const DM &A, const ScalarType beta, const DM &B, DM &dm)

Replace dm with $\\\\alpha A + \\\\beta B$. ";

%feature("docstring")  Anasazi::DenseMatTraits::DMRandom "static void
Anasazi::DenseMatTraits< ScalarType, DM >::DMRandom(DM &dm)

Replace the entries of dm with random numbers. ";

%feature("docstring")  Anasazi::DenseMatTraits::DMInit "static void
Anasazi::DenseMatTraits< ScalarType, DM >::DMInit(DM &dm, const
ScalarType alpha=Teuchos::ScalarTraits< ScalarType >::zero())

Replace each element of mv with alpha. ";

%feature("docstring")  Anasazi::DenseMatTraits::DMPrint "static void
Anasazi::DenseMatTraits< ScalarType, DM >::DMPrint(const DM &dm,
std::ostream &os)

Print the matrix dm to the output stream os. ";


// File: classAnasazi_1_1DirectSolver.xml
%feature("docstring") Anasazi::DirectSolver "

Anasazi's templated abstract base class providing solver capabilities
for projected eigenproblems.

This class provides concrete, templated implementations of solvers for
projected eigenproblems.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziDirectSolver.hpp ";

/*  Constructor/Destructor  */

%feature("docstring")  Anasazi::DirectSolver::DirectSolver "Anasazi::DirectSolver< ScalarType >::DirectSolver()

Basic constructor. ";

%feature("docstring")  Anasazi::DirectSolver::~DirectSolver "virtual
Anasazi::DirectSolver< ScalarType >::~DirectSolver()

Destructor. ";

/*  Eigensolver Projection Methods  */

%feature("docstring")  Anasazi::DirectSolver::directSolver "virtual
int Anasazi::DirectSolver< ScalarType >::directSolver(int size, const
Teuchos::SerialDenseMatrix< int, ScalarType > &KK, const
Teuchos::SerialDenseMatrix< int, ScalarType > *MM,
Teuchos::SerialDenseMatrix< int, ScalarType > &EV, std::vector<
typename Teuchos::ScalarTraits< ScalarType >::magnitudeType > &theta)
const =0

Routine for computing the generalized eigenpairs of the symmetric
pencil (KK, MM)

Parameters:
-----------

size:  [in] Dimension of the eigenproblem (KK, MM)

KK:  [in] Symmetric \"stiffness\" matrix

MM:  [in] Symmetric Positive \"mass\" matrix

EV:  [in] Dense matrix to store the nev eigenvectors

theta:  [out] Array to store the eigenvalues

The code accesses only the upper triangular part of KK and MM.

Integer indicating the number of computed eigenpairs. ";


// File: classAnasazi_1_1Eigenproblem.xml
%feature("docstring") Anasazi::Eigenproblem "

This class defines the interface required by an eigensolver and status
test class to compute solutions to an eigenproblem.

C++ includes: AnasaziEigenproblem.hpp ";

/*  Constructors/Destructor  */

%feature("docstring")  Anasazi::Eigenproblem::Eigenproblem "Anasazi::Eigenproblem< ScalarType, MV, OP >::Eigenproblem()

Empty constructor. ";

%feature("docstring")  Anasazi::Eigenproblem::~Eigenproblem "virtual
Anasazi::Eigenproblem< ScalarType, MV, OP >::~Eigenproblem()

Destructor. ";

/*  Set Methods  */

%feature("docstring")  Anasazi::Eigenproblem::setOperator "virtual
void Anasazi::Eigenproblem< ScalarType, MV, OP >::setOperator(const
Teuchos::RCP< const OP > &Op)=0

Set the operator for which eigenvalues will be computed.

This may be different from the A if a spectral transformation is
employed. For example, this operator may apply the operation
$(A-\\\\sigma I)^{-1}$ if you are looking for eigenvalues of A around
$\\\\sigma$. ";

%feature("docstring")  Anasazi::Eigenproblem::setA "virtual void
Anasazi::Eigenproblem< ScalarType, MV, OP >::setA(const Teuchos::RCP<
const OP > &A)=0

Set the operator A of the eigenvalue problem $Ax=\\\\lambda Mx$. ";

%feature("docstring")  Anasazi::Eigenproblem::setM "virtual void
Anasazi::Eigenproblem< ScalarType, MV, OP >::setM(const Teuchos::RCP<
const OP > &M)=0

Set the operator M of the eigenvalue problem $Ax=\\\\lambda Mx$. ";

%feature("docstring")  Anasazi::Eigenproblem::setPrec "virtual void
Anasazi::Eigenproblem< ScalarType, MV, OP >::setPrec(const
Teuchos::RCP< const OP > &Prec)=0

Set the preconditioner for this eigenvalue problem $Ax=\\\\lambda Mx$.
";

%feature("docstring")  Anasazi::Eigenproblem::setInitVec "virtual
void Anasazi::Eigenproblem< ScalarType, MV, OP >::setInitVec(const
Teuchos::RCP< MV > &InitVec)=0

Set the initial guess.

This multivector should have the same number of columns as the
blocksize. ";

%feature("docstring")  Anasazi::Eigenproblem::setAuxVecs "virtual
void Anasazi::Eigenproblem< ScalarType, MV, OP >::setAuxVecs(const
Teuchos::RCP< const MV > &AuxVecs)=0

Set auxiliary vectors.

This multivector can have any number of columns, and most likely will
contain vectors that will be used by the eigensolver to orthogonalize
against. ";

%feature("docstring")  Anasazi::Eigenproblem::setNEV "virtual void
Anasazi::Eigenproblem< ScalarType, MV, OP >::setNEV(int nev)=0

The number of eigenvalues (NEV) that are requested. ";

%feature("docstring")  Anasazi::Eigenproblem::setHermitian "virtual
void Anasazi::Eigenproblem< ScalarType, MV, OP >::setHermitian(bool
isSym)=0

Specify the symmetry of the eigenproblem.

This knowledge may allow the solver to take advantage of the
eigenproblems' symmetry. Some computational work may be avoided by
setting this properly. ";

%feature("docstring")  Anasazi::Eigenproblem::setProblem "virtual
bool Anasazi::Eigenproblem< ScalarType, MV, OP >::setProblem()=0

Specify that this eigenproblem is fully defined.

This routine serves multiple purpose: sanity check that the
eigenproblem has been fully and consistently defined

opportunity for the eigenproblem to allocate internal storage for
eigenvalues and eigenvectors (to be used by eigensolvers and solver
managers)

The user MUST call this routine before they send the eigenproblem to
any solver or solver manager.

true signifies success, false signifies error. ";

%feature("docstring")  Anasazi::Eigenproblem::setSolution "virtual
void Anasazi::Eigenproblem< ScalarType, MV, OP >::setSolution(const
Eigensolution< ScalarType, MV > &sol)=0

Set the solution to the eigenproblem.

This mechanism allows an Eigensolution struct to be associated with an
Eigenproblem object. setSolution() is usually called by a solver
manager at the end of its SolverManager::solve() routine. ";

/*  Accessor Methods  */

%feature("docstring")  Anasazi::Eigenproblem::getOperator "virtual
Teuchos::RCP<const OP> Anasazi::Eigenproblem< ScalarType, MV, OP
>::getOperator() const =0

Get a pointer to the operator for which eigenvalues will be computed.
";

%feature("docstring")  Anasazi::Eigenproblem::getA "virtual
Teuchos::RCP<const OP> Anasazi::Eigenproblem< ScalarType, MV, OP
>::getA() const =0

Get a pointer to the operator A of the eigenproblem $AX=\\\\lambda
Mx$. ";

%feature("docstring")  Anasazi::Eigenproblem::getM "virtual
Teuchos::RCP<const OP> Anasazi::Eigenproblem< ScalarType, MV, OP
>::getM() const =0

Get a pointer to the operator M of the eigenproblem $AX=\\\\lambda
Mx$. ";

%feature("docstring")  Anasazi::Eigenproblem::getPrec "virtual
Teuchos::RCP<const OP> Anasazi::Eigenproblem< ScalarType, MV, OP
>::getPrec() const =0

Get a pointer to the preconditioner. ";

%feature("docstring")  Anasazi::Eigenproblem::getInitVec "virtual
Teuchos::RCP<const MV> Anasazi::Eigenproblem< ScalarType, MV, OP
>::getInitVec() const =0

Get a pointer to the initial vector. ";

%feature("docstring")  Anasazi::Eigenproblem::getAuxVecs "virtual
Teuchos::RCP<const MV> Anasazi::Eigenproblem< ScalarType, MV, OP
>::getAuxVecs() const =0

Get a pointer to the auxiliary vector. ";

%feature("docstring")  Anasazi::Eigenproblem::getNEV "virtual int
Anasazi::Eigenproblem< ScalarType, MV, OP >::getNEV() const =0

Get the number of eigenvalues (NEV) that are required by this
eigenproblem. ";

%feature("docstring")  Anasazi::Eigenproblem::isHermitian "virtual
bool Anasazi::Eigenproblem< ScalarType, MV, OP >::isHermitian() const
=0

Get the symmetry information for this eigenproblem. ";

%feature("docstring")  Anasazi::Eigenproblem::isProblemSet "virtual
bool Anasazi::Eigenproblem< ScalarType, MV, OP >::isProblemSet() const
=0

If the problem has been set, this method will return true. ";

%feature("docstring")  Anasazi::Eigenproblem::getSolution "virtual
const Eigensolution<ScalarType,MV>& Anasazi::Eigenproblem< ScalarType,
MV, OP >::getSolution() const =0

Get the solution to the eigenproblem.

There is no computation associated with this method. It only provides
a mechanism for associating an Eigensolution with a Eigenproblem. ";


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

/*  Constructors/Destructor  */

%feature("docstring")  Anasazi::Eigensolver::Eigensolver "Anasazi::Eigensolver< ScalarType, MV, OP >::Eigensolver()

Default Constructor. ";

%feature("docstring")  Anasazi::Eigensolver::Eigensolver "Anasazi::Eigensolver< ScalarType, MV, OP >::Eigensolver(const
Teuchos::RCP< Eigenproblem< ScalarType, MV, OP > > &problem, const
Teuchos::RCP< SortManager< ScalarType > > &sorter, const Teuchos::RCP<
OutputManager< ScalarType > > &printer, const Teuchos::RCP<
StatusTest< ScalarType, MV, OP > > &tester, const Teuchos::RCP<
OrthoManager< ScalarType, MV > > &ortho, Teuchos::ParameterList
&params)

Basic Constructor.

This constructor, implemented by all Anasazi eigensolvers, takes an
Anasazi::Eigenproblem, Anasazi::SortManager, Anasazi::OutputManager,
and Teuchos::ParameterList as input. These four arguments are
sufficient enough for constructing any Anasazi::Eigensolver object. ";

%feature("docstring")  Anasazi::Eigensolver::~Eigensolver "virtual
Anasazi::Eigensolver< ScalarType, MV, OP >::~Eigensolver()

Destructor. ";

/*  Solver methods  */

%feature("docstring")  Anasazi::Eigensolver::iterate "virtual void
Anasazi::Eigensolver< ScalarType, MV, OP >::iterate()=0

This method performs eigensolvers iterations until the status test
indicates the need to stop or an error occurs (in which case, an
exception is thrown). ";

%feature("docstring")  Anasazi::Eigensolver::initialize "virtual void
Anasazi::Eigensolver< ScalarType, MV, OP >::initialize()=0

Initialize the solver with the initial vectors from the eigenproblem
or random data. ";

/*  Status methods  */

%feature("docstring")  Anasazi::Eigensolver::getNumIters "virtual int
Anasazi::Eigensolver< ScalarType, MV, OP >::getNumIters() const =0

Get the current iteration count. ";

%feature("docstring")  Anasazi::Eigensolver::resetNumIters "virtual
void Anasazi::Eigensolver< ScalarType, MV, OP >::resetNumIters()=0

Reset the iteration count. ";

%feature("docstring")  Anasazi::Eigensolver::getRitzVectors "virtual
Teuchos::RCP<const MV> Anasazi::Eigensolver< ScalarType, MV, OP
>::getRitzVectors()=0

Get the Ritz vectors from the previous iteration. These are indexed
using getRitzIndex().

For a description of the indexing scheme, see getRitzIndex(). ";

%feature("docstring")  Anasazi::Eigensolver::getRitzValues "virtual
std::vector<Value<ScalarType> > Anasazi::Eigensolver< ScalarType, MV,
OP >::getRitzValues()=0

Get the Ritz values from the previous iteration. ";

%feature("docstring")  Anasazi::Eigensolver::getRitzIndex "virtual
std::vector<int> Anasazi::Eigensolver< ScalarType, MV, OP
>::getRitzIndex()=0

Get the index used for indexing the compressed storage used for Ritz
vectors for real, non-Hermitian problems.

index has length numVecs, where each entry is 0, +1, or -1. These have
the following interpretation: index[i] == 0: signifies that the
corresponding eigenvector is stored as the i column of Evecs. This
will usually be the case when ScalarType is complex, an eigenproblem
is Hermitian, or a real, non- Hermitian eigenproblem has a real
eigenvector.

index[i] == +1: signifies that the corresponding eigenvector is stored
in two vectors: the real part in the i column of Evecs and the
positive imaginary part in the i+1 column of Evecs.

index[i] == -1: signifies that the corresponding eigenvector is stored
in two vectors: the real part in the i-1 column of Evecs and the
negative imaginary part in the i column of Evecs ";

%feature("docstring")  Anasazi::Eigensolver::getResNorms "virtual
std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>
Anasazi::Eigensolver< ScalarType, MV, OP >::getResNorms()=0

Get the current residual norms.

A vector of length blockSize containing the norms of the residuals,
according to the orthogonalization manager norm() method. ";

%feature("docstring")  Anasazi::Eigensolver::getRes2Norms "virtual
std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>
Anasazi::Eigensolver< ScalarType, MV, OP >::getRes2Norms()=0

Get the current residual 2-norms A vector of length blockSize
containing the 2-norms of the residuals. ";

%feature("docstring")  Anasazi::Eigensolver::getRitzRes2Norms "virtual std::vector<typename
Teuchos::ScalarTraits<ScalarType>::magnitudeType>
Anasazi::Eigensolver< ScalarType, MV, OP >::getRitzRes2Norms()=0

Get the 2-norms of the Ritz residuals. A vector of length blockSize
containing the 2-norms of the Ritz residuals. ";

%feature("docstring")  Anasazi::Eigensolver::getCurSubspaceDim "virtual int Anasazi::Eigensolver< ScalarType, MV, OP
>::getCurSubspaceDim() const =0

Get the dimension of the search subspace used to generate the current
eigenvectors and eigenvalues. ";

%feature("docstring")  Anasazi::Eigensolver::getMaxSubspaceDim "virtual int Anasazi::Eigensolver< ScalarType, MV, OP
>::getMaxSubspaceDim() const =0

Get the maximum dimension allocated for the search subspace. ";

/*  Accessor methods  */

%feature("docstring")  Anasazi::Eigensolver::setStatusTest "virtual
void Anasazi::Eigensolver< ScalarType, MV, OP
>::setStatusTest(Teuchos::RCP< StatusTest< ScalarType, MV, OP > >
test)=0

Set a new StatusTest for the solver. ";

%feature("docstring")  Anasazi::Eigensolver::getStatusTest "virtual
Teuchos::RCP<StatusTest<ScalarType,MV,OP> > Anasazi::Eigensolver<
ScalarType, MV, OP >::getStatusTest() const =0

Get the current StatusTest used by the solver. ";

%feature("docstring")  Anasazi::Eigensolver::getProblem "virtual
const Eigenproblem<ScalarType,MV,OP>& Anasazi::Eigensolver<
ScalarType, MV, OP >::getProblem() const =0

Get a constant reference to the eigenvalue problem. ";

%feature("docstring")  Anasazi::Eigensolver::getBlockSize "virtual
int Anasazi::Eigensolver< ScalarType, MV, OP >::getBlockSize() const
=0

Get the blocksize to be used by the iterative solver in solving this
eigenproblem. ";

%feature("docstring")  Anasazi::Eigensolver::setBlockSize "virtual
void Anasazi::Eigensolver< ScalarType, MV, OP >::setBlockSize(int
blockSize)=0

Set the blocksize to be used by the iterative solver in solving this
eigenproblem. ";

%feature("docstring")  Anasazi::Eigensolver::setAuxVecs "virtual void
Anasazi::Eigensolver< ScalarType, MV, OP >::setAuxVecs(const
Teuchos::Array< Teuchos::RCP< const MV > > &auxvecs)=0

Set the auxiliary vectors for the solver. ";

%feature("docstring")  Anasazi::Eigensolver::getAuxVecs "virtual
Teuchos::Array<Teuchos::RCP<const MV> > Anasazi::Eigensolver<
ScalarType, MV, OP >::getAuxVecs() const =0

Get the auxiliary vectors for the solver. ";

%feature("docstring")  Anasazi::Eigensolver::isInitialized "virtual
bool Anasazi::Eigensolver< ScalarType, MV, OP >::isInitialized() const
=0

States whether the solver has been initialized or not. ";

/*  Output methods  */

%feature("docstring")  Anasazi::Eigensolver::currentStatus "virtual
void Anasazi::Eigensolver< ScalarType, MV, OP
>::currentStatus(std::ostream &os)=0

This method requests that the solver print out its current status to
screen. ";


// File: classAnasazi_1_1GenOrthoManager.xml
%feature("docstring") Anasazi::GenOrthoManager "

This class provides an interface for orthogonalization managers to
provide oblique projectors of the form: \\\\[ P_{X,Y} S = S - X
\\\\langle Y, X \\\\rangle^{-1} \\\\langle Y, S \\\\rangle\\\\ . \\\\]
Such a projector modifies the input in the range on $X$ in order to
make the output orthogonal to the range of $Y$.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziGenOrthoManager.hpp ";

/*  Constructor/Destructor  */

%feature("docstring")  Anasazi::GenOrthoManager::GenOrthoManager "Anasazi::GenOrthoManager< ScalarType, MV, OP
>::GenOrthoManager(Teuchos::RCP< const OP > Op=Teuchos::null)

Default constructor. ";

%feature("docstring")  Anasazi::GenOrthoManager::~GenOrthoManager "virtual Anasazi::GenOrthoManager< ScalarType, MV, OP
>::~GenOrthoManager()

Destructor. ";

/*  Orthogonalization methods  */

%feature("docstring")  Anasazi::GenOrthoManager::projectGen "virtual
void Anasazi::GenOrthoManager< ScalarType, MV, OP >::projectGen(MV &S,
Teuchos::Array< Teuchos::RCP< const MV > > X, Teuchos::Array<
Teuchos::RCP< const MV > > Y, bool isBiOrtho, Teuchos::Array<
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > > >
C=Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix< int,
ScalarType > >(Teuchos::null)), Teuchos::RCP< MV > MS=Teuchos::null,
Teuchos::Array< Teuchos::RCP< const MV > >
MX=Teuchos::tuple(Teuchos::RCP< const MV >(Teuchos::null)),
Teuchos::Array< Teuchos::RCP< const MV > >
MY=Teuchos::tuple(Teuchos::RCP< const MV >(Teuchos::null))) const =0

Applies a series of generic projectors.

Given a list of bases X[i] and Y[i] (a projection pair), this method
takes a multivector S and applies the projectors \\\\[ P_{X[i],Y[i]} S
= S - X[i] \\\\langle Y[i], X[i] \\\\rangle^{-1} \\\\langle Y[i], S
\\\\rangle\\\\ . \\\\] This operation projects S onto the space
orthogonal to the Y[i], along the range of the X[i]. The inner product
specified by $\\\\langle \\\\cdot, \\\\cdot \\\\rangle$ is given by
innerProd().

The call is equivalent to the call  The method also returns the
coefficients C[i] associated with each projection pair, so that \\\\[
S_{in} = S_{out} + \\\\sum_i X[i] C[i] \\\\] and therefore \\\\[ C[i]
= \\\\langle Y[i], X[i] \\\\rangle^{-1} \\\\langle Y[i], S
\\\\rangle\\\\ . \\\\]

Lastly, for reasons of efficiency, the user must specify whether the
projection pairs are bi-orthonormal with respect to innerProd(), i.e.,
whether $\\\\langle Y[i], X[i] \\\\rangle = I$. In the case that the
bases are specified to be biorthogonal, the inverse $\\\\langle Y, X
\\\\rangle^{-1}$ will not be computed. Furthermore, the user may
optionally specifiy the image of S and the projection pairs under the
inner product operator getOp().

Parameters:
-----------

S:  [in/out] The multivector to be modified.  On output, the columns
of S will be orthogonal to each Y[i], satisfying \\\\[ \\\\langle
Y[i], S_{out} \\\\rangle = 0 \\\\] Also, \\\\[ S_{in} = S_{out} +
\\\\sum_i X[i] C[i] \\\\]

X:  [in] Multivectors for bases under which $S_{in}$ is modified.

Y:  [in] Multivectors for bases to which $S_{out}$ should be
orthogonal.

isBiortho:  [in] A flag specifying whether the bases X[i] and Y[i] are
biorthonormal, i.e,. whether $\\\\langle Y[i], X[i]\\\\rangle == I$.

C:  [out] Coefficients for reconstructing $S_{in}$ via the bases X[i].
If C[i] is a non-null pointer and C[i] matches the dimensions of S and
X[i], then the coefficients computed during the orthogonalization
routine will be stored in the matrix C[i].  If C[i] points to a
Teuchos::SerialDenseMatrix with size inconsistent with S and  X[i],
then a std::invalid_argument exception will be thrown.  Otherwise, if
C.size() < i or C[i] is a null pointer, the caller will not have
access to the computed coefficients C[i].

MS:  [in/out] If specified by the user, on input MS is required to be
the image of S under the operator getOp(). On output, MS will be
updated to reflect the changes in S.

MX:  [in] If specified by the user, MX[i] is required to be the image
of X[i] under the operator getOp().

MY:  [in] If specified by the user, MY[i] is required to be the image
of Y[i] under the operator getOp().

If X[i] != Teuchos::null or Y[i] != Teuchos::null, then X[i] and Y[i]
are required to have the same number of columns, and each should have
the same number of rows as S.

For any i != j, $\\\\langle Y[i], X[j] \\\\rangle == 0$.

If biOrtho == true, $\\\\langle Y[i], X[i]\\\\rangle == I$

Otherwise, if biOrtho == false, then $\\\\langle Y[i], X[i]\\\\rangle$
should be Hermitian positive-definite.

If X[i] and Y[i] have $xc_i$ columns and S has $sc$ columns, then C[i]
if specified must be $xc_i \\\\times sc$. ";

%feature("docstring")
Anasazi::GenOrthoManager::projectAndNormalizeGen "virtual int
Anasazi::GenOrthoManager< ScalarType, MV, OP
>::projectAndNormalizeGen(MV &S, Teuchos::Array< Teuchos::RCP< const
MV > > X, Teuchos::Array< Teuchos::RCP< const MV > > Y, bool
isBiOrtho, Teuchos::Array< Teuchos::RCP< Teuchos::SerialDenseMatrix<
int, ScalarType > > > C=Teuchos::tuple(Teuchos::RCP<
Teuchos::SerialDenseMatrix< int, ScalarType > >(Teuchos::null)),
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > >
B=Teuchos::null, Teuchos::RCP< MV > MS=Teuchos::null, Teuchos::Array<
Teuchos::RCP< const MV > > MX=Teuchos::tuple(Teuchos::RCP< const MV
>(Teuchos::null)), Teuchos::Array< Teuchos::RCP< const MV > >
MY=Teuchos::tuple(Teuchos::RCP< const MV >(Teuchos::null))) const =0

Applies a series of generic projectors and returns an orthonormal
basis for the residual data.

Given a list of bases X[i] and Y[i] (a projection pair), this method
takes a multivector S and applies the projectors \\\\[ P_{X[i],Y[i]} S
= S - X[i] \\\\langle Y[i], X[i] \\\\rangle^{-1} \\\\langle Y[i], S
\\\\rangle\\\\ . \\\\] These operation project S onto the space
orthogonal to the range of the Y[i], along the range of X[i]. The
inner product specified by $\\\\langle \\\\cdot, \\\\cdot \\\\rangle$
is given by innerProd().

The method returns in S an orthonormal basis for the residual \\\\[
\\\\left( \\\\prod_{i} P_{X[i],Y[i]} \\\\right) S_{in} = S_{out} B\\\\
, \\\\] where B contains the (not necessarily triangular) coefficients
of the residual with respect to the new basis.

The method also returns the coefficients C[i] and B associated with
each projection pair, so that \\\\[ S_{in} = S_{out} B + \\\\sum_i
X[i] C[i] \\\\] and \\\\[ C[i] = \\\\langle Y[i], X[i] \\\\rangle^{-1}
\\\\langle Y[i], S \\\\rangle\\\\ . \\\\]

Lastly, for reasons of efficiency, the user must specify whether the
projection pairs are bi-orthonormal with respect to innerProd(), i.e.,
whether $\\\\langle Y[i], X[i] \\\\rangle = I$. Furthermore, the user
may optionally specifiy the image of S and the projection pairs under
the inner product operator getOp().

Parameters:
-----------

S:  [in/out] The multivector to be modified.  On output, the columns
of S will be orthogonal to each Y[i], satisfying \\\\[ \\\\langle
Y[i], S_{out} \\\\rangle = 0 \\\\] Also, \\\\[ S_{in}(1:m,1:n) =
S_{out}(1:m,1:rank) B(1:rank,1:n) + \\\\sum_i X[i] C[i]\\\\ , \\\\]
where m is the number of rows in S, n is the number of columns in S,
and rank is the value returned from the method.

X:  [in] Multivectors for bases under which $S_{in}$ is modified.

Y:  [in] Multivectors for bases to which $S_{out}$ should be
orthogonal.

isBiortho:  [in] A flag specifying whether the bases X[i] and Y[i] are
biorthonormal, i.e,. whether $\\\\langle Y[i], X[i]\\\\rangle == I$.

C:  [out] Coefficients for reconstructing $S_{in}$ via the bases X[i].
If C[i] is a non-null pointer and C[i] matches the dimensions of X and
Q[i], then the coefficients computed during the orthogonalization
routine will be stored in the matrix C[i].  If C[i] points to a
Teuchos::SerialDenseMatrix with size inconsistent with S and  X[i],
then a std::invalid_argument exception will be thrown.  Otherwise, if
C.size() < i or C[i] is a null pointer, the caller will not have
access to the computed coefficients C[i].

B:  [out] The coefficients of the original S with respect to the
computed basis. If B is a non-null pointer and B matches the
dimensions of B, then the coefficients computed during the
orthogonalization routine will be stored in B, similar to calling If B
points to a Teuchos::SerialDenseMatrix with size inconsistent with S,
then a std::invalid_argument exception will be thrown.  Otherwise, if
B is null, the caller will not have access to the computed
coefficients.

MS:  [in/out] If specified by the user, on input MS is required to be
the image of S under the operator getOp(). On output, MS will be
updated to reflect the changes in S.

MX:  [in] If specified by the user, MX[i] is required to be the image
of X[i] under the operator getOp().

MY:  [in] If specified by the user, MY[i] is required to be the image
of Y[i] under the operator getOp().

The matrix B is not necessarily triangular (as in a QR factorization);
see the documentation of specific orthogonalization managers.

If X[i] != Teuchos::null or Y[i] != Teuchos::null, then X[i] and Y[i]
are required to have the same number of columns, and each should have
the same number of rows as S.

For any i != j, $\\\\langle Y[i], X[j] \\\\rangle == 0$.

If biOrtho == true, $\\\\langle Y[i], X[i]\\\\rangle == I$

Otherwise, if biOrtho == false, then $\\\\langle Y[i], X[i]\\\\rangle$
should be Hermitian positive-definite.

If X[i] and Y[i] have $xc_i$ columns and S has $sc$ columns, then C[i]
if specified must be $xc_i \\\\times sc$.

If S has $sc$ columns, then B if specified must be $sc \\\\times sc $.

Rank of the basis computed by this method. ";


// File: classAnasazi_1_1HelperTraits.xml
%feature("docstring") Anasazi::HelperTraits "

Class which defines basic traits for working with different scalar
types.

An adapter for this traits class must exist for the ScalarType. If
not, this class will produce a compile-time error.

C++ includes: AnasaziHelperTraits.hpp ";


// File: classAnasazi_1_1ICGSOrthoManager.xml
%feature("docstring") Anasazi::ICGSOrthoManager "

An implementation of the Anasazi::GenOrthoManager that performs
orthogonalization using iterated classical Gram- Schmidt.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziICGSOrthoManager.hpp ";

/*  Constructor/Destructor  */

%feature("docstring")  Anasazi::ICGSOrthoManager::ICGSOrthoManager "Anasazi::ICGSOrthoManager< ScalarType, MV, OP
>::ICGSOrthoManager(Teuchos::RCP< const OP > Op=Teuchos::null, int
numIters=2, typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType eps=0.0, typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType tol=0.20)

Constructor specifying the operator defining the inner product as well
as the number of orthogonalization iterations. ";

%feature("docstring")  Anasazi::ICGSOrthoManager::~ICGSOrthoManager "Anasazi::ICGSOrthoManager< ScalarType, MV, OP >::~ICGSOrthoManager()

Destructor. ";

/*  Methods implementing Anasazi::GenOrthoManager  */

%feature("docstring")  Anasazi::ICGSOrthoManager::projectGen "void
Anasazi::ICGSOrthoManager< ScalarType, MV, OP >::projectGen(MV &S,
Teuchos::Array< Teuchos::RCP< const MV > > X, Teuchos::Array<
Teuchos::RCP< const MV > > Y, bool isBiOrtho, Teuchos::Array<
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > > >
C=Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix< int,
ScalarType > >(Teuchos::null)), Teuchos::RCP< MV > MS=Teuchos::null,
Teuchos::Array< Teuchos::RCP< const MV > >
MX=Teuchos::tuple(Teuchos::RCP< const MV >(Teuchos::null)),
Teuchos::Array< Teuchos::RCP< const MV > >
MY=Teuchos::tuple(Teuchos::RCP< const MV >(Teuchos::null))) const

Applies a series of generic projectors.

Given a list of bases X[i] and Y[i] (a projection pair), this method
takes a multivector S and applies the projectors \\\\[ P_{X[i],Y[i]} S
= S - X[i] \\\\langle Y[i], X[i] \\\\rangle^{-1} \\\\langle Y[i], S
\\\\rangle\\\\ . \\\\] This operation projects S onto the space
orthogonal to the Y[i], along the range of the X[i]. The inner product
specified by $\\\\langle \\\\cdot, \\\\cdot \\\\rangle$ is given by
innerProd().

The call is equivalent to the call  The method also returns the
coefficients C[i] associated with each projection pair, so that \\\\[
S_{in} = S_{out} + \\\\sum_i X[i] C[i] \\\\] and therefore \\\\[ C[i]
= \\\\langle Y[i], X[i] \\\\rangle^{-1} \\\\langle Y[i], S
\\\\rangle\\\\ . \\\\]

Lastly, for reasons of efficiency, the user must specify whether the
projection pairs are bi-orthonormal with respect to innerProd(), i.e.,
whether $\\\\langle Y[i], X[i] \\\\rangle = I$. In the case that the
bases are specified to be biorthogonal, the inverse $\\\\langle Y, X
\\\\rangle^{-1}$ will not be computed. Furthermore, the user may
optionally specifiy the image of S and the projection pairs under the
inner product operator getOp().

projectGen() is implemented to apply the projectors via an iterated
Classical Gram-Schmidt, where the iteration is performed getNumIters()
number of times.

Parameters:
-----------

S:  [in/out] The multivector to be modified.  On output, the columns
of S will be orthogonal to each Y[i], satisfying \\\\[ \\\\langle
Y[i], S_{out} \\\\rangle = 0 \\\\] Also, \\\\[ S_{in} = S_{out} +
\\\\sum_i X[i] C[i] \\\\]

X:  [in] Multivectors for bases under which $S_{in}$ is modified.

Y:  [in] Multivectors for bases to which $S_{out}$ should be
orthogonal.

isBiortho:  [in] A flag specifying whether the bases X[i] and Y[i] are
biorthonormal, i.e,. whether $\\\\langle Y[i], X[i]\\\\rangle == I$.

C:  [out] Coefficients for reconstructing $S_{in}$ via the bases X[i].
If C[i] is a non-null pointer and C[i] matches the dimensions of S and
X[i], then the coefficients computed during the orthogonalization
routine will be stored in the matrix C[i].  If C[i] points to a
Teuchos::SerialDenseMatrix with size inconsistent with S and  X[i],
then a std::invalid_argument exception will be thrown.  Otherwise, if
C.size() < i or C[i] is a null pointer, the caller will not have
access to the computed coefficients C[i].

MS:  [in/out] If specified by the user, on input MS is required to be
the image of S under the operator getOp(). On output, MS will be
updated to reflect the changes in S.

MX:  [in] If specified by the user, on MX[i] is required to be the
image of X[i] under the operator getOp().

MY:  [in] If specified by the user, on MY[i] is required to be the
image of Y[i] under the operator getOp().

If X[i] != Teuchos::null or Y[i] != Teuchos::null, then X[i] and Y[i]
are required to have the same number of columns, and each should have
the same number of rows as S.

For any i != j, $\\\\langle Y[i], X[j] \\\\rangle == 0$.

If biOrtho == true, $\\\\langle Y[i], X[i]\\\\rangle == I$

Otherwise, if biOrtho == false, then $\\\\langle Y[i], X[i]\\\\rangle$
should be Hermitian positive-definite.

If X[i] and Y[i] have $xc_i$ columns and S has $sc$ columns, then C[i]
if specified must be $xc_i \\\\times sc$. ";

%feature("docstring")
Anasazi::ICGSOrthoManager::projectAndNormalizeGen "int
Anasazi::ICGSOrthoManager< ScalarType, MV, OP
>::projectAndNormalizeGen(MV &S, Teuchos::Array< Teuchos::RCP< const
MV > > X, Teuchos::Array< Teuchos::RCP< const MV > > Y, bool
isBiOrtho, Teuchos::Array< Teuchos::RCP< Teuchos::SerialDenseMatrix<
int, ScalarType > > > C=Teuchos::tuple(Teuchos::RCP<
Teuchos::SerialDenseMatrix< int, ScalarType > >(Teuchos::null)),
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > >
B=Teuchos::null, Teuchos::RCP< MV > MS=Teuchos::null, Teuchos::Array<
Teuchos::RCP< const MV > > MX=Teuchos::tuple(Teuchos::RCP< const MV
>(Teuchos::null)), Teuchos::Array< Teuchos::RCP< const MV > >
MY=Teuchos::tuple(Teuchos::RCP< const MV >(Teuchos::null))) const

Applies a series of generic projectors and returns an orthonormal
basis for the residual data.

Given a list of bases X[i] and Y[i] (a projection pair), this method
takes a multivector S and applies the projectors \\\\[ P_{X[i],Y[i]} S
= S - X[i] \\\\langle Y[i], X[i] \\\\rangle^{-1} \\\\langle Y[i], S
\\\\rangle\\\\ . \\\\] These operation project S onto the space
orthogonal to the range of the Y[i], along the range of X[i]. The
inner product specified by $\\\\langle \\\\cdot, \\\\cdot \\\\rangle$
is given by innerProd().

The method returns in S an orthonormal basis for the residual \\\\[
\\\\left( \\\\prod_{i} P_{X[i],Y[i]} \\\\right) S_{in} = S_{out} B\\\\
, \\\\] where B contains the (not necessarily triangular) coefficients
of the residual with respect to the new basis.

The method also returns the coefficients C[i] and B associated with
each projection pair, so that \\\\[ S_{in} = S_{out} B + \\\\sum_i
X[i] C[i] \\\\] and \\\\[ C[i] = \\\\langle Y[i], X[i] \\\\rangle^{-1}
\\\\langle Y[i], S \\\\rangle\\\\ . \\\\]

Lastly, for reasons of efficiency, the user must specify whether the
projection pairs are bi-orthonormal with respect to innerProd(), i.e.,
whether $\\\\langle Y[i], X[i] \\\\rangle = I$. Furthermore, the user
may optionally specifiy the image of S and the projection pairs under
the inner product operator getOp().

Parameters:
-----------

S:  [in/out] The multivector to be modified.  On output, the columns
of S will be orthogonal to each Y[i], satisfying \\\\[ \\\\langle
Y[i], S_{out} \\\\rangle = 0 \\\\] Also, \\\\[ S_{in}(1:m,1:n) =
S_{out}(1:m,1:rank) B(1:rank,1:n) + \\\\sum_i X[i] C[i]\\\\ , \\\\]
where m is the number of rows in S, n is the number of columns in S,
and rank is the value returned from the method.

X:  [in] Multivectors for bases under which $S_{in}$ is modified.

Y:  [in] Multivectors for bases to which $S_{out}$ should be
orthogonal.

isBiortho:  [in] A flag specifying whether the bases X[i] and Y[i] are
biorthonormal, i.e,. whether $\\\\langle Y[i], X[i]\\\\rangle == I$.

C:  [out] Coefficients for reconstructing $S_{in}$ via the bases X[i].
If C[i] is a non-null pointer and C[i] matches the dimensions of X and
Q[i], then the coefficients computed during the orthogonalization
routine will be stored in the matrix C[i].  If C[i] points to a
Teuchos::SerialDenseMatrix with size inconsistent with S and  X[i],
then a std::invalid_argument exception will be thrown.  Otherwise, if
C.size() < i or C[i] is a null pointer, the caller will not have
access to the computed coefficients C[i].

B:  [out] The coefficients of the original S with respect to the
computed basis. If B is a non-null pointer and B matches the
dimensions of B, then the coefficients computed during the
orthogonalization routine will be stored in B, similar to calling If B
points to a Teuchos::SerialDenseMatrix with size inconsistent with S,
then a std::invalid_argument exception will be thrown.  Otherwise, if
B is null, the caller will not have access to the computed
coefficients.  The normalization uses classical Gram-Schmidt
iteration, so that B is an upper triangular matrix with positive
diagonal elements.

MS:  [in/out] If specified by the user, on input MS is required to be
the image of S under the operator getOp(). On output, MS will be
updated to reflect the changes in S.

MX:  [in] If specified by the user, on MX[i] is required to be the
image of X[i] under the operator getOp().

MY:  [in] If specified by the user, on MY[i] is required to be the
image of Y[i] under the operator getOp().

If X[i] != Teuchos::null or Y[i] != Teuchos::null, then X[i] and Y[i]
are required to have the same number of columns, and each should have
the same number of rows as S.

For any i != j, $\\\\langle Y[i], X[j] \\\\rangle == 0$.

If biOrtho == true, $\\\\langle Y[i], X[i]\\\\rangle == I$

Otherwise, if biOrtho == false, then $\\\\langle Y[i], X[i]\\\\rangle$
should be Hermitian positive-definite.

If X[i] and Y[i] have $xc_i$ columns and S has $sc$ columns, then C[i]
if specified must be $xc_i \\\\times sc$.

If S has $sc$ columns, then B if specified must be $sc \\\\times sc $.

Rank of the basis computed by this method. ";

/*  Methods implementing Anasazi::MatOrthoManager  */

%feature("docstring")  Anasazi::ICGSOrthoManager::projectMat "void
Anasazi::ICGSOrthoManager< ScalarType, MV, OP >::projectMat(MV &X,
Teuchos::Array< Teuchos::RCP< const MV > > Q, Teuchos::Array<
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > > >
C=Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix< int,
ScalarType > >(Teuchos::null)), Teuchos::RCP< MV > MX=Teuchos::null,
Teuchos::Array< Teuchos::RCP< const MV > >
MQ=Teuchos::tuple(Teuchos::RCP< const MV >(Teuchos::null))) const

Given a list of mutually orthogonal and internally orthonormal bases
Q, this method projects a multivector X onto the space orthogonal to
the individual Q[i], optionally returning the coefficients of X for
the individual Q[i]. All of this is done with respect to the inner
product innerProd().

This method calls projectGen() as follows: See projectGen() for
argument requirements. ";

%feature("docstring")  Anasazi::ICGSOrthoManager::normalizeMat "int
Anasazi::ICGSOrthoManager< ScalarType, MV, OP >::normalizeMat(MV &X,
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > >
B=Teuchos::null, Teuchos::RCP< MV > MX=Teuchos::null) const

This method takes a multivector X and attempts to compute an
orthonormal basis for $colspan(X)$, with respect to innerProd().

This method calls projectAndNormalizeGen() as follows: See
projectAndNormalizeGen() for argument requirements. ";

%feature("docstring")
Anasazi::ICGSOrthoManager::projectAndNormalizeMat "int
Anasazi::ICGSOrthoManager< ScalarType, MV, OP
>::projectAndNormalizeMat(MV &X, Teuchos::Array< Teuchos::RCP< const
MV > > Q, Teuchos::Array< Teuchos::RCP< Teuchos::SerialDenseMatrix<
int, ScalarType > > > C=Teuchos::tuple(Teuchos::RCP<
Teuchos::SerialDenseMatrix< int, ScalarType > >(Teuchos::null)),
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > >
B=Teuchos::null, Teuchos::RCP< MV > MX=Teuchos::null, Teuchos::Array<
Teuchos::RCP< const MV > > MQ=Teuchos::tuple(Teuchos::RCP< const MV
>(Teuchos::null))) const

Given a set of bases Q[i] and a multivector X, this method computes an
orthonormal basis for $colspan(X) - \\\\sum_i colspan(Q[i])$.

This method calls projectAndNormalizeGen() as follows: See
projectAndNormalizeGen() for argument requirements. ";

/*  Error methods  */

%feature("docstring")  Anasazi::ICGSOrthoManager::orthonormErrorMat "Teuchos::ScalarTraits< ScalarType >::magnitudeType
Anasazi::ICGSOrthoManager< ScalarType, MV, OP
>::orthonormErrorMat(const MV &X, Teuchos::RCP< const MV >
MX=Teuchos::null) const

This method computes the error in orthonormality of a multivector,
measured as the Frobenius norm of the difference innerProd(X,Y) - I.
The method has the option of exploiting a caller-provided MX. ";

%feature("docstring")  Anasazi::ICGSOrthoManager::orthogErrorMat "Teuchos::ScalarTraits< ScalarType >::magnitudeType
Anasazi::ICGSOrthoManager< ScalarType, MV, OP >::orthogErrorMat(const
MV &X1, const MV &X2, Teuchos::RCP< const MV > MX1, Teuchos::RCP<
const MV > MX2) const

This method computes the error in orthogonality of two multivectors,
measured as the Frobenius norm of innerProd(X,Y). The method has the
option of exploiting a caller-provided MX. ";

/*  Accessor routines  */

%feature("docstring")  Anasazi::ICGSOrthoManager::setNumIters "void
Anasazi::ICGSOrthoManager< ScalarType, MV, OP >::setNumIters(int
numIters)

Set parameter for re-orthogonalization threshold. ";

%feature("docstring")  Anasazi::ICGSOrthoManager::getNumIters "int
Anasazi::ICGSOrthoManager< ScalarType, MV, OP >::getNumIters() const

Return parameter for re-orthogonalization threshold. ";


// File: classAnasazi_1_1IRTR.xml
%feature("docstring") Anasazi::IRTR "

IRTR is a caching implementation of the Implicit Riemannian Trust-
Region (IRTR) eigensolver.

The solver uses between 10 and 13 blocks of vectors, compared to the
requirements by SIRTR of 6 to 8 blocks of vectors. The base
requirement is 10 blocks of vectors, where a block of vectors contains
a number of vectors equal to the block size specified for the solver
(see RTRBase::getBlockSize()). Additional blocks are required when
solving a generalized eigenvalue problem or when using a
preconditioiner.

For more information, see RTRBase.

Chris Baker

C++ includes: AnasaziIRTR.hpp ";

/*  Constructor/Destructor  */

%feature("docstring")  Anasazi::IRTR::IRTR "Anasazi::IRTR<
ScalarType, MV, OP >::IRTR(const Teuchos::RCP< Eigenproblem<
ScalarType, MV, OP > > &problem, const Teuchos::RCP< SortManager<
typename Teuchos::ScalarTraits< ScalarType >::magnitudeType > >
&sorter, const Teuchos::RCP< OutputManager< ScalarType > > &printer,
const Teuchos::RCP< StatusTest< ScalarType, MV, OP > > &tester, const
Teuchos::RCP< GenOrthoManager< ScalarType, MV, OP > > &ortho,
Teuchos::ParameterList &params)

IRTR constructor with eigenproblem, solver utilities, and parameter
list of solver options.

This constructor takes pointers required by the eigensolver, in
addition to a parameter list of options for the eigensolver. These
options include the following: \"Rho Prime\" - an MagnitudeType
specifying the size of the implicit trust-region radius.

\"Block Size\" - an int specifying the block size used by the
algorithm. This can also be specified using the setBlockSize() method.

\"Leftmost\" - a bool specifying whether the solver is computing the
leftmost (\"SR\") or rightmost (\"LR\") eigenvalues. Default: true.
This must be in accord with the SortManager pass to the constructor.

\"Kappa Convergence\" - a MagnitudeType specifing the rate of
convergence for the linear convergence regime. Default: 0.1

\"Theta Convergence\" - a MagnitudeType specifing the order of
convergence for the linear convergence regime. theta implies a
convergence order of theta+1. Default: 1.0 ";

%feature("docstring")  Anasazi::IRTR::~IRTR "virtual Anasazi::IRTR<
ScalarType, MV, OP >::~IRTR()

IRTR destructor ";

/*  Solver methods  */

%feature("docstring")  Anasazi::IRTR::iterate "void Anasazi::IRTR<
ScalarType, MV, OP >::iterate()

Impemements Eigensolver. The outer IRTR iteration. See
RTRBase::iterate(). ";

/*  Output methods  */

%feature("docstring")  Anasazi::IRTR::currentStatus "void
Anasazi::IRTR< ScalarType, MV, OP >::currentStatus(std::ostream &os)

Impemements Eigensolver. This method requests that the solver print
out its current status to screen. ";


// File: classAnasazi_1_1LOBPCG.xml
%feature("docstring") Anasazi::LOBPCG "

This class provides the Locally Optimal Block Preconditioned Conjugate
Gradient (LOBPCG) iteration, a preconditioned iteration for solving
linear Hermitian eigenproblems.

This implementation is a modification of the one found in A. Knyazev,
\"Toward the optimal preconditioned eigensolver:         Locally
optimal block preconditioner conjugate gradient method\", SIAM J. Sci.
Comput., vol 23, n 2, pp. 517-541.

The modification consists of the orthogonalization steps recommended
in U. Hetmaniuk and R. Lehoucq, \"Basis Selection in LOBPCG\", Journal
of Computational Physics.

These modifcation are referred to as full orthogonalization, and
consist of also conducting the local optimization using an orthonormal
basis.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist

C++ includes: AnasaziLOBPCG.hpp ";

/*  Constructor/Destructor  */

%feature("docstring")  Anasazi::LOBPCG::LOBPCG "Anasazi::LOBPCG<
ScalarType, MV, OP >::LOBPCG(const Teuchos::RCP< Eigenproblem<
ScalarType, MV, OP > > &problem, const Teuchos::RCP< SortManager<
typename Teuchos::ScalarTraits< ScalarType >::magnitudeType > >
&sorter, const Teuchos::RCP< OutputManager< ScalarType > > &printer,
const Teuchos::RCP< StatusTest< ScalarType, MV, OP > > &tester, const
Teuchos::RCP< MatOrthoManager< ScalarType, MV, OP > > &ortho,
Teuchos::ParameterList &params)

LOBPCG constructor with eigenproblem, solver utilities, and parameter
list of solver options.

This constructor takes pointers required by the eigensolver, in
addition to a parameter list of options for the eigensolver. These
options include the following: \"Block Size\" - an int specifying the
block size used by the algorithm. This can also be specified using the
setBlockSize() method.

\"Full Ortho\" - a bool specifying whether the solver should employ a
full orthogonalization technique. This can also be specified using the
setFullOrtho() method. ";

%feature("docstring")  Anasazi::LOBPCG::~LOBPCG "virtual
Anasazi::LOBPCG< ScalarType, MV, OP >::~LOBPCG()

LOBPCG destructor ";

/*  Solver methods  */

%feature("docstring")  Anasazi::LOBPCG::iterate "void
Anasazi::LOBPCG< ScalarType, MV, OP >::iterate()

This method performs LOBPCG iterations until the status test indicates
the need to stop or an error occurs (in which case, an exception is
thrown).

iterate() will first determine whether the solver is initialized; if
not, it will call initialize() using default arguments. After
initialization, the solver performs LOBPCG iterations until the status
test evaluates as Passed, at which point the method returns to the
caller.

The LOBPCG iteration proceeds as follows: The current residual (R) is
preconditioned to form H

H is orthogonalized against the auxiliary vectors and, if full
orthogonalization  is enabled, against X and P.

The basis [X H P] is used to project the problem matrices.

The projected eigenproblem is solved, and the desired eigenvectors and
eigenvalues are selected.

These are used to form the new eigenvector estimates (X) and the
search directions (P).  If full orthogonalization is enabled, these
are generated to be mutually orthogonal and with orthonormal columns.

The new residual (R) is formed.

The status test is queried at the beginning of the iteration.

Possible exceptions thrown include std::logic_error,
std::invalid_argument or one of the LOBPCG-specific exceptions. ";

%feature("docstring")  Anasazi::LOBPCG::initialize "void
Anasazi::LOBPCG< ScalarType, MV, OP >::initialize(LOBPCGState<
ScalarType, MV > newstate)

Initialize the solver to an iterate, optionally providing the Ritz
values, residual, and search direction.

LOBPCGState contains fields V, KV and MV: These are ignored by
initialize()  The LOBPCG eigensolver contains a certain amount of
state relating to the current iterate, including the current residual,
the current search direction, and the images of these spaces under the
eigenproblem operators.

initialize() gives the user the opportunity to manually set these,
although this must be done with caution, abiding by the rules given
below. All notions of orthogonality and orthonormality are derived
from the inner product specified by the orthogonalization manager.

isInitialized() == true (see post-conditions of isInitialize())

If newstate.P != Teuchos::null, hasP() == true.  Otherwise, hasP() ==
false

The user has the option of specifying any component of the state using
initialize(). However, these arguments are assumed to match the post-
conditions specified under isInitialized(). Any component of the state
(i.e., KX) not given to initialize() will be generated. ";

%feature("docstring")  Anasazi::LOBPCG::initialize "void
Anasazi::LOBPCG< ScalarType, MV, OP >::initialize()

Initialize the solver with the initial vectors from the eigenproblem
or random data. ";

%feature("docstring")  Anasazi::LOBPCG::isInitialized "bool
Anasazi::LOBPCG< ScalarType, MV, OP >::isInitialized() const

Indicates whether the solver has been initialized or not.

bool indicating the state of the solver.

If isInitialized() == true: X is orthogonal to auxiliary vectors and
has orthonormal columns

KX == Op*X

MX == M*X if M != Teuchos::null  Otherwise, MX == Teuchos::null

getRitzValues() returns the sorted Ritz values with respect to X

getResNorms(), getRes2Norms(), getRitzResNorms() are correct

If hasP() == true, P orthogonal to auxiliary vectors

If getFullOrtho() == true, P is orthogonal to X and has orthonormal
columns

KP == Op*P

MP == M*P if M != Teuchos::null  Otherwise, MP == Teuchos::null ";

%feature("docstring")  Anasazi::LOBPCG::getState "LOBPCGState<
ScalarType, MV > Anasazi::LOBPCG< ScalarType, MV, OP >::getState()
const

Get the current state of the eigensolver.

The data is only valid if isInitialized() == true. The data for the
search directions P is only meaningful if hasP() == true. Finally, the
data for the preconditioned residual (H) is only meaningful in the
situation where the solver throws an ::LOBPCGRitzFailure exception
during iterate().

An LOBPCGState object containing const views to the current solver
state. ";

/*  Status methods  */

%feature("docstring")  Anasazi::LOBPCG::getNumIters "int
Anasazi::LOBPCG< ScalarType, MV, OP >::getNumIters() const

Get the current iteration count. ";

%feature("docstring")  Anasazi::LOBPCG::resetNumIters "void
Anasazi::LOBPCG< ScalarType, MV, OP >::resetNumIters()

Reset the iteration count. ";

%feature("docstring")  Anasazi::LOBPCG::getRitzVectors "Teuchos::RCP<
const MV > Anasazi::LOBPCG< ScalarType, MV, OP >::getRitzVectors()

Get the Ritz vectors from the previous iteration.

A multivector with getBlockSize() vectors containing the sorted Ritz
vectors corresponding to the most significant Ritz values. The i-th
vector of the return corresponds to the i-th Ritz vector; there is no
need to use getRitzIndex(). ";

%feature("docstring")  Anasazi::LOBPCG::getRitzValues "std::vector<
Value< ScalarType > > Anasazi::LOBPCG< ScalarType, MV, OP
>::getRitzValues()

Get the Ritz values from the previous iteration.

A vector of length getCurSubspaceDim() containing the Ritz values from
the previous projected eigensolve. ";

%feature("docstring")  Anasazi::LOBPCG::getRitzIndex "std::vector<
int > Anasazi::LOBPCG< ScalarType, MV, OP >::getRitzIndex()

Get the index used for extracting Ritz vectors from getRitzVectors().

Because BlockDavidson is a Hermitian solver, all Ritz values are real
and all Ritz vectors can be represented in a single column of a
multivector. Therefore, getRitzIndex() is not needed when using the
output from getRitzVectors().

An int vector of size getCurSubspaceDim() composed of zeros. ";

%feature("docstring")  Anasazi::LOBPCG::getResNorms "std::vector<
typename Teuchos::ScalarTraits< ScalarType >::magnitudeType >
Anasazi::LOBPCG< ScalarType, MV, OP >::getResNorms()

Get the current residual norms.

A vector of length getBlockSize() containing the norms of the
residuals, with respect to the orthogonalization manager norm()
method. ";

%feature("docstring")  Anasazi::LOBPCG::getRes2Norms "std::vector<
typename Teuchos::ScalarTraits< ScalarType >::magnitudeType >
Anasazi::LOBPCG< ScalarType, MV, OP >::getRes2Norms()

Get the current residual 2-norms.

A vector of length getBlockSize() containing the 2-norms of the
residuals. ";

%feature("docstring")  Anasazi::LOBPCG::getRitzRes2Norms "std::vector< typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType > Anasazi::LOBPCG< ScalarType, MV, OP
>::getRitzRes2Norms()

Get the 2-norms of the residuals.

The Ritz residuals are not defined for the LOBPCG iteration. Hence,
this method returns the 2-norms of the direct residuals, and is
equivalent to calling getRes2Norms().

A vector of length getBlockSize() containing the 2-norms of the direct
residuals. ";

%feature("docstring")  Anasazi::LOBPCG::getCurSubspaceDim "int
Anasazi::LOBPCG< ScalarType, MV, OP >::getCurSubspaceDim() const

Get the dimension of the search subspace used to generate the current
eigenvectors and eigenvalues.

LOBPCG employs a sequential subspace iteration, maintaining a fixed-
rank basis, as opposed to an expanding subspace mechanism employed by
Krylov-subspace solvers like BlockKrylovSchur and BlockDavidson.

An integer specifying the rank of the subspace generated by the
eigensolver. If isInitialized() == false, the return is 0. Otherwise,
the return will be 2*getBlockSize() or 3*getBlockSize(). ";

%feature("docstring")  Anasazi::LOBPCG::getMaxSubspaceDim "int
Anasazi::LOBPCG< ScalarType, MV, OP >::getMaxSubspaceDim() const

Get the maximum dimension allocated for the search subspace. For
LOBPCG, this always returns 3*getBlockSize(), the dimension of the
subspace colspan([X H P]). ";

/*  Accessor routines from Eigensolver  */

%feature("docstring")  Anasazi::LOBPCG::setStatusTest "void
Anasazi::LOBPCG< ScalarType, MV, OP >::setStatusTest(Teuchos::RCP<
StatusTest< ScalarType, MV, OP > > test)

Set a new StatusTest for the solver. ";

%feature("docstring")  Anasazi::LOBPCG::getStatusTest "Teuchos::RCP<
StatusTest< ScalarType, MV, OP > > Anasazi::LOBPCG< ScalarType, MV, OP
>::getStatusTest() const

Get the current StatusTest used by the solver. ";

%feature("docstring")  Anasazi::LOBPCG::getProblem "const
Eigenproblem< ScalarType, MV, OP > & Anasazi::LOBPCG< ScalarType, MV,
OP >::getProblem() const

Get a constant reference to the eigenvalue problem. ";

%feature("docstring")  Anasazi::LOBPCG::setBlockSize "void
Anasazi::LOBPCG< ScalarType, MV, OP >::setBlockSize(int blockSize)

Set the blocksize to be used by the iterative solver in solving this
eigenproblem.

If the block size is reduced, then the new iterate (and residual and
search direction) are chosen as the subset of the current iterate
preferred by the sort manager. Otherwise, the solver state is set to
uninitialized. ";

%feature("docstring")  Anasazi::LOBPCG::getBlockSize "int
Anasazi::LOBPCG< ScalarType, MV, OP >::getBlockSize() const

Get the blocksize to be used by the iterative solver in solving this
eigenproblem. ";

%feature("docstring")  Anasazi::LOBPCG::setAuxVecs "void
Anasazi::LOBPCG< ScalarType, MV, OP >::setAuxVecs(const
Teuchos::Array< Teuchos::RCP< const MV > > &auxvecs)

Set the auxiliary vectors for the solver.

Because the current iterate X and search direction P cannot be assumed
orthogonal to the new auxiliary vectors, a call to setAuxVecs() with a
non-empty argument will reset the solver to the uninitialized state.

In order to preserve the current state, the user will need to extract
it from the solver using getState(), orthogonalize it against the new
auxiliary vectors, and manually reinitialize the solver using
initialize(). ";

%feature("docstring")  Anasazi::LOBPCG::getAuxVecs "Teuchos::Array<
Teuchos::RCP< const MV > > Anasazi::LOBPCG< ScalarType, MV, OP
>::getAuxVecs() const

Get the current auxiliary vectors. ";

/*  %LOBPCG-specific accessor routines  */

%feature("docstring")  Anasazi::LOBPCG::setFullOrtho "void
Anasazi::LOBPCG< ScalarType, MV, OP >::setFullOrtho(bool fullOrtho)

Instruct the LOBPCG iteration to use full orthogonality.

If the getFullOrtho() == false and isInitialized() == true and hasP()
== true, then P will be invalidated by setting full orthogonalization
to true. ";

%feature("docstring")  Anasazi::LOBPCG::getFullOrtho "bool
Anasazi::LOBPCG< ScalarType, MV, OP >::getFullOrtho() const

Determine if the LOBPCG iteration is using full orthogonality. ";

%feature("docstring")  Anasazi::LOBPCG::hasP "bool Anasazi::LOBPCG<
ScalarType, MV, OP >::hasP()

Indicates whether the search direction given by getState() is valid.
";

/*  Output methods  */

%feature("docstring")  Anasazi::LOBPCG::currentStatus "void
Anasazi::LOBPCG< ScalarType, MV, OP >::currentStatus(std::ostream &os)

This method requests that the solver print out its current status to
screen. ";


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
return ::Failed; if it returns ::Passed, solve() will throw an
AnasaziError exception.

Much of this behavior is controlled via parameters and options passed
to the solver manager. For more information, see LOBPCGSolMgr().

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist

C++ includes: AnasaziLOBPCGSolMgr.hpp ";

/*  Constructors/Destructor  */

%feature("docstring")  Anasazi::LOBPCGSolMgr::LOBPCGSolMgr "Anasazi::LOBPCGSolMgr< ScalarType, MV, OP >::LOBPCGSolMgr(const
Teuchos::RCP< Eigenproblem< ScalarType, MV, OP > > &problem,
Teuchos::ParameterList &pl)

Basic constructor for LOBPCGSolMgr.

This constructor accepts the Eigenproblem to be solved in addition to
a parameter list of options for the solver manager. These options
include the following: Solver parameters  \"Which\" - a string
specifying the desired eigenvalues: SM, LM, SR or LR. Default: \"SR\"

\"Block Size\" - a int specifying the block size to be used by the
underlying LOBPCG solver. Default: problem->getNEV()

\"Full Ortho\" - a bool specifying whether the underlying solver
should employ the full orthogonalization scheme. Default: true

\"Recover\" - a bool specifying whether the solver manager should
attempt to recover in the case of a LOBPCGRitzFailure when full
orthogonalization is disabled. Default: true

\"Verbosity\" - a sum of MsgType specifying the verbosity. Default:
::Errors

\"Init\" - a LOBPCGState<ScalarType,MV> struct used to initialize the
LOBPCG eigensolver.

Convergence parameters (if using default convergence test; see
setGlobalStatusTest())  \"Maximum Iterations\" - a int specifying the
maximum number of iterations the underlying solver is allowed to
perform. Default: 100

\"Convergence Tolerance\" - a MagnitudeType specifying the level that
residual norms must reach to decide convergence. Default: machine
precision.

\"Relative Convergence Tolerance\" - a bool specifying whether
residuals norms should be scaled by their eigenvalues for the
purposing of deciding convergence. Default: true

\"Convergence Norm\" - a string specifying the norm for convergence
testing: \"2\" or \"M\"

Locking parameters (if using default locking test; see
setLockingStatusTest())  \"Use Locking\" - a bool specifying whether
the algorithm should employ locking of converged eigenpairs. Default:
false

\"Max Locked\" - a int specifying the maximum number of eigenpairs to
be locked. Default: problem->getNEV()

\"Locking Quorum\" - a int specifying the number of eigenpairs that
must meet the locking criteria before locking actually occurs.
Default: 1

\"Locking Tolerance\" - a MagnitudeType specifying the level that
residual norms must reach to decide locking. Default: 0.1*convergence
tolerance

\"Relative Locking Tolerance\" - a bool specifying whether residuals
norms should be scaled by their eigenvalues for the purposing of
deciding locking. Default: true

\"Locking Norm\" - a string specifying the norm for locking testing:
\"2\" or \"M\" ";

%feature("docstring")  Anasazi::LOBPCGSolMgr::~LOBPCGSolMgr "virtual
Anasazi::LOBPCGSolMgr< ScalarType, MV, OP >::~LOBPCGSolMgr()

Destructor. ";

/*  Accessor methods  */

%feature("docstring")  Anasazi::LOBPCGSolMgr::getProblem "const
Eigenproblem<ScalarType,MV,OP>& Anasazi::LOBPCGSolMgr< ScalarType, MV,
OP >::getProblem() const

Return the eigenvalue problem. ";

%feature("docstring")  Anasazi::LOBPCGSolMgr::getNumIters "int
Anasazi::LOBPCGSolMgr< ScalarType, MV, OP >::getNumIters() const

Get the iteration count for the most recent call to  solve(). ";

%feature("docstring")  Anasazi::LOBPCGSolMgr::getTimers "Teuchos::Array<Teuchos::RCP<Teuchos::Time> > Anasazi::LOBPCGSolMgr<
ScalarType, MV, OP >::getTimers() const

Return the timers for this object.

The timers are ordered as follows: time spent in solve() routine

time spent locking converged eigenvectors ";

/*  Solver application methods  */

%feature("docstring")  Anasazi::LOBPCGSolMgr::solve "ReturnType
Anasazi::LOBPCGSolMgr< ScalarType, MV, OP >::solve()

This method performs possibly repeated calls to the underlying
eigensolver's iterate() routine until the problem has been solved (as
decided by the solver manager) or the solver manager decides to quit.

This method calls LOBPCG::iterate(), which will return either because
a specially constructed status test evaluates to ::Passed or an
exception is thrown.

A return from LOBPCG::iterate() signifies one of the following
scenarios: the maximum number of iterations has been exceeded. In this
scenario, the solver manager will place  all converged eigenpairs into
the eigenproblem and return ::Unconverged.

the locking conditions have been met. In this scenario, some of the
current eigenpairs will be removed  from the eigensolver and placed
into auxiliary storage. The eigensolver will be restarted with the
remaining  eigenpairs and some random information to replace the
removed eigenpairs.

global convergence has been met. In this case, the most significant
NEV eigenpairs in the solver and locked storage  have met the
convergence criterion. (Here, NEV refers to the number of eigenpairs
requested by the Eigenproblem.)  In this scenario, the solver manager
will return ::Converged.

an LOBPCGRitzFailure exception has been thrown. If full
orthogonalization is enabled and recovery from this exception  is
requested, the solver manager will attempt to recover from this
exception by gathering the current eigenvectors,  preconditioned
residual, and search directions from the eigensolver,
orthogonormalizing the basis composed of these  three, projecting the
eigenproblem, and restarting the eigensolver with the solution of the
project eigenproblem. Any  additional failure that occurs during this
recovery effort will result in the eigensolver returning
::Unconverged.

::ReturnType specifying: ::Converged: the eigenproblem was solved to
the specification required by the solver manager.

::Unconverged: the eigenproblem was not solved to the specification
desired by the solver manager ";

%feature("docstring")  Anasazi::LOBPCGSolMgr::setGlobalStatusTest "void Anasazi::LOBPCGSolMgr< ScalarType, MV, OP
>::setGlobalStatusTest(const Teuchos::RCP< StatusTest< ScalarType, MV,
OP > > &global)

Set the status test defining global convergence. ";

%feature("docstring")  Anasazi::LOBPCGSolMgr::getGlobalStatusTest "const Teuchos::RCP< StatusTest< ScalarType, MV, OP > > &
Anasazi::LOBPCGSolMgr< ScalarType, MV, OP >::getGlobalStatusTest()
const

Get the status test defining global convergence. ";

%feature("docstring")  Anasazi::LOBPCGSolMgr::setLockingStatusTest "void Anasazi::LOBPCGSolMgr< ScalarType, MV, OP
>::setLockingStatusTest(const Teuchos::RCP< StatusTest< ScalarType,
MV, OP > > &locking)

Set the status test defining locking. ";

%feature("docstring")  Anasazi::LOBPCGSolMgr::getLockingStatusTest "const Teuchos::RCP< StatusTest< ScalarType, MV, OP > > &
Anasazi::LOBPCGSolMgr< ScalarType, MV, OP >::getLockingStatusTest()
const

Get the status test defining locking. ";

%feature("docstring")  Anasazi::LOBPCGSolMgr::setDebugStatusTest "void Anasazi::LOBPCGSolMgr< ScalarType, MV, OP
>::setDebugStatusTest(const Teuchos::RCP< StatusTest< ScalarType, MV,
OP > > &debug)

Set the status test for debugging. ";

%feature("docstring")  Anasazi::LOBPCGSolMgr::getDebugStatusTest "const Teuchos::RCP< StatusTest< ScalarType, MV, OP > > &
Anasazi::LOBPCGSolMgr< ScalarType, MV, OP >::getDebugStatusTest()
const

Get the status test for debugging. ";


// File: structAnasazi_1_1LOBPCGState.xml
%feature("docstring") Anasazi::LOBPCGState "

Structure to contain pointers to Anasazi state variables.

This struct is utilized by LOBPCG::initialize() and
LOBPCG::getState().

C++ includes: AnasaziLOBPCG.hpp ";

%feature("docstring")  Anasazi::LOBPCGState::LOBPCGState "Anasazi::LOBPCGState< ScalarType, MultiVector >::LOBPCGState() ";


// File: classAnasazi_1_1MakePairOp.xml
%feature("docstring") Anasazi::MakePairOp "";


// File: classAnasazi_1_1MatOrthoManager.xml
%feature("docstring") Anasazi::MatOrthoManager "

Anasazi's templated virtual class for providing routines for
orthogonalization and orthonormalization of multivectors using matrix-
based inner products.

This class extends Anasazi::OrthoManager by providing extra calling
arguments to orthogonalization routines, to reduce the cost of
applying the inner product in cases where the user already has the
image of target multivectors under the inner product matrix.

A concrete implementation of this class is necessary. The user can
create their own implementation if those supplied are not suitable for
their needs.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziMatOrthoManager.hpp ";

/*  Constructor/Destructor  */

%feature("docstring")  Anasazi::MatOrthoManager::MatOrthoManager "Anasazi::MatOrthoManager< ScalarType, MV, OP
>::MatOrthoManager(Teuchos::RCP< const OP > Op=Teuchos::null)

Default constructor. ";

%feature("docstring")  Anasazi::MatOrthoManager::~MatOrthoManager "virtual Anasazi::MatOrthoManager< ScalarType, MV, OP
>::~MatOrthoManager()

Destructor. ";

/*  Accessor routines  */

%feature("docstring")  Anasazi::MatOrthoManager::setOp "void
Anasazi::MatOrthoManager< ScalarType, MV, OP >::setOp(Teuchos::RCP<
const OP > Op)

Set operator used for inner product. ";

%feature("docstring")  Anasazi::MatOrthoManager::getOp "Teuchos::RCP<
const OP > Anasazi::MatOrthoManager< ScalarType, MV, OP >::getOp()
const

Get operator used for inner product. ";

%feature("docstring")  Anasazi::MatOrthoManager::getOpCounter "int
Anasazi::MatOrthoManager< ScalarType, MV, OP >::getOpCounter() const

Retrieve operator counter.

This counter returns the number of applications of the operator
specifying the inner product. When the operator is applied to a
multivector, the counter is incremented by the number of vectors in
the multivector. If the operator is not specified, the counter is
never incremented. ";

%feature("docstring")  Anasazi::MatOrthoManager::resetOpCounter "void
Anasazi::MatOrthoManager< ScalarType, MV, OP >::resetOpCounter()

Reset the operator counter to zero.

See getOpCounter() for more details. ";

/*  Matrix-based Orthogonality Methods  */

%feature("docstring")  Anasazi::MatOrthoManager::innerProdMat "void
Anasazi::MatOrthoManager< ScalarType, MV, OP >::innerProdMat(const MV
&X, const MV &Y, Teuchos::SerialDenseMatrix< int, ScalarType > &Z,
Teuchos::RCP< const MV > MX=Teuchos::null, Teuchos::RCP< const MV >
MY=Teuchos::null) const

Provides a matrix-based inner product.

Provides the inner product \\\\[ \\\\langle x, y \\\\rangle = x^H M y
\\\\] Optionally allows the provision of $M y$ and/or $M x$. See
OrthoManager::innerProd() for more details. ";

%feature("docstring")  Anasazi::MatOrthoManager::normMat "void
Anasazi::MatOrthoManager< ScalarType, MV, OP >::normMat(const MV &X,
std::vector< typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType > &normvec, Teuchos::RCP< const MV >
MX=Teuchos::null) const

Provides the norm induced by the matrix-based inner product.

Provides the norm: \\\\[ \\\\|x\\\\|_M = \\\\sqrt{x^H M y} \\\\]
Optionally allows the provision of $M x$. See OrthoManager::norm() for
more details. ";

%feature("docstring")  Anasazi::MatOrthoManager::projectMat "virtual
void Anasazi::MatOrthoManager< ScalarType, MV, OP >::projectMat(MV &X,
Teuchos::Array< Teuchos::RCP< const MV > > Q, Teuchos::Array<
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > > >
C=Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix< int,
ScalarType > >(Teuchos::null)), Teuchos::RCP< MV > MX=Teuchos::null,
Teuchos::Array< Teuchos::RCP< const MV > >
MQ=Teuchos::tuple(Teuchos::RCP< const MV >(Teuchos::null))) const =0

Provides matrix-based projection method.

This method optionally allows the provision of $M X$ and/or the $M
Q[i]$. See OrthoManager::project() for more details.

Parameters:
-----------

X:  Q:  C:  [in/out] As in OrthoManager::project()

MX:  [in/out] If specified by the user, on input MX is required to be
the image of X under the operator getOp(). On output, MX will be
updated to reflect the changes in X.

MQ:  [in] If specified by the user, on MQ[i] is required to be the
image of Q[i] under the operator getOp(). ";

%feature("docstring")  Anasazi::MatOrthoManager::normalizeMat "virtual int Anasazi::MatOrthoManager< ScalarType, MV, OP
>::normalizeMat(MV &X, Teuchos::RCP< Teuchos::SerialDenseMatrix< int,
ScalarType > > B=Teuchos::null, Teuchos::RCP< MV > MX=Teuchos::null)
const =0

Provides matrix-based orthonormalization method.

This method optionally allows the provision of $M X$. See
orthoManager::normalize() for more details.

Parameters:
-----------

X:  B:  [in/out] As in OrthoManager::normalize()

MX:  [in/out] If specified by the user, on input MX is required to be
the image of X under the operator getOp(). On output, MX will be
updated to reflect the changes in X.

Rank of the basis computed by this method, less than or equal to the
number of columns in X. This specifies how many columns in the
returned X and MX and rows in the returned B are valid. ";

%feature("docstring")
Anasazi::MatOrthoManager::projectAndNormalizeMat "virtual int
Anasazi::MatOrthoManager< ScalarType, MV, OP
>::projectAndNormalizeMat(MV &X, Teuchos::Array< Teuchos::RCP< const
MV > > Q, Teuchos::Array< Teuchos::RCP< Teuchos::SerialDenseMatrix<
int, ScalarType > > > C=Teuchos::tuple(Teuchos::RCP<
Teuchos::SerialDenseMatrix< int, ScalarType > >(Teuchos::null)),
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > >
B=Teuchos::null, Teuchos::RCP< MV > MX=Teuchos::null, Teuchos::Array<
Teuchos::RCP< const MV > > MQ=Teuchos::tuple(Teuchos::RCP< const MV
>(Teuchos::null))) const =0

Provides matrix-based projection/orthonormalization method.

This method optionally allows the provision of $M X$ and/or the $M
Q[i]$. See orthoManager::projectAndNormalize() for more details.

Parameters:
-----------

X:  Q:  C:  B:  [in/out] As in OrthoManager::projectAndNormalize()

MX:  [in/out] If specified by the user, on input MX is required to be
the image of X under the operator getOp(). On output, MX will be
updated to reflect the changes in X.

MQ:  [in] If specified by the user, on MQ[i] is required to be the
image of Q[i] under the operator getOp().

Rank of the basis computed by this method, less than or equal to the
number of columns in X. This specifies how many columns in the
returned X and MX and rows in the returned B are valid. ";

%feature("docstring")  Anasazi::MatOrthoManager::orthonormErrorMat "virtual Teuchos::ScalarTraits<ScalarType>::magnitudeType
Anasazi::MatOrthoManager< ScalarType, MV, OP
>::orthonormErrorMat(const MV &X, Teuchos::RCP< const MV >
MX=Teuchos::null) const =0

This method computes the error in orthonormality of a multivector.

This method optionally allows optionally exploits a caller-provided
MX. ";

%feature("docstring")  Anasazi::MatOrthoManager::orthogErrorMat "virtual Teuchos::ScalarTraits<ScalarType>::magnitudeType
Anasazi::MatOrthoManager< ScalarType, MV, OP >::orthogErrorMat(const
MV &X, const MV &Y, Teuchos::RCP< const MV > MX=Teuchos::null,
Teuchos::RCP< const MV > MY=Teuchos::null) const =0

This method computes the error in orthogonality of two multivectors.

This method optionally allows optionally exploits a caller-provided MX
and/or MY. ";

/*  Methods implementing Anasazi::OrthoManager  */

%feature("docstring")  Anasazi::MatOrthoManager::innerProd "void
Anasazi::MatOrthoManager< ScalarType, MV, OP >::innerProd(const MV &X,
const MV &Y, Teuchos::SerialDenseMatrix< int, ScalarType > &Z) const

Implements the interface OrthoManager::innerProd().

This method calls ";

%feature("docstring")  Anasazi::MatOrthoManager::norm "void
Anasazi::MatOrthoManager< ScalarType, MV, OP >::norm(const MV &X,
std::vector< typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType > &normvec) const

Implements the interface OrthoManager::norm().

This method calls ";

%feature("docstring")  Anasazi::MatOrthoManager::project "void
Anasazi::MatOrthoManager< ScalarType, MV, OP >::project(MV &X,
Teuchos::Array< Teuchos::RCP< const MV > > Q, Teuchos::Array<
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > > >
C=Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix< int,
ScalarType > >(Teuchos::null))) const

Implements the interface OrthoManager::project().

This method calls ";

%feature("docstring")  Anasazi::MatOrthoManager::normalize "int
Anasazi::MatOrthoManager< ScalarType, MV, OP >::normalize(MV &X,
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > >
B=Teuchos::null) const

Implements the interface OrthoManager::normalize().

This method calls ";

%feature("docstring")  Anasazi::MatOrthoManager::projectAndNormalize "int Anasazi::MatOrthoManager< ScalarType, MV, OP
>::projectAndNormalize(MV &X, Teuchos::Array< Teuchos::RCP< const MV >
> Q, Teuchos::Array< Teuchos::RCP< Teuchos::SerialDenseMatrix< int,
ScalarType > > > C=Teuchos::tuple(Teuchos::RCP<
Teuchos::SerialDenseMatrix< int, ScalarType > >(Teuchos::null)),
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > >
B=Teuchos::null) const

Implements the interface OrthoManager::projectAndNormalize().

This method calls ";

%feature("docstring")  Anasazi::MatOrthoManager::orthonormError "Teuchos::ScalarTraits< ScalarType >::magnitudeType
Anasazi::MatOrthoManager< ScalarType, MV, OP >::orthonormError(const
MV &X) const

Implements the interface OrthoManager::orthonormError().

This method calls ";

%feature("docstring")  Anasazi::MatOrthoManager::orthogError "Teuchos::ScalarTraits< ScalarType >::magnitudeType
Anasazi::MatOrthoManager< ScalarType, MV, OP >::orthogError(const MV
&X1, const MV &X2) const

Implements the interface OrthoManager::orthogError().

This method calls ";


// File: classAnasazi_1_1MultiVec.xml
%feature("docstring") Anasazi::MultiVec "

Interface for multivectors used by Anasazi's linear solvers.

Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

Parameters:
-----------

ScalarType:  The type of entries of the multivector.

Anasazi accesses multivectors through a traits interface called
MultiVecTraits. If you want to use Anasazi with your own multivector
class MV, you may either specialize MultiVecTraits for MV, or you may
wrap MV in your own class that implements MultiVec. Specializing
MultiVecTraits works via compile-time polymorphism, whereas
implementing the MultiVec interface works via run-time polymorphism.
You may pick whichever option you like. However, specializing
MultiVecTraits is the preferred method. This is because Anasazi's
linear solvers always use a specialization of MultiVecTraits to access
multivector operations. They only use MultiVec through a
specialization of the MultiVecTraits traits class, which is
implemented below in this header file.

If you want your multivector class (or a wrapper thereof) to implement
the MultiVec interface, you should inherit from MultiVec<ScalarType>,
where ScalarType is the type of entries in the multivector. For
example, a multivector with entries of type double would inherit from
MultiVec<double>.

C++ includes: AnasaziMultiVec.hpp ";

/*  Constructor/Destructor  */

%feature("docstring")  Anasazi::MultiVec::MultiVec "Anasazi::MultiVec< ScalarType >::MultiVec()

Default constructor. ";

%feature("docstring")  Anasazi::MultiVec::~MultiVec "virtual
Anasazi::MultiVec< ScalarType >::~MultiVec()

Destructor (virtual for memory safety of derived classes). ";

/*  Creation methods  */

%feature("docstring")  Anasazi::MultiVec::Clone "virtual
MultiVec<ScalarType>* Anasazi::MultiVec< ScalarType >::Clone(const int
numvecs) const =0

Create a new MultiVec with numvecs columns.

Pointer to the new multivector with uninitialized values. ";

%feature("docstring")  Anasazi::MultiVec::CloneCopy "virtual
MultiVec<ScalarType>* Anasazi::MultiVec< ScalarType >::CloneCopy()
const =0

Create a new MultiVec and copy contents of *this into it (deep copy).

Pointer to the new multivector ";

%feature("docstring")  Anasazi::MultiVec::CloneCopy "virtual
MultiVec<ScalarType>* Anasazi::MultiVec< ScalarType >::CloneCopy(const
std::vector< int > &index) const =0

Creates a new Anasazi::MultiVec and copies the selected contents of
*this into the new vector (deep copy). The copied vectors from *this
are indicated by the index.size() indices in index.

Pointer to the new multivector ";

%feature("docstring")  Anasazi::MultiVec::CloneViewNonConst "virtual
MultiVec<ScalarType>* Anasazi::MultiVec< ScalarType
>::CloneViewNonConst(const std::vector< int > &index)=0

Creates a new Anasazi::MultiVec that shares the selected contents of
*this. The index of the numvecs vectors shallow copied from *this are
indicated by the indices given in index.

Pointer to the new multivector ";

%feature("docstring")  Anasazi::MultiVec::CloneView "virtual const
MultiVec<ScalarType>* Anasazi::MultiVec< ScalarType >::CloneView(const
std::vector< int > &index) const =0

Creates a new Anasazi::MultiVec that shares the selected contents of
*this. The index of the numvecs vectors shallow copied from *this are
indicated by the indices given in index.

Pointer to the new multivector ";

/*  Dimension information methods  */

%feature("docstring")  Anasazi::MultiVec::GetVecLength "virtual int
Anasazi::MultiVec< ScalarType >::GetVecLength() const =0

The number of rows in the multivector. ";

%feature("docstring")  Anasazi::MultiVec::GetNumberVecs "virtual int
Anasazi::MultiVec< ScalarType >::GetNumberVecs() const =0

The number of vectors (i.e., columns) in the multivector. ";

/*  Update methods  */

%feature("docstring")  Anasazi::MultiVec::MvTimesMatAddMv "virtual
void Anasazi::MultiVec< ScalarType >::MvTimesMatAddMv(ScalarType
alpha, const MultiVec< ScalarType > &A, const
Teuchos::SerialDenseMatrix< int, ScalarType > &B, ScalarType beta)=0

Update *this with alpha * A * B + beta * ( *this). ";

%feature("docstring")  Anasazi::MultiVec::MvAddMv "virtual void
Anasazi::MultiVec< ScalarType >::MvAddMv(ScalarType alpha, const
MultiVec< ScalarType > &A, ScalarType beta, const MultiVec< ScalarType
> &B)=0

Replace *this with alpha * A + beta * B. ";

%feature("docstring")  Anasazi::MultiVec::MvScale "virtual void
Anasazi::MultiVec< ScalarType >::MvScale(ScalarType alpha)=0

Scale each element of the vectors in *this with alpha. ";

%feature("docstring")  Anasazi::MultiVec::MvScale "virtual void
Anasazi::MultiVec< ScalarType >::MvScale(const std::vector< ScalarType
> &alpha)=0

Scale each element of the i-th vector in *this with alpha[i]. ";

%feature("docstring")  Anasazi::MultiVec::MvTransMv "virtual void
Anasazi::MultiVec< ScalarType >::MvTransMv(ScalarType alpha, const
MultiVec< ScalarType > &A, Teuchos::SerialDenseMatrix< int, ScalarType
> &B) const =0

Compute a dense matrix B through the matrix-matrix multiply alpha *
A^T * ( *this). ";

%feature("docstring")  Anasazi::MultiVec::MvDot "virtual void
Anasazi::MultiVec< ScalarType >::MvDot(const MultiVec< ScalarType >
&A, std::vector< ScalarType > &b) const =0

Compute the dot product of each column of *this with the corresponding
column of A.

Compute a vector b whose entries are the individual dot-products. That
is, b[i] = A[i]^H * (*this)[i] where A[i] is the i-th column of A. ";

/*  Norm method  */

%feature("docstring")  Anasazi::MultiVec::MvNorm "virtual void
Anasazi::MultiVec< ScalarType >::MvNorm(std::vector< typename
Teuchos::ScalarTraits< ScalarType >::magnitudeType > &normvec) const
=0

Compute the 2-norm of each vector in *this.

Parameters:
-----------

normvec:  [out] On output, normvec[i] holds the 2-norm of the i-th
vector of *this. ";

/*  Initialization methods  */

%feature("docstring")  Anasazi::MultiVec::SetBlock "virtual void
Anasazi::MultiVec< ScalarType >::SetBlock(const MultiVec< ScalarType >
&A, const std::vector< int > &index)=0

Copy the vectors in A to a set of vectors in *this.

The numvecs vectors in A are copied to a subset of vectors in *this
indicated by the indices given in index. ";

%feature("docstring")  Anasazi::MultiVec::MvRandom "virtual void
Anasazi::MultiVec< ScalarType >::MvRandom()=0

Fill all the vectors in *this with random numbers. ";

%feature("docstring")  Anasazi::MultiVec::MvInit "virtual void
Anasazi::MultiVec< ScalarType >::MvInit(ScalarType alpha)=0

Replace each element of the vectors in *this with alpha. ";

/*  Print method  */

%feature("docstring")  Anasazi::MultiVec::MvPrint "virtual void
Anasazi::MultiVec< ScalarType >::MvPrint(std::ostream &os) const =0

Print *this multivector to the os output stream. ";


// File: classAnasazi_1_1MultiVecTraits.xml
%feature("docstring") Anasazi::MultiVecTraits "

Traits class which defines basic operations on multivectors.

Parameters:
-----------

ScalarType:  The type of the entries in the multivectors.

MV:  The type of the multivectors themselves.

This traits class tells Anasazi's solvers how to perform multivector
operations for the multivector type MV. These operations include
creating copies or views, finding the number of rows or columns (i.e.,
vectors) in a given multivector, and computing inner products, norms,
and vector sums. (Anasazi's solvers use the OperatorTraits traits
class to apply operators to multivectors.)

Anasazi gives users two different ways to tell its solvers how to
compute with multivectors of a given type MV. The first and preferred
way is for users to specialize MultiVecTraits, this traits class, for
their given MV type. Anasazi provides specializations for MV =
Epetra_MultiVector, Tpetra::MultiVector, and Thyra::MultiVectorBase.
The second way is for users to make their multivector type (or a
wrapper thereof) inherit from MultiVec. This works because Anasazi
provides a specialization of MultiVecTraits for MultiVec. Specializing
MultiVecTraits is more flexible because it does not require a
multivector type to inherit from MultiVec; this is possible even if
you do not have control over the interface of a class.

If you have a different multivector type MV that you would like to use
with Anasazi, and if that type does not inherit from MultiVec, then
you must implement a specialization of MultiVecTraits for MV.
Otherwise, this traits class will report a compile-time error
(relating to UndefinedMultiVecTraits). Specializing MultiVecTraits for
your MV type is not hard. Just look at the examples for
Epetra_MultiVector (in anasazi/epetra/src/AnasaziEpetraAdapter.hpp)
and Tpetra::MultiVector (in
anasazi/tpetra/src/AnasaziTpetraAdapter.hpp).

You do not need to write a specialization of MultiVecTraits if you are
using Epetra, Tpetra, or Thyra multivectors. Anasazi already provides
specializations for these types. Just relax and enjoy using the
solvers!

C++ includes: AnasaziMultiVecTraits.hpp ";

/*  Creation methods  */

%feature("docstring")  Anasazi::MultiVecTraits::Clone "static
Teuchos::RCP<MV> Anasazi::MultiVecTraits< ScalarType, MV
>::Clone(const MV &mv, const int numvecs)

Creates a new empty MV containing numvecs columns.

Reference-counted pointer to the new multivector of type MV. ";

%feature("docstring")  Anasazi::MultiVecTraits::CloneCopy "static
Teuchos::RCP<MV> Anasazi::MultiVecTraits< ScalarType, MV
>::CloneCopy(const MV &mv)

Creates a new MV and copies contents of mv into the new vector (deep
copy).

Reference-counted pointer to the new multivector of type MV. ";

%feature("docstring")  Anasazi::MultiVecTraits::CloneCopy "static
Teuchos::RCP<MV> Anasazi::MultiVecTraits< ScalarType, MV
>::CloneCopy(const MV &mv, const std::vector< int > &index)

Creates a new MV and copies the selected contents of mv into the new
vector (deep copy).

The copied vectors from mv are indicated by the index.size() indices
in index. Reference-counted pointer to the new multivector of type MV.
";

%feature("docstring")  Anasazi::MultiVecTraits::CloneCopy "static
Teuchos::RCP<MV> Anasazi::MultiVecTraits< ScalarType, MV
>::CloneCopy(const MV &mv, const Teuchos::Range1D &index)

Deep copy of specified columns of mv.

Create a new MV, and copy (deep copy) the columns of mv specified by
the given inclusive index range into the new multivector.

Parameters:
-----------

mv:  [in] Multivector to copy

index:  [in] Inclusive index range of columns of mv

Reference-counted pointer to the new multivector of type MV. ";

%feature("docstring")  Anasazi::MultiVecTraits::CloneViewNonConst "static Teuchos::RCP<MV> Anasazi::MultiVecTraits< ScalarType, MV
>::CloneViewNonConst(MV &mv, const std::vector< int > &index)

Creates a new MV that shares the selected contents of mv (shallow
copy).

The index of the numvecs vectors shallow copied from mv are indicated
by the indices given in index. Reference-counted pointer to the new
multivector of type MV. ";

%feature("docstring")  Anasazi::MultiVecTraits::CloneViewNonConst "static Teuchos::RCP<MV> Anasazi::MultiVecTraits< ScalarType, MV
>::CloneViewNonConst(MV &mv, const Teuchos::Range1D &index)

Non-const view of specified columns of mv.

Return a non-const view of the columns of mv specified by the given
inclusive index range.

Parameters:
-----------

mv:  [in] Multivector to view (shallow non-const copy)

index:  [in] Inclusive index range of columns of mv

Reference-counted pointer to the non-const view of specified columns
of mv ";

%feature("docstring")  Anasazi::MultiVecTraits::CloneView "static
Teuchos::RCP<const MV> Anasazi::MultiVecTraits< ScalarType, MV
>::CloneView(const MV &mv, const std::vector< int > &index)

Creates a new const MV that shares the selected contents of mv
(shallow copy).

The index of the numvecs vectors shallow copied from mv are indicated
by the indices given in index. Reference-counted pointer to the new
const multivector of type MV. ";

%feature("docstring")  Anasazi::MultiVecTraits::CloneView "static
Teuchos::RCP<MV> Anasazi::MultiVecTraits< ScalarType, MV
>::CloneView(MV &mv, const Teuchos::Range1D &index)

Const view of specified columns of mv.

Return a const view of the columns of mv specified by the given
inclusive index range.

Parameters:
-----------

mv:  [in] Multivector to view (shallow const copy)

index:  [in] Inclusive index range of columns of mv

Reference-counted pointer to the const view of specified columns of mv
";

/*  Attribute methods  */

%feature("docstring")  Anasazi::MultiVecTraits::GetVecLength "static
int Anasazi::MultiVecTraits< ScalarType, MV >::GetVecLength(const MV
&mv)

Obtain the vector length of mv. ";

%feature("docstring")  Anasazi::MultiVecTraits::GetNumberVecs "static
int Anasazi::MultiVecTraits< ScalarType, MV >::GetNumberVecs(const MV
&mv)

Obtain the number of vectors in mv. ";

/*  Update methods  */

%feature("docstring")  Anasazi::MultiVecTraits::MvTimesMatAddMv "static void Anasazi::MultiVecTraits< ScalarType, MV
>::MvTimesMatAddMv(const ScalarType alpha, const MV &A, const
Teuchos::SerialDenseMatrix< int, ScalarType > &B, const ScalarType
beta, MV &mv)

Update mv with $ \\\\alpha AB + \\\\beta mv $. ";

%feature("docstring")  Anasazi::MultiVecTraits::MvAddMv "static void
Anasazi::MultiVecTraits< ScalarType, MV >::MvAddMv(const ScalarType
alpha, const MV &A, const ScalarType beta, const MV &B, MV &mv)

Replace mv with $\\\\alpha A + \\\\beta B$. ";

%feature("docstring")  Anasazi::MultiVecTraits::MvScale "static void
Anasazi::MultiVecTraits< ScalarType, MV >::MvScale(MV &mv, const
ScalarType alpha)

Scale each element of the vectors in mv with alpha. ";

%feature("docstring")  Anasazi::MultiVecTraits::MvScale "static void
Anasazi::MultiVecTraits< ScalarType, MV >::MvScale(MV &mv, const
std::vector< ScalarType > &alpha)

Scale each element of the i-th vector in mv with alpha[i]. ";

%feature("docstring")  Anasazi::MultiVecTraits::MvTransMv "static
void Anasazi::MultiVecTraits< ScalarType, MV >::MvTransMv(const
ScalarType alpha, const MV &A, const MV &mv,
Teuchos::SerialDenseMatrix< int, ScalarType > &B)

Compute a dense matrix B through the matrix-matrix multiply $
\\\\alpha A^Hmv $. ";

%feature("docstring")  Anasazi::MultiVecTraits::MvDot "static void
Anasazi::MultiVecTraits< ScalarType, MV >::MvDot(const MV &mv, const
MV &A, std::vector< ScalarType > &b)

Compute a vector b where the components are the individual dot-
products of the i-th columns of A and mv, i.e. $b[i] = A[i]^Hmv[i]$.
";

/*  Norm method  */

%feature("docstring")  Anasazi::MultiVecTraits::MvNorm "static void
Anasazi::MultiVecTraits< ScalarType, MV >::MvNorm(const MV &mv,
std::vector< typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType > &normvec)

Compute the 2-norm of each individual vector of mv. Upon return,
normvec[i] holds the value of $||mv_i||_2$, the i-th column of mv. ";

/*  Initialization methods  */

%feature("docstring")  Anasazi::MultiVecTraits::SetBlock "static void
Anasazi::MultiVecTraits< ScalarType, MV >::SetBlock(const MV &A, const
std::vector< int > &index, MV &mv)

Copy the vectors in A to a set of vectors in mv indicated by the
indices given in index.

The numvecs vectors in A are copied to a subset of vectors in mv
indicated by the indices given in index, i.e.  mv[index[i]] = A[i]. ";

%feature("docstring")  Anasazi::MultiVecTraits::SetBlock "static void
Anasazi::MultiVecTraits< ScalarType, MV >::SetBlock(const MV &A, const
Teuchos::Range1D &index, MV &mv)

Deep copy of A into specified columns of mv.

(Deeply) copy the first index.size() columns of A into the columns of
mv specified by the given index range.

Postcondition: mv[i] = A[i - index.lbound()] for all i in
[index.lbound(), index.ubound()]

Parameters:
-----------

A:  [in] Source multivector

index:  [in] Inclusive index range of columns of mv; index set of the
target

mv:  [out] Target multivector ";

%feature("docstring")  Anasazi::MultiVecTraits::Assign "static void
Anasazi::MultiVecTraits< ScalarType, MV >::Assign(const MV &A, MV &mv)

mv := A

Assign (deep copy) A into mv. ";

%feature("docstring")  Anasazi::MultiVecTraits::MvRandom "static void
Anasazi::MultiVecTraits< ScalarType, MV >::MvRandom(MV &mv)

Replace the vectors in mv with random vectors. ";

%feature("docstring")  Anasazi::MultiVecTraits::MvInit "static void
Anasazi::MultiVecTraits< ScalarType, MV >::MvInit(MV &mv, const
ScalarType alpha=Teuchos::ScalarTraits< ScalarType >::zero())

Replace each element of the vectors in mv with alpha. ";

/*  Print method  */

%feature("docstring")  Anasazi::MultiVecTraits::MvPrint "static void
Anasazi::MultiVecTraits< ScalarType, MV >::MvPrint(const MV &mv,
std::ostream &os)

Print the mv multi-vector to the os output stream. ";


// File: classAnasazi_1_1MultiVecTraits_3_01ScalarType_00_01MultiVec_3_01ScalarType_01_4_01_4.xml
%feature("docstring") Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > > "

Specialization of MultiVecTraits for Belos::MultiVec.

Anasazi interfaces to every multivector implementation through a
specialization of MultiVecTraits. Thus, we provide a specialization of
MultiVecTraits for the MultiVec run-time polymorphic interface above.

Parameters:
-----------

ScalarType:  The type of entries in the multivector; the template
parameter of MultiVec.

C++ includes: AnasaziMultiVec.hpp ";

/*  Creation methods  */

%feature("docstring")  Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::Clone " static Teuchos::RCP<MultiVec<ScalarType> >
Anasazi::MultiVecTraits< ScalarType, MultiVec< ScalarType >
>::Clone(const MultiVec< ScalarType > &mv, const int numvecs)

Create a new empty  MultiVec containing numvecs columns.

Reference-counted pointer to the new  MultiVec. ";

%feature("docstring")  Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::CloneCopy " static Teuchos::RCP<MultiVec<ScalarType> >
Anasazi::MultiVecTraits< ScalarType, MultiVec< ScalarType >
>::CloneCopy(const MultiVec< ScalarType > &mv)

Creates a new  Anasazi::MultiVec and copies contents of mv into the
new vector (deep copy).

Reference-counted pointer to the new  Anasazi::MultiVec. ";

%feature("docstring")  Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::CloneCopy " static Teuchos::RCP<MultiVec<ScalarType> >
Anasazi::MultiVecTraits< ScalarType, MultiVec< ScalarType >
>::CloneCopy(const MultiVec< ScalarType > &mv, const std::vector< int
> &index)

Creates a new  Anasazi::MultiVec and copies the selected contents of
mv into the new vector (deep copy).

The copied vectors from mv are indicated by the index.size() indices
in index. Reference-counted pointer to the new  Anasazi::MultiVec. ";

%feature("docstring")  Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::CloneViewNonConst " static
Teuchos::RCP<MultiVec<ScalarType> > Anasazi::MultiVecTraits<
ScalarType, MultiVec< ScalarType > >::CloneViewNonConst(MultiVec<
ScalarType > &mv, const std::vector< int > &index)

Creates a new  Anasazi::MultiVec that shares the selected contents of
mv (shallow copy).

The index of the numvecs vectors shallow copied from mv are indicated
by the indices given in index. Reference-counted pointer to the new
Anasazi::MultiVec. ";

%feature("docstring")  Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::CloneView " static Teuchos::RCP<const
MultiVec<ScalarType> > Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::CloneView(const MultiVec< ScalarType > &mv, const
std::vector< int > &index)

Creates a new const  Anasazi::MultiVec that shares the selected
contents of mv (shallow copy).

The index of the numvecs vectors shallow copied from mv are indicated
by the indices given in index. Reference-counted pointer to the new
const  Anasazi::MultiVec. ";

/*  Attribute methods  */

%feature("docstring")  Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::GetVecLength " static int Anasazi::MultiVecTraits<
ScalarType, MultiVec< ScalarType > >::GetVecLength(const MultiVec<
ScalarType > &mv)

Obtain the vector length of mv. ";

%feature("docstring")  Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::GetNumberVecs " static int Anasazi::MultiVecTraits<
ScalarType, MultiVec< ScalarType > >::GetNumberVecs(const MultiVec<
ScalarType > &mv)

Obtain the number of vectors in mv. ";

/*  Update methods  */

%feature("docstring")  Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::MvTimesMatAddMv " static void Anasazi::MultiVecTraits<
ScalarType, MultiVec< ScalarType > >::MvTimesMatAddMv(ScalarType
alpha, const MultiVec< ScalarType > &A, const
Teuchos::SerialDenseMatrix< int, ScalarType > &B, ScalarType beta,
MultiVec< ScalarType > &mv)

Update mv with $ \\\\alpha AB + \\\\beta mv $. ";

%feature("docstring")  Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::MvAddMv " static void Anasazi::MultiVecTraits<
ScalarType, MultiVec< ScalarType > >::MvAddMv(ScalarType alpha, const
MultiVec< ScalarType > &A, ScalarType beta, const MultiVec< ScalarType
> &B, MultiVec< ScalarType > &mv)

Replace mv with $\\\\alpha A + \\\\beta B$. ";

%feature("docstring")  Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::MvTransMv " static void Anasazi::MultiVecTraits<
ScalarType, MultiVec< ScalarType > >::MvTransMv(ScalarType alpha,
const MultiVec< ScalarType > &A, const MultiVec< ScalarType > &mv,
Teuchos::SerialDenseMatrix< int, ScalarType > &B)

Compute a dense matrix B through the matrix-matrix multiply $
\\\\alpha A^Tmv $. ";

%feature("docstring")  Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::MvDot " static void Anasazi::MultiVecTraits<
ScalarType, MultiVec< ScalarType > >::MvDot(const MultiVec< ScalarType
> &mv, const MultiVec< ScalarType > &A, std::vector< ScalarType > &b)

Compute a vector b where the components are the individual dot-
products of the i-th columns of A and mv, i.e. $b[i] = A[i]^H mv[i]$.
";

%feature("docstring")  Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::MvScale " static void Anasazi::MultiVecTraits<
ScalarType, MultiVec< ScalarType > >::MvScale(MultiVec< ScalarType >
&mv, ScalarType alpha)

Scale each element of the vectors in *this with alpha. ";

%feature("docstring")  Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::MvScale " static void Anasazi::MultiVecTraits<
ScalarType, MultiVec< ScalarType > >::MvScale(MultiVec< ScalarType >
&mv, const std::vector< ScalarType > &alpha)

Scale each element of the i-th vector in *this with alpha[i]. ";

/*  Norm method  */

%feature("docstring")  Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::MvNorm " static void Anasazi::MultiVecTraits<
ScalarType, MultiVec< ScalarType > >::MvNorm(const MultiVec<
ScalarType > &mv, std::vector< typename Teuchos::ScalarTraits<
ScalarType >::magnitudeType > &normvec)

Compute the 2-norm of each individual vector of mv. Upon return,
normvec[i] holds the value of $||mv_i||_2$, the i-th column of mv. ";

/*  Initialization methods  */

%feature("docstring")  Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::SetBlock " static void Anasazi::MultiVecTraits<
ScalarType, MultiVec< ScalarType > >::SetBlock(const MultiVec<
ScalarType > &A, const std::vector< int > &index, MultiVec< ScalarType
> &mv)

Copy the vectors in A to a set of vectors in mv indicated by the
indices given in index.

The numvecs vectors in A are copied to a subset of vectors in mv
indicated by the indices given in index, i.e.  mv[index[i]] = A[i]. ";

%feature("docstring")  Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::MvRandom " static void Anasazi::MultiVecTraits<
ScalarType, MultiVec< ScalarType > >::MvRandom(MultiVec< ScalarType >
&mv)

Replace the vectors in mv with random vectors. ";

%feature("docstring")  Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::MvInit " static void Anasazi::MultiVecTraits<
ScalarType, MultiVec< ScalarType > >::MvInit(MultiVec< ScalarType >
&mv, ScalarType alpha=Teuchos::ScalarTraits< ScalarType >::zero())

Replace each element of the vectors in mv with alpha. ";

/*  Print method  */

%feature("docstring")  Anasazi::MultiVecTraits< ScalarType, MultiVec<
ScalarType > >::MvPrint " static void Anasazi::MultiVecTraits<
ScalarType, MultiVec< ScalarType > >::MvPrint(const MultiVec<
ScalarType > &mv, std::ostream &os)

Print the mv multi-vector to the os output stream. ";


// File: classAnasazi_1_1details_1_1MultiVecTsqrAdapter.xml
%feature("docstring") Anasazi::details::MultiVecTsqrAdapter "

TSQR adapter for MultiVec.

TSQR (Tall Skinny QR factorization) is an orthogonalization kernel
that is as accurate as Householder QR, yet requires only $2 \\\\log P$
messages between $P$ MPI processes, independently of the number of
columns in the multivector.

TSQR works independently of the particular multivector implementation,
and interfaces to the latter via an adapter class. Each multivector
type MV needs its own adapter class. The specialization of
MultiVecTraits for MV refers to its corresponding adapter class as its
tsqr_adaptor_type [sic; sorry about the lack of standard spelling of
\"adapter\"] typedef.

This class is the TSQR adapter for MultiVec. It merely calls
MultiVec's corresponding methods for TSQR functionality.

C++ includes: AnasaziMultiVec.hpp ";

%feature("docstring")
Anasazi::details::MultiVecTsqrAdapter::factorExplicit "void
Anasazi::details::MultiVecTsqrAdapter< ScalarType >::factorExplicit(MV
&A, MV &Q, dense_matrix_type &R, const bool
forceNonnegativeDiagonal=false)

Compute QR factorization A = QR, using TSQR. ";

%feature("docstring")
Anasazi::details::MultiVecTsqrAdapter::revealRank "int
Anasazi::details::MultiVecTsqrAdapter< ScalarType >::revealRank(MV &Q,
dense_matrix_type &R, const magnitude_type &tol)

Compute rank-revealing decomposition using results of
factorExplicit(). ";


// File: classAnasazi_1_1Operator.xml
%feature("docstring") Anasazi::Operator "

Anasazi's templated virtual class for constructing an operator that
can interface with the OperatorTraits class used by the eigensolvers.

A concrete implementation of this class is necessary. The user can
create their own implementation if those supplied are not suitable for
their needs.

Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziOperator.hpp ";

/*  Constructor/Destructor  */

%feature("docstring")  Anasazi::Operator::Operator "Anasazi::Operator< ScalarType >::Operator()

Default constructor. ";

%feature("docstring")  Anasazi::Operator::~Operator "virtual
Anasazi::Operator< ScalarType >::~Operator()

Destructor. ";

/*  Operator application method  */

%feature("docstring")  Anasazi::Operator::Apply "virtual void
Anasazi::Operator< ScalarType >::Apply(const MultiVec< ScalarType >
&x, MultiVec< ScalarType > &y) const =0

This method takes the Anasazi::MultiVec x and applies the operator to
it resulting in the Anasazi::MultiVec y. ";


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

/*  Operator application method.  */

%feature("docstring")  Anasazi::OperatorTraits::Apply "static void
Anasazi::OperatorTraits< ScalarType, MV, OP >::Apply(const OP &Op,
const MV &x, MV &y)

Application method which performs operation y = Op*x. An OperatorError
exception is thrown if there is an error. ";


// File: classAnasazi_1_1OperatorTraits_3_01ScalarType_00_01MultiVec_3_01ScalarType_01_4_00_01Operator_3_01ScalarType_01_4_01_4.xml
%feature("docstring") Anasazi::OperatorTraits< ScalarType, MultiVec<
ScalarType >, Operator< ScalarType > > "

Template specialization of Anasazi::OperatorTraits class using
Anasazi::Operator and Anasazi::MultiVec virtual base classes.

Any class that inherits from Anasazi::Operator will be accepted by the
Anasazi templated solvers due to this interface to the
Anasazi::OperatorTraits class.

C++ includes: AnasaziOperator.hpp ";

/*  Operator application method  */

%feature("docstring")  Anasazi::OperatorTraits< ScalarType, MultiVec<
ScalarType >, Operator< ScalarType > >::Apply " static void
Anasazi::OperatorTraits< ScalarType, MultiVec< ScalarType >, Operator<
ScalarType > >::Apply(const Operator< ScalarType > &Op, const
MultiVec< ScalarType > &x, MultiVec< ScalarType > &y)

This method takes the Anasazi::MultiVec x and applies the
Anasazi::Operator Op to it resulting in the Anasazi::MultiVec y. ";


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

/*  Constructor/Destructor  */

%feature("docstring")  Anasazi::OrthoManager::OrthoManager "Anasazi::OrthoManager< ScalarType, MV >::OrthoManager()

Default constructor. ";

%feature("docstring")  Anasazi::OrthoManager::~OrthoManager "virtual
Anasazi::OrthoManager< ScalarType, MV >::~OrthoManager()

Destructor. ";

/*  Orthogonalization methods  */

%feature("docstring")  Anasazi::OrthoManager::innerProd "virtual void
Anasazi::OrthoManager< ScalarType, MV >::innerProd(const MV &X, const
MV &Y, Teuchos::SerialDenseMatrix< int, ScalarType > &Z) const =0

Provides the inner product defining the orthogonality concepts.

All concepts of orthogonality discussed in this class are defined with
respect to this inner product.

This is potentially different from MultiVecTraits::MvTransMv(). For
example, it is customary in many eigensolvers to exploit a mass matrix
M for the inner product: $x^HMx$.

Parameters:
-----------

Z:  [out] Z(i,j) contains the inner product of X[i] and Y[i]: \\\\[
Z(i,j) = \\\\langle X[i], Y[i] \\\\rangle \\\\] ";

%feature("docstring")  Anasazi::OrthoManager::norm "virtual void
Anasazi::OrthoManager< ScalarType, MV >::norm(const MV &X,
std::vector< typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType > &normvec) const =0

Provides the norm induced by innerProd().

This computes the norm for each column of a multivector. This is the
norm induced by innerProd(): \\\\[ \\\\|x\\\\| = \\\\sqrt{\\\\langle
x, x \\\\rangle} \\\\]

Parameters:
-----------

normvec:  [out] Vector of norms, whose i-th entry corresponds to the
i-th column of X

normvec.size() == GetNumberVecs(X) ";

%feature("docstring")  Anasazi::OrthoManager::project "virtual void
Anasazi::OrthoManager< ScalarType, MV >::project(MV &X,
Teuchos::Array< Teuchos::RCP< const MV > > Q, Teuchos::Array<
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > > >
C=Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix< int,
ScalarType > >(Teuchos::null))) const =0

Given a list of mutually orthogonal and internally orthonormal bases
Q, this method projects a multivector X onto the space orthogonal to
the individual Q[i], optionally returning the coefficients of X for
the individual Q[i]. All of this is done with respect to the inner
product innerProd().

After calling this routine, X will be orthogonal to each of the Q[i].

Parameters:
-----------

X:  [in/out] The multivector to be modified.  On output, the columns
of X will be orthogonal to each Q[i], satisfying \\\\[ \\\\langle
Q[i], X_{out} \\\\rangle = 0 \\\\] Also, \\\\[ X_{out} = X_{in} -
\\\\sum_i Q[i] \\\\langle Q[i], X_{in} \\\\rangle \\\\]

Q:  [in] A list of multivector bases specifying the subspaces to be
orthogonalized against, satisfying \\\\[ \\\\langle Q[i], Q[j]
\\\\rangle = I \\\\quad\\\\textrm{if}\\\\quad i=j \\\\] and \\\\[
\\\\langle Q[i], Q[j] \\\\rangle = 0 \\\\quad\\\\textrm{if}\\\\quad i
\\\\neq j\\\\ . \\\\]

C:  [out] The coefficients of X in the bases Q[i]. If C[i] is a non-
null pointer and C[i] matches the dimensions of X and Q[i], then the
coefficients computed during the orthogonalization routine will be
stored in the matrix C[i], similar to calling If C[i] points to a
Teuchos::SerialDenseMatrix with size inconsistent with X and  Q[i],
then a std::invalid_argument exception will be thrown.  Otherwise, if
C.size() < i or C[i] is a null pointer, the caller will not have
access to the computed coefficients. ";

%feature("docstring")  Anasazi::OrthoManager::normalize "virtual int
Anasazi::OrthoManager< ScalarType, MV >::normalize(MV &X,
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > >
B=Teuchos::null) const =0

This method takes a multivector X and attempts to compute a basis for
$colspan(X)$. This basis is orthonormal with respect to innerProd().

This routine returns an integer rank stating the rank of the computed
basis. If X does not have full rank and the normalize() routine does
not attempt to augment the subspace, then rank may be smaller than the
number of columns in X. In this case, only the first rank columns of
output X and first rank rows of B will be valid.

Parameters:
-----------

X:  [in/out] The multivector to be modified.  On output, the first
rank columns of X satisfy \\\\[ \\\\langle X[i], X[j] \\\\rangle =
\\\\delta_{ij}\\\\ . \\\\] Also, \\\\[ X_{in}(1:m,1:n) =
X_{out}(1:m,1:rank) B(1:rank,1:n)\\\\ , \\\\] where m is the number of
rows in X and n is the number of columns in X.

B:  [out] The coefficients of the original X with respect to the
computed basis. If B is a non-null pointer and B matches the
dimensions of B, then the coefficients computed during the
orthogonalization routine will be stored in B, similar to calling If B
points to a Teuchos::SerialDenseMatrix with size inconsistent with X,
then a std::invalid_argument exception will be thrown.  Otherwise, if
B is null, the caller will not have access to the computed
coefficients.

This matrix is not necessarily triangular (as in a QR factorization);
see the documentation of specific orthogonalization managers.

Rank of the basis computed by this method, less than or equal to the
number of columns in X. This specifies how many columns in the
returned X and rows in the returned B are valid. ";

%feature("docstring")  Anasazi::OrthoManager::projectAndNormalize "virtual int Anasazi::OrthoManager< ScalarType, MV
>::projectAndNormalize(MV &X, Teuchos::Array< Teuchos::RCP< const MV >
> Q, Teuchos::Array< Teuchos::RCP< Teuchos::SerialDenseMatrix< int,
ScalarType > > > C=Teuchos::tuple(Teuchos::RCP<
Teuchos::SerialDenseMatrix< int, ScalarType > >(Teuchos::null)),
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > >
B=Teuchos::null) const =0

Given a set of bases Q[i] and a multivector X, this method computes an
orthonormal basis for $colspan(X) - \\\\sum_i colspan(Q[i])$.

This routine returns an integer rank stating the rank of the computed
basis. If the subspace $colspan(X) - \\\\sum_i colspan(Q[i])$ does not
have dimension as large as the number of columns of X and the
orthogonalization manager does not attempt to augment the subspace,
then rank may be smaller than the number of columns of X. In this
case, only the first rank columns of output X and first rank rows of B
will be valid.

This routine guarantees both the orthogonality of the returned basis
against the Q[i] as well as the orthonormality of the returned basis.
Therefore, this method is not necessarily equivalent to calling
project() followed by a call to normalize(); see the documentation for
specific orthogonalization managers.

Parameters:
-----------

X:  [in/out] On output, the first rank columns of X satisfy \\\\[
\\\\langle X[i], X[j] \\\\rangle = \\\\delta_{ij} \\\\quad
\\\\textrm{and} \\\\quad \\\\langle X, Q[i] \\\\rangle = 0\\\\ . \\\\]
Also, \\\\[ X_{in}(1:m,1:n) = X_{out}(1:m,1:rank) B(1:rank,1:n) +
\\\\sum_i Q[i] C[i] \\\\] where m is the number of rows in X and n is
the number of columns in X.

Q:  [in] A list of multivector bases specifying the subspaces to be
orthogonalized against, satisfying \\\\[ \\\\langle Q[i], Q[j]
\\\\rangle = I \\\\quad\\\\textrm{if}\\\\quad i=j \\\\] and \\\\[
\\\\langle Q[i], Q[j] \\\\rangle = 0 \\\\quad\\\\textrm{if}\\\\quad i
\\\\neq j\\\\ . \\\\]

C:  [out] The coefficients of X in the Q[i]. If C[i] is a non-null
pointer and C[i] matches the dimensions of X and Q[i], then the
coefficients computed during the orthogonalization routine will be
stored in the matrix C[i], similar to calling If C[i] points to a
Teuchos::SerialDenseMatrix with size inconsistent with X and  Q[i],
then a std::invalid_argument exception will be thrown.  Otherwise, if
C.size() < i or C[i] is a null pointer, the caller will not have
access to the computed coefficients.

B:  [out] The coefficients of the original X with respect to the
computed basis. If B is a non-null pointer and B matches the
dimensions of B, then the coefficients computed during the
orthogonalization routine will be stored in B, similar to calling If B
points to a Teuchos::SerialDenseMatrix with size inconsistent with X,
then a std::invalid_argument exception will be thrown.  Otherwise, if
B is null, the caller will not have access to the computed
coefficients.

This matrix is not necessarily triangular (as in a QR factorization);
see the documentation of specific orthogonalization managers.

Rank of the basis computed by this method, less than or equal to the
number of columns in X. This specifies how many columns in the
returned X and rows in the returned B are valid. ";

/*  Error methods  */

%feature("docstring")  Anasazi::OrthoManager::orthonormError "virtual
Teuchos::ScalarTraits< ScalarType >::magnitudeType
Anasazi::OrthoManager< ScalarType, MV >::orthonormError(const MV &X)
const =0

This method computes the error in orthonormality of a multivector.

This method return some measure of $\\\\| \\\\langle X, X \\\\rangle -
I \\\\| $.  See the documentation of specific orthogonalization
managers. ";

%feature("docstring")  Anasazi::OrthoManager::orthogError "virtual
Teuchos::ScalarTraits<ScalarType>::magnitudeType
Anasazi::OrthoManager< ScalarType, MV >::orthogError(const MV &X1,
const MV &X2) const =0

This method computes the error in orthogonality of two multivectors.

This method return some measure of $\\\\| \\\\langle X1, X2 \\\\rangle
\\\\| $.  See the documentation of specific orthogonalization
managers. ";


// File: classAnasazi_1_1OutOfPlaceNormalizerMixin.xml
%feature("docstring") Anasazi::OutOfPlaceNormalizerMixin "

Mixin for out-of-place orthogonalization.

Mark Hoemmen  This class presents an abstract interface for multiple
inheritance (\"mixin\") for special orthogonalization methods that
normalize \"out-of-place.\" OrthoManager and MatOrthoManager both
normalize (and projectAndNormalize) multivectors \"in place,\" meaning
that the input and output multivectors are the same (X, in both
cases). Gram-Schmidt (modified or classical) is an example of an
orthogonalization method that can normalize (and projectAndNormalize)
in place. TSQR (the Tall Skinny QR factorization, see
TsqrOrthoManager.hpp for references) is an orthogonalization method
which cannot normalize (or projectAndNormalize) in place.

Tsqr(Mat) OrthoManager implements (Mat)OrthoManager's normalize() and
projectAndNormalize() methods with scratch space and copying. However,
if you handle Tsqr(Mat) OrthoManager through this mixin, you can
exploit TSQR's unique interface to avoid copying back and forth
between scratch space.

C++ includes: AnasaziTsqrOrthoManager.hpp ";

%feature("docstring")
Anasazi::OutOfPlaceNormalizerMixin::normalizeOutOfPlace "virtual int
Anasazi::OutOfPlaceNormalizerMixin< Scalar, MV
>::normalizeOutOfPlace(MV &X, MV &Q, mat_ptr B) const =0

Normalize X into Q*B.

Parameters:
-----------

X:  [in/out] On input: Multivector to normalize. On output: Possibly
overwritten with invalid values.

Q:  [out] On output: Normalized multivector.

B:  [out] On output: Normalization coefficients.

Rank of the input multivector X. ";

%feature("docstring")
Anasazi::OutOfPlaceNormalizerMixin::projectAndNormalizeOutOfPlace "virtual int Anasazi::OutOfPlaceNormalizerMixin< Scalar, MV
>::projectAndNormalizeOutOfPlace(MV &X_in, MV &X_out, Teuchos::Array<
mat_ptr > C, mat_ptr B, Teuchos::ArrayView< Teuchos::RCP< const MV > >
Q) const =0

Project and normalize X_in into X_out.

Project X_in against Q, storing projection coefficients in C, and
normalize X_in into X_out, storing normalization coefficients in B. On
output, X_out has the resulting orthogonal vectors. X_in may be
overwritten with invalid values.

Parameters:
-----------

X_in:  [in/out] On input: The vectors to project against Q and
normalize. On output: possibly overwritten with invalid values.

X_out:  [out] The normalized input vectors after projection against Q.

C:  [out] Projection coefficients

B:  [out] Normalization coefficients

Q:  [in] The orthogonal basis against which to project

Rank of X_in after projection ";

%feature("docstring")
Anasazi::OutOfPlaceNormalizerMixin::~OutOfPlaceNormalizerMixin "virtual Anasazi::OutOfPlaceNormalizerMixin< Scalar, MV
>::~OutOfPlaceNormalizerMixin()

Trivial virtual destructor, to silence compiler warnings. ";


// File: classAnasazi_1_1OutputManager.xml
%feature("docstring") Anasazi::OutputManager "

Output managers remove the need for the eigensolver to know any
information about the required output. Calling isVerbosity( MsgType
type ) informs the solver if it is supposed to output the information
corresponding to the message type.

Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziOutputManager.hpp ";

/*  Constructors/Destructor  */

%feature("docstring")  Anasazi::OutputManager::OutputManager "Anasazi::OutputManager< ScalarType >::OutputManager(int
vb=Anasazi::Errors)

Default constructor. ";

%feature("docstring")  Anasazi::OutputManager::~OutputManager "virtual Anasazi::OutputManager< ScalarType >::~OutputManager()

Destructor. ";

/*  Set/Get methods  */

%feature("docstring")  Anasazi::OutputManager::setVerbosity "virtual
void Anasazi::OutputManager< ScalarType >::setVerbosity(int vb)

Set the message output types for this manager. ";

%feature("docstring")  Anasazi::OutputManager::getVerbosity "virtual
int Anasazi::OutputManager< ScalarType >::getVerbosity() const

Get the message output types for this manager. ";

/*  Output methods  */

%feature("docstring")  Anasazi::OutputManager::isVerbosity "virtual
bool Anasazi::OutputManager< ScalarType >::isVerbosity(MsgType type)
const =0

Find out whether we need to print out information for this message
type.

This method is used by the solver to determine whether computations
are necessary for this message type. ";

%feature("docstring")  Anasazi::OutputManager::print "virtual void
Anasazi::OutputManager< ScalarType >::print(MsgType type, const
std::string output)=0

Send output to the output manager. ";

%feature("docstring")  Anasazi::OutputManager::stream "virtual
std::ostream& Anasazi::OutputManager< ScalarType >::stream(MsgType
type)=0

Create a stream for outputting to. ";

/*  Undefined methods  */


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


// File: classAnasazi_1_1RTRBase.xml
%feature("docstring") Anasazi::RTRBase "

This class is an abstract base class for Implicit Riemannian Trust-
Region based eigensolvers. The class provides the interfaces shared by
the IRTR solvers (e.g., getState() and initialize()) as well as the
shared implementations (e.g., inner products).

IRTR eigensolvers are capable of solving symmetric/Hermitian
eigenvalue problems. These solvers may be used to compute either the
leftmost (smallest real, \"SR\") or rightmost (largest real, \"LR\")
eigenvalues. For more information, see the publications at theRTR
eigensolvers page.

This class is abstract and objects cannot be instantiated. Instead,
instantiate one of the concrete derived classes: IRTR and SIRTR, the
caching and non-caching implementations of this solver. The main
difference between these solver is the memory allocated by the solvers
in support of the IRTR iteration.

The reduction in memory usage is effected by eliminating the caching
of operator applications. This also results in a reduction in vector
arithmetic required to maintain these caches. The cost is an increase
in the number of operator applications. For inexpensive operator
applications, SIRTR should provide better performance over IRTR. As
the operator applications becomes more expensive, the performance
scale tips towards the IRTR solver. Note, the trajectory of both
solvers is identical in exact arithmetic. However, the effects of
round-off error in the cached results mean that some difference
between the solvers may exist. This effect is seen when a large number
of iterations are required to solve the trust-region subproblem in
solveTRSubproblem(). Also note, the inclusion of auxiliary vectors
increases the memory requirements of these solvers linearly with the
number of auxiliary vectors. The required storage is listed in the
following table:

Number of vectors (bS == blockSize())

Solver

Base requirement

Generalized/B != null

Preconditioned

Generalized and Preconditioned

IRTR

10*bS

11*bS

12*bS

13*bS

SIRTR

6*bS

7*bS

7*bS

8*bS

Chris Baker

C++ includes: AnasaziRTRBase.hpp ";

/*  Constructor/Destructor  */

%feature("docstring")  Anasazi::RTRBase::RTRBase "Anasazi::RTRBase<
ScalarType, MV, OP >::RTRBase(const Teuchos::RCP< Eigenproblem<
ScalarType, MV, OP > > &problem, const Teuchos::RCP< SortManager<
typename Teuchos::ScalarTraits< ScalarType >::magnitudeType > >
&sorter, const Teuchos::RCP< OutputManager< ScalarType > > &printer,
const Teuchos::RCP< StatusTest< ScalarType, MV, OP > > &tester, const
Teuchos::RCP< GenOrthoManager< ScalarType, MV, OP > > &ortho,
Teuchos::ParameterList &params, const std::string &solverLabel, bool
skinnySolver)

RTRBase constructor with eigenproblem, solver utilities, and parameter
list of solver options.

The RTRBase class is abstract and cannot be instantiated; this
constructor is called by derived classes IRTR and RTR. ";

%feature("docstring")  Anasazi::RTRBase::~RTRBase "virtual
Anasazi::RTRBase< ScalarType, MV, OP >::~RTRBase()

RTRBase destructor ";

/*  Solver methods  */

%feature("docstring")  Anasazi::RTRBase::iterate "virtual void
Anasazi::RTRBase< ScalarType, MV, OP >::iterate()=0

This method performs RTR iterations until the status test indicates
the need to stop or an error occurs (in which case, an exception is
thrown).

iterate() will first determine whether the solver is initialized; if
not, it will call initialize() using default arguments. After
initialization, the solver performs RTR iterations until the status
test evaluates as ::Passed, at which point the method returns to the
caller.

The RTR iteration proceeds as follows: the trust-region subproblem at
X is solved for update Eta via a call to solveTRSubproblem()

the new iterate is the Ritz vectors with respect to X+Eta

the eigenproblem residuals are formed with respect to the new iterate

The status test is queried at the beginning of the iteration.

Possible exceptions thrown include std::logic_error,
std::invalid_argument or one of the RTR-specific exceptions. ";

%feature("docstring")  Anasazi::RTRBase::initialize "void
Anasazi::RTRBase< ScalarType, MV, OP >::initialize(RTRState<
ScalarType, MV > newstate)

Initialize the solver to an iterate, optionally providing the Ritz
values and residual.

The RTR eigensolver contains a certain amount of state relating to the
current iterate.

initialize() gives the user the opportunity to manually set these,
although this must be done with caution, abiding by the rules given
below. All notions of orthogonality and orthonormality are derived
from the inner product specified by the orthogonalization manager.

isInitialized() == true (see post-conditions of isInitialize())

The user has the option of specifying any component of the state using
initialize(). However, these arguments are assumed to match the post-
conditions specified under isInitialized(). Any component of the state
(i.e., AX) not given to initialize() will be generated.

If the Ritz values relative to newstate.X are passed in newstate.T,
then newstate.X is assume to contain Ritz vectors, i.e., newstate.T
must be B-orthonormal and it must partially diagonalize A. ";

%feature("docstring")  Anasazi::RTRBase::initialize "void
Anasazi::RTRBase< ScalarType, MV, OP >::initialize()

Initialize the solver with the initial vectors from the eigenproblem
or random data. ";

%feature("docstring")  Anasazi::RTRBase::isInitialized "bool
Anasazi::RTRBase< ScalarType, MV, OP >::isInitialized() const

Indicates whether the solver has been initialized or not.

bool indicating the state of the solver.

If isInitialized() == true: X is orthogonal to auxiliary vectors and
has orthonormal columns

AX == A*X

BX == B*X if B != Teuchos::null  Otherwise, BX == Teuchos::null

getRitzValues() returns the sorted Ritz values with respect to X

getResidualVecs() returns the residual vectors with respect to X ";

%feature("docstring")  Anasazi::RTRBase::getState "RTRState<
ScalarType, MV > Anasazi::RTRBase< ScalarType, MV, OP >::getState()
const

Get the current state of the eigensolver.

The data is only valid if isInitialized() == true.

An RTRState object containing const pointers to the current solver
state. ";

/*  Status methods  */

%feature("docstring")  Anasazi::RTRBase::getNumIters "int
Anasazi::RTRBase< ScalarType, MV, OP >::getNumIters() const

Get the current iteration count. ";

%feature("docstring")  Anasazi::RTRBase::resetNumIters "void
Anasazi::RTRBase< ScalarType, MV, OP >::resetNumIters()

Reset the iteration count. ";

%feature("docstring")  Anasazi::RTRBase::getRitzVectors "Teuchos::RCP< const MV > Anasazi::RTRBase< ScalarType, MV, OP
>::getRitzVectors()

Get the Ritz vectors from the previous iteration.

A multivector with getBlockSize() vectors containing the sorted Ritz
vectors corresponding to the most significant Ritz values. The i-th
vector of the return corresponds to the i-th Ritz vector; there is no
need to use getRitzIndex(). ";

%feature("docstring")  Anasazi::RTRBase::getRitzValues "std::vector<
Value< ScalarType > > Anasazi::RTRBase< ScalarType, MV, OP
>::getRitzValues()

Get the Ritz values from the previous iteration.

A vector of length getCurSubspaceDim() containing the Ritz values from
the previous projected eigensolve. ";

%feature("docstring")  Anasazi::RTRBase::getRitzIndex "std::vector<
int > Anasazi::RTRBase< ScalarType, MV, OP >::getRitzIndex()

Get the index used for extracting Ritz vectors from getRitzVectors().

Because BlockDavidson is a Hermitian solver, all Ritz values are real
and all Ritz vectors can be represented in a single column of a
multivector. Therefore, getRitzIndex() is not needed when using the
output from getRitzVectors().

An int vector of size getCurSubspaceDim() composed of zeros. ";

%feature("docstring")  Anasazi::RTRBase::getResNorms "std::vector<
typename Teuchos::ScalarTraits< ScalarType >::magnitudeType >
Anasazi::RTRBase< ScalarType, MV, OP >::getResNorms()

Get the current residual norms.

A vector of length getCurSubspaceDim() containing the norms of the
residuals, with respect to the orthogonalization manager norm()
method. ";

%feature("docstring")  Anasazi::RTRBase::getRes2Norms "std::vector<
typename Teuchos::ScalarTraits< ScalarType >::magnitudeType >
Anasazi::RTRBase< ScalarType, MV, OP >::getRes2Norms()

Get the current residual 2-norms.

A vector of length getCurSubspaceDim() containing the 2-norms of the
residuals. ";

%feature("docstring")  Anasazi::RTRBase::getRitzRes2Norms "std::vector< typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType > Anasazi::RTRBase< ScalarType, MV, OP
>::getRitzRes2Norms()

Get the 2-norms of the Ritz residuals.

A vector of length getCurSubspaceDim() containing the 2-norms of the
Ritz residuals. ";

%feature("docstring")  Anasazi::RTRBase::getCurSubspaceDim "int
Anasazi::RTRBase< ScalarType, MV, OP >::getCurSubspaceDim() const

Get the dimension of the search subspace used to generate the current
eigenvectors and eigenvalues.

RTR employs a sequential subspace iteration, maintaining a fixed-rank
basis, as opposed to an expanding subspace mechanism employed by
Krylov-subspace solvers like BlockKrylovSchur and BlockDavidson.

An integer specifying the rank of the subspace generated by the
eigensolver. If isInitialized() == false, the return is 0. Otherwise,
the return will be getBlockSize(). ";

%feature("docstring")  Anasazi::RTRBase::getMaxSubspaceDim "int
Anasazi::RTRBase< ScalarType, MV, OP >::getMaxSubspaceDim() const

Get the maximum dimension allocated for the search subspace. For RTR,
this always returns getBlockSize(). ";

/*  Accessor routines from Eigensolver  */

%feature("docstring")  Anasazi::RTRBase::setStatusTest "void
Anasazi::RTRBase< ScalarType, MV, OP >::setStatusTest(Teuchos::RCP<
StatusTest< ScalarType, MV, OP > > test)

Set a new StatusTest for the solver. ";

%feature("docstring")  Anasazi::RTRBase::getStatusTest "Teuchos::RCP<
StatusTest< ScalarType, MV, OP > > Anasazi::RTRBase< ScalarType, MV,
OP >::getStatusTest() const

Get the current StatusTest used by the solver. ";

%feature("docstring")  Anasazi::RTRBase::getProblem "const
Eigenproblem< ScalarType, MV, OP > & Anasazi::RTRBase< ScalarType, MV,
OP >::getProblem() const

Get a constant reference to the eigenvalue problem. ";

%feature("docstring")  Anasazi::RTRBase::setBlockSize "void
Anasazi::RTRBase< ScalarType, MV, OP >::setBlockSize(int blockSize)

Set the blocksize to be used by the iterative solver in solving this
eigenproblem.

If the block size is reduced, then the new iterate (and residual and
search direction) are chosen as the subset of the current iterate
preferred by the sort manager. Otherwise, the solver state is set to
uninitialized. ";

%feature("docstring")  Anasazi::RTRBase::getBlockSize "int
Anasazi::RTRBase< ScalarType, MV, OP >::getBlockSize() const

Get the blocksize to be used by the iterative solver in solving this
eigenproblem. ";

%feature("docstring")  Anasazi::RTRBase::setAuxVecs "void
Anasazi::RTRBase< ScalarType, MV, OP >::setAuxVecs(const
Teuchos::Array< Teuchos::RCP< const MV > > &auxvecs)

Set the auxiliary vectors for the solver.

Because the current iterate X cannot be assumed orthogonal to the new
auxiliary vectors, a call to setAuxVecs() with a non-empty argument
will reset the solver to the uninitialized state.

In order to preserve the current state, the user will need to extract
it from the solver using getState(), orthogonalize it against the new
auxiliary vectors, and manually reinitialize the solver using
initialize().

NOTE: The requirements of the IRTR solvers is such that the auxiliary
vectors must be moved into contiguous storage with the current
iterate. As a result, the multivector data in auxvecs will be copied,
and the multivectors in auxvecs will no longer be referenced. The
(unchanged) internal copies of the auxilliary vectors will be made
available to the caller by the getAuxVecs() routine. This allows the
caller to delete the caller's copies and instead use the copies owned
by the solver, avoiding the duplication of data. This is not
necessary, however. The partitioning of the auxiliary vectors passed
to setAuxVecs() will be preserved. ";

%feature("docstring")  Anasazi::RTRBase::getAuxVecs "Teuchos::Array<
Teuchos::RCP< const MV > > Anasazi::RTRBase< ScalarType, MV, OP
>::getAuxVecs() const

Get the current auxiliary vectors. ";

/*  Output methods  */

%feature("docstring")  Anasazi::RTRBase::currentStatus "void
Anasazi::RTRBase< ScalarType, MV, OP >::currentStatus(std::ostream
&os)

This method requests that the solver print out its current status to
screen. ";


// File: classAnasazi_1_1RTRInitFailure.xml
%feature("docstring") Anasazi::RTRInitFailure "

RTRInitFailure is thrown when the RTR solver is unable to generate an
initial iterate in the RTRBase::initialize() routine.

This exception is thrown from the RTRBase::initialize() method, which
is called by the user or from the RTRBase::iterate() method when
isInitialized() == false.

C++ includes: AnasaziRTRBase.hpp ";

%feature("docstring")  Anasazi::RTRInitFailure::RTRInitFailure "Anasazi::RTRInitFailure::RTRInitFailure(const std::string &what_arg)
";


// File: classAnasazi_1_1RTROrthoFailure.xml
%feature("docstring") Anasazi::RTROrthoFailure "

RTROrthoFailure is thrown when an orthogonalization attempt fails.

This is thrown in one of two scenarios. After preconditioning the
residual, the orthogonalization manager is asked to orthogonalize the
preconditioned residual (H) against the auxiliary vectors. If full
orthogonalization is enabled, H is also orthogonalized against X and P
and normalized.

The second scenario involves the generation of new X and P from the
basis [X H P]. When full orthogonalization is enabled, an attempt is
made to select coefficients for X and P so that they will be mutually
orthogonal and orthonormal.

If either of these attempts fail, the solver throws an RTROrthoFailure
exception.

C++ includes: AnasaziRTRBase.hpp ";

%feature("docstring")  Anasazi::RTROrthoFailure::RTROrthoFailure "Anasazi::RTROrthoFailure::RTROrthoFailure(const std::string &what_arg)
";


// File: classAnasazi_1_1RTRRitzFailure.xml
%feature("docstring") Anasazi::RTRRitzFailure "

RTRRitzFailure is thrown when the RTR solver is unable to continue a
call to RTRBase::iterate() due to a failure of the algorithm.

This signals that the Rayleigh-Ritz analysis of X + Eta detected ill-
conditioning of the projected mass matrix and the inability to
generate a set of orthogonal eigenvectors for the projected problem
(if thrown from iterate()) or that the analysis of the initial iterate
failed in RTRBase::initialize().

After catching this exception, the user can recover the subspace via
RTRBase::getState(). This information can be used to restart the
solver.

C++ includes: AnasaziRTRBase.hpp ";

%feature("docstring")  Anasazi::RTRRitzFailure::RTRRitzFailure "Anasazi::RTRRitzFailure::RTRRitzFailure(const std::string &what_arg)
";


// File: classAnasazi_1_1RTRSolMgr.xml
%feature("docstring") Anasazi::RTRSolMgr "

The Anasazi::RTRSolMgr provides a simple solver manager over the RTR
eigensolver. For more information, see the discussion for RTRBase.

Chris Baker

C++ includes: AnasaziRTRSolMgr.hpp ";

/*  Constructors/Destructor  */

%feature("docstring")  Anasazi::RTRSolMgr::RTRSolMgr "Anasazi::RTRSolMgr< ScalarType, MV, OP >::RTRSolMgr(const
Teuchos::RCP< Eigenproblem< ScalarType, MV, OP > > &problem,
Teuchos::ParameterList &pl)

Basic constructor for RTRSolMgr.

This constructor accepts the Eigenproblem to be solved in addition to
a parameter list of options for the solver manager. These options
include the following: Solver parameters  \"Skinny Solver\" - a bool
specifying whether a non-caching (\"skinny\") solver implementation is
used. Determines whether the underlying solver is IRTR or SIRTR.

\"Which\" - a string specifying the desired eigenvalues: SR or LR,
i.e., smallest or largest algebraic eigenvalues.

\"Block Size\" - a int specifying the block size to be used by the
underlying RTR solver. Default: problem->getNEV()

\"Verbosity\" - a sum of MsgType specifying the verbosity. Default:
::Errors

Convergence parameters  \"Maximum Iterations\" - a int specifying the
maximum number of iterations the underlying solver is allowed to
perform. Default: 100

\"Convergence Tolerance\" - a MagnitudeType specifying the level that
residual norms must reach to decide convergence. Default: machine
precision.

\"Relative Convergence Tolerance\" - a bool specifying whether
residuals norms should be scaled by their eigenvalues for the
purposing of deciding convergence. Default: true

\"Convergence Norm\" - a string specifying the norm for convergence
testing: \"2\" or \"M\" ";

%feature("docstring")  Anasazi::RTRSolMgr::~RTRSolMgr "virtual
Anasazi::RTRSolMgr< ScalarType, MV, OP >::~RTRSolMgr()

Destructor. ";

/*  Accessor methods  */

%feature("docstring")  Anasazi::RTRSolMgr::getProblem "const
Eigenproblem<ScalarType,MV,OP>& Anasazi::RTRSolMgr< ScalarType, MV, OP
>::getProblem() const

Return the eigenvalue problem. ";

%feature("docstring")  Anasazi::RTRSolMgr::getTimers "Teuchos::Array<Teuchos::RCP<Teuchos::Time> > Anasazi::RTRSolMgr<
ScalarType, MV, OP >::getTimers() const

Return the timers for this object.

The timers are ordered as follows: time spent in solve() routine ";

%feature("docstring")  Anasazi::RTRSolMgr::getNumIters "int
Anasazi::RTRSolMgr< ScalarType, MV, OP >::getNumIters() const

Get the iteration count for the most recent call to solve. ";

/*  Solver application methods  */

%feature("docstring")  Anasazi::RTRSolMgr::solve "ReturnType
Anasazi::RTRSolMgr< ScalarType, MV, OP >::solve()

This method performs possibly repeated calls to the underlying
eigensolver's iterate() routine until the problem has been solved (as
decided by the solver manager) or the solver manager decides to quit.

::ReturnType specifying: ::Converged: the eigenproblem was solved to
the specification required by the solver manager.

::Unconverged: the eigenproblem was not solved to the specification
desired by the solver manager. ";


// File: structAnasazi_1_1RTRState.xml
%feature("docstring") Anasazi::RTRState "

Structure to contain pointers to RTR state variables.

This struct is utilized by RTRBase::initialize() and
RTRBase::getState().

C++ includes: AnasaziRTRBase.hpp ";

%feature("docstring")  Anasazi::RTRState::RTRState "Anasazi::RTRState< ScalarType, MV >::RTRState() ";


// File: structAnasazi_1_1BasicSort_1_1sel1st.xml


// File: structAnasazi_1_1BasicSort_1_1sel2nd.xml


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

/*  Constructors/Destructor  */

%feature("docstring")  Anasazi::SimpleLOBPCGSolMgr::SimpleLOBPCGSolMgr
"Anasazi::SimpleLOBPCGSolMgr< ScalarType, MV, OP
>::SimpleLOBPCGSolMgr(const Teuchos::RCP< Eigenproblem< ScalarType,
MV, OP > > &problem, Teuchos::ParameterList &pl)

Basic constructor for SimpleLOBPCGSolMgr.

This constructor accepts the Eigenproblem to be solved in addition to
a parameter list of options for the solver manager. These options
include the following: \"Which\" - a string specifying the desired
eigenvalues: SM, LM, SR or LR. Default: SR

\"Block Size\" - a int specifying the block size to be used by the
underlying LOBPCG solver. Default: problem->getNEV()

\"Maximum Iterations\" - a int specifying the maximum number of
iterations the underlying solver is allowed to perform. Default: 100

\"Verbosity\" - a sum of MsgType specifying the verbosity. Default:
Anasazi::Errors

\"Convergence Tolerance\" - a MagnitudeType specifying the level that
residual norms must reach to decide convergence. Default: machine
precision ";

%feature("docstring")
Anasazi::SimpleLOBPCGSolMgr::~SimpleLOBPCGSolMgr "virtual
Anasazi::SimpleLOBPCGSolMgr< ScalarType, MV, OP
>::~SimpleLOBPCGSolMgr()

Destructor. ";

/*  Accessor methods  */

%feature("docstring")  Anasazi::SimpleLOBPCGSolMgr::getProblem "const
Eigenproblem<ScalarType,MV,OP>& Anasazi::SimpleLOBPCGSolMgr<
ScalarType, MV, OP >::getProblem() const

Return the eigenvalue problem. ";

%feature("docstring")  Anasazi::SimpleLOBPCGSolMgr::getNumIters "int
Anasazi::SimpleLOBPCGSolMgr< ScalarType, MV, OP >::getNumIters() const

Get the iteration count for the most recent call to  solve(). ";

/*  Solver application methods  */

%feature("docstring")  Anasazi::SimpleLOBPCGSolMgr::solve "ReturnType
Anasazi::SimpleLOBPCGSolMgr< ScalarType, MV, OP >::solve()

This method performs possibly repeated calls to the underlying
eigensolver's iterate() routine until the problem has been solved (as
decided by the solver manager) or the solver manager decides to quit.

::ReturnType specifying: ::Converged: the eigenproblem was solved to
the specification required by the solver manager.

::Unconverged: the eigenproblem was not solved to the specification
desired by the solver manager ";


// File: classAnasazi_1_1SIRTR.xml
%feature("docstring") Anasazi::SIRTR "

SIRTR (\"skinny IRTR\") is a non-caching, lower-memory implementation
of the Implicit Riemannian Trust-Region (IRTR) eigensolver.

The solver uses between 6 and 8 blocks of vectors, compared to the
requirements by IRTR of 10 to 13 blocks of vectors. The base
requirement is 6 blocks of vectors, where a block of vectors contains
a number of vectors equal to the block size specified for the solver
(see RTRBase::getBlockSize()). Additional blocks are required when
solving a generalized eigenvalue problem or when using a
preconditioiner.

For more information, see RTRBase.

Chris Baker

C++ includes: AnasaziSIRTR.hpp ";

/*  Constructor/Destructor  */

%feature("docstring")  Anasazi::SIRTR::SIRTR "Anasazi::SIRTR<
ScalarType, MV, OP >::SIRTR(const Teuchos::RCP< Eigenproblem<
ScalarType, MV, OP > > &problem, const Teuchos::RCP< SortManager<
typename Teuchos::ScalarTraits< ScalarType >::magnitudeType > >
&sorter, const Teuchos::RCP< OutputManager< ScalarType > > &printer,
const Teuchos::RCP< StatusTest< ScalarType, MV, OP > > &tester, const
Teuchos::RCP< GenOrthoManager< ScalarType, MV, OP > > &ortho,
Teuchos::ParameterList &params)

SIRTR constructor with eigenproblem, solver utilities, and parameter
list of solver options.

This constructor takes pointers required by the eigensolver, in
addition to a parameter list of options for the eigensolver. These
options include the following: \"Rho Prime\" - an MagnitudeType
specifying the size of the implicit trust-region radius.

\"Block Size\" - an int specifying the block size used by the
algorithm. This can also be specified using the setBlockSize() method.

\"Leftmost\" - a bool specifying whether the solver is computing the
leftmost (\"SR\") or rightmost (\"LR\") eigenvalues. Default: true.
This must be in accord with the SortManager pass to the constructor.

\"Kappa Convergence\" - a MagnitudeType specifing the rate of
convergence for the linear convergence regime. Default: 0.1

\"Theta Convergence\" - a MagnitudeType specifing the order of
convergence for the linear convergence regime. theta implies a
convergence order of theta+1. Default: 1.0 ";

%feature("docstring")  Anasazi::SIRTR::~SIRTR "virtual
Anasazi::SIRTR< ScalarType, MV, OP >::~SIRTR()

SIRTR destructor ";

/*  Solver methods  */

%feature("docstring")  Anasazi::SIRTR::iterate "void Anasazi::SIRTR<
ScalarType, MV, OP >::iterate()

Impemements Eigensolver. The outer IRTR iteration. See
RTRBase::iterate(). ";

/*  Output methods  */

%feature("docstring")  Anasazi::SIRTR::currentStatus "void
Anasazi::SIRTR< ScalarType, MV, OP >::currentStatus(std::ostream &os)

Impemements Eigensolver. This method requests that the solver print
out its current status to screen. ";


// File: classAnasazi_1_1SolverManager.xml
%feature("docstring") Anasazi::SolverManager "

The Anasazi::SolverManager is a templated virtual base class that
defines the basic interface that any solver manager will support.

C++ includes: AnasaziSolverManager.hpp ";

/*  Constructors/Destructor  */

%feature("docstring")  Anasazi::SolverManager::SolverManager "Anasazi::SolverManager< ScalarType, MV, OP >::SolverManager()

Empty constructor. ";

%feature("docstring")  Anasazi::SolverManager::~SolverManager "virtual Anasazi::SolverManager< ScalarType, MV, OP >::~SolverManager()

Destructor. ";

/*  Accessor methods  */

%feature("docstring")  Anasazi::SolverManager::getProblem "virtual
const Eigenproblem<ScalarType,MV,OP>& Anasazi::SolverManager<
ScalarType, MV, OP >::getProblem() const =0

Return the eigenvalue problem. ";

%feature("docstring")  Anasazi::SolverManager::getNumIters "virtual
int Anasazi::SolverManager< ScalarType, MV, OP >::getNumIters() const
=0

Get the iteration count for the most recent call to  solve(). ";

/*  Solver application methods  */

%feature("docstring")  Anasazi::SolverManager::solve "virtual
ReturnType Anasazi::SolverManager< ScalarType, MV, OP >::solve()=0

This method performs possibly repeated calls to the underlying
eigensolver's iterate() routine until the problem has been solved (as
decided by the solver manager) or the solver manager decides to quit.

::ReturnType specifying: ::Converged: the eigenproblem was solved to
the specification required by the solver manager.

::Unconverged: the eigenproblem was not solved to the specification
desired by the solver manager ";


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

/*  Constructor/Destructor  */

%feature("docstring")  Anasazi::SolverUtils::SolverUtils "Anasazi::SolverUtils< ScalarType, MV, OP >::SolverUtils()

Constructor. ";

%feature("docstring")  Anasazi::SolverUtils::~SolverUtils "virtual
Anasazi::SolverUtils< ScalarType, MV, OP >::~SolverUtils()

Destructor. ";

/*  Sorting Methods  */

%feature("docstring")  Anasazi::SolverUtils::permuteVectors "void
Anasazi::SolverUtils< ScalarType, MV, OP >::permuteVectors(const int
n, const std::vector< int > &perm, MV &Q, std::vector< typename
Teuchos::ScalarTraits< ScalarType >::magnitudeType > *resids=0)

Permute the vectors in a multivector according to the permutation
vector perm, and optionally the residual vector resids. ";

%feature("docstring")  Anasazi::SolverUtils::permuteVectors "void
Anasazi::SolverUtils< ScalarType, MV, OP >::permuteVectors(const
std::vector< int > &perm, Teuchos::SerialDenseMatrix< int, ScalarType
> &Q)

Permute the columns of a Teuchos::SerialDenseMatrix according to the
permutation vector perm. ";

/*  Basis update methods  */

%feature("docstring")  Anasazi::SolverUtils::applyHouse "void
Anasazi::SolverUtils< ScalarType, MV, OP >::applyHouse(int k, MV &V,
const Teuchos::SerialDenseMatrix< int, ScalarType > &H, const
std::vector< ScalarType > &tau, Teuchos::RCP< MV >
workMV=Teuchos::null)

Apply a sequence of Householder reflectors (from GEQRF) to a
multivector, using minimal workspace.

Parameters:
-----------

k:  [in] the number of Householder reflectors composing the product

V:  [in/out] the multivector to be modified, with $n$ columns

H:  [in] a $n \\\\times k$ matrix containing the encoded Householder
vectors, as returned from GEQRF (see below)

tau:  [in] the $n$ coefficients for the Householder reflects, as
returned from GEQRF

workMV:  [work] (optional) a multivector used for workspace. it need
contain only a single vector; it if contains more, only the first
vector will be modified.

This routine applies a sequence of Householder reflectors, $H_1 H_2
\\\\cdots H_k$, to a multivector $V$. The reflectors are applied
individually, as rank-one updates to the multivector. The benefit of
this is that the only required workspace is a one-column multivector.
This workspace can be provided by the user. If it is not, it will be
allocated locally on each call to applyHouse.

Each $H_i$ ( $i=1,\\\\ldots,k \\\\leq n$) has the form $ H_i = I -
\\\\tau_i v_i v_i^T $  where $\\\\tau_i$ is a scalar and $v_i$ is a
vector with $v_i(1:i-1) = 0$ and $e_i^T v_i = 1$; $v(i+1:n)$ is stored
below H(i,i) and $\\\\tau_i$ in tau[i-1]. (Note: zero-based indexing
used for data structures H and tau, while one-based indexing used for
mathematic object $v_i$).

If the multivector is $m \\\\times n$ and we apply $k$ Householder
reflectors, the total cost of the method is $4mnk - 2m(k^2-k)$ flops.
For $k=n$, this becomes $2mn^2$, the same as for a matrix-matrix
multiplication by the accumulated Householder reflectors. ";

/*  Eigensolver Projection Methods  */

%feature("docstring")  Anasazi::SolverUtils::directSolver "int
Anasazi::SolverUtils< ScalarType, MV, OP >::directSolver(int size,
const Teuchos::SerialDenseMatrix< int, ScalarType > &KK, Teuchos::RCP<
const Teuchos::SerialDenseMatrix< int, ScalarType > > MM,
Teuchos::SerialDenseMatrix< int, ScalarType > &EV, std::vector<
typename Teuchos::ScalarTraits< ScalarType >::magnitudeType > &theta,
int &nev, int esType=0)

Routine for computing the first NEV generalized eigenpairs of the
Hermitian pencil (KK, MM)

Parameters:
-----------

size:  [in] Dimension of the eigenproblem (KK, MM)

KK:  [in] Hermitian \"stiffness\" matrix

MM:  [in] Hermitian positive-definite \"mass\" matrix

EV:  [in] Dense matrix to store the nev eigenvectors

theta:  [in] Array to store the eigenvalues (Size = nev )

nev:  [in/out] Number of the smallest eigenvalues requested (in) /
computed (out)

esType:  [in] Flag to select the algorithm esType = 0 (default) Uses
LAPACK routine (Cholesky factorization of MM) with deflation of MM to
get orthonormality of eigenvectors ( $S^TMMS = I$)

esType = 1 Uses LAPACK routine (Cholesky factorization of MM) (no
check of orthonormality)

esType = 10 Uses LAPACK routine for simple eigenproblem on KK (MM is
not referenced in this case)

The code accesses only the upper triangular part of KK and MM.

Integer info on the status of the computation Return the integer info
on the status of the computation info = 0 >> Success

info = - 20 >> Failure in LAPACK routine ";

/*  Sanity Checking Methods  */

%feature("docstring")  Anasazi::SolverUtils::errorEquality "Teuchos::ScalarTraits< ScalarType >::magnitudeType
Anasazi::SolverUtils< ScalarType, MV, OP >::errorEquality(const MV &X,
const MV &MX, Teuchos::RCP< const OP > M=Teuchos::null)

Return the maximum coefficient of the matrix $M * X - MX$ scaled by
the maximum coefficient of MX.

When M is not specified, the identity is used. ";

/*  Internal Typedefs  */


// File: classAnasazi_1_1SortManager.xml
%feature("docstring") Anasazi::SortManager "

Anasazi's templated pure virtual class for managing the sorting of
approximate eigenvalues computed by the eigensolver. A concrete
implementation of this class is necessary.

Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist

C++ includes: AnasaziSortManager.hpp ";

%feature("docstring")  Anasazi::SortManager::SortManager "Anasazi::SortManager< MagnitudeType >::SortManager()

Default constructor. ";

%feature("docstring")  Anasazi::SortManager::SortManager "Anasazi::SortManager< MagnitudeType
>::SortManager(Teuchos::ParameterList &pl)

Constructor accepting a Teuchos::ParameterList. This is the default
mode for instantiating a SortManager. ";

%feature("docstring")  Anasazi::SortManager::~SortManager "virtual
Anasazi::SortManager< MagnitudeType >::~SortManager()

Destructor. ";

%feature("docstring")  Anasazi::SortManager::sort "virtual void
Anasazi::SortManager< MagnitudeType >::sort(std::vector< MagnitudeType
> &evals, Teuchos::RCP< std::vector< int > > perm=Teuchos::null, int
n=-1) const =0

Sort real eigenvalues, optionally returning the permutation vector.

Parameters:
-----------

evals:  [in/out] Vector of length at least n containing the
eigenvalues to be sorted.  On output, the first n eigenvalues will be
sorted. The rest will be unchanged.

perm:  [out] Vector of length at least n to store the permutation
index (optional).  If specified, on output the first n eigenvalues
will contain the permutation indices, in the range [0,n-1], such that
evals_out[i] = evals_in[perm[i]]

n:  [in] Number of values in evals to be sorted. If n == -1, all
values will be sorted. ";

%feature("docstring")  Anasazi::SortManager::sort "virtual void
Anasazi::SortManager< MagnitudeType >::sort(std::vector< MagnitudeType
> &r_evals, std::vector< MagnitudeType > &i_evals, Teuchos::RCP<
std::vector< int > > perm=Teuchos::null, int n=-1) const =0

Sort complex eigenvalues, optionally returning the permutation vector.

This routine takes two vectors, one for each part of a complex
eigenvalue. This is helpful for solving real, non-symmetric eigenvalue
problems.

Parameters:
-----------

r_evals:  [in/out] Vector of length at least n containing the real
part of the eigenvalues to be sorted.  On output, the first n
eigenvalues will be sorted. The rest will be unchanged.

i_evals:  [in/out] Vector of length at least n containing the
imaginary part of the eigenvalues to be sorted.  On output, the first
n eigenvalues will be sorted. The rest will be unchanged.

perm:  [out] Vector of length at least n to store the permutation
index (optional).  If specified, on output the first n eigenvalues
will contain the permutation indices, in the range [0,n-1], such that
r_evals_out[i] = r_evals_in[perm[i]] and similarly for i_evals.

n:  [in] Number of values in r_evals, i_evals to be sorted. If n ==
-1, all values will be sorted. ";


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

A pure virtual class for defining the status tests for the ::Anasazi
iterative solvers.

StatusTest is an interface that can be implemented to create
convergence tests for all Anasazi solvers. Almost any kind of test can
be expressed using this mechanism, including composite tests (see
StatusTestCombo).

C++ includes: AnasaziStatusTestDecl.hpp ";

/*  Constructors/destructors  */

%feature("docstring")  Anasazi::StatusTest::StatusTest "Anasazi::StatusTest< ScalarType, MV, OP >::StatusTest()

Constructor. ";

%feature("docstring")  Anasazi::StatusTest::~StatusTest "virtual
Anasazi::StatusTest< ScalarType, MV, OP >::~StatusTest()

Destructor. ";

/*  Status methods  */

%feature("docstring")  Anasazi::StatusTest::checkStatus "virtual
TestStatus Anasazi::StatusTest< ScalarType, MV, OP
>::checkStatus(Eigensolver< ScalarType, MV, OP > *solver)=0

Check status as defined by test.

TestStatus indicating whether the test passed or failed. ";

%feature("docstring")  Anasazi::StatusTest::getStatus "virtual
TestStatus Anasazi::StatusTest< ScalarType, MV, OP >::getStatus()
const =0

Return the result of the most recent checkStatus call, or undefined if
it has not been run. ";

%feature("docstring")  Anasazi::StatusTest::whichVecs "virtual
std::vector<int> Anasazi::StatusTest< ScalarType, MV, OP
>::whichVecs() const =0

Get the indices for the vectors that passed the test. ";

%feature("docstring")  Anasazi::StatusTest::howMany "virtual int
Anasazi::StatusTest< ScalarType, MV, OP >::howMany() const =0

Get the number of vectors that passed the test. ";

/*  Reset methods  */

%feature("docstring")  Anasazi::StatusTest::reset "virtual void
Anasazi::StatusTest< ScalarType, MV, OP >::reset()=0

Informs the status test that it should reset its internal
configuration to the uninitialized state.

This is necessary for the case when the status test is being reused by
another solver or for another eigenvalue problem. The status test may
have information that pertains to a particular problem or solver
state. The internal information will be reset back to the
uninitialized state. The user specified information that the
convergence test uses will remain. ";

%feature("docstring")  Anasazi::StatusTest::clearStatus "virtual void
Anasazi::StatusTest< ScalarType, MV, OP >::clearStatus()=0

Clears the results of the last status test.

This should be distinguished from the reset() method, as it only
clears the cached result from the last status test, so that a call to
getStatus() will return ::Undefined. This is necessary for the SEQOR
and SEQAND tests in the StatusTestCombo class, which may short circuit
and not evaluate all of the StatusTests contained in them. ";

/*  Print methods  */

%feature("docstring")  Anasazi::StatusTest::print "virtual
std::ostream& Anasazi::StatusTest< ScalarType, MV, OP
>::print(std::ostream &os, int indent=0) const =0

Output formatted description of stopping test to output stream. ";


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
evaluates ::Failed. StatusTestCombo::SEQOR run the tests in the order
they were given to the StatusTestCombo class and stops after the first
test that evaluates ::Passed.

C++ includes: AnasaziStatusTestCombo.hpp ";

/*  Constructors/destructors  */

%feature("docstring")  Anasazi::StatusTestCombo::StatusTestCombo "Anasazi::StatusTestCombo< ScalarType, MV, OP >::StatusTestCombo()

Default constructor has no tests and initializes to
StatusTestCombo::ComboType StatusTestCombo::OR.

Constructor ";

%feature("docstring")  Anasazi::StatusTestCombo::StatusTestCombo "Anasazi::StatusTestCombo< ScalarType, MV, OP
>::StatusTestCombo(ComboType type, Teuchos::Array< Teuchos::RCP<
StatusTest< ScalarType, MV, OP > > > tests)

Constructor specifying the StatusTestCombo::ComboType and the tests.

Constructor ";

%feature("docstring")  Anasazi::StatusTestCombo::~StatusTestCombo "virtual Anasazi::StatusTestCombo< ScalarType, MV, OP
>::~StatusTestCombo()

Destructor. ";

/*  Status methods  */

%feature("docstring")  Anasazi::StatusTestCombo::checkStatus "TestStatus Anasazi::StatusTestCombo< ScalarType, MV, OP
>::checkStatus(Eigensolver< ScalarType, MV, OP > *solver)

Check status as defined by test.

TestStatus indicating whether the test passed or failed. ";

%feature("docstring")  Anasazi::StatusTestCombo::getStatus "TestStatus Anasazi::StatusTestCombo< ScalarType, MV, OP >::getStatus()
const

Return the result of the most recent checkStatus call. ";

%feature("docstring")  Anasazi::StatusTestCombo::whichVecs "std::vector<int> Anasazi::StatusTestCombo< ScalarType, MV, OP
>::whichVecs() const

Get the indices for the vectors that passed the test.

This returns some combination of the passing vectors from the tests
comprising the StatusTestCombo. The nature of the combination depends
on the StatusTestCombo::ComboType:  StatusTestCombo::SEQOR,
StatusTestCombo::OR - whichVecs() returns the union of whichVecs()
from all evaluated constituent tests

StatusTestCombo::SEQAND, StatusTestCombo::AND - whichVecs() returns
the intersection of whichVecs() from all evaluated constituent tests
";

%feature("docstring")  Anasazi::StatusTestCombo::howMany "int
Anasazi::StatusTestCombo< ScalarType, MV, OP >::howMany() const

Get the number of vectors that passed the test. ";

/*  Accessor methods  */

%feature("docstring")  Anasazi::StatusTestCombo::setComboType "void
Anasazi::StatusTestCombo< ScalarType, MV, OP >::setComboType(ComboType
type)

Set the maximum number of iterations. This also resets the test status
to ::Undefined. ";

%feature("docstring")  Anasazi::StatusTestCombo::getComboType "ComboType Anasazi::StatusTestCombo< ScalarType, MV, OP
>::getComboType() const

Get the maximum number of iterations. ";

%feature("docstring")  Anasazi::StatusTestCombo::setTests "void
Anasazi::StatusTestCombo< ScalarType, MV, OP
>::setTests(Teuchos::Array< Teuchos::RCP< StatusTest< ScalarType, MV,
OP > > > tests)

Set the tests This also resets the test status to ::Undefined. ";

%feature("docstring")  Anasazi::StatusTestCombo::getTests "Teuchos::Array<Teuchos::RCP<StatusTest<ScalarType,MV,OP> > >
Anasazi::StatusTestCombo< ScalarType, MV, OP >::getTests() const

Get the tests. ";

%feature("docstring")  Anasazi::StatusTestCombo::addTest "void
Anasazi::StatusTestCombo< ScalarType, MV, OP >::addTest(Teuchos::RCP<
StatusTest< ScalarType, MV, OP > > test)

Add a test to the combination.

This also resets the test status to ::Undefined. ";

%feature("docstring")  Anasazi::StatusTestCombo::removeTest "void
Anasazi::StatusTestCombo< ScalarType, MV, OP >::removeTest(const
Teuchos::RCP< StatusTest< ScalarType, MV, OP > > &test)

Removes a test from the combination, if it exists in the tester.

This also resets the test status to ::Undefined, if a test was
removed. ";

/*  Reset methods  */

%feature("docstring")  Anasazi::StatusTestCombo::reset "void
Anasazi::StatusTestCombo< ScalarType, MV, OP >::reset()

Informs the status test that it should reset its internal
configuration to the uninitialized state.

The StatusTestCombo class has no internal state, but children classes
might, so this method will call reset() on all child status tests. It
also resets the test status to ::Undefined. ";

%feature("docstring")  Anasazi::StatusTestCombo::clearStatus "void
Anasazi::StatusTestCombo< ScalarType, MV, OP >::clearStatus()

Clears the results of the last status test.

This should be distinguished from the reset() method, as it only
clears the cached result from the last status test, so that a call to
getStatus() will return ::Undefined. This is necessary for the
StatusTestCombo::SEQOR and StatusTestCombo::SEQAND tests in the
StatusTestCombo class, which may short circuit and not evaluate all of
the StatusTests contained in them. ";

/*  Print methods  */

%feature("docstring")  Anasazi::StatusTestCombo::print "std::ostream
& Anasazi::StatusTestCombo< ScalarType, MV, OP >::print(std::ostream
&os, int indent=0) const

Output formatted description of stopping test to output stream. ";


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

/*  Constructors/destructors  */

%feature("docstring")  Anasazi::StatusTestMaxIters::StatusTestMaxIters
"Anasazi::StatusTestMaxIters< ScalarType, MV, OP
>::StatusTestMaxIters(int maxIter, bool negate=false)

Constructor. ";

%feature("docstring")
Anasazi::StatusTestMaxIters::~StatusTestMaxIters "virtual
Anasazi::StatusTestMaxIters< ScalarType, MV, OP
>::~StatusTestMaxIters()

Destructor. ";

/*  Status methods  */

%feature("docstring")  Anasazi::StatusTestMaxIters::checkStatus "TestStatus Anasazi::StatusTestMaxIters< ScalarType, MV, OP
>::checkStatus(Eigensolver< ScalarType, MV, OP > *solver)

Check status as defined by test.

TestStatus indicating whether the test passed or failed. ";

%feature("docstring")  Anasazi::StatusTestMaxIters::getStatus "TestStatus Anasazi::StatusTestMaxIters< ScalarType, MV, OP
>::getStatus() const

Return the result of the most recent checkStatus call. ";

%feature("docstring")  Anasazi::StatusTestMaxIters::whichVecs "std::vector<int> Anasazi::StatusTestMaxIters< ScalarType, MV, OP
>::whichVecs() const

Get the indices for the vectors that passed the test. ";

%feature("docstring")  Anasazi::StatusTestMaxIters::howMany "int
Anasazi::StatusTestMaxIters< ScalarType, MV, OP >::howMany() const

Get the number of vectors that passed the test. ";

/*  Accessor methods  */

%feature("docstring")  Anasazi::StatusTestMaxIters::setMaxIters "void
Anasazi::StatusTestMaxIters< ScalarType, MV, OP >::setMaxIters(int
maxIters)

Set the maximum number of iterations.

This also resets the test status to ::Undefined. ";

%feature("docstring")  Anasazi::StatusTestMaxIters::getMaxIters "int
Anasazi::StatusTestMaxIters< ScalarType, MV, OP >::getMaxIters()

Get the maximum number of iterations. ";

%feature("docstring")  Anasazi::StatusTestMaxIters::setNegate "void
Anasazi::StatusTestMaxIters< ScalarType, MV, OP >::setNegate(bool
negate)

Set the negation policy for the status test.

This also reset the test status to ::Undefined. ";

%feature("docstring")  Anasazi::StatusTestMaxIters::getNegate "bool
Anasazi::StatusTestMaxIters< ScalarType, MV, OP >::getNegate() const

Get the negation policy for the status test. ";

/*  Reset methods  */

%feature("docstring")  Anasazi::StatusTestMaxIters::reset "void
Anasazi::StatusTestMaxIters< ScalarType, MV, OP >::reset()

Informs the status test that it should reset its internal
configuration to the uninitialized state.

The StatusTestMaxIters class has no internal state, so this call is
equivalent to calling clearStatus(). eigenvalue problem. The status
test may have information that pertains to a particular problem or
solver state. The internal information will be reset back to the
uninitialized state. The user specified information that the
convergence test uses will remain. ";

%feature("docstring")  Anasazi::StatusTestMaxIters::clearStatus "void
Anasazi::StatusTestMaxIters< ScalarType, MV, OP >::clearStatus()

Clears the results of the last status test.

This should be distinguished from the reset() method, as it only
clears the cached result from the last status test, so that a call to
getStatus() will return ::Undefined. This is necessary for the SEQOR
and SEQAND tests in the StatusTestCombo class, which may short circuit
and not evaluate all of the StatusTests contained in them. ";

/*  Print methods  */

%feature("docstring")  Anasazi::StatusTestMaxIters::print "std::ostream& Anasazi::StatusTestMaxIters< ScalarType, MV, OP
>::print(std::ostream &os, int indent=0) const

Output formatted description of stopping test to output stream. ";


// File: classAnasazi_1_1StatusTestOutput.xml
%feature("docstring") Anasazi::StatusTestOutput "

A special StatusTest for printing other status tests.

StatusTestOutput is a wrapper around another StatusTest that calls
StatusTest::print() on the underlying object on calls to
StatusTestOutput::checkStatus(). The frequency and occasion of the
printing can be dictated according to some parameters passed to
StatusTestOutput::StatusTestOutput().

C++ includes: AnasaziStatusTestOutput.hpp ";

/*  Constructors/destructors  */

%feature("docstring")  Anasazi::StatusTestOutput::StatusTestOutput "Anasazi::StatusTestOutput< ScalarType, MV, OP
>::StatusTestOutput(const Teuchos::RCP< OutputManager< ScalarType > >
&printer, Teuchos::RCP< StatusTest< ScalarType, MV, OP > > test, int
mod=1, int printStates=Passed)

Constructor.

The StatusTestOutput requires an OutputManager for printing the
underlying StatusTest on calls to checkStatus(), as well as an
underlying StatusTest.

StatusTestOutput can be initialized with a null pointer for argument
test. However, calling checkStatus() with a null child pointer will
result in a StatusTestError exception being thrown. See checkStatus()
for more information.

The last two parameters, described below, in addition to the verbosity
level of the OutputManager, control when printing is called. When both
the mod criterion and the printStates criterion are satisfied, the
status test will be printed to the OutputManager with ::MsgType of
::StatusTestDetails.

Parameters:
-----------

mod:  A positive number describes how often the output should be
printed. On every call to checkStatus(), an internal counter is
incremented. Printing may only occur when this counter is congruent to
zero modulo mod. Default: 1 (attempt to print on every call to
checkStatus())

printStates:  A combination of ::TestStatus values for which the
output may be printed. Default: ::Passed (attempt to print whenever
checkStatus() will return ::Passed) ";

%feature("docstring")  Anasazi::StatusTestOutput::~StatusTestOutput "virtual Anasazi::StatusTestOutput< ScalarType, MV, OP
>::~StatusTestOutput()

Destructor. ";

/*  Status methods  */

%feature("docstring")  Anasazi::StatusTestOutput::checkStatus "TestStatus Anasazi::StatusTestOutput< ScalarType, MV, OP
>::checkStatus(Eigensolver< ScalarType, MV, OP > *solver)

Check and return status of underlying StatusTest.

This method calls checkStatus() on the StatusTest object passed in the
constructor. If appropriate, the method will follow this call with a
call to print() on the underlying object, using the OutputManager
passed via the constructor with verbosity level ::StatusTestDetails.

The internal counter will be incremented during this call, but only
after performing the tests to decide whether or not to print the
underlying StatusTest. This way, the very first call to checkStatus()
following initialization or reset() will enable the underlying
StatusTest to be printed, regardless of the mod parameter, as the
current number of calls will be zero.

If the specified Teuchos::RCP for the child class is Teuchos::null,
then calling checkStatus() will result in a StatusTestError exception
being thrown.

::TestStatus indicating whether the underlying test passed or failed.
";

%feature("docstring")  Anasazi::StatusTestOutput::getStatus "TestStatus Anasazi::StatusTestOutput< ScalarType, MV, OP
>::getStatus() const

Return the result of the most recent checkStatus call, or undefined if
it has not been run. ";

%feature("docstring")  Anasazi::StatusTestOutput::whichVecs "std::vector<int> Anasazi::StatusTestOutput< ScalarType, MV, OP
>::whichVecs() const

Get the indices for the vectors that passed the test. ";

%feature("docstring")  Anasazi::StatusTestOutput::howMany "int
Anasazi::StatusTestOutput< ScalarType, MV, OP >::howMany() const

Get the number of vectors that passed the test. ";

/*  Accessor methods  */

%feature("docstring")  Anasazi::StatusTestOutput::setChild "void
Anasazi::StatusTestOutput< ScalarType, MV, OP
>::setChild(Teuchos::RCP< StatusTest< ScalarType, MV, OP > > test)

Set child test.

This also resets the test status to ::Undefined. ";

%feature("docstring")  Anasazi::StatusTestOutput::getChild "Teuchos::RCP<StatusTest<ScalarType,MV,OP> > Anasazi::StatusTestOutput<
ScalarType, MV, OP >::getChild() const

Get child test. ";

/*  Reset methods  */

%feature("docstring")  Anasazi::StatusTestOutput::reset "void
Anasazi::StatusTestOutput< ScalarType, MV, OP >::reset()

Informs the status test that it should reset its internal
configuration to the uninitialized state.

This resets the cached state to an ::Undefined state and calls reset()
on the underlying test. It also resets the counter for the number of
calls to checkStatus(). ";

%feature("docstring")  Anasazi::StatusTestOutput::clearStatus "void
Anasazi::StatusTestOutput< ScalarType, MV, OP >::clearStatus()

Clears the results of the last status test. This resets the cached
state to an ::Undefined state and calls clearStatus() on the
underlying test. ";

/*  Print methods  */

%feature("docstring")  Anasazi::StatusTestOutput::print "std::ostream& Anasazi::StatusTestOutput< ScalarType, MV, OP
>::print(std::ostream &os, int indent=0) const

Output formatted description of stopping test to output stream. ";


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
::Passed.

C++ includes: AnasaziStatusTestResNorm.hpp ";

/*  Constructors/destructors  */

%feature("docstring")  Anasazi::StatusTestResNorm::StatusTestResNorm "Anasazi::StatusTestResNorm< ScalarType, MV, OP
>::StatusTestResNorm(typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType tol, int quorum=-1, ResType whichNorm=RES_ORTH, bool
scaled=true, bool throwExceptionOnNaN=true)

Constructor. ";

%feature("docstring")  Anasazi::StatusTestResNorm::~StatusTestResNorm
"virtual Anasazi::StatusTestResNorm< ScalarType, MV, OP
>::~StatusTestResNorm()

Destructor. ";

/*  Status methods  */

%feature("docstring")  Anasazi::StatusTestResNorm::checkStatus "TestStatus Anasazi::StatusTestResNorm< ScalarType, MV, OP
>::checkStatus(Eigensolver< ScalarType, MV, OP > *solver)

Check status as defined by test.

TestStatus indicating whether the test passed or failed. ";

%feature("docstring")  Anasazi::StatusTestResNorm::getStatus "TestStatus Anasazi::StatusTestResNorm< ScalarType, MV, OP
>::getStatus() const

Return the result of the most recent checkStatus call, or undefined if
it has not been run. ";

%feature("docstring")  Anasazi::StatusTestResNorm::whichVecs "std::vector<int> Anasazi::StatusTestResNorm< ScalarType, MV, OP
>::whichVecs() const

Get the indices for the vectors that passed the test. ";

%feature("docstring")  Anasazi::StatusTestResNorm::howMany "int
Anasazi::StatusTestResNorm< ScalarType, MV, OP >::howMany() const

Get the number of vectors that passed the test. ";

/*  Accessor methods  */

%feature("docstring")  Anasazi::StatusTestResNorm::setQuorum "void
Anasazi::StatusTestResNorm< ScalarType, MV, OP >::setQuorum(int
quorum)

Set quorum.

Setting quorum to -1 signifies that all residuals from the solver must
meet the tolerance. This also resets the test status to ::Undefined.
";

%feature("docstring")  Anasazi::StatusTestResNorm::getQuorum "int
Anasazi::StatusTestResNorm< ScalarType, MV, OP >::getQuorum() const

Get quorum. ";

%feature("docstring")  Anasazi::StatusTestResNorm::setTolerance "void
Anasazi::StatusTestResNorm< ScalarType, MV, OP
>::setTolerance(typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType tol)

Set tolerance. This also resets the test status to ::Undefined. ";

%feature("docstring")  Anasazi::StatusTestResNorm::getTolerance "Teuchos::ScalarTraits<ScalarType>::magnitudeType
Anasazi::StatusTestResNorm< ScalarType, MV, OP >::getTolerance()

Get tolerance. ";

%feature("docstring")  Anasazi::StatusTestResNorm::setWhichNorm "void
Anasazi::StatusTestResNorm< ScalarType, MV, OP >::setWhichNorm(ResType
whichNorm)

Set the residual norm to be used by the status test.

This also resets the test status to ::Undefined. ";

%feature("docstring")  Anasazi::StatusTestResNorm::getWhichNorm "ResType Anasazi::StatusTestResNorm< ScalarType, MV, OP
>::getWhichNorm()

Return the residual norm used by the status test. ";

%feature("docstring")  Anasazi::StatusTestResNorm::setScale "void
Anasazi::StatusTestResNorm< ScalarType, MV, OP >::setScale(bool
relscale)

Instruct test to scale norms by eigenvalue estimates (relative scale).
This also resets the test status to ::Undefined. ";

%feature("docstring")  Anasazi::StatusTestResNorm::getScale "bool
Anasazi::StatusTestResNorm< ScalarType, MV, OP >::getScale()

Returns true if the test scales the norms by the eigenvalue estimates
(relative scale). ";

/*  Reset methods  */

%feature("docstring")  Anasazi::StatusTestResNorm::reset "void
Anasazi::StatusTestResNorm< ScalarType, MV, OP >::reset()

Informs the status test that it should reset its internal
configuration to the uninitialized state.

This is necessary for the case when the status test is being reused by
another solver or for another eigenvalue problem. The status test may
have information that pertains to a particular problem or solver
state. The internal information will be reset back to the
uninitialized state. The user specified information that the
convergence test uses will remain. ";

%feature("docstring")  Anasazi::StatusTestResNorm::clearStatus "void
Anasazi::StatusTestResNorm< ScalarType, MV, OP >::clearStatus()

Clears the results of the last status test.

This should be distinguished from the reset() method, as it only
clears the cached result from the last status test, so that a call to
getStatus() will return ::Undefined. This is necessary for the SEQOR
and SEQAND tests in the StatusTestCombo class, which may short circuit
and not evaluate all of the StatusTests contained in them. ";

/*  Print methods  */

%feature("docstring")  Anasazi::StatusTestResNorm::print "std::ostream & Anasazi::StatusTestResNorm< ScalarType, MV, OP
>::print(std::ostream &os, int indent=0) const

Output formatted description of stopping test to output stream. ";


// File: classAnasazi_1_1StatusTestWithOrdering.xml
%feature("docstring") Anasazi::StatusTestWithOrdering "

A status test for testing the norm of the eigenvectors residuals along
with a set of auxiliary eigenvalues.

The test evaluates to ::Passed when then the most significant of the
eigenvalues all have a residual below a certain threshhold. The
purpose of the test is to not only test convergence for some number of
eigenvalues, but to test convergence for the correct ones.

In addition to specifying the tolerance, the user may specify: the
norm to be used: 2-norm or OrthoManager::norm() or getRitzRes2Norms()

the scale: absolute or relative to magnitude of Ritz value

the quorum: the number of vectors required for the test to evaluate as
::Passed.

Finally, the user must specify the Anasazi::SortManager used for
deciding significance.

C++ includes: AnasaziStatusTestWithOrdering.hpp ";

/*  Constructors/destructors  */

%feature("docstring")
Anasazi::StatusTestWithOrdering::StatusTestWithOrdering "Anasazi::StatusTestWithOrdering< ScalarType, MV, OP
>::StatusTestWithOrdering(Teuchos::RCP< StatusTest< ScalarType, MV, OP
> > test, Teuchos::RCP< SortManager< typename Teuchos::ScalarTraits<
ScalarType >::magnitudeType > > sorter, int quorum=-1)

Constructor. ";

%feature("docstring")
Anasazi::StatusTestWithOrdering::~StatusTestWithOrdering "virtual
Anasazi::StatusTestWithOrdering< ScalarType, MV, OP
>::~StatusTestWithOrdering()

Destructor. ";

/*  Status methods  */

%feature("docstring")  Anasazi::StatusTestWithOrdering::checkStatus "TestStatus Anasazi::StatusTestWithOrdering< ScalarType, MV, OP
>::checkStatus(Eigensolver< ScalarType, MV, OP > *solver)

Check status as defined by test.

TestStatus indicating whether the test passed or failed. ";

%feature("docstring")  Anasazi::StatusTestWithOrdering::getStatus "TestStatus Anasazi::StatusTestWithOrdering< ScalarType, MV, OP
>::getStatus() const

Return the result of the most recent checkStatus call, or undefined if
it has not been run. ";

%feature("docstring")  Anasazi::StatusTestWithOrdering::whichVecs "std::vector<int> Anasazi::StatusTestWithOrdering< ScalarType, MV, OP
>::whichVecs() const

Get the indices for the vectors that passed the test.

Non-negative indices correspond to passing vectors from the
constituent status test. Negative entries correspond to auxilliary
values, where the first auxilliary value is indexed by -NumAuxVals,
the second by -NumAuxVals+1, and so forth. ";

%feature("docstring")  Anasazi::StatusTestWithOrdering::howMany "int
Anasazi::StatusTestWithOrdering< ScalarType, MV, OP >::howMany() const

Get the number of vectors that passed the test. ";

/*  Accessor methods  */

%feature("docstring")  Anasazi::StatusTestWithOrdering::setQuorum "void Anasazi::StatusTestWithOrdering< ScalarType, MV, OP
>::setQuorum(int quorum)

Set quorum.

Setting quorum to -1 signifies that all residuals from the solver must
meet the tolerance. This also resets the test status to ::Undefined.
";

%feature("docstring")  Anasazi::StatusTestWithOrdering::getQuorum "int Anasazi::StatusTestWithOrdering< ScalarType, MV, OP >::getQuorum()
const

Get quorum. ";

/*  Reset methods  */

%feature("docstring")  Anasazi::StatusTestWithOrdering::reset "void
Anasazi::StatusTestWithOrdering< ScalarType, MV, OP >::reset()

Informs the status test that it should reset its internal
configuration to the uninitialized state.

This is necessary for the case when the status test is being reused by
another solver or for another eigenvalue problem. The status test may
have information that pertains to a particular problem or solver
state. The internal information will be reset back to the
uninitialized state. The user specified information that the
convergence test uses will remain. ";

%feature("docstring")  Anasazi::StatusTestWithOrdering::clearStatus "void Anasazi::StatusTestWithOrdering< ScalarType, MV, OP
>::clearStatus()

Clears the results of the last status test.

This should be distinguished from the reset() method, as it only
clears the cached result from the last status test, so that a call to
getStatus() will return ::Undefined. This is necessary for the SEQOR
and SEQAND tests in the StatusTestCombo class, which may short circuit
and not evaluate all of the StatusTests contained in them. ";

%feature("docstring")  Anasazi::StatusTestWithOrdering::setAuxVals "void Anasazi::StatusTestWithOrdering< ScalarType, MV, OP
>::setAuxVals(const std::vector< typename Teuchos::ScalarTraits<
ScalarType >::magnitudeType > &vals)

Set the auxiliary eigenvalues.

This routine sets only the real part of the auxiliary eigenvalues; the
imaginary part is set to zero. This routine also resets the state to
::Undefined. ";

%feature("docstring")  Anasazi::StatusTestWithOrdering::setAuxVals "void Anasazi::StatusTestWithOrdering< ScalarType, MV, OP
>::setAuxVals(const std::vector< typename Teuchos::ScalarTraits<
ScalarType >::magnitudeType > &rvals, const std::vector< typename
Teuchos::ScalarTraits< ScalarType >::magnitudeType > &ivals)

Set the auxiliary eigenvalues.

This routine sets both the real and imaginary parts of the auxiliary
eigenvalues. This routine also resets the state to ::Undefined. ";

%feature("docstring")  Anasazi::StatusTestWithOrdering::getAuxVals "void Anasazi::StatusTestWithOrdering< ScalarType, MV, OP
>::getAuxVals(std::vector< typename Teuchos::ScalarTraits< ScalarType
>::magnitudeType > &rvals, std::vector< typename
Teuchos::ScalarTraits< ScalarType >::magnitudeType > &ivals) const

Get the auxiliary eigenvalues.

This routine gets the real and imaginary parts of the auxiliary
eigenvalues. ";

/*  Print methods  */

%feature("docstring")  Anasazi::StatusTestWithOrdering::print "std::ostream & Anasazi::StatusTestWithOrdering< ScalarType, MV, OP
>::print(std::ostream &os, int indent=0) const

Output formatted description of stopping test to output stream. ";


// File: classAnasazi_1_1details_1_1StubTsqrAdapter.xml
%feature("docstring") Anasazi::details::StubTsqrAdapter "

\"Stub\" TSQR adaptor for unsupported multivector types.

TSQR (Tall Skinny QR factorization) is an orthogonalization kernel
that is as accurate as Householder QR, yet requires only $2 \\\\log P$
messages between $P$ MPI processes, independently of the number of
columns in the multivector.

TSQR works independently of the particular multivector implementation,
and interfaces to the latter via an adapter class. Each multivector
type MV needs its own adapter class. The specialization of
MultiVecTraits for MV refers to its corresponding adapter class as its
tsqr_adaptor_type [sic; sorry about the lack of standard spelling of
\"adapter\"] typedef. For examples, please refer to the
Epetra_MultiVector and Tpetra::MultiVector specializations of
Anasazi::MultiVecTraits.

Nevertheless, there may be multivector types for which a TSQR adapter
has not yet been written. This \"stub\" adapter implements the
interface that TSQR adapters must implement, but all of its methods
throw std::logic_error to indicate that this is a stub. Thus, it
allows Anasazi classes like TsqrOrthoManagerImpl to compile
successfully for unsupported MV types. This in turn allows
OrthoManagerFactory to be templated on the MV type.

C++ includes: AnasaziStubTsqrAdapter.hpp ";

%feature("docstring")
Anasazi::details::StubTsqrAdapter::StubTsqrAdapter "Anasazi::details::StubTsqrAdapter< MultiVectorType
>::StubTsqrAdapter(const Teuchos::RCP< Teuchos::ParameterList >
&plist)

Constructor (that accepts a parameter list).

Parameters:
-----------

plist:  [in] List of parameters for configuring TSQR. The specific
parameter keys that are read depend on the TSQR implementation. For
details, call  getValidParameters() and examine the documentation
embedded therein. ";

%feature("docstring")
Anasazi::details::StubTsqrAdapter::StubTsqrAdapter "Anasazi::details::StubTsqrAdapter< MultiVectorType
>::StubTsqrAdapter()

Default constructor (stub; throws std::logic_error). ";

%feature("docstring")
Anasazi::details::StubTsqrAdapter::StubTsqrAdapter "Anasazi::details::StubTsqrAdapter< MultiVectorType
>::StubTsqrAdapter(const StubTsqrAdapter &rhs)

Copy constructor (throws std::logic_error). ";

%feature("docstring")
Anasazi::details::StubTsqrAdapter::getValidParameters "Teuchos::RCP<const Teuchos::ParameterList>
Anasazi::details::StubTsqrAdapter< MultiVectorType
>::getValidParameters() const

Get list of valid default parameters (stub; throws std::logic_error).
";

%feature("docstring")
Anasazi::details::StubTsqrAdapter::setParameterList "void
Anasazi::details::StubTsqrAdapter< MultiVectorType
>::setParameterList(const Teuchos::RCP< Teuchos::ParameterList >
&plist)

Set parameters (stub; throws std::logic_error). ";

%feature("docstring")
Anasazi::details::StubTsqrAdapter::factorExplicit "void
Anasazi::details::StubTsqrAdapter< MultiVectorType
>::factorExplicit(MV &A, MV &Q, dense_matrix_type &R, const bool
forceNonnegativeDiagonal=false)

Compute QR factorization [Q,R] = qr(A,0) (stub; throws
std::logic_error). ";

%feature("docstring")  Anasazi::details::StubTsqrAdapter::revealRank "int Anasazi::details::StubTsqrAdapter< MultiVectorType
>::revealRank(MV &Q, dense_matrix_type &R, const magnitude_type &tol)

Rank-revealing decomposition (stub; does nothing). ";


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

/*  Constructor/Destructor  */

%feature("docstring")  Anasazi::SVQBOrthoManager::SVQBOrthoManager "Anasazi::SVQBOrthoManager< ScalarType, MV, OP
>::SVQBOrthoManager(Teuchos::RCP< const OP > Op=Teuchos::null, bool
debug=false)

Constructor specifying re-orthogonalization tolerance. ";

%feature("docstring")  Anasazi::SVQBOrthoManager::~SVQBOrthoManager "Anasazi::SVQBOrthoManager< ScalarType, MV, OP >::~SVQBOrthoManager()

Destructor. ";

/*  Methods implementing Anasazi::MatOrthoManager  */

%feature("docstring")  Anasazi::SVQBOrthoManager::projectMat "void
Anasazi::SVQBOrthoManager< ScalarType, MV, OP >::projectMat(MV &X,
Teuchos::Array< Teuchos::RCP< const MV > > Q, Teuchos::Array<
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > > >
C=Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix< int,
ScalarType > >(Teuchos::null)), Teuchos::RCP< MV > MX=Teuchos::null,
Teuchos::Array< Teuchos::RCP< const MV > >
MQ=Teuchos::tuple(Teuchos::RCP< const MV >(Teuchos::null))) const

Given a list of mutually orthogonal and internally orthonormal bases
Q, this method projects a multivector X onto the space orthogonal to
the individual Q[i], optionally returning the coefficients of X for
the individual Q[i]. All of this is done with respect to the inner
product innerProd().

After calling this routine, X will be orthogonal to each of the Q[i].

Parameters:
-----------

X:  [in/out] The multivector to be modified.  On output, the columns
of X will be orthogonal to each Q[i], satisfying \\\\[ X_{out} =
X_{in} - \\\\sum_i Q[i] \\\\langle Q[i], X_{in} \\\\rangle \\\\]

MX:  [in/out] The image of X under the inner product operator Op. If $
MX != 0$: On input, this is expected to be consistent with Op  X. On
output, this is updated consistent with updates to X. If $ MX == 0$ or
$ Op == 0$: MX is not referenced.

C:  [out] The coefficients of X in the bases Q[i]. If C[i] is a non-
null pointer and C[i] matches the dimensions of X and Q[i], then the
coefficients computed during the orthogonalization routine will be
stored in the matrix C[i], similar to calling If C[i] points to a
Teuchos::SerialDenseMatrix with size inconsistent with X and  Q[i],
then a std::invalid_argument exception will be thrown. Otherwise, if
C.size() < i or C[i] is a null pointer, the caller will not have
access to the computed coefficients.

Q:  [in] A list of multivector bases specifying the subspaces to be
orthogonalized against, satisfying \\\\[ \\\\langle Q[i], Q[j]
\\\\rangle = I \\\\quad\\\\textrm{if}\\\\quad i=j \\\\] and \\\\[
\\\\langle Q[i], Q[j] \\\\rangle = 0 \\\\quad\\\\textrm{if}\\\\quad i
\\\\neq j\\\\ . \\\\] ";

%feature("docstring")  Anasazi::SVQBOrthoManager::normalizeMat "int
Anasazi::SVQBOrthoManager< ScalarType, MV, OP >::normalizeMat(MV &X,
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > >
B=Teuchos::null, Teuchos::RCP< MV > MX=Teuchos::null) const

This method takes a multivector X and attempts to compute an
orthonormal basis for $colspan(X)$, with respect to innerProd().

This method does not compute an upper triangular coefficient matrix B.

This routine returns an integer rank stating the rank of the computed
basis. If X does not have full rank and the normalize() routine does
not attempt to augment the subspace, then rank may be smaller than the
number of columns in X. In this case, only the first rank columns of
output X and first rank rows of B will be valid.

The method attempts to find a basis with dimension equal to the number
of columns in X. It does this by augmenting linearly dependent vectors
in X with random directions. A finite number of these attempts will be
made; therefore, it is possible that the dimension of the computed
basis is less than the number of vectors in X.

Parameters:
-----------

X:  [in/out] The multivector to be modified.  On output, the first
rank columns of X satisfy \\\\[ \\\\langle X[i], X[j] \\\\rangle =
\\\\delta_{ij}\\\\ . \\\\] Also, \\\\[ X_{in}(1:m,1:n) =
X_{out}(1:m,1:rank) B(1:rank,1:n) \\\\] where m is the number of rows
in X and n is the number of columns in X.

MX:  [in/out] The image of X under the inner product operator Op. If $
MX != 0$: On input, this is expected to be consistent with Op  X. On
output, this is updated consistent with updates to X. If $ MX == 0$ or
$ Op == 0$: MX is not referenced.

B:  [out] The coefficients of the original X with respect to the
computed basis. If B is a non-null pointer and B matches the
dimensions of B, then the coefficients computed during the
orthogonalization routine will be stored in B, similar to calling If B
points to a Teuchos::SerialDenseMatrix with size inconsistent with X,
then a std::invalid_argument exception will be thrown. Otherwise, if B
is null, the caller will not have access to the computed coefficients.
This matrix is not necessarily triangular (as in a QR factorization);
see the documentation of specific orthogonalization managers.  In
general, B has no non-zero structure.

Rank of the basis computed by this method, less than or equal to the
number of columns in X. This specifies how many columns in the
returned X and rows in the returned B are valid. ";

%feature("docstring")
Anasazi::SVQBOrthoManager::projectAndNormalizeMat "int
Anasazi::SVQBOrthoManager< ScalarType, MV, OP
>::projectAndNormalizeMat(MV &X, Teuchos::Array< Teuchos::RCP< const
MV > > Q, Teuchos::Array< Teuchos::RCP< Teuchos::SerialDenseMatrix<
int, ScalarType > > > C=Teuchos::tuple(Teuchos::RCP<
Teuchos::SerialDenseMatrix< int, ScalarType > >(Teuchos::null)),
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > >
B=Teuchos::null, Teuchos::RCP< MV > MX=Teuchos::null, Teuchos::Array<
Teuchos::RCP< const MV > > MQ=Teuchos::tuple(Teuchos::RCP< const MV
>(Teuchos::null))) const

Given a set of bases Q[i] and a multivector X, this method computes an
orthonormal basis for $colspan(X) - \\\\sum_i colspan(Q[i])$.

This routine returns an integer rank stating the rank of the computed
basis. If the subspace $colspan(X) - \\\\sum_i colspan(Q[i])$ does not
have dimension as large as the number of columns of X and the
orthogonalization manager doe not attempt to augment the subspace,
then rank may be smaller than the number of columns of X. In this
case, only the first rank columns of output X and first rank rows of B
will be valid.

The method attempts to find a basis with dimension the same as the
number of columns in X. It does this by augmenting linearly dependent
vectors with random directions. A finite number of these attempts will
be made; therefore, it is possible that the dimension of the computed
basis is less than the number of vectors in X.

Parameters:
-----------

X:  [in/out] The multivector to be modified.  On output, the first
rank columns of X satisfy \\\\[ \\\\langle X[i], X[j] \\\\rangle =
\\\\delta_{ij} \\\\quad \\\\textrm{and} \\\\quad \\\\langle X, Q[i]
\\\\rangle = 0\\\\ . \\\\] Also, \\\\[ X_{in}(1:m,1:n) =
X_{out}(1:m,1:rank) B(1:rank,1:n) + \\\\sum_i Q[i] C[i] \\\\] where m
is the number of rows in X and n is the number of columns in X.

MX:  [in/out] The image of X under the inner product operator Op. If $
MX != 0$: On input, this is expected to be consistent with Op  X. On
output, this is updated consistent with updates to X. If $ MX == 0$ or
$ Op == 0$: MX is not referenced.

C:  [out] The coefficients of X in the Q[i]. If C[i] is a non-null
pointer and C[i] matches the dimensions of X and Q[i], then the
coefficients computed during the orthogonalization routine will be
stored in the matrix C[i], similar to calling If C[i] points to a
Teuchos::SerialDenseMatrix with size inconsistent with X and  Q[i],
then a std::invalid_argument exception will be thrown. Otherwise, if
C.size() < i or C[i] is a null pointer, the caller will not have
access to the computed coefficients.

B:  [out] The coefficients of the original X with respect to the
computed basis. If B is a non-null pointer and B matches the
dimensions of B, then the coefficients computed during the
orthogonalization routine will be stored in B, similar to calling If B
points to a Teuchos::SerialDenseMatrix with size inconsistent with X,
then a std::invalid_argument exception will be thrown. Otherwise, if B
is null, the caller will not have access to the computed coefficients.
This matrix is not necessarily triangular (as in a QR factorization);
see the documentation of specific orthogonalization managers.  In
general, B has no non-zero structure.

Q:  [in] A list of multivector bases specifying the subspaces to be
orthogonalized against, satisfying \\\\[ \\\\langle Q[i], Q[j]
\\\\rangle = I \\\\quad\\\\textrm{if}\\\\quad i=j \\\\] and \\\\[
\\\\langle Q[i], Q[j] \\\\rangle = 0 \\\\quad\\\\textrm{if}\\\\quad i
\\\\neq j\\\\ . \\\\]

Rank of the basis computed by this method, less than or equal to the
number of columns in X. This specifies how many columns in the
returned X and rows in the returned B are valid. ";

/*  Error methods  */

%feature("docstring")  Anasazi::SVQBOrthoManager::orthonormErrorMat "Teuchos::ScalarTraits< ScalarType >::magnitudeType
Anasazi::SVQBOrthoManager< ScalarType, MV, OP
>::orthonormErrorMat(const MV &X, Teuchos::RCP< const MV >
MX=Teuchos::null) const

This method computes the error in orthonormality of a multivector,
measured as the Frobenius norm of the difference innerProd(X,Y) - I.
The method has the option of exploiting a caller-provided MX. ";

%feature("docstring")  Anasazi::SVQBOrthoManager::orthogErrorMat "Teuchos::ScalarTraits< ScalarType >::magnitudeType
Anasazi::SVQBOrthoManager< ScalarType, MV, OP >::orthogErrorMat(const
MV &X, const MV &Y, Teuchos::RCP< const MV > MX=Teuchos::null,
Teuchos::RCP< const MV > MY=Teuchos::null) const

This method computes the error in orthogonality of two multivectors,
measured as the Frobenius norm of innerProd(X,Y). The method has the
option of exploiting a caller-provided MX. ";


// File: classAnasazi_1_1TsqrAdaptor.xml
%feature("docstring") Anasazi::TsqrAdaptor "

Map from multivector class to TSQR adaptor class.

C++ includes: AnasaziTsqrAdaptor.hpp ";


// File: classAnasazi_1_1TsqrMatOrthoManager.xml
%feature("docstring") Anasazi::TsqrMatOrthoManager "

MatOrthoManager subclass using TSQR or SVQB.

When the inner product matrix has not been set, this class uses TSQR +
Block Gram-Schmidt (via  TsqrOrthoManagerImpl). If the inner product
matrix has been set, then this class uses the SVQB algorithm
(Stathopoulos and Wu 2002: CholeskyQR + SVD) for orthogonalization.

TSQR uses multivector scratch space. However, scratch space
initialization is \"lazy,\" so scratch space will not be allocated if
TSQR is not used.

C++ includes: AnasaziTsqrOrthoManager.hpp ";

%feature("docstring")
Anasazi::TsqrMatOrthoManager::TsqrMatOrthoManager "Anasazi::TsqrMatOrthoManager< Scalar, MV, OP
>::TsqrMatOrthoManager(const Teuchos::RCP< Teuchos::ParameterList >
&params, const std::string &label=\"Belos\", Teuchos::RCP< const OP >
Op=Teuchos::null)

Constructor (that sets user-specified parameters).

Parameters:
-----------

params:  [in/out] Configuration parameters, both for this
orthogonalization manager, and for TSQR itself (as the \"TSQR
implementation\" sublist). This can be null, in which case default
parameters will be set for now; you can always call setParameterList()
later to change these.

label:  [in] Label for timers. This only matters if the compile-time
option for enabling timers is set.

Op:  [in] Inner product with respect to which to orthogonalize
vectors. If Teuchos::null, use the Euclidean inner product.

Call  getValidParameters() for default parameters and their
documentation, including TSQR implementation parameters. Call
getFastParameters() to get documented parameters for faster
computation, possibly at the expense of accuracy and robustness. ";

%feature("docstring")
Anasazi::TsqrMatOrthoManager::TsqrMatOrthoManager "Anasazi::TsqrMatOrthoManager< Scalar, MV, OP
>::TsqrMatOrthoManager(const std::string &label=\"Belos\",
Teuchos::RCP< const OP > Op=Teuchos::null)

Constructor (that sets default parameters).

Parameters:
-----------

Op:  [in] Inner product with respect to which to orthogonalize
vectors. If Teuchos::null, use the Euclidean inner product.

label:  [in] Label for timers. This only matters if the compile-time
option for enabling timers is set. ";

%feature("docstring")
Anasazi::TsqrMatOrthoManager::~TsqrMatOrthoManager "virtual
Anasazi::TsqrMatOrthoManager< Scalar, MV, OP >::~TsqrMatOrthoManager()

Destructor (declared virtual for memory safety of derived classes). ";

%feature("docstring")
Anasazi::TsqrMatOrthoManager::getValidParameters "Teuchos::RCP<const
Teuchos::ParameterList> Anasazi::TsqrMatOrthoManager< Scalar, MV, OP
>::getValidParameters() const

Get default parameters for TsqrMatOrthoManager.

Get a (pointer to a) default list of parameters for configuring a
TsqrMatOrthoManager instance.

TSQR implementation configuration options are stored under \"TSQR
implementation\" as a sublist. ";

%feature("docstring")  Anasazi::TsqrMatOrthoManager::getFastParameters
"Teuchos::RCP<const Teuchos::ParameterList>
Anasazi::TsqrMatOrthoManager< Scalar, MV, OP >::getFastParameters()

Get \"fast\" parameters for TsqrMatOrthoManager.

Get a (pointer to a) list of parameters for configuring a
TsqrMatOrthoManager instance for maximum speed, at the cost of
accuracy (no block reorthogonalization) and robustness to rank
deficiency (no randomization of the null space basis).

TSQR implementation configuration options are stored under \"TSQR
implementation\" as a sublist. ";

%feature("docstring")  Anasazi::TsqrMatOrthoManager::setParameterList
"void Anasazi::TsqrMatOrthoManager< Scalar, MV, OP
>::setParameterList(const Teuchos::RCP< Teuchos::ParameterList >
&params) ";

%feature("docstring")  Anasazi::TsqrMatOrthoManager::setOp "virtual
void Anasazi::TsqrMatOrthoManager< Scalar, MV, OP
>::setOp(Teuchos::RCP< const OP > Op)

Set operator used for inner product. ";

%feature("docstring")  Anasazi::TsqrMatOrthoManager::getOp "Teuchos::RCP<const OP> Anasazi::TsqrMatOrthoManager< Scalar, MV, OP
>::getOp() const

Get operator used for inner product. ";

%feature("docstring")  Anasazi::TsqrMatOrthoManager::projectMat "void
Anasazi::TsqrMatOrthoManager< Scalar, MV, OP >::projectMat(MV &X,
Teuchos::Array< Teuchos::RCP< const MV > > Q, Teuchos::Array<
Teuchos::RCP< mat_type > > C=Teuchos::tuple(Teuchos::RCP< mat_type
>(Teuchos::null)), Teuchos::RCP< MV > MX=Teuchos::null,
Teuchos::Array< Teuchos::RCP< const MV > >
MQ=Teuchos::tuple(Teuchos::null)) const

Provides matrix-based projection method.

This method optionally allows the provision of $M X$ and/or the $M
Q[i]$. See OrthoManager::project() for more details.

Parameters:
-----------

X:  Q:  C:  [in/out] As in OrthoManager::project()

MX:  [in/out] If specified by the user, on input MX is required to be
the image of X under the operator getOp(). On output, MX will be
updated to reflect the changes in X.

MQ:  [in] If specified by the user, on MQ[i] is required to be the
image of Q[i] under the operator getOp(). ";

%feature("docstring")  Anasazi::TsqrMatOrthoManager::normalizeMat "int Anasazi::TsqrMatOrthoManager< Scalar, MV, OP >::normalizeMat(MV
&X, mat_ptr B=Teuchos::null, Teuchos::RCP< MV > MX=Teuchos::null)
const ";

%feature("docstring")
Anasazi::TsqrMatOrthoManager::projectAndNormalizeMat "int
Anasazi::TsqrMatOrthoManager< Scalar, MV, OP
>::projectAndNormalizeMat(MV &X, Teuchos::Array< Teuchos::RCP< const
MV > > Q, Teuchos::Array< Teuchos::RCP< mat_type > >
C=Teuchos::tuple(Teuchos::RCP< mat_type >(Teuchos::null)),
Teuchos::RCP< mat_type > B=Teuchos::null, Teuchos::RCP< MV >
MX=Teuchos::null, Teuchos::Array< Teuchos::RCP< const MV > >
MQ=Teuchos::tuple(Teuchos::RCP< const MV >(Teuchos::null))) const

Provides matrix-based projection/orthonormalization method.

This method optionally allows the provision of $M X$ and/or the $M
Q[i]$. See orthoManager::projectAndNormalize() for more details.

Parameters:
-----------

X:  Q:  C:  B:  [in/out] As in OrthoManager::projectAndNormalize()

MX:  [in/out] If specified by the user, on input MX is required to be
the image of X under the operator getOp(). On output, MX will be
updated to reflect the changes in X.

MQ:  [in] If specified by the user, on MQ[i] is required to be the
image of Q[i] under the operator getOp().

Rank of the basis computed by this method, less than or equal to the
number of columns in X. This specifies how many columns in the
returned X and MX and rows in the returned B are valid. ";

%feature("docstring")
Anasazi::TsqrMatOrthoManager::normalizeOutOfPlace "int
Anasazi::TsqrMatOrthoManager< Scalar, MV, OP >::normalizeOutOfPlace(MV
&X, MV &Q, mat_ptr B) const

Normalize X into Q*B.

Parameters:
-----------

X:  [in/out] On input: Multivector to normalize. On output: Possibly
overwritten with invalid values.

Q:  [out] On output: Normalized multivector.

B:  [out] On output: Normalization coefficients.

Rank of the input multivector X. ";

%feature("docstring")
Anasazi::TsqrMatOrthoManager::projectAndNormalizeOutOfPlace "int
Anasazi::TsqrMatOrthoManager< Scalar, MV, OP
>::projectAndNormalizeOutOfPlace(MV &X_in, MV &X_out, Teuchos::Array<
mat_ptr > C, mat_ptr B, Teuchos::ArrayView< Teuchos::RCP< const MV > >
Q) const

Project and normalize X_in into X_out.

Project X_in against Q, storing projection coefficients in C, and
normalize X_in into X_out, storing normalization coefficients in B. On
output, X_out has the resulting orthogonal vectors. X_in may be
overwritten with invalid values.

Parameters:
-----------

X_in:  [in/out] On input: The vectors to project against Q and
normalize. On output: possibly overwritten with invalid values.

X_out:  [out] The normalized input vectors after projection against Q.

C:  [out] Projection coefficients

B:  [out] Normalization coefficients

Q:  [in] The orthogonal basis against which to project

Rank of X_in after projection ";

%feature("docstring")  Anasazi::TsqrMatOrthoManager::orthonormErrorMat
"magnitude_type Anasazi::TsqrMatOrthoManager< Scalar, MV, OP
>::orthonormErrorMat(const MV &X, Teuchos::RCP< const MV >
MX=Teuchos::null) const

This method computes the error in orthonormality of a multivector.

This method optionally allows optionally exploits a caller-provided
MX. ";

%feature("docstring")  Anasazi::TsqrMatOrthoManager::orthogErrorMat "magnitude_type Anasazi::TsqrMatOrthoManager< Scalar, MV, OP
>::orthogErrorMat(const MV &X, const MV &Y, Teuchos::RCP< const MV >
MX=Teuchos::null, Teuchos::RCP< const MV > MY=Teuchos::null) const

This method computes the error in orthogonality of two multivectors.

This method optionally allows optionally exploits a caller-provided MX
and/or MY. ";


// File: classAnasazi_1_1TsqrOrthoError.xml
%feature("docstring") Anasazi::TsqrOrthoError "

TsqrOrthoManager(Impl) error.

Mark Hoemmen

C++ includes: AnasaziTsqrOrthoManagerImpl.hpp ";

%feature("docstring")  Anasazi::TsqrOrthoError::TsqrOrthoError "Anasazi::TsqrOrthoError::TsqrOrthoError(const std::string &what_arg)
";


// File: classAnasazi_1_1TsqrOrthoFault.xml
%feature("docstring") Anasazi::TsqrOrthoFault "

Orthogonalization fault.

Mark Hoemmen  Stewart (SISC 2008) presents a Block Gram-Schmidt (BGS)
algorithm with careful reorthogonalization. He defines an
\"orthogonalization fault\" as happening when the second BGS pass does
not succeed. This is possible in BGS, but not possible in (non-block)
Gram-Schmidt, if you use Stewart's randomization procedure for the
latter. Stewart gives an algorithm for recovering from an
orthogonalization fault, but the algorithm is expensive: it involves
careful reorthogonalization with non-block Gram-Schmidt. If the
\"throwOnReorthogFault\" option is set, we choose instead to report
the orthogonalization fault as an exception.

This is not a (subclass of) TsqrOrthoError, because the latter is a
logic or runtime bug, whereas a TsqrOrthoFault is a property of the
input and admits recovery.

C++ includes: AnasaziTsqrOrthoManagerImpl.hpp ";

%feature("docstring")  Anasazi::TsqrOrthoFault::TsqrOrthoFault "Anasazi::TsqrOrthoFault::TsqrOrthoFault(const std::string &what_arg)
";


// File: classAnasazi_1_1TsqrOrthoManager.xml
%feature("docstring") Anasazi::TsqrOrthoManager "

TSQR-based OrthoManager subclass.

Mark Hoemmen  Subclass of OrthoManager, implemented using
TsqrOrthoManagerImpl (TSQR + Block Gram-Schmidt).

C++ includes: AnasaziTsqrOrthoManager.hpp ";

%feature("docstring")  Anasazi::TsqrOrthoManager::setParameterList "void Anasazi::TsqrOrthoManager< Scalar, MV >::setParameterList(const
Teuchos::RCP< Teuchos::ParameterList > &params) ";

%feature("docstring")
Anasazi::TsqrOrthoManager::getNonconstParameterList "Teuchos::RCP<Teuchos::ParameterList> Anasazi::TsqrOrthoManager<
Scalar, MV >::getNonconstParameterList() ";

%feature("docstring")  Anasazi::TsqrOrthoManager::unsetParameterList "Teuchos::RCP<Teuchos::ParameterList> Anasazi::TsqrOrthoManager<
Scalar, MV >::unsetParameterList() ";

%feature("docstring")  Anasazi::TsqrOrthoManager::getValidParameters "Teuchos::RCP<const Teuchos::ParameterList> Anasazi::TsqrOrthoManager<
Scalar, MV >::getValidParameters() const

Default valid parameter list.

Get a (pointer to a) default list of parameters for configuring a
TsqrOrthoManager instance.

TSQR implementation configuration options are stored under
\"TsqrImpl\" as a sublist. ";

%feature("docstring")  Anasazi::TsqrOrthoManager::getFastParameters "Teuchos::RCP<const Teuchos::ParameterList> Anasazi::TsqrOrthoManager<
Scalar, MV >::getFastParameters() const

Get \"fast\" parameters for TsqrOrthoManager.

Get a (pointer to a) list of parameters for configuring a
TsqrOrthoManager instance for maximum speed, at the cost of accuracy
(no block reorthogonalization) and robustness to rank deficiency (no
randomization of the null space basis).

TSQR implementation configuration options are stored under
\"TsqrImpl\" as a sublist. ";

%feature("docstring")  Anasazi::TsqrOrthoManager::TsqrOrthoManager "Anasazi::TsqrOrthoManager< Scalar, MV >::TsqrOrthoManager(const
Teuchos::RCP< Teuchos::ParameterList > &params, const std::string
&label=\"Anasazi\")

Constructor (that sets user-specified parameters).

Parameters:
-----------

params:  [in/out] Configuration parameters, both for this
orthogonalization manager, and for TSQR itself (as the \"TsqrImpl\"
sublist). This can be null, in which case default parameters will be
set for now; you can always call setParameterList() later to change
these.

label:  [in] Label for timers. This only matters if the compile-time
option for enabling timers is set.

Call  getValidParameters() for default parameters and their
documentation, including TSQR implementation parameters. Call
getFastParameters() to get documented parameters for faster
computation, possibly at the expense of accuracy and robustness. ";

%feature("docstring")  Anasazi::TsqrOrthoManager::TsqrOrthoManager "Anasazi::TsqrOrthoManager< Scalar, MV >::TsqrOrthoManager(const
std::string &label)

Constructor (that sets default parameters).

Parameters:
-----------

label:  [in] Label for timers. This only matters if the compile-time
option for enabling timers is set. ";

%feature("docstring")  Anasazi::TsqrOrthoManager::~TsqrOrthoManager "virtual Anasazi::TsqrOrthoManager< Scalar, MV >::~TsqrOrthoManager()

Destructor, declared virtual for safe inheritance. ";

%feature("docstring")  Anasazi::TsqrOrthoManager::innerProd "void
Anasazi::TsqrOrthoManager< Scalar, MV >::innerProd(const MV &X, const
MV &Y, mat_type &Z) const

Provides the inner product defining the orthogonality concepts.

All concepts of orthogonality discussed in this class are defined with
respect to this inner product.

This is potentially different from MultiVecTraits::MvTransMv(). For
example, it is customary in many eigensolvers to exploit a mass matrix
M for the inner product: $x^HMx$.

Parameters:
-----------

Z:  [out] Z(i,j) contains the inner product of X[i] and Y[i]: \\\\[
Z(i,j) = \\\\langle X[i], Y[i] \\\\rangle \\\\] ";

%feature("docstring")  Anasazi::TsqrOrthoManager::norm "void
Anasazi::TsqrOrthoManager< Scalar, MV >::norm(const MV &X,
std::vector< magnitude_type > &normVec) const ";

%feature("docstring")  Anasazi::TsqrOrthoManager::project "void
Anasazi::TsqrOrthoManager< Scalar, MV >::project(MV &X,
Teuchos::Array< Teuchos::RCP< const MV > > Q, Teuchos::Array<
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, Scalar > > >
C=Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix< int, Scalar
> >(Teuchos::null))) const

Given a list of mutually orthogonal and internally orthonormal bases
Q, this method projects a multivector X onto the space orthogonal to
the individual Q[i], optionally returning the coefficients of X for
the individual Q[i]. All of this is done with respect to the inner
product innerProd().

After calling this routine, X will be orthogonal to each of the Q[i].

Parameters:
-----------

X:  [in/out] The multivector to be modified.  On output, the columns
of X will be orthogonal to each Q[i], satisfying \\\\[ \\\\langle
Q[i], X_{out} \\\\rangle = 0 \\\\] Also, \\\\[ X_{out} = X_{in} -
\\\\sum_i Q[i] \\\\langle Q[i], X_{in} \\\\rangle \\\\]

Q:  [in] A list of multivector bases specifying the subspaces to be
orthogonalized against, satisfying \\\\[ \\\\langle Q[i], Q[j]
\\\\rangle = I \\\\quad\\\\textrm{if}\\\\quad i=j \\\\] and \\\\[
\\\\langle Q[i], Q[j] \\\\rangle = 0 \\\\quad\\\\textrm{if}\\\\quad i
\\\\neq j\\\\ . \\\\]

C:  [out] The coefficients of X in the bases Q[i]. If C[i] is a non-
null pointer and C[i] matches the dimensions of X and Q[i], then the
coefficients computed during the orthogonalization routine will be
stored in the matrix C[i], similar to calling If C[i] points to a
Teuchos::SerialDenseMatrix with size inconsistent with X and  Q[i],
then a std::invalid_argument exception will be thrown.  Otherwise, if
C.size() < i or C[i] is a null pointer, the caller will not have
access to the computed coefficients. ";

%feature("docstring")  Anasazi::TsqrOrthoManager::normalize "int
Anasazi::TsqrOrthoManager< Scalar, MV >::normalize(MV &X, mat_ptr
B=Teuchos::null) const ";

%feature("docstring")  Anasazi::TsqrOrthoManager::projectAndNormalize
"int Anasazi::TsqrOrthoManager< Scalar, MV >::projectAndNormalize(MV
&X, Teuchos::Array< Teuchos::RCP< const MV > > Q, Teuchos::Array<
Teuchos::RCP< Teuchos::SerialDenseMatrix< int, Scalar > > >
C=Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix< int, Scalar
> >(Teuchos::null)), Teuchos::RCP< Teuchos::SerialDenseMatrix< int,
Scalar > > B=Teuchos::null) const

Given a set of bases Q[i] and a multivector X, this method computes an
orthonormal basis for $colspan(X) - \\\\sum_i colspan(Q[i])$.

This routine returns an integer rank stating the rank of the computed
basis. If the subspace $colspan(X) - \\\\sum_i colspan(Q[i])$ does not
have dimension as large as the number of columns of X and the
orthogonalization manager does not attempt to augment the subspace,
then rank may be smaller than the number of columns of X. In this
case, only the first rank columns of output X and first rank rows of B
will be valid.

This routine guarantees both the orthogonality of the returned basis
against the Q[i] as well as the orthonormality of the returned basis.
Therefore, this method is not necessarily equivalent to calling
project() followed by a call to normalize(); see the documentation for
specific orthogonalization managers.

Parameters:
-----------

X:  [in/out] On output, the first rank columns of X satisfy \\\\[
\\\\langle X[i], X[j] \\\\rangle = \\\\delta_{ij} \\\\quad
\\\\textrm{and} \\\\quad \\\\langle X, Q[i] \\\\rangle = 0\\\\ . \\\\]
Also, \\\\[ X_{in}(1:m,1:n) = X_{out}(1:m,1:rank) B(1:rank,1:n) +
\\\\sum_i Q[i] C[i] \\\\] where m is the number of rows in X and n is
the number of columns in X.

Q:  [in] A list of multivector bases specifying the subspaces to be
orthogonalized against, satisfying \\\\[ \\\\langle Q[i], Q[j]
\\\\rangle = I \\\\quad\\\\textrm{if}\\\\quad i=j \\\\] and \\\\[
\\\\langle Q[i], Q[j] \\\\rangle = 0 \\\\quad\\\\textrm{if}\\\\quad i
\\\\neq j\\\\ . \\\\]

C:  [out] The coefficients of X in the Q[i]. If C[i] is a non-null
pointer and C[i] matches the dimensions of X and Q[i], then the
coefficients computed during the orthogonalization routine will be
stored in the matrix C[i], similar to calling If C[i] points to a
Teuchos::SerialDenseMatrix with size inconsistent with X and  Q[i],
then a std::invalid_argument exception will be thrown.  Otherwise, if
C.size() < i or C[i] is a null pointer, the caller will not have
access to the computed coefficients.

B:  [out] The coefficients of the original X with respect to the
computed basis. If B is a non-null pointer and B matches the
dimensions of B, then the coefficients computed during the
orthogonalization routine will be stored in B, similar to calling If B
points to a Teuchos::SerialDenseMatrix with size inconsistent with X,
then a std::invalid_argument exception will be thrown.  Otherwise, if
B is null, the caller will not have access to the computed
coefficients.

This matrix is not necessarily triangular (as in a QR factorization);
see the documentation of specific orthogonalization managers.

Rank of the basis computed by this method, less than or equal to the
number of columns in X. This specifies how many columns in the
returned X and rows in the returned B are valid. ";

%feature("docstring")  Anasazi::TsqrOrthoManager::normalizeOutOfPlace
"int Anasazi::TsqrOrthoManager< Scalar, MV >::normalizeOutOfPlace(MV
&X, MV &Q, mat_ptr B) const

Normalize X into Q*B, overwriting X with invalid values.

We expose this interface to applications because TSQR is not able to
compute an orthogonal basis in place; it needs scratch space.
Applications can exploit this interface to avoid excessive copying of
vectors when using TSQR for orthogonalization.

Parameters:
-----------

X:  [in/out] Input vector(s) to normalize

Q:  [out] Normalized output vector(s)

B:  [out] Normalization coefficients

Rank of X.

Q must have at least as many columns as X. It may have more columns
than X; those columns are ignored. ";

%feature("docstring")
Anasazi::TsqrOrthoManager::projectAndNormalizeOutOfPlace "int
Anasazi::TsqrOrthoManager< Scalar, MV
>::projectAndNormalizeOutOfPlace(MV &X_in, MV &X_out, Teuchos::Array<
mat_ptr > C, mat_ptr B, Teuchos::ArrayView< Teuchos::RCP< const MV > >
Q) const

Project and normalize X_in into X_out; overwrite X_in.

Project X_in against Q, storing projection coefficients in C, and
normalize X_in into X_out, storing normalization coefficients in B. On
output, X_out has the resulting orthogonal vectors and X_in is
overwritten with invalid values.

Parameters:
-----------

X_in:  [in/out] On input: The vectors to project against Q and
normalize. Overwritten with invalid values on output.

X_out:  [out] The normalized input vectors after projection against Q.

C:  [out] Projection coefficients

B:  [out] Normalization coefficients

Q:  [in] The orthogonal basis against which to project

Rank of X_in after projection

We expose this interface to applications for the same reason that we
expose  normalizeOutOfPlace(). ";

%feature("docstring")  Anasazi::TsqrOrthoManager::orthonormError "magnitude_type Anasazi::TsqrOrthoManager< Scalar, MV
>::orthonormError(const MV &X) const

This method computes the error in orthonormality of a multivector.

This method return some measure of $\\\\| \\\\langle X, X \\\\rangle -
I \\\\| $.  See the documentation of specific orthogonalization
managers. ";

%feature("docstring")  Anasazi::TsqrOrthoManager::orthogError "magnitude_type Anasazi::TsqrOrthoManager< Scalar, MV
>::orthogError(const MV &X1, const MV &X2) const

This method computes the error in orthogonality of two multivectors.

This method return some measure of $\\\\| \\\\langle X1, X2 \\\\rangle
\\\\| $.  See the documentation of specific orthogonalization
managers. ";


// File: classAnasazi_1_1TsqrOrthoManagerImpl.xml
%feature("docstring") Anasazi::TsqrOrthoManagerImpl "

TSQR-based OrthoManager subclass implementation.

Mark Hoemmen  TsqrOrthoManagerImpl implements the interface defined by
OrthoManager, as well as the interface defined by
OutOfPlaceNormalizerMixin. We use TsqrOrthoManagerImpl to implement
TsqrOrthoManager and  TsqrMatOrthoManager.

Parameters:
-----------

Scalar:  The type of matrix and (multi)vector entries.

MV:  The type of (multi)vector inputs and outputs.

This class uses a combination of Tall Skinny QR (TSQR) and Block Gram-
Schmidt (BGS) to orthogonalize multivectors. The Block Gram- Schmidt
procedure used here is inspired by that of G. W. Stewart (\"Block
Gram-Schmidt Orthogonalization\", SISC vol 31 #1 pp. 761-- 775, 2008).
The difference is that we use TSQR+SVD instead of Stewart's careful
Gram-Schmidt with reorthogonalization to handle the current block.
\"Orthogonalization faults\" (as defined by Stewart) may still happen,
but we do not handle them by default. Rather, we make one BGS pass, do
TSQR+SVD, check the resulting column norms, and make a second BGS pass
(+ TSQR+SVD) if necessary. If we then detect an orthogonalization
fault, we throw  TsqrOrthoFault.

Despite the \"Impl\" part of the name of this class, we don't actually
use it for the \"pImpl\" C++ idiom. We just separate out the TSQR
implementation to make it easier to implement the OrthoManager and
MatOrthoManager interfaces for the case where the inner product
operator is not the identity matrix.

C++ includes: AnasaziTsqrOrthoManagerImpl.hpp ";

%feature("docstring")
Anasazi::TsqrOrthoManagerImpl::getValidParameters "Teuchos::RCP<
const Teuchos::ParameterList > Anasazi::TsqrOrthoManagerImpl< Scalar,
MV >::getValidParameters() const

Default valid parameter list.

Get a (pointer to a) default list of parameters for configuring a
TsqrOrthoManagerImpl instance.

TSQR implementation configuration options are stored under \"TSQR
implementation\" as a sublist. ";

%feature("docstring")  Anasazi::TsqrOrthoManagerImpl::setParameterList
"void Anasazi::TsqrOrthoManagerImpl< Scalar, MV
>::setParameterList(const Teuchos::RCP< Teuchos::ParameterList >
&params)

Set parameters from the given parameter list. ";

%feature("docstring")
Anasazi::TsqrOrthoManagerImpl::getFastParameters "Teuchos::RCP< const
Teuchos::ParameterList > Anasazi::TsqrOrthoManagerImpl< Scalar, MV
>::getFastParameters()

Get \"fast\" parameters for TsqrOrthoManagerImpl.

Get a (pointer to a) list of parameters for configuring a
TsqrOrthoManager or TsqrMatOrthoManager instance for maximum speed, at
the cost of accuracy (no block reorthogonalization) and robustness to
rank deficiency (no randomization of the null space basis).

TSQR implementation configuration options are stored under \"TSQR
implementation\" as a sublist. ";

%feature("docstring")
Anasazi::TsqrOrthoManagerImpl::TsqrOrthoManagerImpl "Anasazi::TsqrOrthoManagerImpl< Scalar, MV
>::TsqrOrthoManagerImpl(const Teuchos::RCP< Teuchos::ParameterList >
&params, const std::string &label)

Constructor (that sets user-specified parameters).

Parameters:
-----------

params:  [in/out] Configuration parameters, both for this
orthogonalization manager, and for TSQR itself (as the \"TSQR
implementation\" sublist). This can be null, in which case default
parameters will be set for now; you can always call
setParameterList() later to change these.

label:  [in] Label for timers. This only matters if the compile-time
option for enabling timers is set.

Call  getValidParameters() for default parameters and their
documentation, including TSQR implementation parameters. Call
getFastParameters() to get documented parameters for faster
computation, possibly at the expense of accuracy and robustness. ";

%feature("docstring")
Anasazi::TsqrOrthoManagerImpl::TsqrOrthoManagerImpl "Anasazi::TsqrOrthoManagerImpl< Scalar, MV
>::TsqrOrthoManagerImpl(const std::string &label)

Constructor (that sets default parameters).

Parameters:
-----------

label:  [in] Label for timers. This only matters if the compile-time
option for enabling timers is set. ";

%feature("docstring")  Anasazi::TsqrOrthoManagerImpl::setLabel "void
Anasazi::TsqrOrthoManagerImpl< Scalar, MV >::setLabel(const
std::string &label)

Set the label for timers.

This only matters if timers are enabled. If timers are enabled and the
label changes, this method will clear the old timers and replace them
with new ones. The old timers will not appear in the list of timers
shown by Teuchos::TimeMonitor::summarize(). ";

%feature("docstring")  Anasazi::TsqrOrthoManagerImpl::getLabel "const
std::string& Anasazi::TsqrOrthoManagerImpl< Scalar, MV >::getLabel()
const

Get the label for timers (if timers are enabled). ";

%feature("docstring")  Anasazi::TsqrOrthoManagerImpl::innerProd "void
Anasazi::TsqrOrthoManagerImpl< Scalar, MV >::innerProd(const MV &X,
const MV &Y, mat_type &Z) const

Euclidean inner product.

Compute the Euclidean block inner product X^* Y, and store the result
in Z.

Parameters:
-----------

X:  [in]

Y:  [in]

Z:  [out] On output, $X^* Y$ ";

%feature("docstring")  Anasazi::TsqrOrthoManagerImpl::norm "void
Anasazi::TsqrOrthoManagerImpl< Scalar, MV >::norm(const MV &X,
std::vector< magnitude_type > &normvec) const

Compute the 2-norm of each column j of X.

Parameters:
-----------

X:  [in] Multivector for which to compute column norms.

normVec:  [out] On output: normvec[j] is the 2-norm of column j of X.
normVec is resized if necessary so that it has at least as many
entries as there are columns of X.

Performance of this method depends on how MultiVecTraits implements
column norm computation for the given multivector type MV. It may or
may not be the case that a reduction is performed for every column of
X. Furthermore, whether or not the columns of X are contiguous (as
opposed to a view of noncontiguous columns) may also affect
performance. The computed results should be the same regardless,
except perhaps for small rounding differences due to a different order
of operations. ";

%feature("docstring")  Anasazi::TsqrOrthoManagerImpl::project "void
Anasazi::TsqrOrthoManagerImpl< Scalar, MV >::project(MV &X,
Teuchos::Array< mat_ptr > C, Teuchos::ArrayView< Teuchos::RCP< const
MV > > Q)

Compute $C := Q^* X$ and $X := X - Q C$.

Project X against the span of the (Euclidean) orthogonal vectors Q,
and store the resulting coefficients in C.

Parameters:
-----------

X:  [in/out] On input: the vectors to project. On output: $X := X - Q
C$ where $C := Q^* X$.

C:  [out] The projection coefficients $C := Q^* X$

Q:  [in] The orthogonal basis against which to project ";

%feature("docstring")  Anasazi::TsqrOrthoManagerImpl::normalize "int
Anasazi::TsqrOrthoManagerImpl< Scalar, MV >::normalize(MV &X, mat_ptr
B)

Orthogonalize the columns of X in place.

Orthogonalize the columns of X in place, storing the resulting
coefficients in B. Return the rank of X. If X is full rank, then X*B
on output is a QR factorization of X on input. If X is not full rank,
then the first rank columns of X on output form a basis for the column
space of X (on input). Additional options control randomization of the
null space basis.

Parameters:
-----------

X:  [in/out]

B:  [out]

Rank of X ";

%feature("docstring")
Anasazi::TsqrOrthoManagerImpl::normalizeOutOfPlace "int
Anasazi::TsqrOrthoManagerImpl< Scalar, MV >::normalizeOutOfPlace(MV
&X, MV &Q, mat_ptr B)

Normalize X into Q*B, overwriting X.

Normalize X into Q*B, overwriting X with invalid values.

Parameters:
-----------

X:  [in/out] Vector(s) to normalize

Q:  [out] Normalized vector(s)

B:  [out] Normalization coefficients

Rank of X

Q must have at least as many columns as X. It may have more columns
than X; those columns are ignored.

We expose this interface to applications because TSQR is not able to
compute an orthogonal basis in place; it needs scratch space.
Applications can exploit this interface to avoid excessive copying of
vectors when using TSQR for orthogonalization. ";

%feature("docstring")
Anasazi::TsqrOrthoManagerImpl::projectAndNormalize "int
Anasazi::TsqrOrthoManagerImpl< Scalar, MV >::projectAndNormalize(MV
&X, Teuchos::Array< mat_ptr > C, mat_ptr B, Teuchos::ArrayView<
Teuchos::RCP< const MV > > Q)

Project X against Q and normalize X.

This method is equivalent (in exact arithmetic) to project(X,C,Q)
followed by normalize(X,B). However, the interface allows this method
to implement reorthogonalization more efficiently and accurately.

Parameters:
-----------

X:  [in/out] The vectors to project against Q and normalize

C:  [out] The projection coefficients

B:  [out] The normalization coefficients

Q:  [in] The orthogonal basis against which to project

Rank of X after projection ";

%feature("docstring")
Anasazi::TsqrOrthoManagerImpl::projectAndNormalizeOutOfPlace "int
Anasazi::TsqrOrthoManagerImpl< Scalar, MV
>::projectAndNormalizeOutOfPlace(MV &X_in, MV &X_out, Teuchos::Array<
mat_ptr > C, mat_ptr B, Teuchos::ArrayView< Teuchos::RCP< const MV > >
Q)

Project and normalize X_in into X_out; overwrite X_in.

Project X_in against Q, storing projection coefficients in C, and
normalize X_in into X_out, storing normalization coefficients in B. On
output, X_out has the resulting orthogonal vectors and X_in is
overwritten with invalid values.

Parameters:
-----------

X_in:  [in/out] On input: The vectors to project against Q and
normalize. Overwritten with invalid values on output.

X_out:  [out] On output: the normalized input vectors after projection
against Q.

C:  [out] The projection coefficients

B:  [out] The normalization coefficients

Q:  [in] The orthogonal basis against which to project

Rank of X_in after projection

We expose this interface to applications for the same reason that we
expose  normalizeOutOfPlace(). ";

%feature("docstring")  Anasazi::TsqrOrthoManagerImpl::orthonormError "magnitude_type Anasazi::TsqrOrthoManagerImpl< Scalar, MV
>::orthonormError(const MV &X) const

Return $ \\\\| I - X^* \\\\cdot X \\\\|_F $.

Return the Frobenius norm of I - X^* X, which is an absolute measure
of the orthogonality of the columns of X. ";

%feature("docstring")  Anasazi::TsqrOrthoManagerImpl::orthogError "magnitude_type Anasazi::TsqrOrthoManagerImpl< Scalar, MV
>::orthogError(const MV &X1, const MV &X2) const

Return the Frobenius norm of the inner product of X1 with itself. ";

%feature("docstring")
Anasazi::TsqrOrthoManagerImpl::blockReorthogThreshold "magnitude_type
Anasazi::TsqrOrthoManagerImpl< Scalar, MV >::blockReorthogThreshold()
const

Relative tolerance for triggering a block reorthogonalization. If any
column norm in a block decreases by this amount, then we
reorthogonalize. ";

%feature("docstring")
Anasazi::TsqrOrthoManagerImpl::relativeRankTolerance "magnitude_type
Anasazi::TsqrOrthoManagerImpl< Scalar, MV >::relativeRankTolerance()
const

Relative tolerance for determining (via the SVD) whether a block is of
full numerical rank. ";


// File: structAnasazi_1_1UndefinedDenseMatTraits.xml
%feature("docstring") Anasazi::UndefinedDenseMatTraits "

This is the default struct used by DenseMatrixTraits<OrdinalType,
ScalarType> class to produce a compile time error when the
specialization does not exist for dense matrix type DM.

C++ includes: AnasaziDenseMatTraits.hpp ";


// File: structAnasazi_1_1UndefinedMultiVecTraits.xml
%feature("docstring") Anasazi::UndefinedMultiVecTraits "

Used by MultiVecTraits to report lack of a specialization.

MultiVecTraits<ScalarType, MV> uses this struct to produce a compile-
time error when no specialization exists for the scalar type
ScalarType and multivector type MV.

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
/*  Anasazi Enumerations  */

%feature("docstring")  Anasazi::details::Anasazi_Version "std::string
Anasazi::Anasazi_Version() ";

%feature("docstring")  Anasazi::details::TestMultiVecTraits "bool
Anasazi::TestMultiVecTraits(const Teuchos::RCP< OutputManager<
ScalarType > > &om, const Teuchos::RCP< const MV > &A)

This is a function to test the correctness of a MultiVecTraits
specialization and multivector implementation.

Status of the test: true is success, false is error ";

%feature("docstring")  Anasazi::details::TestOperatorTraits "bool
Anasazi::TestOperatorTraits(const Teuchos::RCP< OutputManager<
ScalarType > > &om, const Teuchos::RCP< const MV > &A, const
Teuchos::RCP< const OP > &M)

This function tests the correctness of an operator implementation with
respect to an OperatorTraits specialization.

Status of the test: true is successful, false otherwise. ";


// File: namespaceAnasazi_1_1details.xml


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


// File: AnasaziGenOrthoManager_8hpp.xml


// File: AnasaziHelperTraits_8hpp.xml


// File: AnasaziICGSOrthoManager_8hpp.xml


// File: AnasaziIRTR_8hpp.xml


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


// File: AnasaziRTRBase_8hpp.xml


// File: AnasaziRTRSolMgr_8hpp.xml


// File: AnasaziSimpleLOBPCGSolMgr_8hpp.xml


// File: AnasaziSIRTR_8hpp.xml


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


// File: AnasaziStubTsqrAdapter_8hpp.xml


// File: AnasaziSVQBOrthoManager_8hpp.xml


// File: AnasaziTsqrAdaptor_8hpp.xml


// File: AnasaziTsqrOrthoManager_8hpp.xml


// File: AnasaziTsqrOrthoManagerImpl_8hpp.xml


// File: AnasaziTypes_8hpp.xml


// File: AnasaziVersion_8cpp.xml


// File: dir_cfcf785d330e8f6cf9265ee0cf028086.xml


// File: dir_73432c73b202f8663e3fe783990d1240.xml


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



// File: index.xml

// File: classKomplex__LinearProblem.xml
%feature("docstring") Komplex_LinearProblem "

Komplex_LinearProblem: A class for forming an equivalent real
formulation of a complex valued problem.

The Komplex_LinearProblem class takes a complex linear problem,
separated into real and imaginary parts, and forms an equivalent real
valued system of twice the dimension. The resulting system can then be
solved with any Trilinos solver that understands Epetra objects.

KOMPLEX solves a complex-valued linear system Ax = b by solving an
equivalent real-valued system of twice the dimension. Specifically,
writing in terms of real and imaginary parts, we have

\\\\[ (A_r + i*A_i)*(x_r + i*x_i) = (b_r + i*b_i) \\\\]

or by separating into real and imaginary equations we have

\\\\[ \\\\left( \\\\begin{array}{rr} A_r & -A_i\\\\\\\\ A_i & A_r
\\\\end{array} \\\\right) \\\\left( \\\\begin{array}{r} x_r\\\\\\\\
x_i \\\\end{array} \\\\right) = \\\\left( \\\\begin{array}{r}
b_r\\\\\\\\ b_i \\\\end{array} \\\\right) \\\\] which is a real-valued
system of twice the size. If we find xr and xi, we can form the
solution to the original system as x = xr +i*xi.

KOMPLEX accepts the user linear system as two real-valued matrices
with no assumption about the structure of the matrices, except that
they have compatible RowMap, DomainMap and RangeMap distributions.
Each matrix is multiplied by user-supplied complex constants.

Although formally the system is a 2-by-2 block system, we actually
apply the interleaving at the matrix entry level such that the real
part of the first complex equation is followed by the imaginary part
of the first complex equation, and so on. This approach is documented
in:

David Day and Michael A. Heroux. Solving complex-valued linear systems
via equivalent real formulations. SIAM J. Sci. Comput., 23(2):480498,
2001.

C++ includes: Komplex_LinearProblem.h ";

%feature("docstring")  Komplex_LinearProblem::Komplex_LinearProblem "Komplex_LinearProblem::Komplex_LinearProblem(double c0r, double c0i,
const Epetra_RowMatrix &A0, double c1r, double c1i, const
Epetra_RowMatrix &A1, const Epetra_MultiVector &Xr, const
Epetra_MultiVector &Xi, const Epetra_MultiVector &Br, const
Epetra_MultiVector &Bi)

Komplex_LinearProblem constructor.

Constructs the Komplex operator from the user definition of the
complex-valued matrix C = (c0r+i*c0i)*A0 +(c1r+i*c1i)*A1. Using this
general expression for the complex matrix allows easy formulation of a
variety of common complex problems.

Parameters:
-----------

c0r:  (In) The real part of the complex coefficient multiplying A0.

c0i:  (In) The imag part of the complex coefficient multiplying A0.

A0:  (In) An Epetra_RowMatrix that is one of the matrices used to
define the true complex operator.

c1r:  (In) The real part of the complex coefficient multiplying A1.

c1i:  (In) The imag part of the complex coefficient multiplying A1.

A1:  (In) An Epetra_RowMatrix that is the second of the matrices used
to define the true complex operator.

Xr:  (In) The real part of the complex valued LHS.

Xi:  (In) The imag part of the complex valued LHS.

Br:  (In) The real part of the complex valued RHS.

Bi:  (In) The imag part of the complex valued RHS. ";

%feature("docstring")  Komplex_LinearProblem::~Komplex_LinearProblem "Komplex_LinearProblem::~Komplex_LinearProblem()

Komplex_LinearProblem Destructor. ";

%feature("docstring")  Komplex_LinearProblem::UpdateValues "int
Komplex_LinearProblem::UpdateValues(double c0r, double c0i, const
Epetra_RowMatrix &A0, double c1r, double c1i, const Epetra_RowMatrix
&A1, const Epetra_MultiVector &Xr, const Epetra_MultiVector &Xi, const
Epetra_MultiVector &Br, const Epetra_MultiVector &Bi)

Update the values of the equivalent real valued system.

This method allows the values of an existing Komplex_LinearProblem
object to be updated. Note that the update that there is no change to
the pattern of the matrices.

Error code, set to 0 if no error. ";

%feature("docstring")  Komplex_LinearProblem::ExtractSolution "int
Komplex_LinearProblem::ExtractSolution(Epetra_MultiVector &Xr,
Epetra_MultiVector &Xi)

Extrac a solution for the original complex-valued problem using the
solution of the Komplex problem.

After solving the komplex linear system, this method can be called to
extract the solution of the original problem, assuming the solution
for the komplex system is valid.

Parameters:
-----------

Xr:  (Out) An existing Epetra_MultiVector. On exit it will contain the
real part of the complex valued solution.

Xi:  (Out) An existing Epetra_MultiVector. On exit it will contain the
imag part of the complex valued solution. ";

%feature("docstring")  Komplex_LinearProblem::KomplexProblem "Epetra_LinearProblem* Komplex_LinearProblem::KomplexProblem() const

Returns pointer to the Epetra_LinearProblem object that defines the
Komplex formulation.

The pointer returned from this method will contain the address of a
fully-constructed Epetra_LinearProblem instance that can be used with
any Trilinos preconditioner or solver. ";


// File: Komplex__LinearProblem_8cpp.xml


// File: Komplex__LinearProblem_8h.xml


// File: Komplex__Version_8h.xml
%feature("docstring")  Komplex_Version "string Komplex_Version() ";


// File: dir_5dde9cc40b3bbb9db93e3e148c874a70.xml


// File: dir_6f7b992297991ea7bb945c55cf30835c.xml


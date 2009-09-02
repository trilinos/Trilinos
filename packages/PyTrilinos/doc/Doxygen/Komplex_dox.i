
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

KOMPLEX accept user linear systems in three forms with either global
or local index values.

1) The first form is true complex. The user passes in an MSR or VBR
format matrix where the values are stored like Fortran complex
numbers. Thus, the values array is of type double that is twice as
long as the number of complex values. Each complex entry is stored
with real part followed by imaginary part (as in Fortran).

2) The second form stores real and imaginary parts separately, but the
pattern for each is identical. Thus only the values of the imaginary
part are passed to the creation routines.

3) The third form accepts two real-valued matrices with no assumption
about the structure of the matrices. Each matrix is multiplied by a
user-supplied complex constant. This is the most general form.

Each of the above forms supports a global or local index set. By this
we mean that the index values (stored in bindx) refer to the global
problem indices, or the local indices (for example after calling
AZ_transform).

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

The operator will be explicitly constructed as an Epetra_VbrMatrix
object when the first call to SetKomplexOperator() is made. Subsequent
calls to this method will attempt to reuse the the existing
KomplexVbrMatrix object if possible, rather than reconstructing from
scratch. If this is not possible (typically because the structure has
changed) then a the previous KomplexVbrMatrix object will be deleted
and a new one will be constructed.

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

Bi:  (In) The imag part of the complex valued RHS.

Error code, set to 0 if no error. ";

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


// File: dir_7c33e1574cd31383e666b98827d697d7.xml


// File: dir_f84685b04efc5d362c9b262a426031b6.xml


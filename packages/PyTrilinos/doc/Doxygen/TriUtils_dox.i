
// File: index.xml

// File: classTrilinos__Util_1_1CommandLineParser.xml
%feature("docstring") Trilinos_Util::CommandLineParser "

Trilinos_Util::CommandLineParser: A class for managing the input
arguments and variables.

Using Trilinos_Util::CommandLineParser, it is easy to handle input
line arguments and shell varibles. For instance, the user can write $
./a.out -nx 10 -tol 1e-6 -solver=cg and then easily retrive the value
of nx, tol, and solver.

A simple code using this class is as follows: int main(int argc, char
*argv[])   { Trilinos_Util::CommandLineParser CLP(argc,argv);    int
nx = CLP.GetInt(\"-nx\", 123);    int ny = CLP.GetInt(\"-ny\", 145);
double tol = CLP.GetDouble(\"-tol\", 1e-12);    string solver =
CLP.GetInt(\"-solver\");     cout << \"nx = \" << nx << endl;    cout
<< \"ny = \" << ny << \" (default value)\" << endl;    cout << \"tol =
\" << tol << endl;    cout << \"solver = \" << solver << endl; return
0;      }

Each line option can have a value or not. For options with a value,
the user can specify this values as follows. Let -tolerance be the
name of the option and 1e-12 its value. Both choices are valid:
-option 1e-12 (with one or more spaces)

-option=1e-12 (an `=' sign and no spaces)

Options are indentified with one or more dashes (`-'). Each option
cannot have more than one value.

Note that the user can specify some values without giving them a name.
This can be done as follows:  $ ./a.out value1 value2 value 3 -nx 10
-tol 1e-6 -solver=cg Here, valueX, (X=1,...,9) is stored in the
database entry ARGV_X.

To use this class, the user has to build the database using the
argc,argv input arguments. Then, to retrive the option value, the user
has to use one of the following functions: GetInt

GetDouble

GetString

If option name is not found in the database, a value of 0, 0.0 or an
empty string is returned. If needed, the user can also specify a
default value to return when the option name is not found in the
database. Method HaveOption can be used to query the database for an
option.

The user can modify the database as well, using Set

Add

(GetInt, GetDouble, GetString, Set and Add are derived from the base
class, Trilinos_Util_Map).

Finally, the user can retrive the integer, double or string value of a
shell environmental variable using: GetIntShellVariable

GetDoubleShellVariable

GetStringShellVariable

Marzio Sala, SNL 9214

C++ includes: Trilinos_Util_CommandLineParser.h ";

%feature("docstring")
Trilinos_Util::CommandLineParser::CommandLineParser "Trilinos_Util::CommandLineParser::CommandLineParser(int argc, char
*argv[])

Trilinos_Util_ShellOptions constructor using the options given at the
shell line. ";

%feature("docstring")
Trilinos_Util::CommandLineParser::GetProgramName "string
Trilinos_Util::CommandLineParser::GetProgramName(void)

Returns the name of the program as a C++ string. ";

%feature("docstring")
Trilinos_Util::CommandLineParser::GetIntShellVariable "int
Trilinos_Util::CommandLineParser::GetIntShellVariable(const char *str)

Returns the value of the environmenta variable str as an integer.

This methods returns the value of the environmental variable str. If
the variable does not exists, returns 0. ";

%feature("docstring")
Trilinos_Util::CommandLineParser::GetDoubleShellVariable "double
Trilinos_Util::CommandLineParser::GetDoubleShellVariable(const char
*str)

Returns the value of the environmenta variable str as an double.

This methods returns the value of the environmenta variable str. If
the variable does not exists, returns 0.0. ";

%feature("docstring")
Trilinos_Util::CommandLineParser::GetStringShellVariable "string
Trilinos_Util::CommandLineParser::GetStringShellVariable(const char
*str)

Returns the value of the environmenta variable str as a C++ string.

This methods returns the value of the environmenta variable str. If
the variable does not exists, returns \"\". ";


// File: classTrilinos__Util_1_1CrsMatrixGallery.xml
%feature("docstring") Trilinos_Util::CrsMatrixGallery "";

%feature("docstring")
Trilinos_Util::CrsMatrixGallery::CrsMatrixGallery "Trilinos_Util::CrsMatrixGallery::CrsMatrixGallery(const string name,
const Epetra_Comm &comm)

Triutils_Gallery Constructor.

Creates a Triutils_Gallery instance.

The first parameter is the name of the matrix. We refer to the
Trilinos Tutorial for a detailed description of available matrices.

The matrix name can be empty (\"\"), and set later using, for example,
Set(\"matrix_name\",\"laplace_2d\");  An example of program using this
class is reported below.

Class CommandLineParser can be used as well. In this case, one may
decide to use the following:

Parameters:
-----------

In:  comm - Epetra communicator ";

%feature("docstring")
Trilinos_Util::CrsMatrixGallery::CrsMatrixGallery "Trilinos_Util::CrsMatrixGallery::CrsMatrixGallery(const string name,
const Epetra_Map &map)

Creates an Triutils_Gallery object using a given map.

Create a Triutils_Gallery object using an Epetra_Map. Problem size
must match the elements in map.

Parameters:
-----------

In:  name - definition of the problem to be created.

In:  map - Epetra_Map ";

%feature("docstring")
Trilinos_Util::CrsMatrixGallery::~CrsMatrixGallery "Trilinos_Util::CrsMatrixGallery::~CrsMatrixGallery()

Triutils_Gallery destructor. ";

%feature("docstring")  Trilinos_Util::CrsMatrixGallery::Set "int
Trilinos_Util::CrsMatrixGallery::Set(const string parameter, const int
value)

Sets a gallery options using an interger value. ";

%feature("docstring")  Trilinos_Util::CrsMatrixGallery::Set "int
Trilinos_Util::CrsMatrixGallery::Set(const string parameter, const
string value)

Sets a gallery options using a C++ string . ";

%feature("docstring")  Trilinos_Util::CrsMatrixGallery::Set "int
Trilinos_Util::CrsMatrixGallery::Set(const string parameter, const
double value)

Sets a gallery options using an double value. ";

%feature("docstring")  Trilinos_Util::CrsMatrixGallery::Set "int
Trilinos_Util::CrsMatrixGallery::Set(const string parameter, const
Epetra_Vector &value)

Sets a gallery options using an Epetra_Vector.

Sets a gallery options using an Epetra_Vector. The Epetra_Vector is
copied into internal structures, and freed by the destructor. ";

%feature("docstring")  Trilinos_Util::CrsMatrixGallery::Set "int
Trilinos_Util::CrsMatrixGallery::Set(Trilinos_Util::CommandLineParser
&CLP)

Sets gallery options using values passed from the shell. ";

%feature("docstring")  Trilinos_Util::CrsMatrixGallery::GetMatrix "Epetra_CrsMatrix * Trilinos_Util::CrsMatrixGallery::GetMatrix()

Returns a pointer to the CrsMatrix. ";

%feature("docstring")  Trilinos_Util::CrsMatrixGallery::GetMatrixRef "Epetra_CrsMatrix & Trilinos_Util::CrsMatrixGallery::GetMatrixRef() ";

%feature("docstring")
Trilinos_Util::CrsMatrixGallery::GetExactSolution "Epetra_MultiVector
* Trilinos_Util::CrsMatrixGallery::GetExactSolution()

Returns a pointer to the exact solution.

Returns a pointer to the exact solution.

Some choices are available to define the exact solution, using
Set(\"exact solution\", value). value can be: constant: the exact
solution vector is made up of 1's.

random: a random solution vector

linear: value at node i is defined as alpha*i. The double value alpha
can be set via Set(\"alpha\",DoubleVal). ";

%feature("docstring")
Trilinos_Util::CrsMatrixGallery::GetStartingSolution "Epetra_MultiVector *
Trilinos_Util::CrsMatrixGallery::GetStartingSolution()

Returns a pointer to the starting solution (typically, for HB
problems).

Returns a pointer to the starting solution. This is typically used
while reading a HB problem. However, the user can set a starting
solution using Set(\"starting solution\", \"value\"). Value can be
zero

random ";

%feature("docstring")  Trilinos_Util::CrsMatrixGallery::GetRHS "Epetra_MultiVector * Trilinos_Util::CrsMatrixGallery::GetRHS()

Returns a pointer to the rhs corresponding to the selected exact
solution. ";

%feature("docstring")  Trilinos_Util::CrsMatrixGallery::GetMap "const
Epetra_Map * Trilinos_Util::CrsMatrixGallery::GetMap()

Returns a pointer the internally stored Map. ";

%feature("docstring")  Trilinos_Util::CrsMatrixGallery::GetMapRef "const Epetra_Map & Trilinos_Util::CrsMatrixGallery::GetMapRef() ";

%feature("docstring")
Trilinos_Util::CrsMatrixGallery::GetLinearProblem "Epetra_LinearProblem *
Trilinos_Util::CrsMatrixGallery::GetLinearProblem()

Returns a pointer to Epetra_LinearProblem. ";

%feature("docstring")
Trilinos_Util::CrsMatrixGallery::ComputeResidual "void
Trilinos_Util::CrsMatrixGallery::ComputeResidual(double *residual)

Computes the 2-norm of the residual. ";

%feature("docstring")
Trilinos_Util::CrsMatrixGallery::ComputeDiffBetweenStartingAndExactSolutions
"void
Trilinos_Util::CrsMatrixGallery::ComputeDiffBetweenStartingAndExactSolutions(double
*residual)

Computes the 2-norm of the difference between the starting solution
and the exact solution. ";

%feature("docstring")
Trilinos_Util::CrsMatrixGallery::PrintMatrixAndVectors "void
Trilinos_Util::CrsMatrixGallery::PrintMatrixAndVectors(ostream &os)

Print out matrix and vectors. ";

%feature("docstring")
Trilinos_Util::CrsMatrixGallery::PrintMatrixAndVectors "void
Trilinos_Util::CrsMatrixGallery::PrintMatrixAndVectors() ";

%feature("docstring")
Trilinos_Util::CrsMatrixGallery::GetCartesianCoordinates "void
Trilinos_Util::CrsMatrixGallery::GetCartesianCoordinates(double *&x,
double *&y, double *&z)

Get pointers to double vectors containing coordinates of points. ";

%feature("docstring")  Trilinos_Util::CrsMatrixGallery::WriteMatrix "int Trilinos_Util::CrsMatrixGallery::WriteMatrix(const string
&FileName, const bool UseSparse=true)

Print matrix on file in MATLAB format. ";


// File: classTrilinos__Util_1_1InputFileReader.xml
%feature("docstring") Trilinos_Util::InputFileReader "";

%feature("docstring")  Trilinos_Util::InputFileReader::InputFileReader
"Trilinos_Util::InputFileReader::InputFileReader(const char
FileName[]) ";

%feature("docstring")
Trilinos_Util::InputFileReader::~InputFileReader "Trilinos_Util::InputFileReader::~InputFileReader() ";

%feature("docstring")  Trilinos_Util::InputFileReader::GetFileName "string Trilinos_Util::InputFileReader::GetFileName(void) const ";

%feature("docstring")  Trilinos_Util::InputFileReader::SetCommentChars
"void Trilinos_Util::InputFileReader::SetCommentChars(const string c)
";

%feature("docstring")
Trilinos_Util::InputFileReader::SetSeparationChars "void
Trilinos_Util::InputFileReader::SetSeparationChars(const string c) ";

%feature("docstring")  Trilinos_Util::InputFileReader::ReadFile "int
Trilinos_Util::InputFileReader::ReadFile() ";

%feature("docstring")  Trilinos_Util::InputFileReader::ReadFile "virtual int Trilinos_Util::InputFileReader::ReadFile(const char
FileName[]) ";


// File: structSPBLASMAT__STRUCT.xml
%feature("docstring") SPBLASMAT_STRUCT "";


// File: classTrilinos__Util__Map.xml
%feature("docstring") Trilinos_Util_Map "";

%feature("docstring")  Trilinos_Util_Map::Get "int
Trilinos_Util_Map::Get(const string input, const int def_value)

Gets the value of the specified option as an integer. If not found,
returns the specified default value. ";

%feature("docstring")  Trilinos_Util_Map::Get "double
Trilinos_Util_Map::Get(const string input, const double def_value)

Gets the value of the specified option as a double. If not found,
returns the specified default value. ";

%feature("docstring")  Trilinos_Util_Map::Get "string
Trilinos_Util_Map::Get(const string input, const string def_value)

Gets the value of the specified option as a string. If not found,
returns the specified default value. ";

%feature("docstring")  Trilinos_Util_Map::Set "bool
Trilinos_Util_Map::Set(const string input, const char *value)

Modify the value of a database entry.

This method modifies the value of a database entry. If the entry does
not exist in the database, return false. Otherwise, returns true. ";

%feature("docstring")  Trilinos_Util_Map::Set "bool
Trilinos_Util_Map::Set(const string input, const string value) ";

%feature("docstring")  Trilinos_Util_Map::Set "bool
Trilinos_Util_Map::Set(const string input, const int value) ";

%feature("docstring")  Trilinos_Util_Map::Set "bool
Trilinos_Util_Map::Set(const string input, const double value) ";

%feature("docstring")  Trilinos_Util_Map::Add "bool
Trilinos_Util_Map::Add(const string input, const string value)

Add an entry to the databse.

This method add an entry to the databse. First, it checks that this
entry does not exist. If it exists, the method returns false.
Otherwise, it adds the entry and returns true. ";

%feature("docstring")  Trilinos_Util_Map::SetLabel "bool
Trilinos_Util_Map::SetLabel(string Label) ";

%feature("docstring")  Trilinos_Util_Map::GetLabel "string
Trilinos_Util_Map::GetLabel(string Label) ";

%feature("docstring")  Trilinos_Util_Map::Has "bool
Trilinos_Util_Map::Has(const string input)

Check wheter an option is in the database or not.

This method checks whether option input is in the databse or not. It
returns true if it is, false otherwise. ";

%feature("docstring")  Trilinos_Util_Map::ShowAll "void
Trilinos_Util_Map::ShowAll() const

Show all the databse entries. ";

%feature("docstring")  Trilinos_Util_Map::ShowReallyAll "void
Trilinos_Util_Map::ShowReallyAll() const

Show all the databse entries, including entries beginning with \"_\".
";

%feature("docstring")  Trilinos_Util_Map::Reset "void
Trilinos_Util_Map::Reset(void) ";

%feature("docstring")  Trilinos_Util_Map::Trilinos_Util_Map "Trilinos_Util_Map::Trilinos_Util_Map() ";

%feature("docstring")  Trilinos_Util_Map::~Trilinos_Util_Map "virtual
Trilinos_Util_Map::~Trilinos_Util_Map() ";


// File: classTrilinos__Util_1_1VbrMatrixGallery.xml
%feature("docstring") Trilinos_Util::VbrMatrixGallery "";

%feature("docstring")
Trilinos_Util::VbrMatrixGallery::VbrMatrixGallery "Trilinos_Util::VbrMatrixGallery::VbrMatrixGallery(const string name,
const Epetra_Map &map) ";

%feature("docstring")
Trilinos_Util::VbrMatrixGallery::VbrMatrixGallery "Trilinos_Util::VbrMatrixGallery::VbrMatrixGallery(const string name,
const Epetra_Comm &Comm) ";

%feature("docstring")
Trilinos_Util::VbrMatrixGallery::~VbrMatrixGallery "Trilinos_Util::VbrMatrixGallery::~VbrMatrixGallery() ";

%feature("docstring")  Trilinos_Util::VbrMatrixGallery::GetBlockMap "const Epetra_BlockMap * Trilinos_Util::VbrMatrixGallery::GetBlockMap()

Returns a pointer the internally stored BlockMap. ";

%feature("docstring")  Trilinos_Util::VbrMatrixGallery::GetBlockMapRef
"const Epetra_BlockMap &
Trilinos_Util::VbrMatrixGallery::GetBlockMapRef() ";

%feature("docstring")  Trilinos_Util::VbrMatrixGallery::GetVbrMatrix "Epetra_VbrMatrix * Trilinos_Util::VbrMatrixGallery::GetVbrMatrix(const
int NumPDEEqns)

Returns a VbrMatrix, starting from the CsrMatrix.

Returns a VbrMatrix, starting from the CsrMatrix. This vbr matrix is
formally equivalent to the CrsMatrix returned by GetMatrix(). However,
each node of the CrsMatrix is replicated num_PDE_eqns times (this
value is passed in input, or set via Set(\"num pde eqns\",IntValue)).
";

%feature("docstring")  Trilinos_Util::VbrMatrixGallery::GetVbrMatrix "Epetra_VbrMatrix * Trilinos_Util::VbrMatrixGallery::GetVbrMatrix()

Returns a VbrMatrix, starting from the CsrMatrix. ";

%feature("docstring")
Trilinos_Util::VbrMatrixGallery::GetVbrMatrixRef "Epetra_VbrMatrix &
Trilinos_Util::VbrMatrixGallery::GetVbrMatrixRef() ";

%feature("docstring")  Trilinos_Util::VbrMatrixGallery::GetVbrRHS "Epetra_MultiVector * Trilinos_Util::VbrMatrixGallery::GetVbrRHS()

Returns a pointer to the RHS for the selected Vbr exact solution.

Returns a pointer to the RHS corresponding to the selected exact
solution to the linear systems defined by the Epetra_VbrMatrix. ";

%feature("docstring")
Trilinos_Util::VbrMatrixGallery::GetVbrExactSolution "Epetra_MultiVector *
Trilinos_Util::VbrMatrixGallery::GetVbrExactSolution()

Returns a pointer to the selected Vbr exact solution. ";

%feature("docstring")
Trilinos_Util::VbrMatrixGallery::GetVbrStartingSolution "Epetra_MultiVector *
Trilinos_Util::VbrMatrixGallery::GetVbrStartingSolution()

Returns a pointer to the starting solution for Vbr problems. ";

%feature("docstring")
Trilinos_Util::VbrMatrixGallery::CreateVbrMatrix "void
Trilinos_Util::VbrMatrixGallery::CreateVbrMatrix(void) ";

%feature("docstring")
Trilinos_Util::VbrMatrixGallery::GetVbrLinearProblem "Epetra_LinearProblem *
Trilinos_Util::VbrMatrixGallery::GetVbrLinearProblem()

Returns a pointer to Epetra_LinearProblem for VBR. ";

%feature("docstring")
Trilinos_Util::VbrMatrixGallery::ComputeResidualVbr "void
Trilinos_Util::VbrMatrixGallery::ComputeResidualVbr(double *residual)

Computes the 2-norm of the residual for the VBR problem. ";

%feature("docstring")
Trilinos_Util::VbrMatrixGallery::ComputeDiffBetweenStartingAndExactSolutionsVbr
"void
Trilinos_Util::VbrMatrixGallery::ComputeDiffBetweenStartingAndExactSolutionsVbr(double
*residual)

Computes the 2-norm of the difference between the starting solution
and the exact solution for the VBR problem. ";

%feature("docstring")
Trilinos_Util::VbrMatrixGallery::PrintVbrMatrixAndVectors "void
Trilinos_Util::VbrMatrixGallery::PrintVbrMatrixAndVectors(ostream &os)

Print out Vbr matrix and vectors. ";

%feature("docstring")
Trilinos_Util::VbrMatrixGallery::PrintVbrMatrixAndVectors "void
Trilinos_Util::VbrMatrixGallery::PrintVbrMatrixAndVectors() ";


// File: namespacestd.xml


// File: namespaceTrilinos__Util.xml


// File: Trilinos__Util_8h.xml
%feature("docstring")  Trilinos_Util_read_hb "void
Trilinos_Util_read_hb(char *data_file, int MyPID, int *N_global, int
*n_nonzeros, double **val, int **bindx, double **x, double **b, double
**xexact) ";

%feature("docstring")  Trilinos_Util_read_hb "void
Trilinos_Util_read_hb(char *data_file, int MyPID, int *N_global, int
*n_nonzeros, double **val, int **bindx) ";

%feature("docstring")  Trilinos_Util_read_coo "void
Trilinos_Util_read_coo(char *data_file, int MyPID, int *N_global, int
*n_nonzeros, double **val, int **bindx, double **x, double **b, double
**xexact) ";

%feature("docstring")  Trilinos_Util_ReadHb2Epetra "void
Trilinos_Util_ReadHb2Epetra(char *data_file, const Epetra_Comm &comm,
Epetra_Map *&map, Epetra_CrsMatrix *&A, Epetra_Vector *&x,
Epetra_Vector *&b, Epetra_Vector *&xexact) ";

%feature("docstring")  Trilinos_Util_ReadHpc2Epetra "void
Trilinos_Util_ReadHpc2Epetra(char *data_file, const Epetra_Comm &comm,
Epetra_Map *&map, Epetra_CrsMatrix *&A, Epetra_Vector *&x,
Epetra_Vector *&b, Epetra_Vector *&xexact) ";

%feature("docstring")  Trilinos_Util_ReadHb2EpetraVbr "void
Trilinos_Util_ReadHb2EpetraVbr(char *data_file, char *partitioning,
const Epetra_Comm &comm, Epetra_BlockMap *&map, Epetra_VbrMatrix *&A,
Epetra_Vector *&x, Epetra_Vector *&b, Epetra_Vector *&xexact) ";

%feature("docstring")  Trilinos_Util_distrib_msr_matrix "void
Trilinos_Util_distrib_msr_matrix(const Epetra_Comm &Comm, int
*N_global, int *n_nonzeros, int *N_update, int **update, double **val,
int **bindx, double **x, double **b, double **xexact) ";

%feature("docstring")  Trilinos_Util_distrib_msr_matrix "void
Trilinos_Util_distrib_msr_matrix(const Epetra_Comm &Comm, int
*N_global, int *n_nonzeros, int *N_update, int **update, double **val,
int **bindx) ";

%feature("docstring")  Trilinos_Util_distrib_vbr_matrix "void
Trilinos_Util_distrib_vbr_matrix(const Epetra_Comm &Comm, int
*N_global, int *N_blk_global, int *n_nonzeros, int *n_blk_nonzeros,
int *N_update, int **update, double **val, int **indx, int **rpntr,
int **cpntr, int **bpntr, int **bindx, double **x, double **b, double
**xexact) ";

%feature("docstring")  Trilinos_Util_create_vbr "void
Trilinos_Util_create_vbr(const Epetra_Comm &Comm, char *part_file, int
*N_global, int *N_blk_global, int *n_nonzeros, int *n_blk_nonzeros,
int *N_update, int **update, int *bindx_msr, double *val_msr, double
**val, int **indx, int **rpntr, int **cpntr, int **bpntr, int **bindx)
";

%feature("docstring")  Trilinos_Util_smsrres "double
Trilinos_Util_smsrres(int m, int n, double *val, int *indx, double
*xlocal, double *x, double *b) ";

%feature("docstring")  Trilinos_Util_scscres "double
Trilinos_Util_scscres(int isym, int m, int n, double *val, int *indx,
int *pntr, double *x, double *b) ";

%feature("docstring")  Trilinos_Util_scscmv "void
Trilinos_Util_scscmv(int isym, int m, int n, double *val, int *indx,
int *pntr, double *x, double *b) ";

%feature("docstring")  Trilinos_Util_svbrres "double
Trilinos_Util_svbrres(int m, int n, int m_blk, double *val, int *indx,
int *bindx, int *rpntr, int *cpntr, int *bpntrb, int *bpntre, double
*x, double *b) ";

%feature("docstring")  Trilinos_Util_msr2vbr "void
Trilinos_Util_msr2vbr(double val[], int indx[], int rnptr[], int
cnptr[], int bnptr[], int bindx[], int msr_bindx[], double msr_val[],
int total_blk_rows, int total_blk_cols, int blk_space, int nz_space,
int blk_type) ";

%feature("docstring")  Trilinos_Util_find_block_col "int
Trilinos_Util_find_block_col(int cnptr[], int column, int max_blocks,
int blk_size) ";

%feature("docstring")  Trilinos_Util_find_block_in_row "int
Trilinos_Util_find_block_in_row(int bindx[], int bnptr[], int blk_row,
int blk_col, int indx[], int no_elements, double val[], int blk_space,
int nz_space) ";

%feature("docstring")  Trilinos_Util_add_new_ele "void
Trilinos_Util_add_new_ele(int cnptr[], int col, int blk_row, int
bindx[], int bnptr[], int indx[], double val[], int row, double
new_ele, int maxcols, int blk_space, int nz_space, int blk_type) ";

%feature("docstring")  Trilinos_Util_find_closest_not_larger "int
Trilinos_Util_find_closest_not_larger(int key, int list[], int length)
";

%feature("docstring")  Trilinos_Util_convert_values_to_ptrs "void
Trilinos_Util_convert_values_to_ptrs(int array[], int length, int
start) ";

%feature("docstring")  Trilinos_Util_csrcsc "int
Trilinos_Util_csrcsc(int n, int n2, int job, int ipos, double *a, int
*ja, int *ia, double *ao, int *jao, int *iao) ";

%feature("docstring")  Trilinos_Util_csrmsr "int
Trilinos_Util_csrmsr(int n, double *a, int *ja, int *ia, double *ao,
int *jao, double *wk, int *iwk) ";

%feature("docstring")  Trilinos_Util_ssrcsr "int
Trilinos_Util_ssrcsr(int job, int value2, int nrow, double *a, int
*ja, int *ia, int nzmax, double *ao, int *jao, int *iao, int *indu,
int *iwk) ";

%feature("docstring")  Trilinos_Util_coocsr "int
Trilinos_Util_coocsr(int nrow, int nnz, double *a, int *ir, int *jc,
double *ao, int *jao, int *iao) ";

%feature("docstring")  Trilinos_Util_duscr_vbr "void
Trilinos_Util_duscr_vbr(int n, double *val, int *indx, int *bindx, int
*rpntr, int *cpntr, int *bpntrb, int *bpntre, SPBLASMAT *A) ";

%feature("docstring")  Trilinos_Util_dusmm "void
Trilinos_Util_dusmm(int m, int nrhs, int k, double alpha, SPBLASMAT
*A, double *x, int xstride, double beta, double *b, int bstride) ";

%feature("docstring")  Trilinos_Util_dusds_vbr "void
Trilinos_Util_dusds_vbr(SPBLASMAT *A) ";

%feature("docstring")  Trilinos_Util_GenerateCrsProblem "void
Trilinos_Util_GenerateCrsProblem(int nx, int ny, int npoints, int
*xoff, int *yoff, const Epetra_Comm &comm, Epetra_Map *&map,
Epetra_CrsMatrix *&A, Epetra_Vector *&x, Epetra_Vector *&b,
Epetra_Vector *&xexact, int indexBase=0) ";

%feature("docstring")  Trilinos_Util_GenerateCrsProblem "void
Trilinos_Util_GenerateCrsProblem(int nx, int ny, int npoints, int
*xoff, int *yoff, int nrhs, const Epetra_Comm &comm, Epetra_Map *&map,
Epetra_CrsMatrix *&A, Epetra_MultiVector *&x, Epetra_MultiVector *&b,
Epetra_MultiVector *&xexact, int indexBase=0) ";

%feature("docstring")  Trilinos_Util_GenerateVbrProblem "void
Trilinos_Util_GenerateVbrProblem(int nx, int ny, int npoints, int
*xoff, int *yoff, int nsizes, int *sizes, const Epetra_Comm &comm,
Epetra_BlockMap *&map, Epetra_VbrMatrix *&A, Epetra_Vector *&x,
Epetra_Vector *&b, Epetra_Vector *&xexact) ";

%feature("docstring")  Trilinos_Util_GenerateVbrProblem "void
Trilinos_Util_GenerateVbrProblem(int nx, int ny, int npoints, int
*xoff, int *yoff, int nsizes, int *sizes, int nrhs, const Epetra_Comm
&comm, Epetra_BlockMap *&map, Epetra_VbrMatrix *&A, Epetra_MultiVector
*&x, Epetra_MultiVector *&b, Epetra_MultiVector *&xexact) ";

%feature("docstring")  Trilinos_Util_ReadTriples2Epetra "int
Trilinos_Util_ReadTriples2Epetra(char *data_file, bool symmetric,
const Epetra_Comm &comm, Epetra_Map *&map, Epetra_CrsMatrix *&A,
Epetra_Vector *&x, Epetra_Vector *&b, Epetra_Vector *&xexact, bool
NonUniformMap=false, bool TimDavisHeader=false, bool ZeroBased=false)
";

%feature("docstring")  Trilinos_Util_ReadMatrixMarket2Epetra "int
Trilinos_Util_ReadMatrixMarket2Epetra(char *data_file, const
Epetra_Comm &comm, Epetra_Map *&map, Epetra_CrsMatrix *&A,
Epetra_Vector *&x, Epetra_Vector *&b, Epetra_Vector *&xexact) ";

%feature("docstring")  Trilinos_Util_write_vec "void
Trilinos_Util_write_vec(const char *filename, int n_equations, double
*x) ";

%feature("docstring")  Trilinos_Util_read_vec "void
Trilinos_Util_read_vec(const char *filename, int n_equations, double
*x) ";


// File: Trilinos__Util__CommandLineParser_8cpp.xml


// File: Trilinos__Util__CommandLineParser_8h.xml


// File: Trilinos__Util__coocsr_8cpp.xml
%feature("docstring")  Trilinos_Util_coocsr "int
Trilinos_Util_coocsr(int nrow, int nnz, double *a, int *ir, int *jc,
double *ao, int *jao, int *iao) ";


// File: Trilinos__Util__CountMatrixMarket_8cpp.xml
%feature("docstring")  Trilinos_Util_CountMatrixMarket "void
Trilinos_Util_CountMatrixMarket(const char *data_file, std::vector<
int > &non_zeros, int &N_rows, int &nnz, const Epetra_Comm &comm) ";


// File: Trilinos__Util__CountMatrixMarket_8h.xml
%feature("docstring")  Trilinos_Util_CountMatrixMarket "void
Trilinos_Util_CountMatrixMarket(const char *data_file, std::vector<
int > &non_zeros, int &N_rows, int &nnz, const Epetra_Comm &comm) ";


// File: Trilinos__Util__CountTriples_8cpp.xml
%feature("docstring")  Trilinos_Util_CountTriples "void
Trilinos_Util_CountTriples(const char *data_file, bool symmetric,
std::vector< int > &non_zeros, int &N_rows, int &nnz, const
Epetra_Comm &comm, bool TimDavisHeader=false, bool ZeroBased=false) ";


// File: Trilinos__Util__CountTriples_8h.xml
%feature("docstring")  Trilinos_Util_CountTriples "void
Trilinos_Util_CountTriples(const char *data_file, bool symmetric,
std::vector< int > &non_zeros, int &N_rows, int &nnz, const
Epetra_Comm &comm, bool TimDavisHeader=false, bool ZeroBased=false) ";


// File: Trilinos__Util__create__vbr_8cpp.xml
%feature("docstring")  Trilinos_Util_create_vbr "void
Trilinos_Util_create_vbr(const Epetra_Comm &Comm, char
*partition_file, int *N_global, int *N_blk_global, int *n_nonzeros,
int *n_blk_nonzeros, int *N_update, int **update, int *bindx_msr,
double *val_msr, double **val, int **indx, int **rpntr, int **cpntr,
int **bpntr, int **bindx) ";


// File: Trilinos__Util__CrsMatrixGallery_8cpp.xml


// File: Trilinos__Util__CrsMatrixGallery_8h.xml


// File: Trilinos__Util__csrcsc_8cpp.xml
%feature("docstring")  Trilinos_Util_csrcsc "int
Trilinos_Util_csrcsc(int n, int n2, int job, int ipos, double *a, int
*ja, int *ia, double *ao, int *jao, int *iao) ";


// File: Trilinos__Util__csrmsr_8cpp.xml
%feature("docstring")  Trilinos_Util_csrmsr "int
Trilinos_Util_csrmsr(int n, double *a, int *ja, int *ia, double *ao,
int *jao, double *wk, int *iwk) ";


// File: Trilinos__Util__distrib__msr__matrix_8cpp.xml
%feature("docstring")  Trilinos_Util_distrib_msr_matrix "void
Trilinos_Util_distrib_msr_matrix(const Epetra_Comm &Comm, int
*N_global, int *n_nonzeros, int *N_update, int **update, double **val,
int **bindx, double **x, double **b, double **xexact) ";

%feature("docstring")  Trilinos_Util_distrib_msr_matrix "void
Trilinos_Util_distrib_msr_matrix(const Epetra_Comm &Comm, int
*N_global, int *n_nonzeros, int *N_update, int **update, double **val,
int **bindx) ";


// File: Trilinos__Util__distrib__vbr__matrix_8cpp.xml
%feature("docstring")  Trilinos_Util_distrib_vbr_matrix "void
Trilinos_Util_distrib_vbr_matrix(const Epetra_Comm &Comm, int
*N_global, int *N_blk_global, int *n_nonzeros, int *n_blk_nonzeros,
int *N_update, int **update, double **val, int **indx, int **rpntr,
int **cpntr, int **bpntr, int **bindx, double **x, double **b, double
**xexact) ";


// File: Trilinos__Util__duscr__vbr_8cpp.xml
%feature("docstring")  Trilinos_Util_duscr_vbr "void
Trilinos_Util_duscr_vbr(int n, double *val, int *indx, int *bindx, int
*rpntr, int *cpntr, int *bpntrb, int *bpntre, SPBLASMAT *A) ";


// File: Trilinos__Util__dusds__vbr_8cpp.xml
%feature("docstring")  Trilinos_Util_dusds_vbr "void
Trilinos_Util_dusds_vbr(SPBLASMAT *A) ";


// File: Trilinos__Util__dusmm_8cpp.xml
%feature("docstring")  Trilinos_Util_dusmm "void
Trilinos_Util_dusmm(int m, int nrhs, int k, double alpha, SPBLASMAT
*A, double *x, int xstride, double beta, double *b, int bstride) ";


// File: Trilinos__Util__GenerateCrsProblem_8cpp.xml
%feature("docstring")  Trilinos_Util_GenerateCrsProblem "void
Trilinos_Util_GenerateCrsProblem(int nx, int ny, int npoints, int
*xoff, int *yoff, const Epetra_Comm &comm, Epetra_Map *&map,
Epetra_CrsMatrix *&A, Epetra_Vector *&x, Epetra_Vector *&b,
Epetra_Vector *&xexact, int indexBase) ";

%feature("docstring")  Trilinos_Util_GenerateCrsProblem "void
Trilinos_Util_GenerateCrsProblem(int nx, int ny, int npoints, int
*xoff, int *yoff, int nrhs, const Epetra_Comm &comm, Epetra_Map *&map,
Epetra_CrsMatrix *&A, Epetra_MultiVector *&x, Epetra_MultiVector *&b,
Epetra_MultiVector *&xexact, int indexBase) ";


// File: Trilinos__Util__GenerateVbrProblem_8cpp.xml
%feature("docstring")  Trilinos_Util_GenerateVbrProblem "void
Trilinos_Util_GenerateVbrProblem(int nx, int ny, int npoints, int
*xoff, int *yoff, int nsizes, int *sizes, const Epetra_Comm &comm,
Epetra_BlockMap *&map, Epetra_VbrMatrix *&A, Epetra_Vector *&x,
Epetra_Vector *&b, Epetra_Vector *&xexact) ";

%feature("docstring")  Trilinos_Util_GenerateVbrProblem "void
Trilinos_Util_GenerateVbrProblem(int nx, int ny, int npoints, int
*xoff, int *yoff, int nsizes, int *sizes, int nrhs, const Epetra_Comm
&comm, Epetra_BlockMap *&map, Epetra_VbrMatrix *&A, Epetra_MultiVector
*&x, Epetra_MultiVector *&b, Epetra_MultiVector *&xexact) ";


// File: Trilinos__Util__msr2vbr_8cpp.xml
%feature("docstring")  Trilinos_Util_msr2vbr "void
Trilinos_Util_msr2vbr(double val[], int indx[], int rnptr[], int
cnptr[], int bnptr[], int bindx[], int msr_bindx[], double msr_val[],
int total_blk_rows, int total_blk_cols, int blk_space, int nz_space,
int blk_type) ";

%feature("docstring")  Trilinos_Util_find_block_col "int
Trilinos_Util_find_block_col(int cnptr[], int column, int max_blocks,
int blk_size) ";

%feature("docstring")  Trilinos_Util_find_block_in_row "int
Trilinos_Util_find_block_in_row(int bindx[], int bnptr[], int blk_row,
int blk_col, int indx[], int no_elements, double val[], int blk_space,
int nz_space) ";

%feature("docstring")  Trilinos_Util_add_new_ele "void
Trilinos_Util_add_new_ele(int cnptr[], int col, int blk_row, int
bindx[], int bnptr[], int indx[], double val[], int row, double
new_ele, int maxcols, int blk_space, int nz_space, int blk_type) ";

%feature("docstring")  Trilinos_Util_find_closest_not_larger "int
Trilinos_Util_find_closest_not_larger(int key, int list[], int length)
";

%feature("docstring")  Trilinos_Util_convert_values_to_ptrs "void
Trilinos_Util_convert_values_to_ptrs(int array[], int length, int
start) ";


// File: Trilinos__Util__read__coo_8cpp.xml
%feature("docstring")  Trilinos_Util_read_coo "void
Trilinos_Util_read_coo(char *data_file, int MyPID, int *N_global, int
*n_nonzeros, double **val, int **bindx, double **x, double **b, double
**xexact) ";


// File: Trilinos__Util__read__hb_8cpp.xml
%feature("docstring")  Trilinos_Util_read_hb "void
Trilinos_Util_read_hb(char *data_file, int MyPID, int *N_global, int
*n_nonzeros, double **val, int **bindx, double **x, double **b, double
**xexact) ";

%feature("docstring")  Trilinos_Util_read_hb "void
Trilinos_Util_read_hb(char *data_file, int MyPID, int *N_global, int
*n_nonzeros, double **val, int **bindx) ";


// File: Trilinos__Util__read__vec_8cpp.xml
%feature("docstring")  Trilinos_Util_read_vec "void
Trilinos_Util_read_vec(const char *filename, int n_equations, double
*x) ";


// File: Trilinos__Util__ReadHb2Epetra_8cpp.xml
%feature("docstring")  Trilinos_Util_ReadHb2Epetra "void
Trilinos_Util_ReadHb2Epetra(char *data_file, const Epetra_Comm &comm,
Epetra_Map *&map, Epetra_CrsMatrix *&A, Epetra_Vector *&x,
Epetra_Vector *&b, Epetra_Vector *&xexact) ";


// File: Trilinos__Util__ReadHb2EpetraVbr_8cpp.xml
%feature("docstring")  Trilinos_Util_ReadHb2EpetraVbr "void
Trilinos_Util_ReadHb2EpetraVbr(char *data_file, char *partitioning,
const Epetra_Comm &comm, Epetra_BlockMap *&map, Epetra_VbrMatrix *&A,
Epetra_Vector *&x, Epetra_Vector *&b, Epetra_Vector *&xexact) ";


// File: Trilinos__Util__ReadHpc2Epetra_8cpp.xml
%feature("docstring")  Trilinos_Util_ReadHpc2Epetra "void
Trilinos_Util_ReadHpc2Epetra(char *data_file, const Epetra_Comm &comm,
Epetra_Map *&map, Epetra_CrsMatrix *&A, Epetra_Vector *&x,
Epetra_Vector *&b, Epetra_Vector *&xexact) ";


// File: Trilinos__Util__ReadMatrixMarket2Epetra_8cpp.xml
%feature("docstring")  Trilinos_Util_ReadMatrixMarket2Epetra "int
Trilinos_Util_ReadMatrixMarket2Epetra(char *data_file, const
Epetra_Comm &comm, Epetra_Map *&map, Epetra_CrsMatrix *&A,
Epetra_Vector *&x, Epetra_Vector *&b, Epetra_Vector *&xexact) ";


// File: Trilinos__Util__ReadMatrixMarket2Epetra_8h.xml
%feature("docstring")  Trilinos_Util_ReadMatrixMarket2Epetra "int
Trilinos_Util_ReadMatrixMarket2Epetra(char *data_file, const
Epetra_Comm &comm, Epetra_Map *&map, Epetra_CrsMatrix *&A,
Epetra_Vector *&x, Epetra_Vector *&b, Epetra_Vector *&xexact) ";


// File: Trilinos__Util__ReadTriples2Epetra_8cpp.xml
%feature("docstring")  Trilinos_Util_ReadTriples2Epetra "int
Trilinos_Util_ReadTriples2Epetra(char *data_file, bool symmetric,
const Epetra_Comm &comm, Epetra_Map *&map, Epetra_CrsMatrix *&A,
Epetra_Vector *&x, Epetra_Vector *&b, Epetra_Vector *&xexact, bool
NonUniformMap=false, bool TimDavisHeader=false, bool ZeroBased=false)
";


// File: Trilinos__Util__scscmv_8cpp.xml
%feature("docstring")  Trilinos_Util_scscmv "void
Trilinos_Util_scscmv(int isym, int m, int n, double *val, int *indx,
int *pntr, double *x, double *y) ";


// File: Trilinos__Util__scscres_8cpp.xml
%feature("docstring")  Trilinos_Util_scscres "double
Trilinos_Util_scscres(int isym, int m, int n, double *val, int *indx,
int *pntr, double *x, double *b) ";


// File: Trilinos__Util__smsrres_8cpp.xml
%feature("docstring")  Trilinos_Util_smsrres "double
Trilinos_Util_smsrres(int m, int n, double *val, int *indx, double
*xlocal, double *x, double *b) ";


// File: Trilinos__Util__ssrcsr_8cpp.xml
%feature("docstring")  Trilinos_Util_ssrcsr "int
Trilinos_Util_ssrcsr(int job, int value2, int nrow, double *a, int
*ja, int *ia, int nzmax, double *ao, int *jao, int *iao, int *indu,
int *iwk) ";


// File: Trilinos__Util__svbrres_8cpp.xml
%feature("docstring")  Trilinos_Util_svbrres "double
Trilinos_Util_svbrres(int m, int n, int m_blk, double *val, int *indx,
int *bindx, int *rpntr, int *cpntr, int *bpntrb, int *bpntre, double
*x, double *b) ";


// File: Trilinos__Util__Version_8h.xml
%feature("docstring")  Triutils_Version "string Triutils_Version() ";


// File: Trilinos__Util__write__vec_8cpp.xml
%feature("docstring")  Trilinos_Util_write_vec "void
Trilinos_Util_write_vec(const char *filename, int n_equations, double
*x) ";


// File: dir_a26fd566400729558a414fa6e54cace3.xml


// File: dir_06d92e7e2071799722463955d21cb037.xml


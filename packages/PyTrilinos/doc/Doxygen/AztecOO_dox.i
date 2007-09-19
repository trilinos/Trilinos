
// File: index.xml

// File: classAztecOO.xml
%feature("docstring") AztecOO "

AztecOO: An object-oriented wrapper for Aztec.

Currently it accepts a Petra matrix, initial guess and RHS as separate
arguments, or alternatively, accepts a Epetra_LinearProblem. If
constructed using a Epetra_LinearProblem, AztecOO will infer some
solver/preconditioner, etc., options and parameters. Users may
override these choices and manually choose from among the full set of
Aztec options using the SetAztecOption() and SetAztecParam()
functions.

AztecOO will solve a linear systems of equations: $ AX=B $, using
Epetra objects and the Aztec solver library, where $A$ is an
Epetra_Operator or Epetra_RowMatrix (note that the Epetra_Operator
class is a base class for Epetra_RowMatrix so that Epetra_RowMatrix
isa Epetra_Operator.) $X$ and $B$ are Epetra_MultiVector objects.

WARNING:   AztecOO does not presently support solution of more than
one simultaneous right-hand-side.

C++ includes: AztecOO.h ";

%feature("docstring")  AztecOO::SetPreconditioner "int
AztecOO::SetPreconditioner(AZ_PRECOND *Prec)

AztecOO External Preconditioner Set (object).

Associates an already defined Aztec preconditioner with this solve. ";

%feature("docstring")  AztecOO::SetPreconditioner "int
AztecOO::SetPreconditioner(AZ_PREC_FUN prec_function, void *prec_data)

AztecOO External Preconditioner Set (function and data).

Associates an external function and data pointer with preconditioner
";

%feature("docstring")  AztecOO::SetScaling "int
AztecOO::SetScaling(struct AZ_SCALING *Scaling)

AztecOO External Scaling Set.

Associates an already defined Aztec scaling object with this solve. ";

%feature("docstring")  AztecOO::SetMatrixName "int
AztecOO::SetMatrixName(int label)

AztecOO Label Matrix for Aztec.

This is used to label individual matrices within Aztec. This might be
useful if several Aztec invocations are involved corresponding to
different matrices. ";


// File: structAztecOO_1_1MatrixData.xml
%feature("docstring") AztecOO::MatrixData "";

%feature("docstring")  AztecOO::MatrixData::MatrixData "AztecOO::MatrixData::MatrixData(Epetra_RowMatrix *inA=0, Epetra_Vector
*inX=0, Epetra_Vector *inY=0, Epetra_Vector *inSourceVec=0,
Epetra_Vector *inTargetVec=0) ";

%feature("docstring")  AztecOO::MatrixData::~MatrixData "AztecOO::MatrixData::~MatrixData() ";


// File: structAztecOO_1_1OperatorData.xml
%feature("docstring") AztecOO::OperatorData "";

%feature("docstring")  AztecOO::OperatorData::OperatorData "AztecOO::OperatorData::OperatorData(Epetra_Operator *inA=0,
Epetra_Vector *inX=0, Epetra_Vector *inY=0) ";

%feature("docstring")  AztecOO::OperatorData::~OperatorData "AztecOO::OperatorData::~OperatorData() ";


// File: classAztecOO__Operator.xml
%feature("docstring") AztecOO_Operator "

AztecOO_Operator: An implementation of the Epetra_Operator class.

The AztecOO_Operator class implements Epetra_Operator using a pre-
constructed AztecOO solver object. Once constructed, an
AztecOO_Operator can be used as a preconditioner within another
AztecOO solver object.

C++ includes: AztecOO_Operator.h ";


// File: classAztecOO__StatusTest.xml
%feature("docstring") AztecOO_StatusTest "

AztecOO_StatusTest: A pure virtual class for extending the status
testing capabilities of AztecOO.

C++ includes: AztecOO_StatusTest.h ";

%feature("docstring")  AztecOO_StatusTest::PrintStatus "virtual void
AztecOO_StatusTest::PrintStatus(ostream &os, AztecOO_StatusType type)
const ";


// File: classAztecOO__StatusTestCombo.xml
%feature("docstring") AztecOO_StatusTestCombo "

AztecOO_StatusTestCombo: A class for extending the status testing
capabilities of AztecOO via logical combinations.

AztecOO_StatusTestCombo is an interface that can be implemented to
extend the convergence testing capabilities of AztecOO. This class
supports composite tests. In this situation, two or more existing
AztecOO_StatusTestCombo objects test1 and test2 can be used to create
a new test. For all combinations, if any tests returns Failed or
returns not-a-number (NaN) status, then the combination test returns
Failed. There are three possible combinations: OR combination: If an
OR combination is selected, the status returns Converged if any one of
the subtest returns as Converged.

AND combination: If an AND combination is selected, the status returns
Converged only when all subtests return as Converged.

SEQ combination: SEQ is a form of AND that will perform subtests in
sequence. If the first test returns Unconverged, Failed or NaN, no
other subtests are done, and the status is returned as Unconverged if
the first test was Unconverged, or as Failed if the first test was
Failed or NaN. If the first test returns Converged, the second test is
checked in the same fashion as the first. If the second test is
Converged, the third one is tested, and so on.

The purpose of the SEQ combination is to allow the addition of
expensive but more rigorous convergence tests. For example, we could
define a test that used the implicit residual vector (the one produced
by the iterative method) as the first subtest and define a second test
using the explicitly computed residual vector. Explicitly computing
the residual requires a matrix multiplication with the original matrix
operator, an expensive operation. By using the SEQ combination, we can
avoid the matrix multiplication associated with the explicit residual
calculation until the implicit residual is small.

WARNING:  Presently it is not valid to associate one status test
instance with two different AztecOO objects.

C++ includes: AztecOO_StatusTestCombo.h ";


// File: classAztecOO__StatusTestMaxIters.xml
%feature("docstring") AztecOO_StatusTestMaxIters "

AztecOO_StatusTestMaxIters: An AztecOO_StatusTest class specifying a
maximum number of iterations.

C++ includes: AztecOO_StatusTestMaxIters.h ";


// File: classAztecOO__StatusTestResNorm.xml
%feature("docstring") AztecOO_StatusTestResNorm "

AztecOO_StatusTestResNorm: An implementation of AztecOO_StatusTest
using a family of residual norms.

AztecOO_StatusTestResNorm is an implementation of AztecOO_StatusTest
that allows a user to construct one of a family of residual tests for
use as a status/convergence test for AztecOO. The form of the test is
\\\\[ \\\\frac{\\\\|r\\\\|}{\\\\sigma} \\\\le \\\\tau \\\\] where  $r$
is the residual vector, implicitly or explicitly computed (determined
by enum ResType),

$\\\\|r\\\\|$ is the residual norm determined by the enum NormType
(1-norm, 2-norm or inf-norm),

$\\\\sigma$ is the scale factor that can be passed in as a precomputed
double precision number, or can be selected from by the enum ScaleType
(norm of RHS, norm of initial residual).

$\\\\tau$ is the tolerance that is passed in as a double precision
number to the constructor. The value of $\\\\tau$ can be reset using
the ResetTolerance() method.

WARNING:  Presently it is not valid to associate one status test
instance with two different AztecOO objects.

C++ includes: AztecOO_StatusTestResNorm.h ";


// File: classAztecOOConditionNumber.xml
%feature("docstring") AztecOOConditionNumber "

Condition number estimator using AztecOO.

This object will estimate the condition number of an Epetra_Operator.

C++ includes: AztecOO_ConditionNumber.h ";

%feature("docstring")  AztecOOConditionNumber::AztecOOConditionNumber
"AztecOOConditionNumber::AztecOOConditionNumber()

Constructor. ";

%feature("docstring")  AztecOOConditionNumber::~AztecOOConditionNumber
"AztecOOConditionNumber::~AztecOOConditionNumber()

Destructor. ";

%feature("docstring")  AztecOOConditionNumber::initialize "void
AztecOOConditionNumber::initialize(const Epetra_Operator &op,
SolverType solverType=GMRES_, int krylovSubspaceSize=100, bool
printSolve=false)

Initialization. ";

%feature("docstring")  AztecOOConditionNumber::computeConditionNumber
"int AztecOOConditionNumber::computeConditionNumber(int maxIters,
double tol)

Estimates the condition number. ";

%feature("docstring")  AztecOOConditionNumber::getConditionNumber "double AztecOOConditionNumber::getConditionNumber()

Return condition number computed by last call to
computeConditionNumber. ";


// File: AztecOO_8cpp.xml
%feature("docstring")  AztecOO_uppercase "string
AztecOO_uppercase(const string &s) ";

%feature("docstring")  Epetra_Aztec_matnorminf "double
Epetra_Aztec_matnorminf(AZ_MATRIX *Amat) ";

%feature("docstring")  Epetra_Aztec_matvec "void
Epetra_Aztec_matvec(double x[], double y[], AZ_MATRIX *Amat, int
proc_config[]) ";

%feature("docstring")  Epetra_Aztec_operatornorminf "double
Epetra_Aztec_operatornorminf(AZ_MATRIX *Amat) ";

%feature("docstring")  Epetra_Aztec_operatorvec "void
Epetra_Aztec_operatorvec(double x[], double y[], AZ_MATRIX *Amat, int
proc_config[]) ";

%feature("docstring")  Epetra_Aztec_precond "void
Epetra_Aztec_precond(double x[], int input_options[], int
proc_config[], double input_params[], AZ_MATRIX *Amat, AZ_PRECOND
*prec) ";

%feature("docstring")  Epetra_Aztec_getrow "int
Epetra_Aztec_getrow(int columns[], double values[], int row_lengths[],
AZ_MATRIX *Amat, int N_requested_rows, int requested_rows[], int
allocated_space) ";

%feature("docstring")  Epetra_Aztec_comm_wrapper "int
Epetra_Aztec_comm_wrapper(double vec[], AZ_MATRIX *Amat) ";

%feature("docstring")  AztecOO_StatusTest_wrapper "void
AztecOO_StatusTest_wrapper(void *conv_test_obj, void *res_vector_obj,
int iteration, double *res_vector, int print_info, int sol_updated,
int *converged, int *isnan, double *rnorm, int *r_avail) ";


// File: AztecOO_8h.xml
%feature("docstring")  Epetra_Aztec_matvec "void
Epetra_Aztec_matvec(double x[], double y[], AZ_MATRIX *Amat, int
proc_config[]) ";

%feature("docstring")  Epetra_Aztec_matnorminf "double
Epetra_Aztec_matnorminf(AZ_MATRIX *Amat) ";

%feature("docstring")  Epetra_Aztec_operatorvec "void
Epetra_Aztec_operatorvec(double x[], double y[], AZ_MATRIX *Amat, int
proc_config[]) ";

%feature("docstring")  Epetra_Aztec_operatornorminf "double
Epetra_Aztec_operatornorminf(AZ_MATRIX *Amat) ";

%feature("docstring")  Epetra_Aztec_precond "void
Epetra_Aztec_precond(double x[], int input_options[], int
proc_config[], double input_params[], AZ_MATRIX *Amat, AZ_PRECOND
*prec) ";

%feature("docstring")  Epetra_Aztec_getrow "int
Epetra_Aztec_getrow(int columns[], double values[], int row_lengths[],
AZ_MATRIX *Amat, int N_requested_rows, int requested_rows[], int
allocated_space) ";

%feature("docstring")  Epetra_Aztec_comm_wrapper "int
Epetra_Aztec_comm_wrapper(double vec[], AZ_MATRIX *Amat) ";

%feature("docstring")  AztecOO_StatusTest_wrapper "void
AztecOO_StatusTest_wrapper(void *conv_test_obj, void *res_vector_obj,
int iteration, double *res_vector, int print_info, int sol_updated,
int *converged, int *isnan, double *rnorm, int *r_avail) ";


// File: AztecOO__ConditionNumber_8cpp.xml


// File: AztecOO__ConditionNumber_8h.xml


// File: AztecOO__ConfigDefs_8h.xml


// File: AztecOO__Operator_8cpp.xml


// File: AztecOO__Operator_8h.xml


// File: AztecOO__Scaling_8cpp.xml
%feature("docstring")  AZOO_Scale "int AZOO_Scale(int action,
Epetra_RowMatrix *A, double b[], double x[], int options[], AZ_SCALING
*scaling) ";

%feature("docstring")  AZOO_create_scaling_vector "Epetra_Vector *
AZOO_create_scaling_vector(Epetra_RowMatrix *A, int scaling_type) ";

%feature("docstring")  AztecOO_scale_epetra "int
AztecOO_scale_epetra(int action, AZ_MATRIX *Amat, int options[],
double b[], double x[], int proc_config[], AZ_SCALING *scaling)

A scaling function that can be assigned to the AZ_SCALING.scale
function-pointer, to provide a call-back that Aztec can use to scale
Epetra matrices passed in by AztecOO. ";


// File: AztecOO__Scaling_8h.xml
%feature("docstring")  AztecOO_scale_epetra "int
AztecOO_scale_epetra(int action, AZ_MATRIX *Amat, int options[],
double b[], double x[], int proc_config[], AZ_SCALING *scaling)

A scaling function that can be assigned to the AZ_SCALING.scale
function-pointer, to provide a call-back that Aztec can use to scale
Epetra matrices passed in by AztecOO. ";


// File: AztecOO__StatusTest_8h.xml


// File: AztecOO__StatusTestCombo_8cpp.xml


// File: AztecOO__StatusTestCombo_8h.xml


// File: AztecOO__StatusTestMaxIters_8cpp.xml


// File: AztecOO__StatusTestMaxIters_8h.xml


// File: AztecOO__StatusTestResNorm_8cpp.xml


// File: AztecOO__StatusTestResNorm_8h.xml


// File: AztecOO__StatusType_8h.xml


// File: AztecOO__string__maps_8cpp.xml


// File: AztecOO__string__maps_8h.xml


// File: AztecOO__Version_8h.xml
%feature("docstring")  AztecOO_Version "string AztecOO_Version() ";


// File: dir_4dc5e7d6705411ebb2207ed00e1644b4.xml


// File: dir_c12959d073c8d91317941666aec70eef.xml


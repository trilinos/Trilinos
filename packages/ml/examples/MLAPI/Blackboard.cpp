/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include "ml_config.h"
#include "ml_common.h"
#ifdef HAVE_ML_MLAPI

// This is the easiest way, since all MLAPI include files will 
// automatically be included. However, this also makes the compilation
// longer, so you might want to specify the exact list of include files
// as long as you code gets finalized.
#include "MLAPI.h"

using namespace Teuchos;
using namespace MLAPI;

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  try {

    // Initialize the MLAPI workspace
    Init();

    // all MLAPI objects are based on Space's. A Space defines the number of
    // local and global distribution of elements. Space's can be constructed
    // in several ways. The simplest is to specify the number of global
    // elements:
    
    Space MySpace(4 * GetNumProcs());

    // Class MLAPI::DistributedMatrix is a simple way to define
    // matrices. This class can be used to define distributed matrices. Each
    // processor can set on-processor and off-processor elements. Before any
    // use, the matrix need to be "freezed", so that the communication pattern
    // can be established. This is done by method FillComplete(). After having
    // called this method, elements can no longer be added; however, already
    // inserted elements can be modified using method ReplaceElement().
    //
    // This class is based on the Epetra_FECrsMatrix class. The insertion of
    // elements is somehow slower, but the matrix-vector product is as
    // efficient as those of the Epetra_FECrsmatrix class (with the only
    // overhead of an additional function call).
    
    DistributedMatrix A(MySpace, MySpace);

    // As any processor can set any element, here we fill the entire
    // matrix on processor 0 only. It is of course possible (and preferable)
    // to let each processor fill local rows only.

    if (GetMyPID() == 0) {

      for (int i = 0 ; i < MySpace.GetNumGlobalElements() ; ++i) {

        if (i) 
          A(i, i - 1) = -1.0;
        if (i + 1 != A.NumGlobalCols())
          A(i, i + 1) = -1.0;
        A(i, i) = 2.0;
      }
    }

    A.FillComplete(); 

    // To get the (row, col) value of the matrix, use method
    // value = GetElement(row, col). Note that both `row' and `col'
    // refer to global indices; however row must be a locally hosted row.
    // If (row, col) is not found, value is set to 0.0.

    cout << A;

    // Here we define 3 vectors. The last, z, is empty.
    
    MultiVector x(MySpace);
    MultiVector y(MySpace);
    MultiVector z;

    // It is easy to assign all the elements of x to a constant value,
    // or to populate with random numbers.
    
    x = 1.0;
    y.Random();

    // Work on vector's elements requires the () operator

    for (int i = 0 ; i < y.GetMyLength() ; ++i)
      y(i) = 1.0 * i + x(i);

    // vectors can be manipulated in an intuitive way:
    
    z = x + y;                        // addition
    double norm = sqrt(z * z);        // dot product
    double Anorm = sqrt(z * (A * z)); // dot product with a matrix (not the parenthesis)
    x = x / (x * y);                  // scale x by the dot product between x and y

    if (GetMyPID() == 0) {
      cout << "2-norm of z = " << norm << endl;
      cout << "A-norm of z = " << Anorm << endl;
    }

    // some basic operations on matrices are also supported. For example,
    // one can extract the diagonal of a matrix as a vector, then create a new
    // matrix, containing this vector on the diagonal
    // (that is, B will be a diagonal matrix so that B_{i,i} = A_{i,i})
    
    Operator B = GetDiagonal(GetDiagonal(A));

    // verify that the diagonal of A - B is zero
    
    z = GetDiagonal(B - A);

    // Also operators can be composed intuitively (for compatible
    // operators):

    Operator C = A * B;    // multiplication, A 
    C = A + B;             // addition
    C = 1.0 * A - B * 2.0;             // subtraction

    // use Amesos to apply the inverse of A using LU factorization
    InverseOperator invA(A,"Amesos-KLU");

    // verify that x == inv(A) * A * x
    x = invA * (A * x) - x;
    double NormX = sqrt(x * x);

    if (GetMyPID() == 0)
      cout << "Norm of inv(A) * A * x - x = " << NormX << endl;

    // All MLAPI objects have a Label, which can be set using
    // SetLabel(Label). Also, all objects overload the << operator:

    cout << MySpace;
    cout << x;
    cout << A;
      
    // Function Eig() can be used to compute the eigenvalues of an Operator.
    // This function calls LAPACK, therefore the Operator should be "small".
    // ER will contain the real part of the eigenvalues;
    // EI will contain the imaginary part of the eigenvalues;
    
    MultiVector ER, EI;
    Eig(A, ER, EI);

    for (int i = 0 ; i < ER.GetMyLength() ; ++i)
      cout << "ER(" << MySpace(i) << ") = " << ER(i) << endl;

    // Another nice feature is that objects can be printed in a MATLAB format.
    // To that aim, simply do the following:
    // - set the label of the objects you want to port to MATLAB;
    // - create a MATLABStream object;
    // - use the << operator on MATLABStream
    // - File `Blackboard.m' has been created in your directory. This file
    //   (only one, for serial and parallel runs) can be executed in
    //   MATLAB to reproduce your MLAPI variables.

    MATLABStream matlab("Blackboard.m", false);

    matlab << "% you can add any MATLAB command this way\n";

    x.SetLabel("x");
    matlab << x;

    A.SetLabel("A");
    matlab << A;

    ER.SetLabel("ER");
    EI.SetLabel("EI");

    matlab << ER;
    matlab << EI;
    matlab << "plot(ER, EI, 'o')\n";

    // Finalize the MLAPI work space before leaving the application
    
    Finalize();

  } 
  catch (const int e) {
    cout << "Integer exception, code = " << e << endl;
  } 
  catch (...) {
    cout << "problems here..." << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);

}

#else

#include "ml_include.h"

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("This MLAPI example requires the following configuration options:");
  puts("\t--enable-epetra");
  puts("\t--enable-teuchos");
  puts("\t--enable-ifpack");
  puts("\t--enable-amesos");
  puts("Please check your configure line.");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);
}

#endif // #if defined(HAVE_ML_MLAPI)

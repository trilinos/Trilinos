//@HEADER
// ************************************************************************
// 
//               ML: A Multilevel Preconditioner Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "ml_config.h"
#include "ml_common.h"
#ifdef HAVE_ML_MLAPI

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
    
    Space MySpace(10);

    // MLAPI::Matrix is a very simple and convenient Epetra_RowMatrix derived
    // class. Inserting a new element is just A(row, col) = val.
    // Most of the methods of Epetra_RowMatrix are implemented in MLAPI::Matrix.
    //
    // MLAPI::Matrix can be used for serial computations only. Furthermore, it
    // is *not* meant to be efficient, just easy-to-use. Users should consider
    // other Epetra_RowMatrix derived classes (like Epetra_CrsMatrix or
    // Epetra_VbrMatrix) in order to define parallel and scalable matrices.
    // Function Gallery() returns some parallel matrices, see one of the
    // other MLAPI examples.
    //
    // NOTE: each time A(row,col) is called, a zero element is inserted if
    // not already present in the matrix!
    
    Matrix A_Mat(MySpace, MySpace);
    for (int i = 0 ; i < 10 ; ++i) {
      if (i) A_Mat(i, i - 1) = -1.0;
      A_Mat(i,i) = 2.0;
      if (i + 1 != A_Mat.NumGlobalCols())
        A_Mat(i, i + 1) = -1.0;
    }

    // Note that MLAPI::Matrix cannot be copied or reassigned, and no
    // operators are overloaded on this class. The only way to use an
    // MLAPI::Matrix with other MLAPI objects is to wrap it into an
    // MLAPI::Operator, as done in the following line.
    //
    // NOTE: The last parameter has be to `false' because the Operator 
    // should not delete the A_Mat object.

    Operator A(MySpace, MySpace, &A_Mat, false);

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
    C = A - B;             // subtraction

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
      
    // Function Eig() can be used to compute the eigenvalues of an Operator
    // (for serial runs only). This function calls LAPACK, therefore the
    // Operator should be "small".
    // ER will contain the real part of the eigenvalues;
    // EI will contain the imaginary part of the eigenvalues;
    // V will contain the eigenvalues.
    
    MultiVector ER, EI, V;
    Eig(A, ER, EI, V);

    for (int i = 0 ; i < ER.GetMyLength() ; ++i)
      for (int j = 0 ; j < ER.GetNumVectors() ; ++j)
        cout << "ER(" << i << ", " << j << ") = " << ER(i,j) << endl;
    
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
    V.SetLabel("V");

    matlab << ER;
    matlab << EI;
    matlab << V;
    matlab << "plot(ER, EI, 'o')\n";

    // Finalize the MLAPI work space before leaving the application
    
    Finalize();

  } 
  catch (exception& e) {
    cout << e.what() << endl;
  } 
  catch (int e) {
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

  puts("The ML API requires the following configuration options:");
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

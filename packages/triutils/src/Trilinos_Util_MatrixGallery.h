// @HEADER
// ***********************************************************************
// 
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2001) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#ifndef __TRILINOS_UTILS_GALLERY_H
#define __TRILINOS_UTILS_GALLERY_H

class Epetra_Comm;
class Epetra_Map;
class Epetra_BlockMap;
class Vector;
class Epetra_CrsMatrix;
class Epetra_VbrMatrix;
class Epetra_Export;
class Epetra_LinearProblem;
#include <string>
#include "Trilinos_Util_ShellOptions.h"

class Trilinos_Util_CrsMatrixGallery 
{
public:

  //@{ \name Constructors/Destructor.
  //! Triutils_Gallery Constructor.
  /*! Creates a Triutils_Gallery instance.

  \param In
  name - identity the problem to be generated. At this time, one can

  - \c eye: Creates an identity matrix. The dimension of the problem is set
  using Set("problem size", IntValue), or, alternatively, by Set("nx",
  IntValue).

  - \c diag: Creates a diagonal matrix. The elements on the diagonal can
  set using Set("a",value), where value is either a double, or an
  Epetra_Vector. Problem size set as for "eye."

  - \c tridiag: Creates a tridiagonal matrix. The diagonal element is set
  using Set("a", DiagValue), the subdiagonal using
  Set("b",SubDiagValue), and the superdiagonal using
  Set("c",SupDiagValue).  DiagValue, SubDiagValue, and SupDiagValue can
  be scalar or Epetra vectors.  In the latter case, all vectors must
  have the same number of elements. Elements i will be use to define row
  i. Problem size specified as for "eye."

  - \c laplace_1d: Create the classical tridiagonal matrix with stencil
  [-1, 2, -1].  Problem size specified as for "eye."

  - \c laplace_2d: Creates a matrix corresponding to the stencil of a 2D
  Laplacian operator on a structured cartesian grid. Problem size is
  specified using Set("problem size", IntValue). In this case, IntValue
  must be a square number. Alternatively, one can set the number of
  nodes along the x-axis and y-axis, using Set("nx",IntValue) and
  Set("ny",IntValue).

  - \c cross_stencil_2d: Creates a matrix with the same stencil of "laplace
  2d", but with arbitrary values. The stencil is
  \[
  A = \begin{tabular}{ccc}
    & e &  \\
  b & a & c \\
    & d &   \\
  \end{tabular}
  \]
  and each value can be either a double, or and Epetra_Vector. Problem
  size specified as in "laplace 2d."

  - \c laplace_3d: Creates a matrix corresponding to the stencil of a 3D
  Laplacian operator on a structured cartesian grid. Problem size
  specified using Set("problme size",IntValue). In this case, IntValue
  must be a cube. Alternatively, one can specify the number of nodes
  along the axis, using Set("nx",IntValue), Set("ny",IntValue), and
  Set("nz",IntValue).

  - \c cross_stencil_3d: Similar to the 2D case. The matrix stencil
  correspond to that of a 3D Laplace operator on a structured grid. On a
  given x-y plane, the stencil is a in "laplace 2d". The value on the
  plane below is set using Set("f",F), and in the plane above with
  Set("g",G"). Problem size specifed as in "laplace3d."

  - \c hb: The matrix is read from file. File name is specified by
  Set("file name", FileName). FileName is a C++ string. Problem size is
  automatically read from file.

  - \c lehmer: Returns a symmetric positive definite matrix, such that
  A(i,j) = (i+1)/(j+1) (for j>=i), and (j+1)/(i+1) (for j<j). This
  matrix has three properties:
     - is totally nonnegative;
     - the inverse is tridiagonal and esplicitly known;
     - The condition number is bounded as n <= cond(A) <= 4*n (n is the problem size).

  - \c minij: Returns the symmetric positive definite matrix defined sd
  A(i,j) = min(i+1,j+1). More information can be found in the help of
  MATLAB's gallery function.

  - \c ris: Returns a symmetric Hankel matrix with elements A(i,j) =
  0.5/(n-(i+1)-(j+1)+1.5),
  with n equals problem size. (i and j start from 0.)
  The eigenvalues of A cluster around -pi/2 and pi/2.

  - \c hilbert: This is a famous example of a badly conditioned matrix. The
  elements are defined as 1/((i+1)+(j+1)-1). (i and j start from 0.)

  - \c jordblock: Creates a Jordan block with eigenvalue set via Set("a",DoubleVal).

  - \c cauchy: Creates a particular instance of a Cauchy matrix with
  elements A(i,j) = 1/(i+j) (where i and j range from 1 to
  NumGlobalElements).  Explicit formulas are known for the inverse and
  determinand of a Cauchy matrix. For this particular Cauchy matrix, the
  determinant is nonzero and the matrix is totally positive.

  - \c fiedler: Creates a matrix whose element are A(i,j) = | i - j |.
  The matrix is symmetric, and has a dominant positive eigenvalue, and
  all the other eigenvalues are negative.

  - \c hanowa: Creates a matrix whoe eigenvalues lie on a vertical line
  in the complex plane. The matrix has the 2x2 block structure
  A = [ a * eye(n/2) ; -diag(1:m) ;
        diag(1:m);     a * eye(n/2) ];
  The complex eigenvalues are of the form a +- k*i, for 1 <= k <= n/2.
  The default value for a_ is -1;

  - \c kms: Create the n-by-n Kac-Murdock-Szego Toepliz matrix such that
  A(i,j) = rho^(abs(i-j)) (for real rho only). Default value for rho is
  0.5.  The inverse of this matrix is tridiagonal, and it is positive
  definite if and only if 0 < abs(rho) < 1.
  
  - \c parter: Creates a matrix A(i,j) = 1/(i-j+0.5). A is a Cauchy and
  a Toeplliz matrix. Most of the singular values of A are very close to
  pi.

  - \c pei: Creates the matrix alpha * eye(n) + ones(n). This matrix is
  singular for alpha == 0 or -n.

  - \c vander: Create the Vandermonde matrix whose second last column is
  the vector c. The j-th column is given by A(:,j) = C^(n-j).

  
  An example of program using this class is reported below.

  \code
int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Trilinos_Util_CrsMatrixGallery Gallery("hb", Comm);

  Gallery.Set("matrix name", "bcsstk14.rsa");
  
  Epetra_CrsMatrix * A;
  Epetra_Vector * ExactSolution;
  Epetra_Vector * RHS;
  Epetra_Vector * StartingSolution;
  
  A = Gallery.GetMatrix();
  ExactSolution = Gallery.GetExactSolution();
  RHS = Gallery.GetRHS();
  StartingSolution = Gallery.GetStartingSolution();
  
  // create linear problem
  Epetra_LinearProblem Problem(A,StartingSolution,RHS);
  // create AztecOO instance
  AztecOO Solver(Problem);

  Solver.SetAztecOption( AZ_precond, AZ_dom_decomp );  
  Solver.Iterate(1000,1E-9);

  // compute residual
  double residual;
  
  Gallery.ComputeResidual(residual);
  if( Comm.MyPID()==0 ) cout << "||b-Ax||_2 = " << residual << endl;
  
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(residual);
  if( Comm.MyPID()==0 ) cout << "||x_exact - x||_2 = " << residual << endl;

 #ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

return 0 ;
  } 
  \endcode

  Class Trilinos_Util_CommandLineParser can be used as well. In this case, one may
  decide to use the following:
  \code
  Trilinos_Util_CommandLineParser CLP(argc,argv);
  // set a problem with no matrix name
  Trilinos_Util_CrsMatrixGallery Gallery("", Comm);
  // read parameters and settings from the shell line
  G.Set(CLP);
  // continue with your code...
  \endcode
  
  \param In
  comm - Epetra communicator
  */
  Trilinos_Util_CrsMatrixGallery( const string name, const Epetra_Comm & comm );

  //! Creates an Triutils_Gallery object using a given map.
  /*! Create a Triutils_Gallery object using an Epetra_Map.
    Problem size must match the elements in map.

    \param In
    name - definition of the problem to be created.

    \param In
    map - Epetra_Map
  */
  Trilinos_Util_CrsMatrixGallery( const string name, const Epetra_Map & map );
  
  //! Triutils_Gallery destructor
  ~Trilinos_Util_CrsMatrixGallery();

  //@}

  //@{ \name Setting methods
 
  //! Sets a gallery options using an interger value.
  int Set(const string parameter, const int value);

  //!  Sets a gallery options using a C++ string .
  int Set(const string parameter, const string value );

  //! Sets a gallery options using an double value.
  int Set(const string parameter, const double value);

  //! Sets a gallery options using an Epetra_Vector.
  /*! Sets a gallery options using an Epetra_Vector. The Epetra_Vector
  is copied into internal structures, and freed by the destructor.
  */
  int Set(const string parameter, const Epetra_Vector & value);

  //! Sets gallery options using values passed from the shell
  int Set(Trilinos_Util_CommandLineParser & CLP);

  //@}

  //@{ \name Extraction methods.

  //! Returns a pointer to the CrsMatrix.
  Epetra_CrsMatrix * GetMatrix();

  Epetra_CrsMatrix & GetMatrixRef();

  //! Returns a pointer to the exact solution.
  /*! Returns a pointer to the exact solution.
    
    Some choices are available to define the exact solution, using
    Set("exact solution", value). value can be:
    - constant: the exact solution vector is made up of 1's.
    - random: a random solution vector
    - linear: value at node i is defined as alpha*i. The double value
    alpha can be set via Set("alpha",DoubleVal).
  */
  Epetra_Vector * GetExactSolution();

  //! Returns a pointer to the starting solution (typically, for HB problems).
  /*! Returns a pointer to the starting solution. This is typically used
    while reading a HB problem. However, the user can set a starting
    solution using Set("starting solution", "value"). Value can be
    - zero
    - random 
  */
  Epetra_Vector * GetStartingSolution();
  
  //! Returns a pointer to the rhs corresponding to the selected exact solution.
  Epetra_Vector * GetRHS();

  //! Returns a pointer the internally stored Map.
  const Epetra_Map * GetMap();

  const Epetra_Map & GetMapRef();  

  // ==================== //
  // LINEAR PROBLEM STUFF //
  // ==================== //

  //! Returns a pointer to Epetra_LinearProblem
  Epetra_LinearProblem * GetLinearProblem();

  //! Computes the 2-norm of the residual
  int ComputeResidual(double & residual);

  //! Computes the 2-norm of the difference between the starting solution and the exact solution
  int ComputeDiffBetweenStartingAndExactSolutions(double & residual);

  //! Print out matrix and vectors
  void PrintMatrixAndVectors(ostream & os);

  void PrintMatrixAndVectors();

  //! Get pointers to double vectors containing coordinates of points.
  int GetCartesianCoordinates(double * & x, double * & y, double * & z);
  
  //! Print out detailed information about the problem at hand
  friend ostream & operator << (ostream& os,
				const Trilinos_Util_CrsMatrixGallery & G );
  //! Print matrix on file in MATLAB format
  int WriteMatrix( const string & FileName, const bool UseSparse=true );
  
  //@}

private:

  //@{ \name Creation methods.
  
  //! Creates a map.
  /*! Creates an Epetra_Map. Before calling this function, the problem
  size must have been specified.

  CreateMap() allows some different maps. The type of map is set using
  Set("map",value). Value is a string, defined as:
  - linear: Creates a linear map. Elements are divided into continuous
  chunks among the processors.

  - box: used for problems defined on cartesian grids over a square. The
  nodes is subdivided into mx x my subdomains. mx and my are specified
  via Set("mx",IntValue) and Set("my",IntValue).

  - interlaces: elements are subdivided so that element i is assigned to
  process i%NumProcs.

  - random: assign each node to a random process
  
  - greedy: (only for HB matrices) implements a greedy algorithm to
    decompose the graph of the HB matrix among the processes
    
  */
  int CreateMap();
  
  //! Creates the CrdMatrix.
  int CreateMatrix();

  //! Creates the exact solution.
  int CreateExactSolution();

  //! Creates the starting solution.
  int CreateStartingSolution();

  //! Create the RHS corresponding to the desired exact solution.  
  int CreateRHS();
  
  // Create an identity matrix.
  int CreateEye();

  // Creates a diagonal matrix. Elements on the diagonal are called `a'.
  int CreateMatrixDiag();
    
  // Creates a tridiagonal matrix. Elements on the diagonal are called `a',
  // elements on the sub-diagonal 'b', and on the super-diagonal 'c'.
  int CreateMatrixTriDiag();
  
  // Create a matrix for a Laplacian in 1D
  int CreateMatrixLaplace1d();
  
  int CreateMatrixLaplace1dNeumann();
  
  int CreateMatrixCrossStencil2d();

  int CreateMatrixLaplace2d();

  int CreateMatrixLaplace2d_9pt();

  int CreateMatrixRecirc2d();

  int CreateMatrixLaplace2dNeumann();
  
  int CreateMatrixUniFlow2d();
  
  int CreateMatrixLaplace3d();

  int CreateMatrixCrossStencil3d();

  int CreateMatrixLehmer();

  int CreateMatrixMinij();

  int CreateMatrixRis();

  int CreateMatrixHilbert();

  int CreateMatrixJordblock();

  int CreateMatrixCauchy();

  int CreateMatrixFiedler();

  int CreateMatrixHanowa();

  int CreateMatrixKMS();
  
  int CreateMatrixParter();

  int CreateMatrixPei();

  int CreateMatrixOnes();

  int CreateMatrixVander();
  
  // read an HB matrix. This function requires other Trilinos util files
  int ReadMatrix();

  // Creates a block map, based on map, wich NumPDEEqns equations on each node.
  int CreateBlockMap(void);

  // returns the neighbors of a given node. The node is supposed to be on
  // a 2D Cartesian grid 
  void  GetNeighboursCartesian2d( const int i, const int nx, const int ny,
				  int & left, int & right, 
				  int & lower, int & upper);
  // returns the neighbors of a given node. The node is supposed to be on
  // a 3D Cartesian grid   
  void  GetNeighboursCartesian3d( const int i, const int nx, const int ny, const int nz,
				  int & left, int & right, int & lower, int & upper,
				  int & below, int & above );

  // put to NULL or default values all internal data
  void ZeroOutData();

  //@}
  
  // ======================== //
  // I N T E R N A L  D A T A //
  // ======================== //
  
  const Epetra_Comm * comm_;

  // matrix and vectors (scalar)
  Epetra_CrsMatrix * matrix_;
  Epetra_Vector * ExactSolution_;
  Epetra_Vector * StartingSolution_;
  Epetra_Vector * rhs_;
  Epetra_Map * map_;

  // linear problem
  Epetra_LinearProblem * LinearProblem_;

  // information about the problem to generate
  string name_;
  int NumGlobalElements_;
  int NumMyElements_;
  int * MyGlobalElements_;
  string MapType_;
  string ExactSolutionType_;
  string StartingSolutionType_;
  string ExpandType_;
  
  // parameters
  int nx_, ny_, nz_;
  int mx_, my_, mz_;
  
  int NumPDEEqns_;
  
  Epetra_Vector * VectorA_, * VectorB_, * VectorC_, * VectorD_, * VectorE_, *VectorF_, * VectorG_;
  
  double a_, b_, c_, d_, e_, f_, g_;
  double alpha_, beta_, gamma_, delta_;
  
  string FileName_;

  // others
  string ErrorMsg;
  string OutputMsg;
  bool verbose_;
  
};

#ifdef LATER_ON

class Trilinos_Util_VbrMatrixGallery
{

public:

  Trilinos_Util_VbrMatrixGallery(Trilinos_Util_CrsMatrixGallery & Gallery) :
    Gallery_(Gallery) {} ;

~Trilinos_Util_VbrMatrixGallery() 
{
    // VBR data
  if( VbrLinearProblem_ != NULL ) delete VbrLinearProblem_;
  if( VbrMatrix_ != NULL ) delete VbrMatrix_;  
  if( VbrExactSolution_ != NULL ) delete VbrExactSolution_;
  if( VbrStartingSolution_ != NULL ) delete VbrStartingSolution_;
  if( VbrRhs_ != NULL ) delete VbrRhs_;
  if( BlockMap_ != NULL ) delete BlockMap_;
  
}
  
  // ========= //
  // VBR STUFF //
  // ========= //
  
  //! Returns a pointer the internally stored BlockMap.
  const Epetra_BlockMap * GetBlockMap();

  const Epetra_BlockMap & GetBlockMapRef();
  
  //! Returns a VbrMatrix, starting from the CsrMatrix.
  /*! Returns a VbrMatrix, starting from the CsrMatrix. This vbr matrix
    is formally equivalent to the CrsMatrix returned by
    GetMatrix(). However, each node of the CrsMatrix is replicated
    num_PDE_eqns times (this value is passed in input, or set via Set("num pde
    eqns",IntValue)).
  */
  Epetra_VbrMatrix * GetVbrMatrix(const int NumPDEEqns);

  //! Returns a VbrMatrix, starting from the CsrMatrix.
  Epetra_VbrMatrix * GetVbrMatrix();

  Epetra_VbrMatrix & GetVbrMatrixRef();

  //! Returns a pointer to the RHS for the selected Vbr exact solution
  /*!  Returns a pointer to the RHS  corresponding to the selected exact solution to the linear systems defined by the Epetra_VbrMatrix.
   */
  Epetra_Vector * GetVbrRHS();

  //! Returns a pointer to the selected Vbr exact solution
  Epetra_Vector * GetVbrExactSolution();

  //! Returns a pointer to the starting solution for Vbr problems
  Epetra_Vector * GetVbrStartingSolution();


  // create the Vbr matrix. 
  int CreateVbrMatrix(void);  

  //! Returns a pointer to Epetra_LinearProblem for VBR
  Epetra_LinearProblem * GetVbrLinearProblem();

  //! Computes the 2-norm of the residual for the VBR problem
  int ComputeResidualVbr(double & residual);

  //! Computes the 2-norm of the difference between the starting solution and the exact solution for the VBR problem  
  int ComputeDiffBetweenStartingAndExactSolutionsVbr(double & residual);

  //! Print out Vbr matrix and vectors
  void PrintVbrMatrixAndVectors(ostream & os);

  void PrintVbrMatrixAndVectors();

private:

  
  //! Creates the exact solution for a Epetra_VbrMatrix.
  int CreateVbrExactSolution(void);

  //! Creates the starting solution for Vbr.
  int CreateVbrStartingSolution();
  
  //!  Create the RHS corresponding to the desired exact solution for the Vbr problem.
  int CreateVbrRHS();

  // linear problem  
  Epetra_LinearProblem * VbrLinearProblem_;

  // matrix and vectors (vbr)
  Epetra_VbrMatrix * VbrMatrix_;
  Epetra_Vector * VbrExactSolution_;
  Epetra_Vector * VbrStartingSolution_;
  Epetra_Vector * VbrRhs_;
  Epetra_BlockMap * BlockMap_;
  int MaxBlkSize_;

Trilinos_Util_CrsMatrixGallery * Gallery_;

};
#endif


#endif

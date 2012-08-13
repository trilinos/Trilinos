/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

// Goal of this example if to show how to take advantage of the getrow()
// to define matrix-free operator using the MultiLevelPreconditioner class.
//
// \author Marzio Sala, SNL 9214
// \date Last modified on 25-Mar-05

#include "ml_include.h"

// the following code cannot be compiled without these Trilinos
// packages.
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Import.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Map.h"
#include "AztecOO.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra.h"

using namespace Teuchos;
using namespace ML_Epetra;

// ================================================================== //
// This class is an example of matrix free format.                    //
// It represents a 3D Laplacian on a structured Cartesian grid; the   //
// problem is discretized using a classical 7-point formula.          //
// The computational domain is the unitary cube; Dirichlet boundary   //
// conditions are applied on \partial \Omega. The Dirichlet boundary  //
// nodes have been eliminated from the matrix.                        //
//                                                                    //
// The nonzero elements are queried through a suitable getrow()       //
// function, here given by ExtractMyRowCopy(). To specify the         //
// communication pattern, we proceed as follows:                      //
// - an Epetra_CrsGraph is created, filled, the maps and import       //
//   objects are copied in the class;                                 //
// - the Import() of Epetra_MultiVector is used in method Multiply(). //
// - recall that Aztec and ML require the elements returned by the    //
//   getrow() to be ordered so that external elements come after      //
//   local elements.                                                  //
// ================================================================== //

class Laplace3D : public virtual Epetra_RowMatrix {
      
public:
  //! Constructor.
  /*!
   * \param nx (In) - number of nodes along the X-axis
   *
   * \param ny (In) - number of nodes along the Y-axis
   *
   * \param nz (In) - number of nodes along the Z-axis
   *
   * \param mx (In) - number of subdomains along the X-axis.
   *
   * \param my (In) - number of subdomains along the Y-axis.
   *
   * \param mz (In) - number of subdomains along the Z-axis.
   *
   * \note The global number of unknowns is given by nx * ny * nz.
   *       The total number of subdomains (that is, processes) is given
   *       by mx * my * mz.
   */
  Laplace3D(Epetra_Comm& InputComm,
            const int nx, const int ny, const int nz,
            const int mx, const int my, const int mz,
            const bool BuildCoord = false) :
    Comm_(InputComm),
    nx_(nx),
    ny_(ny),
    nz_(nz),
    mx_(mx),
    my_(my),
    mz_(mz)
  {

    // creates the distribution of rows across the processors,
    // by subdividing the cube into smaller cubes, as specified
    // by the mx, my, mz parameters.
    
    CreateMap();

    // We need to compute the communication pattern. To that aim,
    // we build a graph (using the Epetra_CrsGraph class), call
    // FillComplete() on this class, then extract all the necessary
    // information, including the importer of off-processor
    // elements. Then, we delete the graph.
    
    int left, right, lower, upper, below, above;

    Epetra_CrsGraph Graph(Copy,*RowMap_,7);

    int Indices[7];

    int NumMyElements = RowMap_->NumMyElements();
    int* MyGlobalElements = RowMap_->MyGlobalElements();

    for (int i = 0 ; i < NumMyElements; ++i) 
    {
      Indices[0] = MyGlobalElements[i];
      int NumEntries=1;
      GetNeighboursCartesian3d(MyGlobalElements[i], left, right, 
                               lower, upper, below, above);
      if( left != -1 ) {
        Indices[NumEntries] = left;
        ++NumEntries;
      }
      if( right != -1 ) {
        Indices[NumEntries] = right;
        ++NumEntries;
      }
      if( lower != -1 ) {
        Indices[NumEntries] = lower;
        ++NumEntries;
      }
      if( upper != -1 ) {
        Indices[NumEntries] = upper;
        ++NumEntries;
      }
      if( below != -1 ) {
        Indices[NumEntries] = below;
        ++NumEntries;
      }
      if( above != -1 ) {
        Indices[NumEntries] = above;
        ++NumEntries;
      }

      Graph.InsertGlobalIndices(MyGlobalElements[i], NumEntries, Indices);
    }

    Graph.FillComplete();
    Graph.OptimizeStorage();

    ColMap_ = new Epetra_Map(-1, Graph.ColMap().NumMyElements(),
                             Graph.ColMap().MyGlobalElements(), 0, Comm());

    Importer_ = new Epetra_Import(*ColMap_, *RowMap_);

    // extract all the information from the graph

    NumMyRows_ = Graph.NumMyRows();
    NumMyCols_ = Graph.NumMyCols();
    NumGlobalRows_ = Graph.NumGlobalRows();
    NumGlobalCols_ = Graph.NumGlobalCols();
    NumMyDiagonals_ = Graph.NumMyDiagonals();
    NumGlobalDiagonals_ = Graph.NumGlobalDiagonals();
    NumMyNonzeros_ = Graph.NumMyNonzeros();
    NumGlobalNonzeros_ = Graph.NumGlobalNonzeros();
    MaxNumIndices_ = Graph.MaxNumIndices();

    RowEntries_.resize(NumMyRows());
    RowIndices_.resize(NumMyRows());

    for (int i = 0 ; i < NumMyRows() ; ++i)
    {
      int NumEntries = 0;
      RowEntries_[i] = Graph.NumMyIndices(i);
      if (RowEntries_[i] != 0)
      {
        RowIndices_[i].resize(RowEntries_[i]);
        Graph.ExtractMyRowCopy(i, RowEntries_[i], NumEntries, 
                               &(RowIndices_[i][0]));
      }
    }
    
    // from now on, the Graph object is no longer necessary.
    
    if (BuildCoord)
      CreateCoordinates();
  }

  virtual ~Laplace3D()
  {}

  virtual int NumMyRowEntries(int MyRow, int & NumEntries) const
  {
    NumEntries = RowEntries_[MyRow];
    return(0);
  }

  virtual int MaxNumEntries() const
  {
    return(MaxNumIndices_);
  }

  virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
                               double *Values, int * Indices) const
  {
    ML_RETURN(getrow(MyRow, Length, NumEntries, Values, Indices));
  }

  virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
  {
    for (int i = 0 ; i < Diagonal.MyLength() ; ++i)
      Diagonal[i] = 7.0;
    return(0);
  }

  virtual int Multiply(bool TransA, const Epetra_MultiVector& X, 
                       Epetra_MultiVector& Y) const
  {
    Epetra_MultiVector Xtmp(RowMatrixColMap(), X.NumVectors());
    Xtmp.Import(X, *RowMatrixImporter(), Insert);

    std::vector<int> Indices(MaxNumEntries());
    std::vector<double> Values(MaxNumEntries());

    Y.PutScalar(0.0);

    for (int i = 0 ; i < NumMyRows() ; ++i)
    {
      int NumEntries;
      // use the inlined function
      getrow(i, MaxNumEntries(), NumEntries, &Values[0], &Indices[0]);
      for (int j = 0 ; j < NumEntries ; ++j)
        for (int k = 0 ; k < Y.NumVectors() ; ++k)
          Y[k][i] += Values[j] * Xtmp[k][Indices[j]];
    }
    
    return(0);
  }
  

  virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, 
                    const Epetra_MultiVector& X, 
                    Epetra_MultiVector& Y) const
  {
    ML_RETURN(-1); // not implemented 
  }

  virtual int Apply(const Epetra_MultiVector& X,
                    Epetra_MultiVector& Y) const
  {
    ML_RETURN(Multiply(false,X,Y));
  }

  virtual int ApplyInverse(const Epetra_MultiVector& X,
                           Epetra_MultiVector& Y) const
  {
    ML_RETURN(-1);
  }

  virtual int InvRowSums(Epetra_Vector& x) const
  {
    ML_RETURN(-1); // not implemented
  }

  virtual int LeftScale(const Epetra_Vector& x)
  {
    ML_RETURN(-1); // not implemented
  }

  virtual int InvColSums(Epetra_Vector& x) const
  {
    ML_RETURN(-1); // not implemented
  }

  virtual int RightScale(const Epetra_Vector& x) 
  {
    ML_RETURN(-1); // not implemented
  }

  virtual bool Filled() const
  {
    return(true);
  }

  virtual double NormInf() const
  {
    return(14.0);
  }

  virtual double NormOne() const
  {
    return(14.0);
  }

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
//TODO:CJ:  change data type for class members. check if changing is needed.

  virtual int NumGlobalNonzeros() const
  {
    return(NumGlobalNonzeros_);
  }

  virtual int NumGlobalRows() const
  {
    return(NumGlobalRows_);
  }

  virtual int NumGlobalCols() const
  {
    return(NumGlobalRows_);
  }

  virtual int NumGlobalDiagonals() const
  {
    return(NumGlobalDiagonals_);
  }
#endif

  virtual long long NumGlobalNonzeros64() const
  {
    return(NumGlobalNonzeros_);
  }

  virtual long long NumGlobalRows64() const
  {
    return(NumGlobalRows_);
  }

  virtual long long NumGlobalCols64() const
  {
    return(NumGlobalRows_);
  }

  virtual long long NumGlobalDiagonals64() const
  {
    return(NumGlobalDiagonals_);
  }
    
  virtual int NumMyNonzeros() const
  {
    return(NumMyNonzeros_);
  }

  virtual int NumMyRows() const
  {
    return(NumMyRows_);
  }

  virtual int NumMyCols() const
  {
    return(NumMyCols_);
  }

  virtual int NumMyDiagonals() const
  {
    return(NumMyDiagonals_);
  }

  virtual bool LowerTriangular() const
  {
    return(false);
  }

  virtual bool UpperTriangular() const
  {
    return(false);
  }

  virtual const Epetra_Map & RowMatrixRowMap() const
  {
    return(*RowMap_);
  }

  virtual const Epetra_Map & RowMatrixColMap() const
  {
    return(*ColMap_);
  }

  virtual const Epetra_Import * RowMatrixImporter() const
  {
    return(Importer_);
  }
  
  int SetOwnership(bool ownership)
  {
    return(-1);
  }

  int SetUseTranspose(bool UseTranspose){return(-1);}

  bool UseTranspose() const {return(false);};
  
  bool HasNormInf() const{return(true);};
  
  const Epetra_Comm& Comm() const
  {
    return(Comm_);
  }
  
  const Epetra_Map & OperatorDomainMap() const 
  {
    return(*RowMap_);
  }
  
  const Epetra_Map & OperatorRangeMap() const 
  {
    return(*RowMap_);
  }

  void SetLabel(const char* label) 
  {
    strcpy(Label_,label);
  }

  const char* Label() const{
    return(Label_);
  }

  const Epetra_BlockMap & Map() const
  {
    return(*RowMap_);
  }

  inline int nx() const
  {
    return(nx_);
  }

  inline int ny() const
  {
    return(ny_);
  }

  inline int nz() const
  {
    return(nz_);
  }

  inline int mx() const
  {
    return(mx_);
  }

  inline int my() const
  {
    return(my_);
  }

  inline int mz() const
  {
    return(mz_);
  }

  inline const double* XCoord() const
  {
    return(&x_coord[0]);
  }
    
  inline const double* YCoord() const
  {
    return(&y_coord[0]);
  }
    
  inline const double* ZCoord() const
  {
    return(&z_coord[0]);
  }
    
private:

  void CreateCoordinates()
  {
    x_coord.resize(NumMyRows());
    y_coord.resize(NumMyRows());
    z_coord.resize(NumMyRows());

    double lx = 1.0, ly = 1.0, lz = 1.0;

    double delta_x = lx / (nx() - 1);
    double delta_y = ly / (ny() - 1);
    double delta_z = lz / (nz() - 1);

    for (int i = 0 ; i < NumMyRows() ; i++) 
    {
      int GlobalRow = RowMatrixRowMap().GID(i);
      int ixy = GlobalRow % (nx() * ny());
      int iz = (GlobalRow - ixy) / (nx() * ny());

      int ix = ixy % nx();
      int iy = (ixy - ix) / ny();
      
      x_coord[i] = delta_x * ix;
      y_coord[i] = delta_y * iy;
      z_coord[i] = delta_z * iz;
    }
  }

  //! Interal getrow, to be inlined.
  inline int getrow(int MyRow, int Length, int & NumEntries, 
                    double *Values, int * Indices) const
  {
    if (Length < RowEntries_[MyRow])
      return(-1);

    NumEntries = RowEntries_[MyRow];

    for (int i = 0 ; i < NumEntries ; ++i)
    {
      Indices[i] = RowIndices_[MyRow][i];
      if (Indices[i] == MyRow)
        Values[i] = 7.0;
      else
        Values[i] = -1.0;
    }

    return (0);
  }

  //! Defines the neighbors in a 2D Cartesian grid.
  void GetNeighboursCartesian2d(const int i, int& left, int& right, 
                                int& lower, int& upper) 
  {
    int ix, iy;
    ix = i % nx();
    iy = (i - ix) / nx();

    if (ix == 0) 
      left = -1;
    else 
      left = i-1;
    if (ix == nx() - 1) 
      right = -1;
    else
      right = i + 1;
    if (iy == 0) 
      lower = -1;
    else
      lower = i - nx();
    if (iy == ny() - 1) 
      upper = -1;
    else
      upper = i + nx();

    return;
  }

  //! Defines the neighbors in a 3D Cartesian grid.
  void GetNeighboursCartesian3d(const int i, int& left, int& right, 
                                int& lower, int& upper,
                                int& below, int& above ) 
  {
    int ixy, iz;
    ixy = i % (nx() * ny());

    iz = (i - ixy)/(nx() * ny());

    if (iz == 0) 
      below = -1;
    else 
      below = i - nx() * ny();
    if (iz == nz() - 1) 
      above = -1;
    else
      above = i + nx() * ny();

    GetNeighboursCartesian2d(ixy, left, right, lower, upper);

    if (left  != -1) left  += iz * (nx() * ny());
    if (right != -1) right += iz * (nx() * ny());
    if (lower != -1) lower += iz * (nx() * ny());
    if (upper != -1) upper += iz * (nx() * ny());

    return;
  }

  void CreateMap() 
  {
    int modx = (nx() + (nx() % mx())) / mx();
    int mody = (ny() + (ny() % my())) / my();
    int modz = (nz() + (nz() % mz())) / mz();

    int MyPID = Comm().MyPID(), startx, starty, startz, endx, endy, endz;
    int mxy  = mx() * my();
    int zpid = MyPID / mxy;
    int xpid = (MyPID % mxy) % mx();
    int ypid = (MyPID % mxy) / mx();

    startx = xpid * modx;
    if ((xpid+1)*modx < nx()) endx = (xpid + 1) * modx;
    else endx = nx();
    starty = ypid*mody;
    if ((ypid+1)*mody < ny()) endy = (ypid + 1) * mody;
    else endy = ny();
    startz = zpid*modz;
    if ((zpid+1)*modz < nz()) endz = (zpid + 1) * modz;
    else endz = nz();

    int NumMyElements = (endx - startx) * (endy - starty) * (endz - startz);
    std::vector<int> MyGlobalElements(NumMyElements);
    int count = 0;

    for (int i = startx ; i < endx ; ++i) {
      for (int j = starty ; j < endy ; ++j) {
        for (int k = startz ; k < endz ; ++k) {
          MyGlobalElements[count++] = i + j * nx() + k * (nx() * ny());
        }
      }
    }

    RowMap_ = new Epetra_Map (-1,NumMyElements,&MyGlobalElements[0],0,Comm());
  }

  //! Communicator.
  Epetra_Comm& Comm_;
  //! Total number of nodes along the X-axis.
  int nx_;
  //! Total number of nodes along the Y-axis.
  int ny_;
  //! Total number of nodes along the Z-axis.
  int nz_;
  //! Total number of processors along the X-axis.
  int mx_;
  //! Total number of processors along the Y-axis.
  int my_;
  //! Total number of processors along the Z-axis.
  int mz_;
  //! Number of local rows.
  int NumMyRows_;
  //! Number of local columns.
  int NumMyCols_;
  //! Number of global rows.
  int NumGlobalRows_;
  //! Number of global columns.
  int NumGlobalCols_;
  //! Number of local nonzero diagonals.
  int NumMyDiagonals_;
  //! Number of global nonzero diagonals.
  int NumGlobalDiagonals_;
  //! Maximum number of nonzero elements per row.
  int MaxNumIndices_;
  //! Number of local nonzero elements.
  int NumMyNonzeros_;
  //! Number of global nonzero elements.
  int NumGlobalNonzeros_;
  //! Label of \c this object.
  char Label_[80];
  //! Row map.
  Epetra_Map* RowMap_;
  //! Column map.
  Epetra_Map* ColMap_;
  //! Importer from RowMap() to ColMap().
  Epetra_Import* Importer_;
  //! Contains the already reordered indices for each local row.
  std::vector<std::vector<int> > RowIndices_;
  //! Contains the number of nonzeros in each local row.
  std::vector<int> RowEntries_;
  std::vector<double> x_coord;
  std::vector<double> y_coord;
  std::vector<double> z_coord;

}; // class Laplace3D

// ============================================================ //
// example driver                                               //
//                                                              //
// This example will be run as follows:                         //
//                                                              //
// $ mpirun -np 4 ./ml_MatrixFree.exe --n=<n> --m=<m>           //
//                                                              //
// where <n> * <m> is the number of nodes along each axis and   //
// <m> is the number of subdomains (processes) along each axis. //
// The example requires the number of processors to be a cube   //
// of an integer; however, it is simple to extend the example   //
// to more general cases.                                       //
// ============================================================ //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  CommandLineProcessor CLP;

  int n = 10;
  int m = (int)pow((double)Comm.NumProc(), 0.3334);
  double DampingFactor = 1.333;
  std::string AggregationScheme = "Uncoupled";
  int NPA = 16;
  int MaxLevels = 5;

  CLP.setOption("l", &MaxLevels, "number of levels");
  CLP.setOption("n", &n, "number of nodes along each axis");
  CLP.setOption("damp", &DampingFactor, "prolongator damping factor");
  CLP.setOption("aggr", &AggregationScheme, "aggregation scheme");
  CLP.setOption("npa", &NPA, "nodes per aggregate (if supported by the aggregation scheme)");

  CLP.throwExceptions(false);
  CLP.parse(argc,argv);

  if (m * m * m != Comm.NumProc()) {
    if (Comm.MyPID() == 0) 
    {
      std::cout << "Number of processes must be a perfect cube." << std::endl;
      std::cout << "Please re-run with --help option for details." << std::endl;
    }
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_SUCCESS);
  }
    
  n *= m;

  Laplace3D A(Comm, n, n, n, m, m, m, true);
  Epetra_Vector LHS(A.OperatorDomainMap());
  Epetra_Vector RHS(A.OperatorDomainMap());

  LHS.Random();
  RHS.PutScalar(0.0);

  Epetra_LinearProblem Problem(&A, &LHS, &RHS);
  // Construct a solver object for this problem
  AztecOO solver(Problem);

  // create a parameter list for ML options
  ParameterList MLList;

  // set defaults for classic smoothed aggregation
  ML_Epetra::SetDefaults("SA",MLList);
  
  // overwrite some parameters. Please refer to the user's guide
  // for more information
  MLList.set("max levels",MaxLevels);
  MLList.set("increasing or decreasing","increasing");
  MLList.set("aggregation: type", AggregationScheme);
  MLList.set("aggregation: damping factor", DampingFactor);
  MLList.set("aggregation: nodes per aggregate", NPA);
  MLList.set("smoother: type","symmetric Gauss-Seidel");
  MLList.set("smoother: pre or post", "both");
  MLList.set("coarse: max size", 512);
  MLList.set("coarse: type","Amesos-KLU");
  MLList.set("analyze memory", true);
  MLList.set("repartition: enable", true);
  MLList.set("repartition: max min ratio", 1.1);
  MLList.set("repartition: min per proc", 512);
  MLList.set("low memory usage", true);
  MLList.set("x-coordinates", (double*) A.XCoord());
  MLList.set("y-coordinates", (double*) A.YCoord());
  MLList.set("z-coordinates", (double*) A.ZCoord());
  
  // create the preconditioner object and compute hierarchy
  MultiLevelPreconditioner* MLPrec = new MultiLevelPreconditioner(A, MLList);

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, 1);
  solver.Iterate(500, 1e-10);

  delete MLPrec;
  
  double norm;
  LHS.Norm2(&norm);

  if (Comm.MyPID() == 0) 
    std::cout << "Norm of the error = " << norm << std::endl;

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  exit(EXIT_SUCCESS);
  
}

#else

#include <cstdlib>
#include <cstdio>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with:");
  puts("--enable-epetra");
  puts("--enable-teuchos");
  puts("--enable-aztecoo");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
#endif


//@HEADER
// ************************************************************************
//
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Export.h"
#include "Epetra_Time.h"
#include "Trilinos_Util.h"
#include "Trilinos_Util_MatrixGallery.h"

#include <string>

// ================================================ ====== ==== ==== == =
Trilinos_Util_MatrixGallery::Trilinos_Util_MatrixGallery( string name, Epetra_Comm & comm ) :
  name_(name), comm_(&comm)
{
  ZeroOutData();
  // verbosity level
  if( comm_->MyPID()==0 ) verbose_ = true;
  else verbose_ = false;
  // fix error message
  ErrorMsg = "ERROR (Trilinos_Util_MatrixGallery)";
  OutputMsg = "Triutils_Gallery: ";
  
}

// ================================================ ====== ==== ==== == =
Trilinos_Util_MatrixGallery::Trilinos_Util_MatrixGallery( string name, Epetra_Map & map ) :
  name_(name), comm_(&(map.Comm()))
{
  ZeroOutData();
  // verbosity level
  if( comm_->MyPID()==0 ) verbose_ = true;
  else verbose_ = false;
  // fix error message
  ErrorMsg = "ERROR (Trilinos_Util_MatrixGallery)";
  OutputMsg = "Trilinos_Util_MatrixGallery: ";
  
  map_ = new Epetra_Map(map);
  NumGlobalElements_ = map_->NumGlobalElements();
  NumMyElements_ = map_->NumMyElements();
  MyGlobalElements_ = map_->MyGlobalElements( );
    
}
  
// ================================================ ====== ==== ==== == =
Trilinos_Util_MatrixGallery::~Trilinos_Util_MatrixGallery() 
{
  // Crs data
  if( matrix_ != NULL ) delete matrix_;
  if( ExactSolution_ != NULL ) delete ExactSolution_;
  if( StartingSolution_ != NULL ) delete StartingSolution_;
  if( rhs_ != NULL ) delete rhs_;
  if( map_ != NULL ) delete map_;

  // VBR data
  if( BlockMap_ != NULL ) delete BlockMap_;
  if( VbrMatrix_ != NULL ) delete VbrMatrix_;
  if( VbrExactSolution_ != NULL ) delete VbrExactSolution_;
  if( VbrStartingSolution_ != NULL ) delete VbrStartingSolution_;
  if( VbrRhs_ != NULL ) delete VbrRhs_;
    
  // vectors
  if( VectorA_ != NULL ) delete VectorA_;
  if( VectorB_ != NULL ) delete VectorB_;
  if( VectorC_ != NULL ) delete VectorC_;
  if( VectorD_ != NULL ) delete VectorD_;
  if( VectorE_ != NULL ) delete VectorE_;
  if( VectorF_ != NULL ) delete VectorF_;
  if( VectorG_ != NULL ) delete VectorG_;

  // put to default values
  ZeroOutData();
}
  
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::Set(string parameter, int value)
{
  if( parameter == "problem_size" ) {
    if( value <= 0 ) {
      cerr << ErrorMsg << "problem size must be greater than 1\n";
      return -1;
    }
    if( map_ != NULL ) {
      cerr << ErrorMsg << "map object already set. Continuing with\n"
	   << ErrorMsg << "problemSize = " << NumGlobalElements_ << endl;
      return -2;
    }
    NumGlobalElements_ = value;
    return 0;

  } else if( parameter == "nx" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "nx must be greater than 0\n";
      return -1;
    }

    nx_ = value;
    return 0;
      
  } else if( parameter == "ny" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "ny must be greater than 0\n";
      return -1;
    }

    ny_ = value;
    return 0;
      
  } else if( parameter == "nz" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "nz must be greater than 0\n";
      return -1;
    }

    nz_ = value;
    return 0;
  } else if( parameter == "mx" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "mx must be greater than 0\n";
      return -1;
    }

    mx_ = value;
    return 0;
      
  } else if( parameter == "my" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "my must be greater than 0\n";
      return -1;
    }

    my_ = value;
    return 0;
      
  } else if( parameter == "mz" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "mz must be greater than 0\n";
      return -1;
    }

    mz_ = value;
    return 0;

  } else if( parameter == "num_pde_eqns" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "num pde eqns must be greater than 0\n";
      return -1;
    }

    NumPDEEqns_ = value;
    return 0;
  } 
    
  cerr << ErrorMsg << "input string (" << parameter << ") not valid\n";
  return -2;

}
  
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::Set(string parameter, string value )
{
  if( parameter == "problem_type" ) {
    name_ = value;
  }
  else if( parameter == "map_type" ) {
    MapType_ = value;
  }
  else if( parameter == "exact_solution" ) {
    ExactSolutionType_ = value;
  }
  else if( parameter == "matrix_name" ) {
    FileName_ = value;
  }    
  else if( parameter == "starting_solution" ) {
    StartingSolutionType_ = value;
  }
  else if( parameter == "output" ) {
    if( value == "none" ) verbose_ = false;
    if( value == "proc 0" ) {
      if( comm_->MyPID()==0 ) verbose_ = true;
      else verbose_ = false;
    } else {
      verbose_ = true;
    }
  } else {
    cerr << ErrorMsg << "wrong input parameter (" << parameter << ")\n";
    return -1;
  }
  return 0;
}
  
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::Set(string parameter, double value)
{
  if( parameter == "a" ) {
    a_ = value;
    return 0;
  } else if( parameter == "b" ) {
    b_ = value;
    return 0;
  } else if( parameter == "c" ) {
    c_ = value;
    return 0;
  } else if( parameter == "d" ) {
    d_ = value;
    return 0;
  } else if( parameter == "e" ) {
    e_ = value;
    return 0;
  } else if( parameter == "f" ) {
    f_ = value;
    return 0;
  } else if( parameter == "g" ) {
    g_ = value;
    return 0;
  }

  cerr << ErrorMsg << "input string not valid\n";
  return -2;
}
  
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::Set(string parameter, Epetra_Vector & value)
{

  if( value.Map().SameAs(*map_) == false ) {
    cerr << ErrorMsg << "input vector must have the same map used to\n"
	 << ErrorMsg << "create the Trilinos_Util_MatrixGallery object. Continuing\n";
    return -2;
  }
    
  if( parameter == "a" ) {
    VectorA_ = new Epetra_Vector(value);
    return 0;
  } else if( parameter == "b" ) {
    VectorB_ = new Epetra_Vector(value);
    return 0;
  }
  else if( parameter == "c" ) {
    VectorC_ = new Epetra_Vector(value);
    return 0;
  }
  else if( parameter == "d" ) {
    VectorD_ = new Epetra_Vector(value);
    return 0;
  }
  else if( parameter == "e" ) {
    VectorE_ = new Epetra_Vector(value);
    return 0;
  }
  else if( parameter == "f" ) {
    VectorF_ = new Epetra_Vector(value);
    return 0;
  }
  else if( parameter == "g" ) {
    VectorG_ = new Epetra_Vector(value);
    return 0;
  }

  cerr << ErrorMsg << "input string not valid\n";
  return -3;
}

// ================================================ ====== ==== ==== == =  
#ifdef TRILINOS_UTIL_MATRIX_GALLERY_WITH_SHELL_OPTIONS
#include <set>

int Trilinos_Util_MatrixGallery::Set(Trilinos_Util_ShellOptions & S)
{
  
  // iterate on this set
  set<string>::const_iterator current;

  // all options with strings
  set<string> StrOptions;
  StrOptions.insert("problem_type");
  StrOptions.insert("map_type");
  StrOptions.insert("exact_solution");
  StrOptions.insert("matrix_name");
  StrOptions.insert("starting_solution");
  StrOptions.insert("output");
  
  for( current=StrOptions.begin() ; current != StrOptions.end();
       ++current ) {
    string parameter = "-"+*current;    
    if( S.HaveOption(parameter) == true ) {
      string value = S.GetStringOption(parameter);
      Set(*current,value);
      
    }
  }

  // all options with integers
  set<string>  IntOptions;
  IntOptions.insert("problem_size");
  IntOptions.insert("nx");
  IntOptions.insert("ny");
  IntOptions.insert("nz");
  IntOptions.insert("mx");
  IntOptions.insert("my");
  IntOptions.insert("mz");
  IntOptions.insert("num_pde_eqns");

  for( current=IntOptions.begin() ; current != IntOptions.end();
       ++current ) {
    string parameter = "-"+*current;    
    if( S.HaveOption(parameter) == true ) {
      Set(*current,S.GetIntOption(parameter));
    }
  }
  
  // all options with doubles
  set<string>  DoubleOptions;
  DoubleOptions.insert("a");
  DoubleOptions.insert("b");
  DoubleOptions.insert("c");
  DoubleOptions.insert("d");
  DoubleOptions.insert("e");
  DoubleOptions.insert("f");
  DoubleOptions.insert("g");
  for( current=DoubleOptions.begin() ; current != DoubleOptions.end();
       ++current ) {
    string parameter = "-"+*current;    
    if( S.HaveOption(parameter) == true ) {
      Set(*current,S.GetDoubleOption(parameter));
    }
  }

  return 0;
}
#endif

// ================================================ ====== ==== ==== == =  
Epetra_CrsMatrix * Trilinos_Util_MatrixGallery::GetMatrix() 
{
  if( matrix_ == NULL ) CreateMatrix();
  return matrix_;
}
  
// ================================================ ====== ==== ==== == =
Epetra_Vector * Trilinos_Util_MatrixGallery::GetExactSolution() 
{
  if( ExactSolution_ == NULL ) CreateExactSolution();
  return ExactSolution_;
}

// ================================================ ====== ==== ==== == =
Epetra_Vector * Trilinos_Util_MatrixGallery::GetStartingSolution() 
{
  if( StartingSolution_ == NULL ) CreateStartingSolution();
  return StartingSolution_;
}

// ================================================ ====== ==== ==== == =
Epetra_Vector * Trilinos_Util_MatrixGallery::GetRHS() 
{
  if( rhs_ == NULL ) CreateRHS();
  return rhs_;
}

// ================================================ ====== ==== ==== == =
Epetra_Vector * Trilinos_Util_MatrixGallery::GetVbrRHS() 
{
  if( VbrRhs_ == NULL ) CreateVbrRHS();
  return VbrRhs_;
}

// ================================================ ====== ==== ==== == =
Epetra_Vector * Trilinos_Util_MatrixGallery::GetVbrExactSolution() 
{
  if( VbrExactSolution_ == NULL ) CreateVbrExactSolution();
  return VbrExactSolution_;
}

// ================================================ ====== ==== ==== == =
const Epetra_Map & Trilinos_Util_MatrixGallery::GetMap() 
{
  if( map_ == NULL ) CreateMap();
    
  return *map_;
}

// ================================================ ====== ==== ==== == =
const Epetra_BlockMap & Trilinos_Util_MatrixGallery::GetBlockMap() 
{
  if( BlockMap_ == NULL ) CreateBlockMap();
    
  return *BlockMap_;
}

// ================================================ ====== ==== ==== == =  
Epetra_VbrMatrix * Trilinos_Util_MatrixGallery::GetVbrMatrix(int NumPDEEqns) 
{

  if( NumPDEEqns != NumPDEEqns_ ) {
    if( BlockMap_ != NULL ) {
      delete BlockMap_;
      BlockMap_ = NULL;
    }
    NumPDEEqns_ = NumPDEEqns;
      
  }

  return( GetVbrMatrix() );
    
}

// ================================================ ====== ==== ==== == =
Epetra_VbrMatrix * Trilinos_Util_MatrixGallery::GetVbrMatrix() 
{
    
  if( VbrMatrix_ == NULL ) CreateVbrMatrix();

  return VbrMatrix_;
    
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::ComputeResidual(double & residual)
{

  // create solution and rhs if needed (matrix_ and ExactSolution are
  //  created by CreateRHS is necessary)
  if( rhs_ == NULL ) CreateRHS();

  Epetra_Vector Ax(*map_);
  assert(matrix_->Multiply(false, *ExactSolution_, Ax)==0);

  assert(Ax.Update(1.0, *rhs_, -1.0)==0);

  assert(Ax.Norm2(&residual)==0);

  return 0;
}
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::ComputeDiffBetweenStartingAndExactSolutions(double & residual)
{

  // create solution and rhs if needed (matrix_ and ExactSolution are
  //  created by CreateRHS is necessary)
  if( rhs_ == NULL ) CreateRHS();

  Epetra_Vector temp(*map_);

  assert(temp.Update(1.0, *ExactSolution_, -1.0, *StartingSolution_, 0.0)==0);

  assert(temp.Norm2(&residual)==0);

  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMap() 
{

  Epetra_Time Time(*comm_);
  
  if( verbose_ == true ) {
    cout << OutputMsg << "Creating Map...\n";
  }
  
  if( name_ == "diag" || name_ == "tridiag"  ||
      name_ == "laplace 1d" || name_ == "eye" ||
      name_ == "lehmer" || name_ == "minij" ||
      name_ == "ris" || name_ == "hilbert" ||
      name_ == "jordblock" ) {
    if( NumGlobalElements_ <= 0 ) {
      if( nx_ > 0 ) NumGlobalElements_ = nx_;
      else {
	cerr << ErrorMsg << "problem size not correct\n";
	return -1;
      }
    }
  }
    
  else if( name_ == "laplace_2d" || name_ == "cross_stencil_2d" ) {
    if( NumGlobalElements_ <= 0 ) {  
      if( nx_ > 0 && ny_ > 0 )
	NumGlobalElements_ = nx_*ny_;
      else {
	cerr << ErrorMsg << "problem size not correct\n";
	return -1;
      }
    }
  }
    
  else if( name_ == "laplace_3d" || name_ == "cross_stencil_3d" ) {
    if( NumGlobalElements_ <= 0 ) {
      if( nx_ > 0 && ny_ > 0 && nz_ > 0 )
	NumGlobalElements_ = nx_*ny_*nz_;
      else {
	cerr << ErrorMsg << "problem size not correct\n";
	return -1;
      }
    }
  } else {

    cerr << ErrorMsg << "matrix name is incorrect or not set ("
	 << name_ << ")\n";
    exit( EXIT_FAILURE );

  }
    
  if( MapType_ == "linear" ) {
      
    map_ = new Epetra_Map (NumGlobalElements_,0,*comm_);
      
  } else if( MapType_ == "box" ) {

    if( mx_ == -1 || my_ == -1 ) {
      mx_ = (int)sqrt((double)comm_->NumProc());
      my_ = mx_;
      if( mx_ * my_ != comm_->NumProc() ) {
	cerr << ErrorMsg << "number of processes must be square number\n"
	     << ErrorMsg << "otherwise set mx and my\n";
	return -1;
      }
    }

    if( nx_ == -1 || ny_ == -1 ) {
      nx_ = (int)sqrt((double)NumGlobalElements_);
      ny_ = mx_;
      if( mx_ * my_ != NumGlobalElements_ ) {
	cerr << ErrorMsg << "number of processes must be square number\n"
	     << ErrorMsg << "otherwise set mx and my\n";
	return -1;
      }
    }

    // how to divide the axis

    int modx = (nx_+(nx_%mx_))/mx_;
    int mody = (ny_+(ny_%my_))/my_;

    cout << modx << " " << mody << endl;
      
    int MyPID = comm_->MyPID(), startx, starty, endx, endy;
    int xpid = MyPID/mx_;
    int ypid = MyPID%my_;

    cout << xpid << " " << ypid << endl;

    startx = xpid*modx;
    if( (xpid+1)*modx < nx_ ) endx = (xpid+1)*modx;
    else endx = nx_;
    starty = ypid*mody;
    if( (ypid+1)*mody < ny_ ) endy = (ypid+1)*mody;
    else endy = ny_;

    cout << startx << " - " << endx << " _ " << starty << " _ " << endy << endl;
      
    int NumMyElements = (endx-startx)*(endy-starty);
    int * MyGlobalElements = new int[NumMyElements];
    int count = 0;
      
    for( int i=startx ; i<endx ; ++i ) {
      for( int j=starty ; j<endy ; ++j ) {
	MyGlobalElements[count++] = i+j*nx_;
      }
    }
      
    map_ = new Epetra_Map (NumGlobalElements_,NumMyElements,MyGlobalElements,0,*comm_);

    // I delete this guy to that this case is not different from the
    // others, and I don't have to clean up this mess while
    // destroying the object.
      
    delete MyGlobalElements;
      
  } else if( MapType_ == "interlaced" ) {

    int NumProcs = comm_->NumProc();
    int MyPID = comm_->MyPID();
      
    int NumMyElements = NumGlobalElements_/NumProcs;
    if( MyPID < NumGlobalElements_%NumProcs ) NumMyElements++;
      
    int count = 0;
    int * MyGlobalElements = new int[NumMyElements];
      
    for( int i=0 ; i<NumGlobalElements_ ; ++i ) {
      if( i%NumProcs == MyPID ) 
	MyGlobalElements[count++] = i;
    }

    if( count != NumMyElements ) {
      cerr << ErrorMsg << "something went wrong in CreateMap\n";
      cerr << ErrorMsg << "count = " << count << ", NumMyElements = "
	   << NumMyElements << endl;
      exit( EXIT_FAILURE);
    }
    
    map_ = new Epetra_Map (NumGlobalElements_,NumMyElements,MyGlobalElements,0,*comm_);
    delete MyGlobalElements;

  } else {
    cerr << ErrorMsg << "MapType has an incorrect value (" << MapType_ << ")\n";
    return -1;
  }
    
  // local number of rows
  NumMyElements_ = map_->NumMyElements();
  // get update list
  MyGlobalElements_ = map_->MyGlobalElements( );

  if( verbose_ == true ) {
    cout << OutputMsg << "Time to create Map: "
	 << Time.ElapsedTime() << " (s)\n";
  }
  
  return 0;
}
  
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrix() 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating Matrix...\n";
  }
  
  // HB matrices are different, as their dimension has to be read before.
  // Here the idea is to read the matrix on proc 0, then build the
  // map, then redistribute it linearly

  if( name_ == "hb" ) {

    Epetra_Time Time(*comm_);
    ReadHBMatrix();
    if( verbose_ == true ) {
      cout << OutputMsg << "Time to create matrix: "
	   << Time.ElapsedTime() << " (s)\n";
    }
    
  } else {
      
    if( map_ == NULL ) CreateMap();     

    Epetra_Time Time(*comm_);
    
    if( name_ == "diag" ) CreateMatrixDiag();

    else if( name_ == "eye" ) CreateEye();
      
    else if( name_ == "tridiag" ) CreateMatrixTriDiag();
      
    else if( name_ == "laplace_1d" ) CreateMatrixLaplace1d();
      
    else if( name_ == "laplace_2d" ) CreateMatrixLaplace2d();
      
    else if( name_ == "laplace_3d" ) CreateMatrixLaplace3d();
      
    else if( name_ == "cross_stencil_2d" ) CreateMatrixCrossStencil2d();
      
    else if( name_ == "cross_stencil_3d" ) CreateMatrixCrossStencil3d();

    else if( name_ == "lehmer" ) CreateMatrixLehmer();

    else if( name_ == "minij" ) CreateMatrixMinij();

    else if( name_ == "ris" ) CreateMatrixRis();

    else if( name_ == "hilbert" ) CreateMatrixHilbert();

    else if( name_ == "jordblock" ) CreateMatrixJordblock();
    
    else {
      cerr << ErrorMsg << "matrix name is incorrect or not set ("
	   << name_ << ")\n";
      exit( EXIT_FAILURE );
    }

    if( verbose_ == true ) {
      cout << OutputMsg << "Time to create matrix: "
	   << Time.ElapsedTime() << " (s)\n";
    }
  }

  return 0;    
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util_MatrixGallery::CreateVbrExactSolution(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating exact solution (VBR)...\n";
  }

  // check if already have one; in this case delete and rebuild
  if( VbrExactSolution_ != NULL ) delete VbrExactSolution_;
  // need an exact solution for Crs first
  if( ExactSolution_ == NULL ) CreateExactSolution();
  // need a block map first
  if( BlockMap_ == NULL ) CreateBlockMap();
  // now we can expand to the Vbr format
  VbrExactSolution_ = new Epetra_Vector(*BlockMap_);
  for( int j=0 ; j<NumMyElements_ ; j++ ) 
    for( int i=0 ; i<NumPDEEqns_ ; ++i ) {
      (*VbrExactSolution_)[j*NumPDEEqns_+i] = (*ExactSolution_)[j];
    }
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util_MatrixGallery::CreateExactSolution() 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating exact solution `"
	 << ExactSolutionType_ << "'...\n";
  }
  
  if( map_ == NULL ) CreateMap();

  if( ExactSolution_ == NULL ) {
    ExactSolution_ = new Epetra_Vector(*map_);
    if( ExactSolutionType_ == "random" ) {
      ExactSolution_->Random();
    } else if( ExactSolutionType_ == "constant" ) {
      ExactSolution_->PutScalar(1.0);
    } else if( ExactSolutionType_ == "linear" ) {
      for( int i=0 ; i<NumMyElements_ ; i++ ) 
	(*ExactSolution_)[i] = alpha_*MyGlobalElements_[i];
    } else {
      cerr << ErrorMsg << "exact solution type is not correct : "
	   << ExactSolutionType_ << endl;
    }
  }
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util_MatrixGallery::CreateStartingSolution() 
{
  if( verbose_ == true ) {
    cout << OutputMsg << "Creating starting solution `"
	 << StartingSolutionType_ << "'...\n";
  }

  if( map_ == NULL ) CreateMap();

  if( StartingSolution_ == NULL ) {
    StartingSolution_ = new Epetra_Vector(*map_);
    if( StartingSolutionType_ == "random" ) {
      StartingSolution_->Random();
    } else if( StartingSolutionType_ == "zero" ) {
      StartingSolution_->PutScalar(1.0);
    } else {
      cerr << ErrorMsg << "starting solution type is not correct : "
	   << StartingSolutionType_ << endl;
    }
  }
}
  
// ================================================ ====== ==== ==== == =
void Trilinos_Util_MatrixGallery::CreateRHS() 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating RHS...\n";
  }

  if( map_ == NULL ) CreateMap();
  if( matrix_ == NULL ) CreateMatrix();
  if( ExactSolution_ == NULL )  CreateExactSolution();

  if( rhs_ == NULL ) {
    rhs_ = new Epetra_Vector(*map_);
    matrix_->Multiply(false,*ExactSolution_,*rhs_);
  }
    
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util_MatrixGallery::CreateVbrRHS() 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating RHS (VBR)...\n";
  }

  if( VbrRhs_ != NULL ) {
    delete VbrRhs_;
    VbrRhs_ = NULL;
  }
    
  // need a rhs for crs
  if( rhs_ == NULL ) CreateRHS();
  // need a block map based on map_
  if( BlockMap_ == NULL ) CreateBlockMap();
  // require VbrMatrix to be formed first
  if( VbrMatrix_ == NULL ) CreateVbrMatrix();
  // also need an exact solution
  if( VbrExactSolution_ == NULL )  CreateVbrExactSolution();

  VbrRhs_ = new Epetra_Vector( *BlockMap_);
  VbrMatrix_->Multiply(false,*VbrExactSolution_,*VbrRhs_);
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateEye() 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `eye'...\n";
  }
  
  a_ = 1.0;
  CreateMatrixDiag();
  return 0;    
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixDiag() 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `diag'...\n";
  }
  
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,1);
  double Value;
    
  for( int i=0 ; i<NumMyElements_; ++i ) {

    int Indices = MyGlobalElements_[i];
    if( VectorA_ == NULL ) Value = a_;
    else Value = (*VectorA_)[i];
      
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], 1, &Value, &Indices)==0);
      
  }
    
  assert(matrix_->TransformToLocal()==0);

  return 0;
    
}
  
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixTriDiag() 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `tridiag'...n";
  }
  
  int ierr;

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,3);

  double *Values = new double[2];
  Values[0] = b_; Values[1] = c_;
  int *Indices = new int[2];
  int NumEntries;

  for( int i=0 ; i<NumMyElements_; ++i ) {
    if (MyGlobalElements_[i]==0) {
      Indices[0] = 1;
      NumEntries = 1;
      if( VectorC_ == NULL ) Values[0] = c_;
      else Values[0] = (*VectorC_)[i];
    } else if (MyGlobalElements_[i] == NumGlobalElements_-1) {
      Indices[0] = NumGlobalElements_-2;
      NumEntries = 1;
      if( VectorC_ == NULL ) Values[0] = c_;
      else Values[0] = (*VectorC_)[i];
    } else {
      Indices[0] = MyGlobalElements_[i]-1;
      Indices[1] = MyGlobalElements_[i]+1;
      NumEntries = 2;
    }
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], NumEntries, Values, Indices)==0);
    // Put in the diagonal entry
    if( VectorA_ == NULL ) Values[0] = a_;
    else Values[0] = (*VectorA_)[i];
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], 1, Values, MyGlobalElements_+i)==0);
  }
  
  // Finish up, trasforming the matrix entries into local numbering,
  // to optimize data transfert during matrix-vector products
  assert(matrix_->TransformToLocal()==0);

  delete Values;
  delete Indices;

  return 0;
    
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixLaplace1d() 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `laplace_1d'...\n";
  }

  a_ = 2.0;
  b_ = -1.0;
  c_ = -1.0;

  CreateMatrixTriDiag();
  
  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixCrossStencil2d() 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `cross_stencil_2d'...\n";
  }

  if( nx_ == -1 || ny_ == -1 ) {
    nx_ = (int)sqrt((double)NumGlobalElements_);
    ny_ = nx_;
      
    if( nx_ * ny_ != NumGlobalElements_ ) {
      cerr << ErrorMsg << "The number of global elements must be a square number\n"
	   << ErrorMsg << "(now is " << NumGlobalElements_ << "). Returning...\n";
    }
  }
    
  int left, right, lower, upper;
    
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,5);
    
  // Add  rows one-at-a-time
    
  double Values[4], diag;
  int Indices[4];
  int NumEntries;

  //    e
  //  b a c
  //    d
  for( int i=0 ; i<NumMyElements_; ++i ) {
    int NumEntries=0;
    GetNeighboursCartesian2d(  MyGlobalElements_[i], nx_, ny_, 
			       left, right, lower, upper);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      if( VectorB_ == NULL ) Values[NumEntries] = b_;
      else Values[NumEntries] = (*VectorB_)[i];
      ++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      if( VectorC_ == NULL ) Values[NumEntries] = c_;
      else Values[NumEntries] = (*VectorC_)[i];
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      if( VectorD_ == NULL ) Values[NumEntries] = d_;
      else Values[NumEntries] = (*VectorD_)[i];
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      if( VectorE_ == NULL ) Values[NumEntries] = e_;
      else Values[NumEntries] = (*VectorE_)[i];
      ++NumEntries;
    }
    // put the off-diagonal entries
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], NumEntries, 
				       Values, Indices)==0);
    // Put in the diagonal entry
    if( VectorA_ == NULL ) diag = a_;
    else diag = (*VectorA_)[i];
	
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], 1, 
				       &diag, MyGlobalElements_+i)==0);
  }
  matrix_->TransformToLocal();

  return 0;
      
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixLaplace2d() 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `laplace_2d'...\n";
  }

  a_ = 4.0;
  b_ = -1.0;
  c_ = -1.0;
  d_ = -1.0;
  e_ = -1.0;

  CreateMatrixCrossStencil2d();

  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixLaplace3d()
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `laplace_3d'...\n";
  }
  
  a_ = 6.0;
  b_ = -1.0;
  c_ = -1.0;
  d_ = -1.0;
  e_ = -1.0;
  f_ = -1.0;
  g_ = -1.0;

  CreateMatrixCrossStencil3d();

  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixCrossStencil3d() 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `cross_stencil_3d'...\n";
  }
  
  if( nx_ == -1 || ny_ == -1 || nz_ == -1 ) {
    nx_ = (int)pow((double)NumGlobalElements_,0.333333333);
    ny_ = nx_;
    nz_ = nz_;
      
    if( nx_ * ny_ *nz_ != NumGlobalElements_ ) {
      cerr << ErrorMsg << "The number of global elements must be a cubical number\n"
	   << ErrorMsg << "(now is " << NumGlobalElements_ << "). Returning...\n";
    }
  }
    
  int left, right, lower, upper, below, above;
    
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,7);
    
  // Add  rows one-at-a-time
    
  double Values[6], diag;
  int Indices[6];
  int NumEntries;

  //    e 
  //  b a c
  //    d
  // + f below and g above
    
  for( int i=0 ; i<NumMyElements_; ++i ) {
    int NumEntries=0;
    GetNeighboursCartesian3d(  MyGlobalElements_[i], nx_, ny_, nz_,
			       left, right, lower, upper, below, above);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      if( VectorB_ == NULL ) Values[NumEntries] = b_;
      else Values[NumEntries] = (*VectorB_)[i];
      ++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      if( VectorC_ == NULL ) Values[NumEntries] = c_;
      else Values[NumEntries] = (*VectorC_)[i];
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      if( VectorD_ == NULL ) Values[NumEntries] = d_;
      else Values[NumEntries] = (*VectorD_)[i];
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      if( VectorE_ == NULL ) Values[NumEntries] = e_;
      else Values[NumEntries] = (*VectorE_)[i];
      ++NumEntries;
    }
    if( below != -1 ) {
      Indices[NumEntries] = below;
      if( VectorF_ == NULL ) Values[NumEntries] = f_;
      else Values[NumEntries] = (*VectorF_)[i];
      ++NumEntries;
    }
    if( above != -1 ) {
      Indices[NumEntries] = above;
      if( VectorG_ == NULL ) Values[NumEntries] = g_;
      else Values[NumEntries] = (*VectorG_)[i];
      ++NumEntries;
    }
    // put the off-diagonal entries
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], NumEntries, 
				       Values, Indices)==0);
    // Put in the diagonal entry
    if( VectorA_ == NULL ) diag = a_;
    else diag = (*VectorA_)[i];
	
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], 1, 
				       &diag, MyGlobalElements_+i)==0);
  }

  matrix_->TransformToLocal();
  return 0;
      
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixLehmer() 
{

    if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `lehmer'...\n";
  }

  // this is actually a dense matrix, stored into Crs format
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int * Indices = new int[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  for( int i=0 ; i<NumGlobalElements_ ; ++i ) Indices[i] = i;
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
      if( i>=j ) Values[j] = 1.0*(j+1)/(i+1);
      else       Values[j] = 1.0*(i+1)/(j+1);
    }
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumGlobalElements_, Values, Indices)==0);
      
  }
    
  delete Indices;
  delete Values;

  assert(matrix_->TransformToLocal()==0);
  
  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixMinij() 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `minij'...\n";
  }
  
  // this is actually a dense matrix, stored into Crs format
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int * Indices = new int[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  for( int i=0 ; i<NumGlobalElements_ ; ++i ) Indices[i] = i;
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
      if( i>=j ) Values[j] = 1.0*(j+1);
      else       Values[j] = 1.0*(i+1);
    }
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumGlobalElements_, Values, Indices)==0);
      
  }
    
  delete Indices;
  delete Values;

  assert(matrix_->TransformToLocal()==0);
  
  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixRis() 
{
  
  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `ris'...\n";
  }

  // this is actually a dense matrix, stored into Crs format
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int * Indices = new int[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  for( int i=0 ; i<NumGlobalElements_ ; ++i ) Indices[i] = i;
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
      Values[j] = 0.5/(NumGlobalElements_ -(i+1)-(j+1)+1.5);
      
    }
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumGlobalElements_, Values, Indices)==0);
      
  }
    
  delete Indices;
  delete Values;

  assert(matrix_->TransformToLocal()==0);
  
  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixHilbert() 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `hilbert'...\n";
  }

  // this is actually a dense matrix, stored into Crs format
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int * Indices = new int[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  for( int i=0 ; i<NumGlobalElements_ ; ++i ) Indices[i] = i;
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
      Values[j] = 1.0/((i+1)+(j+1)-1);
    }
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumGlobalElements_, Values, Indices)==0);
      
  }
    
  delete Indices;
  delete Values;

  assert(matrix_->TransformToLocal()==0);
  
  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixJordblock() 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `jordblock'...\n";
  }

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,2);

  int Indices[2];
  double Values[2];
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int NumEntries = 0;
    if( MyGlobalElements_[i] != NumGlobalElements_-1 ) {
      Indices[NumEntries] = MyGlobalElements_[i]+1;
      Values[NumEntries] = 1.0;
      NumEntries++;
    }
    // diagonal contribution
    Indices[NumEntries] = MyGlobalElements_[i];
    if( VectorA_ != NULL ) Values[NumEntries] = (*VectorA_)[i];
    else                   Values[NumEntries] = a_;
    NumEntries++;

    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumEntries, Values, Indices)==0);
      
  }
    
  assert(matrix_->TransformToLocal()==0);
  
  return 0;
}
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::ReadHBMatrix() 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Reading HB matrix `"
	 << FileName_ << "'...\n";
  }

  Epetra_Map * readMap;
  Epetra_CrsMatrix * readA; 
  Epetra_Vector * readx; 
  Epetra_Vector * readb;
  Epetra_Vector * readxexact;
    
  // Call routine to read in HB problem
    
  Trilinos_Util_ReadHb2Epetra((char*)FileName_.c_str(), *comm_, readMap, readA, readx, 
			      readb, readxexact);
    
  NumGlobalElements_ = readMap->NumGlobalElements();

  // Create uniform distributed map
  map_ = new Epetra_Map(NumGlobalElements_, 0, *comm_);

  // Create Exporter to distribute read-in matrix and vectors
  Epetra_Export exporter(*readMap, *map_);
    
  matrix_ = new Epetra_CrsMatrix(Copy, *map_, 0);
  StartingSolution_ = new Epetra_Vector(*map_);
  rhs_ = new Epetra_Vector(*map_);
  ExactSolution_ = new Epetra_Vector(*map_);

  StartingSolution_->Export(*readx, exporter, Add);
  rhs_->Export(*readb, exporter, Add);
  ExactSolution_->Export(*readxexact, exporter, Add);
  matrix_->Export(*readA, exporter, Add);

  assert(matrix_->TransformToLocal()==0);    

  delete readA;
  delete readx;
  delete readb;
  delete readxexact;
  delete readMap;
  
  // local number of rows
  NumMyElements_ = map_->NumMyElements();
  // get update list
  MyGlobalElements_ = map_->MyGlobalElements( );
    
  return 0;

}

// ================================================ ====== ==== ==== == =
void Trilinos_Util_MatrixGallery::CreateBlockMap(void) 
{
        
  if( verbose_ == true ) {
    cout << OutputMsg << "Creating BlockMap...\n";
  }

  // dimension of each block
  Epetra_IntSerialDenseVector ElementSizeList(NumMyElements_);

  MaxBlkSize_ = NumPDEEqns_;
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    ElementSizeList[i] = NumPDEEqns_;
    //      if( ElementSizeList[i] > MaxBlkSize_ ) MaxBlkSize_ = ElementSizeList[i];
  }
  
  BlockMap_ = new Epetra_BlockMap(NumGlobalElements_,NumMyElements_,
				  MyGlobalElements_, 
				  ElementSizeList.Values(),0,*comm_);

}

// ================================================ ====== ==== ==== == =
void Trilinos_Util_MatrixGallery::CreateVbrMatrix(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating VBR matrix...\n";
  }

  if( map_ == NULL ) CreateMap();
  if( matrix_ == NULL ) CreateMatrix();
  if( BlockMap_ == NULL ) CreateBlockMap();

  int MaxNnzPerRow = matrix_->MaxNumEntries();
  if( MaxNnzPerRow == 0 ) {
    cerr << ErrorMsg << "something went wrong in `CreateMatrix'\n"
	 << ErrorMsg << "MaxNnzPerRow == 0 \n";
    exit( EXIT_FAILURE );
  }
    
  // create a VBR matrix based on BlockMap
  VbrMatrix_ = new Epetra_VbrMatrix(Copy, *BlockMap_,MaxNnzPerRow);

  // size of each VBR block
  int MaxBlockSize = MaxBlkSize_*MaxBlkSize_;

  int CrsNumEntries;
  int * CrsIndices;
  double * CrsValues;
    
  int * Indices = new int[MaxNnzPerRow];
  double * Values = new double[MaxBlockSize];
  int BlockRows = NumPDEEqns_;
  int ierr;
    
  cout << MaxNnzPerRow << " nnz per row\n";
    
  // cycle over all the local rows. 
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    
    // get GID of local row
    int GlobalNode = MyGlobalElements_[i];
    // extract Crs row
    ierr = matrix_->ExtractMyRowView(i,CrsNumEntries,
				     CrsValues,CrsIndices);

    cout << ierr;
      
    // with VBR matrices, we have to insert one block at time.
    // This required two more instructions, one to start this
    // process (BeginInsertGlobalValues), and another one to
    // commit the end of submissions (EndSubmitEntries).
    
    VbrMatrix_->BeginInsertGlobalValues(GlobalNode, CrsNumEntries, CrsIndices);

    cout << "---> " << CrsNumEntries << "   " << *CrsIndices << endl;
      
    for( int i=0 ; i<CrsNumEntries ; ++i ) {
	
      for( int k=0 ; k<BlockRows ; ++k ) { // rows
	for( int h=0 ; h<BlockRows ; ++h ) { // cols
	  if( k == h ) Values[k+h*BlockRows] = CrsValues[i];
	  else Values[k+h*BlockRows] = 0.0;
	}
      }
	  
      VbrMatrix_->SubmitBlockEntry(Values,BlockRows,BlockRows,BlockRows);
	
    }

    VbrMatrix_->EndSubmitEntries();
  }
    
  delete Indices;
  delete Values;

  VbrMatrix_->TransformToLocal();
  
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::ComputeResidualVbr(double & residual)
{

  // create solution and rhs if needed (matrix_ and ExactSolution are
  //  created by CreateRHS is necessary)
  if( VbrRhs_ == NULL ) CreateVbrRHS();

  Epetra_Vector Ax(*BlockMap_);
  assert(VbrMatrix_->Multiply(false, *VbrExactSolution_, Ax)==0);

  assert(Ax.Update(1.0, *VbrRhs_, -1.0)==0);

  assert(Ax.Norm2(&residual)==0);

  return 0;
}
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::ComputeDiffBetweenStartingAndExactSolutionsVbr(double & residual)
{

  // create solution and rhs if needed (matrix_ and ExactSolution are
  //  created by CreateRHS is necessary)
  if( VbrRhs_ == NULL ) CreateVbrRHS();

  Epetra_Vector temp(*BlockMap_);

  assert(temp.Update(1.0, *VbrExactSolution_, -1.0, *VbrStartingSolution_, 0.0)==0);

  assert(temp.Norm2(&residual)==0);

  return 0;
}

// ================================================ ====== ==== ==== == =  
void Trilinos_Util_MatrixGallery::GetNeighboursCartesian2d( const int i, const int nx, const int ny,
				int & left, int & right, 
				int & lower, int & upper) 
{

  int ix, iy;
  ix = i%nx;
  iy = (i - ix)/nx;

  if( ix == 0 ) 
    left = -1;
  else 
    left = i-1;
  if( ix == nx-1 ) 
    right = -1;
  else
    right = i+1;
  if( iy == 0 ) 
    lower = -1;
  else
    lower = i-nx;
  if( iy == ny-1 ) 
    upper = -1;
  else
    upper = i+nx;

  return;

}

// ================================================ ====== ==== ==== == =
void Trilinos_Util_MatrixGallery::GetNeighboursCartesian3d( const int i, const int nx, const int ny, const int nz,
				int & left, int & right, int & lower, int & upper,
				int & below, int & above ) 
{

  int ixy, iz;
  ixy = i%(nx*ny);
    
  iz = (i - ixy)/(nx*ny);

  if( iz == 0 ) 
    below = -1;
  else 
    below = i-nx*ny;
  if( iz == nz-1 ) 
    above = -1;
  else
    above = i+nx*ny;

  GetNeighboursCartesian2d( ixy, nx_, ny_, left, right, lower, upper);
    
  return;

}

// ================================================ ====== ==== ==== == =
void Trilinos_Util_MatrixGallery::ZeroOutData() 
{
  NumGlobalElements_ = -1;
  nx_ = -1;    ny_ = -1;     nz_ = -1;
  mx_ = -1;    mx_ = -1;     mz_ = -1;
    
  a_ = 1.0, b_ = 0.0, c_ = 0.0, d_ = 0.0, e_ = 0.0, f_ = 0.0, g_ = 0.0;
  alpha_ = 1.0;
  beta_ = 0.0;
  gamma_ = 0.0;
  delta_ = 0.0;
    
  VectorA_ = NULL;
  VectorB_ = NULL;
  VectorC_ = NULL;
  VectorD_ = NULL;
  VectorE_ = NULL;
  VectorF_ = NULL;
  VectorG_ = NULL;
    
  map_ = NULL;
  matrix_ = NULL;
  BlockMap_ = NULL;
  VbrMatrix_ = NULL;
  ExactSolution_ = NULL;
  StartingSolution_ = NULL;
  rhs_ = NULL;

  VbrExactSolution_ = NULL;
  VbrStartingSolution_ = NULL;
  VbrRhs_ = NULL;
  
  MapType_ = "linear";
  ExactSolutionType_ = "constant";
  StartingSolutionType_ = "zero";
    
  NumPDEEqns_= 1;
    
}

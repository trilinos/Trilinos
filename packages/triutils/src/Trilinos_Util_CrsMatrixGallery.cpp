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
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "Epetra_Util.h"
#include "Trilinos_Util.h"
#include <string>

#include "Trilinos_Util_CommandLineParser.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

const double UNDEF = -99999.87;
const bool Scaling = false;

// ================================================ ====== ==== ==== == =
Trilinos_Util::CrsMatrixGallery::CrsMatrixGallery(const string name, 
							       const Epetra_Comm & comm ) :
  comm_(&comm), name_(name)
{
  ZeroOutData();
  // verbosity level
  if( comm_->MyPID()==0 ) verbose_ = true;
  else verbose_ = false;
  // fix error message
  ErrorMsg = "ERROR [CrsMatrixGallery]: ";
  OutputMsg = "CrsMatrixGallery: ";
  
}

// ================================================ ====== ==== ==== == =
Trilinos_Util::CrsMatrixGallery::CrsMatrixGallery(const string name, 
							       const Epetra_Map & map ) :
  comm_(&(map.Comm())), name_(name)
{
  ZeroOutData();
  // verbosity level
  if( comm_->MyPID()==0 ) verbose_ = true;
  else verbose_ = false;
  // fix error message
  ErrorMsg = "ERROR [Trilinos_Util::CrsMatrixGallery]: ";
  OutputMsg = "Trilinos_Util::CrsMatrixGallery: ";
  
  map_ = new Epetra_Map(map);
  NumGlobalElements_ = map_->NumGlobalElements();
  NumMyElements_ = map_->NumMyElements();
  MyGlobalElements_ = map_->MyGlobalElements( );
    
}
  
// ================================================ ====== ==== ==== == =
Trilinos_Util::CrsMatrixGallery::~CrsMatrixGallery(void) 
{

  // linear problem
  if( LinearProblem_ != NULL ) delete LinearProblem_;

  // Crs data
  if( matrix_ != NULL ) delete matrix_;
  if( ExactSolution_ != NULL ) delete ExactSolution_;
  if( StartingSolution_ != NULL ) delete StartingSolution_;
  if( rhs_ != NULL ) delete rhs_;
  if( map_ != NULL ) delete map_;

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
int Trilinos_Util::CrsMatrixGallery::Set(const string parameter, const int value)
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
int Trilinos_Util::CrsMatrixGallery::Set(const string parameter, const string value )
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
  else if( parameter == "rhs_type" ) {
    RhsType_ = value;
  }
  else if( parameter == "output" ) {
    if( value == "none" ) verbose_ = false;
    if( value == "proc 0" ) {
      if( comm_->MyPID()==0 ) verbose_ = true;
      else verbose_ = false;
    } else {
      verbose_ = true;
    }
  } else if( parameter == "expand_type" ) {
    ExpandType_ = value;
  } else {
    cerr << ErrorMsg << "wrong input parameter (" << parameter << ")\n";
    return -1;
  }

  return 0;
}
  
// ================================================ ====== ==== ==== == =
int Trilinos_Util::CrsMatrixGallery::Set(const string parameter, const double value)
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
  } else if( parameter == "conv" ) {
    conv_ = value;
    return 0;
  } else if( parameter == "diff" ) {
    diff_ = value;
    return 0;
  } else if( parameter == "source" ) {
    source_ = value;
    return 0;
  } else if( parameter == "alpha" ) {
    alpha_ = value;
    return 0;
  }

  cerr << ErrorMsg << "input string not valid\n";
  return -2;
}
  
// ================================================ ====== ==== ==== == =
int Trilinos_Util::CrsMatrixGallery::Set(const string parameter, const Epetra_Vector & value)
{

  if( value.Map().SameAs(*map_) == false ) {
    cerr << ErrorMsg << "input vector must have the same map used to\n"
	 << ErrorMsg << "create the Trilinos_Util::CrsMatrixGallery object. Continuing\n";
    return -2;
  }
    
  if( parameter == "a" ) {
    VectorA_ = new Epetra_Vector(value);
  } else if( parameter == "b" ) {
    VectorB_ = new Epetra_Vector(value);
  }
  else if( parameter == "c" ) {
    VectorC_ = new Epetra_Vector(value);
  }
  else if( parameter == "d" ) {
    VectorD_ = new Epetra_Vector(value);
  }
  else if( parameter == "e" ) {
    VectorE_ = new Epetra_Vector(value);
  }
  else if( parameter == "f" ) {
    VectorF_ = new Epetra_Vector(value);
  }
  else if( parameter == "g" ) {
    VectorG_ = new Epetra_Vector(value);
  } else {
    cerr << ErrorMsg << "input string not valid\n";
    return -3;
  }

  return 0;
}

// ================================================ ====== ==== ==== == =  

int Trilinos_Util::CrsMatrixGallery::Set(Trilinos_Util::CommandLineParser & CLP)
{
  int count;
  
  string Options[15];

  // all options with strings
  count = 0;
  Options[count++] = "problem_type";
  Options[count++] = "map_type";
  Options[count++] = "exact_solution";
  Options[count++] = "matrix_name";
  Options[count++] = "starting_solution";
  Options[count++] = "output";
  Options[count++] = "expand_type";
  Options[count++] = "rhs_type";
  
  for( int i=0 ; i<count ; i++ ) {
    string parameter = "-"+Options[i];    
    if( CLP.Has(parameter) == true ) {
      string value = CLP.Get(parameter,"not-set");
      Set(Options[i],value);
      
    }
  }

  // all options with integers
  Options[0] = "problem_size";
  Options[1] = "nx";
  Options[2] = "ny";
  Options[3] = "nz";
  Options[4] = "mx";
  Options[5] = "my";
  Options[6] = "mz";
  Options[7] = "num_pde_eqns";

  for(  int i=0 ; i<8 ; i++ ) {
    string parameter = "-"+Options[i];   
    if( CLP.Has(parameter) == true ) {
      Set(Options[i],CLP.Get(parameter,(int)1));
    }
  }
  
  // all options with doubles
  Options[0] = "a";
  Options[1] = "b";
  Options[2] = "c";
  Options[3] = "d";
  Options[4] = "e";
  Options[5] = "f";
  Options[6] = "g";
  Options[7] = "conv";
  Options[8] = "diff";
  Options[9] = "source";
  Options[10] = "alpha";
  
  for( int i=0 ; i<11 ; i++ ) {
    string parameter = "-"+Options[i];
    
    if( CLP.Has(parameter) == true ) {
      Set(Options[i],CLP.Get(parameter,1.0));
    }
    
  }

  return 0;
}

// ================================================ ====== ==== ==== == =  
Epetra_CrsMatrix * Trilinos_Util::CrsMatrixGallery::GetMatrix(void) 
{
  if( matrix_ == NULL ) CreateMatrix();
  return( matrix_ );
}

// ================================================ ====== ==== ==== == =  
Epetra_CrsMatrix & Trilinos_Util::CrsMatrixGallery::GetMatrixRef(void) 
{
  if( matrix_ == NULL ) CreateMatrix();
  return( *matrix_ );
}
  
// ================================================ ====== ==== ==== == =
Epetra_Vector * Trilinos_Util::CrsMatrixGallery::GetExactSolution(void) 
{
  if( ExactSolution_ == NULL ) CreateExactSolution();
  return ExactSolution_;
}

// ================================================ ====== ==== ==== == =
Epetra_Vector * Trilinos_Util::CrsMatrixGallery::GetStartingSolution(void)
{
  if( StartingSolution_ == NULL ) CreateStartingSolution();
  return StartingSolution_;
}

// ================================================ ====== ==== ==== == =
Epetra_Vector * Trilinos_Util::CrsMatrixGallery::GetRHS(void)
{
  if( rhs_ == NULL ) CreateRHS();
  return rhs_;
}

// ================================================ ====== ==== ==== == =
const Epetra_Map * Trilinos_Util::CrsMatrixGallery::GetMap(void)
{
  if( map_ == NULL ) CreateMap();
    
  return map_;
}

// ================================================ ====== ==== ==== == =
const Epetra_Map & Trilinos_Util::CrsMatrixGallery::GetMapRef(void)
{
  if( map_ == NULL ) CreateMap();
    
  return *map_;
}

// ================================================ ====== ==== ==== == =
Epetra_LinearProblem * Trilinos_Util::CrsMatrixGallery::GetLinearProblem(void) 
{
  // pointers, not really needed
  Epetra_CrsMatrix * A;
  Epetra_Vector * RHS;
  Epetra_Vector * StartingSolution;

  A = GetMatrix();
  RHS = GetRHS();
  StartingSolution = GetStartingSolution();

  // create linear problem
  if( LinearProblem_ != NULL ) delete LinearProblem_;
  LinearProblem_ = new Epetra_LinearProblem(A,StartingSolution,RHS);

  return LinearProblem_;

}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::ComputeResidual(double & residual)
{

  // create solution and rhs if needed (matrix_ and ExactSolution_ are
  //  created by CreateRHS if necessary)
  if( rhs_ == NULL ) CreateRHS();

  Epetra_Vector Ax(*map_);

  assert(matrix_->Multiply(false, *StartingSolution_, Ax)==0);
  assert(Ax.Update(1.0, *rhs_, -1.0)==0);
  assert(Ax.Norm2(&residual)==0);
  
  return;
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::ComputeDiffBetweenStartingAndExactSolutions(double & residual)
{

  // create solution and rhs if needed (matrix_ and ExactSolution are
  //  created by CreateRHS is necessary)
  if( rhs_ == NULL ) CreateRHS();

  Epetra_Vector temp(*map_);

  assert(temp.Update(1.0, *ExactSolution_, -1.0, *StartingSolution_, 0.0)==0);

  assert(temp.Norm2(&residual)==0);

  return;
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMap(void)
{

  Epetra_Time Time(*comm_);

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating Map `" << MapType_ << "'...\n";
  }

  // first get the problem size. For some problems. the user can
  // specify the problem size using different parameters (e.g.,
  // nx and ny for a 2D Laplace problem). I need the internal
  // variable NumGlobalElements_ properly set before continuing.
  // NOTE: for HB problems, this value has already been set
  
  if( name_ == "diag" || name_ == "tridiag"  ||
      name_ == "laplace_1d" || name_ == "laplace_1d_n" ||
      name_ == "eye" ||
      name_ == "lehmer" || name_ == "minij" ||
      name_ == "ris" || name_ == "hilbert" ||
      name_ == "jordblock" || name_ == "cauchy" ||
      name_ == "fielder" || name_ == "hanowa" ||
      name_ == "kms" || name_ == "parter" ||
      name_ == "pei" || name_ == "ones" ||
      name_ == "vander" ) {
    if( NumGlobalElements_ <= 0 ) {
      if( nx_ > 0 ) NumGlobalElements_ = nx_;
      else {
	cerr << ErrorMsg << "problem size not correct (" << NumGlobalElements_
	     << ")\n";
	exit( EXIT_FAILURE );
      }
    }
  }
    
  else if( name_ == "laplace_2d" || name_ == "laplace_2d_n" 
	   || name_ == "cross_stencil_2d"
	   || name_ == "laplace_2d_9pt" || name_ == "recirc_2d"
	   || name_ == "uni_flow_2d" || name_ == "recirc_2d_divfree" ) {
    if( NumGlobalElements_ <= 0 ) {  
      if( nx_ > 0 && ny_ > 0 )
	NumGlobalElements_ = nx_*ny_;
      else {
	cerr << ErrorMsg << "Problem size not correct (" << NumGlobalElements_
	     << ")" << endl;
	cerr << ErrorMsg << "It should be a perfect square" << endl;
	exit( EXIT_FAILURE );
      }
    }
  }
    
  else if( name_ == "laplace_3d" || name_ == "cross_stencil_3d" ) {
    if( NumGlobalElements_ <= 0 ) {
      if( nx_ > 0 && ny_ > 0 && nz_ > 0 )
	NumGlobalElements_ = nx_*ny_*nz_;
      else {
	cerr << ErrorMsg << "Problem size not correct (" << NumGlobalElements_
	     << ")" << endl;
	cerr << ErrorMsg << "It should be a perfect cube" << endl;
	exit( EXIT_FAILURE );
      }
    }

  } else if( name_ == "hb" || name_ == "matrix_market" ||
	     name_ == "triples_sym" || name_ == "triples_nonsym" ) {
    // The global number of elements has been set in reading matrix
    if( NumGlobalElements_ <= 0 ) {
	cerr << ErrorMsg << "Problem size not correct (" << NumGlobalElements_
	     << ")" << endl;
	exit( EXIT_FAILURE );
    }
    
  } else {

    cerr << ErrorMsg << "matrix name is incorrect or not set ("
	 << name_ << ")\n";
    exit( EXIT_FAILURE );

  }

  // check out whether one is using only one proc or not.
  // If yes, creation of map is straightforward. Then return.
  
  if( comm_->NumProc() == 1 ) {

    map_ = new Epetra_Map(NumGlobalElements_,0,*comm_);

  } else {

    // Here below more than one processor.
  
    if( MapType_ == "linear" ) {
      
      map_ = new Epetra_Map (NumGlobalElements_,0,*comm_);
      
    } else if( MapType_ == "box" ) {

      if( mx_ == -1 || my_ == -1 ) {
	mx_ = (int)sqrt((double)(comm_->NumProc()));
	my_ = mx_;
	if( mx_ * my_ != comm_->NumProc() ) {
	  cerr << ErrorMsg << "number of processes must be perfect square\n"
	       << ErrorMsg << "otherwise set mx and my\n";
	  exit( EXIT_FAILURE );
	}
      }

      SetupCartesianGrid2D();

      // how to divide the axis
      
      int modx = (nx_+(nx_%mx_))/mx_;
      int mody = (ny_+(ny_%my_))/my_;
      
      int MyPID = comm_->MyPID(), startx, starty, endx, endy;
      int xpid = MyPID/mx_;
      int ypid = MyPID%my_;
      
      startx = xpid*modx;
      if( (xpid+1)*modx < nx_ ) endx = (xpid+1)*modx;
      else endx = nx_;
      starty = ypid*mody;
      if( (ypid+1)*mody < ny_ ) endy = (ypid+1)*mody;
      else endy = ny_;
      
      int NumMyElements = (endx-startx)*(endy-starty);
      int * MyGlobalElements = new int[NumMyElements];
      int count = 0;
      
      for( int i=startx ; i<endx ; ++i ) {
	for( int j=starty ; j<endy ; ++j ) {
	  MyGlobalElements[count++] = i+j*nx_;
	}
      }
      
      map_ = new Epetra_Map (NumGlobalElements_,NumMyElements,MyGlobalElements,0,*comm_);
      
      // I delete this guy so that this case is not different from the
      // others, and I don't have to clean up this mess while
      // destroying the object.
      
      delete [] MyGlobalElements;
      
    } else if( MapType_ == "cube" ) {

      if( mx_ == -1 || my_ == -1 || mz_ == -1 ) {
	mx_ = (int)pow((double)(comm_->NumProc()),0.333334);
	my_ = mx_;
	mz_ = mx_;
	
	if( mx_ * my_ * mz_ != comm_->NumProc() ) {
	  cerr << ErrorMsg << "number of processes must be perfect cube\n"
	       << ErrorMsg << "otherwise set mx, my, and mz\n";
	  exit( EXIT_FAILURE );
	}
      }

      SetupCartesianGrid3D();
      
      // how to divide the axis

      int modx = (nx_+(nx_%mx_))/mx_;
      int mody = (ny_+(ny_%my_))/my_;
      int modz = (nz_+(nz_%mz_))/mz_;
      
      int MyPID = comm_->MyPID(), startx, starty, startz, endx, endy, endz;
      int mxy  = mx_*my_;
      int zpid = MyPID/mxy;
      int xpid = (MyPID%mxy)/mx_;
      int ypid = (MyPID%mxy)%my_;

      startx = xpid*modx;
      if( (xpid+1)*modx < nx_ ) endx = (xpid+1)*modx;
      else endx = nx_;
      starty = ypid*mody;
      if( (ypid+1)*mody < ny_ ) endy = (ypid+1)*mody;
      else endy = ny_;
      startz = zpid*modz;
      if( (zpid+1)*modz < nz_ ) endz = (zpid+1)*modz;
      else endz = nz_;
      
      int NumMyElements = (endx-startx)*(endy-starty)*(endz-startz);
      int * MyGlobalElements = new int[NumMyElements];
      int count = 0;
      
      for( int i=startx ; i<endx ; ++i ) {
	for( int j=starty ; j<endy ; ++j ) {
	  for( int k=startz ; k<endz ; ++k ) {
	    MyGlobalElements[count++] = i+j*nx_+k*(nx_*ny_);
	  }
	}
      }
      
      map_ = new Epetra_Map (NumGlobalElements_,NumMyElements,MyGlobalElements,0,*comm_);

      // I delete this guy so that this case is not different from the
      // others, and I don't have to clean up this mess while
      // destroying the object.
      
      delete [] MyGlobalElements;
      
    } else if( MapType_ == "interlaced" ) {
      
      // this is the first funky map. Nodes are assigned so that
      // node 0 is given to proc 0, node 1 to proc 1, and
      // node i to proc i%NumProcs. Probably not the best, but it
      // results in decompositions with lots of boundary nodes.
      
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
      delete [] MyGlobalElements;
      
    } else if( MapType_ == "random" ) {
      
      // this is even funkier. Random decomposition of nodes into procs.
      // It should result in a more ordered decomposition than "interlaced"
      // This is the idea: I create the map on proc 0, then I broadcast
      // it to all procs. This is not very efficient, but saves some
      // MPI calls.
      
      int * part = new int[NumGlobalElements_];
      
      if( comm_->MyPID() == 0 ) {
	Epetra_Util Util;
	
	for( int i=0 ; i<NumGlobalElements_ ; ++i ) {
	  unsigned int r = Util.RandomInt();	
	  part[i] = r%(comm_->NumProc());
	}
      }
      
      comm_->Broadcast(part,NumGlobalElements_,0);
      
      // count the elements assigned to this proc
      int NumMyElements = 0;
      for( int i=0 ; i<NumGlobalElements_ ; ++i ) {
	if( part[i] == comm_->MyPID() ) NumMyElements++;
      }
      
      // get the loc2global list
      int * MyGlobalElements = new int[NumMyElements];
      int count = 0;
      for( int i=0 ; i<NumGlobalElements_ ; ++i ) {
	if( part[i] == comm_->MyPID() ) MyGlobalElements[count++] = i;
      }
      
      map_ = new Epetra_Map (NumGlobalElements_,NumMyElements,MyGlobalElements,
			     0,*comm_);
      
      delete [] MyGlobalElements;
      delete [] part;
      
    } else {
      
      cerr << ErrorMsg << "MapType has an incorrect value (" << MapType_ << ")\n";
      exit( EXIT_FAILURE );
    }
  }

  // local number of rows
  NumMyElements_ = map_->NumMyElements();
  // get update list
  MyGlobalElements_ = map_->MyGlobalElements( );

  if( verbose_ == true ) {
    cout << OutputMsg << "Time to create Map: "
	 << Time.ElapsedTime() << " (s)\n";
  }

  return;
}
  
// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrix(void)
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating Matrix...\n";
  }
  
  // HB matrices are different, as their dimension has to be read before.
  // Here the idea is to read the matrix on proc 0, then build the
  // map, then redistribute it linearly

  if( name_ == "hb" || name_ == "matrix_market" ||
      name_ == "triples_sym" || name_ == "triples_nonsym" ) {

    Epetra_Time Time(*comm_);
    ReadMatrix();
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

    else if( name_ == "laplace_1d_n" ) CreateMatrixLaplace1dNeumann();
      
    else if( name_ == "laplace_2d" ) CreateMatrixLaplace2d();

    else if( name_ == "laplace_2d_n" ) CreateMatrixLaplace2dNeumann();

    else if( name_ == "laplace_2d_9pt" ) CreateMatrixLaplace2d_9pt();

    else if( name_ == "recirc_2d" ) CreateMatrixRecirc2d();

    else if( name_ == "recirc_2d_divfree" ) CreateMatrixRecirc2dDivFree();

    else if( name_ == "uni_flow_2d" ) CreateMatrixUniFlow2d();
      
    else if( name_ == "laplace_3d" ) CreateMatrixLaplace3d();
      
    else if( name_ == "cross_stencil_2d" ) CreateMatrixCrossStencil2d();
      
    else if( name_ == "cross_stencil_3d" ) CreateMatrixCrossStencil3d();

    else if( name_ == "lehmer" ) CreateMatrixLehmer();

    else if( name_ == "minij" ) CreateMatrixMinij();

    else if( name_ == "ris" ) CreateMatrixRis();

    else if( name_ == "hilbert" ) CreateMatrixHilbert();

    else if( name_ == "jordblock" ) CreateMatrixJordblock();

    else if( name_ == "cauchy" ) CreateMatrixCauchy();

    else if( name_ == "fiedler" ) CreateMatrixFiedler();

    else if( name_ == "hanowa" ) CreateMatrixHanowa();

    else if( name_ == "kms" ) CreateMatrixKMS();

    else if( name_ == "parter" ) CreateMatrixParter();

    else if( name_ == "pei" ) CreateMatrixPei();

    else if( name_ == "ones" ) CreateMatrixOnes();

    else if( name_ == "vander" ) CreateMatrixVander();
    
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

  return;    
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateExactSolution(void)
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

    } else if( ExactSolutionType_ == "quad_x" ) {

      // always suppose to have Dirichlet boundary
      // conditions, and those points have already been eliminated
      // from the matrix
      double hx = 1.0/(NumGlobalElements_+1);
      for( int i=0 ; i<NumMyElements_ ; i++ ) {
	double x = (MyGlobalElements_[i]+1)*hx;
	(*ExactSolution_)[i] = x*(1.-x);
      }
      
    } else if( ExactSolutionType_ == "quad_xy" ) {

      SetupCartesianGrid2D();
      
      double hx = 1.0/(nx_+1);
      double hy = 1.0/(ny_+1);
      
      for( int i=0 ; i<NumMyElements_ ; ++i ) {
	int ix, iy;
	ix = (MyGlobalElements_[i])%nx_;
	iy = (MyGlobalElements_[i] - ix)/nx_;
	double x = hx*(ix+1);
	double y = hy*(iy+1);
	double u;
	ExactSolQuadXY(x,y,u);

	(*ExactSolution_)[i] = u;
	
      }
      
      
    } else {
      if( verbose_ ) {
	cerr << ErrorMsg << "exact solution type is not correct : "
	     << ExactSolutionType_ << endl;
	cerr << ErrorMsg << "It should be:\n"
	     << ErrorMsg << "<random> / <constant> / <quad_x> / <quad_xy>" << endl;
      }
      exit( EXIT_FAILURE );
    }
  }

  return;
  
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateStartingSolution(void)
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
      StartingSolution_->PutScalar(0.0);
    } else {
      cerr << ErrorMsg << "starting solution type is not correct : "
	   << StartingSolutionType_ << endl;
      exit( EXIT_FAILURE );
    }
  }

  return;
  
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateRHS(void)
{

  if( map_ == NULL ) CreateMap();
  if( matrix_ == NULL ) CreateMatrix();
  if( ExactSolution_ == NULL )  CreateExactSolution();

  if( rhs_ != NULL ) delete rhs_;
  
  Epetra_Time Time(*comm_);

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating RHS `" << RhsType_ << "' ...\n";
  }
  
  rhs_ = new Epetra_Vector(*map_);

  if( RhsType_ == "from_exact_solution" ) {

    matrix_->Multiply(false,*ExactSolution_,*rhs_);

  } else if( RhsType_ == "exact_rhs_uni_flow_2d" ) {

    // need to set a_ and b_ too 
    if( conv_ == UNDEF ) conv_ = 1;
    if( diff_ == UNDEF ) diff_ = 1e-5;
    if( alpha_ == UNDEF ) alpha_ = 1e-5;

    SetupCartesianGrid2D();
    
    double hx = 1.0/(nx_+1);
    double hy = 1.0/(ny_+1);
    
    for( int i=0 ; i<NumMyElements_ ; ++i ) {
      int ix, iy;
      ix = (MyGlobalElements_[i])%nx_;
      iy = (MyGlobalElements_[i] - ix)/nx_;
      double x = hx*(ix+1);
      double y = hy*(iy+1);
      double u, ux, uy, uxx, uyy;
      ExactSolQuadXY(x,y,u,ux,uy,uxx,uyy);
      
      (*rhs_)[i] = -diff_*( uxx + uyy )  // -b_ \nabla u
	+ conv_*cos(alpha_)*ux               // ux
	+ conv_*sin(alpha_)*uy;              // uy
      
    }
    
  } else if( RhsType_ == "exact_rhs_recirc_2d" ) {

    // need to set a_ and b_ too 
    if( conv_ == UNDEF ) conv_ = 1;
    if( diff_ == UNDEF ) diff_ = 1e-5;

    SetupCartesianGrid2D();
    
    double hx = 1.0/(nx_+1);
    double hy = 1.0/(ny_+1);
    
    for( int i=0 ; i<NumMyElements_ ; ++i ) {
      int ix, iy;
      ix = (MyGlobalElements_[i])%nx_;
      iy = (MyGlobalElements_[i] - ix)/nx_;
      double x = hx*(ix+1);
      double y = hy*(iy+1);
      double u, ux, uy, uxx, uyy;
      ExactSolQuadXY(x,y,u,ux,uy,uxx,uyy);      
      
      (*rhs_)[i] =  -diff_*( uxx + uyy )        // -b_ \nabla u
	+ conv_*4*x*(x-1.)*(1.-2*y)*ux          // ux
	- conv_*4*y*(y-1.)*(1.-2*x)*uy;         // uy
      
    }

  } else if( RhsType_ == "exact_rhs_laplace_2d" ) {

    SetupCartesianGrid2D();
    
    double hx = 1.0/(nx_+1);
    double hy = 1.0/(ny_+1);
    
    for( int i=0 ; i<NumMyElements_ ; ++i ) {
      int ix, iy;
      ix = (MyGlobalElements_[i])%nx_;
      iy = (MyGlobalElements_[i] - ix)/nx_;
      double x = hx*(ix+1);
      double y = hy*(iy+1);
      double u, ux, uy, uxx, uyy;
      ExactSolQuadXY(x,y,u,ux,uy,uxx,uyy);      
      
      (*rhs_)[i] = uxx+uyy;
      
    }
    
  } else {

    cerr << ErrorMsg << "RHS type not correct ("
	 << RhsType_ << ")" << endl;
    exit( EXIT_FAILURE );
  }
  
  if( verbose_ == true ) {
    cout << OutputMsg << "Time to create RHS (matvec): "
         << Time.ElapsedTime() << " (s)\n";
  }
  
  return;
  
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateEye(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `eye'...\n";
  }
  
  a_ = 1.0;
  CreateMatrixDiag();
  return;
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrixDiag(void) 
{

  // default value if not otherwise specified
  if( a_ == UNDEF ) a_ = 1;
  
  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `diag'...\n";
    cout << OutputMsg << "Diagonal element = " << a_ << endl;
  }

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,1);
  double Value;
    
  for( int i=0 ; i<NumMyElements_; ++i ) {

    int Indices = MyGlobalElements_[i];
    Value = a_;
      
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], 1, &Value, &Indices)==0);
      
  }
    
  assert(matrix_->FillComplete()==0);

  return;
    
}
  
// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrixTriDiag(void) 
{

  // default value if not otherwise specified
  if( a_ == UNDEF ) a_ = 2;
  if( b_ == UNDEF ) b_ = 1;
  if( c_ == UNDEF ) c_ = 1;

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `tridiag'...\n";
    cout << OutputMsg << "Row is [" << b_ << ", " << a_ << ", " << c_ << "]\n";
  }
  
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,3);

  double * Values = new double[2];
  int * Indices = new int[2];
  int NumEntries;

  for( int i=0 ; i<NumMyElements_; ++i ) {
    if (MyGlobalElements_[i]==0) {
      // off-diagonal for first row
      Indices[0] = 1;
      NumEntries = 1;
      Values[0] = c_;
    } else if (MyGlobalElements_[i] == NumGlobalElements_-1) {
      // off-diagonal for last row
      Indices[0] = NumGlobalElements_-2;
      NumEntries = 1;
      Values[0] = b_;
    } else {
      // off-diagonal for internal row
      Indices[0] = MyGlobalElements_[i]-1;
      Values[1] = b_;
      Indices[1] = MyGlobalElements_[i]+1;
      Values[0] = c_;
      NumEntries = 2;
    }
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], NumEntries, Values, Indices)==0);
    // Put in the diagonal entry
    Values[0] = a_;
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], 1, Values, MyGlobalElements_+i)==0);
  }
  
  // Finish up, trasforming the matrix entries into local numbering,
  // to optimize data transfert during matrix-vector products
  assert(matrix_->FillComplete()==0);

  delete [] Values;
  delete [] Indices;

  return;
    
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrixLaplace1d(void)
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `laplace_1d'...\n";
  }

  a_ = 2.0;
  b_ = -1.0;
  c_ = -1.0;

  CreateMatrixTriDiag();
  
  return;
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrixLaplace1dNeumann(void)
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `laplace_1d_n'...\n";
  }

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,3);

  double *Values = new double[2];
  int *Indices = new int[2];
  int NumEntries;

  for( int i=0 ; i<NumMyElements_; ++i ) {
    if (MyGlobalElements_[i]==0) {
      Indices[0] = 1;
      NumEntries = 1;
      Values[0] = -1.0;
    } else if (MyGlobalElements_[i] == NumGlobalElements_-1) {
      Indices[0] = NumGlobalElements_-2;
      NumEntries = 1;
      Values[0] = -1.0;
    } else {
      Indices[0] = MyGlobalElements_[i]-1;
      Values[1] = -1.0;
      Indices[1] = MyGlobalElements_[i]+1;
      Values[0] = -1.0;
      NumEntries = 2;
    }
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], NumEntries, Values, Indices)==0);
    // Put in the diagonal entry
    if (MyGlobalElements_[i]==0 || (MyGlobalElements_[i] == NumGlobalElements_-1) )
      Values[0] = 1.0;
    else
      Values[0] = 2.0;
    
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], 1, Values, MyGlobalElements_+i)==0);
  }
  
  // Finish up, trasforming the matrix entries into local numbering,
  // to optimize data transfert during matrix-vector products
  assert(matrix_->FillComplete()==0);

  delete [] Values;
  delete [] Indices;

  return;
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrixCrossStencil2d(void)
{

  // default values if not otherwise specified
  
  if( a_ == UNDEF ) a_ = 4;
  if( b_ == UNDEF ) b_ = 1;
  if( c_ == UNDEF ) c_ = 1;
  if( d_ == UNDEF ) d_ = 1;
  if( e_ == UNDEF ) e_ = 1;

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `cross_stencil_2d'...\n";
    cout << OutputMsg << "with values: a=" << a_ << ", b=" << b_ << ", c=" << c_ 
	 << ", d=" << d_ << ", e=" << e_ << endl;
  }

  
  SetupCartesianGrid2D();
     
  int left, right, lower, upper;
    
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,5);
    
  // Add  rows one-at-a-time
    
  double Values[4], diag;
  int Indices[4];

  //    e
  //  b a c
  //    d
  for( int i=0 ; i<NumMyElements_; ++i ) {
    int NumEntries=0;
    GetNeighboursCartesian2d(  MyGlobalElements_[i], nx_, ny_, 
			       left, right, lower, upper);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      Values[NumEntries] = b_;
      ++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      Values[NumEntries] = c_;
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      Values[NumEntries] = d_;
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      Values[NumEntries] = e_;
      ++NumEntries;
    }
    // put the off-diagonal entries
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], NumEntries, 
				       Values, Indices)==0);
    // Put in the diagonal entry
    diag = a_;
	
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], 1, 
				       &diag, MyGlobalElements_+i)==0);
  }
  matrix_->FillComplete();

  return;
      
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrixCrossStencil2dVector(void)
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `cross_stencil_2d'...\n";
  }

  SetupCartesianGrid2D();
    
  int left, right, lower, upper;
    
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,5);
    
  // Add  rows one-at-a-time
    
  double Values[4], diag;
  int Indices[4];

  //    e
  //  b a c
  //    d
  for( int i=0 ; i<NumMyElements_; ++i ) {

    int NumEntries=0;
    GetNeighboursCartesian2d(  MyGlobalElements_[i], nx_, ny_, 
			       left, right, lower, upper);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      Values[NumEntries] = (*VectorB_)[i];
      ++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      Values[NumEntries] = (*VectorC_)[i];
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      Values[NumEntries] = (*VectorD_)[i];
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      Values[NumEntries] = (*VectorE_)[i];
      ++NumEntries;
    }
    // put the off-diagonal entries
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], NumEntries, 
				       Values, Indices)==0);
    // Put in the diagonal entry
    diag = (*VectorA_)[i];
	
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], 1, 
				       &diag, MyGlobalElements_+i)==0);
  }
  matrix_->FillComplete();

  return;
      
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrixLaplace2dNeumann(void)
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `laplace_2d_n'...\n";
  }

  SetupCartesianGrid2D();  

  int left, right, lower, upper;
    
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,5);
    
  // Add  rows one-at-a-time
    
  double Values[4], diag;
  int Indices[4];
  
  //    e
  //  b a c
  //    d
  for( int i=0 ; i<NumMyElements_; ++i ) {

    bool isBorder = false;

    int NumEntries=0;
    GetNeighboursCartesian2d(  MyGlobalElements_[i], nx_, ny_, 
			       left, right, lower, upper);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    } else isBorder = true;
      
    if( right != -1 ) {
      Indices[NumEntries] = right;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    } else isBorder = true;
    
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    } else isBorder = true;
    
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    } else isBorder = true;
    
    // put the off-diagonal entries
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], NumEntries, 
				       Values, Indices)==0);
    // Put in the diagonal entry
    if( isBorder ) diag = 2;
    else           diag = 4;
	
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], 1, 
				       &diag, MyGlobalElements_+i)==0);
  }
  matrix_->FillComplete();

  return;
      
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrixLaplace2d_9pt(void)
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `laplace_2d_9pt'...\n";
  }

  SetupCartesianGrid2D();
    
  int left, right, lower, upper;
    
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,9);
    
  // Add  rows one-at-a-time
    
  double Values[8], diag;
  for( int i=0 ; i<8 ; ++i ) Values[i] = -1.0;
  int Indices[8];

  diag = 8.0;
  
  //  z3  e  z4
  //   b  a  c
  //  z1  d  z2
  for( int i=0 ; i<NumMyElements_; ++i ) {
    int NumEntries=0;
    GetNeighboursCartesian2d(  MyGlobalElements_[i], nx_, ny_, 
			       left, right, lower, upper);
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
    if( left != -1 && lower != -1 ) {
      Indices[NumEntries] = lower-1;
      ++NumEntries;
    }
    if( right != -1 && lower != -1 ) {
      Indices[NumEntries] = lower+1;
      ++NumEntries;
    }
    if( left != -1 && upper != -1 ) {
      Indices[NumEntries] = upper-1;
      ++NumEntries;
    }
    if( right != -1 && upper != -1 ) {
      Indices[NumEntries] = upper+1;
      ++NumEntries;
    }
    
    // put the off-diagonal entries
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], NumEntries, 
				       Values, Indices)==0);
    // Put in the diagonal entry
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], 1, 
				       &diag, MyGlobalElements_+i)==0);
  }
  matrix_->FillComplete();

  return;
      
}
// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrixLaplace2d(void)
{

  SetupCartesianGrid2D();

  double hx = 1.0/(nx_+1);
  double hy = 1.0/(ny_+1);

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `laplace_2d'...\n";
    if( Scaling ) {
      cout << OutputMsg << "nx = " << nx_ << ", ny = " << ny_ << endl;
      cout << OutputMsg << "hx = " << hx << ", hy = " << hy << endl;
    }
  }

  if( Scaling ) {
    a_ = 2.0/(hx*hx) + 2.0/(hy*hy);
    b_ = -1.0/(hx*hx);
    c_ = -1.0/(hx*hx);
    d_ = -1.0/(hy*hy);
    e_ = -1.0/(hy*hy);
  } else {
    a_ = 4;
    b_ = -1;
    c_ = -1;
    d_ = -1;
    e_ = -1;
  }
  
  CreateMatrixCrossStencil2d();

  return;
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrixRecirc2d(void)
{

  // default values if not specified otherwise

  if( conv_ == UNDEF ) conv_ = 1;
  if( diff_ == UNDEF ) diff_ = 1e-5;

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `recirc_2d'...\n";
    cout << OutputMsg << "with convection = " << conv_ << " and diffusion = " << diff_ << endl;
  }
  
  SetupCartesianGrid2D();
  
  if( VectorA_ ) delete VectorA_;
  if( VectorB_ ) delete VectorB_;
  if( VectorC_ ) delete VectorC_;
  if( VectorD_ ) delete VectorD_;
  if( VectorE_ ) delete VectorE_;
  
  if( VectorA_ == NULL )  VectorA_ = new Epetra_Vector(*map_);
  if( VectorB_ == NULL )  VectorB_ = new Epetra_Vector(*map_);
  if( VectorC_ == NULL )  VectorC_ = new Epetra_Vector(*map_);
  if( VectorD_ == NULL )  VectorD_ = new Epetra_Vector(*map_);
  if( VectorE_ == NULL )  VectorE_ = new Epetra_Vector(*map_);

  assert( VectorA_ != NULL ) ;
  assert( VectorB_ != NULL ) ;
  assert( VectorC_ != NULL ) ;
  assert( VectorD_ != NULL ) ;
  assert( VectorE_ != NULL ) ;
  
  VectorA_->PutScalar(0.0);
  VectorB_->PutScalar(0.0);
  VectorC_->PutScalar(0.0);
  VectorD_->PutScalar(0.0);
  VectorE_->PutScalar(0.0);
  
  double hx = 1.0/(nx_+1);
  double hy = 1.0/(ny_+1);

  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int ix, iy;
    ix = (MyGlobalElements_[i])%nx_;
    iy = (MyGlobalElements_[i] - ix)/nx_;
    double x = hx*(ix+1);
    double y = hy*(iy+1);
    double ConvX = conv_*4*x*(x-1.)*(1.-2*y)/hx;
    double ConvY = -conv_*4*y*(y-1.)*(1.-2*x)/hy;

    // convection part
    
    if( ConvX<0 ) {
      (*VectorC_)[i] += ConvX;
      (*VectorA_)[i] -= ConvX;
    } else {
      (*VectorB_)[i] -= ConvX;
      (*VectorA_)[i] += ConvX;
    }

    if( ConvY<0 ) {
      (*VectorE_)[i] += ConvY;
      (*VectorA_)[i] -= ConvY;
    } else {
      (*VectorD_)[i] -= ConvY;
      (*VectorA_)[i] += ConvY;
    }

    // add diffusion part
    (*VectorA_)[i] += diff_*2./(hx*hx) + diff_*2./(hy*hy);
    (*VectorB_)[i] -= diff_/(hx*hx);
    (*VectorC_)[i] -= diff_/(hx*hx);
    (*VectorD_)[i] -= diff_/(hy*hy);
    (*VectorE_)[i] -= diff_/(hy*hy);
      
    
  }

  CreateMatrixCrossStencil2dVector();

  return;
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrixRecirc2dDivFree(void)
{

  // default values if not specified otherwise

  if( conv_ == UNDEF ) conv_ = 1;
  if( diff_ == UNDEF ) diff_ = 1e-5;

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `recirc_2d_divfree'...\n";
    cout << OutputMsg << "with convection = " << conv_ << " and diffusion = " << diff_ << endl;    
  }


  SetupCartesianGrid2D();
  
  if( VectorA_ ) delete VectorA_;
  if( VectorB_ ) delete VectorB_;
  if( VectorC_ ) delete VectorC_;
  if( VectorD_ ) delete VectorD_;
  if( VectorE_ ) delete VectorE_;

  if( VectorA_ == NULL )  VectorA_ = new Epetra_Vector(*map_);
  if( VectorB_ == NULL )  VectorB_ = new Epetra_Vector(*map_);
  if( VectorC_ == NULL )  VectorC_ = new Epetra_Vector(*map_);
  if( VectorD_ == NULL )  VectorD_ = new Epetra_Vector(*map_);
  if( VectorE_ == NULL )  VectorE_ = new Epetra_Vector(*map_);

  VectorA_->PutScalar(0.0);
  VectorB_->PutScalar(0.0);
  VectorC_->PutScalar(0.0);
  VectorD_->PutScalar(0.0);
  VectorE_->PutScalar(0.0);
  
  double hx = 1.0/(nx_+1);
  double hy = 1.0/(ny_+1);

  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int ix, iy;
    ix = (MyGlobalElements_[i])%nx_;
    iy = (MyGlobalElements_[i] - ix)/nx_;
    double x = hx*(ix+1);
    double y = hy*(iy+1);
    double ConvX = conv_*2*y*(1.-x*x)/hx;
    double ConvY = -conv_*2*x*(1.-y*y)/hy;

    // convection part
    
    if( ConvX<0 ) {
      (*VectorC_)[i] += ConvX;
      (*VectorA_)[i] -= ConvX;
    } else {
      (*VectorB_)[i] -= ConvX;
      (*VectorA_)[i] += ConvX;
    }

    if( ConvY<0 ) {
      (*VectorE_)[i] += ConvY;
      (*VectorA_)[i] -= ConvY;
    } else {
      (*VectorD_)[i] -= ConvY;
      (*VectorA_)[i] += ConvY;
    }

    // add diffusion part
    (*VectorA_)[i] += diff_*2./(hx*hx) + diff_*2./(hy*hy);
    (*VectorB_)[i] -= diff_/(hx*hx);
    (*VectorC_)[i] -= diff_/(hx*hx);
    (*VectorD_)[i] -= diff_/(hy*hy);
    (*VectorE_)[i] -= diff_/(hy*hy);
      
    
  }

  CreateMatrixCrossStencil2d();

  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::CrsMatrixGallery::CreateMatrixUniFlow2d(void)
{

  // default values if not specified otherwise

  if( conv_ == UNDEF  ) conv_ = 1;
  if( diff_ == UNDEF  ) diff_ = 1e-5;
  if( alpha_ == UNDEF ) alpha_ = 0;

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `uni_flow_2d'...\n";
    cout << OutputMsg << "with convection = " << conv_ << ", diffusion = " << diff_ << endl;
    cout << OutputMsg << "and alpha = " << alpha_ << endl;

  }

  SetupCartesianGrid2D();
  
  if( VectorA_ ) delete VectorA_;
  if( VectorB_ ) delete VectorB_;
  if( VectorC_ ) delete VectorC_;
  if( VectorD_ ) delete VectorD_;
  if( VectorE_ ) delete VectorE_;
  
  if( VectorA_ == NULL )  VectorA_ = new Epetra_Vector(*map_);
  if( VectorB_ == NULL )  VectorB_ = new Epetra_Vector(*map_);
  if( VectorC_ == NULL )  VectorC_ = new Epetra_Vector(*map_);
  if( VectorD_ == NULL )  VectorD_ = new Epetra_Vector(*map_);
  if( VectorE_ == NULL )  VectorE_ = new Epetra_Vector(*map_);

  assert( VectorA_ != NULL ) ;
  assert( VectorB_ != NULL ) ;
  assert( VectorC_ != NULL ) ;
  assert( VectorD_ != NULL ) ;
  assert( VectorE_ != NULL ) ;

  VectorA_->PutScalar(0.0);
  VectorB_->PutScalar(0.0);
  VectorC_->PutScalar(0.0);
  VectorD_->PutScalar(0.0);
  VectorE_->PutScalar(0.0);
  
  double hx = 1.0/(nx_+1);
  double hy = 1.0/(ny_+1);

  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int ix, iy;
    ix = (MyGlobalElements_[i])%nx_;
    iy = (MyGlobalElements_[i] - ix)/nx_;

    double ConvX = conv_ * cos(alpha_) / hx;
    double ConvY = conv_ * sin(alpha_) / hy;

    // convection part
    
    if( ConvX<0 ) {
      (*VectorC_)[i] += ConvX;
      (*VectorA_)[i] -= ConvX;
    } else {
      (*VectorB_)[i] -= ConvX;
      (*VectorA_)[i] += ConvX;
    }

    if( ConvY<0 ) {
      (*VectorE_)[i] += ConvY;
      (*VectorA_)[i] -= ConvY;
    } else {
      (*VectorD_)[i] -= ConvY;
      (*VectorA_)[i] += ConvY;
    }

    // add diffusion part
    (*VectorA_)[i] += diff_*2./(hx*hx) + diff_*2./(hy*hy);
    (*VectorB_)[i] -= diff_/(hx*hx);
    (*VectorC_)[i] -= diff_/(hx*hx);
    (*VectorD_)[i] -= diff_/(hy*hy);
    (*VectorE_)[i] -= diff_/(hy*hy);
      
    
  }

  CreateMatrixCrossStencil2dVector();

  return;
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrixLaplace3d(void)
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

  return;
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrixCrossStencil3d(void)
{

  // default values if not otherwise specified

  if( a_ == UNDEF ) a_ = 7;
  if( b_ == UNDEF ) b_ = 1;
  if( c_ == UNDEF ) c_ = 1;
  if( d_ == UNDEF ) d_ = 1;
  if( e_ == UNDEF ) e_ = 1;
  if( f_ == UNDEF ) f_ = 1;
  if( g_ == UNDEF ) g_ = 1;

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `cross_stencil_3d'...\n";
    cout << OutputMsg << "with values: a=" << a_ << ", b=" << b_ << ", c=" << c_ << endl
	 << OutputMsg << "d=" << d_ << ", e=" << e_ << ", f=" << f_
	 << ", g=" << g_ << endl;
  }
  
  // problem size

  SetupCartesianGrid3D();
    
  int left, right, lower, upper, below, above;
    
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,7);
    
  // Add  rows one-at-a-time
    
  double Values[6], diag;
  int Indices[6];

  //    e 
  //  b a c
  //    d
  // + f below and g above
    
  for( int i=0 ; i<NumMyElements_; ++i ) {
    int NumEntries=0;
    GetNeighboursCartesian3d(MyGlobalElements_[i], nx_, ny_, nz_,
			     left, right, lower, upper, below, above);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      Values[NumEntries] = b_;
      ++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      Values[NumEntries] = c_;
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      Values[NumEntries] = d_;
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      Values[NumEntries] = e_;
      ++NumEntries;
    }
    if( below != -1 ) {
      Indices[NumEntries] = below;
      Values[NumEntries] = f_;
      ++NumEntries;
    }
    if( above != -1 ) {
      Indices[NumEntries] = above;
      Values[NumEntries] = g_;
      ++NumEntries;
    }
    // put the off-diagonal entries
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], NumEntries, 
				       Values, Indices)==0);
    // Put in the diagonal entry
    diag = a_;
	
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], 1, 
				       &diag, MyGlobalElements_+i)==0);
  }

  matrix_->FillComplete();
  return;
      
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrixCrossStencil3dVector(void)
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `cross_stencil_3d'...\n";
  }

  // problem size
  
  if( nx_ == -1 || ny_ == -1 || nz_ == -1 ) {
    nx_ = (int)pow(1.0*NumGlobalElements_,0.333334);
    ny_ = nx_;
    nz_ = nx_;
      
    if( nx_ * ny_ *nz_ != NumGlobalElements_ ) {
      cerr << ErrorMsg << "The number of global elements must be a perfect cube\n"
	   << ErrorMsg << "(now is " << NumGlobalElements_ << ")." << endl;
      exit( EXIT_FAILURE );
    }
  }
    
  int left, right, lower, upper, below, above;
    
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,7);
    
  // Add  rows one-at-a-time
    
  double Values[6], diag;
  int Indices[6];

  //    e 
  //  b a c
  //    d
  // + f below and g above
    
  for( int i=0 ; i<NumMyElements_; ++i ) {
    int NumEntries=0;
    GetNeighboursCartesian3d(MyGlobalElements_[i], nx_, ny_, nz_,
			     left, right, lower, upper, below, above);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      Values[NumEntries] = (*VectorB_)[i];
      ++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      Values[NumEntries] = (*VectorC_)[i];
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      Values[NumEntries] = (*VectorD_)[i];
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      Values[NumEntries] = (*VectorE_)[i];
      ++NumEntries;
    }
    if( below != -1 ) {
      Indices[NumEntries] = below;
      Values[NumEntries] = (*VectorF_)[i];
      ++NumEntries;
    }
    if( above != -1 ) {
      Indices[NumEntries] = above;
      Values[NumEntries] = (*VectorG_)[i];
      ++NumEntries;
    }
    // put the off-diagonal entries
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], NumEntries, 
				       Values, Indices)==0);
    // Put in the diagonal entry
    diag = (*VectorA_)[i];
	
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], 1, 
				       &diag, MyGlobalElements_+i)==0);
  }

  matrix_->FillComplete();
  return;
      
}
// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrixLehmer(void)
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
    int iGlobal = MyGlobalElements_[i];
    for( int jGlobal=0 ; jGlobal<NumGlobalElements_ ; ++jGlobal ) {
      if( iGlobal>=jGlobal) Values[jGlobal] = 1.0*(jGlobal+1)/(iGlobal+1);
      else                  Values[jGlobal] = 1.0*(iGlobal+1)/(jGlobal+1);
    }
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumGlobalElements_, Values, Indices)==0);
      
  }
    
  delete [] Indices;
  delete [] Values;

  assert(matrix_->FillComplete()==0);
  
  return;
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrixMinij(void)
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
    int iGlobal = MyGlobalElements_[i];
    for( int jGlobal=0 ; jGlobal<NumGlobalElements_ ; ++jGlobal ) {
      if( iGlobal>=jGlobal ) Values[jGlobal] = 1.0*(jGlobal+1);
      else                   Values[jGlobal] = 1.0*(iGlobal+1);
    }
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumGlobalElements_, Values, Indices)==0);
      
  }
    
  delete [] Indices;
  delete [] Values;

  assert(matrix_->FillComplete()==0);
  
  return;
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrixRis(void)
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
    int iGlobal = MyGlobalElements_[i];
    for( int jGlobal=0 ; jGlobal<NumGlobalElements_ ; ++jGlobal ) {
      Values[jGlobal] = 0.5/(NumGlobalElements_ -(iGlobal+1)-(jGlobal+1)+1.5);
      
    }
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumGlobalElements_, Values, Indices)==0);
      
  }
    
  delete [] Indices;
  delete [] Values;

  assert(matrix_->FillComplete()==0);
  
  return;
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateMatrixHilbert(void) 
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
    int iGlobal = MyGlobalElements_[i];
    for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
      Values[j] = 1.0/((iGlobal+1)+(j+1)-1);
    }
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumGlobalElements_, Values, Indices)==0);
      
  }
    
  delete [] Indices;
  delete [] Values;

  assert(matrix_->FillComplete()==0);
  
  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::CrsMatrixGallery::CreateMatrixJordblock(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `jordblock'...\n";
  }

  // default values is not specified
  if( a_ == UNDEF ) a_ = 0.1;
  
  // create matrix
  
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
    
  assert(matrix_->FillComplete()==0);
  
  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::CrsMatrixGallery::CreateMatrixCauchy(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `cauchy'...\n";
  }

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int * Indices = new int[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int NumEntries = NumGlobalElements_;
    int iGlobal = MyGlobalElements_[i];

    for( int jGlobal=0 ; jGlobal<NumGlobalElements_ ; ++jGlobal ) {
      Indices[jGlobal] = jGlobal;
      Values[jGlobal] = 1.0/(iGlobal+1+jGlobal+1);
    }

    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumEntries, Values, Indices)==0);
      
  }
    
  delete [] Indices;
  delete [] Values;
    
  assert(matrix_->FillComplete()==0);
  
  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::CrsMatrixGallery::CreateMatrixFiedler(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `fiedler'...\n";
  }

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int * Indices = new int[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int NumEntries = NumGlobalElements_;
    int iGlobal = MyGlobalElements_[i];
    for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
      Indices[j] = j;
      Values[j] = (double)abs(iGlobal-j);
    }

    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumEntries, Values, Indices)==0);
      
  }
    
  delete [] Indices;
  delete [] Values;
      
  assert(matrix_->FillComplete()==0);
  
  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::CrsMatrixGallery::CreateMatrixHanowa(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `hanowa'...\n";
  }

  // default values
  
  if( a_ == UNDEF ) a_ = -1;

  // problem size
  
  if( NumGlobalElements_ % 2 ) {
    cerr << ErrorMsg << "`hanowa' matrix requires a even number of points" << endl;
    exit( EXIT_FAILURE );
  }

  int half = NumGlobalElements_/2;
  
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,2);

  int Indices[2];
  double Values[2];
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int NumEntries = 2;
    int Global = MyGlobalElements_[i];
    Indices[0] = Global;
    if( Global < half ) Indices[1] = Global + half;
    else                Indices[1] = Global - half;
    Values[0] = a_;
    // convert from C style to FORTRAN style
    if( Global < half ) Values[1] = (double) -Global - 1;
    else                Values[1] = (double)  (Global - half) + 1;

    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumEntries, Values, Indices)==0);
    
  }
    
  assert(matrix_->FillComplete()==0);
  
  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::CrsMatrixGallery::CreateMatrixKMS(void) 
{
  
  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `kms'...\n";
  }

  // default values
  
  if( a_ == UNDEF ) a_ = 0.5;
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int * Indices = new int[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {

    int NumEntries = NumGlobalElements_;
    int iGlobal = MyGlobalElements_[i];
    for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
      Indices[j] = j;
      Values[j] = pow(a_,abs(iGlobal-j));
    }

    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumEntries, Values, Indices)==0);
    
  }

  delete [] Indices;
  delete [] Values;
      
  assert(matrix_->FillComplete()==0);
  
  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::CrsMatrixGallery::CreateMatrixParter(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `parter'...\n";
  }

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int * Indices = new int[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {

    int NumEntries = NumGlobalElements_;
    int iGlobal = MyGlobalElements_[i];
    for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
      Indices[j] = j;
      Values[j] = 1.0/(iGlobal-j+0.5);
    }

    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumEntries, Values, Indices)==0);
    
  }


  delete [] Indices;
  delete [] Values;
      
  assert(matrix_->FillComplete()==0);
  
  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::CrsMatrixGallery::CreateMatrixPei(void) 
{

  // default values if not specified otherwise
  a_ = 1.0;
  
  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `pei'...\n";
    cout << OutputMsg << "with value a=" << a_ << endl;
  }

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int * Indices = new int[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {

    int NumEntries = NumGlobalElements_;

    for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
      Indices[j] = j;
      Values[j] = 1.0;
      if( j == MyGlobalElements_[i] ) Values[j] += a_;
    }

    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumEntries, Values, Indices)==0);
    
  }

  delete [] Indices;
  delete [] Values;
      
  assert(matrix_->FillComplete()==0);
  
  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::CrsMatrixGallery::CreateMatrixOnes(void) 
{

  // default values if not specified otherwise
  if( a_ == UNDEF ) a_ = 1.0;
  
  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `ones'...\n";
    cout << OutputMsg << "with value a=" << a_ << endl;
  }

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int * Indices = new int[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {

    int NumEntries = NumGlobalElements_;

    for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
      Indices[j] = j;
      Values[j] = a_;
    }

    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumEntries, Values, Indices)==0);
    
  }

  delete [] Indices;
  delete [] Values;
      
  assert(matrix_->FillComplete()==0);
  
  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::CrsMatrixGallery::CreateMatrixVander(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `vander'...\n";
  }

  // create matrix
  
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int * Indices = new int[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {

    int NumEntries = NumGlobalElements_;

    for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
      Indices[j] = j;
      Values[j] = pow((*VectorA_)[i],NumGlobalElements_-j-1);
    }

    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumEntries, Values, Indices)==0);
    
  }

  delete [] Indices;
  delete [] Values;
      
  assert(matrix_->FillComplete()==0);
  
  return;
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::ReadMatrix(void)
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Reading " << name_ << "  matrix `"
	 << FileName_ << "'...\n";
  }

  Epetra_Map * readMap;
  Epetra_CrsMatrix * readA; 
  Epetra_Vector * readx; 
  Epetra_Vector * readb;
  Epetra_Vector * readxexact;
    
  // Call routine to read in problem from file

  if( name_ == "hb" ) 
    Trilinos_Util_ReadHb2Epetra((char*)FileName_.c_str(), *comm_, readMap,
				readA, readx, 
				readb, readxexact);
  else if( name_ == "matrix_market" )
    Trilinos_Util_ReadMatrixMarket2Epetra((char*)FileName_.c_str(), *comm_,
					  readMap, readA, readx, 
					  readb, readxexact );
  else if( name_ == "triples_sym" )
    Trilinos_Util_ReadTriples2Epetra((char*)FileName_.c_str(), false, *comm_,
				     readMap, readA, readx, 
				     readb, readxexact );
  else if( name_ == "triples_nonsym" )
    Trilinos_Util_ReadTriples2Epetra((char*)FileName_.c_str(), true, *comm_,
				     readMap, readA, readx, 
				     readb, readxexact );
  else {
    cerr << ErrorMsg << "problem type not correct (" << name_ << ")\n";
    exit( EXIT_FAILURE );
  }
    
  NumGlobalElements_ = readMap->NumGlobalElements();

  if( map_ != NULL ) delete map_;
  
  // create map for matrix. Use the normal function CreateMap
  // if the user has not specified "greedy" as map type.
  // In this latter case, form on proc 0 a map corresponding to
  // the greedy algorithm. This is some kind of graph decomposition
  // stuff, only cheaper (and that does not required any external
  // library)

  if( MapType_ == "greedy" ) {

    int * part = new int[NumGlobalElements_];

    if( comm_->MyPID() == 0 ) {

      int NumProcs = comm_->NumProc();
      int * ElementsPerDomain = new int[NumProcs];
      int * count = new int[NumProcs];

      // define how many nodes have to be put on each proc
      
      int div = NumGlobalElements_/NumProcs;
      int mod = NumGlobalElements_%NumProcs;

      for( int i=0 ; i<NumProcs ; ++i ) {
	count[i] = 0;
	ElementsPerDomain[i] = div;
	if( i<mod ) ElementsPerDomain[i]++;
      }
      
      for( int i=0 ; i<NumGlobalElements_ ; ++i ) {
	part[i] = -1;
      }
      
      int MaxNnzPerRow = readA->MaxNumEntries();
      if( MaxNnzPerRow == 0 ) {
	cerr << ErrorMsg << "something went wrong in `CreateMatrix'\n"
	     << ErrorMsg << "MaxNnzPerRow == 0 \n";
	exit( EXIT_FAILURE );
      }

      int CrsNumEntries;
      int * CrsIndices;
      double * CrsValues;

      // start from row 0, assigned to domain 0
      int RootNode = 0;
      part[0] = 0;      
      int CurrentDomain = 0;
      
      bool ok = true;
      
      while( ok == true ) {

	readA->ExtractMyRowView(RootNode,CrsNumEntries,
				CrsValues,CrsIndices);

	ok = false;
	
	for( int j=0 ; j<CrsNumEntries ; ++j ) {

	  if( count[CurrentDomain] == ElementsPerDomain[CurrentDomain] ) {
	    CurrentDomain++;
	  }
	  
	  if( part[CrsIndices[j]] == -1 ) {
	    part[CrsIndices[j]] = CurrentDomain;
	    if( ok == false ) {
	      ok = true;
	      RootNode = CrsIndices[j];
	    }
	    count[CurrentDomain]++;
	  }
	}

	// check if some -1 nodes are still available
	if( ok == false ) {
	  for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
	    if( part[j] == -1 ) {
	      RootNode = j;
	      ok = true;
	      break;
	    }
	  }
	}
	      
      }

      delete [] ElementsPerDomain;
      delete [] count;
	
    }

    // now broadcast on all procs. This might be pretty expensive...
    comm_->Broadcast(part,NumGlobalElements_,0);

    for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
      if( part[j] == -1 ) {
	cerr << ErrorMsg << "part[" << j << "] = -1 \n";
      }
    }    

    // count the elements assigned to this proc
    int NumMyElements = 0;
    for( int i=0 ; i<NumGlobalElements_ ; ++i ) {
      if( part[i] == comm_->MyPID() ) NumMyElements++;
    }

    // get the loc2global list
    int * MyGlobalElements = new int[NumMyElements];
    int count = 0;
    for( int i=0 ; i<NumGlobalElements_ ; ++i ) {
      if( part[i] == comm_->MyPID() ) MyGlobalElements[count++] = i;
    }

    map_ = new Epetra_Map (NumGlobalElements_,NumMyElements,MyGlobalElements,
			   0,*comm_);
    
    delete [] MyGlobalElements;
    delete [] part;

  } else {    
    CreateMap();
  }

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
  
  assert(matrix_->FillComplete()==0);    

  delete readA;
  delete readx;
  delete readb;
  delete readxexact;
  delete readMap;
  
  // local number of rows
  NumMyElements_ = map_->NumMyElements();
  // get update list
  MyGlobalElements_ = map_->MyGlobalElements( );
  
  return;

}

// ================================================ ====== ==== ==== == =  
void Trilinos_Util::CrsMatrixGallery::GetNeighboursCartesian2d( const int i, const int nx, const int ny,
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
void Trilinos_Util::CrsMatrixGallery::GetNeighboursCartesian3d( const int i, const int nx, const int ny, const int nz,
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

  GetNeighboursCartesian2d( ixy, nx, ny, left, right, lower, upper);
    
  if( left != -1 ) left += iz*(nx*ny);
  if( right != -1 ) right += iz*(nx*ny);
  if( lower != -1 ) lower += iz*(nx*ny);
  if( upper != -1 ) upper += iz*(nx*ny);

  return;

}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::ZeroOutData() 
{
  NumGlobalElements_ = -1;
  nx_ = -1;    ny_ = -1;     nz_ = -1;
  mx_ = -1;    mx_ = -1;     mz_ = -1;

  a_ = UNDEF, b_ = UNDEF, c_ = UNDEF, d_ = UNDEF, e_ = UNDEF, f_ = UNDEF, g_ = UNDEF;
  alpha_ = UNDEF;
  beta_  = UNDEF;
  gamma_ = UNDEF;
  delta_ = UNDEF;

  conv_ = UNDEF;
  diff_ = UNDEF;
  source_ = UNDEF;
  
  VectorA_ = NULL;
  VectorB_ = NULL;
  VectorC_ = NULL;
  VectorD_ = NULL;
  VectorE_ = NULL;
  VectorF_ = NULL;
  VectorG_ = NULL;
    
  map_ = NULL;
  matrix_ = NULL;
  ExactSolution_ = NULL;
  StartingSolution_ = NULL;
  rhs_ = NULL;

  MapType_ = "linear";
  ExactSolutionType_ = "constant";
  StartingSolutionType_ = "zero";
  ExpandType_ = "zero_off_diagonal";
  RhsType_ = "from_exact_solution";
  
  NumPDEEqns_= 1;

  LinearProblem_ = NULL;
    
}

void Trilinos_Util::CrsMatrixGallery::PrintMatrixAndVectors() 
{
  PrintMatrixAndVectors(cout);
}

void Trilinos_Util::CrsMatrixGallery::PrintMatrixAndVectors(ostream & os) 
{

  if( comm_->MyPID() == 0 ) {
    os << "*** MATRIX ***\n";
  }

  os << *matrix_;

  if( comm_->MyPID() == 0 ) {
    os << "*** RHS ***\n";
  }
  
  os << *rhs_;

  return;

}

namespace Trilinos_Util {

ostream & operator << (ostream& os,
		                      const Trilinos_Util::CrsMatrixGallery & G )
{

  bool verbose = (G.comm_->MyPID() == 0);

  if( verbose ) {
    
    os << " * Solving problem " << G.name_ << endl;
    os << " * Number of global elements : " << G.NumGlobalElements_ << endl;
    os << " * Type of Map : " << G.MapType_ << endl;
    os << " * Number of PDEs : " << G.NumPDEEqns_ << endl;

    // CRS stuff
    if( G.matrix_ != NULL ) {
      os << " * the matrix has been created " << endl;
      os << " * Matrix->OperatorDomainMap().NumGlobalElements() = "
	 << G.matrix_->OperatorDomainMap().NumGlobalElements() << endl;
    }
    if( G.ExactSolution_ != NULL ) 
      os << " * an exact solution (" << G.ExactSolutionType_
	 << ") has been created " << endl;
    if( G.rhs_ != NULL ) 
      os << " * the RHS has been created " << endl;
  }

  //  os << " * On proc " << G.comm_->MyPID() << " there are "
  //     << G.NumMyElements_ << " elements" << endl;

  return os;
  
}
} //namespace Trilinos_Util

void Trilinos_Util::CrsMatrixGallery::GetCartesianCoordinates(double * & x,
							 double * & y,
							 double * & z)
{

  if( map_ == NULL ) CreateMap();
  
  double length = 1.0;
  double delta_x, delta_y, delta_z;

  int ix, iy, iz;
  
  if( name_ == "diag" || name_ == "tridiag"  ||
      name_ == "laplace_1d" || name_ == "eye" ) {
    
    delta_x = length/(nx_-1);

    x = new double[NumMyElements_];
    assert( x != 0 );
    
    for( int i=0 ; i<NumMyElements_ ; ++i ) {

      ix = MyGlobalElements_[i];
      x[i] = delta_x * ix;

    }
    
  } else  if( name_ == "laplace_2d" || name_ == "cross_stencil_2d"
	      || name_ == "laplace_2d_9pt" || name_ == "recirc_2d"
	      || name_ == "laplace_2d_n" || name_ == "uni_flow_2d" ) {
  
    delta_x = length/(nx_-1);
    delta_y = length/(ny_-1);

    x =  new double[NumMyElements_];
    y =  new double[NumMyElements_];
    assert( x != 0 ); assert( y != 0 );
    
    for( int i=0 ; i<NumMyElements_ ; ++i ) {

      ix = MyGlobalElements_[i]%nx_;
      iy = (MyGlobalElements_[i] - ix)/ny_;

      x[i] = delta_x * ix;
      y[i] = delta_y * iy;

    }
    
  } else if( name_ == "laplace_3d" || name_ == "cross_stencil_3d" ) {
  
    delta_x = length/(nx_-1);
    delta_y = length/(ny_-1);
    delta_z = length/(nz_-1);

    x =  new double[NumMyElements_];
    y =  new double[NumMyElements_];
    z =  new double[NumMyElements_];
    assert( x != 0 ); assert( y != 0 ); assert( z != 0 );
    
    for( int i=0 ; i<NumMyElements_ ; i++ ) {

      int ixy = MyGlobalElements_[i]%(nx_*ny_);
      iz = (MyGlobalElements_[i] - ixy)/(nx_*ny_);

      ix = ixy%nx_;
      iy = (ixy - ix)/ny_;
      
      x[i] = delta_x * ix;
      y[i] = delta_y * iy;
      z[i] = delta_z * iz;

    }
    
  } else {

      cerr << ErrorMsg << "You can build Cartesian coordinates" << endl
	   << ErrorMsg << "only with one of the following problem_type:" << endl
	   << ErrorMsg << "<diag> / <tridiag> / <laplace_1d> / <eye>" << endl
	   << ErrorMsg << "<laplace_2d> / <cross_stencil_2d> / <laplace_2d_9pt> / <recirc_2d>" << endl
	   << ErrorMsg << "<laplace_2d_n> / <uni_flow_n>" << endl
	   << ErrorMsg << "<laplace_3d> / <cross_stencil_3d>" << endl;

      exit( EXIT_FAILURE );
  }

  return;
  
}

#include <iostream>
#include <fstream>

// ================================================ ====== ==== ==== == =
int Trilinos_Util::CrsMatrixGallery::WriteMatrix( const string & FileName, const bool UseSparse ) 

{

  // create matrix is not done yet
  
  if( matrix_ == NULL ) CreateMatrix();
  
  int NumMyRows = matrix_->NumMyRows(); // number of rows on this process
  int NumNzRow;   // number of nonzero elements for each row
  int NumEntries; // number of extracted elements for each row
  int NumGlobalRows; // global dimensio of the problem
  int GlobalRow;  // row in global ordering
  int NumGlobalNonzeros; // global number of nonzero elements

  NumGlobalRows = matrix_->NumGlobalRows();
  NumGlobalNonzeros = matrix_->NumGlobalNonzeros();

  // print out on cout if no filename is provided

  int IndexBase = matrix_->IndexBase(); // MATLAB start from 0
  if( IndexBase == 0 ) IndexBase = 1; 

  // write on file the dimension of the matrix

  if( comm_->MyPID() == 0 ) {
    
    ofstream fout(FileName.c_str());
      
    if( UseSparse ) {
      fout << "A = spalloc(";
      fout << NumGlobalRows << ',' << NumGlobalRows;
      fout << ',' << NumGlobalNonzeros << ");\n";
    } else {
      fout << "A = zeros(";
      fout << NumGlobalRows << ',' << NumGlobalRows << ");\n";
    }

    fout.close();
    
  }

  for( int Proc=0 ; Proc<comm_->NumProc() ; ++Proc ) {

    if( comm_->MyPID() == Proc ) {

      ofstream fout(FileName.c_str(),std::ios::app);
      
      fout << "% On proc " << Proc << ": ";
      fout << NumMyRows << " rows and ";
      fout << matrix_->NumMyNonzeros() << " nonzeros\n";

      // cycle over all local rows to find out nonzero elements
      for( int MyRow=0 ; MyRow<NumMyRows ; ++MyRow ) {
	
	GlobalRow = matrix_->GRID(MyRow);
	
	NumNzRow = matrix_->NumMyEntries(MyRow);
	double *Values = new double[NumNzRow];
	int *Indices = new int[NumNzRow];
	
	matrix_->ExtractMyRowCopy(MyRow, NumNzRow, 
			   NumEntries, Values, Indices);
	// print out the elements with MATLAB syntax
	for( int j=0 ; j<NumEntries ; ++j ) {
	  fout << "A(" << GlobalRow  + IndexBase 
	       << "," << matrix_->GCID(Indices[j]) + IndexBase
	       << ") = " << Values[j] << ";\n";
	}
	
	delete Values;
	delete Indices;

      }

      fout.close();

    }
    comm_->Barrier();
  }

  if( comm_->MyPID() == 0 ) {
    ofstream fout(FileName.c_str(),std::ios::app);
    fout << "%End of Matrix Output\n";
    fout.close();
  }
  
  return true;

}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::SetupCartesianGrid2D() 
{
  // needs a square number of nodes or
  // nx and ny set
  if( nx_ == -1 || ny_ == -1 ) {
    nx_ = (int)sqrt((double)NumGlobalElements_);
    ny_ = nx_;
    if( nx_ * ny_ != NumGlobalElements_ ) {
      cerr << ErrorMsg << "The number of global elements must be a perfect square\n"
	   << ErrorMsg << "otherwise set nx and ny. " << endl
	   << ErrorMsg << "(now NumGlobalElements = " << NumGlobalElements_ << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::SetupCartesianGrid3D() 
{
  // needs a cube number of nodes or
  // nx, ny and nz set
  if( nx_ == -1 || ny_ == -1 || nz_ == -1 ) {
    nx_ = (int)pow((double)NumGlobalElements_,0.333334);
    ny_ = nx_;
    nz_ = nx_;
    if( nx_ * ny_ * nz_ != NumGlobalElements_ ) {
      cerr << ErrorMsg << "The number of global elements must be a perfect cube\n"
	   << ErrorMsg << "otherwise set nx, ny, and nz. " << endl
	   << ErrorMsg << "(now NumGlobalElements = " << NumGlobalElements_ << ")" << endl;
      exit( EXIT_FAILURE );
    }
  }
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::ExactSolQuadXY(double x, double y,
						    double & u)
{

  u = x*(1.-x)*y*(1.-y);

  return;
  
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::ExactSolQuadXY(double x, double y,
						    double & u,
						    double & ux, double & uy,
						    double & uxx, double & uyy)
{

  u = x*(1.-x)*y*(1.-y);
  ux = (1-2*x)*y*(1.-y);
  uy = x*(1.-x)*(1.-2*y);
  uxx = -2*(x-x*x);
  uyy = -2*(y-y*y);

  return;

}

// ================================================ ====== ==== ==== == =
Epetra_Vector * Trilinos_Util::VbrMatrixGallery::GetVbrRHS(void)
{
  if( VbrRhs_ == NULL ) CreateVbrRHS();
  return VbrRhs_;
}

// ================================================ ====== ==== ==== == =
Epetra_Vector * Trilinos_Util::VbrMatrixGallery::GetVbrExactSolution(void)
{
  if( VbrExactSolution_ == NULL ) CreateVbrExactSolution();
  return VbrExactSolution_;
}

// ================================================ ====== ==== ==== == =
Epetra_Vector * Trilinos_Util::VbrMatrixGallery::GetVbrStartingSolution(void)
{
  if( VbrStartingSolution_ == NULL ) CreateVbrStartingSolution();
  return VbrStartingSolution_;
}

// ================================================ ====== ==== ==== == =  
Epetra_VbrMatrix * Trilinos_Util::VbrMatrixGallery::GetVbrMatrix(const int NumPDEEqns) 
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
Epetra_VbrMatrix * Trilinos_Util::VbrMatrixGallery::GetVbrMatrix(void)
{
    
  if( VbrMatrix_ == NULL ) CreateVbrMatrix();

  return VbrMatrix_;
    
}

// ================================================ ====== ==== ==== == =
Epetra_VbrMatrix & Trilinos_Util::VbrMatrixGallery::GetVbrMatrixRef(void)
{
    
  if( VbrMatrix_ == NULL ) CreateVbrMatrix();

  return *VbrMatrix_;
    
}

// ================================================ ====== ==== ==== == =

Epetra_LinearProblem * Trilinos_Util::VbrMatrixGallery::GetVbrLinearProblem(void) 
{
  // pointers, not really needed
  Epetra_VbrMatrix * A;
  Epetra_Vector * RHS;
  Epetra_Vector * StartingSolution;

  A = GetVbrMatrix();
  RHS = GetVbrRHS();
  StartingSolution = GetVbrStartingSolution();

  // create linear problem
  if( VbrLinearProblem_ != NULL ) delete VbrLinearProblem_;
  VbrLinearProblem_ = new Epetra_LinearProblem(A,StartingSolution,RHS);

  return VbrLinearProblem_;

}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::VbrMatrixGallery::CreateVbrExactSolution(void) 
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

  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::VbrMatrixGallery::CreateVbrStartingSolution(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating Starting Solution (VBR)...\n";
  }

  if( VbrStartingSolution_ != NULL ) {
    delete VbrStartingSolution_;
    VbrStartingSolution_ = NULL;
  }
    
  // need a rhs for crs
  if( StartingSolution_ == NULL ) CreateStartingSolution();
  // need a block map based on map_
  if( BlockMap_ == NULL ) CreateBlockMap();
  // now we can expand to the Vbr format
  VbrStartingSolution_ = new Epetra_Vector(*BlockMap_);
  for( int j=0 ; j<NumMyElements_ ; j++ ) 
    for( int i=0 ; i<NumPDEEqns_ ; ++i ) {
      (*VbrStartingSolution_)[j*NumPDEEqns_+i] = (*StartingSolution_)[j];
    }

  return;
  
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::VbrMatrixGallery::CreateVbrRHS(void) 
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

  return;
  
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::VbrMatrixGallery::CreateVbrMatrix(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating VBR matrix...\n";
  }

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
    
  int * VbrIndices = new int[MaxNnzPerRow];
  double * VbrValues = new double[MaxBlockSize];
  int BlockRows = NumPDEEqns_;
  int ierr;
    
  // cycle over all the local rows. 
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    
    // get GID of local row
    int GlobalNode = MyGlobalElements_[i];
    // extract Crs row

    ierr = matrix_->ExtractMyRowView(i,CrsNumEntries,
				     CrsValues,CrsIndices);

    // matrix_ is in local form. Need global indices
    for( int kk=0 ; kk<CrsNumEntries ; ++kk) 
      VbrIndices[kk] = matrix_->GCID(CrsIndices[kk]);
    
    // with VBR matrices, we have to insert one block at time.
    // This required two more instructions, one to start this
    // process (BeginInsertGlobalValues), and another one to
    // commit the end of submissions (EndSubmitEntries).
    
    VbrMatrix_->BeginInsertGlobalValues(GlobalNode, CrsNumEntries, VbrIndices);

    int ExpandTypeInt;
    
    if( ExpandType_ == "zero_off_diagonal" ) ExpandTypeInt=0;
    else if( ExpandType_ == "random_off_diagonal" ) ExpandTypeInt=1;
    else {
      cerr << ErrorMsg << "ExpandType not correct (" << ExpandType_ << "\n";
      exit( EXIT_FAILURE );
    }
    Epetra_Util Util;
    double r;
    
    for( int i=0 ; i<CrsNumEntries ; ++i ) {
	
      for( int k=0 ; k<BlockRows ; ++k ) { // rows
	for( int h=0 ; h<BlockRows ; ++h ) { // cols
	  if( k == h ) VbrValues[k+h*BlockRows] = CrsValues[i];
	  else {
	    switch( ExpandTypeInt ) {
	    case 0:
	      r = 0.0;
	      break;
	    case 1:
	      // get a double between -1 and 1
	      r = Util.RandomDouble();
	      // scale it so that the sum of the block off-diagonal
	      // is not greater than the block diangonal
	      r /= (1.5*CrsValues[i]*BlockRows);
	      break;
	    }
	    VbrValues[k+h*BlockRows] = r; 
	  }
	}
      }
	  
      VbrMatrix_->SubmitBlockEntry(VbrValues,BlockRows,BlockRows,BlockRows);
	
    }

    VbrMatrix_->EndSubmitEntries();
  }
    
  delete [] VbrIndices;
  delete [] VbrValues;

  VbrMatrix_->FillComplete();

  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::VbrMatrixGallery::ComputeResidualVbr(double & residual)
{

  // create solution and rhs if needed (matrix_ and ExactSolution are
  //  created by CreateRHS is necessary)
  if( VbrRhs_ == NULL ) CreateVbrRHS();

  Epetra_Vector Ax(*BlockMap_);
  assert(VbrMatrix_->Multiply(false, *VbrStartingSolution_, Ax)==0);

  assert(Ax.Update(1.0, *VbrRhs_, -1.0)==0);

  assert(Ax.Norm2(&residual)==0);

  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::VbrMatrixGallery::ComputeDiffBetweenStartingAndExactSolutionsVbr(double & residual)
{

  // create solution and rhs if needed (matrix_ and ExactSolution are
  //  created by CreateRHS is necessary)
  if( VbrRhs_ == NULL ) CreateVbrRHS();

  Epetra_Vector temp(*BlockMap_);

  assert(temp.Update(1.0, *VbrExactSolution_, -1.0, *VbrStartingSolution_, 0.0)==0);

  assert(temp.Norm2(&residual)==0);

  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::VbrMatrixGallery::PrintVbrMatrixAndVectors(ostream & os)
{

  if( comm_->MyPID() == 0 ) {
    os << "*** MATRIX (VBR) ***\n";
  }

  os << *VbrMatrix_;

  if( comm_->MyPID() == 0 ) {
    os << "*** RHS (VBR) ***\n";
  }

  os << *VbrRhs_;

  return;

}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::VbrMatrixGallery::PrintVbrMatrixAndVectors() 
{
  PrintVbrMatrixAndVectors(cout);
}

// ================================================ ====== ==== ==== == =
const Epetra_BlockMap * Trilinos_Util::VbrMatrixGallery::GetBlockMap(void)
{
  if( BlockMap_ == NULL ) CreateBlockMap();
    
  return BlockMap_;
}

// ================================================ ====== ==== ==== == =
const Epetra_BlockMap & Trilinos_Util::VbrMatrixGallery::GetBlockMapRef(void)
{
  if( BlockMap_ == NULL ) CreateBlockMap();
    
  return *BlockMap_;
}


// ================================================ ====== ==== ==== == =
void Trilinos_Util::VbrMatrixGallery::CreateBlockMap(void) 
{
        
  if( verbose_ == true ) {
    cout << OutputMsg << "Creating BlockMap...\n";
  }

  if( map_ == NULL ) CreateMap();

  Epetra_Time Time(*comm_);
    
  if( NumPDEEqns_ <= 0 ) {
    cerr << ErrorMsg << "NumPDEEqns not correct (" << NumPDEEqns_ << "(\n";
    cerr << ErrorMsg << "Set it to 1\n";
    NumPDEEqns_ = 1;
  }

  MaxBlkSize_ = NumPDEEqns_;
  
  BlockMap_ = new Epetra_BlockMap(NumGlobalElements_,NumMyElements_,
				  MyGlobalElements_, 
				  NumPDEEqns_,0,*comm_);

  if( verbose_ == true ) {
    cout << OutputMsg << "Time to create BlockMap: "
	 << Time.ElapsedTime() << " (s)\n";
  }

  return;
  
}


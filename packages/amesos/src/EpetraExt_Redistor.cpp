// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Import.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "Epetra_Time.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "EpetraExt_Redistor.h"


//=============================================================================

int EpetraExt_Redistor::SourceImport(Epetra_DistObject & Source, Epetra_DistObject & Target,
				     Epetra_CombineMode CombineMode) const
{
  return( Source.Import(Target,*SourceImporter_,CombineMode) );
}

//=============================================================================

int EpetraExt_Redistor::TargetImport(Epetra_DistObject & Source, Epetra_DistObject & Target,
				     Epetra_CombineMode CombineMode) const
{
  return( Target.Import(Source,*TargetImporter_,CombineMode) );
}

//=============================================================================

EpetraExt_Redistor::EpetraExt_Redistor(const Epetra_BlockMap & SourceMap,
				       const int NumTargetProcs )
  
{
  // creates new maps, put View mode ??
  SourceMap_ = new Epetra_BlockMap(SourceMap);
  assert( SourceMap_ != 0 );

  // handle map with constant element size only
  if( SourceMap.ConstantElementSize() == false ) {
    cerr << "Only maps with constant element size are supported" << endl;
    exit( EXIT_FAILURE );
  }
  
  int NumGlobalElements = SourceMap.NumGlobalElements();
  int NumMyElements = 0;
  int ElementsPerTargetProc = NumGlobalElements/NumTargetProcs;
  if( SourceMap.Comm().MyPID() == 0 ) ElementsPerTargetProc+=NumGlobalElements%NumTargetProcs;
  if( SourceMap.Comm().MyPID() < NumTargetProcs ) NumMyElements = ElementsPerTargetProc;

  int IndexBase = SourceMap.IndexBase();
  int ElementSize = SourceMap.ElementSize();
  
  TargetMap_ = new Epetra_BlockMap(NumGlobalElements,NumMyElements,ElementSize,IndexBase,SourceMap.Comm());
  assert( TargetMap_ != 0 );
  
  SourceImporter_ = new Epetra_Import(*SourceMap_,*TargetMap_);
  TargetImporter_ = new Epetra_Import(*TargetMap_,*SourceMap_);

  assert( SourceImporter_ != 0 );
  assert( TargetImporter_ != 0 );
}

//=============================================================================

EpetraExt_Redistor::EpetraExt_Redistor(const Epetra_BlockMap & SourceMap,
				       const Epetra_BlockMap & TargetMap )
{
  // creates new maps, put View mode ??
  SourceMap_ = new Epetra_BlockMap(SourceMap);
  TargetMap_ = new Epetra_BlockMap(TargetMap);

  assert( SourceMap_ != 0 );
  assert( TargetMap_ != 0 );
  
  SourceImporter_ = new Epetra_Import(*TargetMap_,*SourceMap_);
  TargetImporter_ = new Epetra_Import(*SourceMap_,*TargetMap_);

  assert( SourceImporter_ != 0 );
  assert( TargetImporter_ != 0 );
}

//=============================================================================

EpetraExt_Redistor::~EpetraExt_Redistor()
{
  if( SourceImporter_ ) delete SourceImporter_;
  if( TargetImporter_ ) delete TargetImporter_;
  if( SourceMap_ ) delete SourceMap_;
  if( TargetMap_ ) delete TargetMap_;
}

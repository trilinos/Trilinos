/*
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
*/

#ifndef _AMESOS_EPETRA_REDISTRIBUTOR_H_
#define _AMESOS_EPETRA_REDISTRIBUTOR_H_

#include "Amesos_ConfigDefs.h"
#include "Amesos_EpetraInterface.h"

class Amesos_EpetraRedistributor : public Amesos_EpetraInterface
{

public:
  Amesos_EpetraRedistributor(const Epetra_LinearProblem & Problem,
			     const AMESOS::Parameter::List &ParameterList);
  ~Amesos_EpetraRedistributor();

  int SetRedistributor(const int NumProcs);
  
  int CreateTargetMap();
  int CreateImportAndExport();
  
  inline int NumTargetProcs() const
  {
    return( NumTargetProcs_ );
  }
    
  inline Epetra_Import * ImportFromTarget() const
  {
    return( ImportFromTarget_ );
  }
  
  inline Epetra_Import * ImportToTarget() const
  {
    return( ImportToTarget_ );
  }

  inline Epetra_Map * TargetMap() const
  {
    return( TargetMap_ );
  }
  
  inline Epetra_BlockMap * TargetBlockMap() const
  {
    return( TargetBlockMap_ );
  }

  inline Epetra_MultiVector * TargetRHS() const
  {
    return( TargetRHS_ );
  }
  
  int CreateTargetRHS(int nrhs);
     
private:

  int NumTargetProcs_;
     
  Epetra_MultiVector * TargetRHS_;        // MS // MUMPS required rhs (and then solution)
                                          // MS // entirely stored on the host, also for
                                          // MS // ICNTL(18) = 3. 

  Epetra_Import * ImportToTarget_;      // MS // import from distributed to host 
  Epetra_Import * ImportFromTarget_;    // MS // (for rhs and sol) 

  Epetra_Map *TargetMap_ ;               //  Points to a target Map (unused if IsLocal == 1 ) 
  Epetra_BlockMap * TargetBlockMap_;
  
  bool IsTargetMapOK_;
  bool IsImportAndExportOK_;

  
};

#endif

/* Copyright (2001) Sandia Corportation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 *
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

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

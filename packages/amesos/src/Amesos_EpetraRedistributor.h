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
  Amesos_EpetraRedistributor(const Epetra_LinearProblem * LinearProblem);
  ~Amesos_EpetraRedistributor();

  int CreateSerialMap();
  int CreateImportAndExport();

  inline bool IsLocal() const
  {
    return( IsLocal_);
  }

  inline int SetIsLocal(bool flag) 
  {
    IsLocal_ = flag;
    return 0;
  }
  
  inline Epetra_Import * ImportFromProcZero() 
  {
    return( ImportFromProcZero_ );
  }
  
  inline Epetra_Import * ImportToProcZero() 
  {
    return( ImportToProcZero_ );
  }

  inline Epetra_Map * SerialMap() 
  {
    return( SerialMap_ );
  }
  
  inline Epetra_BlockMap * SerialBlockMap() 
  {
    return( SerialBlockMap_ );
  }

  inline Epetra_MultiVector * SerialRHS()
  {
    return( SerialRHS_ );
  }
  
  int CreateSerialRHS(int nrhs);
  
private:
  
   
  Epetra_MultiVector * SerialRHS_;        // MS // MUMPS required rhs (and then solution)
                                          // MS // entirely stored on the host, also for
                                          // MS // ICNTL(18) = 3. 

  Epetra_Import * ImportToProcZero_;      // MS // import from distributed to host 
  Epetra_Import * ImportFromProcZero_;    // MS // (for rhs and sol) 

  Epetra_Map *SerialMap_ ;               //  Points to a Serial Map (unused if IsLocal == 1 ) 
  Epetra_BlockMap * SerialBlockMap_;
  
  bool IsSerialMapOK_;
  bool IsImportAndExportOK_;
  bool IsLocal_;            //  1 if Problem_->GetOperator() is stored entirely on process 0
                           //  Note:  Local Problems do not require redistribution of
                           //  the matrix A or vectors X and B.
  
};

#endif

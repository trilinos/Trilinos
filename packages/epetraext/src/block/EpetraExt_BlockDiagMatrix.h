//@HEADER
/*
************************************************************************

              EpetraExt: Extended Linear Algebra Services Package 
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

#ifndef EPETRAEXT_BLOCKDIAGMATRIX_H
#define EPETRAEXT_BLOCKDIAGMATRIX_H

#include "Epetra_BLAS.h"
#include "Epetra_LAPACK.h"
#include "Epetra_DistObject.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Teuchos_ParameterList.hpp"

class Epetra_Comm;

//! EpetraExt_BlockDiagMatrix: A class for storing distributed block matrices.

/*!
  A dense-block block-diagonal matrix with inversion/factorization
  capabilities.

  This class has a rigid map structure --- the Domain and Range
  maps must be the same. 
*/

//=========================================================================
class EpetraExt_BlockDiagMatrix : public virtual Epetra_Operator, public Epetra_DistObject,public Epetra_BLAS {

 public:
  //! Constructor - This map is the map of the vector this can be applied to
  EpetraExt_BlockDiagMatrix(const Epetra_BlockMap& Map,bool zero_out=true);

  
  //! Copy constructor
  EpetraExt_BlockDiagMatrix(const EpetraExt_BlockDiagMatrix& Source);
    
  //! Destructor
  virtual ~EpetraExt_BlockDiagMatrix();

  //! = Operator.
  /*!
    \param In
           A - EpetraExt_BlockDiagMatrix to copy.

    \return EpetraExt_BlockDiagMatrix.
  */
  EpetraExt_BlockDiagMatrix& operator = (const EpetraExt_BlockDiagMatrix& Source);
  
  //! Block access function.
  /*!
    \return the pointer V[Index].
  */
  double* operator [] (int index) {return &Values_[DataMap_->FirstPointInElement(index)];}
  //! Block access function.
  /*!
    \return the pointer V[Index].
  */
  const double* operator [] (int index) const {return &Values_[DataMap_->FirstPointInElement(index)];}
  //@}


  //! @name Atribute set methods
  //@{ 

  //! SetUseTranspose - not implemented
  virtual int SetUseTranspose(bool UseTranspose){return -1;}

  //! SetParameters
  virtual int SetParameters(Teuchos::ParameterList & List);

  //! Computes the inverse / factorization if such is set on the list
  virtual int Compute();

  //@}

  
  //! @name Attribute access functions
  //@{

  //! Returns a character string describing the operator
  virtual const char * Label() const{return "EpetraExt::BlockDiagMatrix";}//HAQ
  
  //! Returns the current UseTranspose setting.
  virtual bool UseTranspose() const {return false;}

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  virtual bool HasNormInf() const {return false;}

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  virtual const Epetra_Comm & Comm() const {return Map().Comm();}

  //! Returns the Epetra_Map object associated with the domain of this operator.
  virtual const Epetra_Map & OperatorDomainMap() const {return *dynamic_cast<const Epetra_Map*>(&Map());}

  //! Returns the Epetra_Map object associated with the range of this operator.
  virtual const Epetra_Map & OperatorRangeMap() const {return *dynamic_cast<const Epetra_Map*>(&Map());}

  //! Returns the Epetra_BlockMap object associated with the range of this operator.
  virtual const Epetra_BlockMap & BlockMap() const {return Map();} 
  
  //! Returns a pointer to the array containing the blocks.
  double* Values() const {return(Values_);}

  //! Returns the size of the given block
  int BlockSize(int LID) const {return Map().ElementSize(LID);}

  //! Returns the size of the data in the given block
  int DataSize(int LID) const {return DataMap_->ElementSize(LID);}

  //! Returns true if the element size is constant
  bool ConstantBlockSize() const {return Map().ConstantElementSize();}
  
  //! Returns the number of local blocks
  int NumMyBlocks() const {return(Map().NumMyElements());}

  //! Returns the number of global blocks
  int NumGlobalBlocks() const {return(Map().NumGlobalElements());}
  
  //! Returns the number of local unknowns
  int NumMyUnknowns() const {return(Map().NumMyPoints());}

  //! Returns the number of global unknowns
  int NumGlobalUnknowns() const {return(Map().NumGlobalPoints());}

  //! Returns the size of the total Data block
  int NumData() const {return DataMap_->NumMyPoints();}
  
  //! Print method
  virtual void Print(ostream & os) const;
  //@}

  
  //! @name Mathematical functions
  //@{ 

  //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
  /*! 
    \param In
      X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
      Y -A Epetra_MultiVector of dimension NumVectors containing result.    
    \return Integer error code, set to 0 if successful.
  */
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {return -1;}
  
  //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
  /*! 
    \param In
      X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
      Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method must 
              support the case where X and Y are the same object.
  */
  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! NormInf - Not Implemented
  virtual double NormInf() const{return -1;}

  //! PutScalar function
  void PutScalar(double value);

  //@}
  
  
private:
  void Allocate();
  
  int DoCopy(const EpetraExt_BlockDiagMatrix& Source);

  // Routines to implement Epetra_DistObject virtual methods
  // Allows the source and target (\e this) objects to be compared for compatibility, return nonzero if not.
  int CheckSizes(const Epetra_SrcDistObject& Source);
  // Perform ID copies and permutations that are on processor.
  int CopyAndPermute(const Epetra_SrcDistObject& Source,
                     int NumSameIDs, 
                     int NumPermuteIDs,
                     int * PermuteToLIDs,
                     int * PermuteFromLIDs,
                     const Epetra_OffsetIndex * Indexor);

  // Perform any packing or preparation required for call to DoTransfer().
  int PackAndPrepare(const Epetra_SrcDistObject& Source,
                     int NumExportIDs,
                     int* ExportLIDs,
                     int& LenExports,
                     char*& Exports,
                     int& SizeOfPacket,
                     int* Sizes,
                     bool & VarSizes,
                     Epetra_Distributor& Distor);
  
  // Perform any unpacking and combining after call to DoTransfer().
  int UnpackAndCombine(const Epetra_SrcDistObject& Source, 
                       int NumImportIDs,
                       int* ImportLIDs, 
                       int LenImports,
                       char* Imports,
                       int& SizeOfPacket, 
                       Epetra_Distributor& Distor,
                       Epetra_CombineMode CombineMode,
                       const Epetra_OffsetIndex * Indexor);

  Epetra_LAPACK LAPACK;
  Teuchos::ParameterList List_;

  //! Has Computed? Needed for Inverse/Factorization modes
  bool HasComputed_;
  
  //! Which Apply Mode to use
  int ApplyMode_;

  //! Map for the data
  Epetra_BlockMap* DataMap_;
  
  //! Actual Data values
  double *Values_;

  //! Pivots for factorization
  int *Pivots_;

};  /* EPETRAEXT_BLOCKDIAGMATRIX_H */

#endif

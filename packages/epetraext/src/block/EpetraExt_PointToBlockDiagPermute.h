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

#ifndef EPETRAEXT_POINTTOBLOCKDIAGPERMUTE_H
#define EPETRAEXT_POINTTOBLOCKDIAGPERMUTE_H

#include "Epetra_DistObject.h"
#include "Epetra_BlockMap.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_ParameterList.hpp"

class Epetra_Comm;
class Epetra_Import;
class Epetra_Export;
class EpetraExt_BlockDiagMatrix;

//! EpetraExt_PointToBlockDiagPermute: A class for managing point-to-block-diagonal permutations

/*!
  Manages point-to-block-diagonal permutations (and vice versa) as well as service
  routines for extracting block diagonals from matrices.

*/


// NTS: Map() == RowMap is the convention

//=========================================================================
class EpetraExt_PointToBlockDiagPermute : public virtual Epetra_Operator, public Epetra_DistObject {
public:

  //! @ Name Constructors
  //@{
  //! Constructor
  EpetraExt_PointToBlockDiagPermute(const Epetra_CrsMatrix& MAT);
  //@}
  
  //! @name Destructor
  //@{ 
  //! Destructor
  virtual ~EpetraExt_PointToBlockDiagPermute();
  //@}

  
  //! @name Atribute set methods
  //@{ 

  //! Sets the parameter list
  virtual int SetParameters(Teuchos::ParameterList & List);
  
  //! If set true, transpose of this operator will be applied. NOT IMPLEMENTED.
  virtual int SetUseTranspose(bool UseTranspose){return -1;}

  //! Extracts the block-diagonal, builds maps, etc.
  virtual int Compute();
  
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
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

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
  
  //! Returns the infinity norm of the global matrix - NOT IMPLEMENTED
  virtual double NormInf() const {return -1;}
  //@}
  
  //! @name Atribute access functions
  //@{ 

  //! Returns a character string describing the operator
  virtual const char * Label() const{return "Fix Me";}

  //! Returns the current UseTranspose setting.
  virtual bool UseTranspose() const {return false;}

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  virtual bool HasNormInf() const {return false;}

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  virtual const Epetra_Comm & Comm() const {return Map().Comm();}
  
  //! Returns the Epetra_Map object associated with the domain of this
  //operator. 
  virtual const Epetra_Map & OperatorDomainMap() const {return Matrix_->OperatorDomainMap();}

  //! Returns the Epetra_Map object associated with the range of this
  //operator.
  virtual const Epetra_Map & OperatorRangeMap() const {return Matrix_->OperatorRangeMap();}

  
  //! Returns the block matrix. Only call this after calling Compute
  virtual const EpetraExt_BlockDiagMatrix & BlockMatrix(){return *BDMat_;}
  

  //! @name Miscellaneous
  //@{ 
  //! Print method
  virtual void Print(ostream& os) const;
  //@}


  //! @name Import/Export Methods
  //@{ 

  //! Imports an Epetra_DistObject using the Epetra_Import object.
  /*!
    \param In
           Source - Distributed object that will be imported into the "\e this" object.
    \param In
           Importer - A Epetra_Import object specifying the communication required.

    \param In
           CombineMode - A Epetra_CombineMode enumerated type specifying how results should be combined on the 
	   receiving processor.

    \return Integer error code, set to 0 if successful.
  */
  int Import(const Epetra_SrcDistObject& A, const Epetra_Import& Importer, Epetra_CombineMode CombineMode, const Epetra_OffsetIndex * Indexor = 0);

  //! Imports an Epetra_DistObject using the Epetra_Export object.
  /*!
    \param In
           Source - Distributed object that will be imported into the "\e this" object.
    \param In
           Exporter - A Epetra_Export object specifying the communication required.

    \param In
           CombineMode - A Epetra_CombineMode enumerated type specifying how results should be combined on the 
	   receiving processor.

    \return Integer error code, set to 0 if successful.
  */
  int Import(const Epetra_SrcDistObject& A, const Epetra_Export& Exporter, Epetra_CombineMode CombineMode, const Epetra_OffsetIndex * Indexor = 0);

  //! Exports an Epetra_DistObject using the Epetra_Import object.
  /*!
    \param In
           Source - Distributed object that will be exported to the "\e this" object.
    \param In
           Importer - A Epetra_Import object specifying the communication required.

    \param In
           CombineMode - A Epetra_CombineMode enumerated type specifying how results should be combined on the 
	   receiving processor.

    \return Integer error code, set to 0 if successful.
  */
  int Export(const Epetra_SrcDistObject& A, const Epetra_Import & Importer, Epetra_CombineMode CombineMode, const Epetra_OffsetIndex * Indexor = 0);

  //! Exports an Epetra_DistObject using the Epetra_Export object.
  /*!
    \param In
           Source - Distributed object that will be exported to the "\e this" multivector.
    \param In
           Exporter - A Epetra_Export object specifying the communication required.

    \param In
           CombineMode - A Epetra_CombineMode enumerated type specifying how results should be combined on the 
	   receiving processor.

    \return Integer error code, set to 0 if successful.
  */
  int Export(const Epetra_SrcDistObject& A, const Epetra_Export& Exporter, Epetra_CombineMode CombineMode, const Epetra_OffsetIndex * Indexor = 0);
  //@}
  
 protected:


  //! @name Internal utilities  
  //@{ 
  //! Allows the source and target (\e this) objects to be compared for compatibility, return nonzero if not.
  virtual int CheckSizes(const Epetra_SrcDistObject& Source);
  //! Perform ID copies and permutations that are on processor.
  virtual int CopyAndPermute(const Epetra_SrcDistObject& Source,
                             int NumSameIDs, 
                             int NumPermuteIDs,
                             int * PermuteToLIDs,
                             int * PermuteFromLIDs,
                             const Epetra_OffsetIndex * Indexor);

  //! Perform any packing or preparation required for call to DoTransfer().
  virtual int PackAndPrepare(const Epetra_SrcDistObject& Source,
                             int NumExportIDs,
                             int* ExportLIDs,
                             int& LenExports,
                             char*& Exports,
                             int& SizeOfPacket,
                             int* Sizes,
                             bool & VarSizes,
                             Epetra_Distributor& Distor);
  
  //! Perform any unpacking and combining after call to DoTransfer().
  virtual int UnpackAndCombine(const Epetra_SrcDistObject& Source, 
                               int NumImportIDs,
                               int* ImportLIDs, 
                               int LenImports,
                               char* Imports,
                               int& SizeOfPacket, 
                               Epetra_Distributor& Distor,
                               Epetra_CombineMode CombineMode,
                               const Epetra_OffsetIndex * Indexor);

  //@}

private:
  //! Pulls the block diagonal of the matrix and then builds the BDMat_
  int ExtractBlockDiagonal();

  // Copied from Epetra_CrsMatrix
  void UpdateImportVector(int NumVectors) const;
  void UpdateExportVector(int NumVectors) const;

  Teuchos::ParameterList List_;
  const Epetra_CrsMatrix* Matrix_;  
  bool PurelyLocalMode_;  
  int NumBlocks_;
  int *Blockstart_;
  int *Blockids_;
  Epetra_BlockMap *BDMap_;
  Epetra_Map *CompatibleMap_; //A map compatible with BD's block map - used for imports
  EpetraExt_BlockDiagMatrix* BDMat_;
  Epetra_Import *Importer_;
  Epetra_Export *Exporter_;
  mutable Epetra_MultiVector *ImportVector_;
  mutable Epetra_MultiVector *ExportVector_;
  

};  /* EPETRAEXT_POINTTOBLOCKDIAGPERMUTE_H */

#endif

//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef EPETRAEXT_BLOCKDIAGMATRIX_H
#define EPETRAEXT_BLOCKDIAGMATRIX_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

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


  //! @name Attribute set methods
  //@{

  //! SetUseTranspose - not implemented
  virtual int SetUseTranspose(bool /* useTranspose */){return -1;}

  //! SetParameters
  virtual int SetParameters(Teuchos::ParameterList & List);

  //! Computes the inverse / factorization if such is set on the list
  virtual int Compute();

  //@}


  //! @name Attribute access functions
  //@{

  //! Returns a character std::string describing the operator
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
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int NumGlobalBlocks() const {return(Map().NumGlobalElements());}
#endif
  long long NumGlobalBlocks64() const {return(Map().NumGlobalElements64());}

  //! Returns the number of local unknowns
  int NumMyUnknowns() const {return(Map().NumMyPoints());}

  //! Returns the number of global unknowns
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int NumGlobalUnknowns() const {return(Map().NumGlobalPoints());}
#endif
  long long NumGlobalUnknowns64() const {return(Map().NumGlobalPoints64());}

  //! Returns the size of the total Data block
  int NumData() const {return DataMap_->NumMyPoints();}

  //! Gets apply mode info
  int GetApplyMode() {return ApplyMode_;}

  //! Print method
  virtual void Print(std::ostream & os) const;

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
  virtual int Apply(const Epetra_MultiVector& /* X */, Epetra_MultiVector& /* Y */) const {return -1;}

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

  //! Returns the Epetra_BlockMap object with the distribution of underlying values
  virtual const Epetra_BlockMap & DataMap() const {return *DataMap_;}

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
                     const Epetra_OffsetIndex * Indexor,
                     Epetra_CombineMode CombineMode = Zero);

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

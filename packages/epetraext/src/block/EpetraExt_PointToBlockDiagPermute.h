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

#ifndef EPETRAEXT_POINTTOBLOCKDIAGPERMUTE_H
#define EPETRAEXT_POINTTOBLOCKDIAGPERMUTE_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

#include "Epetra_ConfigDefs.h"
#include "Epetra_DistObject.h"
#include "Epetra_BlockMap.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
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


  //! @name Attribute set methods
  //@{

  //! Sets the parameter list
  virtual int SetParameters(Teuchos::ParameterList & List);

  /// \brief Set whether to use the transpose.
  ///
  /// \return 0 if this class can change its transpose state, else nonzero.
  ///
  /// This class does NOT know how to apply its transpose, so this
  /// method always returns an error code.
  virtual int SetUseTranspose (bool /* useTranspose */) {return -1;}

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

  //! @name Attribute access functions
  //@{

  //! Return a string describing the operator
  virtual const char * Label() const{return "Fix Me";}

  //! Return the current UseTranspose setting.
  virtual bool UseTranspose() const {return false;}

  //! Return true if the \e this object can provide an approximate Inf-norm, false otherwise.
  virtual bool HasNormInf() const {return false;}

  //! Return a pointer to the Epetra_Comm communicator associated with this operator.
  virtual const Epetra_Comm & Comm() const {return Map().Comm();}

  //! Return the Epetra_Map object associated with the domain of this operator.
  virtual const Epetra_Map & OperatorDomainMap() const {return Matrix_->OperatorDomainMap();}

  //! Return the Epetra_Map object associated with the range of this operator.
  virtual const Epetra_Map & OperatorRangeMap() const {return Matrix_->OperatorRangeMap();}

  //! Return the block matrix. Only call this after calling Compute.
  virtual const EpetraExt_BlockDiagMatrix & BlockMatrix(){return *BDMat_;}

  /// \brief Create an Epetra_FECrsMatrix from the BlockDiagMatrix.
  ///
  /// This is generally only useful if you want to do a matrix-matrix multiply.
  virtual Epetra_FECrsMatrix * CreateFECrsMatrix();

  //! @name Miscellaneous
  //@{

  //! Print information about this object to the given output stream.
  virtual void Print(std::ostream& os) const;

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
                             const Epetra_OffsetIndex * Indexor,
                             Epetra_CombineMode CombineMode = Zero);

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

  //! Setup/Cleanup for Contiguous Mode
  int SetupContiguousMode();
  int CleanupContiguousMode();

  // Copied from Epetra_CrsMatrix
  void UpdateImportVector(int NumVectors) const;
  void UpdateExportVector(int NumVectors) const;

  Teuchos::ParameterList List_;
  const Epetra_CrsMatrix* Matrix_;
  bool PurelyLocalMode_;

  // For contiguous blocking only
  bool ContiguousBlockMode_;
  int ContiguousBlockSize_;

  int NumBlocks_;
  int *Blockstart_;
  int *Blockids_int_;
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  long long *Blockids_LL_;
#endif
  Epetra_BlockMap *BDMap_;
  Epetra_Map *CompatibleMap_; //A map compatible with BD's block map - used for imports
  EpetraExt_BlockDiagMatrix* BDMat_;
  Epetra_Import *Importer_;
  Epetra_Export *Exporter_;
  mutable Epetra_MultiVector *ImportVector_;
  mutable Epetra_MultiVector *ExportVector_;

  template<typename int_type>
  Epetra_FECrsMatrix * TCreateFECrsMatrix();

  template<typename int_type>
  int TSetupContiguousMode();

  template<typename int_type>
  int TExtractBlockDiagonal();

  template<typename int_type>
  int TSetParameters(Teuchos::ParameterList & List);

  template<typename int_type>
  const int_type* Blockids_const_ptr() const;

  template<typename int_type>
  int_type*& Blockids_ref();

};  /* EPETRAEXT_POINTTOBLOCKDIAGPERMUTE_H */

template<> inline const int* EpetraExt_PointToBlockDiagPermute::Blockids_const_ptr<int>() const { return Blockids_int_; }

template<> inline int*& EpetraExt_PointToBlockDiagPermute::Blockids_ref<int>() { return Blockids_int_; }

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
template<> inline const long long* EpetraExt_PointToBlockDiagPermute::Blockids_const_ptr<long long>() const { return Blockids_LL_; }
template<> inline long long*& EpetraExt_PointToBlockDiagPermute::Blockids_ref<long long>() { return Blockids_LL_; }
#endif

#endif

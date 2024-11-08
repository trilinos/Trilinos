//@HEADER
// ************************************************************************
// 
//              Ifpack_SubdomainFilter
//              Author: Radu Popescu <radu.popescu@epfl.ch>
//              Copyright 2011 EPFL - CMCS
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
// 3. Neither the name of the copyright holder nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY EPFL - CMCS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL EPFL - CMCS OR THE CONTRIBUTORS
// BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
// OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
// BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
// OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
//@HEADER

#ifndef IFPACK_SUBDOMAINFILTER_H
#define IFPACK_SUBDOMAINFILTER_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Ifpack_ConfigDefs.h"

#ifdef HAVE_IFPACK_PARALLEL_SUBDOMAIN_SOLVERS

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_IntVector.h"
#include "Teuchos_RCP.hpp"
#include "Ifpack_OverlappingRowMatrix.h"

class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_Import;
class Epetra_BlockMap;

//! Ifpack_SubdomainFilter a class for the extraction of a distributed subdomain matrix.

/*! 
    This class is used by Ifpack_AdditiveSchwarz when more than 1 MPI process per
    subdomain are desired.

    It implements an Epetra_RowMatrix interface. Given a global distributed matrix,
    this class creates an Epetra MPI subcommunicator on a selected number of processes
    and extracts from the global matrix the rows that belong to the processes in the
    subcommunicator. In these rows, only the column indices that coincide with any
    row indices in the subcommunicator are kept (only couplings between the DOFs of
    processes in the same subdomain are considered, while the couplings between
    subdomains are discarded).

    TODO: currently this class doesn't work with overlapping subdomains in Ifpack
 */ 
class Ifpack_SubdomainFilter : public virtual Epetra_RowMatrix {

public:
  //@{ \name Constructor.
  //! Constructor
  Ifpack_SubdomainFilter(const Teuchos::RCP<const Epetra_RowMatrix>& Matrix, int subdomainId);

  //@}
  //@{ \name Destructor.
  //! Destructor
  ~Ifpack_SubdomainFilter();
  //@}

  //@{ \name Matrix data extraction routines

  //! Returns the number of nonzero entries in MyRow.
  /*! 
    \param 
    MyRow - (In) Local row.
    \param 
    NumEntries - (Out) Number of nonzero values present.

    \return Integer error code, set to 0 if successful.
    */
  virtual int NumMyRowEntries(int MyRow, int & NumEntries) const
  {
    NumEntries = NumEntries_[MyRow];
    return(0);
  }

  //! Returns the maximum of NumMyRowEntries() over all rows.
  virtual int MaxNumEntries() const
  {
    return(MaxNumEntries_);
  }

  //! Returns a copy of the specified local row in user-provided arrays.
  /*! 
    \param
    MyRow - (In) Local row to extract.
    \param
    Length - (In) Length of Values and Indices.
    \param
    NumEntries - (Out) Number of nonzero entries extracted.
    \param
    Values - (Out) Extracted values for this row.
    \param 
    Indices - (Out) Extracted global column indices for the corresponding values.

    \return Integer error code, set to 0 if successful.
    */
  virtual inline int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;

  //! Returns a copy of the main diagonal in a user-provided vector.
  /*! 
    \param
    Diagonal - (Out) Extracted main diagonal.

    \return Integer error code, set to 0 if successful.
    */
  virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const;
  //@}

  //@{ \name Mathematical functions.

  //! Returns the result of a Epetra_RowMatrix multiplied by a Epetra_MultiVector X in Y.
  /*! 
    \param 
    TransA -(In) If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param 
    X - (In) A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param 
    Y -(Out) A Epetra_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
    */
  virtual int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
  {
    if (TransA == true) {
      IFPACK_CHK_ERR(-1);
    }

    IFPACK_CHK_ERR(Apply(X,Y));
    return(0);
  }

  //! Returns result of a local-only solve using a triangular Epetra_RowMatrix with Epetra_MultiVectors X and Y (NOT IMPLEMENTED).
  virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, 
		    Epetra_MultiVector& Y) const
  {
    IFPACK_RETURN(-1); // not implemented 
  }

  virtual int Apply(const Epetra_MultiVector& X,
		    Epetra_MultiVector& Y) const;

  virtual int ApplyInverse(const Epetra_MultiVector& X,
			   Epetra_MultiVector& Y) const;
  //! Computes the sum of absolute values of the rows of the Epetra_RowMatrix, results returned in x (NOT IMPLEMENTED).
  virtual int InvRowSums(Epetra_Vector& x) const
  {
    IFPACK_RETURN(-1); // not implemented
  }

  //! Scales the Epetra_RowMatrix on the left with a Epetra_Vector x (NOT IMPLEMENTED).
  virtual int LeftScale(const Epetra_Vector& x)
  {
    IFPACK_RETURN(-1); // not implemented
  }

  //! Computes the sum of absolute values of the columns of the Epetra_RowMatrix, results returned in x (NOT IMPLEMENTED).
  virtual int InvColSums(Epetra_Vector& x) const
  {
    IFPACK_RETURN(-1); // not implemented
  }


  //! Scales the Epetra_RowMatrix on the right with a Epetra_Vector x (NOT IMPLEMENTED).
  virtual int RightScale(const Epetra_Vector& x) 
  {
    IFPACK_RETURN(-1); // not implemented
  }

  //@}

  //@{ \name Attribute access functions

  //! If FillComplete() has been called, this query returns true, otherwise it returns false.
  virtual bool Filled() const
  {
    return true;
  }

  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\len} \sum_{i=1}^m |a_{ij}| \f].
     */ 
  virtual double NormInf() const
  {
    return(-1.0);
  }

  //! Returns the one norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_1\f$ such that
     \f[\| A \|_1= \max_{1\lej\len} \sum_{j=1}^n |a_{ij}| \f].
     */ 
  virtual double NormOne() const
  {
    IFPACK_RETURN(-1.0);
  }

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  //! Returns the number of nonzero entries in the global matrix.
  virtual int NumGlobalNonzeros() const
  {
    return(NumGlobalNonzeros_);
  }

  //! Returns the number of global matrix rows.
  virtual int NumGlobalRows() const
  {;

    return(NumGlobalRows_);
  }

  //! Returns the number of global matrix columns.
  virtual int NumGlobalCols() const
  {
    return(NumGlobalCols_);
  }

  //! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons.
  virtual int NumGlobalDiagonals() const
  {
    return(NumGlobalRows_);
  }
#endif

  //! Returns the number of nonzero entries in the global matrix.
  virtual long long NumGlobalNonzeros64() const
  {
    return(NumGlobalNonzeros_);
  }

  //! Returns the number of global matrix rows.
  virtual long long NumGlobalRows64() const
  {
    return(NumGlobalRows_);
  }

  //! Returns the number of global matrix columns.
  virtual long long NumGlobalCols64() const
  {
    return(NumGlobalRows_);
  }

  //! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons.
  virtual long long NumGlobalDiagonals64() const
  {
    return(NumGlobalRows_);
  }

  //! Returns the number of nonzero entries in the calling processor's portion of the matrix.
  virtual int NumMyNonzeros() const
  {
    return(NumMyNonzeros_);
  }

  //! Returns the number of matrix rows owned by the calling processor.
  virtual int NumMyRows() const
  {
    return(NumMyRows_);
  }

  //! Returns the number of matrix columns owned by the calling processor.
  virtual int NumMyCols() const
  {
    return(NumMyCols_);
  }

  //! Returns the number of local nonzero diagonal entries, based on global row/column index comparisons.
  virtual int NumMyDiagonals() const
  {
    return(NumMyRows_);
  }

  //! If matrix is lower triangular in local index space, this query returns true, otherwise it returns false.
  virtual bool LowerTriangular() const
  {
    return(Matrix_->LowerTriangular());
  }

  //! If matrix is upper triangular in local index space, this query returns true, otherwise it returns false.
  virtual bool UpperTriangular() const
  {
    return(Matrix_->UpperTriangular());
  }

  //! Returns the Epetra_Map object associated with the rows of this matrix.
  virtual const Epetra_Map & RowMatrixRowMap() const
  {
    return(*Map_);
  }

  //! Returns the Epetra_Map object associated with the columns of this matrix.
  virtual const Epetra_Map & RowMatrixColMap() const
  {
    return(*colMap_);
  }

  //! Returns the Epetra_Import object that contains the import operations for distributed operations.
  virtual const Epetra_Import * RowMatrixImporter() const
  {
    return(&*Importer_);
  }
  //@}

  virtual const Epetra_Import* Importer() const {return(&*Importer_);}

  virtual const Epetra_Export* Exporter() const {return(&*Exporter_);}

  // following functions are required to derive Epetra_RowMatrix objects.

  //! Sets ownership.
  int SetOwnership(bool ownership)
  {
    IFPACK_RETURN(-1);
  }

  //! Sets use transpose (not implemented).
  int SetUseTranspose(bool UseTranspose_in)
  {
    UseTranspose_ = UseTranspose_in;
    return(0);
  }

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const 
  {
    return(UseTranspose_);
  }

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  bool HasNormInf() const
  {
    return(false);
  }

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm & Comm() const
  {
    return(*SubComm_);
  }

  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & OperatorDomainMap() const 
  {
    return(*Map_);
   }

  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & OperatorRangeMap() const 
  {
    return(*Map_);
  }
  //@}

const Epetra_BlockMap& Map() const;

const char* Label() const{
  return(Label_);
};

private:
  void UpdateImportVector(int NumVectors) const;
  void UpdateExportVector(int NumVectors) const;

  //! Pointer to the matrix to be preconditioned.
  Teuchos::RCP<const Epetra_RowMatrix> Matrix_;
  const Ifpack_OverlappingRowMatrix* ovA_;
#ifdef HAVE_MPI
  //! Communicator containing this process only.
  Teuchos::RCP<Epetra_MpiComm> SubComm_;
  MPI_Comm subdomainMPIComm_;
#else
  //! Communicator containing this process only.
  Teuchos::RCP<Epetra_SerialComm> SubComm_;
#endif
  //! Row, domain, and range map, based on SubComm_, containing the local rows only.
  Teuchos::RCP<Epetra_Map> Map_;
  //! Column map based on SubComm_, containing the local rows only.
  Teuchos::RCP<Epetra_Map> colMap_;
  //! Number of rows in the local matrix.
  int NumMyRows_;
  //! Number of cols in the local matrix.
  int NumMyCols_;
  //! Number of nonzeros in the local matrix.
  int NumMyNonzeros_;
  //! Number of rows in the global filtered matrix.
  int NumGlobalRows_;
  //! Number of columns in the global filtered matrix.
  int NumGlobalCols_;
  //! Number of nonzeros in the global filtered matrix.
  int NumGlobalNonzeros_;
  //! Maximum number of nonzero entries in a row for the filtered matrix.
  int MaxNumEntries_;
  //! Maximum number of nonzero entries in a row for Matrix_.
  int MaxNumEntriesA_;
  //! NumEntries_[i] contains the nonzero entries in row `i'.
  std::vector<int> NumEntries_;
  //! Used in ExtractMyRowCopy, to avoid allocation each time.
  mutable std::vector<int> Indices_;
  //! Used in ExtractMyRowCopy, to avoid allocation each time.
  mutable std::vector<double> Values_;
  //! If true, the tranpose of the local matrix will be used.
  bool UseTranspose_;
  //! Label for \c this object.
  char Label_[80];
  Teuchos::RCP<Epetra_Vector> Diagonal_;
  double NormOne_;
  double NormInf_;

  //! Maps to speed LID-LID conversions
  int* Ac_LIDMap_;
  int* Bc_LIDMap_;
  int* Ar_LIDMap_;
  int* Br_LIDMap_;

  //! CrsMatrix pointer, if needed
  const Epetra_CrsMatrix* Acrs_;

  int NumMyRowsA_;
  int NumMyColsA_;
  int NumMyRowsB_;

  //mutable Teuchos::RCP<Epetra_MultiVector> ImportVector_;
  //mutable Teuchos::RCP<Epetra_MultiVector> ExportVector_;
  mutable Epetra_MultiVector* ExportVector_;
  mutable Epetra_MultiVector* ImportVector_;
  Teuchos::RCP<Epetra_Import> Importer_;
  Teuchos::RCP<Epetra_Export> Exporter_;

};
#endif //ifdef HAVE_IFPACK_PARALLEL_SUBDOMAIN_SOLVERS
#endif /* IFPACK_SUBDOMAINFILTER_H */

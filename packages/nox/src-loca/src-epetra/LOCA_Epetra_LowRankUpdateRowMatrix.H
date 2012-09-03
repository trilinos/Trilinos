// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
                                                                                
#ifndef LOCA_EPETRA_LOWRANKUPDATEROWMATRIX_H
#define LOCA_EPETRA_LOWRANKUPDATEROWMATRIX_H

#include "Teuchos_RCP.hpp"

#include "LOCA_Epetra_LowRankUpdateOp.H" // base class
#include "Epetra_RowMatrix.h"	    // base class

namespace LOCA {

  namespace Epetra {
    
    /*! 
     * \brief An Epetra row matrix for implementing the operator
     * \f$P = J + U V^T\f$.
     */
    /*!
     * This class implements the Epetra_RowMatrix interface for 
     * \f$P = J + U V^T\f$ where \f$J\f$ is an Epetra_RowMatrix and 
     * \f$U\f$ and \f$V\f$ are Epetra_MultiVectors.  It is derived
     * from LOCA::Epetra::LowRankUpdateOp to implement the Epetra_Operator
     * interface.  The interface here implements the Epetra_RowMatrix
     * interface when the matrix \f$J\f$ is itself a row matrix.  This 
     * allows preconditioners to be computed and scaling in linear systems
     * to be performed when using this operator.  The implementation
     * here merely adds the corresponding entries for \f$U V^T\f$ to the
     * rows of \f$J\f$.  Note however this is only an approximation to the
     * true matrix \f$J + U V^T\f$.
     *
     * This class assumes \f$U\f$ and \f$V\f$ have the same distribution
     * as the rows of \f$J\f$.
     */
    class LowRankUpdateRowMatrix : 
      public LOCA::Epetra::LowRankUpdateOp,
      public virtual Epetra_RowMatrix {
      
    public:

      //! Constructor 
      /*!
       * \param global_data [in] The global data object
       * \param jacRowMatrix [in] Jacobian operator J as a row matrix
       * \param U_multiVec [in] Multivector representing U
       * \param V_multiVec [in] Multivector representing V
       * \param setup_for_solve [in] Setup data structures for ApplyInverse()
       * \param include_UV_terms [in] Include \f$U V^T\f$ terms in RowMatrix
       *        routines ExtractRowCopy(), ExtactDiagonalCopy(), InvRowSums(),
       *        InvColSums(), NormInf() and NormOne().
       */
      LowRankUpdateRowMatrix(
        const Teuchos::RCP<LOCA::GlobalData>& global_data,
	const Teuchos::RCP<Epetra_RowMatrix>& jacRowMatrix, 
	const Teuchos::RCP<Epetra_MultiVector>& U_multiVec, 
	const Teuchos::RCP<Epetra_MultiVector>& V_multiVec,
	bool setup_for_solve,
	bool include_UV_terms);

      //! Destructor
      virtual ~LowRankUpdateRowMatrix();

      //! Returns a reference to the Epetra_BlockMap for this object.
      virtual const Epetra_BlockMap & Map() const;

      //@{ \name Matrix data extraction routines

      //! Returns the number of nonzero entries in MyRow.
      virtual int NumMyRowEntries(int MyRow, int & NumEntries) const;


      //! Returns the maximum of NumMyRowEntries() over all rows.
      virtual int MaxNumEntries() const;

      //! Returns a copy of the specified local row in user-provided arrays.
      virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
				   double *Values, int * Indices) const;

      //! Returns a copy of the main diagonal in a user-provided vector.
      virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const;

      //@}
  
      //@{ \name Mathematical functions.

      /*! 
       * \brief Returns the result of a Epetra_RowMatrix multiplied by a 
       * Epetra_MultiVector X in Y.
       */
      virtual int Multiply(bool TransA, const Epetra_MultiVector& X, 
			   Epetra_MultiVector& Y) const;

      /*! 
       * \brief Returns result of a local-only solve using a triangular 
       * Epetra_RowMatrix with Epetra_MultiVectors X and Y.
       */
      virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, 
			const Epetra_MultiVector& X, 
			Epetra_MultiVector& Y) const;

      /*! 
       * \brief Computes the sum of absolute values of the rows of the 
       * Epetra_RowMatrix, results returned in x.
       */
      virtual int InvRowSums(Epetra_Vector& x) const;

      //! Scales the Epetra_RowMatrix on the left with a Epetra_Vector x.
      virtual int LeftScale(const Epetra_Vector& x);

      /*! 
       * \brief Computes the sum of absolute values of the columns of the 
       * Epetra_RowMatrix, results returned in x.
       */
      virtual int InvColSums(Epetra_Vector& x) const;

      //! Scales the Epetra_RowMatrix on the right with a Epetra_Vector x.
      virtual int RightScale(const Epetra_Vector& x);

      //@}
  
      //@{ \name Attribute access functions
      
      /*! 
       * \brief If FillComplete() has been called, this query returns true, 
       * otherwise it returns false.
       */
      virtual bool Filled() const;

      //! Returns the infinity norm of the global matrix.
      virtual double NormInf() const;

      //! Returns the one norm of the global matrix.
      virtual double NormOne() const;

      //! Returns the number of nonzero entries in the global matrix.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
      virtual int NumGlobalNonzeros() const;
#endif
      virtual long long NumGlobalNonzeros64() const;

      //! Returns the number of global matrix rows.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
      virtual int NumGlobalRows() const;
#endif
      virtual long long NumGlobalRows64() const;

      //! Returns the number of global matrix columns.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
      virtual int NumGlobalCols() const;
#endif
      virtual long long NumGlobalCols64() const;

      /*! 
       * \brief Returns the number of global nonzero diagonal entries, 
       * based on global row/column index comparisons.
       */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
      virtual int NumGlobalDiagonals() const;
#endif
      virtual long long NumGlobalDiagonals64() const;
    
      /*! 
       * \brief Returns the number of nonzero entries in the calling 
       * processor's portion of the matrix.
       */
      virtual int NumMyNonzeros() const;

      //! Returns the number of matrix rows owned by the calling processor.
      virtual int NumMyRows() const;

      //! Returns the number of matrix columns owned by the calling processor.
      virtual int NumMyCols() const;

      /*! 
       * \brief Returns the number of local nonzero diagonal entries, based 
       * on global row/column index comparisons.
       */
      virtual int NumMyDiagonals() const;

      /*! 
       * \brief If matrix is lower triangular in local index space, 
       * this query returns true, otherwise it returns false.
       */
      virtual bool LowerTriangular() const;

      /*! 
       * \brief If matrix is upper triangular in local index space, this 
       * query returns true, otherwise it returns false.
       */
      virtual bool UpperTriangular() const;

      /*! 
       * \brief Returns the Epetra_Map object associated with the rows of 
       * this matrix.
       */
      virtual const Epetra_Map & RowMatrixRowMap() const;

      /*! 
       * \brief Returns the Epetra_Map object associated with the columns of 
       * this matrix.
       */
      virtual const Epetra_Map & RowMatrixColMap() const;

      /*! 
       * \brief Returns the Epetra_Import object that contains the import 
       * operations for distributed operations.
       */
      virtual const Epetra_Import * RowMatrixImporter() const;

      //@}

    protected:

      //! Compute \c MyRow, \c MyCol entry of \f$U V^T\f$.
      double computeUV(int MyRow, int MyCol) const;

    private:

      //! Private to prohibit copying
      LowRankUpdateRowMatrix(const LowRankUpdateRowMatrix&);

      //! Private to prohibit copying
      LowRankUpdateRowMatrix& operator=(const LowRankUpdateRowMatrix&);

    protected:
    
      //! Stores row matrix representing J
      Teuchos::RCP<Epetra_RowMatrix> J_rowMatrix;

      //! Stores pointer to non-const U
      Teuchos::RCP<Epetra_MultiVector> nonconst_U;

      //! Stores pointer to non-const V
      Teuchos::RCP<Epetra_MultiVector> nonconst_V;

      //! Flag indicating whether to include U*V^T terms
      bool includeUV;

      //! Number of columns in U and V
      int m;

      //! Map for U
      const Epetra_BlockMap& U_map;

      //! Map for V
      const Epetra_BlockMap& V_map;

      //! Row map for J
      const Epetra_BlockMap& row_map;

    };
  } // namespace Epetra
} // namespace LOCA

#endif 

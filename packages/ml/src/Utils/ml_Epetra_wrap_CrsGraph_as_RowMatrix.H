/*
#@HEADER
# ************************************************************************
#
#               ML: A Multilevel Preconditioner Package
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#if defined(HAVE_ML_EPETRA)

#ifndef ML_EPETRA_CRSGRAPHWRAP_H
#define ML_EPETRA_CRSGRAPHWRAP_H

#include "Epetra_Import.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_RowMatrix.h"
class Epetra_Map;
class Epetra_Comm;
class Epetra_Import;
class Epetra_Export;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_SrcDistObject;

namespace ML_Epetra 
{

class CrsGraphWrapper: public virtual Epetra_RowMatrix 
{
      
 public:
  //@{ \name Constructor.
    //! Constructor
    CrsGraphWrapper(const Epetra_CrsGraph& graph, Epetra_Map& dm,Epetra_Map& rm,Epetra_Comm& comm) : 
    graph_(graph),
    DomainMap_(dm),
    RangeMap_(rm),
    comm_(comm)
    {
      name_  = "ML_Epetra::CrsGraphWrapper";
      return;
    }

  //@{ \name Destructor.
    //! Destructor
    virtual ~CrsGraphWrapper() {};

  //@}
  
  //@{ \name Matrix data extraction routines

    //! Returns the number of nonzero entries in MyRow.
    /*! 
    \param In
           MyRow - Local row.
    \param Out
	   NumEntries - Number of nonzero values present.
	  
    \return Integer error code, set to 0 if successful.
  */
    int NumMyRowEntries(int MyRow, int & NumEntries) const
    {
      NumEntries = graph_.NumMyIndices(MyRow);
      return(0);
    }


    //! Returns the maximum of NumMyRowEntries() over all rows.
    int MaxNumEntries() const
    {
      return(graph_.MaxNumIndices());
    }

    //! Returns a copy of the specified local row in user-provided arrays.
    /*! 
    \param In
           MyRow - Local row to extract.
    \param In
	   Length - Length of Values and Indices.
    \param Out
	   NumEntries - Number of nonzero entries extracted.
    \param Out
	   Values - Extracted values for this row.
    \param Out
	   Indices - Extracted global column indices for the corresponding values.
	  
    \return Integer error code, set to 0 if successful.
  */
    int ExtractMyRowCopy(int MyRow, int Length, int& NumEntries, double *Values, int * Indices) const
    {
      int err = graph_.ExtractMyRowCopy(MyRow,Length,NumEntries,Indices);
      if (!err)
      {
        int i,min;
        if (NumEntries<Length) min = NumEntries;
        else                   min = Length;
        for (i=0; i<min; i++) 
          Values[i] = 1.0;
      }
      return (err);
    }

    //! Returns a copy of the main diagonal in a user-provided vector.
    /*! 
    \param Out
	   Diagonal - Extracted main diagonal.

    \return Integer error code, set to 0 if successful.
  */
    int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.ExtractDiagonalCopy(...) is not impl!\n";
       throw -1;
       return(-1);
    }
  //@}
  
  //@{ \name Mathematical functions.

    //! Not implmented.
    /*! 
    \param In
	   TransA -If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
  */
    int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.Multiply(...) is not impl!\n";
       throw -1;
       return(-1);
    }

    //! Not implemented
    /*! This method will perform a triangular solve independently on each processor of the parallel machine.
        No communication is performed.
    \param In
	   Upper -If true, solve Ux = y, otherwise solve Lx = y.
    \param In
	   Trans -If true, solve transpose problem.
    \param In
	   UnitDiagonal -If true, assume diagonal is unit (whether it's stored or not).
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
    int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, 
	      Epetra_MultiVector& Y) const
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.Solve(...) is not impl!\n";
       throw -1;
       return(-1);
    }

    //! Not implemented
    /*! The vector x will return such that x[i] will contain the inverse of sum of the absolute values of the 
        \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.  Using the resulting vector from this function as input to LeftScale()
	will make the infinity norm of the resulting matrix exactly 1.
    \param Out
	   x -A Epetra_Vector containing the row sums of the \e this matrix. 
	   \warning It is assumed that the distribution of x is the same as the rows of \e this.

    \return Integer error code, set to 0 if successful.
  */
    int InvRowSums(Epetra_Vector& x) const
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.InvRowSums(...) is not impl!\n";
       throw -1;
       return(-1);
    }

    //! Not implemented
    /*! The \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the row number of A
        and j denotes the column number of A.
    \param In
	   x -A Epetra_Vector to solve for.

    \return Integer error code, set to 0 if successful.
  */
    int LeftScale(const Epetra_Vector& x) 
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.LeftScale(...) is not impl!\n";
       throw -1;
       return(-1);
    }

    //! Not implemented
    /*! The vector x will return such that x[j] will contain the inverse of sum of the absolute values of the 
        \e this matrix will be sca such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.  Using the resulting vector from this function as input to 
	RighttScale() will make the one norm of the resulting matrix exactly 1.
    \param Out
	   x -A Epetra_Vector containing the column sums of the \e this matrix. 
	   \warning It is assumed that the distribution of x is the same as the rows of \e this.

    \return Integer error code, set to 0 if successful.
  */
    int InvColSums(Epetra_Vector& x) const
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.InvColSums(...) is not impl!\n";
       throw -1;
       return(-1);
    }

    //! Not implemented
    /*! The \e this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.
    \param In
	   x -The Epetra_Vector used for scaling \e this.

    \return Integer error code, set to 0 if successful.
  */
    int RightScale(const Epetra_Vector& x) 
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.RightScale(...) is not impl!\n";
       throw -1;
       return(-1);
    }
  //@}
  
  //@{ \name Atribute access functions

    //! If FillComplete() has been called, this query returns true, otherwise it returns false.
    bool Filled() const 
    {
       return(graph_.Filled());
    }

    //! Not implemented
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
       \f[\| A \|_\infty = \max_{1\lei\len} \sum_{i=1}^m |a_{ij}| \f].
    */ 
    double NormInf() const
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.NormInf(...) is not impl!\n";
       throw -1;
       return(-1.0);
    }

    //! Not implemented
    /* Returns the quantity \f$ \| A \|_1\f$ such that
       \f[\| A \|_1= \max_{1\lej\len} \sum_{j=1}^n |a_{ij}| \f].
    */ 
    double NormOne() const
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.NormOne(...) is not impl!\n";
       throw -1;
       return(-1.0);
    }

    //! Returns the number of nonzero entries in the global matrix.
    int NumGlobalNonzeros() const {return(graph_.NumGlobalNonzeros()); }

    //! Returns the number of global matrix rows.
    int NumGlobalRows() const {return(graph_.NumGlobalRows()); }

    //! Returns the number of global matrix columns.
    int NumGlobalCols() const {return(graph_.NumGlobalCols()); }

    //! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons.
    int NumGlobalDiagonals() const {return(graph_.NumGlobalDiagonals()); }
    
    //! Returns the number of nonzero entries in the calling processor's portion of the matrix.
    int NumMyNonzeros() const {return(graph_.NumMyNonzeros()); }

    //! Returns the number of matrix rows owned by the calling processor.
    int NumMyRows() const {return(graph_.NumMyRows()); }

    //! Returns the number of matrix columns owned by the calling processor.
    int NumMyCols() const {return(graph_.NumMyCols()); }

    //! Returns the number of local nonzero diagonal entries, based on global row/column index comparisons.
    int NumMyDiagonals() const {return(graph_.NumMyDiagonals()); }

    //! If matrix is lower triangular in local index space, this query returns true, otherwise it returns false.
    bool LowerTriangular() const {return(graph_.LowerTriangular()); }

    //! If matrix is upper triangular in local index space, this query returns true, otherwise it returns false.
    bool UpperTriangular() const {return(graph_.UpperTriangular()); }

    //! Returns the Epetra_Map object associated with the rows of this matrix.
    const Epetra_Map & RowMatrixRowMap() const
    {
      return(dynamic_cast<const Epetra_Map&>(graph_.RowMap()));
    }

    //! Returns the Epetra_Map object associated with the columns of this matrix.
    const Epetra_Map & RowMatrixColMap() const 
    {
      return(dynamic_cast<const Epetra_Map &>(graph_.ColMap()));
    }

    //! Returns the Epetra_Import object that contains the import operations for distributed operations.
    const Epetra_Import * RowMatrixImporter() const 
    {
      return(graph_.Importer());
    }
    
    //! derived from Epetra_SrcDistObject
    const Epetra_BlockMap& Map() const
    {
      return(graph_.RowMap());
    }
    
    // derived from Epetra_Operator
    int SetUseTranspose(bool UseTranspose)
    {
      if (UseTranspose==true) usetranspose_=true;
      else                    usetranspose_=false;
      return(-1); //  the implementation does not support use of transpose
    }
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
    { 
       cout << "ML_Epetra::CrsGraphWrapper::Apply() is not implemented\n"; throw -1;
       return(-1);
    }
    
    //! Not implemented
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.ApplyInverse(...) is not impl!\n";
       throw -1;
       return(-1);
    }
    char * Label() const                         {return(name_);}
    bool UseTranspose() const                    {return(usetranspose_);}
    bool HasNormInf() const                      {return(false);}
    const Epetra_Comm & Comm() const             {return(comm_);}
    const Epetra_Map & OperatorDomainMap() const {return(DomainMap_);}
    const Epetra_Map & OperatorRangeMap() const  {return(RangeMap_);}
    
  //@}
 private:
    
    const Epetra_CrsGraph& graph_;
    bool                   usetranspose_;
    char*                  name_;
    Epetra_Comm&           comm_;
    Epetra_Map             DomainMap_;
    Epetra_Map             RangeMap_;
};

} // namespace ML_Epetra

#endif // defined(HAVE_ML_EPETRA)
#endif // ML_EPETRA_CRSGRAPHWRAP_H

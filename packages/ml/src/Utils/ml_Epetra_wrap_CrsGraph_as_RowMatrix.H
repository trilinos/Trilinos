/*!
 * \file ml_Epetra_wrap_CrsGraph_as_RowMatrix.H
 *
 * \class ML_Epetra::CrsGraphWrapper
 *
 * \brief Wraps a Epetra_CrsGraph as an Epetra_RowMatrix
 *
 * \date Last update to Doxygen: 01-Apr-05
 *
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

/*!
 
   \brief ML_Epetra::CrsGraphWrapper: a class to wrap an Epetra_CrsGraph as Epetra_RowMatrix

   Class ML_Epetra::CrsGraphWrapper takes an Epetra_CrsGraph object and wraps it as an 
   Epetra_RowMatrix. It  can then be used as input
   to ML to generate a plain aggegration MG hierarchy. Note that
   the resulting coarse grid operators are rubbish and need to be replaced before
   actually applying this hierarchy to something.
   This class' destructor does not destoy any of the passed objects, so they need to be 
   destoyed separately.
   
  
    To make use of the full functionality of this class it 
    requires ML to be configured with the following options:
    - \c --enable-epetra
    
    \author Michael Gee, SNL 9214
*/  
class CrsGraphWrapper: public virtual Epetra_RowMatrix 
{
      
 public:
  //@{ \name Constructor.
  //! Constructs a CrsGraphWrapper class.
  /*! 
   Constructs a CrsGraphWrapper class.
   As this wrapper implements an Epetra_RowMatrix, it can then be used as input
   to ML to generate a plain aggegration MG hierarchy. Note that
   the resulting coarse grid operators are rubbish and need to be replaced before
   actually applying this hierarchy to something.
   \param graph (In) : ref to Epetra_CrsGraph
   \param dm        (In) : ref to the DomainMap
   \param rm        (In) : ref to the RangeMap
   \param comm      (In) : the Epetra_Comm to be used
   */
    CrsGraphWrapper(const Epetra_CrsGraph& graph, const Epetra_Map& dm,const Epetra_Map& rm,const Epetra_Comm& comm) : 
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
    /*! As this class wraps a graph and implements an Epetra_RowMatrix,
        int-values are taken from the graph and converted to double values.
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

    //! not implemented, method will throw error and terminate execution.
    int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.ExtractDiagonalCopy(...) is not impl!\n";
       throw -1;
       return(-1);
    }
  //@}
  
  //@{ \name Mathematical functions.

    //! not implemented, method will throw error and terminate execution.
    int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.Multiply(...) is not impl!\n";
       throw -1;
       return(-1);
    }

    //! not implemented, method will throw error and terminate execution.
    int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, 
	      Epetra_MultiVector& Y) const
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.Solve(...) is not impl!\n";
       throw -1;
       return(-1);
    }

    //! not implemented, method will throw error and terminate execution.
    int InvRowSums(Epetra_Vector& x) const
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.InvRowSums(...) is not impl!\n";
       throw -1;
       return(-1);
    }

    //! not implemented, method will throw error and terminate execution.
    int LeftScale(const Epetra_Vector& x) 
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.LeftScale(...) is not impl!\n";
       throw -1;
       return(-1);
    }

    //! not implemented, method will throw error and terminate execution.
    int InvColSums(Epetra_Vector& x) const
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.InvColSums(...) is not impl!\n";
       throw -1;
       return(-1);
    }

    //! not implemented, method will throw error and terminate execution.
    int RightScale(const Epetra_Vector& x) 
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.RightScale(...) is not impl!\n";
       throw -1;
       return(-1);
    }
  //@}
  
  //@{ \name Attribute access functions

    //! If FillComplete() has been called, this query returns true, otherwise it returns false.
    bool Filled() const 
    {
       return(graph_.Filled());
    }

    //! not implemented, method will throw error and terminate execution.
    double NormInf() const
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.NormInf(...) is not impl!\n";
       throw -1;
       return(-1.0);
    }

    //! not implemented, method will throw error and terminate execution.
    double NormOne() const
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.NormOne(...) is not impl!\n";
       throw -1;
       return(-1.0);
    }

    //! Returns the number of nonzero entries in the global matrix.
    long long NumGlobalNonzeros() const {return(graph_.NumGlobalNonzeros()); }

    //! Returns the number of global matrix rows.
    long long NumGlobalRows() const {return(graph_.NumGlobalRows()); }

    //! Returns the number of global matrix columns.
    long long NumGlobalCols() const {return(graph_.NumGlobalCols()); }

    //! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons.
    long long NumGlobalDiagonals() const {return(graph_.NumGlobalDiagonals()); }
    
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
    
    //! Returns the Epetra_Map object associated with the rows of this matrix, derived from Epetra_Operator
    const Epetra_BlockMap& Map() const
    {
      return(graph_.RowMap());
    }
    
    //! not implemented, method will always return -1
    int SetUseTranspose(bool UseTranspose)
    {
      usetranspose_ = UseTranspose;
      return(-1); //  the implementation does not support use of transpose
    }
    
    //! not implemented, method will throw error and terminate execution.
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
    { 
       cout << "ML_Epetra::CrsGraphWrapper::Apply() is not implemented\n"; throw -1;
       return(-1);
    }
    
    //! not implemented, method will throw error and terminate execution.
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
       cout << "**ERR** ML_Epetra::CrsGraphWrapper.ApplyInverse(...) is not impl!\n";
       throw -1;
       return(-1);
    }

    //! returns the label of this class
    const char * Label() const                   {return(name_);}

    //! returns the usetranspose_ flag
    bool UseTranspose() const                    {return(usetranspose_);}

    //! always returns false
    bool HasNormInf() const                      {return(false);}

    //! returns ref to the Epetra_Comm associated with this class
    const Epetra_Comm & Comm() const             {return(comm_);}

    //! returns ref to the OperatorDomainMap associated with this class
    const Epetra_Map & OperatorDomainMap() const {return(DomainMap_);}

    //! returns ref to the OperatorRangeMap associated with this class
    const Epetra_Map & OperatorRangeMap() const  {return(RangeMap_);}
    
  //@}
 private:
    
    const Epetra_CrsGraph& graph_;
    bool                   usetranspose_;
    const char*            name_;
    Epetra_Map             DomainMap_;
    Epetra_Map             RangeMap_;
    const Epetra_Comm&     comm_;
};

} // namespace ML_Epetra

#endif // defined(HAVE_ML_EPETRA)
#endif // ML_EPETRA_CRSGRAPHWRAP_H

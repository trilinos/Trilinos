
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

#ifndef _EPETRA_CRSSINGLETONFILTER_H_
#define _EPETRA_CRSSINGLETONFILTER_H_

#include "Epetra_Object.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MapColoring.h"
class Epetra_LinearProblem;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Import;
class Epetra_Export;
class Epetra_IntVector;

//! Epetra_CrsSingletonFilter: A class for explicitly eliminating matrix rows and columns.

/*! The Epetra_CrsSingletonFilter class takes an existing Epetra_LinearProblem object, analyzes
    it structure and explicitly eliminates rows and columns from the matrix based on density
    of nonzero entries.
*/    

class Epetra_CrsSingletonFilter {
      
 public:

  //@{ \name Constructors/Destructor.
  //! Epetra_CrsSingletonFilter default constructor.
  Epetra_CrsSingletonFilter();

  //! Epetra_CrsSingletonFilter Destructor
  virtual ~Epetra_CrsSingletonFilter();
  //@}
  //@{ \name Analyze methods.
  //! Analyze the input matrix, removing row/column pairs that have singletons.
  /*! Analyzes the user's input matrix to determine rows and columns that should be explicitly
      eliminated to create the reduced system.  Look for rows and columns that have single entries.  
      These rows/columns
      can easily be removed from the problem.  
      The results of calling this method are two MapColoring objects accessible via RowMapColors() and 
      ColMapColors() accessor methods.  All rows/columns that would be eliminated in the reduced system
      have a color of 1 in the corresponding RowMapColors/ColMapColors object.  All kept rows/cols have a 
      color of 0.
  */
  int Analyze(Epetra_RowMatrix * FullMatrix);

  //! Returns true if singletons were detected in this matrix (must be called after Analyze() to be effective).
  bool SingletonsDetected() const {if (!AnalysisDone_) return(false); else return(RowMapColors_->MaxNumColors()>1);};

  //! Analyze the input matrix using density thresholds and detection of singletons to determine rows/cols to eliminate.
  /*! Analyzes the user's input matrix to determine rows and columns that should be explicitly
      eliminated to create the reduced system.  First, look for rows and columns that have single entries.  
      These rows/columns
      can easily be removed from the problem.  Second, if one or both of AbsoluteThreshold and RelativeThreshold are
      nonzero, rows/columns above the threshold will be eliminated.
      The results of calling this method are two Epetra_Maps, RowEliminateMap and ColEliminateMap,
      that contain the list of global IDs of the row and columns to be eliminated.  These two maps can be access via
      RowEliminateMap() and ColEliminateMap() accessor methods.
    
    \param In 
           AbsoluteThreshold - If any row or column has a count above this value, the row/column
	   pair will be formally eliminated from the reduced matrix.  Setting this value to 0 will
	   turn off this test.
    \param In 
           RelativeThreshold - If the ratio of the nonzero count of any row or column over the average
	   row or column count is above the relative threshold, the row/column
	   pair will be formally eliminated from the reduced matrix.  Setting this value to 0 will
	   turn off this test.
  */
  int Analyze(int AbsoluteThreshold, double RelativeThreshold);

  //! Print statistics about the reduction analysis.
  int Statistics() const;
  //@}

  //@{ \name Reduce methods.
  //! Return a reduced linear problem based on results of Analyze().
  /*! Creates a new Epetra_LinearProblem object based on the results of the Analyze phase.  A pointer
      to the reduced problem is obtained via a call to ReducedProblem().  
    	   
    \return Error code, set to 0 if no error.
  */
  int ConstructReducedProblem(Epetra_LinearProblem * Problem);

  //! Update a reduced linear problem using new values.
  /*! Updates an existing Epetra_LinearProblem object using new matrix, LHS and RHS values.  The matrix
      structure must be \e identical to the matrix that was used to construct the original reduced problem.  
    	   
    \return Error code, set to 0 if no error.
  */
  int UpdateReducedProblem(Epetra_LinearProblem * Problem);

  //! Return a reduced linear problem using user-defined elimination maps.
  /*! Creates a new Epetra_LinearProblem object based on the maps that the user
      provides.  A pointer
      to the reduced problem is obtained via a call to ReducedProblem().  
    
    \param In 
           RowEliminateMap - An Epetra_Map specifying the global IDs of rows and columns that
	   should be explicitly eliminated from the matrix.
    \param In 
           RowEliminateMap - An Epetra_Map specifying the global IDs of rows and columns that
	   should be explicitly eliminated from the matrix.
	   
    \return Error code, set to 0 if no error.
  */
  int ConstructReducedProblem(const Epetra_Map & RowEliminateMap, const Epetra_Map & ColEliminateMap);

  //@}
  //@{ \name Methods to construct Full System Solution.
  //! Compute a solution for the full problem using the solution of the reduced problem, put in LHS of FullProblem().
  /*! After solving the reduced linear system, this method can be called to compute the
      solution to the original problem, assuming the solution for the reduced system is valid. The solution of the 
      unreduced, original problem will be in the LHS of the original Epetra_LinearProblem.
    
  */
  int ComputeFullSolution();
  //@}
  //@{ \name Attribute Access Methods.

  //! Returns pointer to the original unreduced Epetra_LinearProblem.
  Epetra_LinearProblem * FullProblem() const {return(FullProblem_);};

  //! Returns pointer to the derived reduced Epetra_LinearProblem.
  Epetra_LinearProblem * ReducedProblem() const {return(ReducedProblem_);};

  //! Returns pointer to Epetra_CrsMatrix from full problem.
  Epetra_RowMatrix * FullMatrix() const {return(FullMatrix_);};

  //! Returns pointer to Epetra_CrsMatrix from full problem.
  Epetra_CrsMatrix * FullCrsMatrix() const {return(FullCrsMatrix_);};

  //! Returns pointer to Epetra_CrsMatrix from full problem.
  Epetra_CrsMatrix * ReducedMatrix() const {return(ReducedMatrix_);};

  //! Returns Epetra_Map containing global row IDs of full matrix.
  const Epetra_Map & FullMatrixRowMap() const {return(FullMatrix()->RowMatrixRowMap());};

  //! Returns Epetra_Map containing column map of full matrix.
  const Epetra_Map & FullMatrixColMap() const {return(FullMatrix()->RowMatrixColMap());};

  //! Returns Epetra_Map containing domain map of full matrix.
  const Epetra_Map & FullMatrixDomainMap() const {return((FullMatrix()->OperatorDomainMap()));};

  //! Returns Epetra_Map containing range map of full matrix.
  const Epetra_Map & FullMatrixRangeMap() const {return((FullMatrix()->OperatorRangeMap()));};

  //! Returns pointer to Epetra_MapColoring object: color 0 rows are part of reduced system.
  Epetra_MapColoring * RowMapColors() const {return(RowMapColors_);};

  //! Returns pointer to Epetra_MapColoring object: color 0 columns are part of reduced system.
  Epetra_MapColoring * ColMapColors() const {return(ColMapColors_);};

  //! Returns pointer to Epetra_Map describing the reduced system row distribution.
  Epetra_Map * ReducedMatrixRowMap() const {return(ReducedMatrixRowMap_);};

  //! Returns pointer to Epetra_Map describing the reduced system column distribution.
  Epetra_Map * ReducedMatrixColMap() const {return(ReducedMatrixColMap_);};

  //! Returns pointer to Epetra_Map describing the domain map for the reduced system.
  Epetra_Map * ReducedMatrixDomainMap() const {return(ReducedMatrixDomainMap_);};

  //! Returns pointer to Epetra_Map describing the range map for the reduced system.
  Epetra_Map * ReducedMatrixRangeMap() const {return(ReducedMatrixRangeMap_);};
  //@}

 protected:

    void InitializeDefaults();
    int ComputeEliminateMaps();
    int Setup(Epetra_LinearProblem * Problem);
    int InitFullMatrixAccess();
    int GetRow(int Row, int & NumIndices, int * & Indices);
    int GetRowGCIDs(int Row, int & NumIndices, double * & Values, int * & GlobalIndices);
    int GetRow(int Row, int & NumIndices, double * & Values, int * & Indices);
    int CreatePostSolveArrays(int * RowIDs,
			      const Epetra_MapColoring & RowMapColors,
			      const Epetra_IntVector & ColProfiles,
			      const Epetra_IntVector & NewColProfiles,
			      const Epetra_IntVector & ColHasRowWithSingleton);

    int ConstructRedistributeExporter(Epetra_Map * SourceMap, Epetra_Map * TargetMap,
				      Epetra_Export * & RedistributeExporter,
				      Epetra_Map * & RedistributeMap);

    Epetra_LinearProblem * FullProblem_;
    Epetra_LinearProblem * ReducedProblem_;
    Epetra_RowMatrix * FullMatrix_;
    Epetra_CrsMatrix * FullCrsMatrix_;
    Epetra_CrsMatrix * ReducedMatrix_;
    Epetra_MultiVector * ReducedRHS_;
    Epetra_MultiVector * ReducedLHS_;

    Epetra_Map * ReducedMatrixRowMap_;
    Epetra_Map * ReducedMatrixColMap_;
    Epetra_Map * ReducedMatrixDomainMap_;
    Epetra_Map * ReducedMatrixRangeMap_;
    Epetra_Map * OrigReducedMatrixDomainMap_;
    Epetra_Import * Full2ReducedRHSImporter_;
    Epetra_Import * Full2ReducedLHSImporter_;
    Epetra_Export * RedistributeDomainExporter_;
    
    int * ColSingletonRowLIDs_;
    int * ColSingletonColLIDs_;
    int * ColSingletonPivotLIDs_;
    double * ColSingletonPivots_;

    
    int AbsoluteThreshold_;
    double RelativeThreshold_;

    int NumRowSingletons_;
    int NumColSingletons_;

    bool HaveReducedProblem_;
    bool UserDefinedEliminateMaps_;
    bool AnalysisDone_;
    bool SymmetricElimination_;

    Epetra_MultiVector * tempExportX_;
    Epetra_MultiVector * tempX_;
    Epetra_MultiVector * tempB_;
    Epetra_MultiVector * RedistributeReducedLHS_;
    int * Indices_;
    double * Values_;

    Epetra_MapColoring * RowMapColors_;
    Epetra_MapColoring * ColMapColors_;
    bool FullMatrixIsCrsMatrix_;
    int MaxNumMyEntries_;
	
    
 private:
 //! Copy constructor (defined as private so it is unavailable to user).
  Epetra_CrsSingletonFilter(const Epetra_CrsSingletonFilter & Problem){};
};
#endif /* _EPETRA_CRSSINGLETONFILTER_H_ */

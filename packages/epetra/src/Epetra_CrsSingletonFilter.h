
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
class Epetra_LinearProblem;
class Epetra_Map;
class Epetra_MultiVector;

//! Epetra_CrsSingletonFilter: A class for explicitly eliminating matrix rows and columns.

/*! The Epetra_CrsSingletonFilter class takes an existing Epetra_LinearProblem object, analyzes
    it structure and explicitly eliminates rows and columns from the matrix based on density
    of nonzero entries.
*/    

class Epetra_CrsSingletonFilter {
      
 public:

  //@{ \name Constructors/Destructor.
  //! Epetra_CrsSingletonFilter default constructor.
  Epetra_CrsSingletonFilter(Epetra_LinearProblem * Problem);

  //! Epetra_CrsSingletonFilter Destructor
  virtual ~Epetra_CrsSingletonFilter();
  //@}
  //@{ \name Analyze methods.
  //! Analyze the input matrix, removing row/column pairs that have singletons.
  /*! Analyzes the user's input matrix to determine rows and columns that should be explicitly
      eliminated to create the reduced system.  Look for rows and columns that have single entries.  
      These rows/columns
      can easily be removed from the problem.  
      The results of calling this method are four Epetra_Maps, RowEliminateMap and ColEliminateMap, 
      ReducedMatrixRowMap and ReducedMatrixDomainMap
      that contain the list of global IDs of the row and columns to be eliminated and kept, respectively.  
      These  maps can be access via
      RowEliminateMap(), ColEliminateMap(), ReducedMatrixRowMap() and ReducedMatrixDomainMap() accessor methods.
  */
  int Analyze();

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
  int Statistics();
  //@}

  //@{ \name Reduce methods.
  //! Return a reduced linear problem based on results of Analyze().
  /*! Creates a new Epetra_LinearProblem object based on the results of the Analyze phase.  A pointer
      to the reduced problem is obtained via a call to ReducedProblem().  
    	   
    \return Error code, set to 0 if no error.
  */
  int ConstructReducedProblem();

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
  Epetra_LinearProblem * FullProblem(){return(FullProblem_);};
  //! Returns pointer to the derived reduced Epetra_LinearProblem.
  Epetra_LinearProblem * ReducedProblem(){return(ReducedProblem_);};
  //! Returns pointer to Epetra_CrsMatrix from full problem.
  Epetra_CrsMatrix * FullMatrix(){return(FullMatrix_);};
  //! Returns pointer to Epetra_CrsMatrix from full problem.
  Epetra_CrsMatrix * ReducedMatrix(){return(ReducedMatrix_);};
  //! Returns Epetra_Map containing Global Row IDs of full matrix.
  const Epetra_Map & FullMatrixRowMap(){return(FullMatrix()->RowMap());};
  //! Returns Epetra_Map containing domain map of full matrix.
  const Epetra_Map & FullMatrixDomainMap(){return(dynamic_cast<const Epetra_Map &>(FullMatrix()->DomainMap()));};
  //! Returns Epetra_Map containing import map of full matrix.
  const Epetra_Map & FullMatrixImportMap(){return(FullMatrix()->ImportMap());};
  //! Returns pointer to Epetra_Map containing Global Row IDs that are eliminated from full problem.
  Epetra_Map * RowEliminateMap(){return(RowEliminateMap_);};
  //! Returns pointer to Epetra_Map containing Global Column IDs that are eliminated from full problem.
  Epetra_Map * ColEliminateMap(){return(ColEliminateMap_);};
  //! Returns pointer to Epetra_Map describing the reduced system row distribution.
  Epetra_Map * ReducedMatrixRowMap(){return(ReducedMatrixRowMap_);};
  //! Returns pointer to Epetra_Map describing the domain map for the reduced system.
  Epetra_Map * ReducedMatrixDomainMap(){return(ReducedMatrixDomainMap_);};
  //@}

 protected:

    void InitializeDefaults();
    int ComputeEliminateMaps();
    int Setup(Epetra_LinearProblem * Problem);

    Epetra_LinearProblem * FullProblem_;
    Epetra_LinearProblem * ReducedProblem_;
    Epetra_CrsMatrix * FullMatrix_;
    Epetra_CrsMatrix * ReducedMatrix_;
    Epetra_MultiVector * ReducedRHS_;
    Epetra_MultiVector * ReducedLHS_;

    Epetra_Map * RowEliminateMap_;
    Epetra_Map * ColEliminateMap_;
    Epetra_Map * ReducedMatrixRowMap_;
    Epetra_Map * ReducedMatrixDomainMap_;
    Epetra_Import * Full2ReducedRHSImporter_;
    Epetra_Import * Full2ReducedLHSImporter_;
    
    int * ColSingletonRowLIDs_;
    int * ColSingletonColLIDs_;
    double * ColSingletonPivots_;

    
    int AbsoluteThreshold_;
    double RelativeThreshold_;

    int NumRowSingletons_;
    int NumColSingletons_;

    bool HaveReducedProblem_;
    bool UserDefinedEliminateMaps_;
    bool AnalysisDone_;
    
 private:
 //! Copy constructor (defined as private so it is unavailable to user).
  Epetra_CrsSingletonFilter(const Epetra_CrsSingletonFilter & Problem){};
};
#endif /* _EPETRA_CRSSINGLETONFILTER_H_ */

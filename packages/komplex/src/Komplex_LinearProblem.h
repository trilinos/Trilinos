
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

#ifndef KOMPLEX_LINEARPROBLEM_H
#define KOMPLEX_LINEARPROBLEM_H

#include "Epetra_Object.h"
#include "Epetra_CrsMatrix.h"
class Epetra_LinearProblem;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Import;
class Epetra_Export;
class Epetra_MapColoring;
class Epetra_IntVector;

//! Komplex_LinearProblem: A class for explicitly eliminating matrix rows and columns.

/*! The Komplex_LinearProblem class takes an existing Epetra_LinearProblem object, analyzes
    it structure and explicitly eliminates rows and columns from the matrix based on density
    of nonzero entries.
*/    

class Komplex_LinearProblem {
      
 public:

  //@{ \name Constructors/Destructor.
  //! Komplex_LinearProblem default constructor.
  Komplex_LinearProblem();

  //! Komplex_LinearProblem Destructor
  virtual ~Komplex_LinearProblem();
  //@}
  //@{ \name Analyze methods.
  //! Analyze the input matrices tat are in Vbr format to determine best Komplex formulation (Unimplemented at this time, defaults to K1 form).
  /*! This method will analyze the real and imaginary parts of the complex-valued problem to determine the best Komplex
      formulation.  It is designed to work specifically with Vbr Matrix objects.
      It is unimplemented at this time.  The K1 formulation is used by default.
  */
  int AnalyzeVbrMatrices(const Epetra_VbrMatrix & Ar, const Epetra_VbrMatrix & Ai);
  //! Analyze the input matrices that implement the Epetra_RowMatrix interface to determine best Komplex formulation (Unimplemented at this time, defaults to K1 form).
  /*! This method will analyze the real and imaginary parts of the complex-valued problem to determine the best Komplex
      formulation.  It will work with any matrix that implements the Epetra_RowMatrix interface (including Epetra_VbrMatrix
      objects).  If the user matrix is an Epetra_VbrMatrix object, the AnalyzeVbrMatrices method is usually a better approach.
      It is unimplemented at this time.  The K1 formulation is used by default.
  */
  int AnalyzeRowMatrices(const Epetra_RowMatrix & Ar, const Epetra_RowMatrix & Ai);

  //! Print statistics about the reduction analysis (Unimplemented at this time).
  int Statistics() const{return(0);};
  //@}

  //@{ \name Construction/update methods.
  //! Construct a komplex linear problem from VbrMatrix objects based on results of Analyze().
  /*! Creates a new Epetra_LinearProblem object containing the komplex (equivalent real formulation)
       based on the results of the Analyze phase.  A pointer
      to the komplex problem is obtained via a call to KomplexProblem().  
    	   
    \return Error code, set to 0 if no error.
  */
  int ConstructKomplexProblemFromVbrMatrices(const Epetra_VbrMatrix & Ar, const Epetra_VbrMatrix & Ai, 
				     const Epetra_MultiVector & Xr, const Epetra_MultiVector & Xi,
				     const Epetra_MultiVector & Br, const Epetra_MultiVector & Bi,
				     );
  //! Construct a komplex linear problem from RowMatrix objects based on results of Analyze().
  /*! Creates a new Epetra_LinearProblem object containing the komplex (equivalent real formulation)
       based on the results of the Analyze phase.  A pointer
      to the komplex problem is obtained via a call to KomplexProblem().  
    	   
    \return Error code, set to 0 if no error.
  */
  int ConstructKomplexProblemFromRowMatrices(const Epetra_RowMatrix & Ar, const Epetra_RowMatrix & Ai, 
				     const Epetra_MultiVector & Xr, const Epetra_MultiVector & Xi,
				     const Epetra_MultiVector & Br, const Epetra_MultiVector & Bi,
				     );

  //! Update an existing Komplex linear problem using new values.
  /*! Updates an existing Komplex_LinearProblem object using new matrix, LHS and RHS values.  The matrix
      structure must be \e identical to the matrix that was used to construct the original reduced problem.  
    	   
    \return Error code, set to 0 if no error.
  */
  int UpdateKomplexProblemFromVbrMatrices(const Epetra_VbrMatrix & Ar, const Epetra_VbrMatrix & Ai, 
					  const Epetra_MultiVector & Xr, const Epetra_MultiVector & Xi,
					  const Epetra_MultiVector & Br, const Epetra_MultiVector & Bi,
					  );

  //! Update an existing Komplex linear problem using new values.
  /*! Updates an existing Komplex_LinearProblem object using new matrix, LHS and RHS values.  The matrix
      structure must be \e identical to the matrix that was used to construct the original reduced problem.  
    	   
    \return Error code, set to 0 if no error.
  */
  int UpdateKomplexProblemFromRowMatrices(const Epetra_VbrMatrix & Ar, const Epetra_VbrMatrix & Ai, 
					  const Epetra_MultiVector & Xr, const Epetra_MultiVector & Xi,
					  const Epetra_MultiVector & Br, const Epetra_MultiVector & Bi,
					  );
  //@}
  //@{ \name Methods to extract complex system solution.
  //! Extrac a solution for the original complex-valued problem using the solution of the Komplex problem.
  /*! After solving the komplex linear system, this method can be called to extract the
      solution of the original problem, assuming the solution for the komplex system is valid.
    
  */
  int ExtractSolution(Epetra_MultiVector * Xr, Epetra_MultiVector * Xi);
  //@}
  //@{ \name Attribute Access Methods.

  //! Returns pointer to the Epetra_LinearProblem object that defines the Komplex formulation.
  Epetra_LinearProblem * KomplexProblem() const {return(KomplexProblem_);};

  //! Returns pointer to Epetra_VbrMatrix from Komplex problem.
  Epetra_VbrMatrix * KomplexVbrMatrix() const {return(KomplexVbrMatrix_);};
  //@}

 protected:

    void InitializeDefaults();
    int ComputeKomplexMaps();
    int Setup(Epetra_LinearProblem * Problem);
    int InitFullMatrixAccess();
    int GetRow(int Row, int & NumIndices, int * & Indices);
    int GetRowGCIDs(int Row, int & NumIndices, double * & Values, int * & GlobalIndices);
    int GetRow(int Row, int & NumIndices, double * & Values, int * & Indices);

    Epetra_LinearProblem * KomplexProblem_;
    Epetra_VbrMatrix * KomplexMatrix_;
    Epetra_MultiVector * KomplexRHS_;
    Epetra_MultiVector * KomplexLHS_;

    Epetra_BlockMap * KomplexMatrixRowMap_;
    Epetra_BlockMap * KomplexMatrixColMap_;
    Epetra_BlockMap * KomplexMatrixDomainMap_;
    Epetra_BlockMap * KomplexMatrixRangeMap_;
    

    bool HaveKomplexProblem_;
    bool AnalysisDone_;

    int * Indices_;
    double * Values_;

    bool UserMatrixIsCrsMatrix_;
    bool UserMatrixIsVbrMatrix_;
    int MaxNumMyEntries_;
	
    
 private:
 //! Copy constructor (defined as private so it is unavailable to user).
  Komplex_LinearProblem(const Komplex_LinearProblem & Problem){};
};
#endif /* KOMPLEX_LINEARPROBLEM_H */


//@HEADER
// ***********************************************************************
// 
//                Komplex: Complex Linear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#ifndef KOMPLEX_LINEARPROBLEM_H
#define KOMPLEX_LINEARPROBLEM_H

class Epetra_LinearProblem;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Import;
class Epetra_Export;
class Epetra_MapColoring;
class Epetra_IntVector;
class Epetra_VbrMatrix;

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
  //@{ \name Set methods.
  //! Construct the komplex linear operator from VbrMatrix objects.
  /*! Constructs the Komplex operator from the user definition 
      of the complex-valued matrix C = (c0r+i*c0i)*A0 +(c1r+i*c1i)*A1.
      Using this general expression for the complex matrix allows easy formulation of a variety of common
      complex problems.

      The operator will be explicitly constructed as an Epetra_VbrMatrix object when the first call to
      SetKomplexOperator() is made.  Subsequent calls to this method will attempt to reuse the the existing
      KomplexVbrMatrix object if possible, rather than reconstructing from scratch.  If this is not possible (typically
      because the structure has changed) then a the previous KomplexVbrMatrix object will be deleted and a new one will be 
      constructed.

      \param c0r (In) The real part of the complex coefficient multiplying A0.
      \param c0i (In) The imag part of the complex coefficient multiplying A0.
      \param A0 (In) An Epetra_RowMatrix that is one of the matrices used to define the true complex operator.
      \param c1r (In) The real part of the complex coefficient multiplying A1.
      \param c1i (In) The imag part of the complex coefficient multiplying A1.
      \param A1 (In) An Epetra_RowMatrix that is the second of the matrices used to define the true complex operator.

    \return Error code, set to 0 if no error.
  */
  int SetKomplexOperator(double c0r, double c0i, const Epetra_RowMatrix & A0,
			 double c1r, double c1i, const Epetra_RowMatrix & Ai);

  //! Set the left hand side of a Komplex linear problem.
  /*! 
    \param Xr (In) The real part of the complex valued LHS.  
    \param Xi (In) The imag part of the complex valued LHS.  
    	   
    \return Error code, set to 0 if no error.
  */
  int SetKomplexLHS(const Epetra_MultiVector & Xr, const Epetra_MultiVector & Xi);

  //! Set the right hand side of a Komplex linear problem.
  /*! 
    \param Br (In) The real part of the complex valued RHS.  
    \param Bi (In) The imag part of the complex valued RHS.  
    	   
    \return Error code, set to 0 if no error.
  */
  int SetKomplexRHS(const Epetra_MultiVector & Br, const Epetra_MultiVector & Bi);
  //@}
  //@{ \name Methods to extract complex system solution.
  //! Extrac a solution for the original complex-valued problem using the solution of the Komplex problem.
  /*! After solving the komplex linear system, this method can be called to extract the
      solution of the original problem, assuming the solution for the komplex system is valid.
    \param Xr (In) The real part of the complex valued solution.  
    \param Xi (In) The imag part of the complex valued solution. 
    
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
    SetKomplexOperatorVbr(double c0r, double c0i, const Epetra_VbrMatrix & A0,
			  double c1r, double c1i, const Epetra_VbrMatrix & Ai);
    SetKomplexOperatorRow(double c0r, double c0i, const Epetra_RowMatrix & A0,
			  double c1r, double c1i, const Epetra_RowMatrix & Ai);
    int ConstructKomplexMaps(const Epetra_BlockMap & A0DomainMap, const Epetra_BlockMap & A0RangeMap, 
			     const Epetra_BlockMap & A0RowMap, const Epetra_BlockMap & A0ColMap,
			     const Epetra_BlockMap & A1DomainMap, const Epetra_BlockMap & A1RangeMap, 
			     const Epetra_BlockMap & A1RowMap, const Epetra_BlockMap & A1ColMap);
    int MakeKomplexMap(const Epetra_BlockMap & Map, Epetra_BlockMap * & KMap); 
    int InitFullMatrixAccess();
    int GetRow(int Row, int & NumIndices, int * & Indices);
    int GetRowGCIDs(int Row, int & NumIndices, double * & Values, int * & GlobalIndices);
    int GetRow(int Row, int & NumIndices, double * & Values, int * & Indices);

    Epetra_LinearProblem * KomplexProblem_;
    Epetra_VbrMatrix * KomplexMatrix_;
    Epetra_CrsGraph * KomplexGraph_;
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

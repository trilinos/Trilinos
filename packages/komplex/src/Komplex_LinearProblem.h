
//@HEADER
// ***********************************************************************
// 
//                Komplex: Complex Linear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#ifndef KOMPLEX_LINEARPROBLEM_H
#define KOMPLEX_LINEARPROBLEM_H

#include "Teuchos_RCP.hpp"
#include <vector>
class Epetra_LinearProblem;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Import;
class Epetra_Export;
class Epetra_MapColoring;
class Epetra_IntVector;
class Epetra_VbrMatrix;

//! Komplex_LinearProblem: A class for forming an equivalent real formulation of a complex valued problem.

/*! The Komplex_LinearProblem class takes a complex linear problem, separated into real and imaginary parts,
    and forms an equivalent real valued system of twice the dimension.  The resulting system can then be
    solved with any Trilinos solver that understands Epetra objects.

KOMPLEX solves a complex-valued linear system Ax = b by solving
an equivalent real-valued system of twice the dimension.  Specifically,
writing in terms of real and imaginary parts, we have

 \f[ (A_r + i*A_i)*(x_r + i*x_i) = (b_r + i*b_i) \f]

  or by separating into real and imaginary equations we have

\f[
  \left( \begin{array}{rr}
                                    A_r & -A_i\\
                                    A_i &  A_r
                             \end{array}
   \right)
   \left( \begin{array}{r}
                                    x_r\\
                                    x_i
                             \end{array}
   \right)
   =
   \left( \begin{array}{r}
                                    b_r\\
                                    b_i
                             \end{array}
   \right)
\f]
  which is a real-valued system of twice the size.  If we find xr and xi, we
  can form the solution to the original system as x = xr +i*xi.


KOMPLEX accepts the user linear system as two real-valued matrices with no assumption
about the structure of the matrices, except that they have compatible RowMap, DomainMap 
and RangeMap distributions.  Each matrix is multiplied by
user-supplied complex constants.
 
Although formally the system is a 2-by-2 block system, we actually apply the interleaving at the matrix entry level
such that the real part of the first complex equation is followed by the imaginary part of the first complex equation,
and so on.  This approach is documented in:
 
 David Day and Michael A. Heroux. Solving complex-valued linear systems via equivalent real formulations. 
 SIAM J. Sci. Comput., 23(2):480–498, 2001.


*/    

class Komplex_LinearProblem {
      
 public:

  //@{ \name Constructors/Destructor.
  //! Komplex_LinearProblem constructor.
  /*! Constructs the Komplex operator from the user definition 
      of the complex-valued matrix C = (c0r+i*c0i)*A0 +(c1r+i*c1i)*A1.
      Using this general expression for the complex matrix allows easy formulation of a variety of common
      complex problems.

      \param c0r (In) The real part of the complex coefficient multiplying A0.
      \param c0i (In) The imag part of the complex coefficient multiplying A0.
      \param A0 (In) An Epetra_RowMatrix that is one of the matrices used to define the true complex operator.
      \param c1r (In) The real part of the complex coefficient multiplying A1.
      \param c1i (In) The imag part of the complex coefficient multiplying A1.
      \param A1 (In) An Epetra_RowMatrix that is the second of the matrices used to define the true complex operator.
      \param Xr (In) The real part of the complex valued LHS.  
      \param Xi (In) The imag part of the complex valued LHS.  
      \param Br (In) The real part of the complex valued RHS.  
      \param Bi (In) The imag part of the complex valued RHS.  
      
    */
  Komplex_LinearProblem(double c0r, double c0i, const Epetra_RowMatrix & A0,
			double c1r, double c1i, const Epetra_RowMatrix & A1,
			const Epetra_MultiVector & Xr, const Epetra_MultiVector & Xi,
			const Epetra_MultiVector & Br, const Epetra_MultiVector & Bi);

  //! Komplex_LinearProblem Destructor
  virtual ~Komplex_LinearProblem();
  //@}
  //@{ \name Set methods.
  //! Update the values of the equivalent real valued system.
  /*! This method allows the values of an existing Komplex_LinearProblem object to be updated.  Note that
      the update that there is no change to the pattern of the matrices. 
    	   
    \return Error code, set to 0 if no error.
  */
  int UpdateValues(double c0r, double c0i, const Epetra_RowMatrix & A0,
		   double c1r, double c1i, const Epetra_RowMatrix & A1,
		   const Epetra_MultiVector & Xr, const Epetra_MultiVector & Xi,
		   const Epetra_MultiVector & Br, const Epetra_MultiVector & Bi);
  //@}
  //@{ \name Methods to extract complex system solution.
  //! Extrac a solution for the original complex-valued problem using the solution of the Komplex problem.
  /*! After solving the komplex linear system, this method can be called to extract the
      solution of the original problem, assuming the solution for the komplex system is valid.
    \param Xr (Out) An existing Epetra_MultiVector.  On exit it will contain the real part of the complex valued solution.  
    \param Xi (Out) An existing Epetra_MultiVector.  On exit it will contain the imag part of the complex valued solution. 
    
  */
  int ExtractSolution(Epetra_MultiVector & Xr, Epetra_MultiVector & Xi);
  //@}
  //@{ \name Attribute Access Methods.

  //! Returns pointer to the Epetra_LinearProblem object that defines the Komplex formulation.
  /*! The pointer returned from this method will contain the address of a fully-constructed Epetra_LinearProblem
      instance that can be used with any Trilinos preconditioner or solver.
  */
  Epetra_LinearProblem * KomplexProblem() const {return(KomplexProblem_.get());}

  //@}

 protected:
  int ProcessValues(double c0r, double c0i, const Epetra_RowMatrix & A0,
		    double c1r, double c1i, const Epetra_RowMatrix & A1,
		    const Epetra_MultiVector & Xr, const Epetra_MultiVector & Xi,
		    const Epetra_MultiVector & Br, const Epetra_MultiVector & Bi,
		    bool firstTime);
  
  int TestMaps (const Epetra_RowMatrix & A0, const Epetra_RowMatrix & A1,
		const Epetra_MultiVector & Xr, const Epetra_MultiVector & Xi,
		const Epetra_MultiVector & Br, const Epetra_MultiVector & Bi);
  
  int ConstructKomplexMaps(const Epetra_Map & A0DomainMap, const Epetra_Map & A0RangeMap,
			   const Epetra_Map & A0RowMap);
  
  int MakeKomplexMap(const Epetra_Map & Map, Teuchos::RCP<Epetra_Map> & KMap);
  
  int InitMatrixAccess(const Epetra_RowMatrix & A0, const Epetra_RowMatrix & A1);
 
  int GetRow(int Row, const Epetra_RowMatrix & A0, const Epetra_RowMatrix & A1,
	     int & NumIndices0, double * & Values0, int * & Indices0,
	     int & NumIndices1, double * & Values1, int * & Indices1);

  int PutRow(int Row, int & NumIndices, double * Values, int * Indices, bool firstTime); 

  Teuchos::RCP<Epetra_LinearProblem> KomplexProblem_;
  Teuchos::RCP<Epetra_CrsMatrix> KomplexMatrix_;
  Teuchos::RCP<Epetra_MultiVector> KomplexRHS_;
  Teuchos::RCP<Epetra_MultiVector> KomplexLHS_;
  
  Teuchos::RCP<Epetra_Map> KomplexMatrixRowMap_;
  Teuchos::RCP<Epetra_Map> KomplexMatrixColMap_;
  Teuchos::RCP<Epetra_Map> KomplexMatrixDomainMap_;
  Teuchos::RCP<Epetra_Map> KomplexMatrixRangeMap_;
  
  const Epetra_CrsMatrix * CrsA0_;
  const Epetra_CrsMatrix * CrsA1_;
  bool A0A1AreCrs_;
  
  std::vector<int> Indices0_;
  std::vector<double> Values0_;
  int MaxNumMyEntries0_;
  std::vector<int> Indices1_;
  std::vector<double> Values1_;
  int MaxNumMyEntries1_;	
  std::vector<int> IndicesK_;
  std::vector<double> ValuesK_;
  int MaxNumMyEntriesK_;	
  
 private:
  //! Copy constructor (defined as private so it is unavailable to user).
  Komplex_LinearProblem(const Komplex_LinearProblem & Problem){};
};
#endif /* KOMPLEX_LINEARPROBLEM_H */

//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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
#include <EpetraExt_ConfigDefs.h>
class Epetra_Comm;
class Epetra_CrsMatrix;
class Epetra_Map;
namespace EpetraExt {
 
  //! Constructs an Epetra_CrsMatrix object from a Matrix Market format file, only row map specified.
  /*! This function constructs an Epetra_CrsMatrix object by reading a Matrix Market file.

      \param filename (In) A filename, including path if desired.  The matrix to be read should be in this file in 
                           Matrix Market coordinate format.

      \param rowMap (In) An Epetra_Map object describing the desired row distribution of the matrix.

      \param A (Out) An Epetra_CrsMatrix object constructed from file contents.  
      \warning User must delete!!.

      \return Returns 0 if no error, -1 if any problems with file system.

  */
  int MatrixMarketFileToCrsMatrix( const char *filename, const Epetra_Map & rowMap, Epetra_CrsMatrix * & A);
 
  //! Constructs an Epetra_CrsMatrix object from a Matrix Market format file, both row and column map specified.
  /*! This function constructs an Epetra_CrsMatrix object by reading a Matrix Market file.

      \param filename (In) A filename, including path if desired.  The matrix to be read should be in this file in 
                           Matrix Market coordinate format.

      \param rowMap (In) An Epetra_Map object describing the desired row distribution of the matrix.
      \param colMap (In) An Epetra_Map object describing the desired column distribution of the matrix.

      \param A (Out) An Epetra_CrsMatrix object constructed from file contents.  
      \warning User must delete!!.

      \return Returns 0 if no error, -1 if any problems with file system.

  */
  int MatrixMarketFileToCrsMatrix( const char *filename, const Epetra_Map & rowMap, const Epetra_Map & colMap, Epetra_CrsMatrix * & A);

  /** Reads an Epetra_CrsMatrix object from a matrix-market file, but
     uses the specified maps for constructing and 'FillComplete()'ing the
     matrix. Successfully creates rectangular matrices.
  */
  int MatrixMarketFileToCrsMatrix(const char *filename,
                                const Epetra_Map & rowMap,
                                const Epetra_Map & colMap,
                                const Epetra_Map& rangeMap,
                                const Epetra_Map& domainMap,
                                Epetra_CrsMatrix * & A);

  // Internal function
  int MatrixMarketFileToCrsMatrixHandle( const char *filename,
                                         Epetra_CrsMatrix * A,
                                         const Epetra_Map* rangeMap = NULL,
                                         const Epetra_Map* domainMap = NULL);

} // namespace EpetraExt

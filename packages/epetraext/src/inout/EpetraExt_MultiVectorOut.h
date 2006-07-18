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
class Epetra_MultiVector;
namespace EpetraExt {
 
  //! Writes an Epetra_MultiVector object to a Matrix Market format file
  /*! This function takes an Epetra_MultiVector object and writes it
      to the specified file.  The multivector can be distributed or serial.  The user can provide
      a strings containing the object name, a description, and specify that header information
      should or should not be printed to the file.

      \param filename (In) A filename, including path if desired.  If a file with this name already exists,
                      it will be deleted.  On exit, this file will contained any requested header information
		      followed by the matrix coefficients.  The file will contain a row for each entry.  All entries
		      for a column are listed before going to the next column.
      \param A (In) An Epetra_MultiVector Object containing the user matrix to be dumped to file.
      \param matrixName (In) A C-style string pointer to a name that will be stored in the comment field of the file.
                         This is not a required argument.  Note that it is possible to pass in the method A.Label().
      \param matrixDescription (In) A C-style string pointer to a matrix description that will be stored in the comment 
                                    field of the file.
      \param writeHeader (In) If true, the header will be written, otherwise only the matrix entries will be written.

      \return Returns 0 if no error, -1 if any problems with file system.

  */
  int MultiVectorToMatrixMarketFile( const char *filename, const Epetra_MultiVector & A, 
				   const char * matrixName=0,
				   const char *matrixDescription=0, 
				   bool writeHeader=true);

  //! Writes an Epetra_MultiVector object to a file that is compatible with Matlab.
  /*! This function takes any matrix that implements the Epetra_MultiVector interface and writes it
      to the specified file.  The matrix can be distributed or serial.  This function is a convenience wrapper 
      around MultiVectorToMatrixMarketFile.  The following Matlab commands can be used to read the resulting file
      and convert to it to a Matlab sparse matrix:
      <ol>
      <li> load \e filename;
      </ol>
      For example:
      <ol>
      <li> load A.dat;
      </ol>
      The above produces a dense matrix A with each vector in the multivector as a column in A.

      \param filename (In) A filename, including path if desired.  If a file with this name already exists,
                      it will be deleted.  On exit, this file will contain a row for each row of the multivector.
      \param A (In) An Epetra_MultiVector Object containing the user matrix to be dumped to file. 

      \return Returns 0 if no error, -1 if any problems with file system.

  */
  int MultiVectorToMatlabFile( const char *filename, const Epetra_MultiVector & A);
   

  //! Writes an Epetra_MultiVector object that is compatible with Matrix Market array format to a file handle.
  /*! This function takes an Epetra_MultiVector and writes it
      to the specified file handle.

      \param handle (In) A C-style file handle, already opened.  On exit, the file associated with this handle will
                      have appended to it a row for each multivector row.
      \param A (In) An Epetra_MultiVector Object containing the user object to be dumped to file.

      \return Returns 0 if no error, -1 if any problems with file system.

  */
  int MultiVectorToMatrixMarketHandle(FILE * handle, const Epetra_MultiVector & A);

  //! Writes an Epetra_MultiVector object that is compatible with Matlab to a file handle.
  /*! This function takes an Epetra_MultiVector and writes it
      to the specified file handle.

      \param handle (In) A C-style file handle, already opened.  On exit, the file associated with this handle will
                      have appended to it a row for each multivector row.
      \param A (In) An Epetra_MultiVector Object containing the user object to be dumped to file.

      \return Returns 0 if no error, -1 if any problems with file system.

  */
  int MultiVectorToMatlabHandle(FILE * handle, const Epetra_MultiVector & A);

  // Internal functions
  int MultiVectorToHandle(FILE * handle, const Epetra_MultiVector & A, bool mmFormat);
  int writeMultiVector(FILE * handle, const Epetra_MultiVector & A, bool mmFormat);

} // namespace EpetraExt

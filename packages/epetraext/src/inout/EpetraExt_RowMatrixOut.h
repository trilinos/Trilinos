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
#ifndef EPETRAEXT_ROWMATRIXOUT_H
#define EPETRAEXT_ROWMATRIXOUT_H
#include <EpetraExt_ConfigDefs.h>
class Epetra_RowMatrix;

namespace EpetraExt {
 
  //! Writes an Epetra_RowMatrix object to a Matrix Market format file
  /*! This function takes any matrix that implements the Epetra_RowMatrix interface and writes it
      to the specified file.  The matrix can be distributed or serial.  The user can provide
      a strings containing the matrix name, a matrix description, and specify that header information
      should or should not be printed to the file.

      \param filename (In) A filename, including path if desired.  If a file with this name already exists,
                      it will be deleted.  On exit, this file will contained any requested header information
		      followed by the matrix coefficients.  The file will contain a row for each matrix entry
		      The first column is the global row index, using base 1, the second column is the global
		      column index of the entry, the third value is the matrix value for that entry.
      \param A (In) An Epetra_RowMatrix Object containing the user matrix to be dumped to file.  Any object
                    that implements the Epetra_RowMatrix interface can be passed in.  In particular, the 
		    Epetra_CrsMatrix, Epetra_VbrMatrix, Epetra_FECrsMatrix, Epetra_FEVbrMatrix and Epetra_MsrMatrix
		    classes are compatible with this interface.
      \param matrixName (In) A C-style std::string pointer to a name that will be stored in the comment field of the file.
                         This is not a required argument.  Note that it is possible to pass in the method A.Label() if 
			 the matrix is one of the four types: Epetra_CrsMatrix, Epetra_VbrMatrix, Epetra_FECrsMatrix, 
			 Epetra_FEVbrMatrix.
      \param matrixDescription (In) A C-style std::string pointer to a matrix description that will be stored in the comment 
                                    field of the file.
      \param writeHeader (In) If true, the header will be written, otherwise only the matrix entries will be written.

      \return Returns 0 if no error, -1 if any problems with file system.

  */
  int RowMatrixToMatrixMarketFile( const char *filename, const Epetra_RowMatrix & A, 
				   const char * matrixName=0,
				   const char *matrixDescription=0, 
				   bool writeHeader=true);

  //! Writes an Epetra_RowMatrix object to a file that is compatible with Matlab.
  /*! This function takes any matrix that implements the Epetra_RowMatrix interface and writes it
      to the specified file.  The matrix can be distributed or serial.  This function is a convenience wrapper 
      around RowMatrixToMatrixMarketFile.  The following Matlab commands can be used to read the resulting file
      and convert to it to a Matlab sparse matrix:
      <ol>
      <li> load \e filename;
      <li> matrix_name = spconvert(filename_root);
      </ol>
      For example:
      <ol>
      <li> load A.dat;
      <li> A = spconvert(filename_root);
      </ol>
      The above produces a sparse matrix A.

      \param filename (In) A filename, including path if desired.  If a file with this name already exists,
                      it will be deleted.  On exit, this file will contain a row for each matrix entry
		      The first column is the global row index, using base 1, the second column is the global
		      column index of the entry, the third value is the matrix value for that entry.
      \param A (In) An Epetra_RowMatrix Object containing the user matrix to be dumped to file.  Any object
                    that implements the Epetra_RowMatrix interface can be passed in.  In particular, the 
		    Epetra_CrsMatrix, Epetra_VbrMatrix, Epetra_FECrsMatrix, Epetra_FEVbrMatrix and Epetra_MsrMatrix
		    classes are compatible with this interface.

      \return Returns 0 if no error, -1 if any problems with file system.

  */
  int RowMatrixToMatlabFile( const char *filename, const Epetra_RowMatrix & A);
   

  //! Writes an Epetra_RowMatrix object to a format file that is compatible with Matlab.
  /*! This function takes any matrix that implements the Epetra_RowMatrix interface and writes it
      to the specified file handle.  The matrix can be distributed or serial.  This function is a convenience wrapper 
      around RowMatrixToMatrixMarketFile.

      \param handle (In) A C-style file handle, already opened.  On exit, the file associated with this handle will
                      have appended to it a row for each matrix entry
		      The first column is the global row index, using base 1, the second column is the global
		      column index of the entry, the third value is the matrix value for that entry.
      \param A (In) An Epetra_RowMatrix Object containing the user matrix to be dumped to file.  Any object
                    that implements the Epetra_RowMatrix interface can be passed in.  In particular, the 
		    Epetra_CrsMatrix, Epetra_VbrMatrix, Epetra_FECrsMatrix, Epetra_FEVbrMatrix and Epetra_MsrMatrix
		    classes are compatible with this interface.

      \return Returns 0 if no error, -1 if any problems with file system.

  */
  int RowMatrixToHandle(std::FILE * handle, const Epetra_RowMatrix & A);
  int writeRowMatrix(std::FILE * handle, const Epetra_RowMatrix & A);

} // namespace EpetraExt
#endif /* EPETRAEXT_ROWMATRIXOUT_H */

//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
#ifndef EPETRAEXT_OPERATOROUT_H
#define EPETRAEXT_OPERATOROUT_H
#include <EpetraExt_ConfigDefs.h>
class Epetra_Operator;
class Epetra_MultiVector;
class Epetra_Map;

namespace EpetraExt {
 
  //! Writes an Epetra_Operator object to a Matrix Market format file, forming the coefficients by applying the operator to the e_j vectors.
  /*! This function takes any linear operator that implements the Epetra_Operator interface and writes it
      to the specified file.  The operator can be distributed or serial.  The user can provide
      a strings containing the matrix name, a matrix description, and specify that header information
      should or should not be printed to the file.

      The coeffients are formed by applying the operator to the canonical vectors 
      \f[ e_j = (0, \ldots, 0, 1, 0, \ldots, 0) \f]
      where the value 1 appears in the the jth entry.  The number of canonical vectors used is determined by the size of the 
      OperatorDomainMap() and the lengths by the size of OperatorRangeMap().

      \param filename (In) A filename, including path if desired.  If a file with this name already exists,
                      it will be deleted.  On exit, this file will contained any requested header information
		      followed by the matrix coefficients.  The file will contain a row for each matrix entry
		      The first column is the global row index, using base 1, the second column is the global
		      column index of the entry, the third value is the matrix value for that entry.
      \param A (In) An Epetra_Operator Object containing the implicit user matrix to be dumped to file.  Any object
                    that implements the Epetra_Operator interface can be passed in.  In particular, the 
		    Epetra_CrsMatrix, Epetra_VbrMatrix, Epetra_FECrsMatrix, Epetra_FEVbrMatrix and Epetra_MsrMatrix
		    classes are compatible with this interface, as is AztecOO_Operator, all Ifpack and ML preconditioners.
      \param matrixName (In) A C-style string pointer to a name that will be stored in the comment field of the file.
                         This is not a required argument.  Note that it is possible to pass in the method A.Label().
      \param matrixDescription (In) A C-style string pointer to a matrix description that will be stored in the comment 
                                    field of the file.
      \param writeHeader (In) If true, the header will be written, otherwise only the matrix entries will be written.

      \return Returns 0 if no error, -1 if any problems with file system.

  */
  int OperatorToMatrixMarketFile( const char *filename, const Epetra_Operator & A, 
				   const char * matrixName=0,
				   const char *matrixDescription=0, 
				   bool writeHeader=true);

  //! Writes an Epetra_Operator object to a file that is compatible with Matlab.
  /*! This function takes any matrix that implements the Epetra_Operator interface and writes it
      to the specified file.  The matrix can be distributed or serial.  This function is a convenience wrapper 
      around OperatorToMatrixMarketFile.  The following Matlab commands can be used to read the resulting file
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
      \param A (In) An Epetra_Operator Object containing the implicit user matrix to be dumped to file.  Any object
                    that implements the Epetra_Operator interface can be passed in.  In particular, the 
		    Epetra_CrsMatrix, Epetra_VbrMatrix, Epetra_FECrsMatrix, Epetra_FEVbrMatrix and Epetra_MsrMatrix
		    classes are compatible with this interface, as is AztecOO_Operator, all Ifpack and ML preconditioners.

      \return Returns 0 if no error, -1 if any problems with file system.

  */
  int OperatorToMatlabFile( const char *filename, const Epetra_Operator & A);
   

  //! Writes an Epetra_Operator object to a format file that is compatible with Matlab.
  /*! This function takes any matrix that implements the Epetra_Operator interface and writes it
      to the specified file handle.  The matrix can be distributed or serial.  This function is a convenience wrapper 
      around OperatorToMatrixMarketFile.

      \param handle (In) A C-style file handle, already opened.  On exit, the file associated with this handle will
                      have appended to it a row for each matrix entry
		      The first column is the global row index, using base 1, the second column is the global
		      column index of the entry, the third value is the matrix value for that entry.
      \param A (In) An Epetra_Operator Object containing the user matrix to be dumped to file.  Any object
                    that implements the Epetra_Operator interface can be passed in.  In particular, the 
		    Epetra_CrsMatrix, Epetra_VbrMatrix, Epetra_FECrsMatrix, Epetra_FEVbrMatrix and Epetra_MsrMatrix
		    classes are compatible with this interface, as is AztecOO_Operator, all Ifpack and ML preconditioners.

      \return Returns 0 if no error, -1 if any problems with file system.

  */
  int OperatorToHandle(std::FILE * handle, const Epetra_Operator & A);
  int writeOperatorStrip(std::FILE * handle, const Epetra_MultiVector & y, const Epetra_Map & rootDomainMap, const Epetra_Map & rootRangeMap, int startColumn);
  int get_nz(const Epetra_Operator & A, int & nz);

} // namespace EpetraExt
#endif /* EPETRAEXT_OPERATOROUT_H */

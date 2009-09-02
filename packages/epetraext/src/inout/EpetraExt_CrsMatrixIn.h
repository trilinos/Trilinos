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
#ifndef EPETRAEXT_CRSMATRIXIN_H
#define EPETRAEXT_CRSMATRIXIN_H
#include <EpetraExt_ConfigDefs.h>
class Epetra_Comm;
class Epetra_CrsMatrix;
class Epetra_Map;
namespace EpetraExt {

  //! Constructs an Epetra_CrsMatrix object from a Matlab format file, distributes rows evenly across processors.
  /*! This function constructs an Epetra_CrsMatrix object by reading a Matlab (i,j,value) format file.

      \param filename (In) A filename, including path if desired.  The matrix to be read should be in this file in 
                           Matlab coordinate format.
      \param comm (In) An Epetra_Comm object.  The matrix will have its rows distributed evenly by row-count across the parallel machine.
      \param A (Out) An Epetra_CrsMatrix object constructed from file contents.  The input matrix can be any dimension, square or rectangular.
      \warning User must delete matrix A!!.

      \return Returns 0 if no error, -1 if any problems with file system.

      Notes: 
<ol>
<li> The file will be read twice: first to get the maximum row and column dimensions.  Next to insert values.
<li> The global row and column dimensions will be determined by the maximum row and column index, respectively, 
     contained in the file.  If some rows or columns are empty they will still be present in the matrix.
</ol>

     The format expected for the input file is a list of nonzero entries with one entry per row.  Each row will have 
     the row index, column index and value listed with space in between each item.  The number of lines in the file should
     be exactly the number of entries of the matrix.  For example, consider the following matrix where only the nonzero values are stored:

\f[
\left[\begin{array}{cccc}
5 & 7 & 0 & 0 \\
3 & 2 & 0 & 1 \\
0 & 0 & 0 & 4 \\
\end{array}\right].
\f]

A Matlab format file for this matrix would be:
\verbatim
1 1 5.0
1 2 7.0
2 1 3.0
2 2 2.0
2 4 1.0
4 4 4.0
\endverbatim

Note that the entries can be listed in any order and that the matrix does not need to be square.  Values in the first and second columns
must be integer values and in those in the third column must be floating point format.


\htmlonly
(See the <a href="http://www.mathworks.com">Matlab</a> home page for details.)
\endhtmlonly
  */
  int MatlabFileToCrsMatrix( const char *filename, const Epetra_Comm & comm, Epetra_CrsMatrix * & A);

  //! Constructs an Epetra_CrsMatrix object from a Matrix Market format file, simplest version: requires matrix to be square, distributes rows evenly across processors.
  /*! This function constructs an Epetra_CrsMatrix object by reading a Matrix Market file.

      \param filename (In) A filename, including path if desired.  The matrix to be read should be in this file in 
                           Matrix Market coordinate format.
      \param comm (In) An Epetra_Comm object.
      \param transpose (In) A boolean value indicating whether the reader should transpose the matrix as it is read into matrix A.  (default = 0).
      \param verbose (In) A boolean value indicating whether the reader should print diagnostic statements to stdout.  (default = 0).
      \param A (Out) An Epetra_CrsMatrix object constructed from file contents.  
      \warning User must delete!!.

      \return Returns 0 if no error, -1 if any problems with file system.

\htmlonly
(See the <a href="http://math.nist.gov/MatrixMarket">Matrix Market</a> home page for details.)
\endhtmlonly

  */
  int MatrixMarketFileToCrsMatrix( const char *filename, const Epetra_Comm & comm, Epetra_CrsMatrix * & A, const bool transpose=0, const bool verbose=0);
 
  //! Constructs an Epetra_CrsMatrix object from a Matrix Market format file, row, range and domain map specified; typically used for rectangular matrices.
  /*! Reads an Epetra_CrsMatrix object from a matrix-market file, but
     uses the specified maps for constructing and 'FillComplete()'ing the
     matrix. Successfully creates rectangular matrices.

      \param filename (In) A filename, including path if desired.  The matrix to be read should be in this file in 
                           Matrix Market coordinate format.

      \param rowMap (In) An Epetra_Map object describing the desired row distribution of the matrix.
      \param rangeMap (In) An Epetra_Map object describing the distribution of range vectors that will be used with this matrix, must be 1-to-1.
      \param domainMap (In) An Epetra_Map object describing the distribution of domain vectors that will be used with this matrix, must be 1-to-1.

      \param transpose (In) A boolean value indicating whether the reader should transpose the matrix as it is read into matrix A.  (default = 0).
      \param verbose (In) A boolean value indicating whether the reader should print diagnostic statements to stdout.  (default = 0).
      \param A (Out) An Epetra_CrsMatrix object constructed from file contents.  
      \warning User must delete!!.

      \return Returns 0 if no error, -1 if any problems with file system.
\htmlonly
(See the <a href="http://math.nist.gov/MatrixMarket">Matrix Market</a> home page for details.)
\endhtmlonly
  */
  int MatrixMarketFileToCrsMatrix(const char *filename,const Epetra_Map & rowMap, 
				  const Epetra_Map& rangeMap, const Epetra_Map& domainMap, Epetra_CrsMatrix * & A, const bool transpose=0, const bool verbose=0);

  //! Constructs an Epetra_CrsMatrix object from a Matrix Market format file, only row map specified; allows user defined distribution of matrix rows, requires square matrix.
  /*! This function constructs an Epetra_CrsMatrix object by reading a Matrix Market file.

      \param filename (In) A filename, including path if desired.  The matrix to be read should be in this file in 
                           Matrix Market coordinate format.
      \param rowMap (In) An Epetra_Map object describing the desired row distribution of the matrix.
      \param transpose (In) A boolean value indicating whether the reader should transpose the matrix as it is read into matrix A.  (default = 0).
      \param verbose (In) A boolean value indicating whether the reader should print diagnostic statements to stdout.  (default = 0).
      \param A (Out) An Epetra_CrsMatrix object constructed from file contents.  
      \warning User must delete!!.

      \return Returns 0 if no error, -1 if any problems with file system.

\htmlonly
(See the <a href="http://math.nist.gov/MatrixMarket">Matrix Market</a> home page for details.)
\endhtmlonly

  */
  int MatrixMarketFileToCrsMatrix( const char *filename, const Epetra_Map & rowMap, Epetra_CrsMatrix * & A, const bool transpose=0, const bool verbose=0);
 
  //! Constructs an Epetra_CrsMatrix object from a Matrix Market format file, both row and column map specified; this version is seldom used unless you want explicit control over column map.
  /*! This function constructs an Epetra_CrsMatrix object by reading a Matrix Market file.

      \param filename (In) A filename, including path if desired.  The matrix to be read should be in this file in 
                           Matrix Market coordinate format.

      \param rowMap (In) An Epetra_Map object describing the desired row distribution of the matrix.
      \param colMap (In) An Epetra_Map object describing the desired column distribution of the matrix.

      \param transpose (In) A boolean value indicating whether the reader should transpose the matrix as it is read into matrix A.  (default = 0).
      \param verbose (In) A boolean value indicating whether the reader should print diagnostic statements to stdout.  (default = 0).
      \param A (Out) An Epetra_CrsMatrix object constructed from file contents.  
      \warning User must delete!!.

      \return Returns 0 if no error, -1 if any problems with file system.

\htmlonly
(See the <a href="http://math.nist.gov/MatrixMarket">Matrix Market</a> home page for details.)
\endhtmlonly
  */
  int MatrixMarketFileToCrsMatrix( const char *filename, const Epetra_Map & rowMap, const Epetra_Map & colMap, Epetra_CrsMatrix * & A, const bool transpose=0, const bool verbose=0);


  //! Constructs an Epetra_CrsMatrix object from a Matrix Market format file, row, column, range and domain map specified; this version is seldom required unless you want explicit control over column map.
  /*! Reads an Epetra_CrsMatrix object from a matrix-market file, but
     uses the specified maps for constructing and 'FillComplete()'ing the
     matrix. Successfully creates rectangular matrices.

      \param filename (In) A filename, including path if desired.  The matrix to be read should be in this file in 
                           Matrix Market coordinate format.

      \param rowMap (In) An Epetra_Map object describing the desired row distribution of the matrix.
      \param colMap (In) An Epetra_Map object describing the desired column distribution of the matrix.
      \param rangeMap (In) An Epetra_Map object describing the distribution of range vectors that will be used with this matrix, must be 1-to-1.
      \param domainMap (In) An Epetra_Map object describing the distribution of domain vectors that will be used with this matrix, must be 1-to-1.
      \param transpose (In) A boolean value indicating whether the reader should transpose the matrix as it is read into matrix A.  (default = 0).
      \param verbose (In) A boolean value indicating whether the reader should print diagnostic statements to stdout.  (default = 0).

      \param A (Out) An Epetra_CrsMatrix object constructed from file contents.  
      \warning User must delete!!.

      \return Returns 0 if no error, -1 if any problems with file system.
\htmlonly
(See the <a href="http://math.nist.gov/MatrixMarket">Matrix Market</a> home page for details.)
\endhtmlonly
  */
  int MatrixMarketFileToCrsMatrix(const char *filename, const Epetra_Map & rowMap, const Epetra_Map & colMap,
				  const Epetra_Map& rangeMap, const Epetra_Map& domainMap, Epetra_CrsMatrix * & A, const bool transpose=0, const bool verbose=0);

  //! Constructs an Epetra_CrsMatrix object from a Hypre Matrix Print command, the row map is specified.
  /*! Reads an Epetra_CrsMatrix object from a Hypre Matrix Printout, the matrix should be square.

      \param filename (In) A filename not including the processor id extension, including path if desired.

      \param comm (In) An Epetra_Comm object describing the communication among processors.

      \param A (Out) An Epetra_CrsMatrix object constructed from file contents.  
      \warning User must delete!!.

      \return Returns 0 if no error, -1 if any problems with file system.
  */
  int HypreFileToCrsMatrix(const char *filename, const Epetra_Comm &comm, Epetra_CrsMatrix *&A);
  // Internal function
  int MatrixMarketFileToCrsMatrixHandle( const char *filename,
					 const Epetra_Comm & comm,
                                         Epetra_CrsMatrix * & A,
					 const Epetra_Map * rowMap = 0,
					 const Epetra_Map * colMap = 0,
					 const Epetra_Map * rangeMap = 0,
					 const Epetra_Map * domainMap = 0,
                                         const bool transpose = 0, 
                                         const bool verbose=0);

} // namespace EpetraExt
#endif /* EPETRAEXT_CRSMATRIXIN_H */

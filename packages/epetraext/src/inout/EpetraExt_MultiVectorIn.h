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
#ifndef EPETRAEXT_MULTIVECTORIN_H
#define EPETRAEXT_MULTIVECTORIN_H
#include <EpetraExt_ConfigDefs.h>
class Epetra_MultiVector;
class Epetra_BlockMap;
namespace EpetraExt {
 
  //! Constructs an Epetra_MultiVector object from a Matrix Market format file.
  /*! This function constructs an Epetra_MultiVector object by reading a Matrix Market file.

      \param filename (In) A filename, including path if desired.  The multivector to be read should be in this file in 
                           Matrix Market array format.

      \param map (In) An Epetra_Map or Epetra_BlockMap object describing the desired distribution of the multivector.

      \param A (Out) An Epetra_MultiVector object constructed from file contents.  
      \warning User must delete!!.

      \return Returns 0 if no error, -1 if any problems with file system.

  */
  int MatrixMarketFileToMultiVector( const char *filename, const Epetra_BlockMap & map, Epetra_MultiVector * & A);

} // namespace EpetraExt
#endif /* EPETRAEXT_MULTIVECTORIN_H */

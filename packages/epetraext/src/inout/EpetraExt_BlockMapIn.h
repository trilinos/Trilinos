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
#ifndef EPETRAEXT_BLOCKMAPIN_H
#define EPETRAEXT_BLOCKMAPIN_H

#include <EpetraExt_ConfigDefs.h>
class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_Map;
  /*! \namespace EpetraExt
     \brief The EpetraExt namespace contains a variety of useful functions and class that extend Epetra capabilities.

     All EpetraExt capabilities are available from the EpetraExt package in Trilinos.  
  */

namespace EpetraExt {
 
  //! Constructs an Epetra_BlockMap object from a Matrix Market format file
  /*! This function constructs an Epetra_BlockMap or Epetra_Map object by reading a Matrix Market file.
      If the file was created using the EpetraExt::BlockMapOut functions,
      special information was encoded in the comment field of this map that allows for identical reproduction
      of the map, including distribution across processors and element size information.
      If the same of processors is being used to create the object as were used to write it, the object will be an exact
      reproduction of the original.  Otherwise, a uniform distribution of the GIDs will be created.

      The first column of the input file will must be the list of GIDs in the map.
      If the block map has non-uniform sizes, a second column must contain the element sizes.

      \param filename (In) A filename, including path if desired.  If a file with this name already exists,
                      it will be deleted.  On exit, this file will contained any requested header information
		      followed by the map GIDs.  A second column may be present if the BlockMap has nonuniform sizes.

      \param comm (In) An Epetra_Comm object describing the parallel machine.

      \param blockMap (Out) An Epetra_BlockMap object constructed from file contents.  
      \warning User must delete!!.

      \return Returns 0 if no error, -1 if any problems with file system, returns 1 if number of processors differs from file creator.

  */
  int MatrixMarketFileToBlockMap( const char *filename, const Epetra_Comm & comm, Epetra_BlockMap * & blockMap);

  //! Constructs an Epetra_BlockMap object from a Matrix Market format file
  /*! This function constructs an Epetra_BlockMap or Epetra_Map object by reading a Matrix Market file.
      If the file was created using the EpetraExt::BlockMapOut functions,
      special information was encoded in the comment field of this map that allows for identical reproduction
      of the map, including distribution across processors and element size information.
      If the same of processors is being used to create the object as were used to write it, the object will be an exact
      reproduction of the original.  Otherwise, a uniform distribution of the GIDs will be created.

      The first column of the input file will must be the list of GIDs in the map.
      If the block map has non-uniform sizes, a second column must contain the element sizes.

      \param filename (In) A filename, including path if desired.  If a file with this name already exists,
                      it will be deleted.  On exit, this file will contained any requested header information
		      followed by the map GIDs.  A second column may be present if the BlockMap has nonuniform sizes.

      \param comm (In) An Epetra_Comm object describing the parallel machine.

      \param map (Out) An Epetra_Map object constructed from file contents.

      \warning User must delete!!.


      \return Returns 0 if no error, -1 if any problems with file system, -2 if file contained nontrivial Epetra_BlockMap, 1 if number of processors differs from file creator.

  */
  int MatrixMarketFileToMap( const char *filename, const Epetra_Comm & comm, Epetra_Map * & map);


  /** Constructs row,col,range and domain maps from a matrix-market matrix file.
  */
  int MatrixMarketFileToBlockMaps(const char* filename,
                                  const Epetra_Comm& comm,
                                  Epetra_BlockMap*& rowmap,
                                  Epetra_BlockMap*& colmap,
                                  Epetra_BlockMap*& rangemap,
                                  Epetra_BlockMap*& domainmap);
} // namespace EpetraExt
#endif /* EPETRAEXT_BLOCKMAPIN_H */

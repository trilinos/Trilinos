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

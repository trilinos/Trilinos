//@HEADER
// ************************************************************************
//
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
                                                                                                    
#ifndef EpetraExt_LINEARPROBLEM_REINDEX_H
#define EpetraExt_LINEARPROBLEM_REINDEX_H

#include <EpetraExt_Transform.h>

class Epetra_Map;
class Epetra_LinearProblem;

namespace EpetraExt {

class CrsMatrix_Reindex;
class MultiVector_Reindex;

///
/** Given and input Epetra_LinearProblem, a "reindexed" version will be returned
 *  using the given NewRowMap.  If a null map is given, a lexigraphically indexed
 *  LP will be returned.  The data in the new E_LP is a "reindexed" view of the 
 *  original.
 */
class LinearProblem_Reindex : public ViewTransform<Epetra_LinearProblem>
{
  CrsMatrix_Reindex * MatTrans_;
  MultiVector_Reindex * LHSTrans_;
  MultiVector_Reindex * RHSTrans_;

  Epetra_Map * NewRowMap_;

  bool NewRowMapOwned_;

 public:

  ///
  /** Destructor
   */
  ~LinearProblem_Reindex();

  ///
  /** Constructor
   */
  LinearProblem_Reindex( Epetra_Map * NewRowMap )
  : NewRowMap_(NewRowMap),
    MatTrans_(0),
    LHSTrans_(0),
    RHSTrans_(0),
    NewRowMapOwned_(false)
  {}

  ///
  /** Constructs a new view the original LP, "reindexed" using the given NewRowMap.
   */
  NewTypeRef operator()( OriginalTypeRef orig );

};

} //namespace EpetraExt

#endif //EpetraExt_LINEARPROBLEM_REINDEX_H


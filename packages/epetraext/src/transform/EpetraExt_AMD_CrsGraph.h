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

#ifndef EpetraExt_CRSGRAPH_AMD_H
#define EpetraExt_CRSGRAPH_AMD_H

#include <EpetraExt_Transform.h>

class Epetra_CrsGraph;
class Epetra_Map;

namespace EpetraExt {

///
/** Approximate Minimum Degree Reordering using Tim Daley's AMD Algorithm.
 */
class CrsGraph_AMD : public StructuralSameTypeTransform<Epetra_CrsGraph> {

 public:

  ///
  /** Destructor
   */
  ~CrsGraph_AMD();

  ///
  /** Constructor
   */
  CrsGraph_AMD( bool verbose = false )
  : NewMap_(0),
    NewGraph_(0),
    verbose_(verbose)
  {}

  ///
  /** Constructs AMD ordered Epetra_CrsGraph from <tt>orig</tt> object.
   */
  NewTypeRef operator()( OriginalTypeRef orig );

 private:

  Epetra_Map * NewMap_;
  Epetra_CrsGraph * NewGraph_;

  const bool verbose_;
};

} //namespace EpetraExt

#endif //EpetraExt_CRSGRAPH_AMD_H

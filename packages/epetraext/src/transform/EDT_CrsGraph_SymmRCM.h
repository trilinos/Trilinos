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

#ifndef EDT_CRSGRAPH_SYMMRCM_H
#define EDT_CRSGRAPH_SYMMRCM_H

#include <vector>

#include <Epetra_Transform.h>

class Epetra_Map;
class Epetra_CrsGraph;

namespace EpetraExt {

struct CrsGraph_SymmRCM : public StructuralSameTypeTransform<Epetra_CrsGraph> {

 public:

  ~CrsGraph_SymmRCM();

  CrsGraph_SymmRCM( int testLeafWidth = 5 )
  : testLeafWidth_(testLeafWidth),
    RCMMap_(0)
  {}

  NewTypeRef operator()( OriginalTypeRef orig );

 private:

  Epetra_Map * RCMMap_;
  int testLeafWidth_;

  class BFT {
    
   public:

     BFT( const std::vector< std::vector<int> > & adjlist,
          int root,
          int max_width,
          bool & failed );

     int Width() { return width_; }
     int Depth() { return depth_; }

     void NonNeighborLeaves( std::vector<int> & leaves, int count );
     void ReverseVector( std::vector<int> & ordered );

   private:

     bool failed_;
     int width_;
     int depth_;
     int nodes_;

     std::vector< std::vector<int> > levelSets_;
     std::vector< std::vector<int> > adjList_;

  };

};

} //namespace EpetraExt

#endif //EDT_CRSGRAPH_SYMMRCM_H

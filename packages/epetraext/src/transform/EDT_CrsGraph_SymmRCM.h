
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

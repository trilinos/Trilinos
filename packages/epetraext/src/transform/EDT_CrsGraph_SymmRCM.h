
#ifndef EDT_CRSGRAPH_SYMMRCM_H
#define EDT_CRSGRAPH_SYMMRCM_H

#include <Epetra_Transform.h>

class Epetra_CrsGraph;

namespace EpetraExt {

struct CrsGraph_SymmRCM : public StructuralSameTypeTransform<Epetra_CrsGraph> {

 public:

  ~CrsGraph_SymmRCM() {}

  CrsGraph_SymmRCM() {}

  NewTypePtr operator()( OriginalTypeRef original );

 private:

  class BFT {
    
   public:

     BFT( const vector< vector<int> > & adjlist,
          int root,
          int max_width,
          bool & failed );

     int Width() { return width_; }
     int Depth() { return depth_; }

     void NonNeighborLeaves( vector<int> & leaves, int count );
     void ReverseVector( vector<int> & ordered );

   private:

     int width_;
     int depth_;
     int nodes_;

     vector< vector<int> > levelSets_;
     vector< vector<int> > adjList_;

  };

};

} //namespace EpetraExt

#endif //EDT_CRSGRAPH_SYMMRCM_H

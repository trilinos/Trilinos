#include "Tpetra_Core.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_FEMultiVector.hpp"
//#include "Zoltan2_IceSheetGraph.h"

#include <string>
#include <sstream>
#include <iostream>
#include <queue>

#ifndef ZOLTAN2_VTXLABEL_
#define ZOLTAN2_VTXLABEL_

namespace Zoltan2{
        
        /*
 *
struct graph {
  int n;
  unsigned m;
  int* out_array;
  unsigned* out_degree_list;
  int max_degree_vert;
  double avg_out_degree;
};

#define out_degree(g, n) (g->out_degree_list[n+1] - g->out_degree_list[n])
#define out_vertices(g, n) (&g->out_array[g->out_degree_list[n]])

 *      */
        template<typename lno_t>
	class icePropGraph {
          public:
            lno_t n;
            lno_t m;
            lno_t* out_array;
            lno_t* out_degree_list;
            lno_t max_degree_vert;
            double avg_out_degree;
            lno_t out_degree(lno_t vert){
              return (out_degree_list[vert+1] - out_degree_list[vert]);
            }
            lno_t* out_vertices(lno_t vert){
              return (&out_array[out_degree_list[vert]]);
            }
        }; 
        
        
	std::queue<int> icePropArtQueue;
	std::queue<int> icePropRegQueue;


	enum IcePropGrounding_Status {ICEPROPGS_FULL=2, ICEPROPGS_HALF=1, ICEPROPGS_NONE = 0};
	// Struct representing a vertex label.
	// We define our own "addition" for these labels.
	// Later, we'll create a Tpetra::FEMultiVector of these labels.
	class IcePropVtxLabel {
	public:
	  //int gid;
          int id;
	  int first_label;
	  int first_sender;
	  int second_label;
	  int second_sender;
          int bcc_name;
	  bool is_art;
	  // Constructors
	  IcePropVtxLabel(int idx_, int first_ = -1, int first_sender_ = -1, int second_ = -1, int second_sender_ = -1,bool art_ = false, int bcc_name_ = -1) { 
	    id = idx_;
            //gid = idx;
	    if(id == -1){
		std::cout<<"A label's ID is -1\n";
	    }
	    first_label = first_;
	    first_sender = first_sender_;
	    second_label = second_;
	    second_sender = second_sender_; 
	    is_art = art_;
            bcc_name = bcc_name_;
	  }
	  IcePropVtxLabel() {
	    id = -1;
	    first_label = -1;
	    first_sender = -1;
	    second_label = -1;
	    second_sender = -1;
	    is_art = false; 
            bcc_name = -1;
	  }

	  IcePropVtxLabel(volatile const IcePropVtxLabel& other){
	    id = other.id;
            //gid = other.gid;
	    first_label = other.first_label;
	    first_sender = other.first_sender;
	    second_label = other.second_label;
	    second_sender = other.second_sender;
	    is_art = other.is_art;
            bcc_name = other.bcc_name;
	  }
	  IcePropVtxLabel(const IcePropVtxLabel& other){
	    id = other.id;
            //gid = other.gid;
	    first_label = other.first_label;
	    first_sender = other.first_sender;
	    second_label = other.second_label;
	    second_sender = other.second_sender;
	    is_art = other.is_art;
            bcc_name = other.bcc_name;
	  }
	  // IcePropVtxLabel assignment
          volatile IcePropVtxLabel operator=(const IcePropVtxLabel& other) volatile{
	    id = other.id;
            //gid = other.gid;
	    first_label = other.first_label;
	    first_sender = other.first_sender;
	    second_label = other.second_label;
	    second_sender = other.second_sender;
	    is_art = other.is_art;
            bcc_name = other.bcc_name;
	    return *this; 
	  }

          IcePropVtxLabel& operator=(const IcePropVtxLabel& other) {
	    id = other.id;
            //gid = other.gid;
	    first_label = other.first_label;
	    first_sender = other.first_sender;
	    second_label = other.second_label;
	    second_sender = other.second_sender;
	    is_art = other.is_art;
            bcc_name = other.bcc_name;
	    return *this; 
	  } 

	  // int assignment
	  IcePropVtxLabel& operator=(const int& other) { 
	    first_label = other;
	    first_sender = other;
	    second_label = -1;
	    second_sender = -1;
	    return *this;
	  }
	  // += overload
	  // for communicating copy's labels over processor boundaries.
          IcePropVtxLabel& operator+=(const IcePropVtxLabel& copy) {
            IcePropGrounding_Status owned_gs = getGroundingStatus();
            IcePropGrounding_Status copy_gs = copy.getGroundingStatus();
            //The only cases we care about are 
            //owned	copy	
            //NONE  <   HALF
            //NONE  <	FULL
            //HALF  ==	HALF
            //HALF  <	FULL
          
	    //handles NONE < HALF, HALF < FULL
            if(owned_gs < copy_gs){
              first_label = copy.first_label;
              first_sender = copy.first_sender;
              second_label = copy.second_label;
              second_sender = copy.second_sender;
              bcc_name = copy.bcc_name;
            //handles HALF == HALF
            } else if(owned_gs == copy_gs && owned_gs == ICEPROPGS_HALF){
              if(copy.first_label != first_label){
                second_label = copy.first_label;
                second_sender = copy.first_sender;
              }
            }
            
            if(getGroundingStatus() != owned_gs){
              if(!is_art) icePropRegQueue.push(id);
              else icePropArtQueue.push(id);
            }
            
	    return *this;
	  }
	  
	  // IcePropVtxLabel equality overload
	  friend bool operator==(const IcePropVtxLabel& lhs, const IcePropVtxLabel& rhs) {
	    return ((lhs.first_label == rhs.first_label)&&(lhs.first_sender == rhs.first_sender)&&(lhs.second_label == rhs.second_label)&&(lhs.second_sender == rhs.second_sender));
	  }
	  // int equality overload
	  friend bool operator==(const IcePropVtxLabel& lhs, const int& rhs) {
	    return ((lhs.first_label == rhs)&&(lhs.first_sender == rhs));
	  }
	  // output stream overload
	  friend std::ostream& operator<<(std::ostream& os, const IcePropVtxLabel& a) {
	    os<<a.id<<": "<< a.first_label<<", "<<a.first_sender<<"; "<<a.second_label<<", "<<a.second_sender;
	    return os;
	  }
	  IcePropGrounding_Status getGroundingStatus() const {
	    return (IcePropGrounding_Status)((first_label != -1) + (second_label != -1));
	  }
	};

}//end namespace iceProp
        
/////////////////////////////////////////////////////////////////////////
// ArithTraits -- arithmetic traits needed for struct IcePropVtxLabel
// Needed so that Tpetra compiles.
// Not all functions were needed; this is a subset of ArithTraits' traits.
// Modified from kokkos-kernels/src/Kokkos_ArithTraits.hpp's 
// <int> specialization
namespace Kokkos {
  namespace Details {

    template<>
    class ArithTraits<Zoltan2::IcePropVtxLabel> {  // specialized for IcePropVtxLabel struct
    public:
      typedef Zoltan2::IcePropVtxLabel val_type;
      typedef int mag_type;
    
      static const bool is_specialized = true;
      static const bool is_signed = true;
      static const bool is_integer = true;
      static const bool is_exact = true;
      static const bool is_complex = false;
    
      static KOKKOS_FORCEINLINE_FUNCTION bool isInf(const val_type &) {
	return false;
      }
      static KOKKOS_FORCEINLINE_FUNCTION bool isNan(const val_type &) {
	return false;
      }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type abs(const val_type &x) {
	return (x.first_label >= 0 ? x.first_label : -(x.first_label));
      }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type zero() { return 0; }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type one() { return 1; }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type min() { return INT_MIN; }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type max() { return INT_MAX; }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type nan() { return -1; }
    
      // Backwards compatibility with Teuchos::ScalarTraits.
      typedef mag_type magnitudeType;
      static const bool isComplex = false;
      static const bool isOrdinal = true;
      static const bool isComparable = true;
      static const bool hasMachineParameters = false;
      static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude(
	const val_type &x) 
      {
	return abs(x);
      }
      static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf(const val_type &) {
	return false;
      }
      static std::string name() { return "Zoltan2::IcePropVtxLabel"; }
    };
  }
}

/////////////////////////////////////////////////////////////////////////////
// Teuchos::SerializationTraits are needed to copy vtxLabels into MPI buffers
// Because sizeof(vtxLabel) works for struct vtxLabel, we'll use a 
// provided serialization of vtxLabel into char*.
namespace Teuchos{
template<typename Ordinal>
struct SerializationTraits<Ordinal, Zoltan2::IcePropVtxLabel> :
       public Teuchos::DirectSerializationTraits<Ordinal, Zoltan2::IcePropVtxLabel>
{};
}//end namespace Teuchos

namespace Zoltan2{
template<typename MAP>
class iceSheetPropagation {
public:
  //typedefs for the FEMultiVector
  //typedef Tpetra::Map<> map_t;
  typedef typename MAP::local_ordinal_type lno_t;
  typedef typename MAP::global_ordinal_type gno_t;
  typedef IcePropVtxLabel scalar_t;
  typedef Tpetra::FEMultiVector<scalar_t,lno_t, gno_t> femv_t;	
  
  
 	
  //Constructor assigns vertices to processors and builds maps with and
  //without copies ICE SHEET VERSION
  iceSheetPropagation(const Teuchos::RCP<const Teuchos::Comm<int> > &comm_, Teuchos::RCP<const MAP> mapOwned_, Teuchos::RCP<const MAP> mapWithCopies_, icePropGraph<lno_t>* g_,int* boundary_flags, bool* grounding_flags,int localOwned,int localCopy):
    me(comm_->getRank()), np(comm_->getSize()),
    nLocalOwned(localOwned), nLocalCopy(localCopy),
    nVec(1), comm(comm_),g(g_),mapOwned(mapOwned_),
    mapWithCopies(mapWithCopies_)
  {
    typedef Tpetra::Import<lno_t, gno_t> import_t;
    Teuchos::RCP<import_t> importer = rcp(new import_t(mapOwned, mapWithCopies));
    femv = rcp(new femv_t(mapOwned, importer, nVec, true)); 
    femv->beginFill();
    //set the member variable that stores femv->getData(0);
    femvData = femv->getData(0);
    for(lno_t i = 0; i < nLocalOwned + nLocalCopy; i++){
      gno_t gid = mapWithCopies->getGlobalElement(i);
      IcePropVtxLabel label(i);
      if(boundary_flags[i] > 2 && i < nLocalOwned) {
	label.is_art = true;
      }
      if(grounding_flags[i]){
        label.first_label = gid;
        label.first_sender = gid;
        icePropRegQueue.push(label.id);
        //if(label.is_art) icePropArtQueue.push(label.id);
        //else icePropRegQueue.push(label.id);
      }
      femv->replaceLocalValue(i,0,label);
      //femv->replaceGlobalValue(gid,1,me);
    }
    //IcePropVtxLabel label1 = femv->getData(0)[0];
    //IcePropVtxLabel label2 = femv->getData(0)[1];
    //label1 += label2;

    //femv->replaceGlobalValue(0,0,label1);
    //printFEMV("BeforeFill");   
 
    femv->endFill(); 
    std::cout<<me<<": finished constructing successfully\n";
    //printFEMV("AfterFill");
  }
  void printFEMV(const char *msg){
    for (int v = 0; v < 1; v++){
      std::cout << me << " OWNED " << msg << " FEMV[" << v << "] Owned: ";
      auto value = femv->getData(v);
      for(lno_t i = 0; i < nLocalOwned; i++) std::cout<<value[i]<<" ";
      std::cout<<std::endl;
    }
    femv->switchActiveMultiVector();
    for(int v = 0; v < 1; v++){
      std::cout << me << " WITHCOPIES " <<msg<<" FEMV[" <<v<< "] Owned: ";
      auto value = femv->getData(v);
      for(lno_t i = 0; i < nLocalOwned; i++) std::cout << value[i] <<" ";
      std::cout<<" Copies: ";
      for(lno_t i = nLocalOwned; i < nLocalOwned+nLocalCopy; i++)
        std::cout << value[i] <<" ";
      std::cout << std::endl;
    }
    femv->switchActiveMultiVector();

  }
  //propagation functions
  //returns an array of vertices to remove
  int* propagate(void){ 
    //run bfs_prop
    //std::cout<<me<<": starting propagation\n"; 
    bfs_prop();
    //std::cout<<me<<": done with initial propagation\n";
    //check for potentially false articulation points
    while(true){
      femv->switchActiveMultiVector(); 
      for(int i = 0; i < g->n; i++){
	IcePropVtxLabel curr_node = femvData[i];
        if(curr_node.is_art && curr_node.getGroundingStatus() == ICEPROPGS_FULL){
          lno_t out_degree = g->out_degree(curr_node.id);//out_degree(g, curr_node.id);
          lno_t* outs = g->out_vertices(curr_node.id);//out_vertices(g, curr_node.id);
          for(int j = 0; j < out_degree; j++){
            IcePropVtxLabel neighbor = femvData[outs[j]];
            if(neighbor.getGroundingStatus() == ICEPROPGS_HALF && neighbor.first_label != mapWithCopies->getGlobalElement(curr_node.id) && neighbor.first_sender == mapWithCopies->getGlobalElement(curr_node.id)){
              icePropRegQueue.push(curr_node.id);
            }
          }
        }
      }
      int local_done = icePropRegQueue.empty();
      int done = 0;
      Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_MIN,1, &local_done,&done);
      
      if(done) break;
      
      //add all articulation points to the icePropRegQueue
      for(int i = 0; i < g->n; i++){
	IcePropVtxLabel curr_node = femvData[i];
        if(curr_node.getGroundingStatus() == ICEPROPGS_HALF){
	  IcePropVtxLabel cleared(curr_node.id);
          cleared.is_art = curr_node.is_art;
          femv->replaceLocalValue(i,0,cleared);
        }
        if(curr_node.getGroundingStatus() == ICEPROPGS_FULL && curr_node.is_art){
          icePropRegQueue.push(curr_node.id);
        }
      }
      //std::cout<<me<<": Running BFS-prop again\n";
      //re-run bfs_prop until incomplete propagation is fixed
      bfs_prop();     
    }
    //check for nodes that are less than full.
    //return flags for each node, -2 for keep, -1 for remove, <vtxID> for singly grounded nodes.
    lno_t* removed = new lno_t[g->n];
    for(int i = 0; i < g->n; i++){
      IcePropVtxLabel curr_node = femvData[i];
      IcePropGrounding_Status gs = curr_node.getGroundingStatus();
      if(gs == ICEPROPGS_FULL) removed[i] = -2;
      else if(gs == ICEPROPGS_HALF) removed[i] = curr_node.first_label;
      else removed[i] = -1;
    }
    femv->switchActiveMultiVector();
    //std::cout<<me<<": returning answer\n";
    return removed;
  }

  //performs one bfs_prop iteration (does not check for incomplete propagation)
  void bfs_prop(){
    //do
    //  femv->beginFill()
    //  //call giveLabels in a while loop
    //  femv->endFill()
    //  reduceAll()
    //while(!done);
    int done = 0;
    while(!done){
      //std::cout<<me<<": Propagating...\n";
      femv->beginFill();
      //visit every node that changed
      std::queue<int>* curr;
      if(icePropRegQueue.empty()) curr = &icePropArtQueue;
      else curr = &icePropRegQueue;
      while(!curr->empty()){
	IcePropVtxLabel curr_node = femvData[curr->front()];
        curr->pop();
        
        //if the current node is a copy, it shouldn't propagate out to its neighbors.
        if(curr_node.id >= nLocalOwned) continue;

        lno_t out_degree = g->out_degree(curr_node.id);//out_degree(g, curr_node.id);
        lno_t* outs = g->out_vertices(curr_node.id);//out_vertices(g, curr_node.id);
        for(int i = 0; i < out_degree; i++){
	  IcePropVtxLabel neighbor = femvData[outs[i]];
	  IcePropGrounding_Status old_gs = neighbor.getGroundingStatus();
          
          //give curr_node's neighbor some more labels
          giveLabels(curr_node, neighbor);
          
          if(old_gs != neighbor.getGroundingStatus()){
            femv->replaceLocalValue(outs[i],0,neighbor);
            if(neighbor.is_art) icePropArtQueue.push(neighbor.id);
            else icePropRegQueue.push(neighbor.id);
          }
        }
        if(curr->empty()){
          if(curr == &icePropRegQueue) curr = &icePropArtQueue;
          else curr = &icePropRegQueue;
        }
        
      }
      //std::cout<<me<<": art queue front = "<<icePropArtQueue.front()<<"\n";
      //std::cout<<me<<": art queue size = "<<icePropArtQueue.size()<<"\n";
      femv->endFill();
      //printFEMV("Before Communicating");
      femv->doOwnedToOwnedPlusShared(Tpetra::ADD);
      femv->doOwnedPlusSharedToOwned(Tpetra::ADD);
      //std::cout<<me<<": reg queue front = "<<icePropRegQueue.front()<<"\n";
      //std::cout<<me<<": reg queue size = "<<icePropRegQueue.size()<<"\n";
      int local_done = icePropRegQueue.empty() && icePropArtQueue.empty();
      //printFEMV("After Communicating");
      //this call makes sure that if any inter-processor communication changed labels
      //we catch the changes and keep propagating them.
      Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_MIN,1, &local_done,&done);
    }
    
    
  }
  
  //function that exchanges labels between two nodes
  //curr_node gives its labels to neighbor.
  void giveLabels(IcePropVtxLabel& curr_node, IcePropVtxLabel& neighbor){
    IcePropGrounding_Status curr_gs = curr_node.getGroundingStatus();
    IcePropGrounding_Status nbor_gs = neighbor.getGroundingStatus();
    int curr_node_gid = mapWithCopies->getGlobalElement(curr_node.id);
    //if the neighbor is full, we don't need to pass labels
    if(nbor_gs == ICEPROPGS_FULL) return;
    //if the current nod is empty (shouldn't happen) we can't pass any labels
    if(curr_gs == ICEPROPGS_NONE) return;
    //if the current node is full (and not an articulation point), pass both labels on
    if(curr_gs == ICEPROPGS_FULL && !curr_node.is_art){
      neighbor.first_label = curr_node.first_label;
      neighbor.first_sender = curr_node_gid;
      neighbor.second_label = curr_node.second_label;
      neighbor.second_sender = curr_node_gid;
      neighbor.bcc_name = curr_node.bcc_name;
      return;
    } else if (curr_gs == ICEPROPGS_FULL) {
      //if it is an articulation point, and it hasn't sent to this neighbor
      if(neighbor.first_sender != curr_node_gid){
        //send itself as a label
        if(nbor_gs == ICEPROPGS_NONE){
          neighbor.first_label = curr_node_gid;
          neighbor.first_sender = curr_node_gid;
          neighbor.bcc_name = curr_node.bcc_name;
        } else if(nbor_gs == ICEPROPGS_HALF){
          if(neighbor.first_label != curr_node_gid){
            neighbor.second_label = curr_node_gid;
            neighbor.second_sender = curr_node_gid;
          }
        }
      }
      return;
    }
    //if the current node has only one label
    if(curr_gs == ICEPROPGS_HALF){
      //pass that on appropriately
      if(nbor_gs == ICEPROPGS_NONE){
        neighbor.first_label = curr_node.first_label;
        neighbor.first_sender = curr_node_gid;
        neighbor.bcc_name = curr_node.bcc_name;
      } else if(nbor_gs == ICEPROPGS_HALF){
        //make sure you aren't giving a duplicate label, and that
        //you haven't sent a label to this neighbor before.
        if(neighbor.first_label != curr_node.first_label && neighbor.first_sender != curr_node_gid){
          neighbor.second_label = curr_node.first_label;
          neighbor.second_sender = curr_node_gid;
        }
      }
    }
    
    if(nbor_gs != neighbor.getGroundingStatus()){
      if(neighbor.is_art) icePropArtQueue.push(neighbor.id);
      else icePropRegQueue.push(neighbor.id);
    } 
  }
  
  Teuchos::ArrayRCP<const scalar_t> bccPropagate(){
    int done = 0;
    std::queue<int> art_queue;
    int bcc_count = 0;
    int* articulation_point_flags = new int[nLocalOwned];
    for(int i = 0; i < nLocalOwned; i++){
      articulation_point_flags[i] = 0;
    }
    //while there are empty vertices
    while(!done){
      
      //see how many articulation points all processors know about
      int globalArtPtCount = 0;
      int localArtPtCount = art_queue.size();
      Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_SUM,1, &localArtPtCount,&globalArtPtCount);
       
      //if none, no one is making progress, so ground two empty neighbors.
      if(globalArtPtCount == 0){
        //search for a pair of empty vertices where one is ghosted on a processor
        int foundGhostPair = 0;
        int ownedVtx = -1, ghostVtx = -1;
        for(int i = 0; i < nLocalOwned; i++){
          if(femvData[i].getGroundingStatus() == ICEPROPGS_NONE){
            lno_t out_degree = g->out_degree(i);//out_degree(g, i);
            lno_t* outs = g->out_vertices(i);//out_vertices(g, i);
            for(int j = 0; j < out_degree; j++){
              if(outs[j] >= nLocalOwned){
                //we're dealing with a ghosted vertex
                if(femvData[outs[j]].getGroundingStatus() == ICEPROPGS_NONE){
                  ownedVtx = i;
                  ghostVtx = outs[j];
                  foundGhostPair = 1;
                  break;
                }
              }
            }
          }
          if(foundGhostPair) break;
        }
        femv->beginFill();
        //if you found a ghost pair, send your own procID, otherwise send -1
        int neighborProc = -1;
        int neighborSend = -1;
        if(foundGhostPair) neighborSend = me;
        Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_MAX, 1, &neighborSend, &neighborProc);
        
        //if neighborProc is me, I have to ground the neighbors.
        if(neighborProc == me){
          //std::cout<<me<<": Grounding empty neighbors,"<<mapWithCopies->getGlobalElement(ownedVtx)<<" and "<<mapWithCopies->getGlobalElement(ghostVtx)<<"\n";
          //replace local value with self-grounded vertex with new bcc_name
          IcePropVtxLabel firstNeighbor = femvData[ownedVtx];
          IcePropVtxLabel secondNeighbor = femvData[ghostVtx];
          gno_t firstNeighbor_gid = mapWithCopies->getGlobalElement(firstNeighbor.id);          
	  gno_t secondNeighbor_gid = mapWithCopies->getGlobalElement(secondNeighbor.id);
          firstNeighbor.first_label = firstNeighbor_gid;
          firstNeighbor.first_sender = firstNeighbor_gid;
          firstNeighbor.bcc_name = bcc_count*np + me;
          secondNeighbor.first_label = secondNeighbor_gid;
          secondNeighbor.first_sender = secondNeighbor_gid;
          secondNeighbor.bcc_name = bcc_count*np + me;
          femv->replaceLocalValue(ownedVtx,0, firstNeighbor);
          femv->replaceLocalValue(ghostVtx,0, secondNeighbor);
          //push the neighbors on the reg queue
          icePropRegQueue.push(ownedVtx);
          icePropRegQueue.push(ghostVtx);
        } else if(neighborProc == -1){
          int foundEmptyPair = 0;
          lno_t vtx1 = -1, vtx2 = -1;
          //if none are found, find any pair of empty vertices. (similar procedure)
          for(int i = 0; i < nLocalOwned; i++){
            if(femvData[i].getGroundingStatus() == ICEPROPGS_NONE){
              lno_t out_degree =g->out_degree(i);// out_degree(g, i);
              lno_t* outs = g->out_vertices(i);//out_vertices(g, i);
              for(int j = 0; j < out_degree; j++){
                if(femvData[outs[j]].getGroundingStatus() == ICEPROPGS_NONE){
                  foundEmptyPair = 1;
                  vtx1 = i;
                  vtx2 = outs[j];
                  break;
                }
              }
            }
            if(foundEmptyPair) break;
          }
          
          int emptyProc = -1;
          int emptySend = -1;
          if(foundEmptyPair) emptySend = me;
          Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_MAX, 1, &emptySend, &emptyProc);
          
          //if emptyProc is -1, no processor has a pair of empty neighboring vertices, so we can't do anything
          
          if(emptyProc == -1){
            //std::cout<<me<<": Couldn't find two empty neighbors\n";
            femv->endFill();
            break;
          }
          else if(emptyProc == me){
            //std::cout<<me<<": Grounding "<<mapWithCopies->getGlobalElement(vtx1)<<" and neighbor "<<mapWithCopies->getGlobalElement(vtx2)<<"\n";
            //this processor will ground two random, empty neighbors.
            IcePropVtxLabel firstNeighbor = femvData[vtx1];
            IcePropVtxLabel secondNeighbor = femvData[vtx2];
            gno_t firstNeighbor_gid = mapWithCopies->getGlobalElement(firstNeighbor.id);
            gno_t secondNeighbor_gid = mapWithCopies->getGlobalElement(secondNeighbor.id);
            firstNeighbor.first_label = firstNeighbor_gid;
            firstNeighbor.first_sender = firstNeighbor_gid;
            firstNeighbor.bcc_name = bcc_count*np + me;
            secondNeighbor.first_label = secondNeighbor_gid;
            secondNeighbor.first_sender = secondNeighbor_gid;
            secondNeighbor.bcc_name = bcc_count*np + me;
            femv->replaceLocalValue(vtx1, 0, firstNeighbor);
            femv->replaceLocalValue(vtx2, 0, secondNeighbor);
            //put the neighbors on the regular queue
            icePropRegQueue.push(vtx1);
            icePropRegQueue.push(vtx2);
          }
        }
        femv->endFill();
      } else {
        femv->beginFill();
        
        //if this processor knows about articulation points
        if(!art_queue.empty()){
          //look at the front, and ground a neighbor.
          lno_t art_pt = art_queue.front();
          lno_t out_degree = g->out_degree(art_pt);//out_degree(g, art_pt);
          lno_t* outs = g->out_vertices(art_pt);//out_vertices(g, art_pt);
          for(int i = 0;i < out_degree; i++){
            if(femvData[outs[i]].getGroundingStatus() == ICEPROPGS_NONE){
              //std::cout<<me<<": Grounding "<<mapWithCopies->getGlobalElement(art_queue.front())<<" and neighbor "<<mapWithCopies->getGlobalElement(outs[i])<<"\n";
              IcePropVtxLabel neighbor = femvData[outs[i]];
              IcePropVtxLabel artvtx = femvData[art_pt];
              gno_t neighbor_gid = mapWithCopies->getGlobalElement(neighbor.id);
              neighbor.first_label = neighbor_gid;
              neighbor.first_sender = neighbor_gid;
              neighbor.bcc_name = bcc_count*np+me;
              artvtx.bcc_name = bcc_count*np+me;
              femv->replaceLocalValue(art_pt,0,artvtx);
              femv->replaceLocalValue(outs[i],0, neighbor);
              icePropRegQueue.push(art_pt);
              icePropRegQueue.push(outs[i]);
              break;
            }
          }
          
        }
        femv->endFill();
      }
      
      //call this->propagate, which, per-processor, finds one bcc at a time.
      //std::cout<<me<<": Calling propagate\n";
      int* t = propagate();
      //printFEMV("AfterPropagation");
      bcc_count++;
      delete [] t;
      //std::cout<<me<<": Checking for articulation points\n";
      //check for OWNED articulation points
      for(int i = 0; i < nLocalOwned; i++){
        if(femvData[i].getGroundingStatus() == ICEPROPGS_FULL){
          lno_t out_degree = g->out_degree(i);//out_degree(g, i);
          lno_t* outs = g->out_vertices(i);//out_vertices(g, i);
          for(int j = 0; j < out_degree; j++){
            if(femvData[outs[j]].getGroundingStatus() < ICEPROPGS_FULL){
              art_queue.push(i);
              articulation_point_flags[i] = 1;
              break;
            }
          }
        }
      }
      //std::cout<<me<<": Clearing half labels\n";
      //clear half labels
      for(int i = 0; i < nLocalOwned+nLocalCopy; i++){
        IcePropVtxLabel label = femvData[i];
        if(label.getGroundingStatus() == ICEPROPGS_HALF){
          label.first_label = -1;
          label.first_sender = -1;
          label.bcc_name = -1;
          femv->switchActiveMultiVector();
          femv->replaceLocalValue(i,0,label);
          femv->switchActiveMultiVector();
        }
      }
      //pop articulation points off of the art_queue if necessary.
      //std::cout<<me<<": Removing fully explored articulation points\n";
      if(!art_queue.empty()){
        bool pop_art = true;
        while(pop_art && !art_queue.empty()){
          //std::cout<<me<<": Checking artpt "<<art_queue.front()<<"\n";
          lno_t top_art_degree = g->out_degree(art_queue.front());//out_degree(g,art_queue.front());
          lno_t* art_outs = g->out_vertices(g,art_queue.front());//out_vertices(g,art_queue.front());
        
          for(int i = 0; i < top_art_degree; i++){
            if(femvData[art_outs[i]].getGroundingStatus() < ICEPROPGS_FULL){
              pop_art = false;
            }
          }
          if(pop_art) art_queue.pop();
        }
      }
      //std::cout<<me<<": Starting over\n";
    }
    //std::cout<<me<<": found "<<bcc_count<<" biconnected components\n";
    for(int i = 0; i < nLocalOwned; i++){
      if(!articulation_point_flags[i]){
        IcePropVtxLabel label = femvData[i];
        label.is_art = false;
        femv->replaceLocalValue(i,0,label);
      }
    }
    return femv->getData(0);
  }
  //int vtxLabelUnitTest();

private:
  int me; 	    //my processor rank
  int np;	    //number of processors
  int nLocalOwned;  //number of vertices owned by this processor
  int nLocalCopy;   //number of copies of off-processor vertices on this proc
  int nVec;	    //number of vectors in multivector
  
  const Teuchos::RCP<const Teuchos::Comm<int> > comm; //MPI communicator
  icePropGraph<lno_t>* g;	    //csr representation of vertices on this processor

  Teuchos::RCP<const MAP> mapOwned;       //Tpetra::Map including only owned
  Teuchos::RCP<const MAP> mapWithCopies;  //Tpetra::Map including owned
                                            //vertices and copies

  Teuchos::RCP<femv_t> femv;
  
  Teuchos::ArrayRCP<const scalar_t> femvData;
};
}//end namespace iceProp

//  Unit test for the vtxLabel struct
//  Make sure vtxLabel's overloaded operators compile and work as expected.
//
//  Tests for coverage of += operator	done
//  	empty non-art += empty non-art    x
//  	empty non-art += half non-art     x
//  	empty non-art += full non-art	  x
//  	empty non-art += empty art	  x
//  	empty non-art += half art         x
//  	empty non-art += full art         x
//  	half non-art += empty non-art     x
//  	half non-art += half non-art   	  x
//	half non-art += full non-art      x
//	half non-art += empty art	  x
//	half non-art += half art	  x
//	half non-art += full art	  x
//	full non-art += empty non-art     x
//	full non-art += half non-art      x
//	full non-art += full non-art	  x
//	full non-art += empty art 	  x
//	full non-art += half art          x
//	full non-art += full art          x
//  
//  Note: suspected articulation points (SAPs?) receive labels in the 
//  same way that regular vertices do, so the only cases we need to test
//  for exercising the articulation point logic are the cases where
//  suspected articulation points are on the right hand side of the +=.
//
/*template<typename MAP>
int  iceProp::iceSheetPropagation<MAP>::vtxLabelUnitTest()
{
  int ierr = 0;
  iceProp::vtxLabel a(0);
  iceProp::vtxLabel a_same(0,-1,-1,-1,-1,false);
  
  //equality test 
  if(!(a == a_same)){
    std::cout<<"UnitTest Error: default constructor arguments not as expected "<<a<<" != "<<a_same<<"\n";
    ierr++;
  }

//  	empty non-art += empty non-art
  iceProp::vtxLabel ena1(0);
  iceProp::vtxLabel ena2(2);
  iceProp::vtxLabel result(0);
  giveLabels(ena2,ena1);
  
  if(!(ena1==result)){
    std::cout<<"UnitTest Error: empty non-art += empty non-art not as expected: "<<ena1<<" != "<<result<<"\n";
    ierr++;
  }

//  	empty non-art += half non-art 
  iceProp::vtxLabel hna1(0,1,1,-1,-1,false);
  ena1 = iceProp::vtxLabel(2);
  
  if(!iceProp::reg.empty()){
    std::cout<<"UnitTest Error: regular queue is not empty, but it should be\n";
    ierr++;
  }

  giveLabels(hna1,ena1);
  result = iceProp::vtxLabel(2,1,0,-1,-1,false);
  
  if(!(ena1 == result)){
    std::cout<<"UnitTest Error: empty non-art += half non-art not as expected "<<ena1<<" != "<<result<<"\n";
    ierr++;
  }

  if(iceProp::reg.front() != 2){
    std::cout<<"UnitTest Error: += did not update the regular queue: size is "<<iceProp::reg.size()<<", front is "<<iceProp::reg.front()<<"\n";
    ierr++;
  }
  
//  	empty non-art += full non-art
  ena1 = iceProp::vtxLabel(2);
  iceProp::vtxLabel fna1(10,3,5,7,9,false);
  result = iceProp::vtxLabel(2,3,10,7,10,false);

  //ena1 += fna1;
  giveLabels(fna1,ena1);
  if(!(ena1 == result)){
    std::cout<<"UnitTest Error: empty non-art += full non-art not as expected: "<<ena1<<" != "<<result<<"\n";
    ierr++;
  }

//  	empty non-art += empty art
  ena1 = iceProp::vtxLabel(2);
  iceProp::vtxLabel ea1(10,-1,-1,-1,-1,true);
  result = iceProp::vtxLabel(2);

  //ena1 += ea1;
  giveLabels(ea1,ena1);
  if(!(result == ena1)){
    std::cout<<"UnitTest Error: empty non-art += empty art not as expected: "<<ena1<<" != "<<result<<"\n";
    ierr++;
  }
  
//  	empty non-art += half art
  ena1 = iceProp::vtxLabel(2);
  iceProp::vtxLabel ha1(10,3,4,-1,-1,true);
  result = iceProp::vtxLabel(2,3,10,-1,-1,false);

  //ena1 += ha1;
  giveLabels(ha1,ena1);
  
  if(!(ena1 == result)){
    std::cout<<"UnitTest Error: empty non-art += half art not as expected: "<<ena1<<" != "<<result<<"\n";
    ierr++;
  }

//  	empty non-art += full art
  ena1 = iceProp::vtxLabel(2);
  iceProp::vtxLabel fa1(10,3,4,7,8,true);
  result = iceProp::vtxLabel(2,10,10,-1,-1,false);
 
  //ena1 += fa1;
  giveLabels(fa1,ena1);
  
  if(!(ena1 == result)){
    std::cout<<"UnitTest Error: empty non-art += full art not as expected: "<<ena1<<" != "<<result<<"\n";
    ierr++;
  }
  
//  	half non-art += empty non-art
  hna1 = iceProp::vtxLabel(2,3,5,-1,-1,false);
  ena1 = iceProp::vtxLabel(10);
  result = iceProp::vtxLabel(2,3,5,-1,-1,false);

  //hna1 += ena1;
  giveLabels(ena1,hna1);
  
  if(!(hna1 == result)){
    std::cout<<"UnitTest Error: half non-art += empty non-art not as expected: "<<hna1<<" != "<<result<<"\n";
    ierr++;
  }

//  	half non-art += half non-art   	  
  //same labels (senders don't matter in this case)
  hna1 = iceProp::vtxLabel(2,3,5,-1,-1,false);
  iceProp::vtxLabel hna2(4,3,9,-1,-1,false);
  result = iceProp::vtxLabel(2,3,5,-1,-1,false);
  
  //hna1+=hna2;
  giveLabels(hna2, hna1);

  if(!(hna1 == result)){
    std::cout<<"UnitTest Error: half non-art += half non-art with same labels not as expected: "<<hna1<<" != "<<result<<"\n";
    ierr++;
  }

  //hna2 has already sent to hna1
  hna1 = iceProp::vtxLabel(2,3,5,-1,-1,false);
  hna2 = iceProp::vtxLabel(5,3,9,-1,-1,false);
  result = iceProp::vtxLabel(2,3,5,-1,-1,false);
  
  //hna1 += hna2;
  giveLabels(hna2,hna1);

  if(!(result == hna1)){
    std::cout<<"UnitTest Error: half non-art += half non-art where rhs already sent to lhs not as expected: "<<hna1<<" != "<<result<<"\n";
    ierr++;
  }
  
  //different labels
  hna1 = iceProp::vtxLabel(2,3,5,-1,-1,false);
  hna2 = iceProp::vtxLabel(4,7,9,-1,-1,false);
  result = iceProp::vtxLabel(2,3,5,7,4,false);

  //hna1 += hna2;
  giveLabels(hna2,hna1);

  if(!(hna1 == result)){
    std::cout<<"UnitTest Error: half non-art += half non-art with different labels not as expected: "<<hna1<<" != "<<result<<"\n";
    ierr++;
  }
    
//	half non-art += full non-art
  //different labels, different senders
  hna1 = iceProp::vtxLabel(2,3,5,-1,-1,false);
  fna1 = iceProp::vtxLabel(6,4,9,1,0,false);
  result = iceProp::vtxLabel(2,4,6,1,6,false);

  //hna1 += fna1;
  giveLabels(fna1,hna1);
  
  if(!(result == hna1)){
    std::cout<<"UnitTest Error: half non-art += full non-art with different labels not as expected: "<<hna1<<" != "<<result<<"\n";
    ierr++;
  }

  //full sent the half's only label
  hna1 = iceProp::vtxLabel(2,3,5,-1,-1,false);
  fna1 = iceProp::vtxLabel(5,3,3,6,0,false);
  result = iceProp::vtxLabel(2,3,5,6,5,false);
  
  //hna1 += fna1;
  giveLabels(fna1,hna1);  

  if(!(result == hna1)){
    std::cout<<"UnitTest Error: half non-art += full non-art where full already sent half one label not as expected: "<<hna1<<" != "<<result<<"\n";
    ierr++;
  } 
  //full and half have one label in common
  hna1 = iceProp::vtxLabel(2,3,5,-1,-1,false);
  fna1 = iceProp::vtxLabel(1,3,6,8,0,false);
  result = iceProp::vtxLabel(2,3,1,8,1,false);
  
  //hna1 += fna1;
  giveLabels(fna1,hna1);  

  if(!(result == hna1)){
    std::cout<<"UnitTest Error: half non-art += full non-art where full & half share one label not as expected: "<<hna1 <<" !=  " <<result<<"\n";
    ierr++;
  }
  
//	half non-art += empty art
  hna1 = iceProp::vtxLabel(2,3,5,-1,-1,false);
  ea1 = iceProp::vtxLabel(29,-1,-1,-1,-1,false);
  result = hna1;
  
  //hna1 += ea1;
  giveLabels(ea1,hna1);

  if(!(result == hna1)){
    std::cout<<"UnitTest Error: half non-art += empty art not as expected: "<<hna1<<" != "<<result<<"\n";
    ierr++;
  }

//	half non-art += half art
  //same label
  hna1 = iceProp::vtxLabel(2,4,5,-1,-1,false);
  ha1 = iceProp::vtxLabel(10,4,3,-1,-1,true);
  result = hna1;

  //hna1 += ha1;
  giveLabels(ha1,hna1);

  if(!(result == hna1)){
    std::cout<<"UnitTest Error: half non-art += half art same labels not as expected: "<<hna1<<" != "<<result<<"\n";
    ierr++;
  }

  //ha1 sent hna1 its only label
  hna1 = iceProp::vtxLabel(2,4,5,-1,-1,false);
  ha1 = iceProp::vtxLabel(5,4,3,-1,-1,true);
  result = hna1;

  //hna1 += ha1;
  giveLabels(ha1,hna1);  

  if(!(result == hna1)){
    std::cout<<"UnitTest Error: half non-art += half art where art sent to non-art already not as expected: "<<hna1<<" != "<<result<<"\n";
    ierr++;
  }

  //different labels
  hna1 = iceProp::vtxLabel(2,3,4,-1,-1,false);
  ha1 = iceProp::vtxLabel(5,6,7,-1,-1,true);
  result = iceProp::vtxLabel(2,3,4,6,5,false);

  //hna1 += ha1;
  giveLabels(ha1,hna1);  

  if(!(result == hna1)){
    std::cout<<"UnitTest Error: half non-art += half art with different labels not as expected: "<<hna1<<" != "<<result<<"\n";
    ierr++;
  }

//	half non-art += full art
  //sharing a label
  hna1 = iceProp::vtxLabel(2,4,5,-1,-1,false);
  fa1 = iceProp::vtxLabel(6,4,3,9,0,true);
  result = iceProp::vtxLabel(2,4,5,6,6,false);

  //hna1 += fa1;
  giveLabels(fa1,hna1); 

  if(!(result == hna1)){
    std::cout<<"UnitTest Error: half non-art += full art sharing a label not as expected: "<<hna1<<" != "<< result<<"\n";
    ierr++;
  }

  //full sent a label to half
  hna1 = iceProp::vtxLabel(2,4,5,-1,-1,false);
  fa1 = iceProp::vtxLabel(5,4,3,6,7,true);
  result = hna1;

  //hna1 += fa1;
  giveLabels(fa1, hna1);

  if(!(hna1 == result)){
    std::cout<<"UnitTest Error: half non-art += full art full already sent a label not as expected: "<<hna1 <<" != "<<result<<"\n";
    ierr++;
  }

  //different labels
  hna1 = iceProp::vtxLabel(2,3,4,-1,-1,false);
  fa1 = iceProp::vtxLabel(1,5,6,7,8,true);
  result = iceProp::vtxLabel(2,3,4,1,1,false);

  //hna1 += fa1;
  giveLabels(fa1,hna1);

  if(!(result == hna1)){
    std::cout<<"UnitTest Error half non-art += full art different labels not as expected: "<<hna1<<" != "<<result<<"\n";
    ierr++;
  }

//	full non-art += empty non-art
  fna1 = iceProp::vtxLabel(2,3,4,5,6,false);
  ena1 = iceProp::vtxLabel(7);
  result = fna1;

  //fna1 += ena1;
  giveLabels(ena1, fna1);

  if(!(result == fna1)){
    std::cout<<"UnitTest Error full non-art += empty non-art not as expected: "<<fna1<<" != "<<result<<"\n";
    ierr++;
  }
  
//	full non-art += half non-art
  fna1 = iceProp::vtxLabel(2,3,4,5,6,false);
  hna1 = iceProp::vtxLabel(9,1,3,-1,-1,false);
  result = fna1;
  
  //fna1 += hna1;
  giveLabels(hna1, fna1);  

  if(!(result == fna1)){
    std::cout<<"UnitTest Error full non-art += half non-art not as expected: "<<fna1<<" != "<<result<<"\n";
    ierr++;
  } 


//	full non-art += full non-art
  fna1 = iceProp::vtxLabel(2,3,4,5,6,false);
  iceProp::vtxLabel fna2(6,7,8,9,10,false);
  result = fna1;

  //fna1 += fna2;
  giveLabels(fna2, fna1);

  if(!(result == fna1)){
    std::cout<<"UnitTest Error full non-art += full non-art not as expected: "<<fna1<<" != "<<result<<"\n";
    ierr++;
  }
 
//	full non-art += empty art
  fna1 = iceProp::vtxLabel(2,4,5,6,7,false);
  ea1 = iceProp::vtxLabel(1,-1,-1,-1,-1,true);
  result = fna1;
  
  //fna1 += ea1;
  giveLabels(ea1, fna1);

  if(!(result == fna1)){
    std::cout<<"UnitTest Error full non-art += empty art not as expected: "<<fna1<<" != "<<result<<"\n";
    ierr++;
  }

//	full non-art += half art
  fna1 = iceProp::vtxLabel(2,3,4,5,6,false);
  ha1 = iceProp::vtxLabel(1,10,11,-1,-1,true);
  result = fna1;

  //fna1 += ha1;
  giveLabels(ha1, fna1);

  if(!(result == fna1)){
    std::cout<<"UnitTest Error full non-art += half art not as expected: "<<fna1<<" != "<<result<<"\n";
    ierr++;
  }

//	full non-art += full art
  fna1 = iceProp::vtxLabel(2,3,4,5,6,false);
  fa1 = iceProp::vtxLabel(10,11,12,13,14,true);
  result = fna1;
  
  //fna1 += fa1;
  giveLabels(fa1,fna1);  

  if(!(result == fna1)) {
    std::cout<<"UnitTest Error full non-art += full art not as expected: "<<fna1<<" != "<<result<<"\n";
    ierr++;
  }
  while(!iceProp::reg.empty())iceProp::reg.pop();
  if(ierr == 0 ){
    std::cout<<"iceSheetPropagation::giveLabels OK\n";
  }
  return ierr;
}*/


#endif

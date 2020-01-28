#include "Tpetra_Core.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_FEMultiVector.hpp"

#include <string>
#include <sstream>
#include <iostream>
#include <queue>

#ifndef ZOLTAN2_VTXLABEL_
#define ZOLTAN2_VTXLABEL_

namespace Zoltan2{
        
	//Graph datastructure that represents the local
	//graph in ice sheet propagation problems
        template<typename lno_t>
	class icePropGraph {
          public:
            //n represents the number of vertices in the local graph
            size_t n;
            //m represents the number of edges in the local graph
            size_t m;
            //out_array is the adjacency array for each vertex
            lno_t* out_array;
            //out_degree_list is the offset array for indexing the out_array
            lno_t* out_degree_list;

            //member function to return the degree of a given vertex.
            lno_t out_degree(lno_t vert){
              return (out_degree_list[vert+1] - out_degree_list[vert]);
            }
            //member function to return the neighbors for a given vertex
            lno_t* out_vertices(lno_t vert){
              return (&out_array[out_degree_list[vert]]);
            }
        }; 
        
        //global queues to allow easy propagation from the
	//+= operator of IcePropVtxLabel.
	std::queue<int> icePropArtQueue;
	std::queue<int> icePropRegQueue;

        //Enum that denotes
	enum IcePropGrounding_Status {ICEPROPGS_FULL=2, ICEPROPGS_HALF=1, ICEPROPGS_NONE = 0};
	// Struct representing a vertex label.
	// We define our own "addition" for these labels.
	// Later, we'll create a Tpetra::FEMultiVector of these labels.
	template<typename lno_t, typename gno_t>
	class IcePropVtxLabel {
	public: 
          //The local ID for the vertex represented by this label
          lno_t id;
          
          //the global ID of a grounded vertex
	  gno_t first_label;
          
          //the global ID of the vertex that sent the label
	  gno_t first_sender;
          
          //this field indicates whether or not the first label is currently set
          //(this is necessary for unsigned global types)
          bool  first_used;
          
          //the global ID of a second grounded vertex
	  gno_t second_label;

          //the global ID of the vertex that sent the second identifier
	  gno_t second_sender;

          //this field indicates whether of not the second label is currently set
          //(this is necessary for unsigned global types)
          bool  second_used;

          //this field may potentially be used to find biconnected components
          int bcc_name;

          //this flag indicates whether or not this vertex is a potential articulation point
	  bool is_art;

	  // Constructors
	  IcePropVtxLabel(int idx_, int first_ = -1, int first_sender_ = -1, bool first_used_ = false, 
                                    int second_ = -1, int second_sender_ = -1, bool second_used_ = false,
                                    bool art_ = false, int bcc_name_ = -1) { 
	    id = idx_;
	    if(id == -1){
		std::cout<<"A label's ID is -1\n";
	    }
	    first_label = first_;
	    first_sender = first_sender_;
            first_used = first_used_;
	    second_label = second_;
	    second_sender = second_sender_; 
            second_used = second_used_;
	    is_art = art_;
            bcc_name = bcc_name_;
	  }
	  IcePropVtxLabel() {
	    id = -1;
	    first_label = -1;
	    first_sender = -1;
            first_used = false;
	    second_label = -1;
	    second_sender = -1;
            second_used = false;
	    is_art = false; 
            bcc_name = -1;
	  }

	  IcePropVtxLabel(volatile const IcePropVtxLabel& other){
	    id = other.id;
	    first_label = other.first_label;
	    first_sender = other.first_sender;
            first_used = other.first_used;
	    second_label = other.second_label;
	    second_sender = other.second_sender;
            second_used = other.second_used;
	    is_art = other.is_art;
            bcc_name = other.bcc_name;
	  }
	  IcePropVtxLabel(const IcePropVtxLabel& other){
	    id = other.id;
	    first_label = other.first_label;
	    first_sender = other.first_sender;
            first_used = other.first_used;
	    second_label = other.second_label;
	    second_sender = other.second_sender;
            second_used = other.second_used;
	    is_art = other.is_art;
            bcc_name = other.bcc_name;
	  }
	  // IcePropVtxLabel assignment
          volatile IcePropVtxLabel operator=(const IcePropVtxLabel& other) volatile{
	    id = other.id;
	    first_label = other.first_label;
	    first_sender = other.first_sender;
            first_used = other.first_used;
	    second_label = other.second_label;
	    second_sender = other.second_sender;
            second_used = other.second_used;
	    is_art = other.is_art;
            bcc_name = other.bcc_name;
	    return *this; 
	  }

          IcePropVtxLabel& operator=(const IcePropVtxLabel& other) {
	    id = other.id;
	    first_label = other.first_label;
	    first_sender = other.first_sender;
            first_used = other.first_used;
	    second_label = other.second_label;
	    second_sender = other.second_sender;
            second_used = other.second_used;
	    is_art = other.is_art;
            bcc_name = other.bcc_name;
	    return *this; 
	  } 

	  // int assignment
	  IcePropVtxLabel& operator=(const int& other) { 
	    first_label = other;
	    first_sender = other;
            first_used = true;
	    second_label = -1;
	    second_sender = -1;
            second_used = false;
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
              first_used = copy.first_used;
              second_label = copy.second_label;
              second_sender = copy.second_sender;
              second_used = copy.second_used;
              bcc_name = copy.bcc_name;
            //handles HALF == HALF
            } else if(owned_gs == copy_gs && owned_gs == ICEPROPGS_HALF){
              if(copy.first_label != first_label){
                second_label = copy.first_label;
                second_sender = copy.first_sender;
                second_used = copy.first_used;
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
	    return (IcePropGrounding_Status)((first_used) + (second_used));
	  }
	};

}//end namespace Zoltan2
        
/////////////////////////////////////////////////////////////////////////
// ArithTraits -- arithmetic traits needed for struct IcePropVtxLabel
// Needed so that Tpetra compiles.
// Not all functions were needed; this is a subset of ArithTraits' traits.
// Modified from kokkos-kernels/src/Kokkos_ArithTraits.hpp's 
// <int> specialization
namespace Kokkos {
  namespace Details {

    template<>
    class ArithTraits<Zoltan2::IcePropVtxLabel<Tpetra::Map<>::local_ordinal_type, Tpetra::Map<>::global_ordinal_type> > {  // specialized for IcePropVtxLabel struct
    public:
      typedef Zoltan2::IcePropVtxLabel<Tpetra::Map<>::local_ordinal_type,Tpetra::Map<>::global_ordinal_type> val_type;
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
struct SerializationTraits<Ordinal, Zoltan2::IcePropVtxLabel<Tpetra::Map<>::local_ordinal_type, Tpetra::Map<>::global_ordinal_type> >:
       public Teuchos::DirectSerializationTraits<Ordinal, Zoltan2::IcePropVtxLabel<Tpetra::Map<>::local_ordinal_type, Tpetra::Map<>::global_ordinal_type> >
{};
}//end namespace Teuchos

namespace Zoltan2{
template<typename MAP>
class iceSheetPropagation {
public:

  //typedefs for the FEMultiVector
  typedef typename MAP::local_ordinal_type lno_t;
  typedef typename MAP::global_ordinal_type gno_t;
  typedef IcePropVtxLabel<lno_t,gno_t> scalar_t;
  typedef IcePropVtxLabel<lno_t,gno_t> label_t;
  typedef Tpetra::FEMultiVector<scalar_t,lno_t, gno_t> femv_t;	
  
  
 	
  //Constructor creates FEMultiVector of IcePropVertexLabels,
  //and initializes the labels according to the input grounding information.
  //It also sets the is_art flag on each label if the vertex is a potential articulation point.
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
      label_t label(i);
      if(boundary_flags[i] > 2 && i < nLocalOwned) {
	label.is_art = true;
      }
      if(grounding_flags[i]){
        label.first_label = gid;
        label.first_sender = gid;
        label.first_used=true;
        icePropRegQueue.push(label.id);
      }
      femv->replaceLocalValue(i,0,label);
    }
 
    femv->endFill(); 
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
  //returns a flag for each vertex:
  // -2 means keep the vertex
  // -1 means the vertex is completely floating
  // any flag > -1 means the vertex is a part of a hinge, and the flag is the ID of the top-level hinge vertex
  const Teuchos::ArrayRCP<const scalar_t>& propagate(void){ 
    //run bfs_prop
    bfs_prop();
    //check for potentially false articulation points
    while(true){
      femv->switchActiveMultiVector(); 
      for(size_t i = 0; i < g->n; i++){
	label_t curr_node = femvData[i];
        if(curr_node.is_art && curr_node.getGroundingStatus() == ICEPROPGS_FULL){
          lno_t out_degree = g->out_degree(curr_node.id);//out_degree(g, curr_node.id);
          lno_t* outs = g->out_vertices(curr_node.id);//out_vertices(g, curr_node.id);
          for(int j = 0; j < out_degree; j++){
            label_t neighbor = femvData[outs[j]];
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
      for(size_t i = 0; i < g->n; i++){
	label_t curr_node = femvData[i];
        if(curr_node.getGroundingStatus() == ICEPROPGS_HALF){
	  label_t cleared(curr_node.id);
          cleared.is_art = curr_node.is_art;
          femv->replaceLocalValue(i,0,cleared);
        }
        if(curr_node.getGroundingStatus() == ICEPROPGS_FULL && curr_node.is_art){
          icePropRegQueue.push(curr_node.id);
        }
      }
      //re-run bfs_prop until incomplete propagation is fixed
      bfs_prop();     
    }
    //check for nodes that are less than full.
    //return flags for each node, -2 for keep, -1 for remove, <vtxID> for singly grounded nodes.
    /*lno_t* removed = new lno_t[g->n];
    for(size_t i = 0; i < g->n; i++){
      IcePropVtxLabel curr_node = femvData[i];
      IcePropGrounding_Status gs = curr_node.getGroundingStatus();
      if(gs == ICEPROPGS_FULL) removed[i] = -2;
      else if(gs == ICEPROPGS_HALF) removed[i] = curr_node.first_label;
      else removed[i] = -1;
    }
    femv->switchActiveMultiVector();
    return removed;*/
    return femvData;
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
      femv->beginFill();
      //visit every node that changed
      std::queue<int>* curr;
      if(icePropRegQueue.empty()) curr = &icePropArtQueue;
      else curr = &icePropRegQueue;
      while(!curr->empty()){
	label_t curr_node = femvData[curr->front()];
        curr->pop();
        
        //if the current node is a copy, it shouldn't propagate out to its neighbors.
        if(curr_node.id >= nLocalOwned) continue;

        lno_t out_degree = g->out_degree(curr_node.id);
        lno_t* outs = g->out_vertices(curr_node.id);
        for(int i = 0; i < out_degree; i++){
	  label_t neighbor = femvData[outs[i]];
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
      femv->endFill();
      femv->doOwnedToOwnedPlusShared(Tpetra::ADD);
      femv->doOwnedPlusSharedToOwned(Tpetra::ADD);
      int local_done = icePropRegQueue.empty() && icePropArtQueue.empty();
      //this call makes sure that if any inter-processor communication changed labels
      //we catch the changes and keep propagating them.
      Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_MIN,1, &local_done,&done);
    }
    
    
  }
  
  //function that exchanges labels between two nodes
  //curr_node gives its labels to neighbor.
  void giveLabels(IcePropVtxLabel<lno_t,gno_t>& curr_node, IcePropVtxLabel<lno_t,gno_t>& neighbor){
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
      neighbor.first_used = curr_node.first_used;
      neighbor.second_label = curr_node.second_label;
      neighbor.second_sender = curr_node_gid;
      neighbor.second_used = curr_node.second_used;
      neighbor.bcc_name = curr_node.bcc_name;
      return;
    } else if (curr_gs == ICEPROPGS_FULL) {
      //if it is an articulation point, and it hasn't sent to this neighbor
      if(neighbor.first_sender != curr_node_gid){
        //send itself as a label
        if(nbor_gs == ICEPROPGS_NONE){
          neighbor.first_label = curr_node_gid;
          neighbor.first_sender = curr_node_gid;
          neighbor.first_used = true;
          neighbor.bcc_name = curr_node.bcc_name;
        } else if(nbor_gs == ICEPROPGS_HALF){
          if(neighbor.first_label != curr_node_gid){
            neighbor.second_label = curr_node_gid;
            neighbor.second_sender = curr_node_gid;
            neighbor.second_used = true;
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
        neighbor.first_used = curr_node.first_used;
        neighbor.bcc_name = curr_node.bcc_name;
      } else if(nbor_gs == ICEPROPGS_HALF){
        //make sure you aren't giving a duplicate label, and that
        //you haven't sent a label to this neighbor before.
        if(neighbor.first_label != curr_node.first_label && neighbor.first_sender != curr_node_gid){
          neighbor.second_label = curr_node.first_label;
          neighbor.second_sender = curr_node_gid;
          neighbor.second_used = true;
        }
      }
    }
    
    if(nbor_gs != neighbor.getGroundingStatus()){
      if(neighbor.is_art) icePropArtQueue.push(neighbor.id);
      else icePropRegQueue.push(neighbor.id);
    } 
  }
  
  //propagation rules to find biconnected components 
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
            lno_t out_degree = g->out_degree(i);
            lno_t* outs = g->out_vertices(i);
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
          //replace local value with self-grounded vertex with new bcc_name
          label_t firstNeighbor = femvData[ownedVtx];
          label_t secondNeighbor = femvData[ghostVtx];
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
              lno_t out_degree =g->out_degree(i);
              lno_t* outs = g->out_vertices(i);
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
            femv->endFill();
            break;
          }
          else if(emptyProc == me){
            //this processor will ground two random, empty neighbors.
            label_t firstNeighbor = femvData[vtx1];
            label_t secondNeighbor = femvData[vtx2];
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
          lno_t out_degree = g->out_degree(art_pt);
          lno_t* outs = g->out_vertices(art_pt);
          for(int i = 0;i < out_degree; i++){
            if(femvData[outs[i]].getGroundingStatus() == ICEPROPGS_NONE){
              label_t neighbor = femvData[outs[i]];
              label_t artvtx = femvData[art_pt];
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
      int* t = propagate();
      bcc_count++;
      delete [] t;
      //check for OWNED articulation points
      for(int i = 0; i < nLocalOwned; i++){
        if(femvData[i].getGroundingStatus() == ICEPROPGS_FULL){
          lno_t out_degree = g->out_degree(i);
          lno_t* outs = g->out_vertices(i);
          for(int j = 0; j < out_degree; j++){
            if(femvData[outs[j]].getGroundingStatus() < ICEPROPGS_FULL){
              art_queue.push(i);
              articulation_point_flags[i] = 1;
              break;
            }
          }
        }
      }
      //clear half labels
      for(int i = 0; i < nLocalOwned+nLocalCopy; i++){
        label_t label = femvData[i];
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
      if(!art_queue.empty()){
        bool pop_art = true;
        while(pop_art && !art_queue.empty()){
          lno_t top_art_degree = g->out_degree(art_queue.front());
          lno_t* art_outs = g->out_vertices(g,art_queue.front());
        
          for(int i = 0; i < top_art_degree; i++){
            if(femvData[art_outs[i]].getGroundingStatus() < ICEPROPGS_FULL){
              pop_art = false;
            }
          }
          if(pop_art) art_queue.pop();
        }
      }
    }
    for(int i = 0; i < nLocalOwned; i++){
      if(!articulation_point_flags[i]){
        label_t label = femvData[i];
        label.is_art = false;
        femv->replaceLocalValue(i,0,label);
      }
    }
    return femv->getData(0);
  }

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
}//end namespace Zoltan2

#endif

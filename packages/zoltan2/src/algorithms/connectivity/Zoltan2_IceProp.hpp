#include "Tpetra_Core.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_FEMultiVector.hpp"
#include "Tpetra_MultiVector_decl.hpp"
#include "Tpetra_MultiVector_def.hpp"
#include "Zoltan2_VtxLabel.hpp"

#include <string>
#include <sstream>
#include <iostream>
#include <queue>
#include <vector>

#ifndef ZOLTAN2_ICEPROP_
#define ZOLTAN2_ICEPROP_

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
  
  std::queue<lno_t> artQueue;
  std::queue<lno_t> regQueue;
 	
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
        regQueue.push(label.id);
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
          lno_t out_degree = g->out_degree(curr_node.id);
          lno_t* outs = g->out_vertices(curr_node.id);
          for(int j = 0; j < out_degree; j++){
            label_t neighbor = femvData[outs[j]];
            if(neighbor.getGroundingStatus() == ICEPROPGS_HALF && neighbor.first_label != mapWithCopies->getGlobalElement(curr_node.id) && neighbor.first_sender == mapWithCopies->getGlobalElement(curr_node.id)){
              regQueue.push(curr_node.id);
            }
          }
        }
      }
      int local_done = regQueue.empty();
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
          regQueue.push(curr_node.id);
        }
      }
      //re-run bfs_prop until incomplete propagation is fixed
      bfs_prop();     
    }
    //return the final labels for each owned vertex
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
      if(regQueue.empty()) curr = &artQueue;
      else curr = &regQueue;
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
            if(neighbor.is_art) artQueue.push(neighbor.id);
            else regQueue.push(neighbor.id);
          }
        }
        if(curr->empty()){
          if(curr == &regQueue) curr = &artQueue;
          else curr = &regQueue;
        }
        
      }
      
      //store old vertex statuses for detecting changes
      std::vector<IcePropGrounding_Status> old_statuses(nLocalOwned,ICEPROPGS_NONE);      
      for(size_t i = 0; i < nLocalOwned; i ++){
        old_statuses[i] = femvData[i].getGroundingStatus();
      }
      
      femv->endFill();
      femv->doOwnedToOwnedPlusShared(Tpetra::ADD);
      femv->doOwnedPlusSharedToOwned(Tpetra::ADD);

      //detect any change in status due to communication
      for(size_t i = 0; i < nLocalOwned; i++){
        if(femvData[i].getGroundingStatus() != old_statuses[i]){
          if(femvData[i].is_art){
            artQueue.push(i);
          } else {
            regQueue.push(i);
          }
        }
      }
      int local_done = regQueue.empty() && artQueue.empty();
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
      if(neighbor.is_art) artQueue.push(neighbor.id);
      else regQueue.push(neighbor.id);
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
          regQueue.push(ownedVtx);
          regQueue.push(ghostVtx);
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
            regQueue.push(vtx1);
            regQueue.push(vtx2);
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
              regQueue.push(art_pt);
              regQueue.push(outs[i]);
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

//std::queue<iceSheetPropagation<Tpetra::Map<>>::lno_t> iceSheetPropagation<Tpetra::Map<>> ::regQueue;
//std::queue<iceSheetPropagation<Tpetra::Map<>>::lno_t> iceSheetPropagation<Tpetra::Map<>> ::artQueue;


}//end namespace Zoltan2

#endif

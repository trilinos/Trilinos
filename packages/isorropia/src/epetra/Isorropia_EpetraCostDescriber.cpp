//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Alan Williams (william@sandia.gov)
                or Erik Boman    (egboman@sandia.gov)

************************************************************************
*/
//@HEADER

#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_Exception.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_BlockMap.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

namespace Isorropia {
namespace Epetra {

CostDescriber::CostDescriber()
  : vertex_weights_(),
    graph_edge_weights_(),
    hg_edge_gids_(NULL),
    hg_edge_weights_(NULL),
    num_hg_edge_weights_(0),
    numGlobalVertexWeights_(0),
    numGlobalGraphEdgeWeights_(0),
    numGlobalHypergraphEdgeWeights_(0),
    paramlist_()
{
}

CostDescriber::~CostDescriber()
{
  free_hg_edge_weights_();
}

void CostDescriber::setParameters(const Teuchos::ParameterList& paramlist)
{
  paramlist_ = paramlist;
}

void CostDescriber::setVertexWeights(Teuchos::RefCountPtr<Epetra_Vector> vwts)
{
  if (vertex_weights_.get() != 0){
    vertex_weights_.release();
  }
  vertex_weights_ = vwts;
}

void
CostDescriber::setGraphEdgeWeights(Teuchos::RefCountPtr<Epetra_CrsMatrix> gewts)
{
  if (graph_edge_weights_.get() != 0){
    graph_edge_weights_.release();
    graph_self_edges_.clear();
  }
  graph_edge_weights_ = gewts;

  if (gewts->NumMyDiagonals() > 0){

    // Save list of self edges - we omit them in the Zoltan query functions

    const Epetra_Map &rowmap = gewts->RowMap();

    Epetra_Vector diag(rowmap);
    
    gewts->ExtractDiagonalCopy(diag);

    int nvals = gewts->NumMyRows();
    double *entry;
    diag.ExtractView(&entry);
    for (int i=0; i<nvals; i++){
      if (entry[i] != 0){
        graph_self_edges_.insert( rowmap.GID(i));
      }
    }
  }
}

void
CostDescriber::setHypergraphEdgeWeights(Teuchos::RefCountPtr<Epetra_Vector> hgewts)
{
  free_hg_edge_weights_();
  const Epetra_BlockMap& map = hgewts->Map();

  int numWeights = map.NumMyElements();

  if (numWeights > 0){
    allocate_hg_edge_weights_(numWeights);
    map.MyGlobalElements(hg_edge_gids_);
    double *v;
    int stride;
    hgewts->ExtractView(&v, &stride);
    for (int i=0; i<numWeights; i++){
      hg_edge_weights_[i] = (float)v[i];
    }
  }
}

void
CostDescriber::setHypergraphEdgeWeights(int numHGedges, int *hgGIDs, float *hgEwgts)
{
  free_hg_edge_weights_();
  if (numHGedges > 0){
    allocate_hg_edge_weights_(numHGedges);
    for (int i=0; i<numHGedges; i++){
      hg_edge_weights_[i] = hgEwgts[i];
      hg_edge_gids_[i] = hgGIDs[i];
    }
  }
}

void
CostDescriber::setHypergraphEdgeWeights(int numHGedges, int *hgGIDs, double *hgEwgts)
{
  free_hg_edge_weights_();
  if (numHGedges > 0){
    allocate_hg_edge_weights_(numHGedges);
    for (int i=0; i<numHGedges; i++){
      hg_edge_weights_[i] = (float)hgEwgts[i];
      hg_edge_gids_[i] = hgGIDs[i];
    }
  }
}

bool CostDescriber::haveVertexWeights() const
{
  const int n = getNumVertices();
  return( n > 0);
}

int CostDescriber::getNumVertices() const
{
  return( vertex_weights_.get()==0 ? 0 : vertex_weights_->MyLength() );
}

int CostDescriber::getVertexWeights(std::map<int, float> &wgtMap) const
{
  double *wgts;

  const Epetra_BlockMap& map = vertex_weights_->Map();
  const int length = map.NumMyElements();

  if (length < 1) return 0;

  int *global_ids = map.MyGlobalElements();
  vertex_weights_->ExtractView(&wgts);

  for(int i=0; i<length; ++i) {
    wgtMap[global_ids[i]] = (float)wgts[i];
  }
  return length;
}

void CostDescriber::getVertexWeights(int numVertices,
                                     int* global_ids,
                                     float* weights) const
{
  if (getNumVertices() == 0){
    return;
  }

  const Epetra_BlockMap& map = vertex_weights_->Map();

  if (numVertices != map.NumMyElements()) {
    throw Isorropia::Exception("CostDescriber::getVertexWeights: wrong numVertices");
  }

  map.MyGlobalElements(global_ids);

  double* vals = vertex_weights_->Values();
  for(int i=0; i<numVertices; ++i) {
    weights[i] = vals[i];
  }
}

bool CostDescriber::haveGraphEdgeWeights() const
{
  int n = 0;
  if (graph_edge_weights_.get()){
    n = graph_edge_weights_->NumMyNonzeros();
  }
  return( n > 0);
}

int CostDescriber::getNumGraphEdges(int vertex_global_id) const
{
  int n = 0;
  if (graph_edge_weights_.get() != 0) {
    int lrid = graph_edge_weights_->LRID(vertex_global_id);
    if (lrid >= 0){   
      n = graph_edge_weights_->NumMyEntries(lrid);

      if (graph_self_edges_.size() > 0){
        std::set<int>::iterator it = graph_self_edges_.find(vertex_global_id);
        if (it != graph_self_edges_.end()){
          n--;     // don't count self edges
        }
      }

    }
  }

  return n;
}

int CostDescriber::getEdges(int vertexGID, int len, int *nborGID, float *weights) const
{
  const Epetra_Map &colmap = graph_edge_weights_->ColMap();
  const Epetra_Map &rowmap = graph_edge_weights_->RowMap();

  int vertexLID = rowmap.LID(vertexGID);
  int numRealEdges = getNumGraphEdges(vertexGID);  //excluding self edges

  if (numRealEdges < 1){
    return 0;
  }

  if (len < numRealEdges ){
    throw Isorropia::Exception("CostDescriber::getEdges: length of allocated arrays");
  }

  int self_edge = 0;
  std::set<int>::iterator it = graph_self_edges_.find(vertexGID);
  if (it != graph_self_edges_.end()){
    self_edge = 1;
  }

  int *viewIds;
  double *viewWgts;
  int numedges;         // including self edges

  int rc = graph_edge_weights_->ExtractMyRowView(vertexLID, numedges,
                                                 viewWgts, viewIds);

  if (rc){
    throw Isorropia::Exception("CostDescriber::getEdges: Extract matrix row view");
  }

  if (numedges != (numRealEdges + self_edge)){
    throw Isorropia::Exception("CostDescriber::getEdges: Extract matrix count");
  }

  int nextID = 0;

  for (int j=0; j < numedges; j++){
    int gid = colmap.GID(viewIds[j]);
    if (gid == vertexGID) continue;   // skip the self edges

    nborGID[nextID] = gid;
    weights[nextID] = (float)viewWgts[j];

    nextID++;
  }

  return nextID;
}

int CostDescriber::getGraphEdgeWeights(int vertex_global_id, std::map<int, float> &wgtMap) const
{
  int rowlen = getNumGraphEdges(vertex_global_id);

  if (rowlen < 1){
    return 0;
  }

  float *wgt = new float [rowlen];
  int *nborGID = new int [rowlen];

  int numEdges = getEdges(vertex_global_id, rowlen, nborGID, wgt);

  for (int i=0; i<numEdges; i++){
    wgtMap[nborGID[i]] = wgt[i];
  }

  if (rowlen > 0){
    delete [] nborGID;
    delete [] wgt;
  }
  return numEdges;
}

void CostDescriber::getGraphEdgeWeights(int vertex_global_id,
                                       int num_neighbors,
                                       int* neighbor_global_ids,
                                       float* weights) const
{
  int rowlen = getNumGraphEdges(vertex_global_id);

  if (rowlen < 1){
    return;
  }

  if (rowlen > num_neighbors) {
    throw Isorropia::Exception("CostDescriber::getGraphEdgeWeights: wrong num_neighbors");
  }

  getEdges(vertex_global_id, num_neighbors, neighbor_global_ids, weights);
}

bool CostDescriber::haveHypergraphEdgeWeights() const
{
  return(num_hg_edge_weights_ > 0);
}

int CostDescriber::getNumHypergraphEdgeWeights() const
{
  return num_hg_edge_weights_;
}

void CostDescriber::getHypergraphEdgeWeights(int numEdges,
					     int* global_ids,
					     float* weights) const
{
  if (numEdges != num_hg_edge_weights_) {
    throw Isorropia::Exception("CostDescriber::getHypergraphEdgeWeights: wrong numEdges");
  }

  for(int i=0; i<numEdges; ++i) {
    weights[i] = hg_edge_weights_[i];
    global_ids[i] = hg_edge_gids_[i];
  }
}

bool CostDescriber::haveGlobalVertexWeights() const
{
  return (numGlobalVertexWeights_ > 0);
}
void CostDescriber::setNumGlobalVertexWeights(int num)
{
  numGlobalVertexWeights_ = num;
}
bool CostDescriber::haveGlobalGraphEdgeWeights() const
{
  return (numGlobalGraphEdgeWeights_ > 0);
}
void CostDescriber::setNumGlobalGraphEdgeWeights(int num)
{
  numGlobalGraphEdgeWeights_ = num;
}
bool CostDescriber::haveGlobalHypergraphEdgeWeights() const
{
  return (numGlobalHypergraphEdgeWeights_ > 0);
}
void CostDescriber::setNumGlobalHypergraphEdgeWeights(int num)
{
  numGlobalHypergraphEdgeWeights_ = num;
}



void CostDescriber::allocate_hg_edge_weights_(int n)
{
  free_hg_edge_weights_();
  if (n > 0){
    hg_edge_gids_ = new int [n];
    hg_edge_weights_ = new float [n];
    num_hg_edge_weights_ = n;
  }
}
void CostDescriber::free_hg_edge_weights_()
{
  if (hg_edge_gids_){
    delete [] hg_edge_gids_;
    delete [] hg_edge_weights_;
    hg_edge_gids_ = NULL;
    hg_edge_weights_ = NULL;
    num_hg_edge_weights_ = 0;
  }
}


}//namespace Epetra
}//namespace Isorropia
#endif


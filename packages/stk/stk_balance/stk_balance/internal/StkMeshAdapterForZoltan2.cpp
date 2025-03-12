#include "StkMeshAdapterForZoltan2.hpp"
#include "Zoltan2_MeshAdapter.hpp"      // for MeshEntityType, etc
#include <stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp> // for LocalIdMapper
#include <stk_balance/balanceUtils.hpp>

StkMeshZoltanAdapter::StkMeshZoltanAdapter(const Zoltan2ParallelGraph &graph)
  : mGraph(graph)
{ }

size_t StkMeshZoltanAdapter::getLocalNumOf(Zoltan2::MeshEntityType etype) const
{
  size_t localNum = 0;
  if(etype == Zoltan2::MESH_REGION)
  {
    localNum = mGraph.num_vertices();
  }
  return localNum;
}

size_t StkMeshZoltanAdapter::getGlobalNumOf(Zoltan2::MeshEntityType etype) const
{
  size_t globalNum = 0;
  if(etype == Zoltan2::MESH_REGION)
  {
    globalNum = mGraph.get_num_global_elements();
  }
  return globalNum;
}

void StkMeshZoltanAdapter::getIDsViewOf(Zoltan2::MeshEntityType etype, BalanceGlobalNumber const *&Ids) const
{
  Ids = NULL;
  if(etype == Zoltan2::MESH_REGION)
  {
    if(mGraph.num_vertices() != 0)
    {
      Ids = mGraph.get_vertex_ids().data();
    }
  }
}

int StkMeshZoltanAdapter::getDimension() const
{
  return mGraph.get_spatial_dim();
}

void StkMeshZoltanAdapter::debuggingInfo(int proc_id, std::ofstream& out) const
{
  int stride = 0;
  const scalar_t** coords = new const scalar_t*[this->getDimension()];

  for(int coordim = 0; coordim < this->getDimension(); ++coordim)
  {
    this->getCoordinatesViewOf(getPrimaryEntityType(), coords[coordim], stride, coordim);
  }

  // vertex weights
  const scalar_t** weights = new const scalar_t*[getNumWeightsPerOf(getPrimaryEntityType())];
  int weightstride = 0;
  for(int weightdim = 0; weightdim < getNumWeightsPerOf(getPrimaryEntityType()); ++weightdim)
  {
    this->getWeightsViewOf(getPrimaryEntityType(), weights[weightdim], weightstride, weightdim);
  }

  std::ostringstream os;
  os << "Output of coordinates on processor " << proc_id << std::endl;
  for(size_t i=0;i<getLocalNumOf(getPrimaryEntityType());++i)
  {
    os << "[" << proc_id << "] Vertex " << i << " has coordinates: ";
    for(int coordim = 0; coordim < this->getDimension(); ++coordim)
    {
      os << coords[coordim][stride*i] << "\t";
    }

    os << " with weights: ";
    for(int weightdim = 0; weightdim < getNumWeightsPerOf(getPrimaryEntityType()); ++weightdim)
    {
      os << weights[weightdim][weightstride*i] << "\t";
    }
    os << std::endl;
  }

  out << os.str();

  delete [] weights;
  delete [] coords;
}

void StkMeshZoltanAdapter::getCoordinatesViewOf(Zoltan2::MeshEntityType etype, const scalar_t *&coords, int &stride, int coordDim) const
{
  coords = NULL;
  stride = getDimension();
  if(etype == Zoltan2::MESH_REGION)
  {
    if(!mGraph.get_vertex_coords().empty())
    {
      coords = mGraph.get_vertex_coords().data() + coordDim;
    }
  }
}

int StkMeshZoltanAdapter::getNumWeightsPerOf(Zoltan2::MeshEntityType etype) const
{
  int numWeightsPerVertex = 0;
  if(etype == Zoltan2::MESH_REGION)
  {
    numWeightsPerVertex = mGraph.get_num_field_criteria();
  }
  return numWeightsPerVertex;
}

void StkMeshZoltanAdapter::getWeightsViewOf(Zoltan2::MeshEntityType etype, const scalar_t *&weights, int &stride, int idx) const
{
  weights = NULL;
  stride = mGraph.get_num_field_criteria();
  if(etype == Zoltan2::MESH_REGION)
  {
    if(!mGraph.get_vertex_weights().empty())
    {
      weights = mGraph.get_vertex_weights().data() + idx;
    }
  }
}

bool StkMeshZoltanAdapter::avail2ndAdjs(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through) const
{
  bool isAvailable = false;
  if(sourcetarget == Zoltan2::MESH_REGION && through == Zoltan2::MESH_FACE)
  {
    isAvailable = true;
  }
  return isAvailable;
}

size_t StkMeshZoltanAdapter::getLocalNum2ndAdjs(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through) const
{
  size_t numAdjacencies = 0u;
  if(sourcetarget == Zoltan2::MESH_REGION && through == Zoltan2::MESH_FACE)
  {
    numAdjacencies = mGraph.get_adjacency().size();
  }
  return numAdjacencies;
}

void StkMeshZoltanAdapter::get2ndAdjsView(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through, const BalanceLocalNumber *&offsets, const BalanceGlobalNumber *&adjacencyIds) const
{
  offsets = NULL;
  adjacencyIds = NULL;
  if(sourcetarget == Zoltan2::MESH_REGION && through == Zoltan2::MESH_FACE)
  {
    offsets = mGraph.get_offsets().data();
    if(!mGraph.get_adjacency().empty())
    {
      adjacencyIds = mGraph.get_adjacency().data();
    }
  }
}

int StkMeshZoltanAdapter::getNumWeightsPer2ndAdj(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through) const
{
  size_t numWeightsPerEdge = 0u;
  if(sourcetarget == Zoltan2::MESH_REGION && through == Zoltan2::MESH_FACE)
  {
    numWeightsPerEdge = 1;
  }
  return numWeightsPerEdge;
}

void StkMeshZoltanAdapter::get2ndAdjWeightsView(Zoltan2::MeshEntityType sourcetarget,
                                                Zoltan2::MeshEntityType through,
                                                const scalar_t *&weights,
                                                int &stride, int /*idx*/) const
{
  weights = NULL;
  stride = 1;
  if(sourcetarget == Zoltan2::MESH_REGION && through == Zoltan2::MESH_FACE)
  {
    if(!mGraph.get_edge_weights().empty())
    {
      weights = mGraph.get_edge_weights().data();
    }
  }
}



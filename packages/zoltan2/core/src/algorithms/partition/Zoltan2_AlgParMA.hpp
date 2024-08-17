// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_ALGPARMA_HPP_
#define _ZOLTAN2_ALGPARMA_HPP_

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Util.hpp>
#include <Zoltan2_TPLTraits.hpp>

//////////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgParMA.hpp
//! \brief interface to the ParMA library
//  
//  This design creates an apf mesh to run the ParMA algorithms on. The
//  final solution is determined by changes from beginning to end of the mesh.
//  This approach allows development closer to that of PUMI setup but at the 
//  cost of creating an extra mesh representation.
// 
//  Available ParMA algorithms are given by setting the parma_method parameter 
//  of the sublist parma_paramaters to one of the following:
//  Vertex       - Balances targeting vertex imbalance
//  Element      - Balances targeting element imbalance
//  VtxElm       - Balances targeting vertex and element imbalance
//  VtxEdgeElm   - Balances targeting vertex, edge, and element imbalance
//  Ghost        - Balances using ghost element aware diffusion      
//  Shape        - Optimizes shape of parts by increasing the size of small part boundaries
//  Centroid     - Balances using centroid diffusion
//////////////////////////////////////////////////////////////////////////////

#ifndef HAVE_ZOLTAN2_PARMA

// Error handling for when ParMA is requested
// but Zoltan2 not built with ParMA.

namespace Zoltan2 {
template <typename Adapter>
class AlgParMA : public Algorithm<Adapter>
{
public:
  typedef typename Adapter::user_t user_t;

  AlgParMA(const RCP<const Environment> &/* env */,
           const RCP<const Comm<int> > &/* problemComm */,
           const RCP<const BaseAdapter<user_t> > &/* adapter */)
  {
    throw std::runtime_error(
          "BUILD ERROR:  ParMA requested but not compiled into Zoltan2.\n"
          "Please set CMake flag Zoltan2_ENABLE_ParMA:BOOL=ON.");
  }
};
}

#endif

#ifdef HAVE_ZOLTAN2_PARMA


#include <apf.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <PCU.h>
#include <parma.h>
#include <apfConvert.h>
#include <apfShape.h>
#include <map>
#include <cassert>

namespace Zoltan2 {

template <typename Adapter>
class AlgParMA : public Algorithm<Adapter>
{

private:

  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::offset_t offset_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::userCoord_t userCoord_t;

  const RCP<const Environment> env;
  const RCP<const Comm<int> > problemComm;
  const RCP<const MeshAdapter<user_t> > adapter;
  
  apf::Mesh2* m;
  apf::Numbering* gids;
  apf::Numbering* origin_part_ids;
  std::map<gno_t, lno_t> mapping_elm_gids_index;

  MPI_Comm mpicomm;
  bool pcu_outside;

  void setMPIComm(const RCP<const Comm<int> > &problemComm__) {
#   ifdef HAVE_ZOLTAN2_MPI
      mpicomm = Teuchos::getRawMpiComm(*problemComm__);
#   else
      mpicomm = MPI_COMM_WORLD;  // taken from siMPI
#   endif
  }
  // provides conversion from an APF entity dimension to a Zoltan2 entity type
  enum MeshEntityType entityAPFtoZ2(int dimension) const {return static_cast<MeshEntityType>(dimension);}

  //provides a conversion from the Zoltan2 topology type to and APF type
  //  throws an error on topology types not supported by APF
  enum apf::Mesh::Type topologyZ2toAPF(enum EntityTopologyType ttype) const {
    if (ttype==POINT)
      return apf::Mesh::VERTEX;
    else if (ttype==LINE_SEGMENT)
      return apf::Mesh::EDGE;
    else if (ttype==TRIANGLE)
      return apf::Mesh::TRIANGLE;
    else if (ttype==QUADRILATERAL)
      return apf::Mesh::QUAD;
    else if (ttype==TETRAHEDRON) 
      return apf::Mesh::TET;
    else if (ttype==HEXAHEDRON)
      return apf::Mesh::HEX;
    else if (ttype==PRISM)
      return apf::Mesh::PRISM;
    else if (ttype==PYRAMID)
      return apf::Mesh::PYRAMID;
    else 
      throw std::runtime_error("APF does not support this topology type");
    
  }

  //Sets the weights of each entity in dimension 'dim' to those provided by the mesh adapter
  //sets all weights in the mesh adapter but currently only one is considered by ParMA
  void setEntWeights(int dim, apf::MeshTag* tag) {
    MeshEntityType etype = entityAPFtoZ2(dim);
    for (int i=0;i<m->getTagSize(tag);i++) {
      apf::MeshIterator* itr = m->begin(dim);
      apf::MeshEntity* ent;
      const scalar_t* ws=NULL;
      int stride;
      if (i<adapter->getNumWeightsPerOf(etype)) 
        adapter->getWeightsViewOf(etype,ws,stride,i);
      int j=0;
      while ((ent= m->iterate(itr)))  {
        double w = 1.0;
        if (ws!=NULL)
          w = static_cast<double>(ws[j]);
        m->setDoubleTag(ent,tag,&w);
        j++;
      }
      m->end(itr);
    }
  }
  
  //Helper function to set the weights of each dimension needed by the specific parma algorithm
  apf::MeshTag* setWeights(bool vtx, bool edge, bool face, bool elm) {
    int num_ws=1;
    if (vtx)
      num_ws = std::max(num_ws,adapter->getNumWeightsPerOf(MESH_VERTEX));
    if (edge)
      num_ws = std::max(num_ws,adapter->getNumWeightsPerOf(MESH_EDGE));
    if (face)
      num_ws = std::max(num_ws,adapter->getNumWeightsPerOf(MESH_FACE));
    if (elm) 
      num_ws = std::max(num_ws,adapter->getNumWeightsPerOf(entityAPFtoZ2(m->getDimension())));
    apf::MeshTag* tag = m->createDoubleTag("parma_weight",num_ws);
    if (vtx)
      setEntWeights(0,tag);
    if (edge)
      setEntWeights(1,tag);
    if (face)
      setEntWeights(2,tag);
    if (elm) {
      setEntWeights(m->getDimension(),tag);
    }
    return tag;
  }


  //APF Mesh construction helper functions modified and placed here to support arbitrary entity types
  void constructElements(const gno_t* conn, lno_t nelem, const offset_t* offsets,
                         const EntityTopologyType* tops, apf::GlobalToVert& globalToVert)
  {
    apf::ModelEntity* interior = m->findModelEntity(m->getDimension(), 0);
    for (lno_t i = 0; i < nelem; ++i) {
      apf::Mesh::Type etype = topologyZ2toAPF(tops[i]);
      apf::Downward verts;
      for (offset_t j = offsets[i]; j < offsets[i+1]; ++j)
	verts[j-offsets[i]] = globalToVert[conn[j]];
      buildElement(m, interior, etype, verts);
    }
  }
  int getMax(const apf::GlobalToVert& globalToVert)
  {
    int max = -1;
    APF_CONST_ITERATE(apf::GlobalToVert, globalToVert, it)
      max = std::max(max, it->first);
    PCU_Max_Ints(&max, 1); // this is type-dependent
    return max;
  }
  void constructResidence(apf::GlobalToVert& globalToVert)
  {
    int max = getMax(globalToVert);
    int total = max + 1;
    int peers = PCU_Comm_Peers();
    int quotient = total / peers;
    int remainder = total % peers;
    int mySize = quotient;
    int self = PCU_Comm_Self();
    if (self == (peers - 1))
      mySize += remainder;
    typedef std::vector< std::vector<int> > TmpParts;
    TmpParts tmpParts(mySize);
    /* if we have a vertex, send its global id to the
       broker for that global id */
    PCU_Comm_Begin();
    APF_ITERATE(apf::GlobalToVert, globalToVert, it) {
      int gid = it->first;
      int to = std::min(peers - 1, gid / quotient);
      PCU_COMM_PACK(to, gid);
    }
    PCU_Comm_Send();
    int myOffset = self * quotient;
    /* brokers store all the part ids that sent messages
       for each global id */
    while (PCU_Comm_Receive()) {
      int gid;
      PCU_COMM_UNPACK(gid);
      int from = PCU_Comm_Sender();
      tmpParts.at(gid - myOffset).push_back(from);
    }
    /* for each global id, send all associated part ids
       to all associated parts */
    PCU_Comm_Begin();
    for (int i = 0; i < mySize; ++i) {
      std::vector<int>& parts = tmpParts[i];
      for (size_t j = 0; j < parts.size(); ++j) {
	int to = parts[j];
	int gid = i + myOffset;
	int nparts = parts.size();
	PCU_COMM_PACK(to, gid);
	PCU_COMM_PACK(to, nparts);
	for (size_t k = 0; k < parts.size(); ++k)
	  PCU_COMM_PACK(to, parts[k]);
      }
    }
    PCU_Comm_Send();
    /* receiving a global id and associated parts,
     lookup the vertex and classify it on the partition
     model entity for that set of parts */
    while (PCU_Comm_Receive()) {
      int gid;
      PCU_COMM_UNPACK(gid);
      int nparts;
      PCU_COMM_UNPACK(nparts);
      apf::Parts residence;
      for (int i = 0; i < nparts; ++i) {
	int part;
	PCU_COMM_UNPACK(part);
	residence.insert(part);
      }
      apf::MeshEntity* vert = globalToVert[gid];
      m->setResidence(vert, residence);
    }
  }

  /* given correct residence from the above algorithm,
   negotiate remote copies by exchanging (gid,pointer)
   pairs with parts in the residence of the vertex */
  void constructRemotes(apf::GlobalToVert& globalToVert)
  {
    int self = PCU_Comm_Self();
    PCU_Comm_Begin();
    APF_ITERATE(apf::GlobalToVert, globalToVert, it) {
      int gid = it->first;
      apf::MeshEntity* vert = it->second;
      apf::Parts residence;
      m->getResidence(vert, residence);
      APF_ITERATE(apf::Parts, residence, rit)
	if (*rit != self) {
	  PCU_COMM_PACK(*rit, gid);
	  PCU_COMM_PACK(*rit, vert);
	}
    }
    PCU_Comm_Send();
    while (PCU_Comm_Receive()) {
      int gid;
      PCU_COMM_UNPACK(gid);
      apf::MeshEntity* remote;
      PCU_COMM_UNPACK(remote);
      int from = PCU_Comm_Sender();
      apf::MeshEntity* vert = globalToVert[gid];
      m->addRemote(vert, from, remote);
    }
  }

public:

  /*! ParMA constructor
   *  \param env  parameters for the problem and library configuration
   *  \param problemComm  the communicator for the problem
   *  \param adapter  the user's input adapter (MeshAdapter)
   */
  AlgParMA(const RCP<const Environment> &env__,
            const RCP<const Comm<int> > &problemComm__,
            const RCP<const IdentifierAdapter<user_t> > &adapter__)
  { 
    throw std::runtime_error("ParMA needs a MeshAdapter but you haven't given it one");
  }

  AlgParMA(const RCP<const Environment> &env__,
            const RCP<const Comm<int> > &problemComm__,
            const RCP<const VectorAdapter<user_t> > &adapter__) 
  { 
    throw std::runtime_error("ParMA needs a MeshAdapter but you haven't given it one");
  }

  AlgParMA(const RCP<const Environment> &env__,
            const RCP<const Comm<int> > &problemComm__,
            const RCP<const GraphAdapter<user_t,userCoord_t> > &adapter__) 
  { 
    throw std::runtime_error("ParMA needs a MeshAdapter but you haven't given it one");
  }

  AlgParMA(const RCP<const Environment> &env__,
            const RCP<const Comm<int> > &problemComm__,
            const RCP<const MatrixAdapter<user_t,userCoord_t> > &adapter__)
  { 
    throw std::runtime_error("ParMA needs a MeshAdapter but you haven't given it one");
    
  }

  AlgParMA(const RCP<const Environment> &env__,
            const RCP<const Comm<int> > &problemComm__,
            const RCP<const MeshAdapter<user_t> > &adapter__) :
    env(env__), problemComm(problemComm__), adapter(adapter__)
  { 
    setMPIComm(problemComm__);
   
    //Setup PCU communications
    //If PCU was already initialized outside (EX: for the APFMeshAdapter) 
    // we don't initialize it again.
    pcu_outside=false;
    if (!PCU_Comm_Initialized())
      PCU_Comm_Init();
    else
      pcu_outside=true;
    PCU_Switch_Comm(mpicomm);

    //Find the mesh dimension based on if there are any regions or faces in the part
    // an all reduce is needed in case one part is empty (Ex: after hypergraph partitioning)
    int dim;
    if (adapter->getLocalNumOf(MESH_REGION)>0)
      dim=3;
    else if (adapter->getLocalNumOf(MESH_FACE)>0)
      dim=2;
    else
      dim=0;
    PCU_Max_Ints(&dim,1);
    if (dim<2)
      throw std::runtime_error("ParMA neeeds faces or region information");
    
    //GFD Currently not allowing ParMA to balance non element primary types
    if (dim!=adapter->getPrimaryEntityType())
      throw std::runtime_error("ParMA only supports balancing primary type==mesh element");
    
    //Create empty apf mesh
    gmi_register_null();
    gmi_model* g = gmi_load(".null");
    enum MeshEntityType primary_type = entityAPFtoZ2(dim);
    m = apf::makeEmptyMdsMesh(g,dim,false);

    //Get entity topology types
    const EntityTopologyType* tops;
    try {
      adapter->getTopologyViewOf(primary_type,tops);
    }
    Z2_FORWARD_EXCEPTIONS
    
    //Get element global ids and part ids
    const gno_t* element_gids;
    const part_t* part_ids;
    adapter->getIDsViewOf(primary_type,element_gids);
    adapter->getPartsView(part_ids);
    for (size_t i =0;i<adapter->getLocalNumOf(primary_type);i++)
      mapping_elm_gids_index[element_gids[i]] = i;
    
    //get vertex global ids
    const gno_t* vertex_gids;
    adapter->getIDsViewOf(MESH_VERTEX,vertex_gids);
    
    //Get vertex coordinates
    int c_dim = adapter->getDimension();
    const scalar_t ** vertex_coords = new const scalar_t*[c_dim];
    int* strides = new int[c_dim];
    for (int i=0;i<c_dim;i++)
      adapter->getCoordinatesViewOf(MESH_VERTEX,vertex_coords[i],strides[i],i);

    //Get first adjacencies from elements to vertices
    if (!adapter->availAdjs(primary_type,MESH_VERTEX))
      throw "APF needs adjacency information from elements to vertices";
    const offset_t* offsets;
    const gno_t* adjacent_vertex_gids;
    adapter->getAdjsView(primary_type, MESH_VERTEX,offsets,adjacent_vertex_gids);
    
    //build the apf mesh
    apf::GlobalToVert vertex_mapping;
    apf::ModelEntity* interior = m->findModelEntity(m->getDimension(), 0);
    for (size_t i=0;i<adapter->getLocalNumOf(MESH_VERTEX);i++) {
      apf::MeshEntity* vtx = m->createVert_(interior);
      scalar_t temp_coords[3];
      for (int k=0;k<c_dim&&k<3;k++) 
	temp_coords[k] = vertex_coords[k][i*strides[k]];

      for (int k=c_dim;k<3;k++)
	temp_coords[k] = 0;  
      apf::Vector3 point(temp_coords[0],temp_coords[1],temp_coords[2]);    
      m->setPoint(vtx,0,point);
      vertex_mapping[vertex_gids[i]] = vtx;
    }
    //Call modified helper functions to build the mesh from element to vertex adjacency
    constructElements(adjacent_vertex_gids, adapter->getLocalNumOf(primary_type), offsets, tops, vertex_mapping);
    constructResidence(vertex_mapping);
    constructRemotes(vertex_mapping);
    stitchMesh(m);
    m->acceptChanges();
    
    
    //Setup numberings of global ids and original part ids
    // for use after ParMA is run
    apf::FieldShape* s = apf::getConstant(dim);
    gids = apf::createNumbering(m,"global_ids",s,1);
    origin_part_ids = apf::createNumbering(m,"origin",s,1);

    //number the global ids and original part ids
    apf::MeshIterator* itr = m->begin(dim);
    apf::MeshEntity* ent;
    int i = 0;
    while ((ent = m->iterate(itr))) {
      apf::number(gids,ent,0,0,element_gids[i]);
      apf::number(origin_part_ids,ent,0,0,PCU_Comm_Self());
      i++; 
    }
    m->end(itr);

    //final setup for apf mesh
    apf::alignMdsRemotes(m);
    apf::deriveMdsModel(m);
    m->acceptChanges();
    m->verify();

    //cleanup temp storage
    delete [] vertex_coords;
    delete [] strides;
  }
  void partition(const RCP<PartitioningSolution<Adapter> > &solution);
  
};

/////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void AlgParMA<Adapter>::partition(
  const RCP<PartitioningSolution<Adapter> > &solution
)
{
  //Get parameters
  std::string alg_name = "VtxElm";
  double imbalance = 1.1;
  double step = .5;
  int ghost_layers=3;
  int ghost_bridge=m->getDimension()-1;
 
  //Get the parameters for ParMA
  const Teuchos::ParameterList &pl = env->getParameters();
  try {
    const Teuchos::ParameterList &ppl = pl.sublist("parma_parameters");
    for (ParameterList::ConstIterator iter =  ppl.begin();
	 iter != ppl.end(); iter++) {
      const std::string &zname = pl.name(iter);
      if (zname == "parma_method") {
	std::string zval = pl.entry(iter).getValue(&zval);
	alg_name = zval;
      }
      else if (zname == "step_size") {
	double zval = pl.entry(iter).getValue(&zval);
	step = zval;
      }
      else if (zname=="ghost_layers" || zname=="ghost_bridge") {
	int zval = pl.entry(iter).getValue(&zval);
	if (zname=="ghost_layers")
	  ghost_layers = zval;
	else
	  ghost_bridge = zval;
      }
    }
  }
  catch (std::exception &e) {
    //No parma_parameters sublist found
  }

  const Teuchos::ParameterEntry *pe2 = pl.getEntryPtr("imbalance_tolerance");
  if (pe2){
    imbalance = pe2->getValue<double>(&imbalance);
  }

  //booleans for which dimensions need weights
  bool weightVertex,weightEdge,weightFace,weightElement;
  weightVertex=weightEdge=weightFace=weightElement=false;

  //Build the selected balancer
  apf::Balancer* balancer;
  const int verbose = 1;
  if (alg_name=="Vertex") {
    balancer = Parma_MakeVtxBalancer(m, step, verbose);
    weightVertex = true;
  }
  else if (alg_name=="Element") {
    balancer = Parma_MakeElmBalancer(m, step, verbose);
    weightElement=true;
  }
  else if (alg_name=="VtxElm") {
    balancer = Parma_MakeVtxElmBalancer(m,step,verbose);
    weightVertex = weightElement=true;
  }
  else if (alg_name=="VtxEdgeElm") {
    balancer = Parma_MakeVtxEdgeElmBalancer(m, step, verbose);
    weightVertex=weightEdge=weightElement=true;    
  }
  else if (alg_name=="Ghost") {
    balancer = Parma_MakeGhostDiffuser(m, ghost_layers, ghost_bridge, step, verbose);
    weightVertex=weightEdge=weightFace=true;
    if (3 == m->getDimension()) {
      weightElement=true;
    }
  }
  else if (alg_name=="Shape") {
    balancer = Parma_MakeShapeOptimizer(m,step,verbose);
    weightElement=true;
  }
  else if (alg_name=="Centroid") {
    balancer = Parma_MakeCentroidDiffuser(m,step,verbose);
    weightElement=true;
  }
  else  {
    throw std::runtime_error("No such parma method defined");
  }
  
  //build the weights
  apf::MeshTag* weights = setWeights(weightVertex,weightEdge,weightFace,weightElement);

  //balance the apf mesh
  balancer->balance(weights, imbalance);
  delete balancer;

  // Load answer into the solution.
  int num_local = adapter->getLocalNumOf(adapter->getPrimaryEntityType());
  ArrayRCP<part_t> partList(new part_t[num_local], 0, num_local, true);

  //Setup for communication
  PCU_Comm_Begin();
  apf::MeshEntity* ent;
  apf::MeshIterator* itr = m->begin(m->getDimension());
  //Pack information back to each elements original owner
  while ((ent=m->iterate(itr))) {
    if (m->isOwned(ent)) {
      part_t target_part_id = apf::getNumber(origin_part_ids,ent,0,0);
      gno_t element_gid = apf::getNumber(gids,ent,0,0);
      PCU_COMM_PACK(target_part_id,element_gid);
    }
  }
  m->end(itr);

  //Send information off
  PCU_Comm_Send();

  //Unpack information and set new part ids
  while (PCU_Comm_Receive()) {
    gno_t global_id;
    PCU_COMM_UNPACK(global_id);
    lno_t local_id = mapping_elm_gids_index[global_id];
    part_t new_part_id = PCU_Comm_Sender();
    partList[local_id] = new_part_id;
  }
  //construct partition solution
  solution->setParts(partList);
  
  // Clean up
  apf::destroyNumbering(gids);
  apf::destroyNumbering(origin_part_ids);
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  m->destroyNative();
  apf::destroyMesh(m);
  //only free PCU if it isn't being used outside
  if (!pcu_outside)
    PCU_Comm_Free();
}

} // namespace Zoltan2

#endif // HAVE_ZOLTAN2_PARMA

#endif

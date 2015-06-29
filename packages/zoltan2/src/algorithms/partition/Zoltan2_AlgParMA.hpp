// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
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
//  This first design creates an apf mesh to run the ParMA algorithms on. The
//  final solution is determined by changes from beginning to end of the mesh.
//  This approach allows development closer to that of PUMI setup but at the 
//  cost of creating an extra mesh representation.
// 
//!  Available ParMA algorithms are given by setting the parma_method parameter 
//!  of the sublist parma_paramaters to one of the following:
//!  Vertex       - Balances targeting vertex imbalance
//!  Element      - Balances targeting element imbalance
//!  VtxElm       - Balances targeting vertex and element imbalance
//!  VtxEdgeElm   - Balances targeting vertex, edge, and element imbalance
//!  Ghost        - Balances using ghost element aware diffusion      
//!  Shape        - Optimizes shape of parts by increasing the size of small part boundaries
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

  AlgParMA(const RCP<const Environment> &env,
              const RCP<const Comm<int> > &problemComm,
           const RCP<const BaseAdapter<user_t> > &adapter)
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
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::userCoord_t userCoord_t;
  typedef typename Adapter::zgid_t zgid_t;

  const RCP<const Environment> env;
  const RCP<const Comm<int> > problemComm;
  const RCP<const MeshAdapter<user_t> > adapter;
  
  apf::Mesh2* m;
  apf::Numbering* gids;
  apf::Numbering* origin_part_ids;
  std::map<zgid_t, lno_t> mapping_elm_gids_index;

  MPI_Comm mpicomm;
  bool pcu_outside;

  void setMPIComm(const RCP<const Comm<int> > &problemComm__) {
#   ifdef HAVE_ZOLTAN2_MPI
      mpicomm = TeuchosConst2MPI(problemComm__);
#   else
      mpicomm = MPI_COMM_WORLD;  // taken from siMPI
#   endif
  }

  enum MeshEntityType entityAPFtoZ2(int dimension) const {return static_cast<MeshEntityType>(dimension);}

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

  void setEntWeights(int dim, apf::MeshTag* tag) {
    apf::MeshIterator* itr = m->begin(dim);
    apf::MeshEntity* ent;
    double w = 1.0;
    while ((ent= m->iterate(itr)))  {
      m->setDoubleTag(ent,tag,&w);
      assert(m->hasTag(ent,tag));
    }
    m->end(itr);
   
  }
  
  apf::MeshTag* setWeights(bool vtx, bool edge, bool elm) {
    apf::MeshTag* tag = m->createDoubleTag("parma_weight",1);
    if (vtx)
      setEntWeights(0,tag);
    if (edge)
      setEntWeights(1,tag);
    if (elm) {
      setEntWeights(m->getDimension(),tag);
    }
    return tag;
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
    pcu_outside=false;
    if (!PCU_Comm_Initialized())
      PCU_Comm_Init();
    else
      pcu_outside=true;
    PCU_Switch_Comm(mpicomm);

    int dim;
    //Get region topology types if its NULL then we have a 2D mesh
    const EntityTopologyType* tops;
    adapter->getTopologyViewOf(MESH_REGION,tops);
    if (tops==NULL)
      dim=2;
    else 
      dim=3;

    //Create empty apf mesh
    gmi_register_null();
    gmi_model* g = gmi_load(".null");
    enum MeshEntityType primary_type = entityAPFtoZ2(dim);
    m = apf::makeEmptyMdsMesh(g,dim,false);

    //Get entity topology types
    adapter->getTopologyViewOf(primary_type,tops);
    
    //Get element global ids and part ids
    const zgid_t* element_gids;
    const part_t* part_ids;
    adapter->getIDsViewOf(primary_type,element_gids);
    adapter->getPartsView(part_ids);
    for (size_t i =0;i<adapter->getLocalNumOf(primary_type);i++)
      mapping_elm_gids_index[element_gids[i]] = i;
    
    //get vertex global ids
    const zgid_t* vertex_gids;
    adapter->getIDsViewOf(MESH_VERTEX,vertex_gids);
    
    
    //Get vertex coordinates
    const scalar_t ** vertex_coords = new const scalar_t*[dim];
    int* strides = new int[dim];
    for (int i=0;i<dim;i++)
      adapter->getCoordinatesViewOf(MESH_VERTEX,vertex_coords[i],strides[i],i);

    //Get first adjacencies from elements to vertices
    if (!adapter->availAdjs(primary_type,MESH_VERTEX))
      throw "APF needs adjacency information from elements to vertices";
    const lno_t* offsets;
    const zgid_t* adjacent_vertex_gids;
    adapter->getAdjsView(primary_type, MESH_VERTEX,offsets,adjacent_vertex_gids);
    
    //build the apf mesh
    apf::GlobalToVert vertex_mapping;
    for (size_t i=0;i<adapter->getLocalNumOf(MESH_VERTEX);i++) {
      apf::MeshEntity* vtx = m->createVert(0);
      scalar_t temp_coords[3];
      for (int k=0;k<dim;k++) 
	temp_coords[k] = vertex_coords[k][i*strides[k]];

      for (int k=dim;k<3;k++)
	temp_coords[k] = 0;  
      apf::Vector3 point(temp_coords[0],temp_coords[1],temp_coords[2]);    
      m->setPoint(vtx,0,point);
      vertex_mapping[vertex_gids[i]] = vtx;
    }
    apf::Mesh::Type t = topologyZ2toAPF(tops[0]);
    apf::construct(m,adjacent_vertex_gids,adapter->getLocalNumOf(primary_type),t,vertex_mapping);
    
    
    //Setup numberings
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
    //cleanup
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
  //Get paramters
  std::string alg_name = "VtxElm";
  double imbalance = 1.1;
  double step = .1;
  int ghost_layers=3;
  int ghost_bridge=m->getDimension()-1;
 
  const Teuchos::ParameterList &pl = env->getParameters();
  try {
    const Teuchos::ParameterList &ppl = pl.sublist("parma_parameters");
    for (ParameterList::ConstIterator iter =  ppl.begin();
	 iter != ppl.end(); iter++) {
      const std::string &zname = pl.name(iter);
      // Convert the value to a string to pass to Zoltan
      if (zname == "parma_method") {
	std::string &zval = pl.entry(iter).getValue(&zval);
	alg_name = zval;
      }
      else if (zname == "step_size") {
	double &zval = pl.entry(iter).getValue(&zval);
	step = zval;
      }
      else if (zname=="ghost_layers" || zname=="ghost_bridge") {
	int &zval = pl.entry(iter).getValue(&zval);
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

  bool weightVertex,weightEdge,weightElement;
  weightVertex=weightEdge=weightElement=false;

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
    weightVertex = true;
  }
  else if (alg_name=="Shape") {
    balancer = Parma_MakeShapeOptimizer(m,step,verbose);
    weightElement=true;
  }
  else  {
    //Should be caught by the validator
    throw std::runtime_error("No such parma method defined");
  }
  apf::MeshTag* weights = setWeights(weightVertex,weightEdge,weightElement);

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
      zgid_t element_gid = apf::getNumber(gids,ent,0,0);
      PCU_COMM_PACK(target_part_id,element_gid);
    }
  }
  m->end(itr);
  //Send information off
  PCU_Comm_Send();
  //Unpack information and set new part ids
  while (PCU_Comm_Receive()) {
    zgid_t global_id;
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
  if (!pcu_outside)
    PCU_Comm_Free();
}

} // namespace Zoltan2

#endif // HAVE_ZOLTAN2_PARMA

#endif

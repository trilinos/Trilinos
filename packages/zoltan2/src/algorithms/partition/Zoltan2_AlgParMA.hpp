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

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Util.hpp>
#include <Zoltan2_TPLTraits.hpp>

#include <Zoltan2_AlgZoltanCallbacks.hpp>
#include <zoltan_cpp.h>

//////////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgParMA.hpp
//! \brief interface to the ParMA library
//  
//  This first design creates an apf mesh to run the ParMA algorithms on. The
//  final solution is determined by changes from beginning to end of the mesh.
//  This approach allows development closer to that of PUMI setup but at the 
//  cost of creating an extra mesh representation.
// 
//  Another approach might be to provide callbacks to modify the adapter 
//  in each iteration of a ParMA algorithm. This may be able to remove the 
//  need for an intermediate structure.
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

#include <apfMesh.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfMDS.h>
#include <PCU.h>
#include <parma.h>
#include <gmi_null.h>
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
  std::map<zgid_t, lno_t> mapping_gids_index;

  MPI_Comm mpicomm;
  
  void setMPIComm(const RCP<const Comm<int> > &problemComm__) {
#   ifdef HAVE_ZOLTAN2_MPI
      mpicomm = TeuchosConst2MPI(problemComm__);
#   else
      mpicomm = MPI_COMM_WORLD;  // taken from siMPI
#   endif
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
    throw("ParMA needs a MeshAdapter but you haven't given it one");
  }

  AlgParMA(const RCP<const Environment> &env__,
            const RCP<const Comm<int> > &problemComm__,
            const RCP<const VectorAdapter<user_t> > &adapter__) 
  { 
    throw("ParMA needs a MeshAdapter but you haven't given it one");
  }

  AlgParMA(const RCP<const Environment> &env__,
            const RCP<const Comm<int> > &problemComm__,
            const RCP<const GraphAdapter<user_t,userCoord_t> > &adapter__) 
  { 
    throw("ParMA needs a MeshAdapter but you haven't given it one");
  }

  AlgParMA(const RCP<const Environment> &env__,
            const RCP<const Comm<int> > &problemComm__,
            const RCP<const MatrixAdapter<user_t,userCoord_t> > &adapter__)
  { 
    throw("ParMA needs a MeshAdapter but you haven't given it one");
    
  }
  

  AlgParMA(const RCP<const Environment> &env__,
            const RCP<const Comm<int> > &problemComm__,
            const RCP<const MeshAdapter<user_t> > &adapter__) :
    env(env__), problemComm(problemComm__), adapter(adapter__)
  { 
    setMPIComm(problemComm__);
    
    //Setup numberings
    gids = apf::createNumbering(m,"global_ids",m->getShape(),1);
    origin_part_ids = apf::createNumbering(m,"origin",m->getShape(),1);

    //build the mesh
    gmi_register_null();
    gmi_model* g = gmi_load(".null");
    enum MeshEntityType primary_type = adapter->getPrimaryEntityType();
    m = apf::makeEmptyMdsMesh(g,adapter->getDimension(),false);

    const zgid_t* element_gids;
    const part_t* part_ids;
    adapter->getIDsViewOf(primary_type,element_gids);
    adapter->getPartsView(part_ids);

    const zgid_t* vertex_gids;
    adapter->getIDsViewOf(MESH_VERTEX,vertex_gids);
    
    for (size_t i =0;i<adapter->getLocalNumOf(MESH_VERTEX);i++)
      mapping_gids_index[vertex_gids[i]] = i;

    const scalar_t ** vertex_coords = new const scalar_t*[adapter->getDimension()];
    int stride;
    for (int i=0;i<adapter->getDimension();i++)
      adapter->getCoordinatesViewOf(MESH_VERTEX,vertex_coords[i],stride,i);

    const lno_t* offsets;
    const zgid_t* adjacent_vertex_gids;
    adapter->getAdjsView(primary_type, MESH_VERTEX,offsets,adjacent_vertex_gids);

    for (size_t i=0;i<adapter->getLocalNumOf(primary_type);i++) {
      lno_t num_verts = offsets[i+1]-offsets[i];
      apf::MeshEntity** vertices = new apf::MeshEntity*[num_verts];
      for (lno_t j=0; j<num_verts;j++) {
	lno_t vertex_lid = mapping_gids_index[adjacent_vertex_gids[offsets[i]+j]];
	scalar_t temp_coords[3];
	for (int k=0;k<adapter->getDimension();k++)
	  temp_coords[k] = vertex_coords[k][vertex_lid];
	for (int k=adapter->getDimension();k<3;k++)
	  temp_coords[k] = 0;
	    
	apf::Vector3 point(temp_coords[0],temp_coords[1],temp_coords[2]);
	vertices[j] = m->createVert(0);
	m->setPoint(vertices[j],0,point);
      }
      apf::MeshEntity* element = apf::buildElement(m, 0, apf::Mesh::TET, vertices);
      apf::number(gids,element,0,0,element_gids[i]);
      apf::number(origin_part_ids,element,0,0,part_ids[i]);
      
      delete [] vertices;
    }
    apf::deriveMdsModel(m);
    m->acceptChanges();
    m->verify();
    
  }
  void partition(const RCP<PartitioningSolution<Adapter> > &solution);
  
};

/////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void AlgParMA<Adapter>::partition(
  const RCP<PartitioningSolution<Adapter> > &solution
)
{
  //TODO See if we need to change the comm of PCU
  PCU_Comm_Init();
  apf::Balancer* balancer;
  apf::MeshTag* weights= NULL;
  if (true) {
    //weights = setWeights(m);
    const double step = 0.1; const int verbose = 1;
    balancer = Parma_MakeElmBalancer(m, step, verbose);
  }


  balancer->balance(weights, 1.05);
  delete balancer;

  /*
  env->globalInputAssertion(__FILE__, __LINE__, "Zoltan LB_Partition", 
    (ierr==ZOLTAN_OK || ierr==ZOLTAN_WARN), BASIC_ASSERTION, problemComm);
  */

  // Load answer into the solution.
  int num_local = adapter->getLocalNumOf(adapter->getPrimaryEntityType());
  ArrayRCP<part_t> partList(new part_t[num_local], 0, num_local, true);

  PCU_Comm_Begin();
  apf::MeshEntity* ent;
  apf::MeshIterator* itr = m->begin(m->getDimension());
  while ((ent=m->iterate(itr))) {
    if (m->isOwned(ent)) {
      part_t target_part_id = apf::getNumber(origin_part_ids,ent,0,0);
      zgid_t element_gid = apf::getNumber(gids,ent,0,0);
      PCU_COMM_PACK(target_part_id,element_gid);
    }
  }

  PCU_Comm_Send();
  while (PCU_Comm_Listen()) {
    zgid_t global_id;
    PCU_COMM_UNPACK(global_id);
    lno_t local_id = mapping_gids_index[global_id];
    part_t new_part_id = PCU_Comm_Sender();
    partList[local_id] = new_part_id;
  }

  solution->setParts(partList);
  

  // Clean up
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);

  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
}

} // namespace Zoltan2

#endif // HAVE_ZOLTAN2_PARMA

#endif

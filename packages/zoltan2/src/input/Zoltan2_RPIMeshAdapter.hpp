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

/*! \file Zoltan2_RPIMeshAdapter.hpp
    \brief Defines the RPIMeshAdapter class.
*/

#ifndef _ZOLTAN2_RPIMESHADAPTER_HPP_
#define _ZOLTAN2_RPIMESHADAPTER_HPP_

#include <Zoltan2_MeshAdapter.hpp>
#include <Zoltan2_StridedData.hpp>
#include <map>
#include <vector>
#include <string>
#include <cassert>

#ifndef HAVE_ZOLTAN2_PARMA

namespace apf {
  class Mesh;
}
namespace Zoltan2 {
template <typename User>
class RPIMeshAdapter : public MeshAdapter<User>
{
public:
  
  
  RPIMeshAdapter(const Comm<int> &comm, apf::Mesh* m)
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
#include <apfDynamicArray.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <PCU.h>
namespace Zoltan2 {

/*! \brief This class represents a mesh.
 *
 *  A mesh can be a collection of global Identifiers
 *           and their associated weights, if any.
 *
 *  The user supplies the identifiers and weights by way of pointers
 *    to arrays.  
 *
    The template parameter (\c User) is a C++ class type which provides the
    actual data types with which the Zoltan2 library will be compiled, through
    a Traits mechanism.  \c User may be the
    actual class used by application to represent coordinates, or it may be
    the empty helper class \c BasicUserTypes with which a Zoltan2 user
    can easily supply the data types for the library.

    The \c scalar_t type, representing use data such as matrix values, is
    used by Zoltan2 for weights, coordinates, part sizes and
    quality metrics.
    Some User types (like Tpetra::CrsMatrix) have an inherent scalar type,
    and some
    (like Tpetra::CrsGraph) do not.  For such objects, the scalar type is
    set by Zoltan2 to \c float.  If you wish to change it to double, set
    the second template parameter to \c double.

 */

template <typename User>
  class RPIMeshAdapter: public MeshAdapter<User> {

public:

  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::zgid_t    zgid_t;
  typedef typename InputTraits<User>::part_t   part_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef MeshAdapter<User>       base_adapter_t;
  typedef User user_t;

  /*! \brief Constructor for mesh with identifiers but no coordinates or edges
   *  \param etype is the mesh entity type of the identifiers
   *
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this InputAdapter.
   */

  RPIMeshAdapter(const Comm<int> &comm, apf::Mesh* m);
  void destroy();
  void print(int);
  template <typename Adapter>
  void applyPartitioningSolution(const User &in, User *&out,
    const PartitioningSolution<Adapter> &solution) const{

    apf::Migration* plan = new apf::Migration(*out);
    const part_t* new_part_ids = solution.getPartListView();

    apf::MeshIterator* itr = (*out)->begin(m_dimension);
    apf::MeshEntity* ent;
    int i=0;
    while ((ent=(*out)->iterate(itr)))  {
      plan->send(ent,new_part_ids[i]);
      i++;
    }
    (*out)->migrate(plan);
  }

  ////////////////////////////////////////////////////////////////
  // The MeshAdapter interface.
  // This is the interface that would be called by a model or a problem .
  ////////////////////////////////////////////////////////////////

  size_t getGlobalNumOf(MeshEntityType etype) const
  {
    int dim = entityZ2toAPF(etype);
    if (dim<=m_dimension&&dim>=0)
      return num_global[dim];
    return 0;
  }

  size_t getLocalNumOf(MeshEntityType etype) const
  {
    int dim = entityZ2toAPF(etype);
    if (dim<=m_dimension&&dim>=0)
      return num_local[dim];
    return 0;
  }
   
  void getIDsViewOf(MeshEntityType etype, const zgid_t *&Ids) const
  {
    int dim = entityZ2toAPF(etype);
    if (dim<=m_dimension&&dim>=0)
      Ids = gid_mapping[dim];
    else
      Ids = NULL;
  }

  void getTopologyViewOf(MeshEntityType etype,
                         enum EntityTopologyType const *&Types) const {
    int dim = entityZ2toAPF(etype);
    if (dim<=m_dimension&&dim>=0)
      Types = topologies[dim];
    else
      Types = NULL;
  }

  void getWeightsViewOf(MeshEntityType etype, const scalar_t *&weights,
                        int &stride, int idx = 0) const
  {
    weights = NULL;
    stride = 0;
  }

  int getDimension() const { return coord_dimension; }

  void getCoordinatesViewOf(MeshEntityType etype, const scalar_t *&coords,
                            int &stride, int coordDim) const {
    if (coordDim>=0 && coordDim<3) {
      int dim = entityZ2toAPF(etype);
      if (dim<=m_dimension&&dim>=0) {
        coords = ent_coords[dim]+coordDim;
        stride = 3;
      }
    }
    else {
      coords = NULL;
      stride = 0;
    }
  }

  bool availAdjs(MeshEntityType source, MeshEntityType target) const {
    int dim_source = entityZ2toAPF(source);
    int dim_target = entityZ2toAPF(target);
    return dim_source<=m_dimension && dim_source>=0 &&
      dim_target<=m_dimension && dim_target>=0 &&
      dim_target!=dim_source;
  }

  size_t getLocalNumAdjs(MeshEntityType source, MeshEntityType target) const
  {
    int dim_source = entityZ2toAPF(source);
    int dim_target = entityZ2toAPF(target);
    if (availAdjs(source,target)) 
      return adj_gids[dim_source][dim_target].size();
    return 0;
  }

  void getAdjsView(MeshEntityType source, MeshEntityType target,
		   const lno_t *&offsets, const zgid_t *& adjacencyIds) const
  {
    int dim_source = entityZ2toAPF(source);
    int dim_target = entityZ2toAPF(target);
    if (availAdjs(source,target)) {
      offsets = adj_offsets[dim_source][dim_target];
      adjacencyIds = &(adj_gids[dim_source][dim_target][0]);
    } 
    else {
      offsets=NULL;
      adjacencyIds = NULL;
    }
  }
  //#define USE_MESH_ADAPTER
#ifndef USE_MESH_ADAPTER
  bool avail2ndAdjs(MeshEntityType sourcetarget, MeshEntityType through) const
  {
    int dim_source = entityZ2toAPF(sourcetarget);
    int dim_target = entityZ2toAPF(through);
    return dim_source<=m_dimension && dim_source>=0 &&
      dim_target<=m_dimension && dim_target>=0 &&
      dim_target!=dim_source;
  }

  size_t getLocalNum2ndAdjs(MeshEntityType sourcetarget, 
			    MeshEntityType through) const
  {
    int dim_source = entityZ2toAPF(sourcetarget);
    int dim_target = entityZ2toAPF(through);
    if (avail2ndAdjs(sourcetarget,through)) 
      return adj2_gids[dim_source][dim_target].size();
    return 0;

  }

  void get2ndAdjsView(MeshEntityType sourcetarget, MeshEntityType through, 
		      const lno_t *&offsets, const zgid_t *&adjacencyIds) const
  {
    int dim_source = entityZ2toAPF(sourcetarget);
    int dim_target = entityZ2toAPF(through);
    if (avail2ndAdjs(sourcetarget,through)) {
      offsets=adj2_offsets[dim_source][dim_target];
      adjacencyIds=&(adj2_gids[dim_source][dim_target][0]);
    }
    
  }
#endif

private:
  int entityZ2toAPF(enum MeshEntityType etype) const {return static_cast<int>(etype);}
  enum EntityTopologyType topologyAPFtoZ2(enum apf::Mesh::Type ttype) const {
    if (ttype==apf::Mesh::VERTEX)
      return POINT;
    else if (ttype==apf::Mesh::EDGE)
      return LINE_SEGMENT;
    else if (ttype==apf::Mesh::TRIANGLE)
      return TRIANGLE;
    else if (ttype==apf::Mesh::QUAD)
      return QUADRILATERAL;
    else if (ttype==apf::Mesh::TET)
      return TETRAHEDRON;
    else if (ttype==apf::Mesh::HEX)
      return HEXAHEDRON;
    else if (ttype==apf::Mesh::PRISM)
      return PRISM;
    else if (ttype==apf::Mesh::PYRAMID)
      return PYRAMID;
    else 
      throw "No such APF topology type";
    
  }
  

  int m_dimension;  //Dimension of the mesh
  apf::Numbering** lids; //[dimension] array of local id numbers
  apf::GlobalNumbering** gids;//[dimension] array of global id numbers
  zgid_t** gid_mapping; //[dimension][lid] corresponding global id numbers
  size_t* num_local; //[dimension] number of local entities
  size_t* num_global; //[dimension] number of global entities
  EntityTopologyType** topologies; //[dimension] topologies for each entity
  lno_t*** adj_offsets; //[first_dimension][second_dimension] array of offsets
  std::vector<zgid_t>** adj_gids; //[first_dimension][second_dimension] global_ids of first adjacencies
  lno_t*** adj2_offsets; //[first_dimension][second_dimension] array of offsets for second adjacencies
  std::vector<zgid_t>** adj2_gids; //[first_dimension][second_dimension] global_ids of second adjacencies
  int coord_dimension; //dimension of coordinates (always 3 for APF)
  scalar_t** ent_coords; //[dimension] array of coordinates [xs ys zs]
};

////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////

template <typename User>
RPIMeshAdapter<User>::RPIMeshAdapter(const Comm<int> &comm,
                                     apf::Mesh* m) {
  
  
  //mesh dimension
  m_dimension = m->getDimension();

  if (m_dimension==3)
    this->setEntityTypes("region","face","face");
  else
    this->setEntityTypes("face","edge","edge");
  //count the local and global numbers as well as assign ids and map local to global
  lids = new apf::Numbering*[m_dimension+1];
  gids = new apf::GlobalNumbering*[m_dimension+1];
  gid_mapping = new zgid_t*[m_dimension+1];
  std::map<zgid_t,lno_t>* lid_mapping = new std::map<zgid_t,lno_t>[m_dimension+1];
  num_local = new size_t[m_dimension+1];
  num_global = new size_t[m_dimension+1];
  topologies = new EntityTopologyType*[m_dimension+1];

  for (int i=0;i<=m_dimension;i++) {  
    //number of local and global entities
    long global_count = num_local[i] = m->count(i);
    PCU_Add_Longs(&global_count,1);
    num_global[i] = global_count;
    
    //Number each entity with local and global numbers
    char lids_name[15];
    sprintf(lids_name,"lids%d",i);
    char gids_name[15];
    sprintf(gids_name,"ids%d",i);
    apf::FieldShape* shape = apf::getConstant(i);
    lids[i] = apf::createNumbering(m,lids_name,shape,1);
    apf::Numbering* tmp = apf::numberOwnedDimension(m,gids_name,i);
    gids[i] = apf::makeGlobal(tmp);
    apf::synchronize(gids[i]);
    apf::MeshIterator* itr = m->begin(i);
    apf::MeshEntity* ent;
    unsigned int num=0;
    while ((ent=m->iterate(itr))) {
      apf::number(lids[i],ent,0,0,num); 
      lid_mapping[i][apf::getNumber(gids[i],apf::Node(ent,0))]=num;
      num++;
    }
    m->end(itr);
    assert(num==num_local[i]);
    
    //Make a mapping from local to global
    //While we are at it take the topology types
    gid_mapping[i] = new zgid_t[num_local[i]];
    topologies[i] = new EntityTopologyType[num_local[i]];
    apf::DynamicArray<apf::Node> nodes;
    itr = m->begin(i);
    num=0;
    while((ent=m->iterate(itr)))  {
      zgid_t gid = apf::getNumber(gids[i],apf::Node(ent,0));
      gid_mapping[i][ apf::getNumber(lids[i],ent,0,0)] = gid;
      topologies[i][num] = topologyAPFtoZ2(static_cast<apf::Mesh::Type>(m->getType(ent)));
      num++;
    }
    m->end(itr);
   
    
  }
  //First Adjacency and Second Adjacency data
  adj_gids = new std::vector<zgid_t>*[m_dimension+1];
  adj_offsets = new lno_t**[m_dimension+1];
  adj2_gids = new std::vector<zgid_t>*[m_dimension+1];
  adj2_offsets = new lno_t**[m_dimension+1];
      
  for (int i=0;i<=m_dimension;i++) {
    adj_gids[i] = new std::vector<zgid_t>[m_dimension+1];
    adj_offsets[i] = new lno_t*[m_dimension+1];
    adj2_gids[i] = new std::vector<zgid_t>[m_dimension+1];
    adj2_offsets[i] = new lno_t*[m_dimension+1];
    for (int j=0;j<=m_dimension;j++) {
      
      if (i==j) {
        adj_offsets[i][j]=NULL;
        adj2_offsets[i][j]=NULL;
        continue;
      }
      
      //Loop through each entity
      apf::MeshIterator* itr = m->begin(i);
      apf::MeshEntity* ent;
      
      adj_offsets[i][j] = new lno_t[num_local[i]+1];
      adj_offsets[i][j][0] =0;
      adj2_offsets[i][j] = new lno_t[num_local[i]+1];
      adj2_offsets[i][j][0] =0;
      int k=1;
      
      //We need communication for second adjacency
      PCU_Comm_Begin();
      std::map<zgid_t,apf::MeshEntity*> part_boundary_mapping;
      while ((ent=m->iterate(itr))) {
        std::set<zgid_t> temp_adjs; //temp storage for second adjacency
        //Get First Adjacency
        apf::Adjacent adj;
        m->getAdjacent(ent,j,adj);
        for (unsigned int l=0;l<adj.getSize();l++) {
          adj_gids[i][j].push_back(apf::getNumber(gids[j],apf::Node(adj[l],0)));
          //Now look at Second Adjacency
          apf::Adjacent adj2;
          m->getAdjacent(adj[l],i,adj2);
          for (unsigned int o=0;o<adj2.getSize();o++)
            temp_adjs.insert(apf::getNumber(gids[i],apf::Node(adj2[o],0)));
	  if (i==m_dimension) {
	    apf::Parts res;
	    m->getResidence(adj[l],res);
	     
	    part_boundary_mapping[apf::getNumber(gids[j],apf::Node(adj[l],0))] = adj[l];
	    for (apf::Parts::iterator it=res.begin();it!=res.end();it++) {
	      zgid_t send_vals[2];
	      send_vals[1]=apf::getNumber(gids[i],apf::Node(ent,0));
	      send_vals[0]=apf::getNumber(gids[j],apf::Node(adj[l],0));
	      
	      PCU_Comm_Pack(*it,send_vals,2*sizeof(zgid_t));
	    }
	  }
        }
        adj_offsets[i][j][k] = adj_gids[i][j].size();
        k++;
        //Copy over local second adjacencies to copies
	if (i!=m_dimension) {
	  apf::Parts res;
	  m->getResidence(ent,res);
	  typename std::set<zgid_t>::iterator adj_itr;
	  for (adj_itr=temp_adjs.begin();adj_itr!=temp_adjs.end();adj_itr++) {
	    for (apf::Parts::iterator it=res.begin();it!=res.end();it++) {
	      zgid_t send_vals[2];
	      send_vals[0]=apf::getNumber(gids[i],apf::Node(ent,0));
	      send_vals[1] = *adj_itr;
	      if (send_vals[0]!=send_vals[1])
		PCU_Comm_Pack(*it,send_vals,2*sizeof(zgid_t));
	    }
          }
        }
      }
      m->end(itr);
      //Now capture mesh wide second adjacency locally
      PCU_Comm_Send();
      std::set<zgid_t>* adjs2 = new std::set<zgid_t>[num_local[i]];
      while (PCU_Comm_Receive()) {
        zgid_t adj2[2];
        PCU_Comm_Unpack(adj2,2*sizeof(zgid_t));
	
	if (i==m_dimension) {
	  apf::MeshEntity* through = part_boundary_mapping[adj2[0]];
	  apf::Adjacent adj;
	  m->getAdjacent(through,i,adj);
	  for (unsigned int l=0;l<adj.getSize();l++) {
	    if (apf::getNumber(gids[i],apf::Node(adj[l],0))!=adj2[1])
	      adjs2[apf::getNumber(lids[i],adj[l],0,0)].insert(adj2[1]);
	  }
	}
	else {
	  lno_t index = lid_mapping[i][adj2[0]];
	  adjs2[index].insert(adj2[1]);
	}
      }      
      //And finally convert the second adjacency to a vector to be returned to user
      for (size_t l=0;l<num_local[i];l++) {
        for (typename std::set<zgid_t>::iterator sitr = adjs2[l].begin();sitr!=adjs2[l].end();sitr++) {
          adj2_gids[i][j].push_back(*sitr);
        }
        adj2_offsets[i][j][l+1]=adj2_gids[i][j].size();
      }
    }
  }
  //Coordinates
  coord_dimension = 3;
  ent_coords = new scalar_t*[m_dimension+1];
  for (int i=0;i<=m_dimension;i++) {
    apf::MeshIterator* itr = m->begin(i);
    apf::MeshEntity* ent;
    ent_coords[i] = new scalar_t[3*num_local[i]];
    int j=0;
    while((ent=m->iterate(itr))) {
      apf::Vector3 point;
      if (i==0) {
        m->getPoint(ent,0,point);
      }
      else {
        point = apf::getLinearCentroid(m,ent);
      }
      for (int k=0;k<3;k++)  
        ent_coords[i][j*3+k] = point[k];
      j++;
    }
    m->end(itr);
  }
  delete [] lid_mapping;
}
template <typename User>
void RPIMeshAdapter<User>::destroy() {
  for (int i=0;i<=m_dimension;i++) {
    delete [] ent_coords[i];
    delete [] adj_gids[i];
    delete [] adj2_gids[i];
    for (int j=0;j<=m_dimension;j++) {
      if (i!=j) {
        delete [] adj_offsets[i][j];
	delete [] adj2_offsets[i][j];
      }
    }
    delete [] adj2_offsets[i];
    delete [] adj_offsets[i];
    delete [] gid_mapping[i];
    apf::destroyGlobalNumbering(gids[i]);
    apf::destroyNumbering(lids[i]);
  }
  delete [] ent_coords;
  delete [] adj_gids;
  delete [] adj_offsets;
  delete [] adj2_gids;
  delete [] adj2_offsets;
  delete [] gid_mapping;
  delete [] gids;
  delete [] lids;
  delete [] num_local;
  delete [] num_global;
  
     
}  
template <typename User>
void RPIMeshAdapter<User>::print(int me)
{
  std::string fn("RPIMesh ");
  std::cout << me << fn 
            << " dimension = " << m_dimension
            << std::endl;
  for (int i=0;i<=m_dimension;i++)
    std::cout<<me<<" Number of dimension " << i<< " = " <<num_local[i] <<std::endl;
}

 

#endif //HAVE_ZOLTAN2_PARMA
  
}  //namespace Zoltan2
  
#endif

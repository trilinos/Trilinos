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
//
// Basic testing of Zoltan2::APFMeshAdapter

#include <Zoltan2_APFMeshAdapter.hpp>

#ifdef HAVE_ZOLTAN2_PARMA
#include <apf.h>
#include <apfMesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <gmi_mesh.h>
#endif

// Teuchos includes
#include "Teuchos_XMLParameterListHelpers.hpp"

using Teuchos::RCP;

/*********************************************************/
/*                     Typedefs                          */
/*********************************************************/
//Tpetra typedefs
typedef Tpetra::MultiVector<double, int, int>     tMVector_t;

//Topology type helper function
enum Zoltan2::EntityTopologyType topologyAPFtoZ2(enum apf::Mesh::Type ttype) {
  if (ttype==apf::Mesh::VERTEX)
    return Zoltan2::POINT;
  else if (ttype==apf::Mesh::EDGE)
    return Zoltan2::LINE_SEGMENT;
  else if (ttype==apf::Mesh::TRIANGLE)
    return Zoltan2::TRIANGLE;
  else if (ttype==apf::Mesh::QUAD)
    return Zoltan2::QUADRILATERAL;
  else if (ttype==apf::Mesh::TET)
    return Zoltan2::TETRAHEDRON;
  else if (ttype==apf::Mesh::HEX)
    return Zoltan2::HEXAHEDRON;
  else if (ttype==apf::Mesh::PRISM)
    return Zoltan2::PRISM;
  else if (ttype==apf::Mesh::PYRAMID)
    return Zoltan2::PYRAMID;
  else
    throw "No such APF topology type";

}

/*****************************************************************************/
/******************************** MAIN ***************************************/
/*****************************************************************************/

int main(int narg, char *arg[]) {

  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > CommT = Tpetra::getDefaultComm();

#ifdef HAVE_ZOLTAN2_PARMA
  //Open up PCU for communication required by all APF operations
  PCU_Comm_Init();

  //Open up the mesh
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh("pumiTri14/plate.dmg","pumiTri14/2/");
  apf::verify(m);

  int dim = m->getDimension();

  //Contruct the MeshAdapter
  typedef Zoltan2::APFMeshAdapter<apf::Mesh2*> Adapter;
  typedef Adapter::lno_t lno_t;
  typedef Adapter::gno_t gno_t;
  typedef Adapter::scalar_t scalar_t;
  typedef Adapter::offset_t offset_t;

  std::string pri = "face";
  std::string adj = "vertex";
  if (dim==3) {
    adj=pri;
    pri="region";
  }

  Zoltan2::APFMeshAdapter<apf::Mesh2*> ia(*CommT,m,pri,adj,true);

  Zoltan2::MeshEntityType ents[4] = {Zoltan2::MESH_VERTEX,
                                     Zoltan2::MESH_EDGE,
                                     Zoltan2::MESH_FACE,
                                     Zoltan2::MESH_REGION};

  //Check the local number of each entity
  bool* has = new bool[dim+1];
  for (int i=0;i<=dim;i++) {

    has[i]=true;
    if (ia.getLocalNumOf(ents[i])==0) {
      has[i]=false;
      continue;
    }
    if (ia.getLocalNumOf(ents[i])!=m->count(i)) {
      std::cerr<<"Local number of entities does not match in dimension "<<i<<"\n";
      return 1;
    }

  }


  //Check the coordinate dimension
  apf::GlobalNumbering** gnums = new apf::GlobalNumbering*[dim];
  apf::Numbering** lnums = new apf::Numbering*[dim];
  int sub=0;
  for (int i=0;i<=dim;i++) {
    if (!has[i]) {
      sub++;
      continue;
    }
    gnums[i] = m->getGlobalNumbering(i-sub);
    lnums[i] = m->getNumbering(i-sub);
  }

  for (int i=0;i<=dim;i++) {
    if (!has[i])
      continue;
    const gno_t* gids;
    const Zoltan2::EntityTopologyType* topTypes;
    const scalar_t* x_coords;
    const scalar_t* y_coords;
    const scalar_t* z_coords;
    int x_stride;
    int y_stride;
    int z_stride;
    const offset_t** offsets = new const offset_t*[dim];
    const gno_t** adj_gids = new const gno_t*[dim];

    ia.getIDsViewOf(ents[i],gids);
    ia.getTopologyViewOf(ents[i],topTypes);
    ia.getCoordinatesViewOf(ents[i],x_coords,x_stride,0);
    ia.getCoordinatesViewOf(ents[i],y_coords,y_stride,1);
    ia.getCoordinatesViewOf(ents[i],z_coords,z_stride,2);
    //Check availability of First Adjacency
    for (int j=0;j<=dim;j++) {
      if (!has[j])
        continue;
      if (ia.availAdjs(ents[i],ents[j])!=(i!=j)) {
        std::cerr<<"First Adjacency does not exist from "<<i<<" to "<< j<<"\n";
        return 5;
      }
      ia.getAdjsView(ents[i],ents[j],offsets[j],adj_gids[j]);
    }
    int j=0;
    apf::MeshIterator* itr = m->begin(i);
    apf::MeshEntity* ent;
    size_t* numAdjs = new size_t[dim+1];
    for (int k=0;k<=dim;k++)
      numAdjs[k]=0;
    while ((ent=m->iterate(itr))) {
      //Check Local ID numbers
      if (apf::getNumber(lnums[i],ent,0,0)!=j) {
        std::cerr<<"Local numbering does not match in dimension "<<i<<" on entity "<<j<<"\n";
        return 2;
      }
      //Check Global Id numbers
      if (apf::getNumber(gnums[i],apf::Node(ent,0))!=gids[j]) {
        std::cerr<<"Global numbering does not match in dimension "<<i<<" on entity "<<j<<"\n";
        return 2;
      }

      //Check Topology Types
      if (topologyAPFtoZ2(m->getType(ent))!=topTypes[j]) {
        std::cerr<<"Topology types do not match in dimension "<<i<<" on entity "<<j<<"\n";
        return 3;
      }

      //Check the coordinates
      apf::Vector3 pnt;
      if (i==0)
        m->getPoint(ent,0,pnt);
      else
        pnt = apf::getLinearCentroid(m,ent);
      float eps=.00001;
      if (fabs(pnt[0] - x_coords[j*x_stride])>eps) {
        std::cerr<<"X coordinate do not match in dimension "<<i<<" on entity "<<j<<"\n";
        return 4;
      }
      if (fabs(pnt[1] - y_coords[j*y_stride])>eps) {
        std::cerr<<"Y coordinate do not match in dimension "<<i<<" on entity "<<j<<"\n";
        return 4;
      }
      if (fabs(pnt[2] - z_coords[j*z_stride])>eps) {
        std::cerr<<"Z coordinate do not match in dimension "<<i<<" on entity "<<j<<"\n";
        return 4;
      }

      //Check first adjacencies
      for (int k=0;k<=dim;k++) {
        if (!has[k])
          continue;
        if (i==k) {
          //Check first adjacency to self is set to NULL
          if (offsets[k]!=NULL || adj_gids[k]!=NULL) {
            std::cerr<<"[WARNING] First adjacency to self is not set to NULL in dimension "<<i
                     <<" to dimension "<<k<<"\n";
          }

          continue;
        }
        apf::Adjacent adjs;
        m->getAdjacent(ent,k,adjs);
        offset_t ind = offsets[k][j];
        for (unsigned int l=0;l<adjs.getSize();l++) {
          if (apf::getNumber(gnums[k],apf::Node(adjs[l],0))!=adj_gids[k][ind]) {
            std::cerr<<"First adjacency does not match in dimension " <<i<<" to dimension "<<k
                     <<" on entity "<<j<<"\n";
            return 7;
          }
          ind++;
        }
        if (ind!=offsets[k][j+1]) {
          std::cerr<<"First adjacency length does not match in dimension "<<i<<" to dimension "
                   <<k<<" on entity "<<j<<"\n";
          return 8;
        }
        numAdjs[k]+=adjs.getSize();

      }
      j++;
    }
    m->end(itr);
    delete [] offsets;
    delete [] adj_gids;
    //Check the number of first adjacency
    for (int k=0;k<=dim;k++) {
      if (!has[k])
        continue;
      if (ia.getLocalNumAdjs(ents[i],ents[k])!=numAdjs[k]) {
        std::cerr<<"Local number of first adjacencies do not match in dimension "<<i
                 <<" through dimension "<<k<<"\n";
        return 6;
      }
    }
    delete [] numAdjs;

  }
  delete [] lnums;
  delete [] gnums;
  const Adapter::scalar_t arr[] = {1,2,1,3,1,5,1,2,1,6,1,7,1,8};
  ia.setWeights(Zoltan2::MESH_FACE,arr,2);

  if (ia.getNumWeightsPerOf(Zoltan2::MESH_FACE)!=1) {
    std::cerr<<"Number of weights incorrect\n";
    return 9;

  }


  const Adapter::scalar_t* arr2;
  int stride;
  ia.getWeightsViewOf(Zoltan2::MESH_FACE,arr2,stride);
  for (int i=0;i<7;i++) {
    if (arr[i*stride]!=arr2[i*stride]) {
      std::cerr<<"Weights do not match the original input\n";
      return 10;

    }
  }
  bool ifIdontfeellikeit = false;
  apf::MeshTag* weights = m->createDoubleTag("weights",2);
  apf::MeshIterator* itr = m->begin(0);
  apf::MeshEntity* ent;
  while ((ent=m->iterate(itr))) {
    if (!ifIdontfeellikeit||PCU_Comm_Self()) {
      double w[]={1.0,PCU_Comm_Self()+1.0};
      m->setDoubleTag(ent,weights,w);
    }
    ifIdontfeellikeit=!ifIdontfeellikeit;
  }
  m->end(itr);



  ia.setWeights(Zoltan2::MESH_VERTEX,m,weights);
  if (ia.getNumWeightsPerOf(Zoltan2::MESH_VERTEX)!=2) {
    std::cerr<<"Number of weights incorrect\n";
    return 9;

  }

  ia.getWeightsViewOf(Zoltan2::MESH_VERTEX,arr2,stride,0);

  itr = m->begin(0);

  int i=0;
  while ((ent=m->iterate(itr))) {
    double w=1;
    if (w!=arr2[i*stride]) {
      std::cerr<<"Weights do not match the original input\n";
      return 10;

    }
    i++;
  }
  m->end(itr);

  ia.getWeightsViewOf(Zoltan2::MESH_VERTEX,arr2,stride,1);

  itr = m->begin(0);
  i=0;
  while ((ent=m->iterate(itr))) {
    double w[2];
    if (PCU_Comm_Self())
      m->getDoubleTag(ent,weights,w);
    else
      w[1]=1;
    if (w[1]!=arr2[i*stride]) {
      std::cerr<<"Weights do not match the original input\n";
      return 10;

    }
    i++;
  }
  m->end(itr);

  //ia.print(PCU_Comm_Self(),5);

  //Delete the MeshAdapter
  ia.destroy();

  //Delete the APF Mesh
  m->destroyNative();
  apf::destroyMesh(m);

  //End PCU communications
  PCU_Comm_Free();
#endif
  std::cout<<"PASS\n";

  return 0;
}
/*****************************************************************************/
/********************************* END MAIN **********************************/
/*****************************************************************************/

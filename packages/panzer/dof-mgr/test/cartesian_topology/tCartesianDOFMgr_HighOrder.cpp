// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include "Kokkos_Core.hpp"
#include "Kokkos_DynRankView.hpp"

#include "Intrepid2_config.h"
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"

#include "PanzerCore_config.hpp"

#include "Panzer_ConnManager.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_UniqueGlobalIndexer_Utilities.hpp"

#include "CartesianConnManager.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {
namespace unit_test {

typedef CartesianConnManager<int,Ordinal64>::Triplet<Ordinal64> Triplet;

RCP<const panzer::FieldPattern> buildFieldPattern(RCP<Intrepid2::Basis<PHX::Device,double,double> > basis)
{
  // build a geometric pattern from a single basis
  RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
  return pattern;
}

std::string getElementBlock(const Triplet & element,
                                    const CartesianConnManager<int,Ordinal64> & connManager)
                                    
{
  int localElmtId = connManager.computeLocalElementIndex(element); 
  return connManager.getBlockId(localElmtId);
}

TEUCHOS_UNIT_TEST(tCartesianDOFMgr_HighOrder, ho_gid_values)
{
  typedef CartesianConnManager<int,Ordinal64> CCM;
  typedef panzer::DOFManager<int,Ordinal64> DOFManager;

  // build global (or serial communicator)
  #ifdef HAVE_MPI
    Teuchos::MpiComm<int> comm(MPI_COMM_WORLD);
  #else
    THIS_REALLY_DOES_NOT_WORK
  #endif

  int np = comm.getSize(); // number of processors
  //int rank = comm.getRank(); // processor rank

  // mesh description
  Ordinal64 nx = 8, ny = 4;//, nz = 4;
  //Ordinal64 nx = 4, ny = 3;//, nz = 4;
  int px = np, py = 1;//, pz = 1; // npx1 processor grids
  int bx =  1, by = 1;//, bz = 1; // 1x2 blocks

  const int poly_U = 4;
  const int poly_P = 3;
  RCP<const panzer::FieldPattern> pattern_U = buildFieldPattern( rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>( poly_U )) );
  RCP<const panzer::FieldPattern> pattern_P = buildFieldPattern( rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>( poly_P )) );

  // build the topology
  RCP<CCM> connManager = rcp(new CCM);
  connManager->initialize(comm,nx,ny,px,py,bx,by);

  // build the dof manager, and assocaite with the topology
  RCP<DOFManager> dofManager = rcp(new DOFManager);
  dofManager->setConnManager(connManager,*comm.getRawMpiComm());

  // add velocity (U) and PRESSURE fields to the MHD element block
  dofManager->addField("eblock-0_0","U",pattern_U);
  dofManager->addField("eblock-0_0","P",pattern_P);

  // build global unknowns (useful comment!)
  dofManager->buildGlobalUnknowns();

  auto myOffset   = connManager->getMyOffsetTriplet();
  auto myElements = connManager->getMyElementsTriplet();

  // check sharing locally on this processor
  {
    // choose an element in the middle of the day
    Triplet element;
    element.x = myOffset.x + myElements.x/2-1;
    element.y = myOffset.y + myElements.y/2-1;
    element.z = 0; 

    out << "Root element = " << element.x << " " << element.y << " " << element.z << std::endl;

    int localElmtId    = connManager->computeLocalElementIndex(element);

    TEST_ASSERT(localElmtId>=0);

    std::string eblock    = getElementBlock(element,*connManager);

    std::vector<Ordinal64> gids;
    dofManager->getElementGIDs(localElmtId,   gids);

    std::set<Ordinal64> s_gids;
    s_gids.insert(gids.begin(),gids.end());
 
    // ensure that the expected number of GIDs are produced
    TEST_EQUALITY(s_gids.size(),gids.size());
  }

  std::vector<Ordinal64> indices;
  dofManager->getOwnedIndices(indices);
  std::set<Ordinal64> s_indices;
  s_indices.insert(indices.begin(),indices.end());
  TEST_EQUALITY(s_indices.size(),indices.size()); // these should be the same

  int count = Teuchos::as<int>(s_indices.size());
  int totalCount = 0;
  Teuchos::reduceAll(comm,Teuchos::REDUCE_SUM,1,&count,&totalCount);

  TEST_EQUALITY(totalCount,(nx*poly_U+1)*(ny*poly_U+1)+(nx*poly_P+1)*(ny*poly_P+1))

  std::vector<Ordinal64> ghosted_indices;
  dofManager->getOwnedAndGhostedIndices(ghosted_indices);
  std::set<Ordinal64> s_ghosted_indices;
  s_ghosted_indices.insert(ghosted_indices.begin(),ghosted_indices.end());

  TEST_ASSERT(s_ghosted_indices.size()>=s_indices.size()); // should have more ghosted indices then owned indices
  TEST_EQUALITY(s_ghosted_indices.size(),ghosted_indices.size()); // these should be the same

}

TEUCHOS_UNIT_TEST(tCartesianDOFMgr_HighOrder, gid_values)
{
  typedef CartesianConnManager<int,Ordinal64> CCM;
  typedef panzer::DOFManager<int,Ordinal64> DOFManager;

  // build global (or serial communicator)
  #ifdef HAVE_MPI
    Teuchos::MpiComm<int> comm(MPI_COMM_WORLD);
  #else
    THIS_REALLY_DOES_NOT_WORK
  #endif

  int np = comm.getSize(); // number of processors
  //int rank = comm.getRank(); // processor rank

  // mesh description
  Ordinal64 nx = 8, ny = 4;//, nz = 4;
  //Ordinal64 nx = 4, ny = 3;//, nz = 4;
  int px = np, py = 1;//, pz = 1; // npx1 processor grids
  int bx =  1, by = 2;//, bz = 1; // 1x2 blocks

  // const int poly_U = 4, poly_P = 1, poly_T = 3;
  const int poly_U = 1;

  RCP<const panzer::FieldPattern> pattern_U = buildFieldPattern( rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>( poly_U )) );

  // build the topology
  RCP<CCM> connManager = rcp(new CCM);
  connManager->initialize(comm,nx,ny,px,py,bx,by);

  // build the dof manager, and assocaite with the topology
  RCP<DOFManager> dofManager = rcp(new DOFManager);
  dofManager->setConnManager(connManager,*comm.getRawMpiComm());

  // add velocity (U) and PRESSURE fields to the MHD element block
  dofManager->addField("eblock-0_0","UX",pattern_U);
  dofManager->addField("eblock-0_0","P",pattern_U);
  // dofManager->addField("eblock-0_1","UX",pattern_U);
  dofManager->addField("eblock-0_1","P",pattern_U);

  // build global unknowns (useful comment!)
  dofManager->buildGlobalUnknowns();

  auto myOffset   = connManager->getMyOffsetTriplet();
  auto myElements = connManager->getMyElementsTriplet();

  // check sharing locally on this processor
  {
    // choose an element in the middle of the day
    Triplet element;
    element.x = myOffset.x + myElements.x/2-1;
    element.y = myOffset.y + myElements.y/2-1;
    element.z = 0; 

    out << "Root element = " << element.x << " " << element.y << " " << element.z << std::endl;

    int localElmtId    = connManager->computeLocalElementIndex(element);
    int localElmtId_px = connManager->computeLocalElementIndex(Triplet(element.x+1,element.y,element.z));
    int localElmtId_py = connManager->computeLocalElementIndex(Triplet(element.x,element.y+1,element.z));

    TEST_ASSERT(localElmtId>=0);
    TEST_ASSERT(localElmtId_px>=0);
    TEST_ASSERT(localElmtId_py>=0);

    std::string eblock    = getElementBlock(element,*connManager);
    std::string eblock_px = getElementBlock(Triplet(element.x+1,element.y,element.z),*connManager);
    std::string eblock_py = getElementBlock(Triplet(element.x,element.y+1,element.z),*connManager);

    std::vector<Ordinal64> gids, gids_px, gids_py;

    dofManager->getElementGIDs(localElmtId,   gids);
    dofManager->getElementGIDs(localElmtId_px,gids_px);
    dofManager->getElementGIDs(localElmtId_py,gids_py);

    out << "gids = ";
    for(std::size_t i=0;i<gids.size();i++)
      out << gids[i] << " ";
    out << std::endl;

    out << "gids_px = ";
    for(std::size_t i=0;i<gids_px.size();i++)
      out << gids_px[i] << " ";
    out << std::endl;

    out << "gids_py = ";
    for(std::size_t i=0;i<gids_py.size();i++)
      out << gids_py[i] << " ";
    out << std::endl;

    std::sort(gids.begin(),gids.end());
    std::sort(gids_px.begin(),gids_px.end());
    std::sort(gids_py.begin(),gids_py.end());

    TEST_ASSERT(gids[0]>=0);
    TEST_ASSERT(gids_px[0]>=0);
    TEST_ASSERT(gids_py[0]>=0);
  }
}

TEUCHOS_UNIT_TEST(tCartesianDOFMgr_HighOrder, quad2d)
{
  typedef CartesianConnManager<int,Ordinal64> CCM;
  typedef panzer::DOFManager<int,Ordinal64> DOFManager;

  // build global (or serial communicator)
  #ifdef HAVE_MPI
    Teuchos::MpiComm<int> comm(MPI_COMM_WORLD);
  #else
    THIS_REALLY_DOES_NOT_WORK
  #endif

  int np = comm.getSize(); // number of processors
  int rank = comm.getRank(); // processor rank

  // mesh description
  Ordinal64 nx = 10, ny = 7;//, nz = 4;
  //Ordinal64 nx = 4, ny = 3;//, nz = 4;
  int px = np, py = 1;//, pz = 1; // npx1 processor grids
  int bx =  1, by = 2;//, bz = 1; // 1x2 blocks

  // build velocity, temperature and pressure fields
  const int poly_U = 4, poly_P = 1, poly_T = 3;

  RCP<const panzer::FieldPattern> pattern_U = buildFieldPattern( rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>( poly_U )) );
  RCP<const panzer::FieldPattern> pattern_P = buildFieldPattern( rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>( poly_P )) );
  RCP<const panzer::FieldPattern> pattern_T = buildFieldPattern( rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>( poly_T )) );
  
  // build the topology
  RCP<CCM> connManager = rcp(new CCM);
  connManager->initialize(comm,nx,ny,px,py,bx,by);

  // build the dof manager, and assocaite with the topology
  RCP<DOFManager> dofManager = rcp(new DOFManager);
  dofManager->setConnManager(connManager,*comm.getRawMpiComm());

  // add TEMPERATURE field to all element blocks (MHD and solid)
  dofManager->addField("TEMPERATURE",pattern_T);

  // add velocity (U) and PRESSURE fields to the MHD element block
  dofManager->addField("eblock-0_0","UX",pattern_U);
  dofManager->addField("eblock-0_0","UY",pattern_U);
  dofManager->addField("eblock-0_0","PRESSURE",pattern_P);

  // add velocity (U) fields to the solid element block
  dofManager->addField("eblock-0_1","UX",pattern_U);
  dofManager->addField("eblock-0_1","UY",pattern_U);

  // int temp_num = dofManager->getFieldNum("TEMPERATURE");
  int p_num   = dofManager->getFieldNum("PRESSURE");
  int t_num   = dofManager->getFieldNum("TEMPERATURE");
  int ux_num   = dofManager->getFieldNum("UX");
  int uy_num   = dofManager->getFieldNum("UY");

  out << "\n\nP, T, UX, UY = " << p_num << ", " << t_num << ", " 
                           << ux_num << ", " << uy_num << std::endl;

  // build global unknowns (useful comment!)
  dofManager->buildGlobalUnknowns();

  TEST_EQUALITY(dofManager->getElementBlockGIDCount("eblock-0_0"),
                pattern_U->numberIds() + 
                pattern_U->numberIds() + 
                pattern_P->numberIds() +
                pattern_T->numberIds() );

  TEST_EQUALITY(dofManager->getElementBlockGIDCount("eblock-0_1"),
                pattern_U->numberIds() + 
                pattern_U->numberIds() + 
                pattern_T->numberIds() );
  
  // print out some diagnostic information 
  ///////////////////////////////////////////////////////////

  dofManager->printFieldInformation(out); 

  out << std::endl << "Load balancing: " << printUGILoadBalancingInformation(*dofManager) << std::endl;

  out << std::endl << "Mesh Topology: " << std::endl;
  printMeshTopology(out,*dofManager);

  auto myOffset   = connManager->getMyOffsetTriplet();
  auto myElements = connManager->getMyElementsTriplet();

  out << "My Offset   = " << myOffset.x << " " << myOffset.y << " " << myOffset.z << std::endl;
  out << "My myElements = " << myElements.x << " " << myElements.y << " " << myElements.z << std::endl;

  // check sharing locally on this processor
  {
    // choose an element in the middle of the day
    Triplet element;
    element.x = myOffset.x + myElements.x/2-1;
    element.y = myOffset.y + myElements.y/2-1;
    element.z = 0; 

    out << "Root element = " << element.x << " " << element.y << " " << element.z << std::endl;

    int localElmtId    = connManager->computeLocalElementIndex(element);
    int localElmtId_px = connManager->computeLocalElementIndex(Triplet(element.x+1,element.y,element.z));
    int localElmtId_py = connManager->computeLocalElementIndex(Triplet(element.x,element.y+1,element.z));

    TEST_ASSERT(localElmtId>=0);
    TEST_ASSERT(localElmtId_px>=0);
    TEST_ASSERT(localElmtId_py>=0);

    std::string eblock    = getElementBlock(element,*connManager);
    std::string eblock_px = getElementBlock(Triplet(element.x+1,element.y,element.z),*connManager);
    std::string eblock_py = getElementBlock(Triplet(element.x,element.y+1,element.z),*connManager);

    std::vector<Ordinal64> gids, gids_px, gids_py;

    dofManager->getElementGIDs(localElmtId,   gids);
    dofManager->getElementGIDs(localElmtId_px,gids_px);
    dofManager->getElementGIDs(localElmtId_py,gids_py);

    // check that GIDs are all unique within an element
    {
      std::set<Ordinal64> s_gids, s_gids_px, s_gids_py;
      s_gids.insert(gids.begin(),gids.end());
      s_gids_px.insert(gids_px.begin(),gids_px.end());
      s_gids_py.insert(gids_py.begin(),gids_py.end());
 
      // ensure that the expected number of GIDs are produced
      TEST_EQUALITY(s_gids.size(),gids.size());
      TEST_EQUALITY(s_gids_px.size(),gids_px.size());
      TEST_EQUALITY(s_gids_py.size(),gids_py.size());
    }

    {
      out << "Elements " << localElmtId << " " << localElmtId_px << std::endl;
      auto offsets   = dofManager->getGIDFieldOffsets_closure(eblock,   ux_num,1,1); // +x
      auto offsets_n = dofManager->getGIDFieldOffsets_closure(eblock_px,ux_num,1,3); // -x

      TEST_EQUALITY(offsets.first.size(),offsets_n.first.size());

      std::vector<Ordinal64> gid_sub, gid_sub_px;
      for(std::size_t i=0;i<offsets.first.size();i++) {
        gid_sub.push_back(gids[offsets.first[i]]);
        gid_sub_px.push_back(gids_px[offsets_n.first[i]]);
      }

      std::sort(gid_sub.begin(),gid_sub.end());
      std::sort(gid_sub_px.begin(),gid_sub_px.end());

      // index comparison
      for(std::size_t i=0;i<gid_sub.size();i++)
        TEST_EQUALITY(gid_sub[i],gid_sub_px[i]);
    }

    {
      out << "Elements " << localElmtId << " " << localElmtId_py << std::endl;
      auto offsets   = dofManager->getGIDFieldOffsets_closure(   eblock,ux_num,1,2); // +y
      auto offsets_n = dofManager->getGIDFieldOffsets_closure(eblock_py,ux_num,1,0); // -y

      TEST_EQUALITY(offsets.first.size(),offsets_n.first.size());

      std::vector<Ordinal64> gid_sub, gid_sub_py;
      for(std::size_t i=0;i<offsets.first.size();i++) {
        gid_sub.push_back(gids[offsets.first[i]]);
        gid_sub_py.push_back(gids_py[offsets_n.first[i]]);
      }

      std::sort(gid_sub.begin(),gid_sub.end());
      std::sort(gid_sub_py.begin(),gid_sub_py.end());

      for(std::size_t i=0;i<gid_sub.size();i++)
        TEST_EQUALITY(gid_sub[i],gid_sub_py[i]);
    }
  }

  // assuming a 1d partition, check shared boundaries between processors
  {
    Triplet element_l;
    element_l.x = myOffset.x;
    element_l.y = myOffset.y + myElements.y/2;
    element_l.z = 0;

    out << "left element = " << element_l.x << " " << element_l.y << " " << element_l.z << std::endl;

    Triplet element_r;
    element_r.x = myOffset.x + myElements.x-1;
    element_r.y = myOffset.y + myElements.y/2;
    element_r.z = 0;

    out << "right element = " << element_r.x << " " << element_r.y << " " << element_r.z << std::endl;

    int localElmtId_l    = connManager->computeLocalElementIndex(element_l);
    int localElmtId_r    = connManager->computeLocalElementIndex(element_r);

    TEST_ASSERT(localElmtId_l>=0);
    TEST_ASSERT(localElmtId_r>=0);

    std::vector<Ordinal64> gids_l, gids_r;
    dofManager->getElementGIDs(   localElmtId_l,   gids_l);
    dofManager->getElementGIDs(   localElmtId_r,   gids_r);

    // check that GIDs are all unique within an element
    {
      std::set<Ordinal64> s_gids_l, s_gids_r;
      s_gids_l.insert(gids_l.begin(),gids_l.end());
      s_gids_r.insert(gids_r.begin(),gids_r.end());
 
      // ensure that the expected number of GIDs are produced
      TEST_EQUALITY(s_gids_l.size(),gids_l.size());
      TEST_EQUALITY(s_gids_r.size(),gids_r.size());
    }

    std::string eblock_l = getElementBlock(element_l,*connManager);
    std::string eblock_r = getElementBlock(element_r,*connManager);

    auto offsets_l   = dofManager->getGIDFieldOffsets_closure(   eblock_l,uy_num,1,3); // -x
    auto offsets_r   = dofManager->getGIDFieldOffsets_closure(   eblock_r,uy_num,1,1); // +x

    TEST_EQUALITY(offsets_l.first.size(),offsets_r.first.size());

    out << "Elements L/R " << localElmtId_l << " " << localElmtId_r << std::endl;
    std::vector<Ordinal64> gid_sub_l, gid_sub_r;
    for(std::size_t i=0;i<offsets_l.first.size();i++) {
      gid_sub_l.push_back(gids_l[offsets_l.first[i]]);
      gid_sub_r.push_back(gids_r[offsets_r.first[i]]);
    }

    std::sort(gid_sub_l.begin(),gid_sub_l.end());
    std::sort(gid_sub_r.begin(),gid_sub_r.end());

    // send left
    if(rank!=0) {
      Teuchos::send(comm,Teuchos::as<int>(gid_sub_l.size()),&gid_sub_l[0],rank-1);
    }

    // recieve right, check 
    if(rank!=np-1) {
      std::vector<Ordinal64> gid_remote(gid_sub_r.size(),-1);
      Teuchos::receive(comm,rank+1,Teuchos::as<int>(gid_sub_r.size()),&gid_remote[0]);

      for(std::size_t i=0;i<gid_sub_r.size();i++)
        TEST_EQUALITY(gid_sub_r[i],gid_remote[i]);
    }
  }
    
}

} // end unit test
} // end panzer

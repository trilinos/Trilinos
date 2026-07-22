// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#include "Intrepid2_HGRAD_HEX_Cn_FEM.hpp"
#include "Intrepid2_HCURL_HEX_In_FEM.hpp"
#include "Intrepid2_HDIV_HEX_In_FEM.hpp"

#include "PanzerCore_config.hpp"

#include "Panzer_ConnManager.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_GlobalIndexer_Utilities.hpp"

#include "CartesianConnManager.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace std {
  template<typename T>
  std::ostream& operator<<(ostream& os, const std::vector<T>& v)
  {
    os << "(";
    for (std::size_t i=0; i < v.size(); ++i) {
      if (i !=0)
        os << ",";
      os << v[i];
    }
    os << ")";
    return os;
  }
}

namespace panzer::unit_test {

using Triplet = CartesianConnManager::Triplet<panzer::GlobalOrdinal>;

std::string getElementBlock(
  const Triplet & element,
  const CartesianConnManager & connManager
) {
  int localElmtId = connManager.computeLocalBrickElementIndex(element);
  return connManager.getBlockId(localElmtId);
}

TEUCHOS_UNIT_TEST(tCartesianDOFMgr_DG, basic)
{
  Teuchos::MpiComm<int> comm(MPI_COMM_WORLD);    
  int np = comm.getSize();
  int rank = comm.getRank();

  const panzer::GlobalOrdinal nx = 12, ny = 5, nz = 3;
  const int px = np, py = 1, pz = 1;
  const int bx =  1, by = 1, bz = 1;

  // We need at least two elements in each direction on each process for
  // testing below.
  TEUCHOS_ASSERT((nx/np) >= 2); // parallel only in x dim
  TEUCHOS_ASSERT(ny >= 2);
  TEUCHOS_ASSERT(nz >= 2);

  // build the topology
  using CCM = CartesianConnManager;
  const auto connManager = Teuchos::make_rcp<CCM>();
  connManager->initialize(comm,nx,ny,nz,px,py,pz,bx,by,bz);

  // build the dof manager, and assocaite with the topology
  using DOFManager = panzer::DOFManager;
  const auto dofManager = Teuchos::make_rcp<DOFManager>();
  dofManager->setConnManager(connManager,*comm.getRawMpiComm());

  using Basis = Intrepid2::Basis<PHX::Device,double,double>;
  
  RCP<Basis> bhgrad2 = Teuchos::make_rcp<Intrepid2::Basis_HGRAD_HEX_Cn_FEM<PHX::Device, double, double>>(2);
  RCP<const FieldPattern> hgrad2 = Teuchos::make_rcp<Intrepid2FieldPattern>(bhgrad2);
  out << "HGRAD2\n" << *hgrad2 << std::endl;

  RCP<Basis> bhgrad1 = Teuchos::make_rcp<Intrepid2::Basis_HGRAD_HEX_Cn_FEM<PHX::Device, double, double>>(1);
  RCP<const FieldPattern> hgrad1 = Teuchos::make_rcp<Intrepid2FieldPattern>(bhgrad1);
  out << "HGRAD1\n" << *hgrad1 << std::endl;
  
  RCP<Basis> bhcurl = Teuchos::make_rcp<Intrepid2::Basis_HCURL_HEX_In_FEM<PHX::Device, double, double>>(1);
  RCP<const FieldPattern> hcurl = Teuchos::make_rcp<Intrepid2FieldPattern>(bhcurl);
  out << "HCURL2\n" << *hcurl << std::endl;
  
  RCP<Basis> bhdiv = Teuchos::make_rcp<Intrepid2::Basis_HDIV_HEX_In_FEM<PHX::Device, double, double>>(1);
  RCP<const FieldPattern> hdiv = Teuchos::make_rcp<Intrepid2FieldPattern>(bhdiv);
  out << "HDIV2\n" << *hdiv << std::endl;

  // Add fields. Use Hgrad2 for testing DG so that we have unknowns on
  // all subcells.
  dofManager->addField("DG",hgrad2,FieldType::DG);
  // the rest are CG. Testing includes mixing all FieldTypes
  dofManager->addField("UX",hgrad2);
  dofManager->addField("UY",hgrad2);
  dofManager->addField("UZ",hgrad2);
  dofManager->addField("P",hgrad1);
  dofManager->addField("T",hgrad2);
  dofManager->addField("E",hcurl);
  dofManager->addField("B",hdiv);

  // build global unknowns (useful comment!)
  dofManager->buildGlobalUnknowns();
 
  // print out some diagnostic information 
  ///////////////////////////////////////////////////////////

  out.setShowProcRank(true); // show process rank

  dofManager->printFieldInformation(out); 
  out << std::endl << "Load balancing: " << printUGILoadBalancingInformation(*dofManager) << std::endl;
  out << std::endl << "Mesh Topology: " << std::endl;
  printMeshTopology(out,*dofManager);

  // ***********************
  // Check that the total number of DOFs are consistent
  // ***********************
  {
    using ord_t = long long; // for MPI communication
    const ord_t numHgrad2FieldsCG = 4;
    const ord_t numHgrad1FieldsCG = 1;
    const ord_t numHcurlFieldsCG = 1;
    const ord_t numHdivFieldsCG = 1;
    const ord_t numHgrad2FieldsDG = 1;
    const ord_t numHgrad2UnknownsCG = (2*nx+1)*(2*ny+1)*(2*nz+1);
    const ord_t numHgrad1UnknownsCG = (nx+1)*(ny+1)*(nz+1);
    const ord_t numHcurlUnknownsCG = (nx+1) * ( (ny+1)*nz + (nz+1)*ny ) + (ny+1)*(nz+1)*nx;
    const ord_t numHdivUnknownsCG = (nx+1)*(ny*nz) + (ny+1)*(nx*nz) + (nz+1)*(nx*ny);
    const ord_t numCells = nx * ny * nz;
    const ord_t numHgrad2UnknownsDG = numCells * hgrad2->numberIds();
    
    const ord_t expectedTotalDOFs =
      numHgrad2FieldsCG * numHgrad2UnknownsCG 
      + numHgrad1FieldsCG * numHgrad1UnknownsCG
      + numHcurlFieldsCG * numHcurlUnknownsCG
      + numHdivFieldsCG * numHdivUnknownsCG
      + numHgrad2FieldsDG * numHgrad2UnknownsDG;
    
    ord_t myNumDOFs = (ord_t) dofManager->getNumOwned();
    ord_t numDOFs = 0;
    MPI_Reduce(&myNumDOFs,&numDOFs,1,MPI_LONG_LONG_INT,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Bcast(&numDOFs,1,MPI_LONG_LONG_INT,0,MPI_COMM_WORLD);
    TEST_EQUALITY(numDOFs, expectedTotalDOFs);
  }

  // ***********************
  // Check that field DOFs on the face of adjoining elements are the
  // same for CG but different for DG
  // ***********************
  auto myOffset   = connManager->getMyBrickOffsetTriplet();
  auto myNumElements = connManager->getMyBrickElementsTriplet();

  // out.setOutputToRootOnly(rank);
  out << "My Offset(" << myOffset.x << " " << myOffset.y << " " << myOffset.z << ")" << std::endl;
  out << "My myElements(" << myNumElements.x << "," << myNumElements.y << "," << myNumElements.z << ")"<< std::endl;
  // out.setOutputToRootOnly(0);

  int uxId = dofManager->getFieldNum("UX");
  int uyId = dofManager->getFieldNum("UY");
  int dgId = dofManager->getFieldNum("DG");

  {
    // choose an element in this MPI process workset that has at least one adjoining element local to this process.
    Triplet element;
    element.x = myOffset.x + myNumElements.x/2-1;
    element.y = myOffset.y + myNumElements.y/2-1;
    element.z = myOffset.z + myNumElements.z/2-1;

    out << "Root element = " << element.x << " " << element.y << " " << element.z << std::endl;

    int localElmtId    = connManager->computeLocalBrickElementIndex(element);
    int localElmtId_px = connManager->computeLocalBrickElementIndex(Triplet(element.x+1,element.y,element.z));
    int localElmtId_py = connManager->computeLocalBrickElementIndex(Triplet(element.x,element.y+1,element.z));
    int localElmtId_pz = connManager->computeLocalBrickElementIndex(Triplet(element.x,element.y,element.z+1));

    TEST_ASSERT(localElmtId>=0);
    TEST_ASSERT(localElmtId_px>=0);
    TEST_ASSERT(localElmtId_py>=0);
    TEST_ASSERT(localElmtId_pz>=0);

    std::string eblock    = getElementBlock(element,*connManager);
    std::string eblock_px = getElementBlock(Triplet(element.x+1,element.y,element.z),*connManager);
    std::string eblock_py = getElementBlock(Triplet(element.x,element.y+1,element.z),*connManager);
    std::string eblock_pz = getElementBlock(Triplet(element.x,element.y,element.z+1),*connManager);

    std::vector<panzer::GlobalOrdinal> gids, gids_px, gids_py, gids_pz;

    dofManager->getElementGIDs(   localElmtId,   gids);
    dofManager->getElementGIDs(localElmtId_px,gids_px);
    dofManager->getElementGIDs(localElmtId_py,gids_py);
    dofManager->getElementGIDs(localElmtId_pz,gids_pz);

    // CG in x
    {
      out << "\nCG in x" << std::endl; 
      out << "Elements " << localElmtId << " " << localElmtId_px << std::endl;
      auto offsets   = dofManager->getGIDFieldOffsets_closure(   eblock,uxId,2,1); // +x
      auto offsets_n = dofManager->getGIDFieldOffsets_closure(eblock_px,uxId,2,3); // -x

      TEST_EQUALITY(offsets.first.size(),offsets_n.first.size());
      TEST_ASSERT(offsets.first.size() > 0);

      std::vector<panzer::GlobalOrdinal> gid_sub, gid_sub_px;
      for(std::size_t i=0;i<offsets.first.size();i++) {
        gid_sub.push_back(gids[offsets.first[i]]);
        gid_sub_px.push_back(gids_px[offsets_n.first[i]]);
      }

      std::sort(gid_sub.begin(),gid_sub.end());
      std::sort(gid_sub_px.begin(),gid_sub_px.end());

      for(std::size_t i=0;i<gid_sub.size();i++)
        TEST_EQUALITY(gid_sub[i],gid_sub_px[i]);
    }

    // DG in x (CG matched, use the exact same procedure, but make sure DG is different)
    {
      out << "\nDG in x" << std::endl; 
      out << "Elements " << localElmtId << " " << localElmtId_px << std::endl;
      auto offsets   = dofManager->getGIDFieldOffsets_closure(   eblock,dgId,2,1); // +x
      auto offsets_n = dofManager->getGIDFieldOffsets_closure(eblock_px,dgId,2,3); // -x

      TEST_EQUALITY(offsets.first.size(),offsets_n.first.size());

      std::vector<panzer::GlobalOrdinal> gid_sub, gid_sub_px;
      for(std::size_t i=0;i<offsets.first.size();i++) {
        gid_sub.push_back(gids[offsets.first[i]]);
        gid_sub_px.push_back(gids_px[offsets_n.first[i]]);
      }

      std::sort(gid_sub.begin(),gid_sub.end());
      std::sort(gid_sub_px.begin(),gid_sub_px.end());

      for(std::size_t i=0;i<gid_sub.size();i++) {
        TEST_INEQUALITY(gid_sub[i],gid_sub_px[i]);
      }
      
      // Now check that the closure gids are in the volume gid list
      const auto& offsets_vol = dofManager->getGIDFieldOffsets(eblock,dgId);
      std::vector<panzer::GlobalOrdinal> gid_vol;
      for(std::size_t i=0;i<offsets_vol.size();i++)
        gid_vol.push_back(gids[offsets_vol[i]]);
      TEST_EQUALITY(gid_vol.size(), 27);
      TEST_EQUALITY(gid_sub.size(), 9);
      for (auto gid : gid_sub) {
        auto search = std::find(gid_vol.cbegin(),gid_vol.cend(),gid);
        TEST_ASSERT(search != gid_vol.cend());
      }
    }

    // CG in y
    {
      out << "\nCG in y" << std::endl; 
      out << "Elements " << localElmtId << " " << localElmtId_py << std::endl;
      auto offsets   = dofManager->getGIDFieldOffsets_closure(   eblock,uxId,2,2); // +y
      auto offsets_n = dofManager->getGIDFieldOffsets_closure(eblock_py,uxId,2,0); // -y

      TEST_EQUALITY(offsets.first.size(),offsets_n.first.size());
      TEST_ASSERT(offsets.first.size() > 0);

      std::vector<panzer::GlobalOrdinal> gid_sub, gid_sub_py;
      for(std::size_t i=0;i<offsets.first.size();i++) {
        gid_sub.push_back(gids[offsets.first[i]]);
        gid_sub_py.push_back(gids_py[offsets_n.first[i]]);
      }

      std::sort(gid_sub.begin(),gid_sub.end());
      std::sort(gid_sub_py.begin(),gid_sub_py.end());

      for(std::size_t i=0;i<gid_sub.size();i++)
        TEST_EQUALITY(gid_sub[i],gid_sub_py[i]);
    }

    // DG in y (CG matched, use the exact same procedure, but make sure DG is different)
    {
      out << "\nDG in y" << std::endl; 
      out << "Elements " << localElmtId << " " << localElmtId_py << std::endl;
      auto offsets   = dofManager->getGIDFieldOffsets_closure(   eblock,dgId,2,2); // +y
      auto offsets_n = dofManager->getGIDFieldOffsets_closure(eblock_py,dgId,2,0); // -y

      TEST_EQUALITY(offsets.first.size(),offsets_n.first.size());
      TEST_ASSERT(offsets.first.size() > 0);

      std::vector<panzer::GlobalOrdinal> gid_sub, gid_sub_py;
      for(std::size_t i=0;i<offsets.first.size();i++) {
        gid_sub.push_back(gids[offsets.first[i]]);
        gid_sub_py.push_back(gids_py[offsets_n.first[i]]);
      }

      std::sort(gid_sub.begin(),gid_sub.end());
      std::sort(gid_sub_py.begin(),gid_sub_py.end());

      for(std::size_t i=0;i<gid_sub.size();i++)
        TEST_INEQUALITY(gid_sub[i],gid_sub_py[i]);

      // Now check that the closure gids are in the volume gid list
      const auto& offsets_vol = dofManager->getGIDFieldOffsets(eblock,dgId);
      std::vector<panzer::GlobalOrdinal> gid_vol;
      for(std::size_t i=0;i<offsets_vol.size();i++)
        gid_vol.push_back(gids[offsets_vol[i]]);
      TEST_EQUALITY(gid_vol.size(), 27);
      TEST_EQUALITY(gid_sub.size(), 9);
      for (auto gid : gid_sub) {
        auto search = std::find(gid_vol.cbegin(),gid_vol.cend(),gid);
        TEST_ASSERT(search != gid_vol.cend());
      }
    }

    // CG in z
    {
      out << "\nCG in z" << std::endl; 
      out << "Elements " << localElmtId << " " << localElmtId_pz << std::endl;
      auto offsets   = dofManager->getGIDFieldOffsets_closure(   eblock,uxId,2,5); // +z
      auto offsets_n = dofManager->getGIDFieldOffsets_closure(eblock_pz,uxId,2,4); // -z

      TEST_EQUALITY(offsets.first.size(),offsets_n.first.size());
      TEST_ASSERT(offsets.first.size() > 0);

      std::vector<panzer::GlobalOrdinal> gid_sub, gid_sub_pz;
      for(std::size_t i=0;i<offsets.first.size();i++) {
        gid_sub.push_back(gids[offsets.first[i]]);
        gid_sub_pz.push_back(gids_pz[offsets_n.first[i]]);
      }

      std::sort(gid_sub.begin(),gid_sub.end());
      std::sort(gid_sub_pz.begin(),gid_sub_pz.end());

      for(std::size_t i=0;i<gid_sub.size();i++)
        TEST_EQUALITY(gid_sub[i],gid_sub_pz[i]);
    }

    // DG in z (CG matched, use the exact same procedure, but make
    // sure DG is different)
    {
      out << "\nDG in z" << std::endl; 
      out << "Elements " << localElmtId << " " << localElmtId_pz << std::endl;
      auto offsets   = dofManager->getGIDFieldOffsets_closure(   eblock,dgId,2,5); // +z
      auto offsets_n = dofManager->getGIDFieldOffsets_closure(eblock_pz,dgId,2,4); // -z

      TEST_EQUALITY(offsets.first.size(),offsets_n.first.size());
      TEST_ASSERT(offsets.first.size() > 0);

      std::vector<panzer::GlobalOrdinal> gid_sub, gid_sub_pz;
      for(std::size_t i=0;i<offsets.first.size();i++) {
        gid_sub.push_back(gids[offsets.first[i]]);
        gid_sub_pz.push_back(gids_pz[offsets_n.first[i]]);
      }

      std::sort(gid_sub.begin(),gid_sub.end());
      std::sort(gid_sub_pz.begin(),gid_sub_pz.end());

      for(std::size_t i=0;i<gid_sub.size();i++)
        TEST_INEQUALITY(gid_sub[i],gid_sub_pz[i]);

      // Now check that the closure gids are in the volume gid list
      const auto& offsets_vol = dofManager->getGIDFieldOffsets(eblock,dgId);
      std::vector<panzer::GlobalOrdinal> gid_vol;
      for(std::size_t i=0;i<offsets_vol.size();i++)
        gid_vol.push_back(gids[offsets_vol[i]]);
      TEST_EQUALITY(gid_vol.size(), 27);
      TEST_EQUALITY(gid_sub.size(), 9);
      for (auto gid : gid_sub) {
        auto search = std::find(gid_vol.cbegin(),gid_vol.cend(),gid);
        TEST_ASSERT(search != gid_vol.cend());
      }
    }

  }

  // assuming a 1d MPI partition, check shared boundaries between processors
  {
    Triplet element_l;
    element_l.x = myOffset.x;
    element_l.y = myOffset.y + myNumElements.y/2;
    element_l.z = myOffset.z + myNumElements.z/2;

    Triplet element_r;
    element_r.x = myOffset.x + myNumElements.x-1;
    element_r.y = myOffset.y + myNumElements.y/2;
    element_r.z = myOffset.z + myNumElements.z/2;

    int localElmtId_l    = connManager->computeLocalBrickElementIndex(element_l);
    int localElmtId_r    = connManager->computeLocalBrickElementIndex(element_r);

    TEST_ASSERT(localElmtId_l>=0);
    TEST_ASSERT(localElmtId_r>=0);

    std::vector<panzer::GlobalOrdinal> gids_l, gids_r;
    dofManager->getElementGIDs(   localElmtId_l,   gids_l);
    dofManager->getElementGIDs(   localElmtId_r,   gids_r);

    std::string eblock_l = getElementBlock(element_l,*connManager);
    std::string eblock_r = getElementBlock(element_r,*connManager);

    auto offsets_l   = dofManager->getGIDFieldOffsets_closure(   eblock_l,uyId,2,3); // -x
    auto offsets_r   = dofManager->getGIDFieldOffsets_closure(   eblock_r,uyId,2,1); // +x

    TEST_EQUALITY(offsets_l.first.size(),offsets_r.first.size());

    out << "Elements L/R " << localElmtId_l << " " << localElmtId_r << std::endl;
    std::vector<panzer::GlobalOrdinal> gid_sub_l, gid_sub_r;
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
      std::vector<panzer::GlobalOrdinal> gid_remote(gid_sub_r.size(),-1);
      Teuchos::receive(comm,rank+1,Teuchos::as<int>(gid_sub_r.size()),&gid_remote[0]);

      for(std::size_t i=0;i<gid_sub_r.size();i++)
        TEST_EQUALITY(gid_sub_r[i],gid_remote[i]);
    }
  }
  
}

} // namespace panzer::unit test

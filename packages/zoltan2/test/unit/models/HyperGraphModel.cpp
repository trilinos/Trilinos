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
// Testing of HyperGraphModel built from APF mesh adapters.
//

/*! \brief Test of HyperGraphModel interface.
 *
 *  \todo test all methods of HyperGraphModel
 *  \todo add functionality to choose action from command line
 */

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_HyperGraphModel.hpp>
#include <Zoltan2_APFMeshAdapter.hpp>
#include <Zoltan2_Environment.hpp>

// Teuchos includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// SCOREC includes
#ifdef HAVE_ZOLTAN2_PARMA
#include <parma.h>
#include <apf.h>
#include <apfMesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <gmi_mesh.h>
#endif

using Teuchos::ParameterList;
using Teuchos::RCP;

int main(int narg, char *arg[]) {

  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > CommT = Tpetra::getDefaultComm();

#ifdef HAVE_ZOLTAN2_PARMA
  //Setup for SCOREC
  PCU_Comm_Init();

  // Generate mesh with MDS
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh("../partition/pumiTri14/plate.dmg","../partition/pumiTri14/2/");

  typedef Zoltan2::APFMeshAdapter<apf::Mesh2*> inputAdapter_t;
  typedef Zoltan2::MeshAdapter<apf::Mesh2*> baseMeshAdapter_t;
  Teuchos::ParameterList params("test params");
  params.set("timer_output_stream" , "std::cout");
  params.set("debug_level", "verbose_detailed_status");
  params.set("hypergraph_model_type","ghosting");

  RCP<Zoltan2::Environment> env;
  try{
    env = rcp(new Zoltan2::Environment(params, Tpetra::getDefaultComm()));
  }
  Z2_FORWARD_EXCEPTIONS

  RCP<const Zoltan2::Environment> envConst = Teuchos::rcp_const_cast<const Zoltan2::Environment>(env);


  inputAdapter_t* ia = new inputAdapter_t(*CommT, m,"vertex","edge",false);
  inputAdapter_t::scalar_t* arr = new inputAdapter_t::scalar_t[ia->getLocalNumOf(ia->getPrimaryEntityType())];
  for (size_t i=0;i<ia->getLocalNumOf(ia->getPrimaryEntityType());i++) {
    arr[i]=PCU_Comm_Self();
  }
  const inputAdapter_t::scalar_t* weights=arr;
  ia->setWeights(ia->getPrimaryEntityType(),weights,1);

  const baseMeshAdapter_t *base_ia = dynamic_cast<const baseMeshAdapter_t*>(ia);
  Zoltan2::modelFlag_t graphFlags_;
  RCP<const baseMeshAdapter_t> baseInputAdapter_(base_ia,false);

  Zoltan2::HyperGraphModel<inputAdapter_t> model(baseInputAdapter_,envConst,CommT,
                                                 graphFlags_,Zoltan2::HYPEREDGE_CENTRIC);
  ia->destroy();
  delete ia;

  //Delete APF Mesh;
  m->destroyNative();
  apf::destroyMesh(m);
  //End communications
  PCU_Comm_Free();

#endif

  std::cout<<"PASS\n";
  return 0;
}

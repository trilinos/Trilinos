// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "PanzerSTK_UnitTest_BuildMesh.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"

#include "Panzer_STK_Interface.hpp"

#include "PanzerSTK_UnitTest_STKInterfaceGenerator.hpp"

namespace panzer_stk
{

Teuchos::RCP<panzer_stk::STK_Interface>
buildMesh(const std::vector<int> & N,
          const std::vector<int> & B,
          const std::vector<double> & L)
{

  const int num_dims = N.size();

  TEUCHOS_ASSERT(L.size() == static_cast<std::size_t>(num_dims));

  Teuchos::ParameterList mesh_parameters("Mesh");

  mesh_parameters.set<double>("X0",-L[0]/2.);
  mesh_parameters.set<double>("Xf", L[0]/2.);
  mesh_parameters.set<int>("X Elements", N[0]);
  mesh_parameters.set<int>("X Blocks", B[0]);
  if(num_dims>1){
    mesh_parameters.set<double>("Y0",-L[1]/2.);
    mesh_parameters.set<double>("Yf", L[1]/2.);
    mesh_parameters.set<int>("Y Elements", N[1]);
    mesh_parameters.set<int>("Y Blocks", B[1]);
  }
  if(num_dims>2){
    mesh_parameters.set<double>("Z0",-L[2]/2.);
    mesh_parameters.set<double>("Zf", L[2]/2.);
    mesh_parameters.set<int>("Z Elements", N[2]);
    mesh_parameters.set<int>("Z Blocks", B[2]);
  }

  if(num_dims == 1){
    mesh_parameters.set<std::string>("Mesh Type","Line");
  } else if(num_dims == 2){
    mesh_parameters.set<std::string>("Mesh Type","Quad");
  } else if(num_dims == 3){
    mesh_parameters.set<std::string>("Mesh Type","Hex");
  } else {
    TEUCHOS_ASSERT(false);
  }

  return generateMesh(mesh_parameters);
}


Teuchos::RCP<panzer_stk::STK_Interface>
buildParallelMesh(const std::vector<int> & N,
                  const std::vector<int> & B,
                  const std::vector<int> & P,
                  const std::vector<double> & L,
                  const std::vector<int> & periodic_dims)
{

  const size_t num_dims = N.size();

  TEUCHOS_ASSERT(L.size() == num_dims);
  TEUCHOS_ASSERT(P.size() == num_dims);

  Teuchos::ParameterList mesh_parameters("Mesh");

  mesh_parameters.set<double>("X0",-L[0]/2.);
  mesh_parameters.set<double>("Xf", L[0]/2.);
  mesh_parameters.set<int>("X Elements", N[0]);
  mesh_parameters.set<int>("X Blocks", B[0]);

  // The line mesh factory doesn't allow this
  if(num_dims>1)
    mesh_parameters.set<int>("X Procs", P[0]);

  if(num_dims>1){
    mesh_parameters.set<double>("Y0",-L[1]/2.);
    mesh_parameters.set<double>("Yf", L[1]/2.);
    mesh_parameters.set<int>("Y Elements", N[1]);
    mesh_parameters.set<int>("Y Blocks", B[1]);
    mesh_parameters.set<int>("Y Procs", P[1]);
  }
  if(num_dims>2){
    mesh_parameters.set<double>("Z0",-L[2]/2.);
    mesh_parameters.set<double>("Zf", L[2]/2.);
    mesh_parameters.set<int>("Z Elements", N[2]);
    mesh_parameters.set<int>("Z Blocks", B[2]);
    mesh_parameters.set<int>("Z Procs", P[2]);
  }

  if(num_dims == 1){
    mesh_parameters.set<std::string>("Mesh Type","Line");
  } else if(num_dims == 2){
    mesh_parameters.set<std::string>("Mesh Type","Quad");
  } else if(num_dims == 3){
    mesh_parameters.set<std::string>("Mesh Type","Hex");
  } else {
    TEUCHOS_ASSERT(false);
  }

  if(periodic_dims.size() > 0){
    auto & periodic_params = mesh_parameters.sublist("Periodic BCs");
    int arg=1;
    for(const auto & periodic_dim : periodic_dims){
      if(periodic_dim == 0){
        if(num_dims==1)
          periodic_params.set("Periodic Condition " +std::to_string(arg++), "y-coord 1.e-8: left;right");
        else if(num_dims==2){
          periodic_params.set("Periodic Condition " +std::to_string(arg++), "y-coord 1.e-8: left;right");
          periodic_params.set("Periodic Condition " +std::to_string(arg++), "y-edge 1.e-8: left;right");
        } else if(num_dims==3){
          periodic_params.set("Periodic Condition " +std::to_string(arg++), "yz-coord 1.e-8: left;right");
          periodic_params.set("Periodic Condition " +std::to_string(arg++), "yz-edge 1.e-8: left;right");
          periodic_params.set("Periodic Condition " +std::to_string(arg++), "yz-face 1.e-8: left;right");
        }
      } else if(periodic_dim == 1){
        if(num_dims==2){
          periodic_params.set("Periodic Condition " +std::to_string(arg++), "x-coord 1.e-8: top;bottom");
          periodic_params.set("Periodic Condition " +std::to_string(arg++), "x-edge 1.e-8: top;bottom");
        } else if(num_dims==3){
          periodic_params.set("Periodic Condition " +std::to_string(arg++), "xz-coord 1.e-8: top;bottom");
          periodic_params.set("Periodic Condition " +std::to_string(arg++), "xz-edge 1.e-8: top;bottom");
          periodic_params.set("Periodic Condition " +std::to_string(arg++), "xz-face 1.e-8: top;bottom");
        }
      } else if(periodic_dim == 2){
        if(num_dims==3){
          periodic_params.set("Periodic Condition " +std::to_string(arg++), "xy-coord 1.e-8: front;back");
          periodic_params.set("Periodic Condition " +std::to_string(arg++), "xy-edge 1.e-8: front;back");
          periodic_params.set("Periodic Condition " +std::to_string(arg++), "xy-face 1.e-8: front;back");
        }
      }
    }
    periodic_params.set("Count",arg-1);
  }

  return generateMesh(mesh_parameters);
}

}

// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.



/* ========================================================================================== */
/*  This is generated code, do not edit - see SmootherMetricGen.mhpp and *.m, *.nb files      */
/* ========================================================================================== */


#ifndef ScalingMatricesGen_hpp
#define ScalingMatricesGen_hpp

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#ifdef __GNUC__
# if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wsign-compare"
#endif // GCC_VERSION
#endif // __GNUC__

#ifdef __INTEL_COMPILER
#pragma warning disable 1599
#pragma warning disable 1478
#endif // __INTEL_COMPILER

#include <percept/mesh/mod/smoother/SmootherMetric.hpp>
#include <percept/mesh/mod/smoother/MeshSmoother.hpp>

namespace percept {

  static const int indices_wedge[8][4] = {{0, 1, 2, 3}, {1, 2, 0, 4},
                                       {2, 0, 1, 5}, {3, 5, 4, 0},
                                       {4, 3, 5, 1}, {5, 4, 3, 2}};
  static const int indices_hex[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5},
                                     {2, 3, 1, 6}, {3, 0, 2, 7},
                                     {4, 7, 5, 0}, {5, 4, 6, 1},
                                     {6, 5, 7, 2}, {7, 6, 4, 3}};

  const int indices_tri[8][4] = {{0,1,2,0},{0,1,2,0},{0,1,2,0},{0,1,2,0}};  // last slot just to keep all the code below easy and each corner has 3 connections
  const int indices_tet[8][4] = {{0,1,2,3},{0,1,2,3},{0,1,2,3},{0,1,2,3}};  // jacobian is the same for all nodes

  static int indices_pyr[4][4];

  class ScalingMatrices {
  public:
    static ScalingMatrices s_scalingMatrices;
    enum types {
      Tri,
      Quad,
      Hex,
      Tet,
      Pyr,
      Wedge,
      END_elem_types
    };

    MDArray SM[stk::topology::END_TOPOLOGY+1];

    ScalingMatrices() {
      Teuchos::Array<int> dim33(3,3);
      static double sTet[3][3] = {{1, -isqrt3, -isqrt6}, {0, 2*isqrt3, -isqrt6}, {0, 0, 3*isqrt6}} ;
      SM[stk::topology::TET_4] = MDArray(dim33, &sTet[0][0], true, true);
      static double sTri[3][3] = {{1, -isqrt3, 0}, {0, 2*isqrt3, 0}, {0, 0, 1}} ;
      SM[stk::topology::TRI_3] = MDArray(dim33, &sTri[0][0], true, true);
      SM[stk::topology::TRI_3_2D] = MDArray(dim33, &sTri[0][0], true, true);
      static double sPyr[3][3] = {{1, 0, -0.5}, {0, 1, -0.5}, {0, 0, 1.}} ;
      SM[stk::topology::PYRAMID_5] = MDArray(dim33, &sPyr[0][0], true, true);
      static double sQuad[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}} ;
      SM[stk::topology::QUAD_4] = MDArray(dim33, &sQuad[0][0], true, true);
      SM[stk::topology::QUAD_4_2D] = MDArray(dim33, &sQuad[0][0], true, true);
      static double sHex[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}} ;
      SM[stk::topology::HEX_8] = MDArray(dim33, &sHex[0][0], true, true);
      static double sWedge[3][3] = {{1, -isqrt3, 0}, {0, 2*isqrt3, 0}, {0, 0, 1}} ;
      SM[stk::topology::WEDGE_6] = MDArray(dim33, &sWedge[0][0], true, true);
    }

    inline const MDArray& get(stk::topology::topology_t type)
    {
      return SM[type];
    }
  };

  class Indices {
    const int *s_allIndices[stk::topology::END_TOPOLOGY+1][8];
  public:
    static Indices s_indices;
    Indices() {
      for (int i = 0; i < 4; ++i) {
        indices_pyr[i][0] = i;
        indices_pyr[i][1] = (i+1)%4;
        indices_pyr[i][2] = (i+3)%4;
        indices_pyr[i][3] = 4;
      }
      for (int i = 0; i < 8; ++i) {
        s_allIndices[stk::topology::TRI_3][i] = indices_tri[i];
        s_allIndices[stk::topology::TRI_3_2D][i] = indices_tri[i];
        s_allIndices[stk::topology::QUAD_4][i] = indices_hex[i];
        s_allIndices[stk::topology::QUAD_4_2D][i] = indices_hex[i];
        s_allIndices[stk::topology::HEX_8][i] = indices_hex[i];
        s_allIndices[stk::topology::TET_4][i] = indices_tet[i];
        s_allIndices[stk::topology::WEDGE_6][i] = indices_wedge[i];
        s_allIndices[stk::topology::PYRAMID_5][i] = indices_pyr[i];
      }
    }

    inline const int *get_indices(stk::topology::topology_t type, int i)
    {
#ifndef NDEBUG
      switch(type) {
      case stk::topology::TRI_3:
      case stk::topology::TRI_3_2D:
      case stk::topology::QUAD_4:
      case stk::topology::QUAD_4_2D:
      case stk::topology::HEX_8:
      case stk::topology::TET_4:
      case stk::topology::WEDGE_6:
      case stk::topology::PYRAMID_5:
        break;
      default:
        {
          stk::topology tt;
          tt.m_value = type;
          std::cout << "topology not handled= " << tt << std::endl;
          throw std::runtime_error("topology not handled = " + toString(type));
        }
      }
#endif
      return s_allIndices[type][i];
    }
  };


}


#ifdef __GNUC__
# if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406
#pragma GCC diagnostic pop
#endif // GCC_VERSION
#endif // __GNUC__

#endif
#endif

// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef MeshSmoother_hpp
#define MeshSmoother_hpp

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#include <percept/PerceptMesh.hpp>
#include <percept/MeshType.hpp>

#if defined(STK_PERCEPT_HAS_GEOMETRY)
#include <percept/mesh/geometry/kernel/MeshGeometry.hpp>
#endif

#undef USE_CALLGRIND_MESH_SMOOTHER
//#define USE_CALLGRIND_MESH_SMOOTHER
#ifdef USE_CALLGRIND_MESH_SMOOTHER
#include "/usr/netpub/valgrind-3.8.1/include/valgrind/callgrind.h"
#endif
#include <percept/mesh/mod/smoother/GenericAlgorithm_total_element_metric.hpp>


  namespace percept {

#if !defined(STK_PERCEPT_HAS_GEOMETRY)
    class MeshGeometry {};
#endif

    /// Abstract base class smoother
    template<typename MeshType>
    class MeshSmootherImpl
    {

    protected:
      PerceptMesh *m_eMesh;
      int innerIter;
      double gradNorm;
      int parallelIterations;
      STKMesh::MTSelector *m_stk_boundarySelector;
    public:
      typename MeshType::MTMeshGeometry *m_meshGeometry;

    public:

      MeshSmootherImpl(PerceptMesh *eMesh,
                   STKMesh::MTSelector *stk_select=0,
                   typename MeshType::MTMeshGeometry *meshGeometry=0,
                       int innerIter=100, double gradNorm = 1.e-8, int parallelIterations=20);

      STKMesh::MTSelector * get_stkmesh_select() const {return m_stk_boundarySelector;}

      virtual ~MeshSmootherImpl() {}

      static size_t parallel_count_invalid_elements(PerceptMesh *eMesh);

      void run();
      virtual void run_algorithm() = 0;

      static bool select_bucket(typename MeshType::MTBucket& bucket, PerceptMesh *eMesh);
      std::pair<bool,int> get_fixed_flag(typename MeshType::MTNode node_ptr);
      int classify_node(typename MeshType::MTNode node, GeometryHandle curveOrSurfaceEvaluator /*size_t& curveOrSurfaceEvaluator*/) const;
      void project_delta_to_tangent_plane(typename MeshType::MTNode node, double *delta, double *norm=0);
      void enforce_tangent_plane(typename MeshType::MTNode node, double rhs[3], double lhs[3][3], double *norm=0);

      /// if reset is true, don't actually modify the node's coordinates and only
      /// return the snapped position in @param coordinate
      void snap_to(typename MeshType::MTNode node_ptr,  double *coordinate, bool reset=false) const;
    };

    using MeshSmoother = MeshSmootherImpl<STKMesh>;
  }


#endif
#endif

// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef GeometryKernelGregoryPatch_hpp
#define GeometryKernelGregoryPatch_hpp

#include <percept/PerceptMesh.hpp>
#include "GeometryKernel.hpp"
#include "xfer/GPSTKMeshTransferSetup.hpp"

namespace percept {
class GeometryKernelGregoryPatch : public GeometryKernel
{
public:
  typedef std::set<stk::mesh::Entity> EntitySet;
  typedef std::vector<stk::mesh::Entity> EntityVector;
  typedef std::map<stk::mesh::Entity, EntitySet > EntityToEntitySetMap;

  GeometryKernelGregoryPatch(percept::PerceptMesh& meshWithNodes, bool debug=false) :
    m_nodeMesh(meshWithNodes), m_geometryMesh(0), m_geometryMeshActiveParts(), m_debug(debug)
  {}
  virtual ~GeometryKernelGregoryPatch();

  virtual bool read_file(const std::string& file_name,
                         std::vector<GeometryHandle>& geometry_entities) override;

  virtual void pre_process(percept::PerceptMesh *eMeshWithNode, const stk::mesh::PartVector& parts) override;

  virtual bool debug_dump_file(const std::string& file_name) override;

  virtual void snap_to(KernelPoint& point, GeometryHandle geom,
                       double *converged_tolerance = NULL,
                       double *uvw_computed = NULL,
                       double *uvw_hint = NULL,
                       void *extra_hint = NULL) override;
  virtual void normal_at(KernelPoint& point, GeometryHandle geom, std::vector<double>& normal, void *extra_hint = NULL) override;

  bool in_face(stk::mesh::Entity face, double *uv) const;

  percept::PerceptMesh *get_mesh_with_nodes() { return &m_nodeMesh; }
  percept::PerceptMesh *get_mesh_with_geometry() { return m_geometryMesh; }

private:

  // given a set of faces (@param neighbors), and an input point, find the closest
  // valid point (valid = point when projected is in the parameter space of the face)
  bool findClosestPointInternal(const double *point,
                                std::set<stk::mesh::Entity>& neighbors,
                                stk::mesh::Entity node_hint,
                                double *closest_point,
                                double *uvw_computed,
                                bool allowError,  // ignore errors and report closest point
                                bool debug);

  percept::PerceptMesh& m_nodeMesh;
  percept::PerceptMesh *m_geometryMesh; // the mesh that holds the geometry
  stk::mesh::PartVector m_geometryMeshActiveParts;
  stk::mesh::PartVector m_nodeMeshActiveParts;
  bool m_debug;
  std::shared_ptr<percept::GPSTKMeshTransfer> m_meshTransfer;
  stk::mesh::Entity m_found_face;

};

}

#endif // GeometryKernelGregoryPatch_hpp

// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef MeshType_hpp_
#define MeshType_hpp_

#include <percept/Percept.hpp>

#include <array>
#include <memory>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>

#include <Shards_CellTopologyData.h>
#include <Shards_CellTopology.hpp>
#include <percept/function/MDArray.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

namespace percept {


enum NodeClassifyType {
    MS_VERTEX,
    MS_CURVE,
    MS_SURFACE,
    MS_VOLUME,
    MS_ON_BOUNDARY,
    MS_NOT_ON_BOUNDARY
};

  typedef long double Double;

  class PerceptMesh;

  class MeshGeometry;

  struct STKMesh {
    typedef stk::mesh::Selector MTSelector;
    typedef MeshGeometry MTMeshGeometry;
    typedef stk::mesh::Entity MTNode;
    typedef stk::mesh::Entity MTElement;
    typedef stk::mesh::Bucket MTBucket;
    typedef stk::mesh::FieldBase MTField;
    typedef CellTopologyData MTCellTopology;
    enum {NELEM_TYPES = 10 };
  };

  template<typename MeshType>
  unsigned get_num_nodes(PerceptMesh *eMesh, typename MeshType::MTElement element);

  template<typename MeshType>
  const typename MeshType::MTNode *get_nodes(PerceptMesh *eMesh, typename MeshType::MTElement element, std::vector<typename MeshType::MTNode> *nodes );

  template<typename MeshType>
  bool MTisGhostNode(PerceptMesh *m_eMesh, typename MeshType::MTNode node);

  template<typename MeshType>
  bool MTnode_locally_owned(PerceptMesh *m_eMesh, typename MeshType::MTNode node);

  template<typename MeshType>
  void MTcommFields(std::vector<const typename MeshType::MTField*>& fields, PerceptMesh *m_eMesh);

  template<typename MeshType>
  void MTsum_fields(std::vector<const typename MeshType::MTField*>& fields, PerceptMesh *m_eMesh);

  /// gets @param field data from @param node into @param fld
  template<typename MeshType>
  void get_field(double *fld, unsigned size, PerceptMesh *eMesh, typename MeshType::MTField *field, typename MeshType::MTNode node);

  template<typename MeshType>
  void get_field_new(double *fld, unsigned size, PerceptMesh *eMesh, typename MeshType::MTField *field, typename MeshType::MTNode node);

  /// gets @param field data from @param node into @param fld[@param index]
  template<typename MeshType>
  void get_field(double *fld, unsigned size, int index, PerceptMesh *eMesh, typename MeshType::MTField *field, typename MeshType::MTNode node);

  /// sets @param field data from @param fld into @param node
  template<typename MeshType>
  void set_field(const double * fld, unsigned size, PerceptMesh *eMesh, typename MeshType::MTField *field, typename MeshType::MTNode node);

  /// sets @param field data from @param fld[@param index] into @param node
  template<typename MeshType>
  void set_field(const double * fld, unsigned size, int index, PerceptMesh *eMesh, typename MeshType::MTField *field, typename MeshType::MTNode node);

}
#endif

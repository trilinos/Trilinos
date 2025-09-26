// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/PerceptMesh.hpp>
#include <percept/MeshType.hpp>

namespace percept {

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // specializations
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template<>
  unsigned get_num_nodes<STKMesh>(PerceptMesh *eMesh, typename STKMesh::MTElement element)
  {
    return eMesh->get_bulk_data()->num_nodes(element);
  }

  template<>
  const typename STKMesh::MTNode * get_nodes<STKMesh>(PerceptMesh *eMesh, typename STKMesh::MTElement element, std::vector<typename STKMesh::MTNode > */*nodes*/ )
  {
    return eMesh->get_bulk_data()->begin_nodes(element);
  }

  template<>
  bool MTisGhostNode<STKMesh>(PerceptMesh *m_eMesh, typename STKMesh::MTNode node)
  {
    return m_eMesh->aura(node);
  }

  template<>
  bool MTnode_locally_owned<STKMesh>(PerceptMesh *m_eMesh, typename STKMesh::MTNode node)
  {
    return (m_eMesh->get_rank() == m_eMesh->owner_rank(node));
  }

  template<>
  void MTcommFields<STKMesh>(std::vector<const typename STKMesh::MTField*>& fields, PerceptMesh *m_eMesh)
  {
    stk::mesh::copy_owned_to_shared(*m_eMesh->get_bulk_data(), fields);
    stk::mesh::communicate_field_data(m_eMesh->get_bulk_data()->aura_ghosting(), fields);
  }

  template<>
  void MTsum_fields<STKMesh>(std::vector<const typename STKMesh::MTField*>& fields, PerceptMesh *m_eMesh)
  {
    stk::mesh::parallel_sum(*m_eMesh->get_bulk_data(), fields);
  }

  template<>
  void get_field<STKMesh>(double * fld, unsigned size, PerceptMesh *eMesh, typename STKMesh::MTField *field, typename STKMesh::MTNode node)
  {
    double *data = eMesh->field_data(field, node);
    for (unsigned i=0; i < size; ++i)
      {
        fld[i] = data[i];
      }
  }

  /// gets @param field data from @param node into @param fld[@param index]
  template<>
  void get_field<STKMesh>(double *fld, unsigned /*size*/, int index, PerceptMesh *eMesh, typename STKMesh::MTField *field, typename STKMesh::MTNode node)
  {
    double *data = eMesh->field_data(field, node);
    fld[index] = data[index];
  }

  /// sets @param field data from @param fld into @param node
  template<>
  void set_field<STKMesh>(const double *fld, unsigned size, PerceptMesh *eMesh, typename STKMesh::MTField *field, typename STKMesh::MTNode node)
  {
    double *data = eMesh->field_data(field, node);
    for (unsigned i=0; i < size; ++i)
      {
        data[i] = fld[i];
      }
  }

  /// sets @param field data from @param fld[@param index] into @param node
  template<>
  void set_field<STKMesh>(const double *fld, unsigned /*size*/, int index, PerceptMesh *eMesh, typename STKMesh::MTField *field, typename STKMesh::MTNode node)
  {
    double *data = eMesh->field_data(field, node);
    data[index] = fld[index];
  }
}

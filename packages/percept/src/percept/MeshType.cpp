// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/PerceptMesh.hpp>
#include <percept/MeshType.hpp>
#include <percept/structured/BlockStructuredGrid.hpp>

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
  const typename STKMesh::MTNode * get_nodes<STKMesh>(PerceptMesh *eMesh, typename STKMesh::MTElement element, std::vector<typename STKMesh::MTNode > *nodes )
  {
    return eMesh->get_bulk_data()->begin_nodes(element);
  }

  template<>
  unsigned get_num_nodes<StructuredGrid>(PerceptMesh *eMesh, typename StructuredGrid::MTElement element)
  {
    return 8;
  }

  template<>
  const typename StructuredGrid::MTNode * get_nodes<StructuredGrid>(PerceptMesh *eMesh, typename StructuredGrid::MTElement element, std::vector<typename StructuredGrid::MTNode> *nodes )
  {
    VERIFY_OP_ON(nodes, !=, 0, "bad nodes");
    nodes->resize(0);
    for (unsigned i2 = element[2]; i2 <= element[2]+1; ++i2)
      {
        for (unsigned i1 = element[1]; i1 <= element[1]+1; ++i1)
          {
            for (unsigned i0 = element[0]; i0 <= element[0]+1; ++i0)
              {
                nodes->push_back(StructuredCellIndex{{i0,i1,i2,element[3]}});
              }
          }
      }
    VERIFY_OP_ON(nodes->size(), ==, 8, "bad nodes size");
    return &(*nodes)[0];
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
  bool MTisGhostNode<StructuredGrid>(PerceptMesh *m_eMesh, typename StructuredGrid::MTNode node)
  {
    //VERIFY_MSG("not impl");
    //FIXME
    return false;
  }

  template<>
  bool MTnode_locally_owned<StructuredGrid>(PerceptMesh *m_eMesh, typename StructuredGrid::MTNode node)
  {
    //VERIFY_MSG("not impl");
    //FIXME
    return true;
  }

  template<>
  void MTcommFields<STKMesh>(std::vector<const typename STKMesh::MTField*>& fields, PerceptMesh *m_eMesh)
  {
    stk::mesh::copy_owned_to_shared(*m_eMesh->get_bulk_data(), fields);
    stk::mesh::communicate_field_data(m_eMesh->get_bulk_data()->aura_ghosting(), fields);
  }

  template<>
  void MTcommFields<StructuredGrid>(std::vector<const typename StructuredGrid::MTField*>& fields, PerceptMesh *m_eMesh)
  {
    m_eMesh->get_block_structured_grid()->comm_fields(fields);
  }

  template<>
  void MTsum_fields<STKMesh>(std::vector<const typename STKMesh::MTField*>& fields, PerceptMesh *m_eMesh)
  {
    stk::mesh::parallel_sum(*m_eMesh->get_bulk_data(), fields);
  }

  template<>
  void MTsum_fields<StructuredGrid>(std::vector<const typename StructuredGrid::MTField*>& fields, PerceptMesh *m_eMesh)
  {
    m_eMesh->get_block_structured_grid()->sum_fields(fields);
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

  template<>
  void get_field<StructuredGrid>(double *fld, unsigned size, PerceptMesh *eMesh, typename StructuredGrid::MTField *field, typename StructuredGrid::MTNode node)
  {
    unsigned iblock = node[3];
    std::shared_ptr<StructuredBlock> sgrid = eMesh->get_block_structured_grid()->m_sblocks[iblock];
    const unsigned A0 = sgrid->m_access_ordering[0], A1 = sgrid->m_access_ordering[1], A2 = sgrid->m_access_ordering[2];

    for (unsigned i=0; i < size; ++i)
      {
        fld[i] = (*field->m_block_fields[iblock])(node[A0], node[A1], node[A2], i);
      }
  }

  /// gets @param field data from @param node into @param fld[@param index]
  template<>
  void get_field<STKMesh>(double *fld, unsigned size, int index, PerceptMesh *eMesh, typename STKMesh::MTField *field, typename STKMesh::MTNode node)
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
  void set_field<STKMesh>(const double *fld, unsigned size, int index, PerceptMesh *eMesh, typename STKMesh::MTField *field, typename STKMesh::MTNode node)
  {
    double *data = eMesh->field_data(field, node);
    data[index] = fld[index];
  }

  template<>
  void set_field<StructuredGrid>(const double *fld, unsigned size, PerceptMesh *eMesh, typename StructuredGrid::MTField *field, typename StructuredGrid::MTNode node)
  {
    unsigned iblock = node[3];
    std::shared_ptr<StructuredBlock> sgrid = eMesh->get_block_structured_grid()->m_sblocks[iblock];
    const unsigned A0 = sgrid->m_access_ordering[0], A1 = sgrid->m_access_ordering[1], A2 = sgrid->m_access_ordering[2];

    for (unsigned i=0; i < size; ++i)
      {
        (*field->m_block_fields[iblock])(node[A0], node[A1], node[A2], i) = fld[i];
      }
  }

  template<>
  void set_field<StructuredGrid>(const double *fld, unsigned size, int index, PerceptMesh *eMesh, typename StructuredGrid::MTField *field, typename StructuredGrid::MTNode node)
  {
    unsigned iblock = node[3];
    std::shared_ptr<StructuredBlock> sgrid = eMesh->get_block_structured_grid()->m_sblocks[iblock];
    const unsigned A0 = sgrid->m_access_ordering[0], A1 = sgrid->m_access_ordering[1], A2 = sgrid->m_access_ordering[2];

    (*field->m_block_fields[iblock])(node[A0], node[A1], node[A2], index) = fld[index];

  }

}

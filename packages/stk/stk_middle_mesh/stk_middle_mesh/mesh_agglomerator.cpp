#include "mesh_agglomerator.hpp"
#include "stk_util/util/ReportHandler.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

std::vector<int> MeshAgglomerator::get_group_idxs(SetType<MeshEntityPtr>& vertsIn, const int nthres)
{
  std::vector<int> idxs;
  for (int i = 0; i < get_num_groups(); ++i)
  {
    int nfound   = 0;
    int nvertsRemaining = vertsIn.size();
    for (auto& v : vertsIn)
    {
      if (contains_entity_sorted(m_verts[i], v))
      {
        nfound += 1;
      }

      if (nfound == nthres)
      {
        idxs.push_back(i);
        break;
      }

      nvertsRemaining--;
      if (nfound + nvertsRemaining < nthres)
      {
        break;
      }
    }
  }

  return idxs;
}

void MeshAgglomerator::get_verts(std::shared_ptr<Mesh> meshIn, const VectorType<MeshEntityPtr>& els, VectorType<MeshEntityPtr>& vertVector)
{
  using Bool = int_least8_t;
  FieldPtr<Bool> seenEntitiesFieldPtr = create_field<Bool>(meshIn, FieldShape(1, 0, 0), 1, false);
  auto& seenEntitiesField = *seenEntitiesFieldPtr;

  std::array<MeshEntityPtr, MAX_DOWN> elVerts;
  for (auto& el : els)
  {
    int nverts = get_downward(el, 0, elVerts.data());
    for (int i = 0; i < nverts; ++i)
    {
      MeshEntityPtr vert = elVerts[i];
      if (!seenEntitiesField(vert, 0, 0))
      {
        vertVector.push_back(vert);
        seenEntitiesField(vert, 0, 0) = true;
      }
    }
  }  
}

bool MeshAgglomerator::contains_entity_sorted(const VectorType<MeshEntityPtr>& entities, MeshEntityPtr entity)
{
  STK_ThrowAssertMsg(std::is_sorted(entities.begin(), entities.end(), MeshEntityCompare()), "Vector must be sorted by pointer value");
  return std::binary_search(entities.begin(), entities.end(), entity, MeshEntityCompare());
}


std::ostream& operator<<(std::ostream& os, const MeshAgglomerator& agg)
{
  os << "MeshAgglomerator with " << agg.get_num_groups() << " groups:" << std::endl;
  for (int i = 0; i < agg.get_num_groups(); ++i)
  {
    os << "  group " << i << " verts: ";
    const auto& verts = agg.get_group_verts(i);
    for (auto& v : verts)
      os << v->get_id() << " ";
    os << std::endl;
  }

  return os;
}

} // namespace impl
} // namespace mesh
} // namespace middle_mesh
} // namespace stk

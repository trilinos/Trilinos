#ifndef stk_mesh_DefaultFEM_hpp
#define stk_mesh_DefaultFEM_hpp

#include <stdexcept>
#include <vector>
#include <map>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/fem/FEMInterface.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

namespace stk {
namespace mesh {
// namespace fem {

class DefaultFEM : public fem::FEMInterface
{
public:
  typedef std::map<fem::CellTopology, std::pair<Part *, EntityRank> > CellTopologyPartEntityRankMap;
  typedef std::vector<fem::CellTopology> PartCellTopologyVector;

  DefaultFEM(MetaData &meta_data);
  
  DefaultFEM(MetaData &meta_data, size_t spatial_dimension);

  virtual size_t get_spatial_dimension() const {
    return m_spatialDimension;
  }

  virtual void set_spatial_dimension(size_t spatial_dimension);

  virtual void register_cell_topology(const fem::CellTopology cell_topology, EntityRank entity_rank);
  
  virtual void set_cell_topology(Part &part, fem::CellTopology cell_topology);


  virtual fem::CellTopology get_cell_topology(const Part &part) const;

  virtual EntityRank get_entity_rank(const fem::CellTopology cell_topology) const;
  
  virtual Part &get_part(const fem::CellTopology cell_topology) const;

private:
  void initialize(size_t spatial_dimension);

private:
  MetaData &                    m_metaData;
  size_t                        m_spatialDimension;
  CellTopologyPartEntityRankMap m_cellTopologyPartEntityRankMap;
  PartCellTopologyVector        m_partCellTopologyVector;
};

} // namespace mesh
} // namespace stk

#endif // stk_mesh_DefaultFEM_hpp

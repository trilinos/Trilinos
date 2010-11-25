#ifndef stk_mesh_DefaultFEM_hpp
#define stk_mesh_DefaultFEM_hpp

#include <stdexcept>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/fem/FEMInterface.hpp>

namespace stk {
namespace mesh {
// namespace fem {

class DefaultFEM : public fem::FEMInterface
{
public:
  typedef std::vector< std::pair< fem::CellTopology , EntityRank > > TopologyEntityRankMap ;
  typedef std::vector< std::pair< unsigned , fem::CellTopology > > PartCellTopologyMap ;

  DefaultFEM(MetaData &meta_data);
  
  DefaultFEM(MetaData &meta_data, size_t spatial_dimension);

  virtual size_t get_spatial_dimension() const {
    return m_spatialDimension;
  }

  virtual void set_spatial_dimension(size_t spatial_dimension);

  
  virtual void set_cell_topology(const Part &part, fem::CellTopology cell_topology);

  virtual fem::CellTopology get_cell_topology(const Part &part) const;


  virtual void set_entity_rank(const fem::CellTopology cell_topology, EntityRank entity_rank);
  
  virtual EntityRank get_entity_rank(fem::CellTopology cell_topology) const;
  
private:
  void initialize(size_t spatial_dimension);

private:
  size_t                m_spatialDimension;
  TopologyEntityRankMap m_topEntityRank;
  PartCellTopologyMap   m_partCellTopologyMap;
};

namespace fem {

std::vector<std::string> entity_rank_names(size_t spatial_dimension);

} // namespace fem

} // namespace mesh
} // namespace stk

#endif // stk_mesh_DefaultFEM_hpp

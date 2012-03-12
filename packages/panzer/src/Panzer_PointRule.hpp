#ifndef PANZER_POINT_RULE_HPP
#define PANZER_POINT_RULE_HPP

#include "Teuchos_ArrayRCP.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Shards_CellTopology.hpp"

namespace panzer {

  class CellData;

  /** Base class useful for constructing data layouts
    * for points on a reference cell.
    */
  class PointRule {
  public:
    
    /** if side = -1 then we use the cell as an reference frame
      *
      * \param[in] ptName Name of the point rule.
      * \param[in] np Number of points per cell
      * \param[in] cell_data Description of the cell
      */
    PointRule(const std::string & ptName,int np, const panzer::CellData& cell_data);

    void setup(const std::string & ptName,int np, const panzer::CellData& cell_data);
  
    // Returns true if this point rule is for a sideset
    bool isSide();

    /** Get the name of this point rule.
      */
    const std::string & getName() const;

    Teuchos::RCP<const shards::CellTopology> topology;
    
    Teuchos::RCP<shards::CellTopology> side_topology;
    
    //! Data layout for scalar fields
    Teuchos::RCP<PHX::DataLayout> dl_scalar;
    //! Data layout for vector fields
    Teuchos::RCP<PHX::DataLayout> dl_vector;
    //! Data layout for rank-2 tensor fields
    Teuchos::RCP<PHX::DataLayout> dl_tensor;
    
    int spatial_dimension;
    int workset_size;
    int num_points;

    //! Defaults to -1 if this is volume and not sideset
    int side;

    //! print information about the integration rule
    virtual void print(std::ostream & os);
  
  protected:
    PointRule() : side(-1) {}

    /** Look up side topology for a cell_data object. Returns null if
      * cell data does not correspond to a side object.
      */
    static Teuchos::RCP<shards::CellTopology> getSideTopology(const CellData & cell_data);

  private:
    std::string point_name;
  };

}

#endif

// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_POINT_RULE_HPP
#define PANZER_POINT_RULE_HPP

#include "Teuchos_ArrayRCP.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Shards_CellTopology.hpp"

#include "Panzer_PointDescriptor.hpp"

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


    /** \brief
      *
      * \param[in] point_rule_name Name for point rule
      * \param[in] num_cells Number of cells
      * \param[in] num_points_per_cell Number of points in each cell
      * \param[in] num_face Number of faces
      * \param[in] num_points_per_face Number of points on each face
      * \param[in] cell_topology Cell topology
      */
    PointRule(const std::string & point_rule_name,
              const int num_cells,
              const int num_points_per_cell,
              const int num_faces,
              const int num_points_per_face,
              const Teuchos::RCP<const shards::CellTopology> & cell_topology);

    /** Constructor from a point description.
      */
    PointRule(const panzer::PointDescriptor& description,
              const Teuchos::RCP<const shards::CellTopology> & cell_topology,
              const int num_cells);

    //! Destructor (Satisfying the compiler)
    virtual ~PointRule() {}

    void setup(const std::string & ptName,int np, const panzer::CellData& cell_data);
  
    // Returns true if this point rule is for a sideset
    bool isSide() const;

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
    
    //! Data layout for vector fields - full (x,y,z)
    Teuchos::RCP<PHX::DataLayout> dl_vector3;

    //! Data layout for vector fields - full ((xx,xy,xz),(yx,yy,yz),(zx,zy,zz)) (or transpose?)
    Teuchos::RCP<PHX::DataLayout> dl_tensor3x3;

    int spatial_dimension;
    int workset_size;
    int num_points;

    //! Defaults to -1 if this is volume and not sideset
    int side;

    //! print information about the integration rule
    virtual void print(std::ostream & os);

    // TODO: These need to be moved to a DataLayoutGenerator

    Teuchos::RCP<PHX::DataLayout> getCellDataLayout() const;
    Teuchos::RCP<PHX::DataLayout> getCellDataLayout(const int dim0) const;
    Teuchos::RCP<PHX::DataLayout> getCellDataLayout(const int dim0, const int dim1) const;

    Teuchos::RCP<PHX::DataLayout> getCellPointDataLayout() const;
    Teuchos::RCP<PHX::DataLayout> getCellPointDataLayout(const int dim0) const;
    Teuchos::RCP<PHX::DataLayout> getCellPointDataLayout(const int dim0, const int dim1) const;

    Teuchos::RCP<PHX::DataLayout> getFaceDataLayout() const;
    Teuchos::RCP<PHX::DataLayout> getFaceDataLayout(const int dim0) const;
    Teuchos::RCP<PHX::DataLayout> getFaceDataLayout(const int dim0, const int dim1) const;

    Teuchos::RCP<PHX::DataLayout> getFacePointDataLayout() const;
    Teuchos::RCP<PHX::DataLayout> getFacePointDataLayout(const int dim0) const;
    Teuchos::RCP<PHX::DataLayout> getFacePointDataLayout(const int dim0, const int dim1) const;
  
  protected:
    PointRule() : side(-1) {}

    /** Look up side topology for a cell_data object. Returns null if
      * cell data does not correspond to a side object.
      */
    static Teuchos::RCP<shards::CellTopology> getSideTopology(const CellData & cell_data);


    void setup(const std::string & point_rule_name,
               const int num_cells,
               const int num_points_per_cell,
               const int num_faces,
               const int num_points_per_face,
               const Teuchos::RCP<const shards::CellTopology> & cell_topology);

    int _num_faces;
    int _num_points_per_face;

  private:
    std::string point_name;
  };

}

#endif

#include <stk_percept/Percept.hpp>

#include <stk_percept/function/internal/IsInElement.hpp>


#include <Intrepid_CellTools.hpp>
#include <Intrepid_FunctionSpaceTools.hpp>

#include <stk_percept/norm/IntrepidManager.hpp>

namespace stk
{
  namespace percept
  {

    IsInElement::IsInElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates) :
      m_found_it(false), m_input_phy_points(input_phy_points), m_found_parametric_coordinates(found_parametric_coordinates),
      m_foundElement(0)
    {}

    void IsInElement::init_elementOp()
    {
      m_found_it=false;
    }
    void IsInElement::fini_elementOp()
    {
    }

    bool IsInElement::operator()(const stk::mesh::Entity& element, stk::mesh::FieldBase *field,  const mesh::BulkData& bulkData)
    {

      unsigned found_it;
      isInElement(m_input_phy_points, m_found_parametric_coordinates, found_it, element, bulkData);
      //if (EXTRA_PRINT) std::cout << "IsInElement::operator() found_it = " << found_it << std::endl;
      // break out of enclosing element loop
      if (found_it)
        {
          m_found_it = true;
          m_foundElement = &element;
          return true;
        }
      else
        return false;
    }

    /**
     *  Dimensions of input_phy_points = ([P]=1, [D])
     *  Dimensions of found_parametric_coordinates = ([P]=1, [D])
     */
    void IsInElement::isInElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates, unsigned& found_it, const mesh::Entity& element,
                                  const mesh::BulkData& bulkData)
    {
      IntrepidManager::isInElement(input_phy_points, found_parametric_coordinates, found_it, element, bulkData);
#if 0
      found_it = 0;

      // FIXME consider caching the coords_field in FieldFunction
      const mesh::fem::FEMMetaData& metaData = stk::mesh::fem::FEMMetaData::get(bulkData);
      VectorFieldType *coords_field = metaData.get_field<VectorFieldType >("coordinates");

      const mesh::Bucket & bucket = element.bucket();
      const CellTopologyData * const bucket_cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(bucket);

      unsigned numCells = 1; // FIXME

      shards::CellTopology topo(bucket_cell_topo_data);
      unsigned numNodes = topo.getNodeCount();
      unsigned cellDim  = topo.getDimension();
      MDArray cellWorkset(numCells, numNodes, cellDim);

      /// FIXME -- fill cellWorkset
      const mesh::PairIterRelation elem_nodes = element.relations( mesh::Node );

      for (unsigned iCell = 0; iCell < numCells; iCell++)
        {
          for (unsigned iNode = 0; iNode < numNodes; iNode++)
            {
              mesh::Entity& node = *elem_nodes[iNode].entity();
              double * node_coord_data = stk::mesh::field_data( *coords_field , node);
              for (unsigned iDim=0; iDim < cellDim; iDim++)
                {
                  cellWorkset(iCell, iNode, iDim) = node_coord_data[iDim];
                }
            }
        }

      // FIXME for multiple points
      if (input_phy_points.rank() == 1)
        {
          VERIFY_1("IsInElement::isInElement bad rank of input_phy_points");
        }
      VERIFY_OP(input_phy_points.dimension(0), == , 1, "IsInElement::isInElement bad input_phy_points 1st dim");
      VERIFY_OP(input_phy_points.dimension(1), >= , (int)cellDim, "IsInElement::isInElement bad input_phy_points 2nd dim");

      if (found_parametric_coordinates.rank() == 1)
        {
          VERIFY_1("IsInElement::isInElement bad rank of found_parametric_coordinates");
        }
      VERIFY_OP(found_parametric_coordinates.dimension(0), == , 1, "IsInElement::isInElement bad found_parametric_coordinates 1st dim");
      VERIFY_OP(found_parametric_coordinates.dimension(1), == , (int)cellDim,
                "IsInElement::isInElement bad found_parametric_coordinates 2nd dim");

      unsigned cellOrd = 0;  // FIXME
      Intrepid::CellTools<double>::mapToReferenceFrame(found_parametric_coordinates, input_phy_points, cellWorkset, topo, cellOrd);
      MDArrayUInt inclusion_results(1);  // FIXME
      Intrepid::CellTools<double>::checkPointwiseInclusion(inclusion_results, found_parametric_coordinates, topo);
      found_it = inclusion_results(0);
      if (found_it)
        {
          // for testing only
          if (0)
            {
              FieldContainer<double> images(1, cellDim );
              //Intrepid::CellTools<double>::mapToPhysicalFrame(images, preImages, triNodes, triangle_3, whichCell);
              Intrepid::CellTools<double>::mapToPhysicalFrame(images, found_parametric_coordinates, cellWorkset, topo, cellOrd);
            }
        }
#endif
    }


  }
}

#ifndef stk_percept_BuildBoundingBoxes_hpp
#define stk_percept_BuildBoundingBoxes_hpp

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/diag/IdentProc.hpp>

#include <Shards_CellTopology.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_percept/function/ElementOp.hpp>

#define EXTRA_PRINT 0

namespace stk_classic
{
  namespace percept
  {

    typedef mesh::Field<double>                     ScalarFieldType ;
    typedef mesh::Field<double, mesh::Cartesian>    VectorFieldType ;

    template<unsigned SpatialDim>
    class BuildBoundingBoxes : public ElementOp
    {

    public:
      typedef stk_classic::search::ident::IdentProc<uint64_t,unsigned> IdentProc;
      typedef stk_classic::search::box::PointBoundingBox<IdentProc,double,SpatialDim> BoundingPoint;
      typedef stk_classic::search::box::AxisAlignedBoundingBox<IdentProc,double,SpatialDim> AABoundingBox;

      std::vector<AABoundingBox>& m_boxes;
      VectorFieldType *m_coords_field;
      bool m_notInitialized;
    public:
      BuildBoundingBoxes(std::vector<AABoundingBox>& boxes, VectorFieldType *coords_field) :  m_boxes(boxes), m_coords_field(coords_field),
                                                                                              m_notInitialized(false)
      {
      }

      void init_elementOp()
      {
      }
      void fini_elementOp()
      {
        m_notInitialized=true;  // force this object to be used only once 
      }
      bool operator()(const stk_classic::mesh::Entity& element, stk_classic::mesh::FieldBase *field,  const mesh::BulkData& bulkData);

      AABoundingBox getBoundingBox(const stk_classic::mesh::Entity& element, const mesh::BulkData& bulkData)
      {
        double bbox[2*SpatialDim];
        const mesh::PairIterRelation elem_nodes = element.relations( stk_classic::mesh::fem::FEMMetaData::NODE_RANK );
        unsigned numNodes = elem_nodes.size();
        for (unsigned iNode = 0; iNode < numNodes; iNode++)
          {
            mesh::Entity& node = *elem_nodes[iNode].entity();
            double * coord_data = mesh::field_data( *m_coords_field, node);
            if (iNode == 0)
              {
                for (unsigned iDim = 0; iDim < SpatialDim; iDim++)
                  {
                    bbox[iDim]              = coord_data[iDim];
                    bbox[iDim + SpatialDim] = coord_data[iDim];
                  }
              }
            else
              {
                for (unsigned iDim = 0; iDim < SpatialDim; iDim++)
                  {
                    bbox[iDim]              = std::min(bbox[iDim],              coord_data[iDim]);
                    bbox[iDim + SpatialDim] = std::max(bbox[iDim + SpatialDim], coord_data[iDim]);
                  }
              }
          }

        AABoundingBox bb;
        bb.key.ident = element.identifier();
        bb.set_box(bbox);

        return bb;
      }
    };


    // FIXME
    //template<unsigned SpatialDim>
    //std::ostream &operator<<(std::ostream &out, const typename BuildBoundingBoxes<SpatialDim>::BoundingBox &bbox)
    std::ostream &operator<<(std::ostream &out, const BuildBoundingBoxes<3>::AABoundingBox &bbox);

    //template<unsigned SpatialDim>
    //std::ostream &operator<<(std::ostream &out, const typename BuildBoundingBoxes<SpatialDim>::BoundingBox &bbox)
    std::ostream &operator<<(std::ostream &out, const BuildBoundingBoxes<3>::BoundingPoint &bbox);

    std::ostream &operator<<(std::ostream &out, const BuildBoundingBoxes<2>::AABoundingBox &bbox);

    //template<unsigned SpatialDim>
    //std::ostream &operator<<(std::ostream &out, const typename BuildBoundingBoxes<SpatialDim>::BoundingBox &bbox)
    std::ostream &operator<<(std::ostream &out, const BuildBoundingBoxes<2>::BoundingPoint &bbox);

#if 0
    template<unsigned SpatialDim>
    std::ostream &operator<<(std::ostream &out, const typename BuildBoundingBoxes<SpatialDim>::AABoundingBox &bbox)
    {
      out << "bbox_min = { " 
          << bbox.lower(0) << ",  "
          << bbox.lower(1) << ",  "
          << bbox.lower(2) << "}  "
          << "bbox_max = { " 
          << bbox.upper(0) << ",  "
          << bbox.upper(1) << ",  "
          << bbox.upper(2) << "}  ";
      return out;
    }
#endif

#if 0
    template<typename BBox_loc>
    std::ostream &operator<<(std::ostream &out, const BBox_loc &bbox)
    {
      out << "bbox_min = { " 
          << bbox.lower(0) << ",  "
          << bbox.lower(1) << ",  "
          << bbox.lower(2) << "}  "
          << "bbox_max = { " 
          << bbox.upper(0) << ",  "
          << bbox.upper(1) << ",  "
          << bbox.upper(2) << "}  ";
      return out;
    }
#endif

    template<unsigned SpatialDim>
    bool BuildBoundingBoxes<SpatialDim>::operator()(const stk_classic::mesh::Entity& element, stk_classic::mesh::FieldBase *field,  const mesh::BulkData& bulkData)
    {
      if (m_notInitialized)
        throw std::runtime_error("BuildBoundingBoxes::operator(): you must re-construct this object before reusing it");

      AABoundingBox bb = getBoundingBox(element, bulkData);
      if (0 || EXTRA_PRINT) std::cout << "bb = " << bb << std::endl;
      m_boxes.push_back(bb);
      return false;  // never break out of the enclosing loop
    }


  }
}

#endif

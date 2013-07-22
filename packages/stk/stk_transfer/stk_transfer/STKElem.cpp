

#include <limits>

#include <boost/shared_ptr.hpp>

#include <Intrepid_Utils.hpp>
#include <Intrepid_FunctionSpaceTools.hpp>
#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_CellTools.hpp>
#include <Intrepid_Basis.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

#include <stk_util/environment/ReportHandler.hpp>
#include <stk_search/BoundingBox.hpp>

namespace stk {
namespace transfer {
namespace STKElemUtil {

typedef Intrepid::FieldContainer<double>   MDArray;
typedef Intrepid::FieldContainer<unsigned> MDArrayUInt;

typedef Intrepid::Basis<double, MDArray> BasisType;
typedef boost::shared_ptr<BasisType>     BasisTypeRCP;
typedef std::map<unsigned,BasisTypeRCP>  BasisTable;

template<unsigned DIM> 
unsigned parametric(std::vector<double> &para_coords, 
                    const double *to, 
                    const mesh::Entity element,
                    const mesh::FieldBase &coords_field,
                    const mesh::BulkData& bulkData)
{   

  const unsigned dimension  = DIM; 
  const unsigned numCells = 1; // FIXME

  unsigned found_it = 0;

  MDArray input_phy_points (1,dimension);
  for (unsigned i=0; i<dimension; ++i) input_phy_points(0,i) = to[i];

  const mesh::Bucket & bucket = bulkData.bucket(element);
  const CellTopologyData * const bucket_cell_topo_data = mesh::get_cell_topology(bucket).getCellTopologyData();
  ThrowRequireMsg (bucket_cell_topo_data, __FILE__<<":"<<__LINE__<<" parametric::bogus topology");

  shards::CellTopology topo(bucket_cell_topo_data);
  const unsigned numNodes = topo.getNodeCount();

  mesh::Entity const* elem_node_rels = bulkData.begin_nodes(element);
  const unsigned num_nodes = bulkData.num_nodes(element);

  ThrowRequireMsg (topo.getDimension() == dimension,__FILE__<<":"<<__LINE__<<" Wrong spatical dimension"
    <<" for topology. Expected "<<dimension<<" found "<<topo.getDimension());
  ThrowRequireMsg (numNodes == num_nodes ,
    __FILE__<<":"<<__LINE__<<" Expected "<<numNodes<<" nodes but found "<<num_nodes);

  /// FIXME -- fill cellWorkset
  MDArray cellWorkset(numCells, numNodes, dimension);
  for (unsigned iCell = 0; iCell < numCells; iCell++) {   
    for (unsigned iNode = 0; iNode < numNodes; iNode++) {   
      const mesh::Entity node = elem_node_rels[iNode];
      const double * coords = static_cast<const double*>(bulkData.field_data(coords_field, node));
      for (unsigned iDim=0; iDim < dimension; iDim++) cellWorkset(iCell, iNode, iDim) = coords[iDim];
    }   
  }   

  MDArray parametric_coordinates(1,dimension);
  const unsigned cellOrd = 0;  // FIXME
  Intrepid::CellTools<double>::mapToReferenceFrame(parametric_coordinates, 
                                                   input_phy_points, 
                                                   cellWorkset, 
                                                   topo, 
                                                   cellOrd);
  MDArrayUInt inclusion_results(1);  // FIXME
  const double threshold = 1.e-4; // (INTREPID_THRESHOLD default = 10*double_eps ~ 20e-16)
  Intrepid::CellTools<double>::checkPointwiseInclusion(inclusion_results,
                                                       parametric_coordinates,
                                                       topo,
                                                       threshold);
  found_it = inclusion_results(0);
  if (found_it) {
    // for testing only
    if (0) {
      MDArray images(1, dimension );
      //Intrepid::CellTools<double>::mapToPhysicalFrame(images, preImages, triNodes, triangle_3, whichCell);
      Intrepid::CellTools<double>::mapToPhysicalFrame(images,
                                                      parametric_coordinates,
                                                      cellWorkset,
                                                      topo,
                                                      cellOrd);
    }
  }
  para_coords.resize(dimension);
  for (unsigned i=0; i<dimension; ++i) para_coords[i] = parametric_coordinates(0,i);
  return found_it;
}

template unsigned parametric<1> (std::vector<double> &para_coords, 
                     const double *to, 
                     const mesh::Entity element,
                     const mesh::FieldBase &coords_field,
                     const mesh::BulkData& bulkData);
template unsigned parametric<2> (std::vector<double> &para_coords, 
                     const double *to, 
                     const mesh::Entity element,
                     const mesh::FieldBase &coords_field,
                     const mesh::BulkData& bulkData);
template unsigned parametric<3> (std::vector<double> &para_coords, 
                     const double *to, 
                     const mesh::Entity element,
                     const mesh::FieldBase &coords_field,
                     const mesh::BulkData& bulkData);

BasisTable setupBasisTable() {
  BasisTable basisTable; 
  basisTable[shards::getCellTopologyData<shards::Line<2> >()-> key]               .reset( new Intrepid::Basis_HGRAD_LINE_C1_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::Triangle<3> >()-> key]           .reset( new Intrepid::Basis_HGRAD_TRI_C1_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::Triangle<6> >()-> key]           .reset( new Intrepid::Basis_HGRAD_TRI_C2_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::Quadrilateral<4> >()-> key]      .reset( new Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::Quadrilateral<9> >()-> key]      .reset( new Intrepid::Basis_HGRAD_QUAD_C2_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::Hexahedron<8> >()-> key]         .reset( new Intrepid::Basis_HGRAD_HEX_C1_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::Hexahedron<27> >()-> key]        .reset( new Intrepid::Basis_HGRAD_HEX_C2_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::Tetrahedron<4> >()-> key]        .reset( new Intrepid::Basis_HGRAD_TET_C1_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::Tetrahedron<10> >()-> key]       .reset( new Intrepid::Basis_HGRAD_TET_C2_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::Wedge<6> >()-> key]              .reset( new Intrepid::Basis_HGRAD_WEDGE_C1_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::ShellTriangle<3> >()-> key]      .reset( new Intrepid::Basis_HGRAD_TRI_C1_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::ShellTriangle<6> >()-> key]      .reset( new Intrepid::Basis_HGRAD_TRI_C2_FEM<double, MDArray >() );
  basisTable[shards::getCellTopologyData<shards::ShellQuadrilateral<4> >()-> key] .reset( new Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, MDArray >() );
  return basisTable;
}

const BasisTypeRCP getBasis(const shards::CellTopology& topo) {    
  static const BasisTable basisTable(setupBasisTable());

  const unsigned key = topo.getKey();
  BasisTable::const_iterator b = basisTable.find(key);
  ThrowRequireMsg( (b != basisTable.end()), "No basis available for this topology");

  const BasisTypeRCP basis = b->second;
  return basis;
}  

void fill_ref_vals(MDArray &refVals, const MDArray &refPoints, const shards::CellTopology &topo) {
  const BasisTypeRCP basis = getBasis(topo);
  basis->getValues(refVals, refPoints, Intrepid::OPERATOR_VALUE);
}


template<unsigned DIM> 
void parametric(std::vector<std::vector<double> > &val,
          const std::vector<double>               &para_coords, 
          const mesh::Entity                       element,
          const std::vector<mesh::FieldBase*>     &values_field,
          const mesh::BulkData&                   bulkData) {

  typedef Intrepid::FieldContainer<double>   MDArray;

  const unsigned dimension  = DIM; 
  const unsigned numIntrp   = 1; // FIXME

  mesh::Entity const* elem_node_rels = bulkData.begin_nodes(element);
  const unsigned num_nodes = bulkData.num_nodes(element);

  const mesh::Bucket & elem_bucket = bulkData.bucket(element);
  const CellTopologyData * const bucket_cell_topo_data = mesh::get_cell_topology(elem_bucket).getCellTopologyData();
  ThrowRequireMsg (bucket_cell_topo_data, __FILE__<<":"<<__LINE__<<" parametric::bogus topology");

  const shards::CellTopology topo(bucket_cell_topo_data);
  const unsigned numNodes = topo.getNodeCount();

  ThrowRequireMsg (topo.getDimension() == dimension,__FILE__<<":"<<__LINE__<<" Wrong spatical dimension"
    <<" for topology. Expected "<<dimension<<" found "<<topo.getDimension());
  ThrowRequireMsg (numNodes == num_nodes ,
    __FILE__<<":"<<__LINE__<<" Expected "<<numNodes<<" nodes but found "<<num_nodes);


  MDArray refPoints (numIntrp,  dimension);
  MDArray refVals   (           num_nodes, numIntrp);


  for (unsigned intrp = 0; intrp < numIntrp; intrp++)    
    for (unsigned i=0; i<dimension; ++i) refPoints(intrp,i) = para_coords[i];

  fill_ref_vals(refVals, refPoints, topo);

  for (unsigned intrp = 0; intrp < numIntrp; intrp++) {   
    const unsigned num_values = values_field.size();
    val.resize(num_values);
    for (unsigned ival=0; ival < num_values; ++ival) {

      const mesh::FieldBase &field = *values_field[ival];
      const mesh::Bucket & node_bucket = bulkData.bucket(elem_node_rels[0]);
      const unsigned bytes = bulkData.field_data_size_per_entity(field, node_bucket);
      const unsigned bytes_per_entry = field.data_traits().size_of;
      const unsigned num_entry = bytes/bytes_per_entry;

      ThrowRequireMsg (bytes == num_entry * bytes_per_entry,
         __FILE__<<":"<<__LINE__<<" Error:" <<"  bytes:" <<bytes<<"  num_entry:" <<num_entry
              <<"  bytes_per_entry:" <<bytes_per_entry);

      val[ival].resize(num_entry);
      MDArray basis     (num_entry,  num_nodes, numIntrp);
      MDArray dataVals  (num_entry,  num_nodes);
      MDArray fieldVals (num_entry,  numIntrp);
      fieldVals.initialize(0.0);

      // transfer reference basis values to physical frame values
      Intrepid::FunctionSpaceTools::HGRADtransformVALUE<double>(basis, refVals);
      for (unsigned iNode = 0; iNode < num_nodes; iNode++) {   
        const mesh::Entity node = elem_node_rels[iNode];
        const double * values = static_cast<const double*>(bulkData.field_data(field, node));
        for (unsigned e = 0; e < num_entry; ++e) dataVals(e, iNode) = values[e];
      }
      // evaluate function at specified points
      Intrepid::FunctionSpaceTools::evaluate<double>(fieldVals, dataVals, basis);
      for (unsigned e = 0; e < num_entry; ++e) val[ival][e] = fieldVals(e, intrp);
    }   
  }   
}

template void parametric<1>(std::vector<std::vector<double> > &val,
                    const std::vector<double> &para_coords, 
                    const mesh::Entity element,
                    const std::vector<mesh::FieldBase*> &values_field,
                    const mesh::BulkData& bulkData) ;
template void parametric<2>(std::vector<std::vector<double> > &val,
                    const std::vector<double> &para_coords, 
                    const mesh::Entity element,
                    const std::vector<mesh::FieldBase*> &values_field,
                    const mesh::BulkData& bulkData) ;
template void parametric<3>(std::vector<std::vector<double> > &val,
                    const std::vector<double> &para_coords, 
                    const mesh::Entity element,
                    const std::vector<mesh::FieldBase*> &values_field,
                    const mesh::BulkData& bulkData) ;
}
}
}

/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <vector>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/diag/EntityKey.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>

#include <stk_io/MeshReadWriteUtils.hpp>
#include <init/Ionit_Initializer.h>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>
#include <stk_util/diag/Writer.hpp>
#include <stk_util/diag/WriterExt.hpp>

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/diag/IdentProc.hpp>

#include <stk_search_util/stk_mesh/CreateBoundingBox.hpp>
#include <stk_search_util/stk_mesh/PrintBoundingBox.hpp>

using namespace stk_classic::search;
using namespace use_case;

using namespace stk_classic::diag;

typedef stk_classic::mesh::Field<double>                        ScalarField ;
typedef stk_classic::mesh::Field<double, stk_classic::mesh::Cartesian>  VectorField ;

static const size_t spatial_dimension = 3;

void
use_case_4_driver(stk_classic::ParallelMachine  comm,
		  const std::string &working_directory,
		  const std::string &range_mesh_filename,
		  const std::string &range_mesh_type,
		  const std::string &range_entity,
		  const std::string &domain_mesh_filename,
		  const std::string &domain_mesh_type,
		  const std::string &domain_entity)
{
  stk_classic::diag::WriterThrowSafe _write_throw_safe(dw());

  dw().m(LOG_SEARCH) << "Use case 4" << stk_classic::diag::push << stk_classic::diag::dendl;

  // ========================================================================
  // Initialize IO system.  Registers all element types and storage
  // types and the exodusII default database type.
  Ioss::Init::Initializer init_db;

  // Define range mesh...
  stk_classic::mesh::fem::FEMMetaData range_meta_data( spatial_dimension );
  stk_classic::io::MeshData range_mesh_data;
  std::string filename = working_directory + range_mesh_filename;
  stk_classic::io::create_input_mesh(range_mesh_type, filename, comm,
			     range_meta_data, range_mesh_data);
  range_meta_data.commit();

  stk_classic::mesh::BulkData range_bulk_data(range_meta_data.get_meta_data(range_meta_data) , comm);
  stk_classic::io::populate_bulk_data(range_bulk_data, range_mesh_data);


  // Define domain mesh...
  stk_classic::mesh::fem::FEMMetaData domain_meta_data( spatial_dimension );
  stk_classic::io::MeshData domain_mesh_data;
  filename = working_directory + domain_mesh_filename;
  stk_classic::io::create_input_mesh(domain_mesh_type, domain_mesh_filename, comm,
			     domain_meta_data, domain_mesh_data);
  domain_meta_data.commit();

  stk_classic::mesh::BulkData domain_bulk_data(domain_meta_data.get_meta_data(domain_meta_data) , comm);
  stk_classic::io::populate_bulk_data(domain_bulk_data, domain_mesh_data);

  // ========================================================================
  // For this use case, the domain consists of an axis-aligned
  // bounding box for each 'domain_entity' in the mesh.  The range is also an
  // axis-aligned bounding box of each 'range_entity'.

  VectorField *range_coord_field = range_meta_data.get_field<VectorField>("coordinates");
  std::vector<AxisAlignedBoundingBox3D> range_vector;
  stk_classic::search_util::build_axis_aligned_bbox(range_bulk_data,
                                            range_meta_data.entity_rank(range_entity),
                                            range_coord_field, range_vector);

  VectorField *domain_coord_field = domain_meta_data.get_field<VectorField>("coordinates");
  std::vector<AxisAlignedBoundingBox3D> domain_vector;
  stk_classic::search_util::build_axis_aligned_bbox(domain_bulk_data,
					domain_meta_data.entity_rank(domain_entity),
					domain_coord_field, domain_vector);

  if (range_vector.size() <= 100)
     dw().m(LOG_SEARCH) << "range  " << range_vector << dendl;
  else
    dw().m(LOG_SEARCH) << "range vector size  = " << range_vector.size() << dendl;

  if (domain_vector.size() <= 100)
    dw().m(LOG_SEARCH) << "domain " << domain_vector << dendl;
  else
    dw().m(LOG_SEARCH) << "domain vector size = " << domain_vector.size() << dendl;

  FactoryOrder order;
  order.m_communicator = comm;
  //order.m_algorithm = FactoryOrder::BIHTREE;

  dw() << "Search algorithm " << order.m_algorithm << dendl;
//  dw().m(LOG_SEARCH) << "Search tree " << *range_search << dendl;

  IdentProcRelation relation;

  stk_classic::search::coarse_search(relation, range_vector, domain_vector,order);

  if (relation.size() <= 100)
    dw().m(LOG_SEARCH) << "relation " << relation << dendl;
  else
    dw().m(LOG_SEARCH) << "relation size = " << relation.size() << dendl;

  dw().m(LOG_SEARCH) << stk_classic::diag::pop;
}


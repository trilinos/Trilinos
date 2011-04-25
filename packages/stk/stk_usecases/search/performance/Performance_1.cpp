/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <vector>
#include <fstream>

#include <assert.h>

#include <search/performance/Performance.hpp>
#include <stk_io/MeshReadWriteUtils.hpp>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>
#include <Ioss_ParallelUtils.h>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>
#include <stk_util/diag/Writer.hpp>
#include <stk_util/diag/WriterExt.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/SkinMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/diag/EntityKey.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>

#include <stk_search/SearchTypes.hpp>
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/diag/IdentProc.hpp>

#include <stk_search_util/stk_mesh/PrintEntityProc.hpp>
#include <stk_search_util/stk_mesh/CreateBoundingBox.hpp>
#include <stk_search_util/stk_mesh/PrintBoundingBox.hpp>

using namespace stk::diag;
using namespace use_case;

typedef stk::mesh::Field<double, stk::mesh::Cartesian>       VectorField ;

namespace {
  void get_domain_range_used(const IdentProcRelation &relation,
			     std::vector<stk::mesh::EntityKey> &domain_used,
			     std::vector<stk::mesh::EntityKey> &range_used,
			     int &on_proc, int &off_proc,
			     unsigned int my_rank,
			     bool same_mesh);

  template<typename T>
  std::ostream& operator<<(std::ostream& ostr, const std::vector<T> &v){
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(ostr, ", "));
    return ostr;
  }
}


void
performance_driver(stk::ParallelMachine  comm,
		   const std::string &working_directory,
		   const std::string &search_type,
		   const stk::search::Options &range,
		   const stk::search::Options &domain,
		   bool performance)
{
  static const size_t spatial_dimension = 3;
  
  stk::diag::WriterThrowSafe _write_throw_safe(dw());

  stk::diag::Timer timer("SearchPerformance", use_case::TIMER_SEARCH, use_case::timer());
  stk::diag::Timer timer_search_performance(domain.mesh_filename + " to " + range.mesh_filename, timer);
  stk::diag::Timer timer_range_bb("Range bounding box", timer_search_performance);
  stk::diag::Timer timer_domain_bb("Domain bounding box", timer_search_performance);
  stk::diag::Timer timer_range_search("Search", timer_search_performance);

  stk::diag::TimeBlock timer_search_performance_block(timer_search_performance);

  bool do_output = !performance;

  bool same_mesh = (range.mesh_filename == domain.mesh_filename);

  // For this use case, if 'node' is the range.entity_rank, then we
  // want to use the universal set nodes, not a nodeset.
  bool range_use_universal_set = true;
  bool domain_use_universal_set = true;
  if (range.entity == "nodeset") {
    range_use_universal_set = false;
  }
  if (domain.entity == "nodeset") {
    domain_use_universal_set = false;
  }

  // If 'faceset' is specified, then use non-skin faceset faces;
  // If 'face' is specified, then skin the model and use the faces in the 'skin' part.
  if (range.entity == "faceset") {
    range_use_universal_set = false;
  }
  if (domain.entity == "faceset") {
    domain_use_universal_set = false;
  }

  stk::CommAll comm_all( comm );

  dw() << "Performance use case: " << search_type   << stk::diag::dendl;
  dw() << "Range  Entity Type =  " << range.entity  << stk::diag::dendl;
  dw() << "Domain Entity Type =  " << domain.entity << stk::diag::dendl;

  // Initialize IO system.  Registers all element types and storage
  // types and the exodusII default database type.
  Ioss::Init::Initializer init_db;

  stk::mesh::fem::FEMMetaData range_meta_data( spatial_dimension );
  stk::io::MeshData range_mesh_data;

  dw() << "Build range metadata...\n";
  std::string filename = working_directory + range.mesh_filename;
  stk::io::create_input_mesh(range.mesh_type, filename, comm, range_meta_data, range_mesh_data);

  stk::mesh::Part *range_skin_part = NULL;
  stk::mesh::Part *domain_skin_part = NULL;
  stk::mesh::fem::CellTopology skin_top(shards::getCellTopologyData<shards::Quadrilateral<4> >());

  if ((range.entity == "face" && range_use_universal_set) ||
      (same_mesh && domain.entity == "face" && domain_use_universal_set)) {
    range_skin_part = &range_meta_data.declare_part("skin", skin_top);
  }

  range_meta_data.commit();

  dw() << "Build range bulkdata...\n";
  stk::mesh::BulkData range_bulk_data(range_meta_data.get_meta_data(range_meta_data) , comm);
  stk::io::populate_bulk_data(range_bulk_data, range_mesh_data);

  if ((range.entity == "face" && range_use_universal_set) ||
      (same_mesh && domain.entity == "face" && domain_use_universal_set)) {
    stk::mesh::skin_mesh( range_bulk_data, range_meta_data.element_rank(), range_skin_part);
  }

  stk::mesh::fem::FEMMetaData *domain_meta_data = NULL;
  stk::mesh::BulkData *domain_bulk_data = NULL;

  if (!same_mesh) {
    stk::io::MeshData domain_mesh_data;
    dw() << "Build domain metadata...\n";
    domain_meta_data = new stk::mesh::fem::FEMMetaData( spatial_dimension );
    filename = working_directory + domain.mesh_filename;
    stk::io::create_input_mesh(domain.mesh_type, filename, comm, *domain_meta_data, domain_mesh_data);

    if (domain.entity == "face" && domain_use_universal_set) {
      domain_skin_part = &domain_meta_data->declare_part("skin", skin_top);
    }
    domain_meta_data->commit();

    dw() << "Build domain bulkdata...\n";
    domain_bulk_data = new stk::mesh::BulkData(domain_meta_data->get_meta_data(*domain_meta_data) , comm);
    stk::io::populate_bulk_data(*domain_bulk_data, domain_mesh_data);

    if (domain.entity == "face" && domain_use_universal_set) {
      stk::mesh::skin_mesh( *domain_bulk_data, domain_meta_data->element_rank(), domain_skin_part);
    }
  } else {
    dw() << "Domain shares metadata and bulkdata with range...\n";
    domain_meta_data = &range_meta_data;
    domain_bulk_data = &range_bulk_data;
    domain_skin_part = range_skin_part;
  }

  // ========================================================================
  // End of mesh/geometry generation. The remainder is related to
  // searching on thie geometry...
  // ========================================================================


  stk::search::FactoryOrder order;
  order.m_communicator = comm;

  VectorField *range_coord_field = range_meta_data.get_field<VectorField>("coordinates");
  VectorField *domain_coord_field = domain_meta_data->get_field<VectorField>("coordinates");

  IdentProcRelation relation;
  size_t domain_vector_size = 0;
  size_t range_vector_size = 0;

  if (search_type == "point_in_box") {
    dw() << "For the " << search_type << " use case:\n"
	 << "\t* domain consists of an axis-aligned bounding box for each " << domain.entity << " in the mesh.\n"
	 << "\t* range is a PointBoundingBox3D at the centroid of each " << range.entity << " in the mesh.\n\n";

    std::vector<PointBoundingBox3D> range_vector;

    {
      dw() << "Build range bbox...\n";
      stk::diag::TimeBlock __timer_range_bb(timer_range_bb);
      stk::search_util::build_centroid_bbox(range_bulk_data,
					    range_meta_data.entity_rank(range.entity),
					    range_coord_field, range_vector,
					    range_use_universal_set);
    }

    std::vector<AxisAlignedBoundingBox3D> domain_vector;
    {
      dw() << "Build domain bbox...\n";
      stk::diag::TimeBlock __timer_domain_bb(timer_domain_bb);
      stk::search_util::build_axis_aligned_bbox(*domain_bulk_data,
						domain_meta_data->entity_rank(domain.entity),
						domain_coord_field, domain_vector,
						domain_use_universal_set,
						stk::search_util::OffsetScaleOp(domain.scale,
										domain.offset));
    }

    {
      stk::diag::TimeBlock __timer_range_search(timer_range_search);
      dw() << "Build search object...\n";

      dw() << "Perform search...\n";
      stk::search::coarse_search(relation, range_vector, domain_vector,order);
    }

    if (do_output) {
      dw() << "Search algorithm " << order.m_algorithm << dendl;
      dw() << "range  " << range_vector  << stk::diag::dendl;
      dw() << "domain " << domain_vector << stk::diag::dendl;
      stk::search_util::print_stk_mesh_relation_map(dw(), stk::mesh::fem::entity_rank_names(spatial_dimension), relation);
    }
    range_vector_size = range_vector.size();
    domain_vector_size = domain_vector.size();
  }

  else if (search_type == "overlap") {  // NEED BETTER NAMES FOR SEARCHES?
    dw() << "For the " << search_type << " use case:\n"
	 << "\t* domain consists of an axis-aligned bounding box for each " << domain.entity << " in the mesh.\n"
	 << "\t* range is also an axis-aligned bounding box for each " << range.entity << " in the mesh.\n\n";

    std::vector<AxisAlignedBoundingBox3D> range_vector;
    {
      stk::diag::TimeBlock __timer_range_bb(timer_range_bb);
      stk::search_util::build_axis_aligned_bbox(range_bulk_data,
						range_meta_data.entity_rank(range.entity),
						range_coord_field, range_vector,
						range_use_universal_set,
						stk::search_util::OffsetScaleOp(range.scale,
										range.offset));
    }

    std::vector<AxisAlignedBoundingBox3D> domain_vector;
    {
      stk::diag::TimeBlock __timer_domain_bb(timer_domain_bb);
      stk::search_util::build_axis_aligned_bbox(*domain_bulk_data,
						domain_meta_data->entity_rank(domain.entity),
						domain_coord_field, domain_vector,
						domain_use_universal_set,
						stk::search_util::OffsetScaleOp(domain.scale,
										domain.offset));
    }

    {
      stk::diag::TimeBlock __timer_range_search(timer_range_search);
      stk::search::coarse_search(relation, range_vector, domain_vector, order);
    }

    if (do_output) {
      dw() << "Search algorithm " << order.m_algorithm << dendl;
      dw() << "range  " << range_vector  << stk::diag::dendl;
      dw() << "domain " << domain_vector << stk::diag::dendl;
      stk::search_util::print_stk_mesh_relation_map(dw(), stk::mesh::fem::entity_rank_names(spatial_dimension), relation);
    }
    range_vector_size = range_vector.size();
    domain_vector_size = domain_vector.size();
  }
  else {
    std::cerr << "ERROR: Invalid search option '" << search_type << "'.\n";
    return;
  }

  // Output metrics/timing information...
  // -- Size of domain and range
  // -- Percentage of domain objects matched with at least one range object
  // -- Percentage of range objects matched with at least one domain object
  // -- Need idea of false positives -- matches found in coarse that
  //    are filtered by the detailed... Also, cost of a detailed
  //    search/entity...
  //
  std::vector<stk::mesh::EntityKey> domain_used;
  std::vector<stk::mesh::EntityKey> range_used;
  int lsizes[] = {0,0,0,0,0,0,0};
  int gsizes[] = {0,0,0,0,0,0,0};

  get_domain_range_used(relation, domain_used, range_used, lsizes[4], lsizes[5],
			comm_all.parallel_rank(), same_mesh);

  lsizes[0] = static_cast<int>(domain_used.size());
  lsizes[1] = static_cast<int>(range_used.size());
  lsizes[2] = static_cast<int>(domain_vector_size);
  lsizes[3] = static_cast<int>(range_vector_size);
  lsizes[6] = static_cast<int>(relation.size());

  MPI_Allreduce(lsizes, gsizes, sizeof(lsizes)/sizeof(lsizes[0]), MPI_INT, MPI_SUM, comm);

  dw() << "\nSearch type = " << search_type << dendl
       << "Domain Size = " << gsizes[2] << " " << domain.entity
       << "  Used = "      << gsizes[0] << "  (" << (100*static_cast<long>(gsizes[0]))/(gsizes[2]) << "%)\n"
       << "Range  Size = " << gsizes[3] << " " << range.entity
       << "  Used = "      << gsizes[1] << "  (" << (100*static_cast<long>(gsizes[1]))/(gsizes[3]) << "%)\n"
       << "Total matches: " << gsizes[4] + gsizes[5];
  if (gsizes[4] + gsizes[5] > 0) {
    dw() << "  On processor: "  << gsizes[4] << "  (" << (101*static_cast<long>(gsizes[4])-1)/(gsizes[4] + gsizes[5]) << "%)"
         << "  Off processor: " << gsizes[5] << "  (" << (101*static_cast<long>(gsizes[5])-1)/(gsizes[4] + gsizes[5]) << "%)";
  }
  dw() << std::endl;

  // Different search types:
  // -- expect one-to-one match - each domain object matches exactly
  //    one range object.
  // -- expect few matches
  // -- large domain/small range or vice-versa

  if (!same_mesh) {
    delete domain_bulk_data;
    delete domain_meta_data;
  }
}

namespace {
  void uniqify(std::vector<stk::mesh::EntityKey> &vec)
  {
    std::sort(vec.begin(), vec.end());
    std::vector<stk::mesh::EntityKey>::iterator i = std::unique(vec.begin(), vec.end());
    vec.erase(i, vec.end());
  }

  void get_domain_range_used(const IdentProcRelation &relation,
			     std::vector<stk::mesh::EntityKey> &domain_used,
			     std::vector<stk::mesh::EntityKey> &range_used,
			     int &on_proc, int &off_proc, unsigned int my_rank, bool same_mesh)
  {
    IdentProcRelation::const_iterator i=relation.begin(), end=relation.end();
    for (;i!=end; ++i) {
      if (same_mesh && i->first.ident == i->second.ident)
	continue;

      std::size_t domain_owning_proc = i->first.proc;
      std::size_t range_owning_proc = i->second.proc;
      if (domain_owning_proc == my_rank) {
	domain_used.push_back(i->first.ident);
	if (range_owning_proc == my_rank) {
	  on_proc++;
	} else {
	  off_proc++;
	}
     }

      if (range_owning_proc == my_rank) {
	range_used.push_back(i->second.ident);
	if (domain_owning_proc != my_rank) {
	  off_proc++;
	}
      }
    }
    uniqify(domain_used);
    uniqify(range_used);
  }
}


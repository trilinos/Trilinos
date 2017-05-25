#include <Ioss_CodeTypes.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_StructuredBlock.h>
#include <Ioss_Utils.h>
#include <cgns/Iocgns_DecompositionData.h>
#include <cgns/Iocgns_Utils.h>

#include <cgnsconfig.h>
#if CG_BUILD_PARALLEL
#include <pcgnslib.h>
#else
#include <cgnslib.h>
#endif

#include <Ioss_TerminalColor.h>

#include <algorithm>
#include <numeric>

#include <assert.h>

#define IOSS_DEBUG_OUTPUT 0
namespace {
  int rank = 0;
#define OUTPUT                                                                                     \
  if (rank == 0)                                                                                   \
  std::cerr

// ZOLTAN Callback functions...

#if !defined(NO_ZOLTAN_SUPPORT)
  int zoltan_num_dim(void *data, int *ierr)
  {
    // Return dimensionality of coordinate data.
    Iocgns::DecompositionDataBase *zdata = (Iocgns::DecompositionDataBase *)(data);

    *ierr = ZOLTAN_OK;
    return zdata->spatial_dimension();
  }

  int zoltan_num_obj(void *data, int *ierr)
  {
    // Return number of objects (element count) on this processor...
    Iocgns::DecompositionDataBase *zdata = (Iocgns::DecompositionDataBase *)(data);

    *ierr = ZOLTAN_OK;
    return zdata->decomp_elem_count();
  }

  void zoltan_obj_list(void *data, int ngid_ent, int nlid_ent, ZOLTAN_ID_PTR gids,
                       ZOLTAN_ID_PTR lids, int wdim, float *wgts, int *ierr)
  {
    // Return list of object IDs, both local and global.
    Iocgns::DecompositionDataBase *zdata = (Iocgns::DecompositionDataBase *)(data);

    // At the time this is called, we don't have much information
    // These routines are the ones that are developing that
    // information...
    size_t element_count  = zdata->decomp_elem_count();
    size_t element_offset = zdata->decomp_elem_offset();

    *ierr = ZOLTAN_OK;

    if (lids) {
      std::iota(lids, lids + element_count, 0);
    }

    if (wdim) {
      std::fill(wgts, wgts + element_count, 1.0);
    }

    if (ngid_ent == 1) {
      std::iota(gids, gids + element_count, element_offset);
    }
    else if (ngid_ent == 2) {
      int64_t *global_ids = (int64_t *)gids;
      std::iota(global_ids, global_ids + element_count, element_offset);
    }
    else {
      *ierr = ZOLTAN_FATAL;
    }
    return;
  }

  void zoltan_geom(void *data, int ngid_ent, int nlid_ent, int nobj, ZOLTAN_ID_PTR gids,
                   ZOLTAN_ID_PTR lids, int ndim, double *geom, int *ierr)
  {
    // Return coordinates for objects.
    Iocgns::DecompositionDataBase *zdata = (Iocgns::DecompositionDataBase *)(data);

    std::copy(zdata->centroids().begin(), zdata->centroids().end(), &geom[0]);

    *ierr = ZOLTAN_OK;
    return;
  }
#endif

  // These are used for structured parallel decomposition...
  void create_zone_data(int cgnsFilePtr, std::vector<Iocgns::StructuredZoneData *> &zones)
  {
    int base      = 1;
    int num_zones = 0;
    cg_nzones(cgnsFilePtr, base, &num_zones);

    std::map<std::string, int> zone_name_map;

    for (cgsize_t zone = 1; zone <= num_zones; zone++) {
      cgsize_t size[9];
      char     zone_name[33];
      cg_zone_read(cgnsFilePtr, base, zone, zone_name, size);
      zone_name_map[zone_name] = zone;

      assert(size[0] - 1 == size[3]);
      assert(size[1] - 1 == size[4]);
      assert(size[2] - 1 == size[5]);

      assert(size[6] == 0);
      assert(size[7] == 0);
      assert(size[8] == 0);

      cgsize_t index_dim = 0;
      cg_index_dim(cgnsFilePtr, base, zone, &index_dim);

      auto *zone_data = new Iocgns::StructuredZoneData(zone_name, zone, size[3], size[4], size[5]);
      zone_data->m_adam = zone_data;
      zones.push_back(zone_data);

      // Handle zone-grid-connectivity...
      int nconn = 0;
      cg_n1to1(cgnsFilePtr, base, zone, &nconn);
      for (int i = 0; i < nconn; i++) {
        char connectname[33];
        char donorname[33];
        std::array<cgsize_t, 6> range;
        std::array<cgsize_t, 6> donor_range;
        Ioss::IJK_t transform;

        cg_1to1_read(cgnsFilePtr, base, zone, i + 1, connectname, donorname, range.data(),
                     donor_range.data(), transform.data());

        // Get number of nodes shared with other "previous" zones...
        // A "previous" zone will have a lower zone number this this zone...
        int  donor_zone = -1;
        auto donor_iter = zone_name_map.find(donorname);
        if (donor_iter != zone_name_map.end()) {
          donor_zone = (*donor_iter).second;
        }
        Ioss::IJK_t range_beg{{range[0], range[1], range[2]}};
        Ioss::IJK_t range_end{{range[3], range[4], range[5]}};
        Ioss::IJK_t donor_beg{{donor_range[0], donor_range[1], donor_range[2]}};
        Ioss::IJK_t donor_end{{donor_range[3], donor_range[4], donor_range[5]}};

#if defined(IOSS_DEBUG_OUTPUT)
        OUTPUT << "Adding zgc " << connectname << " to " << zone_name << " donor: " << donorname
               << "\n";
#endif
        bool owns_nodes = zone < donor_zone || donor_zone == -1;
        zone_data->m_zoneConnectivity.emplace_back(connectname, zone, donorname, donor_zone,
                                                   transform, range_beg, range_end, donor_beg,
                                                   donor_end, owns_nodes);
      }
    }

    // If there are any Structured blocks, need to iterate them and their 1-to-1 connections
    // and update the donor_zone id for zones that had not yet been processed at the time of
    // definition...
    for (auto &zone : zones) {
      for (auto &conn : zone->m_zoneConnectivity) {
        if (conn.m_donorZone < 0) {
          auto donor_iter = zone_name_map.find(conn.m_donorName);
          assert(donor_iter != zone_name_map.end());
          conn.m_donorZone = (*donor_iter).second;
        }
      }
    }
  }

  ssize_t proc_with_minimum_work(const std::vector<size_t> &work, ssize_t exclude_proc = -1)
  {
    size_t  min_work = std::numeric_limits<size_t>::max();
    ssize_t min_proc = -1;
    for (ssize_t i = 0; i < (ssize_t)work.size(); i++) {
      if (work[i] < min_work && i != exclude_proc) {
        min_work = work[i];
        min_proc = i;
        if (min_work == 0) {
          break;
        }
      }
    }
    return min_proc;
  }
}

namespace Iocgns {
  template DecompositionData<int>::DecompositionData(const Ioss::PropertyManager &props,
                                                     MPI_Comm                     communicator);
  template DecompositionData<int64_t>::DecompositionData(const Ioss::PropertyManager &props,
                                                         MPI_Comm                     communicator);

  template <typename INT>
  DecompositionData<INT>::DecompositionData(const Ioss::PropertyManager &props,
                                            MPI_Comm                     communicator)
      : DecompositionDataBase(communicator), m_loadBalanceThreshold(1.4),
        m_decomposition(props, communicator)
  {
    rank = m_decomposition.m_processor;

    if (props.exists("LOAD_BALANCE_THRESHOLD")) {
      if (props.get("LOAD_BALANCE_THRESHOLD").get_type() == Ioss::Property::STRING) {
        std::string lb_thresh  = props.get("LOAD_BALANCE_THRESHOLD").get_string();
        m_loadBalanceThreshold = std::strtod(lb_thresh.c_str(), nullptr);
      }
      else if (props.get("LOAD_BALANCE_THRESHOLD").get_type() == Ioss::Property::REAL) {
        m_loadBalanceThreshold = props.get("LOAD_BALANCE_THRESHOLD").get_real();
      }
    }
  }

  template <typename INT>
  void DecompositionData<INT>::decompose_model(int filePtr, CG_ZoneType_t common_zone_type)
  {
    if (common_zone_type == CG_Unstructured) {
      decompose_unstructured(filePtr);
    }
    else if (common_zone_type == CG_Structured) {
      decompose_structured(filePtr);
    }
    else {
      std::ostringstream errmsg;
      errmsg << "ERROR: CGNS: The common zone type is not of type Unstructured or Structured "
                "which are the only types currently supported";
      IOSS_ERROR(errmsg);
    }
  }

  template <typename INT> void DecompositionData<INT>::decompose_structured(int filePtr)
  {
    create_zone_data(filePtr, m_structuredZones);
    if (m_structuredZones.empty()) {
      return;
    }

    size_t work = 0;
    for (const auto &z : m_structuredZones) {
      work += z->work();
      assert(z->is_active());
    }

    size_t new_zone_id = m_structuredZones.size() + 1;
    size_t px          = 0;
    size_t num_split   = 0;
    bool   split       = false;
    double avg_work    = (double)work / m_decomposition.m_processorCount;

    auto num_active = m_structuredZones.size();
    OUTPUT << "Decomposing structured mesh with " << num_active << " zones for "
           << m_decomposition.m_processorCount << " processors.\nAverage workload is " << avg_work
           << ", Load Balance Threshold is " << m_loadBalanceThreshold << ", Work range "
           << avg_work / m_loadBalanceThreshold << " to " << avg_work * m_loadBalanceThreshold
           << "\n";

    if (avg_work < 1.0) {
      OUTPUT << "ERROR: Model size too small to distribute over "
             << m_decomposition.m_processorCount << " processors.\n";
      std::exit(EXIT_FAILURE);
    }

#if defined(IOSS_DEBUG_OUTPUT)
    OUTPUT << "========================================================================\n";
    OUTPUT << "Pre-Splitting:\n";
#endif
    // Split all blocks where block->work() > avg_work * m_loadBalanceThreshold
    do {
      auto zone_new(m_structuredZones);
      split = false;
      for (auto zone : m_structuredZones) {
        if (zone->is_active() && zone->work() > avg_work * m_loadBalanceThreshold) {
          auto children = zone->split(new_zone_id, zone->work() / avg_work);
          if (children.first != nullptr && children.second != nullptr) {
            zone_new.push_back(children.first);
            zone_new.push_back(children.second);
            split = true;
            new_zone_id += 2;
          }
        }
      }
      std::swap(zone_new, m_structuredZones);
    } while (split);
#if defined(IOSS_DEBUG_OUTPUT)
    OUTPUT << "========================================================================\n";
#endif
    do {
      // Sort zones based on work.  Most work first..
      // TODO: Possibly filter 'zones' down to only active zones to
      // reduce sort and iteration time.
      std::sort(m_structuredZones.begin(), m_structuredZones.end(),
                [](Iocgns::StructuredZoneData *a, Iocgns::StructuredZoneData *b) {
                  return a->work() > b->work();
                });

      std::vector<size_t> work_vector(m_decomposition.m_processorCount);

      auto zone_new(m_structuredZones);
      for (auto &zone : m_structuredZones) {
        zone->m_proc = -1;
        if (zone->is_active()) {
          // Assign zone to processor with minimum work...
          ssize_t proc = proc_with_minimum_work(work_vector);

          // See if any other zone on this processor has the same adam zone...
          // TODO: Currently only do one "re-search".  Need to do something
          // better to make sure; or be able to handle this condition correctly.
          for (auto &pzone : m_structuredZones) {
            if (pzone->is_active() && pzone->m_proc == proc) {
              if (pzone->m_adam == zone->m_adam) {
                proc = proc_with_minimum_work(work_vector, proc);
                break;
              }
            }
          }

          zone->m_proc = proc;
          work_vector[proc] += zone->work();
#if defined(IOSS_DEBUG_OUTPUT)
          OUTPUT << "Assigning zone " << zone->m_zone << " with work " << zone->work()
                 << " to processor " << proc << "\n";
#endif
        }
      }

      // Calculate workload ratio for each processor...
      px = 0; // Number of processors where workload ratio exceeds threshold.
      std::vector<bool> exceeds(m_decomposition.m_processorCount);
      for (size_t i = 0; i < work_vector.size(); i++) {
        double workload_ratio = double(work_vector[i]) / double(avg_work);
#if defined(IOSS_DEBUG_OUTPUT)
        OUTPUT << "Processor " << i << " work: " << work_vector[i]
               << ", workload ratio: " << workload_ratio << "\n";
#endif
        if (workload_ratio > m_loadBalanceThreshold) {
          exceeds[i] = true;
          px++;
        }
      }
#if defined(IOSS_DEBUG_OUTPUT)
      OUTPUT << "Workload threshold exceeded on " << px << " processors.\n";
#endif
      num_split = 0;
      if (px > 0) {
        for (auto zone : m_structuredZones) {
          if (zone->is_active() && exceeds[zone->m_proc]) {
            // Since 'zones' is sorted from most work to least,
            // we just iterate zones and check whether the zone
            // is on a proc where the threshold was exceeded.
            // if so, split the block and set exceeds[proc] to false;
            // Exit the loop when num_split >= px.
            auto children = zone->split(new_zone_id);
            if (children.first != nullptr && children.second != nullptr) {
              zone_new.push_back(children.first);
              zone_new.push_back(children.second);

              new_zone_id += 2;
              exceeds[zone->m_proc] = false;
              num_split++;
              if (num_split >= px) {
                break;
              }
            }
          }
        }
        std::swap(zone_new, m_structuredZones);
      }
      auto active = std::count_if(m_structuredZones.begin(), m_structuredZones.end(),
                                  [](Iocgns::StructuredZoneData *a) { return a->is_active(); });
#if defined(IOSS_DEBUG_OUTPUT)
      OUTPUT << "Number of active zones = " << active << ", average work = " << avg_work << "\n";
      OUTPUT << "========================================================================\n";
#endif
    } while (px > 0 && num_split > 0);

    std::sort(m_structuredZones.begin(), m_structuredZones.end(),
              [](Iocgns::StructuredZoneData *a, Iocgns::StructuredZoneData *b) {
                return a->m_zone < b->m_zone;
              });

    for (auto zone : m_structuredZones) {
      if (zone->is_active()) {
        zone->resolve_zgc_split_donor(m_structuredZones);
      }
    }

    // Update and Output the processor assignments
    for (auto &zone : m_structuredZones) {
      if (zone->is_active()) {
        zone->update_zgc_processor(m_structuredZones);
#if defined(IOSS_DEBUG_OUTPUT)
        OUTPUT << "Zone " << zone->m_zone << " assigned to processor " << zone->m_proc
               << ", Adam zone = " << zone->m_adam->m_zone << "\n";
        auto zgcs = zone->m_zoneConnectivity;
        for (auto &zgc : zgcs) {
          OUTPUT << zgc << "\n";
        }
#endif
      }
    }

    for (auto &zone : m_structuredZones) {
      if (!zone->is_active()) {
        zone->m_proc = -1;
      }
    }

#if defined(IOSS_DEBUG_OUTPUT)
    OUTPUT << Ioss::trmclr::green << "Returning from decomposition\n" << Ioss::trmclr::normal;
#endif
  }

  template <typename INT> void DecompositionData<INT>::decompose_unstructured(int filePtr)
  {
    // Initial decomposition is linear where processor #p contains
    // elements from (#p * #element/#proc) to (#p+1 * #element/#proc)

    // ========================================================================
    // Get the number of zones (element blocks) in the mesh...
    int num_zones = 0;
    int base      = 1; // Only single base supported so far.

    {
      cgsize_t cell_dimension = 0;
      cgsize_t phys_dimension = 0;
      char     base_name[33];
      cg_base_read(filePtr, base, base_name, &cell_dimension, &phys_dimension);
      m_decomposition.m_spatialDimension = phys_dimension;
    }

    cg_nzones(filePtr, base, &num_zones);
    m_zones.resize(num_zones + 1); // Use 1-based zones.

    size_t global_cell_node_count = 0;
    size_t global_element_count   = 0;
    for (int zone = 1; zone <= num_zones; zone++) {
      // All zones are "Unstructured" since this was checked prior to
      // calling this function...
      cgsize_t size[3];
      char     zone_name[33];
      cg_zone_read(filePtr, base, zone, zone_name, size);

      INT total_block_nodes = size[0];
      INT total_block_elem  = size[1];

      m_zones[zone].m_nodeCount     = total_block_nodes;
      m_zones[zone].m_nodeOffset    = global_cell_node_count;
      m_zones[zone].m_name          = zone_name;
      m_zones[zone].m_elementOffset = global_element_count;
      global_cell_node_count += total_block_nodes;
      global_element_count += total_block_elem;
    }

    // Generate element_dist/node_dist --  size m_decomposition.m_processorCount + 1
    // processor p contains all elements/nodes from X_dist[p] .. X_dist[p+1]
    m_decomposition.generate_entity_distributions(global_cell_node_count, global_element_count);

    generate_adjacency_list(filePtr, m_decomposition);

    // Get min and max node used on this processor...
    auto min_max =
        std::minmax_element(m_decomposition.m_adjacency.begin(), m_decomposition.m_adjacency.end());
    INT min_node = *(min_max.first);
    INT max_node = *(min_max.second);
    generate_zone_shared_nodes(filePtr, min_node, max_node);

    // Now iterate adjacency list and update any "zone_shared_node" nodes
    // with their "sharee"
    if (!m_zoneSharedMap.empty()) {
      for (auto &node : m_decomposition.m_adjacency) {
        auto alias = m_zoneSharedMap.find(node);
        if (alias != m_zoneSharedMap.end()) {
          node = (*alias).second;
        }
      }
    }

#if IOSS_DEBUG_OUTPUT
    OUTPUT << "Processor " << m_decomposition.m_processor << " has " << decomp_elem_count()
           << " elements; offset = " << decomp_elem_offset() << "\n";
    OUTPUT << "Processor " << m_decomposition.m_processor << " has " << decomp_node_count()
           << " nodes; offset = " << decomp_node_offset() << ".\n";
#endif

    if (m_decomposition.needs_centroids()) {
      // Get my coordinate data using direct cgns calls
      std::vector<double> x(decomp_node_count());
      std::vector<double> y;
      std::vector<double> z;

      get_file_node_coordinates(filePtr, 0, TOPTR(x));
      if (m_decomposition.m_spatialDimension > 1) {
        y.resize(decomp_node_count());
        get_file_node_coordinates(filePtr, 1, TOPTR(y));
      }
      if (m_decomposition.m_spatialDimension > 2) {
        z.resize(decomp_node_count());
        get_file_node_coordinates(filePtr, 2, TOPTR(z));
      }

      m_decomposition.calculate_element_centroids(x, y, z);
    }

#if !defined(NO_ZOLTAN_SUPPORT)
    float version = 0.0;
    Zoltan_Initialize(0, nullptr, &version);

    Zoltan zz(m_decomposition.m_comm);

    // Register Zoltan Callback functions...
    zz.Set_Num_Obj_Fn(zoltan_num_obj, this);
    zz.Set_Obj_List_Fn(zoltan_obj_list, this);
    zz.Set_Num_Geom_Fn(zoltan_num_dim, this);
    zz.Set_Geom_Multi_Fn(zoltan_geom, this);
#endif

    m_decomposition.decompose_model(
#if !defined(NO_ZOLTAN_SUPPORT)
        zz,
#endif
        m_elementBlocks);

    if (!m_sideSets.empty()) {
      // Create elemGTL map which is used for sidesets (also element sets)
      build_global_to_local_elem_map();
    }

    get_sideset_data(filePtr);

    // Have all the decomposition data needed (except for boundary
    // conditions...)
    // Can now populate the Ioss metadata...
  }

  template <typename INT>
  void DecompositionData<INT>::generate_zone_shared_nodes(int filePtr, INT min_node, INT max_node)
  {
    // Begin of Zone-Shared node information

    // Modify adjacency list based on shared nodes between zones...
    // Need the map from "global" to "global-shared"
    // * This is not necessarily nodes only on my processor since connectivity can include
    //   nodes other than those I own.
    // * Potentially large number of shared nodes; practically small(?)

    // * Maintain hash map from old id to new (if any)
    // * TODO: Determine whether the node is used on this processor...

    int base = 1; // Only single base supported so far.

    // Donor zone is always lower numbered, so zone 1 has no donor zone. Start at zone 2.
    for (cgsize_t zone = 2; zone < (cgsize_t)m_zones.size(); zone++) {

      // Determine number of "shared" nodes (shared with other zones)
      int nconn = 0;
      cg_nconns(filePtr, base, zone, &nconn);
      for (int i = 0; i < nconn; i++) {
        char                      connectname[33];
        CG_GridLocation_t         location;
        CG_GridConnectivityType_t connect_type;
        CG_PointSetType_t         ptset_type;
        cgsize_t                  npnts = 0;
        char                      donorname[33];
        CG_ZoneType_t             donor_zonetype;
        CG_PointSetType_t         donor_ptset_type;
        CG_DataType_t             donor_datatype;
        cgsize_t                  ndata_donor;

        cg_conn_info(filePtr, base, zone, i + 1, connectname, &location, &connect_type, &ptset_type,
                     &npnts, donorname, &donor_zonetype, &donor_ptset_type, &donor_datatype,
                     &ndata_donor);

        if (connect_type != CG_Abutting1to1 || ptset_type != CG_PointList ||
            donor_ptset_type != CG_PointListDonor) {
          std::ostringstream errmsg;
          errmsg << "ERROR: CGNS: Zone " << zone
                 << " adjacency data is not correct type. Require Abutting1to1 and PointList."
                 << connect_type << "\t" << ptset_type << "\t" << donor_ptset_type;
          IOSS_ERROR(errmsg);
        }

        // Verify data consistency...
        if (npnts != ndata_donor) {
          std::ostringstream errmsg;
          errmsg << "ERROR: CGNS: Zone " << zone << " point count (" << npnts
                 << ") does not match donor point count (" << ndata_donor << ").";
          IOSS_ERROR(errmsg);
        }

        // Get number of nodes shared with other "previous" zones...
        // A "previous" zone will have a lower zone number this this zone...
        std::string dz_name(donorname);
        int         dz = 1;
        for (; dz < zone; dz++) {
          if (m_zones[dz].m_name == dz_name)
            break;
        }

        if (dz != zone) {
          std::cout << "Zone " << zone << " shares " << npnts << " nodes with " << donorname
                    << "\n";

          std::vector<cgsize_t> points(npnts);
          std::vector<cgsize_t> donors(npnts);

          cg_conn_read(filePtr, base, zone, i + 1, TOPTR(points), donor_datatype, TOPTR(donors));

          for (int j = 0; j < npnts; j++) {
            cgsize_t point = points[j] - 1 + m_zones[zone].m_nodeOffset;
            if (point >= min_node && point <= max_node) {
              cgsize_t donor = donors[j] - 1 + m_zones[dz].m_nodeOffset;

              // See if 'donor' is mapped to a different node already
              auto donor_map = m_zoneSharedMap.find(donor);
              if (donor_map != m_zoneSharedMap.end()) {
                donor = (*donor_map).second;
              }
              assert(m_zoneSharedMap.find(point) == m_zoneSharedMap.end());
              m_zoneSharedMap.insert({point, donor});
            }
          }
        }
      }
    }
  }

  template <typename INT>
  void DecompositionData<INT>::generate_adjacency_list(int                       filePtr,
                                                       Ioss::Decomposition<INT> &decomposition)
  {
    int base = 1; // Only single base supported so far.

    // Range of elements currently handled by this processor [)
    size_t p_start = decomp_elem_offset();
    size_t p_end   = p_start + decomp_elem_count();

    assert(sizeof(INT) == sizeof(cgsize_t));
    size_t sum    = 0; // Size of adjacency vector.
    size_t offset = 0;

    int num_zones        = 0;
    INT zone_node_offset = 0;

    cg_nzones(filePtr, base, &num_zones);
    for (int zone = 1; zone <= num_zones; zone++) {
      cgsize_t size[3];
      char     zone_name[33];
      cg_zone_read(filePtr, base, zone, zone_name, size);

      INT total_elements = size[1];
      // NOTE: A Zone will have a single set of nodes, but can have
      //       multiple sections each with their own element type...
      //       Keep treating sections as element blocks until we
      //       have handled 'size[1]' number of elements; the remaining
      //       sections are then the boundary faces (?)
      int num_sections = 0;
      cg_nsections(filePtr, base, zone, &num_sections);

      size_t last_blk_location = 0;
      for (int is = 1; is <= num_sections; is++) {
        char             section_name[33];
        CG_ElementType_t e_type;
        cgsize_t         el_start    = 0;
        cgsize_t         el_end      = 0;
        int              num_bndry   = 0;
        int              parent_flag = 0;

        // Get the type of elements in this section...
        cg_section_read(filePtr, base, zone, is, section_name, &e_type, &el_start, &el_end,
                        &num_bndry, &parent_flag);

        INT num_entity = el_end - el_start + 1;

        if (parent_flag == 0 && total_elements > 0) {
          total_elements -= num_entity;

          // Range of elements in element block b [)
          size_t b_start = offset; // offset is index of first element in this block...
          offset += num_entity;
          size_t b_end = offset;

          int element_nodes;
          cg_npe(e_type, &element_nodes);

          if (b_start < p_end && p_start < b_end) {
            // Some of this blocks elements are on this processor...
            size_t overlap = std::min(b_end, p_end) - std::max(b_start, p_start);
            sum += overlap * element_nodes;
          }

          Ioss::BlockDecompositionData block;
          block.zone_          = zone;
          block.section_       = is;
          block.name_          = zone_name;
          block.topologyType   = Utils::map_cgns_to_topology_type(e_type);
          block.nodesPerEntity = element_nodes;
          block.fileCount      = num_entity;
          block.zoneNodeOffset = zone_node_offset;

          last_blk_location = m_elementBlocks.size();
          m_elementBlocks.push_back(block);
        }
        else {
          // This is a boundary-condition -- sideset (?)
          std::string ss_name(section_name);

          Ioss::SetDecompositionData sset;
          sset.zone_            = zone;
          sset.section_         = is;
          sset.name_            = ss_name;
          sset.fileCount        = num_entity;
          sset.topologyType     = Utils::map_cgns_to_topology_type(e_type);
          sset.parentBlockIndex = last_blk_location;
          m_sideSets.push_back(sset);
        }
      }
      zone_node_offset += size[0];
    }
    int block_count = (int)m_elementBlocks.size();

    // Get the global element block index list at this time also.
    // The global element at index 'I' (0-based) is on block B
    // if global_block_index[B] <= I && global_block_index[B+1] < I
    // allocate and TODO: Fill
    m_decomposition.m_fileBlockIndex.reserve(block_count + 1);
    for (auto block : m_elementBlocks) {
      m_decomposition.m_fileBlockIndex.push_back(block.file_count());
    }
    m_decomposition.m_fileBlockIndex.push_back(0);
    Ioss::Utils::generate_index(m_decomposition.m_fileBlockIndex);

    // Make sure 'sum' can fit in INT...
    INT tmp_sum = (INT)sum;
    if ((size_t)tmp_sum != sum) {
      std::ostringstream errmsg;
      errmsg << "ERROR: The decomposition of this mesh requires 64-bit integers, but is being\n"
             << "       run with 32-bit integer code. Please rerun with the property "
                "INTEGER_SIZE_API\n"
             << "       set to 8. The details of how to do this vary with the code that is being "
                "run.\n"
             << "       Contact gdsjaar@sandia.gov for more details.\n";
      OUTPUT << errmsg.str();
      exit(EXIT_FAILURE);
    }

    // Now, populate the vectors...
    decomposition.m_pointer.reserve(decomp_elem_count() + 1);
    decomposition.m_adjacency.reserve(sum);
    offset = 0;
    sum    = 0; // Size of adjacency vector.

    for (auto &block : m_elementBlocks) {
      // Range of elements in element block b [)
      size_t b_start = offset; // offset is index of first element in this block...
      offset += block.file_count();
      size_t b_end = b_start + block.file_count();

      ssize_t overlap      = std::min(b_end, p_end) - std::max(b_start, p_start);
      overlap              = std::max(overlap, (ssize_t)0);
      block.fileCount      = overlap;
      size_t element_nodes = block.nodesPerEntity;
      int    zone          = block.zone_;
      int    section       = block.section_;

      // Get the connectivity (raw) for this portion of elements...
      std::vector<cgsize_t> connectivity(overlap * element_nodes);
      INT                   blk_start = std::max(b_start, p_start) - b_start + 1;
      INT                   blk_end   = blk_start + overlap - 1;
      blk_start                       = blk_start < 0 ? 0 : blk_start;
      blk_end                         = blk_end < 0 ? 0 : blk_end;
#if IOSS_DEBUG_OUTPUT
      OUTPUT << "Processor " << m_decomposition.m_processor << " has " << overlap
             << " elements on element block " << block.name() << "\t(" << blk_start << " to "
             << blk_end << ")\n";
#endif
      block.fileSectionOffset = blk_start;
#if CG_BUILD_PARALLEL
      cgp_elements_read_data(filePtr, base, zone, section, blk_start, blk_end, TOPTR(connectivity));
#else
      cg_elements_partial_read(filePtr, base, zone, section, blk_start, blk_end,
                               TOPTR(connectivity), nullptr);
#endif
      size_t el          = 0;
      INT    zone_offset = block.zoneNodeOffset;

      for (ssize_t elem = 0; elem < overlap; elem++) {
        decomposition.m_pointer.push_back(decomposition.m_adjacency.size());
        for (size_t k = 0; k < element_nodes; k++) {
          INT node = connectivity[el++] - 1 + zone_offset; // 0-based node
          decomposition.m_adjacency.push_back(node);
        }
      }
      sum += overlap * element_nodes;
    }
    decomposition.m_pointer.push_back(decomposition.m_adjacency.size());
  }

  template <typename INT> void DecompositionData<INT>::get_sideset_data(int filePtr)
  {
    int base = 1; // Only single base supported so far.

    // Get total length of sideset elemlists...
    size_t elemlist_size = 0;
    for (auto &sset : m_sideSets) {
      elemlist_size += sset.file_count();
    }

    // Calculate the max "buffer" size usable for storing sideset
    // elemlists. This is basically the space used to store the file
    // decomposition nodal coordinates. The "decomp_node_count()/2*2" is to
    // equalize the decomp_node_count() among processors since some procs have 1
    // more node than others. For small models, assume we can handle
    // at least 10000 nodes.
    //    size_t max_size = std::max(10000, (decomp_node_count() / 2) * 2 * 3 *sizeof(double) /
    //    sizeof(cgsize_t));

    bool subsetting = false; // elemlist_size > max_size;

    if (subsetting) {
      assert(1 == 0);
    }
    else {
      // Can handle reading all sideset elem lists on a single
      // processor simultaneously.
      std::vector<cgsize_t> elemlist(elemlist_size);

      size_t offset = 0;
#if !CG_BUILD_PARALLEL
      // Read the elemlists on root processor.
      // If CG_BUILD_PARALLEL, then all processors must call the cg_elements_read
      int root = 0; // Root processor that reads all sideset bulk data (nodelists)
      if (m_decomposition.m_processor == root) {
#endif
        for (auto &sset : m_sideSets) {

          // TODO? Possibly rewrite using cgi_read_int_data so can skip reading element connectivity
          int                   nodes_per_face = 4; // FIXME: sb->topology()->number_nodes();
          std::vector<cgsize_t> elements(nodes_per_face *
                                         sset.file_count()); // Not needed, but can't skip

          // We get:
          // *  num_to_get parent elements,
          // *  num_to_get zeros (other parent element for face, but on boundary so 0)
          // *  num_to_get face_on_element
          // *  num_to_get zeros (face on other parent element)
          std::vector<cgsize_t> parent(4 * sset.file_count());

          int ierr = cg_elements_read(filePtr, base, sset.zone(), sset.section(), TOPTR(elements),
                                      TOPTR(parent));
          if (ierr < 0) {
            Utils::cgns_error(filePtr, __FILE__, __func__, __LINE__, m_decomposition.m_processor);
          }

          // Move from 'parent' to 'elementlist'
          size_t zone_element_id_offset = m_zones[sset.zone()].m_elementOffset;
          for (size_t i = 0; i < sset.file_count(); i++) {
            elemlist[offset++] = parent[i] + zone_element_id_offset;
          }
        }
        assert(offset == elemlist_size);
#if !CG_BUILD_PARALLEL
      }

      // Broadcast this data to all other processors...
      MPI_Bcast(TOPTR(elemlist), sizeof(cgsize_t) * elemlist.size(), MPI_BYTE, root,
                m_decomposition.m_comm);
#endif

      // Each processor now has a complete list of all elems in all
      // sidesets.
      // Determine which of these are owned by the current
      // processor...
      {
        offset = 0;
        for (auto &sset : m_sideSets) {
          size_t ss_beg = offset;
          size_t ss_end = ss_beg + sset.file_count();

          for (size_t n = ss_beg; n < ss_end; n++) {
            cgsize_t elem = elemlist[n];
            // See if elem owned by this processor...
            if (i_own_elem(elem)) {
              // Save elem in this processors elemlist for this set.
              // The saved data is this elems location in the global
              // elemlist for this set.
              sset.entitylist_map.push_back(n - offset);
            }
          }
          offset = ss_end;
        }
      }

      // Each processor knows how many of the sideset elems it owns;
      // broadcast that information (the count) to the other
      // processors. The first processor with non-zero elem count is
      // the "root" for this sideset.
      {
        std::vector<int> has_elems_local(m_sideSets.size());
        for (size_t i = 0; i < m_sideSets.size(); i++) {
          has_elems_local[i] = m_sideSets[i].entitylist_map.empty() ? 0 : 1;
        }

        std::vector<int> has_elems(m_sideSets.size() * m_decomposition.m_processorCount);
        MPI_Allgather(TOPTR(has_elems_local), has_elems_local.size(), MPI_INT, TOPTR(has_elems),
                      has_elems_local.size(), MPI_INT, m_decomposition.m_comm);

        for (size_t i = 0; i < m_sideSets.size(); i++) {
          m_sideSets[i].hasEntities.resize(m_decomposition.m_processorCount);
          m_sideSets[i].root_ = m_decomposition.m_processorCount;
          for (int p = 0; p < m_decomposition.m_processorCount; p++) {
            if (p < m_sideSets[i].root_ && has_elems[p * m_sideSets.size() + i] != 0) {
              m_sideSets[i].root_ = p;
            }
            m_sideSets[i].hasEntities[p] = has_elems[p * m_sideSets.size() + i];
          }
        }
      }
    }
  }

  template <typename INT>
  void DecompositionData<INT>::get_file_node_coordinates(int filePtr, int direction,
                                                         double *data) const
  {
#if !CG_BUILD_PARALLEL
    const std::string coord_name[] = {"CoordinateX", "CoordinateY", "CoordinateZ"};
#endif

    int      base        = 1; // Only single base supported so far.
    cgsize_t beg         = 0;
    cgsize_t end         = 0;
    cgsize_t offset      = 0;
    cgsize_t node_count  = decomp_node_count();
    cgsize_t node_offset = decomp_node_offset();

    int num_zones = (int)m_zones.size() - 1;
    for (int zone = 1; zone <= num_zones; zone++) {
      end += m_zones[zone].m_nodeCount;

      cgsize_t start  = std::max(node_offset, beg);
      cgsize_t finish = std::min(end, node_offset + node_count);
      cgsize_t count  = (finish > start) ? finish - start : 0;

      // Now adjust start for 1-based node numbering and the start of this zone...
      start  = start - beg + 1;
      finish = finish - beg;
      if (count == 0) {
        start  = 0;
        finish = 0;
      }
#if defined(IOSS_DEBUG_OUTPUT)
      OUTPUT << m_decomposition.m_processor << ": reading " << count << " nodes from zone " << zone
             << " starting at " << start << " with an offset of " << offset << " ending at "
             << finish << "\n";
#endif
#if CG_BUILD_PARALLEL
      double *coords = nullptr;
      if (count > 0) {
        coords = &data[offset];
      }
      int ierr = cgp_coord_read_data(filePtr, base, zone, direction + 1, &start, &finish, coords);
#else
      int ierr = 0;
      if (count > 0) {
        ierr = cg_coord_read(filePtr, base, zone, coord_name[direction].c_str(), CG_RealDouble,
                             &start, &finish, &data[offset]);
      }
#endif
      if (ierr < 0) {
        Utils::cgns_error(filePtr, __FILE__, __func__, __LINE__, m_decomposition.m_processor);
      }
      offset += count;
      beg = end;
    }
  }

  template <typename INT>
  void DecompositionData<INT>::get_node_coordinates(int filePtr, double *ioss_data,
                                                    const Ioss::Field &field) const
  {
    std::vector<double> tmp(decomp_node_count());
    if (field.get_name() == "mesh_model_coordinates_x") {
      get_file_node_coordinates(filePtr, 0, TOPTR(tmp));
      communicate_node_data(TOPTR(tmp), ioss_data, 1);
    }

    else if (field.get_name() == "mesh_model_coordinates_y") {
      get_file_node_coordinates(filePtr, 1, TOPTR(tmp));
      communicate_node_data(TOPTR(tmp), ioss_data, 1);
    }

    else if (field.get_name() == "mesh_model_coordinates_z") {
      get_file_node_coordinates(filePtr, 2, TOPTR(tmp));
      communicate_node_data(TOPTR(tmp), ioss_data, 1);
    }

    else if (field.get_name() == "mesh_model_coordinates") {
      // Data required by upper classes store x0, y0, z0, ... xn,
      // yn, zn. Data stored in cgns file is x0, ..., xn, y0,
      // ..., yn, z0, ..., zn so we have to allocate some scratch
      // memory to read in the data and then map into supplied
      // 'data'

      std::vector<double> ioss_tmp(ioss_node_count());

      // This implementation trades off extra communication for
      // reduced memory overhead.
      // * This method uses 'ioss_node_count' extra memory; 3
      // reads; and 3 communicate_node_data calls.
      //
      // * Other method uses 6*ioss_node_count extra memory; 3 reads;
      // and 1 communicate_node_data call.
      //
      for (int d = 0; d < m_decomposition.m_spatialDimension; d++) {
        get_file_node_coordinates(filePtr, d, TOPTR(tmp));
        communicate_node_data(TOPTR(tmp), TOPTR(ioss_tmp), 1);

        size_t index = d;
        for (size_t i = 0; i < ioss_node_count(); i++) {
          ioss_data[index] = ioss_tmp[i];
          index += m_decomposition.m_spatialDimension;
        }
      }
    }
  }

  template void DecompositionData<int>::get_sideset_element_side(
      int filePtr, const Ioss::SetDecompositionData &sset, int *data) const;
  template void DecompositionData<int64_t>::get_sideset_element_side(
      int filePtr, const Ioss::SetDecompositionData &sset, int64_t *data) const;
  template <typename INT>
  void DecompositionData<INT>::get_sideset_element_side(int                               filePtr,
                                                        const Ioss::SetDecompositionData &sset,
                                                        INT *ioss_data) const
  {
    std::vector<INT> element_side;
#if !CG_BUILD_PARALLEL
    // Cannot do this due to parallel read -- hangs in hdf5
    if (m_decomposition.m_processor == sset.root_) {
#endif
      int base = 1;

      int                   nodes_per_face = 4; // FIXME: sb->topology()->number_nodes();
      std::vector<cgsize_t> nodes(nodes_per_face * sset.file_count());

      // TODO? Possibly rewrite using cgi_read_int_data so can skip reading element connectivity

      // We get:
      // *  num_to_get parent elements,
      // *  num_to_get zeros (other parent element for face, but on boundary so 0)
      // *  num_to_get face_on_element
      // *  num_to_get zeros (face on other parent element)
      std::vector<cgsize_t> parent(4 * sset.file_count());

      int ierr =
          cg_elements_read(filePtr, base, sset.zone(), sset.section(), TOPTR(nodes), TOPTR(parent));
      // Get rid of 'nodes' list -- not used.
      nodes.resize(0);
      nodes.shrink_to_fit();

      if (ierr < 0) {
        Utils::cgns_error(filePtr, __FILE__, __func__, __LINE__, m_decomposition.m_processor);
      }

      // Move from 'parent' to 'element_side' and interleave. element, side, element, side, ...
      element_side.reserve(sset.file_count() * 2);
      size_t zone_element_id_offset = m_zones[sset.zone()].m_elementOffset;
      for (size_t i = 0; i < sset.file_count(); i++) {
        element_side.push_back(parent[0 * sset.file_count() + i] + zone_element_id_offset);
        element_side.push_back(parent[2 * sset.file_count() + i]);
      }
#if !CG_BUILD_PARALLEL
    }
#endif
    // The above was all on root processor for this side set, now need to send data to other
    // processors that own any of the elements in the sideset.
    communicate_set_data(TOPTR(element_side), ioss_data, sset, 2);
  }

  template void DecompositionData<int>::get_block_connectivity(int filePtr, int *data,
                                                               int blk_seq) const;
  template void DecompositionData<int64_t>::get_block_connectivity(int filePtr, int64_t *data,
                                                                   int blk_seq) const;

  template <typename INT>
  void DecompositionData<INT>::get_block_connectivity(int filePtr, INT *data, int blk_seq) const
  {
    auto                  blk = m_elementBlocks[blk_seq];
    std::vector<cgsize_t> file_conn(blk.file_count() * blk.nodesPerEntity);
    int                   base = 1;
#if CG_BUILD_PARALLEL
    cgp_elements_read_data(filePtr, base, blk.zone(), blk.section(), blk.fileSectionOffset,
                           blk.fileSectionOffset + blk.file_count() - 1, TOPTR(file_conn));
#else
    cg_elements_partial_read(filePtr, base, blk.zone(), blk.section(), blk.fileSectionOffset,
                             blk.fileSectionOffset + blk.file_count() - 1, TOPTR(file_conn),
                             nullptr);
#endif
    // Map from zone-local node numbers to global implicit
    for (auto &node : file_conn) {
      node += blk.zoneNodeOffset;
    }

    if (!m_zoneSharedMap.empty()) {
      for (auto &node : file_conn) {
        std::unordered_map<cgsize_t, cgsize_t>::const_iterator alias =
            m_zoneSharedMap.find(node - 1);
        if (alias != m_zoneSharedMap.end()) {
          node = (*alias).second + 1;
        }
      }
    }

    communicate_block_data(TOPTR(file_conn), data, blk, blk.nodesPerEntity);
  }

  DecompositionDataBase::~DecompositionDataBase()
  {
    for (auto &zone : m_structuredZones) {
      delete zone;
    }
  }

  template void DecompositionDataBase::communicate_node_data(int *file_data, int *ioss_data,
                                                             size_t comp_count) const;
  template void DecompositionDataBase::communicate_node_data(int64_t *file_data, int64_t *ioss_data,
                                                             size_t comp_count) const;
  template void DecompositionDataBase::communicate_node_data(double *file_data, double *ioss_data,
                                                             size_t comp_count) const;

  template <typename T>
  void DecompositionDataBase::communicate_node_data(T *file_data, T *ioss_data,
                                                    size_t comp_count) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int> *>(this);
      Ioss::Utils::check_dynamic_cast(this32);
      this32->communicate_node_data(file_data, ioss_data, comp_count);
    }
    else {
      const DecompositionData<int64_t> *this64 =
          dynamic_cast<const DecompositionData<int64_t> *>(this);
      Ioss::Utils::check_dynamic_cast(this64);
      this64->communicate_node_data(file_data, ioss_data, comp_count);
    }
  }

  template void DecompositionDataBase::communicate_element_data(int *file_data, int *ioss_data,
                                                                size_t comp_count) const;
  template void DecompositionDataBase::communicate_element_data(int64_t *file_data,
                                                                int64_t *ioss_data,
                                                                size_t   comp_count) const;
  template void DecompositionDataBase::communicate_element_data(double *file_data,
                                                                double *ioss_data,
                                                                size_t  comp_count) const;

  template <typename T>
  void DecompositionDataBase::communicate_element_data(T *file_data, T *ioss_data,
                                                       size_t comp_count) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int> *>(this);
      Ioss::Utils::check_dynamic_cast(this32);
      this32->communicate_element_data(file_data, ioss_data, comp_count);
    }
    else {
      const DecompositionData<int64_t> *this64 =
          dynamic_cast<const DecompositionData<int64_t> *>(this);
      Ioss::Utils::check_dynamic_cast(this64);
      this64->communicate_element_data(file_data, ioss_data, comp_count);
    }
  }

  void DecompositionDataBase::get_node_entity_proc_data(void *                    entity_proc,
                                                        const Ioss::MapContainer &node_map,
                                                        bool                      do_map) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int> *>(this);
      Ioss::Utils::check_dynamic_cast(this32);
      this32->m_decomposition.get_node_entity_proc_data((int *)entity_proc, node_map, do_map);
    }
    else {
      const DecompositionData<int64_t> *this64 =
          dynamic_cast<const DecompositionData<int64_t> *>(this);
      Ioss::Utils::check_dynamic_cast(this64);
      this64->m_decomposition.get_node_entity_proc_data((int64_t *)entity_proc, node_map, do_map);
    }
  }

  void DecompositionDataBase::get_block_connectivity(int filePtr, void *data, int blk_seq) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int> *>(this);
      Ioss::Utils::check_dynamic_cast(this32);
      this32->get_block_connectivity(filePtr, (int *)data, blk_seq);
    }
    else {
      const DecompositionData<int64_t> *this64 =
          dynamic_cast<const DecompositionData<int64_t> *>(this);
      Ioss::Utils::check_dynamic_cast(this64);
      this64->get_block_connectivity(filePtr, (int64_t *)data, blk_seq);
    }
  }

  void DecompositionDataBase::get_sideset_element_side(int                               filePtr,
                                                       const Ioss::SetDecompositionData &sset,
                                                       void *                            data) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int> *>(this);
      Ioss::Utils::check_dynamic_cast(this32);
      this32->get_sideset_element_side(filePtr, sset, (int *)data);
    }
    else {
      const DecompositionData<int64_t> *this64 =
          dynamic_cast<const DecompositionData<int64_t> *>(this);
      Ioss::Utils::check_dynamic_cast(this64);
      this64->get_sideset_element_side(filePtr, sset, (int64_t *)data);
    }
  }
}

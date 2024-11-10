#include <gtest/gtest.h>
#include <stk_unit_test_utils/getOption.h>
#include <string>
#include <memory>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include "stk_mesh/baseImpl/EquivalentEntityBlocks.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_io/FillMesh.hpp"
#include <stk_balance/internal/StkBalanceUtils.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>

#include <stk_search/SearchMethod.hpp>
#include <stk_search/CoarseSearch.hpp>

#include <stk_balance/internal/privateDeclarations.hpp>
#include "stk_mesh/base/SkinMeshUtil.hpp"

#include <stk_balance/balance.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include "stk_balance/internal/LastStepFieldWriter.hpp"

#include "stk_balance/search_tolerance_algs/SecondShortestEdgeFaceSearchTolerance.hpp"
#include "test_utils/StkBalanceRunner.hpp"

namespace
{
using stk::unit_test_util::build_mesh;

TEST(Stkbalance, DISABLED_Ticket15830)
{
  std::string filename = stk::unit_test_util::get_option("-i", "rs1.rsout");

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;

  EXPECT_NO_THROW(stk::io::fill_mesh_with_auto_decomp(filename, bulk));
}

void getMaxMinForNodes(stk::mesh::BulkData& bulk, stk::mesh::EntityVector& nodes, std::vector<double>& minCoord, std::vector<double>& maxCoord, const stk::mesh::FieldBase & coordField)
{
  for(stk::mesh::Entity node : nodes)
  {
    double *coord = static_cast<double*>(stk::mesh::field_data(coordField, node));
    for(unsigned j=0;j<bulk.mesh_meta_data().spatial_dimension();++j)
    {
      minCoord[j] = std::min(minCoord[j], coord[j]);
      maxCoord[j] = std::max(maxCoord[j], coord[j]);
    }
  }
}

stk::balance::internal::SearchBoxIdentProcs fill_bounding_boxes(const std::vector<stk::mesh::SideSetEntry>& skinnedSideSet,
                                                                stk::mesh::BulkData& bulk,
                                                                std::vector<double>& coordMinOnProc,
                                                                std::vector<double>& coordMaxOnProc,
                                                                double searchTol)
{
  stk::balance::internal::SearchBoxIdentProcs local_domain;
  const stk::mesh::FieldBase * coordField = bulk.mesh_meta_data().get_field(stk::topology::NODE_RANK, "Coordinates" );

  int boxCounter = 0;
  for(stk::mesh::SideSetEntry sidesetEntry : skinnedSideSet)
  {
    stk::mesh::Entity sidesetElement = sidesetEntry.element;
    stk::mesh::ConnectivityOrdinal sidesetSide = sidesetEntry.side;
    stk::mesh::EntityVector sideNodes;
    stk::mesh::get_subcell_nodes(bulk, sidesetElement, bulk.mesh_meta_data().side_rank(), sidesetSide, sideNodes);
    std::vector<double> minCoord(3, std::numeric_limits<double>::max());
    std::vector<double> maxCoord(3, std::numeric_limits<double>::lowest());
    getMaxMinForNodes(bulk, sideNodes, minCoord, maxCoord, *coordField);
    for(size_t i = 0; i < minCoord.size(); ++i)
    {
      coordMinOnProc[i] = std::min(coordMinOnProc[i], minCoord[i]);
      coordMaxOnProc[i] = std::max(coordMaxOnProc[i], maxCoord[i]);
    }
    stk::balance::internal::StkBox box(minCoord[0] - searchTol,
        minCoord[1] - searchTol,
        minCoord[2] - searchTol,
        maxCoord[0] + searchTol,
        maxCoord[1] + searchTol,
        maxCoord[2] + searchTol);
    stk::balance::internal::SearchIdentProc id(boxCounter, bulk.parallel_rank());
    boxCounter++;
    local_domain.push_back(std::make_pair(box, id));
  }
  return local_domain;
}

void write_metrics_for_overlapping_BB(int myId, const stk::balance::internal::SearchBoxIdentProcs &local_domain, const std::vector<double>& coordMinOnProc, const std::vector<double>& coordMaxOnProc)
{
  stk::balance::internal::SearchBoxIdentProcs local_range = local_domain;

  stk::balance::internal::SearchElemPairs searchResults;
  stk::search::coarse_search(local_domain, local_range, stk::search::KDTREE, MPI_COMM_WORLD, searchResults);

  std::ostringstream os;
  double dx = coordMaxOnProc[0] - coordMinOnProc[0];
  double dy = coordMaxOnProc[1] - coordMinOnProc[1];
  double dz = coordMaxOnProc[2] - coordMinOnProc[2];

  os << "For processor " << myId << ", Number of interactions: " << searchResults.size() << ", and num boxes = " << local_range.size() <<  ", and volume = " << dx*dy*dz;

  size_t counter = 0;
  for(size_t i=0;i<searchResults.size();++i)
  {
    if( searchResults[i].first.proc() != searchResults[i].second.proc() )
    {
      if(searchResults[i].first.id() != searchResults[i].second.id())
      {
        counter++;
      }
    }
  }

  os << ", Number unique overlaps of faces between procs:, " << counter;

  size_t globalCounts = 0;
  stk::all_reduce_sum(MPI_COMM_WORLD, &counter, &globalCounts, 1);
  os << ", Global counts ," << globalCounts << std::endl;

  std::cerr << os.str();
}

// This is a useful diagnostic tool in unit test form, so it doesn't
// actually check anything.
TEST(Stkbalance, NumOverlappingBB)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 3) return;

  const std::string dummyFileName("ARefLA.e");
  std::string filename = stk::unit_test_util::get_option("-i", dummyFileName);

  std::vector<double> coordMinOnProc(3, std::numeric_limits<double>::max());
  std::vector<double> coordMaxOnProc(3, std::numeric_limits<double>::lowest());

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::io::fill_mesh(filename, bulk);

  std::vector<stk::mesh::SideSetEntry> skinnedSideSet = stk::mesh::SkinMeshUtil::get_skinned_sideset(bulk, bulk.mesh_meta_data().locally_owned_part());

  double searchTol = 1e-4;
  stk::balance::internal::SearchBoxIdentProcs local_domain = fill_bounding_boxes(skinnedSideSet, bulk, coordMinOnProc, coordMaxOnProc, searchTol);

  int myId = bulk.parallel_rank();
  write_metrics_for_overlapping_BB(myId, local_domain, coordMinOnProc, coordMaxOnProc);
}

void addComponentToSet(std::set<stk::mesh::Entity> &component, stk::mesh::ElemElemGraph& graph, stk::mesh::Entity startingElement)
{
  bool done = false;
  stk::mesh::EntityVector elementsToConsider;
  stk::mesh::EntityVector newElements;
  newElements.push_back(startingElement);

  while(!done)
  {
    elementsToConsider = newElements;
    newElements.clear();
    for(stk::mesh::Entity entity : elementsToConsider)
    {
      size_t numConn = graph.get_num_connected_elems(entity);
      for(size_t i=0;i<numConn;++i)
      {
        if(graph.is_connected_elem_locally_owned(entity, i))
        {
          stk::mesh::impl::ElementViaSidePair elemSide = graph.get_connected_element_and_via_side(entity, i);
          stk::mesh::Entity otherElem = elemSide.element;
          auto iter=component.find(otherElem);
          if(iter==component.end())
          {
            newElements.push_back(otherElem);
            component.insert(otherElem);
          }
        }
      }
    }

    done = newElements.size()==0;
  }
}

// This is a useful diagnostic tool in unit test form, so it doesn't
// actually check anything.
TEST(Stkbalance, modifyMeshIfNeeded)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 3) return;

  std::string filename = stk::unit_test_util::get_option("-i", "ARefLA.e");

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;

  EXPECT_NO_THROW(stk::io::fill_mesh(filename, bulk));

  stk::mesh::EntityVector elements;
  const bool sortById = true;
  stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, bulk.mesh_meta_data().locally_owned_part(), elements, sortById);

  std::vector<std::set<stk::mesh::Entity>> components;

  stk::mesh::ElemElemGraph &graph = bulk.get_face_adjacent_element_graph();

  stk::mesh::Entity startingElement = elements[0];
  size_t memoryIndex = 1;

  bool done = false;
  while(!done)
  {
    std::set<stk::mesh::Entity> tmp;
    tmp.insert(startingElement);
    components.push_back(tmp);

    addComponentToSet(components.back(), graph, startingElement);

    size_t totalSize = 0;
    for(size_t i=0;i<components.size();++i)
      totalSize += components[i].size();

    done = totalSize == elements.size();
    if(!done)
    {
      for(size_t i=memoryIndex;i<elements.size();++i)
      {
        stk::mesh::Entity e = elements[i];
        bool uniqueElement = true;
        for(size_t j=0;j<components.size();++j)
        {
          auto iter=components[j].find(e);
          if(iter!=components[j].end())
          {
            uniqueElement = false;
            break;
          }
        }
        if(uniqueElement)
        {
          startingElement = elements[i];
          memoryIndex = i+1;
          break;
        }
      }
    }
  }

  done = false;
  std::vector<int> indicesToErase;
  while(!done)
  {
    indicesToErase.clear();
    for(size_t i=1;i<components.size();++i)
    {
      stk::mesh::EntityVector diff(elements.size());
      auto iter=std::set_intersection(components[i].begin(), components[i].end(), components[i-1].begin(), components[i-1].end(), diff.begin());
      size_t numDiff = iter-diff.begin();
      if(numDiff > 0)
      {
        components[i-1].insert(components[i].begin(), components[i].end());
        indicesToErase.push_back(i);
      }
    }

    done = indicesToErase.size() == 0;

    if(!done)
    {
      std::sort(indicesToErase.begin(), indicesToErase.end(), std::greater<int>());
      for(int index : indicesToErase)
        components.erase(components.begin()+index);
    }
  }

  std::ostringstream os;
  os << "Processor " << bulk.parallel_rank() << " has " << components.size() << " components.\n";
  std::cerr << os.str();
}

// This is a useful diagnostic tool in unit test form, so it doesn't
// actually check anything.
TEST(Stkbalance, checkForDegenerateElements)
{
  std::string filename = stk::unit_test_util::get_option("-i", "ZDZ.e");

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  EXPECT_NO_THROW(stk::io::fill_mesh_with_auto_decomp(filename, bulk));

  stk::mesh::EntityVector elements;
  const bool sortById = true;
  stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, meta.locally_owned_part(), elements, sortById);

  std::ostringstream os;
  std::vector<stk::mesh::PartOrdinal> ords;

  for(stk::mesh::Entity element : elements)
  {
    stk::mesh::EntityVector nodes(bulk.begin_nodes(element), bulk.end_nodes(element));
    size_t sizeBefore = nodes.size();
    stk::util::sort_and_unique(nodes);
    if(nodes.size() != sizeBefore)
    {
      os << "[" << bulk.parallel_rank() << "] for element " << bulk.identifier(element) << " has degeneracy in element block ";
      stk::mesh::impl::get_element_block_part_ordinals(element, bulk, ords);
      os << meta.get_parts()[ords[0]]->name() << " with " << nodes.size() << " unique nodes\n";
    }
  }

  std::cerr << os.str();
}

}

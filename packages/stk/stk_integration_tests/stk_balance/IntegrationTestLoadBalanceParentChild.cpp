#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_unit_test_utils/getOption.h"
#include "stk_mesh/base/GetEntities.hpp"
#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>

namespace {

// This is meant to *minimally* represent the kind of parent/child relationships
// among elements refined during CDFEM simulations using Krino and adaptivity
// using Percept.
//

//==============================================================================
class ParentChildManager
{
public:
  ParentChildManager() {};
  virtual ~ParentChildManager() {};

  virtual void add_child_element(stk::mesh::Entity parentElement, stk::mesh::Entity childElement) = 0;
  virtual stk::mesh::EntityVector get_child_elements(stk::mesh::Entity parentElement) const = 0;
  virtual stk::mesh::Entity get_parent_element(stk::mesh::Entity childElement) const = 0;
  virtual bool has_children(stk::mesh::Entity parentElement) const = 0;
  virtual void move_related_entities_with_parent_element(stk::balance::DecompositionChangeList & decomp, stk::mesh::Entity parentElement, const int destination)
  {
    if (has_children(parentElement)) {
      set_child_element_destination_from_parent_element(decomp, parentElement, destination);
    }
  }

protected:
  void set_entity_destination(stk::balance::DecompositionChangeList & decomp, stk::mesh::Entity entity, const int destination) const
  {
    EXPECT_TRUE(decomp.get_bulk().is_valid(entity));
    decomp.set_entity_destination(entity, destination);
  }

  void set_child_element_destination_from_parent_element(stk::balance::DecompositionChangeList & decomp, stk::mesh::Entity parentElement, const int destination) const
  {
    stk::mesh::EntityVector childElements = get_child_elements(parentElement);
    for (const stk::mesh::Entity & childElement : childElements) {
      set_entity_destination(decomp, childElement, destination);
    }
  }
};

//==============================================================================
class MapParentChildManager : public ParentChildManager
{
public:
  MapParentChildManager() {};
  virtual ~MapParentChildManager() {};

  void add_child_element(stk::mesh::Entity parentElement, stk::mesh::Entity childElement)
  {
    m_childMap[parentElement].push_back(childElement);
    m_parentMap[childElement] = parentElement;
  }

  stk::mesh::EntityVector get_child_elements(stk::mesh::Entity parentElement) const
  {
    return m_childMap.at(parentElement);
  }

  stk::mesh::Entity get_parent_element(stk::mesh::Entity childElement) const
  {
    return m_parentMap.at(childElement);
  }

  bool has_children(stk::mesh::Entity parentElement) const
  {
    return (m_childMap.find(parentElement) != m_childMap.end());
  }

private:
  std::map<stk::mesh::Entity, stk::mesh::EntityVector> m_childMap;
  std::map<stk::mesh::Entity, stk::mesh::Entity> m_parentMap;
};

//==============================================================================
class RelationParentChildManager : public ParentChildManager
{
public:
  RelationParentChildManager(stk::mesh::BulkData & bulkData)
    : m_stkMeshBulkData(bulkData) {};
  virtual ~RelationParentChildManager() {};

  virtual void move_related_entities_with_parent_element(stk::balance::DecompositionChangeList & decomp, stk::mesh::Entity parentElement, const int destination)
  {
    ParentChildManager::move_related_entities_with_parent_element(decomp, parentElement, destination);
    stk::mesh::Entity familyTree = get_family_tree(parentElement);

    if (decomp.get_bulk().is_valid(familyTree)) {
      set_entity_destination(decomp, familyTree, destination);
    }
  }

  void add_child_element(stk::mesh::Entity parentElement, stk::mesh::Entity childElement)
  {
    stk::mesh::Entity familyTree = get_or_create_family_tree(parentElement);
    unsigned newOrdinal = count_relations(familyTree, stk::topology::ELEMENT_RANK);
    m_stkMeshBulkData.declare_relation(familyTree, childElement, newOrdinal);
  }

  stk::mesh::EntityVector get_child_elements(stk::mesh::Entity parentElement) const
  {
    stk::mesh::EntityVector childElements;

    stk::mesh::Entity familyTree = get_family_tree(parentElement);
    if (familyTree != stk::mesh::Entity::InvalidEntity) {
      gather_child_elements_for_parent(parentElement, familyTree, childElements);
    }

    return childElements;
  }

  stk::mesh::Entity get_parent_element(stk::mesh::Entity childElement) const
  {
    stk::mesh::Entity familyTree = get_family_tree(childElement);
    STK_ThrowRequire(familyTree != stk::mesh::Entity::InvalidEntity);
    return get_parent_from_family_tree(familyTree);
  }

  bool has_children(stk::mesh::Entity parentElement) const
  {
    return (get_child_elements(parentElement).size() > 0);
  }

private:
  unsigned count_relations(stk::mesh::Entity entity, stk::topology::rank_t rank) const
  {
    stk::mesh::Entity const * entityStart = m_stkMeshBulkData.begin(entity, rank);
    stk::mesh::Entity const * entityEnd   = m_stkMeshBulkData.end  (entity, rank);

    const unsigned numRelatedEntities = std::distance(entityStart, entityEnd);
    return numRelatedEntities;
  }

  stk::mesh::Entity get_family_tree(stk::mesh::Entity entity) const
  {
    stk::mesh::Entity familyTree;

    const unsigned numConstraints = count_relations(entity, stk::topology::CONSTRAINT_RANK);
    if (1 == numConstraints) {
      familyTree = *(m_stkMeshBulkData.begin(entity, stk::topology::CONSTRAINT_RANK));
    }
    return familyTree;
  }

  stk::mesh::Entity create_family_tree(stk::mesh::Entity entity) const
  {
    stk::mesh::EntityId id = m_stkMeshBulkData.identifier(entity);
    stk::mesh::Entity familyTree = m_stkMeshBulkData.declare_constraint(id);
    m_stkMeshBulkData.declare_relation(familyTree, entity, 0);
    return familyTree;
  }

  stk::mesh::Entity get_or_create_family_tree(stk::mesh::Entity entity)
  {
    stk::mesh::Entity familyTree = get_family_tree(entity);
    if (familyTree == stk::mesh::Entity::InvalidEntity) {
      familyTree = create_family_tree(entity);
    }
    return familyTree;
  }

  void gather_child_elements_for_parent(stk::mesh::Entity parentElement, stk::mesh::Entity familyTree, stk::mesh::EntityVector& childElements) const {
    const stk::mesh::Entity* entity    = m_stkMeshBulkData.begin(familyTree, stk::topology::ELEMENT_RANK);
    const stk::mesh::Entity* entityEnd = m_stkMeshBulkData.end  (familyTree, stk::topology::ELEMENT_RANK);
    for ( ; entity != entityEnd; ++entity) {
      if (*entity != parentElement) {
        childElements.push_back(*entity);
      }
    }
  }

  stk::mesh::Entity get_parent_from_family_tree(stk::mesh::Entity familyTree) const {
    const stk::mesh::Entity* entityStart = m_stkMeshBulkData.begin(familyTree, stk::topology::ELEMENT_RANK);
    const stk::mesh::Entity* entityEnd   = m_stkMeshBulkData.end  (familyTree, stk::topology::ELEMENT_RANK);
    STK_ThrowRequire(entityStart != entityEnd);
    return (*entityStart);
  }

private:
  stk::mesh::BulkData & m_stkMeshBulkData;
};


//==============================================================================
class StkFieldWeightRebalance : public stk::balance::FieldVertexWeightSettings
{
public:
  StkFieldWeightRebalance(stk::mesh::BulkData &stkMeshBulkData,
                          const stk::balance::DoubleFieldType &weightField)
    : FieldVertexWeightSettings(stkMeshBulkData, weightField)
  {
  }

  virtual bool allowModificationOfVertexWeightsForSmallMeshes() const { return false; }

  virtual ~StkFieldWeightRebalance() = default;

  virtual std::string getCoordinateFieldName() const { return "model_coordinates"; }

protected:
  StkFieldWeightRebalance() = delete;
  StkFieldWeightRebalance(const StkFieldWeightRebalance&) = delete;
  StkFieldWeightRebalance& operator=(const StkFieldWeightRebalance&) = delete;
};

//==============================================================================
class StkParentRebalance : public StkFieldWeightRebalance
{
public:
  StkParentRebalance(ParentChildManager & parentChildManager,
                     const stk::mesh::Selector & selector,
                     stk::mesh::BulkData &stkMeshBulkData,
                     const stk::balance::DoubleFieldType &weightField)
    : StkFieldWeightRebalance(stkMeshBulkData, weightField),
      m_parentChildManager(parentChildManager),
      m_parentElementSelector(selector)
  { }

  virtual ~StkParentRebalance() = default;

  virtual std::string getCoordinateFieldName() const { return "model_coordinates"; }

  virtual void modifyDecomposition(stk::balance::DecompositionChangeList & decomp) const
  {
    delete_child_elements_from_decomposition(decomp);
    move_related_entities_with_parent_element(decomp);
  };

  virtual bool shouldFixMechanisms() const { return false; }

private:

  stk::mesh::EntityVector get_elements_from_selector(stk::mesh::BulkData & stkMeshBulkData, stk::mesh::Selector selector) const
  {
    stk::mesh::EntityVector elements;
    const bool sortById = true;
    stk::mesh::get_entities(stkMeshBulkData, stk::topology::ELEM_RANK, selector, elements, sortById);
    return elements;
  }

  void move_related_entities_with_parent_element(stk::balance::DecompositionChangeList & decomp) const
  {
    stk::mesh::EntityProcVec changedEntities = decomp.get_all_partition_changes();
    for (stk::mesh::EntityProc & entityProcs : changedEntities) {
      m_parentChildManager.move_related_entities_with_parent_element(decomp, entityProcs.first, entityProcs.second);
    }
  }

  void delete_child_elements_from_decomposition(stk::balance::DecompositionChangeList & decomp) const
  {
    stk::mesh::Selector localChildElementSelector = (!m_parentElementSelector) & decomp.get_bulk().mesh_meta_data().locally_owned_part();
    stk::mesh::EntityVector childElements = get_elements_from_selector(decomp.get_bulk(), localChildElementSelector);
    for (const auto & childElement : childElements) {
      decomp.delete_entity(childElement);
    }
  }

protected:
  StkParentRebalance() = delete;
  StkParentRebalance(const StkParentRebalance&) = delete;
  StkParentRebalance& operator=(const StkParentRebalance&) = delete;

  ParentChildManager & m_parentChildManager;
  const stk::mesh::Selector & m_parentElementSelector;
};

//
//
//   Base mesh:                      | <-- Level set boundary requiring element cut
//      [4.0]             [3.0]      |      [6.0]             [8.1]            [10.1]
//        O-----------------O-----------------O-----------------O-----------------O
//        |\                |\       |        |\                |\                |
//        |  \        2.0   |  \     | 4.0    |  \       6.1    |  \       8.1    |
//        |    \            |    \   |        |    \            |    \            |
//        |      \          |      \ |        |      \          |      \          |
//        |        \        |        \        |        \        |        \        |
//        |          \      |        | \      |          \      |          \      |
//        |            \    |        |   \    |            \    |            \    |
//        |    1.0       \  |    3.0 |     \  |    5.1       \  |   7.1        \  |
//        |                \|        |       \|                \|                \|
//        O-----------------O-----------------O-----------------O-----------------O
//      [1.0]             [2.0]      |      [5.0]             [7.1]             [9.1]
//
//
//   Decomposed mesh with child elements:
//
//      [4.0]             [3.0]   [13.0]    [6.0]             [8.1]            [10.1]
//        O-----------------O========o========O-----------------O-----------------O
//        |\                ║\\(14.1)|(13.1) /║\                |\                |
//        |  \       2.0    ║  \\    | 4.1 /  ║  \       6.1    |  \       8.1    |
//        |    \            ║    \\  |   /    ║    \            |    \            |
//        |      \          ║      \\| /(12.1)║      \          |      \          |
//        |        \        ║(9.0)   o [12.0] ║        \        |        \        |
//        |          \      ║      / |\\      ║          \      |          \      |
//        |            \    ║    /   |  \\    ║            \    |            \    |
//        |    1.0       \  ║  / 3.0 |    \\  ║    5.1       \  |   7.1        \  |
//        |                \║/ (10.0)|(11.0)\\║                \|                \|
//        O-----------------O========o========O-----------------O-----------------O
//      [1.0]             [2.0]   [11.0]    [5.0]             [7.1]             [9.1]
//
//  Key:  X.Y  = ID.proc for base elements
//       (X.Y) = ID.proc for child elements
//       [X.Y] = ID.proc for nodes
//           O = Base mesh nodes
//           o = Refined mesh nodes (for child elements only)
//
//
class RebalanceParentChildMesh : public ::testing::Test
{
public:
  RebalanceParentChildMesh()
    : m_communicator(MPI_COMM_WORLD),
      m_metaData(nullptr),
      m_bulkData(nullptr),
      m_mapManager(),
      m_relationManager(nullptr),
      m_elementWeightField(nullptr),
      m_parentPart(nullptr)
  {
  }


  virtual ~RebalanceParentChildMesh()
  {
    delete m_relationManager;
  }

  void setup_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    allocate_bulk(auraOption);
    create_parts();
    register_fields();
    create_coarse_mesh();

    m_relationManager = new RelationParentChildManager(get_bulk());
  }

  void create_child_elements(ParentChildManager &parentChildManager)
  {
    std::vector<ParentChildRelationship> parentChildRels = create_child_mesh();

    get_bulk().modification_begin();
    connect_parents_and_children(parentChildRels, parentChildManager);
    get_bulk().modification_end();
  }

  MPI_Comm get_comm()
  {
    return m_communicator;
  }

  stk::mesh::MetaData& get_meta()
  {
    STK_ThrowRequireMsg(m_metaData!=nullptr, "Unit test error. Trying to get meta data before it has been initialized.");
    return *m_metaData;
  }

  stk::mesh::BulkData& get_bulk()
  {
    STK_ThrowRequireMsg(m_bulkData!=nullptr, "Unit test error. Trying to get bulk data before it has been initialized.");
    return *m_bulkData;
  }

protected:
  struct ParentChildRelationship {
    ParentChildRelationship(stk::mesh::EntityId parentElementArg, stk::mesh::EntityIdVector childElementsArg)
      : parentElement(parentElementArg), childElements(childElementsArg) {}

    stk::mesh::EntityId parentElement;
    stk::mesh::EntityIdVector childElements;
  };

  struct NodeSharingData {
    NodeSharingData(stk::mesh::EntityId nodeIdArg, std::vector<int> sharingProcsArg)
      : nodeId(nodeIdArg), sharingProcs(sharingProcsArg) {}

    stk::mesh::EntityId nodeId;
    std::vector<int> sharingProcs;
  };


  void create_parts()
  {
    m_parentPart = & get_meta().declare_part_with_topology("Parent Elements", stk::topology::TRI_3_2D, true);
  }

  void register_fields()
  {
    m_elementWeightField = & get_meta().declare_field<double>(stk::topology::ELEM_RANK, "Element Weights");
    stk::mesh::put_field_on_mesh(*m_elementWeightField, get_meta().universal_part(), nullptr);

    stk::mesh::FieldBase & coordinateField = get_meta().declare_field<double>(stk::topology::NODE_RANK, "model_coordinates");
    stk::mesh::put_field_on_mesh(coordinateField, get_meta().universal_part(), get_meta().spatial_dimension(), nullptr);
  }

  void allocate_bulk(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    std::vector< std::string > names;

    names.reserve( 5 );
    names.push_back(std::string("NODE"));
    names.push_back(std::string("EDGE"));
    names.push_back(std::string("FACE"));
    names.push_back(std::string("ELEMENT"));
    names.push_back(std::string("CONSTRAINTS"));

    stk::mesh::MeshBuilder builder(m_communicator);
    builder.set_spatial_dimension(2);
    builder.set_entity_rank_names(names);
    builder.set_aura_option(auraOption);
    m_bulkData = builder.create();
    m_metaData = &(m_bulkData->mesh_meta_data());
  }

  void create_coarse_mesh()
  {
    std::vector<stk::mesh::EntityIdVector> elemIDsForEachProc {
      {1, 2, 3, 4},
      {5, 6, 7, 8}
    };

    std::map<stk::mesh::EntityId, stk::mesh::EntityIdVector> nodeIDsForEachElem {
      {1, {1,  2, 4}},
      {2, {2,  3, 4}},
      {3, {2,  5, 3}},
      {4, {5,  6, 3}},
      {5, {5,  7, 6}},
      {6, {7,  8, 6}},
      {7, {7,  9, 8}},
      {8, {9, 10, 8}}
    };

    std::vector<NodeSharingData> nodeSharing {
      NodeSharingData(5, {0, 1}),
          NodeSharingData(6, {0, 1})
    };

    build_elements_and_set_up_sharing(elemIDsForEachProc, nodeIDsForEachElem, nodeSharing, *m_parentPart);
  }

  std::vector<ParentChildRelationship> create_child_mesh()
  {
    std::vector<ParentChildRelationship> parentChildRels {
      ParentChildRelationship(3, { 9, 10, 11}),
          ParentChildRelationship(4, {12, 13, 14})
    };

    std::vector<stk::mesh::EntityIdVector> elemIDsForEachProc {
      {9, 10, 11, 12, 13, 14},
      {                     }
    };

    std::map<stk::mesh::EntityId, stk::mesh::EntityIdVector> nodeIDsForEachElem {
      { 9, { 2, 12,  3}},
      {10, { 2, 11, 12}},
      {11, {11,  5, 12}},
      {12, { 5,  6, 12}},
      {13, {12,  6, 13}},
      {14, {12, 13,  3}}
    };

    std::vector<NodeSharingData> nodeSharing;  // No sharing

    stk::mesh::Part & triRootPart = get_meta().get_topology_root_part(stk::topology::TRI_3_2D);

    build_elements_and_set_up_sharing(elemIDsForEachProc, nodeIDsForEachElem, nodeSharing, triRootPart);

    return parentChildRels;
  }

  void build_elements_and_set_up_sharing(const std::vector<stk::mesh::EntityIdVector> & elemIDsForEachProc,
                                         const std::map<stk::mesh::EntityId, stk::mesh::EntityIdVector> & nodeIDsForEachElem,
                                         const std::vector<NodeSharingData> & nodeSharing,
                                         stk::mesh::Part & partForNewElements)
  {
    get_bulk().modification_begin();
    build_elements(elemIDsForEachProc, nodeIDsForEachElem, partForNewElements);
    set_up_node_sharing(nodeSharing);
    get_bulk().modification_end();
  }

  void set_element_weights(ParentChildManager &parentChildManager)
  {
    stk::mesh::EntityVector parentElements;
    stk::mesh::Selector parentElementSelector = (*m_parentPart) & get_meta().locally_owned_part();
    const bool sortById = true;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, parentElementSelector, parentElements, sortById);

    set_parent_element_weights(parentElements, parentChildManager);
  }

  void set_parent_element_weights(const stk::mesh::EntityVector& parentElements, ParentChildManager& parentChildManager)
  {
    for (stk::mesh::Entity parentElement : parentElements) {
      double * weight = stk::mesh::field_data(*m_elementWeightField, parentElement);
      EXPECT_TRUE(nullptr != weight);

      *weight = compute_weight_factor(parentElement, parentChildManager);
    }
  }

  double compute_weight_factor(stk::mesh::Entity parentElement, ParentChildManager& parentChildManager)
  {
    double weightFactor = 1.0;
    if (parentChildManager.has_children(parentElement)) {
      weightFactor = parentChildManager.get_child_elements(parentElement).size();
    }
    return weightFactor;
  }

  void rebalance_all_elements()
  {
    StkFieldWeightRebalance graphSettings(get_bulk(), *m_elementWeightField);
    stk::balance::balanceStkMesh(graphSettings, get_bulk());
  }

  void rebalance_parent_elements_with_manager(ParentChildManager &parentChildManager)
  {
    stk::mesh::Selector selector = (*m_parentPart);
    StkParentRebalance graphSettings(parentChildManager, selector, get_bulk(), *m_elementWeightField);
    stk::balance::balanceStkMesh(graphSettings, get_bulk(), {selector});
  }

  double get_total_weight_for_these_elements(const stk::mesh::EntityVector & solidElements)
  {
    double totalWeightTheseElements = 0.0;
    for (const stk::mesh::Entity element : solidElements) {
      double* elementWeight = stk::mesh::field_data(*m_elementWeightField, element);
      totalWeightTheseElements += (*elementWeight);
    }
    return totalWeightTheseElements;
  }

  double get_total_element_weight_for_this_proc()
  {
    stk::mesh::EntityVector parentElements;
    stk::mesh::Selector parentSelector = (*m_parentPart) & get_meta().locally_owned_part();
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, parentSelector, parentElements);
    return get_total_weight_for_these_elements(parentElements);
  }

  bool is_valid_and_owned(stk::mesh::Entity entity)
  {
    return (get_bulk().is_valid(entity) && get_bulk().bucket(entity).owned());
  }

  struct FamilyRelationships {
    FamilyRelationships(stk::mesh::Entity inParent, stk::mesh::EntityVector inChildren)
      : parent(inParent),
        children(inChildren) {};

    stk::mesh::Entity parent;
    stk::mesh::EntityVector children;
  };

  void check_all_children_are_with_parent(const FamilyRelationships & family)
  {
    if (is_valid_and_owned(family.parent)) {
      for (const stk::mesh::Entity & child : family.children) {
        EXPECT_TRUE(is_valid_and_owned(child));
      }
    }
  }

  void check_expected_parent_with_id_3()
  {
    FamilyRelationships family( get_bulk().get_entity(stk::topology::ELEM_RANK,  3),
    { get_bulk().get_entity(stk::topology::ELEM_RANK,  9),
      get_bulk().get_entity(stk::topology::ELEM_RANK, 10),
      get_bulk().get_entity(stk::topology::ELEM_RANK, 11)
                                } );

    check_all_children_are_with_parent(family);
  }

  void check_expected_parent_with_id_4()
  {
    FamilyRelationships family( get_bulk().get_entity(stk::topology::ELEM_RANK,  4),
    { get_bulk().get_entity(stk::topology::ELEM_RANK, 12),
      get_bulk().get_entity(stk::topology::ELEM_RANK, 13),
      get_bulk().get_entity(stk::topology::ELEM_RANK, 14)
                                } );

    check_all_children_are_with_parent(family);
  }

  bool check_if_all_children_are_with_parent(const FamilyRelationships & family)
  {
    bool childrenAreWithParent = true;
    if (is_valid_and_owned(family.parent)) {
      for (const stk::mesh::Entity & child : family.children) {
        childrenAreWithParent = (childrenAreWithParent && is_valid_and_owned(child));
      }
    }
    return childrenAreWithParent;
  }

  bool check_if_no_children_are_without_parent(const FamilyRelationships & family)
  {
    bool childrenAreWithoutParent = false;
    if (!is_valid_and_owned(family.parent)) {
      for (const stk::mesh::Entity & child : family.children) {
        childrenAreWithoutParent = (childrenAreWithoutParent || is_valid_and_owned(child));
      }
    }
    return (!childrenAreWithoutParent);
  }


  bool check_if_family_is_together(const FamilyRelationships & family)
  {
    bool familyIsTogether = true;
    familyIsTogether = (familyIsTogether && check_if_all_children_are_with_parent(family));
    familyIsTogether = (familyIsTogether && check_if_no_children_are_without_parent(family));
    return familyIsTogether;
  }

  bool get_failure_confirmation_for_parent_with_id_3()
  {
    FamilyRelationships family( get_bulk().get_entity(stk::topology::ELEM_RANK,  3),
    { get_bulk().get_entity(stk::topology::ELEM_RANK,  9),
      get_bulk().get_entity(stk::topology::ELEM_RANK, 10),
      get_bulk().get_entity(stk::topology::ELEM_RANK, 11)
                                } );

    return check_if_family_is_together(family);
  }

  bool get_failure_confirmation_for_parent_with_id_4()
  {
    FamilyRelationships family( get_bulk().get_entity(stk::topology::ELEM_RANK,  4),
    { get_bulk().get_entity(stk::topology::ELEM_RANK, 12),
      get_bulk().get_entity(stk::topology::ELEM_RANK, 13),
      get_bulk().get_entity(stk::topology::ELEM_RANK, 14)
                                } );

    return check_if_family_is_together(family);
  }

  void check_failed_element_partitioning()
  {
    bool failureForParent3 = get_failure_confirmation_for_parent_with_id_3();
    bool failureForParent4 = get_failure_confirmation_for_parent_with_id_4();

    EXPECT_TRUE(failureForParent3 || failureForParent4);
  }

  void check_expected_element_partitioning()
  {
    double totalElementWeight = get_total_element_weight_for_this_proc();
    if (get_bulk().parallel_rank() == 0) {
      EXPECT_EQ(5.0, totalElementWeight);
    }
    if (get_bulk().parallel_rank() == 1) {
      EXPECT_EQ(7.0, totalElementWeight);
    }
    check_expected_parent_with_id_3();
    check_expected_parent_with_id_4();
  }

  void check_number_of_child_elements(ParentChildManager & parentChildManager, stk::mesh::Entity parentElement, size_t numChildrenExpected)
  {
    stk::mesh::EntityVector childElements = parentChildManager.get_child_elements(parentElement);
    EXPECT_TRUE(numChildrenExpected == childElements.size());
  }

  void check_expected_parent_of_children(ParentChildManager & parentChildManager, stk::mesh::Entity parentElement)
  {
    stk::mesh::EntityVector childElements = parentChildManager.get_child_elements(parentElement);
    for (stk::mesh::Entity childElement : childElements) {
      stk::mesh::Entity myParent = parentChildManager.get_parent_element(childElement);
      EXPECT_EQ(myParent, parentElement);
    }
  }

  void verify_parent_child_connectivity(const stk::mesh::EntityVector& parentElements, ParentChildManager& parentChildManager) {
    for (stk::mesh::Entity parentElement : parentElements) {
      if (parentChildManager.has_children(parentElement)) {
        check_number_of_child_elements(parentChildManager, parentElement, 3);
        check_expected_parent_of_children(parentChildManager, parentElement);
      }
    }
  }
  void check_initial_parent_child_relation(ParentChildManager &parentChildManager)
  {
    stk::mesh::EntityVector parentElements;
    stk::mesh::Selector parentSelector = (*m_parentPart) & get_meta().locally_owned_part();
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, parentSelector, parentElements);
    verify_parent_child_connectivity(parentElements, parentChildManager);
  }

  void run_rebalance_all_elements_with_map_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_mesh(auraOption);
    create_child_elements(m_mapManager);
    set_element_weights(m_mapManager);
    check_initial_parent_child_relation(m_mapManager);
    rebalance_all_elements();
    check_failed_element_partitioning();
  }

  void run_rebalance_parent_elements_with_map_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_mesh(auraOption);
    create_child_elements(m_mapManager);
    set_element_weights(m_mapManager);
    check_initial_parent_child_relation(m_mapManager);
    rebalance_parent_elements_with_manager(m_mapManager);
    check_expected_element_partitioning();
  }

  void run_rebalance_parent_elements_with_relation_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_mesh(auraOption);
    create_child_elements(*m_relationManager);
    set_element_weights(*m_relationManager);
    check_initial_parent_child_relation(*m_relationManager);
    rebalance_parent_elements_with_manager(*m_relationManager);
    check_expected_element_partitioning();
  }

protected:
  MPI_Comm m_communicator;
  stk::mesh::MetaData *m_metaData;
  std::shared_ptr<stk::mesh::BulkData> m_bulkData;

  MapParentChildManager      m_mapManager;
  RelationParentChildManager *m_relationManager;

  stk::balance::DoubleFieldType * m_elementWeightField;
  stk::mesh::Part * m_parentPart;

private:
  void build_elements(const std::vector<stk::mesh::EntityIdVector> & elemIDsForEachProc,
                      const std::map<stk::mesh::EntityId, stk::mesh::EntityIdVector> & nodeIDsForEachElem,
                      stk::mesh::Part & partForNewElements)
  {
    const int p_rank = get_bulk().parallel_rank();
    stk::mesh::EntityIdVector elemIDsForThisProc = elemIDsForEachProc[p_rank];
    for (stk::mesh::EntityId elemId : elemIDsForThisProc) {
      stk::mesh::declare_element(get_bulk(), partForNewElements, elemId, nodeIDsForEachElem.at(elemId));
    }
  }

  void add_all_sharing_for_this_node(int myRank, const NodeSharingData & nodeSharing)
  {
    stk::mesh::Entity node = get_bulk().get_entity(stk::topology::NODE_RANK, nodeSharing.nodeId);
    for (const int proc : nodeSharing.sharingProcs) {
      if (proc != get_bulk().parallel_rank()) {
        get_bulk().add_node_sharing(node, proc);
      }
    }
  }

  void set_up_node_sharing(const std::vector<NodeSharingData> & nodeSharing)
  {
    for (size_t i = 0; i < nodeSharing.size(); ++i) {
      for (const int proc : nodeSharing[i].sharingProcs) {
        if (proc == get_bulk().parallel_rank()) {
          add_all_sharing_for_this_node(proc, nodeSharing[i]);
        }
      }
    }
  }

  void connect_parent_with_children(stk::mesh::Entity parentElement, const stk::mesh::EntityIdVector & childElements, ParentChildManager & parentChildManager) {
    for (const stk::mesh::EntityId childElemId : childElements) {
      stk::mesh::Entity childElement = get_bulk().get_entity(stk::topology::ELEM_RANK, childElemId);
      EXPECT_TRUE(is_valid_and_owned(childElement));
      parentChildManager.add_child_element(parentElement, childElement);
    }
  }

  void connect_parents_and_children(const std::vector<ParentChildRelationship>& parentChildRels, ParentChildManager& parentChildManager) {
    for (const ParentChildRelationship& parentChildRel : parentChildRels) {
      stk::mesh::Entity parentElement = get_bulk().get_entity(stk::topology::ELEM_RANK, parentChildRel.parentElement);
      if (is_valid_and_owned(parentElement)) {
        connect_parent_with_children(parentElement, parentChildRel.childElements, parentChildManager);
      }
    }
  }

};

TEST_F(RebalanceParentChildMesh, KrinoBalanceParentsAndChildren2ProcWithAura)
{
  if (2 == stk::parallel_machine_size(get_comm())) {
    run_rebalance_all_elements_with_map_test(stk::mesh::BulkData::AutomaticAuraOption::AUTO_AURA);
  }
}

TEST_F(RebalanceParentChildMesh, KrinoBalanceParentsAndChildren2ProcWithoutAura)
{
  if (2 == stk::parallel_machine_size(get_comm())) {
    run_rebalance_all_elements_with_map_test(stk::mesh::BulkData::AutomaticAuraOption::NO_AUTO_AURA);
  }
}

TEST_F(RebalanceParentChildMesh, KrinoBalanceParents2ProcWithAura)
{
  if (2 == stk::parallel_machine_size(get_comm())) {
    run_rebalance_parent_elements_with_map_test(stk::mesh::BulkData::AutomaticAuraOption::AUTO_AURA);
  }
}

TEST_F(RebalanceParentChildMesh, KrinoBalanceParents2ProcWithoutAura)
{
  if (2 == stk::parallel_machine_size(get_comm())) {
    run_rebalance_parent_elements_with_map_test(stk::mesh::BulkData::AutomaticAuraOption::NO_AUTO_AURA);
  }
}

TEST_F(RebalanceParentChildMesh, EncoreBalanceParents2ProcWithAura)
{
  if (2 == stk::parallel_machine_size(get_comm())) {
    run_rebalance_parent_elements_with_relation_test(stk::mesh::BulkData::AutomaticAuraOption::AUTO_AURA);
  }
}

TEST_F(RebalanceParentChildMesh, EncoreBalanceParents2ProcWithoutAura)
{
  if (2 == stk::parallel_machine_size(get_comm())) {
    run_rebalance_parent_elements_with_relation_test(stk::mesh::BulkData::AutomaticAuraOption::NO_AUTO_AURA);
  }
}
}

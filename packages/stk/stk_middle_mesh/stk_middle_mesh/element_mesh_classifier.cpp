#include "element_mesh_classifier.hpp"

#include "mesh_io.hpp"
#include <map>
#include <set>

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

using predicates::impl::PointClassification;

void ElementMeshClassifier::classify(ElementMeshData& elementMeshData, int numConstraintEdges, const std::vector<mesh::MeshEntityPtr>& mesh2Els)
{
  m_elementMeshData = elementMeshData;
  m_mesh2Els = mesh2Els;
  get_mesh2_element_classification(numConstraintEdges);
}

void ElementMeshClassifier::get_mesh2_element_classification(int numConstraintEdges)
{
  if (M_OUTPUT)
  {
    int myrank = utils::impl::comm_rank(MPI_COMM_WORLD); // CHECK: ALLOW MPI_COMM_WORLD
    mesh::impl::print_vert_edges(std::string("mesh1_el_") + std::to_string(myrank) , {m_elementMeshData.el1});
    mesh::impl::print_vert_edges(std::string("mesh2_els_") + std::to_string(myrank), m_mesh2Els);
    mesh::impl::print_vert_edges(std::string("meshIn_els_") + std::to_string(myrank), m_elementMeshData.elementMeshIn);
  }

  // printTIN("mesh_class", class_i.mesh_in);

  // handle special case
  if (m_mesh2Els.size() == 1)
  {
    classify_all(m_mesh2Els[0]);
    return;
  }

  mesh::impl::MeshAgglomerator agg = setup_agglomerations(numConstraintEdges);

  if (!need_numerical_classification())
  {
    std::vector<mesh::MeshEntityPtr> els2List;
    mesh::MeshEntityPtr el22 = classification_pass1(m_mesh2Els, agg, els2List);

    if (el22)
    {
      if (M_OUTPUT)
        std::cout << "only possible el2 = " << el22;
      classify_all(el22);
      return;
    }

    //------------------------------------------------------------------
    // process elements that might overlap el1 but might not
    // Any elements with two vertices not already classified must be
    // part of el2
    classification_pass2(els2List, agg);
  } else
  {
    //-----------------------------------------------------------------
    // if there are any remaining mesh_in elements not classified, fall back
    // to using PointClassifier to determine the status of the third vertex
    classification_pass3(m_mesh2Els, agg);
  }

  // std::cout << "\nfinished classification" << std::endl;

#ifndef NDEBUG
  auto& elementMeshEntitiesToMeshInEntities = *m_elementMeshData.elementMeshEntitiesToMeshInEntities;
  auto& meshInElementsToMesh2Elements       = *(m_relationalData->meshInElementsToMesh2Elements);
  for (auto& elementMeshEl : m_elementMeshData.elementMeshIn->get_elements())
    if (elementMeshEl)
    {
      mesh::MeshEntityPtr meshInEl = elementMeshEntitiesToMeshInEntities(elementMeshEl, 0, 0);
      // std::cout << "element_mesh_el = " << (void*)element_mesh_el << std::endl;
      // std::cout << "mesh_in_el = " << (void*)meshInEl << std::endl;

      if (!meshInElementsToMesh2Elements(meshInEl, 0, 0))
      {
        std::cout << "element mesh el id = " << elementMeshEl->get_id() << std::endl;
        std::cout << "mesh1 el id = " << m_elementMeshData.el1->get_id() << std::endl;
        std::cout << "mesh1 el = " << m_elementMeshData.el1 << std::endl;
        std::cout << "meshIn El = " << meshInEl << std::endl;
        std::cout << "number of element mesh verts = " << m_elementMeshData.elementMeshIn->get_vertices().size()
                  << std::endl;
        // std::cout << "element mesh vert id 96 maps to vert_in id " <<
        // element_mesh_entities_to_mesh_in_entities(m_elementMeshData.element_mesh_in->get_vertices()[96], 0,
        // 0)->get_id() << std::endl; std::cout << "element mesh vert id 99 maps to vert_in id " <<
        // element_mesh_entities_to_mesh_in_entities(m_elementMeshData.element_mesh_in->get_vertices()[99], 0,
        // 0)->get_id() << std::endl;

        std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
        get_downward(meshInEl, 0, verts.data());
        std::cout << "mesh_in element id = " << meshInEl->get_id()
                  << ", element_mesh element id = " << elementMeshEl->get_id() << std::endl;
        std::cout << "mesh_in vert ids = " << verts[0]->get_id() << ", " << verts[1]->get_id() << ", "
                  << verts[2]->get_id() << std::endl;

        get_downward(elementMeshEl, 0, verts.data());
        std::cout << "element mesh vert ids = " << verts[0]->get_id() << ", " << verts[1]->get_id() << ", "
                  << verts[2]->get_id() << std::endl;

        std::cout << "first vert pt = " << verts[0]->get_point_orig(0) << std::endl;
      }
      assert(meshInElementsToMesh2Elements(meshInEl, 0, 0));
    }
#endif
}

mesh::impl::MeshAgglomerator ElementMeshClassifier::setup_agglomerations(int numConstraintEdges)
{
  auto entityRemap = [&](mesh::MeshEntityPtr entity) {
    auto& elementMeshEntitiesToMeshInEntities = *m_elementMeshData.elementMeshEntitiesToMeshInEntities;
    return elementMeshEntitiesToMeshInEntities(entity, 0, 0);
  };
  auto f = [numConstraintEdges](mesh::MeshEntityPtr edge) { return edge->get_id() < numConstraintEdges; };
  mesh::impl::MeshAgglomerator agg(m_elementMeshData.elementMeshIn, f, entityRemap);

  m_agglomerationsClassified.resize(agg.get_num_groups());
  std::fill(m_agglomerationsClassified.begin(), m_agglomerationsClassified.end(), false);
  // checkGroupsUnique(agg);  //TODO: DEBUGGING

  return agg;
}

void ElementMeshClassifier::get_mesh2_edges(std::vector<mesh::MeshEntityPtr>& edges2)
{
  edges2.clear();
  for (auto& el2 : m_mesh2Els)
  {
    for (int j = 0; j < el2->count_down(); ++j)
      edges2.push_back(el2->get_down(j));
  }

  std::sort(edges2.begin(), edges2.end(), mesh::is_less);
  auto it = std::unique(edges2.begin(), edges2.end());
  edges2.erase(it, edges2.end());
}

bool ElementMeshClassifier::is_classified_on_element(const predicates::impl::PointRecord& record)
{
  assert(record.type != PointClassification::Exterior);
  if (m_elementMeshData.el1 == record.el)
    return true;
  else if (record.type == PointClassification::Interior && m_elementMeshData.el1 != record.el)
    return false;
  else
  {
    assert(record.type == PointClassification::Vert || record.type == PointClassification::Edge);

    mesh::MeshEntityPtr entity = get_entity(record);

    int dim = record.type == PointClassification::Vert ? 0 : 1;
    std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> entities;
    int nentities = get_downward(m_elementMeshData.el1, dim, entities.data());

    for (int i = 0; i < nentities; ++i)
      if (entities[i] == entity)
        return true;

    return false;
  }
}

bool ElementMeshClassifier::need_numerical_classification()
{
  // TODO: bail out early in easy case: number of element_mesh_in verts > 6?

  auto& edges2ToFakeVertsIn    = *(m_relationalData->edges2ToFakeVertsIn);
  auto& fakeVertsToRealVertsIn = m_relationalData->fakeVertsToVertsIn;
  auto& vertsInClassOnMesh1    = *(m_relationalData->vertsInClassOnMesh1);
  std::vector<mesh::MeshEntityPtr> edges2;
  get_mesh2_edges(edges2);

  SetType<mesh::MeshEntityPtr> vertsFromEdges2;
  for (auto& edge2 : edges2)
  {
    for (int i = 0; i < edges2ToFakeVertsIn.get_num_comp(edge2, 0); ++i)
    {
      FakeVert fv                                 = edges2ToFakeVertsIn(edge2, 0, i).vert;
      mesh::MeshEntityPtr vertIn                  = fakeVertsToRealVertsIn[fv.id];
      if (vertIn)
      {
        const predicates::impl::PointRecord& record = vertsInClassOnMesh1(vertIn, 0, 0);
        if (is_classified_on_element(record))
        {
          vertsFromEdges2.insert(vertIn);
        }
      }
    }
  }

  return vertsFromEdges2.size() <= 2;
}

mesh::MeshEntityPtr ElementMeshClassifier::classification_pass1(const std::vector<mesh::MeshEntityPtr>& mesh2Els,
                                                                mesh::impl::MeshAgglomerator& agg,
                                                                std::vector<mesh::MeshEntityPtr>& els2List)
{
  if (M_OUTPUT)
    std::cout << "\nEntered classificationPass1" << std::endl;

  auto& meshInElementsToMesh2Elements = *(m_relationalData->meshInElementsToMesh2Elements);

  std::vector<mesh::MeshEntityPtr> elsIn;
  SetType<mesh::MeshEntityPtr> vertsIn;
  // std::vector<mesh::MeshEntityPtr> els2_list;
  els2List.clear();
  int nelemTwoVerts               = 0;       // number of elements with at least two vertices on mesh_in
  mesh::MeshEntityPtr el2TwoVerts = nullptr; // if nelem_2 == 1 at the end of the loop, this
                                             // is the element with only 2 vertices
  for (auto& el2 : mesh2Els)
  {
    if (M_OUTPUT)
    {
      std::cout << "el2 = " << el2->get_id() << ", " << el2 << std::endl;

      std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> el2Verts;
      int ndown = mesh::get_downward(el2, 0, el2Verts.data());
      std::cout << "  el2 vert ids = ";
      for (int i=0; i < ndown; ++i)
        std::cout << el2Verts[i]->get_id() << ", ";
      std::cout << std::endl;
    }
    vertsIn.clear();

    get_verts_in(m_elementMeshData.el1, el2, vertsIn);
    if (M_OUTPUT)
    {
      std::cout << "verts_in:" << std::endl;
      for (auto& v : vertsIn)
      {
        std::cout << v << ", id = " << v->get_id()
                  << ", element mesh vert id = " << get_element_mesh_entity(v)->get_id() << std::endl;
      }
    }

    if (vertsIn.size() <= 1) // el2 is exterior element
      continue;
    else if (vertsIn.size() == 2) // el2 could cut two edges of el1, or
    {                             // could have an edge overlap,
      els2List.push_back(el2);
      nelemTwoVerts++;
      el2TwoVerts = el2;
      continue;
    } else
    {
      nelemTwoVerts++;
      el2TwoVerts = el2;
    }

    if (are_all_verts_on_same_edge(vertsIn, m_relationalData->vertsInClassOnMesh1))
      continue;

    // else: 3 or more vertices: elements can be definitively classified
    //       one el2
    elsIn.clear();
    // get_commonElements(verts_in, 3, els_in);
    auto idxs = agg.get_group_idxs(vertsIn, 3);
    if (M_OUTPUT)
    {
      std::cout << "idxs = ";
      for (auto idx : idxs)
        std::cout << idx << ", ";
      std::cout << std::endl;
    }
    assert(idxs.size() == 1);
    
    // std::cout << "classifying group " << idxs[0] << std::endl;

    for (auto& elIn : agg.get_group_elements(idxs[0]))
    {
      if (M_OUTPUT && elIn->get_id() == 2660)
        std::cout << "classifying el_in id " << elIn->get_id() << std::endl;

      if (M_OUTPUT)
      {
        std::cout << "classifying el_in id " << elIn->get_id() << std::endl;
        // std::cout << "classifying element with centroid " << compute_centroid_3d(el_in) << " on el2" << std::endl;
        // mesh::MeshEntityPtr _verts_in[mesh::MAX_DOWN]; get_downward(el_in, 0, _verts_in);
        // std::cout << "vert_ids = " << _verts_in[0]->get_id() << ", " << _verts_in[1]->get_id() << ", " <<
        // _verts_in[2]->get_id() << std::endl;
      }
      assert(!(meshInElementsToMesh2Elements(elIn, 0, 0)));
      meshInElementsToMesh2Elements(elIn, 0, 0) = el2;
      m_agglomerationsClassified[idxs[0]]       = true;
    }
  }

  if (nelemTwoVerts == 1)
    return el2TwoVerts;
  else
    return nullptr;
}

// TODO: DEBUGGING
mesh::MeshEntityPtr ElementMeshClassifier::get_element_mesh_entity(mesh::MeshEntityPtr meshInEntity)
{
  auto& elementMeshEntitesToMeshInEntities = *m_elementMeshData.elementMeshEntitiesToMeshInEntities;
  auto& elementMeshEntities =
      m_elementMeshData.elementMeshIn->get_mesh_entities(get_type_dimension(meshInEntity->get_type()));

  for (auto& elementMeshEntity : elementMeshEntities)
    if (elementMeshEntity && elementMeshEntitesToMeshInEntities(elementMeshEntity, 0, 0) == meshInEntity)
      return elementMeshEntity;

  throw std::runtime_error("unable to find entity");
}

void ElementMeshClassifier::classification_pass2(const std::vector<mesh::MeshEntityPtr>& els2List,
                                                 mesh::impl::MeshAgglomerator& agg)
{
  // the logic of this function uses process of elimination to determine
  // that certain mesh_in elements are on a given mesh2 element.  The
  // process of elimination only works if we excluded certain cases
  // where there is no way to figure out which mesh2 element
  // a mesh_in element is inside of.
  assert(!need_numerical_classification());
  std::vector<mesh::MeshEntityPtr> elsIn;
  SetType<mesh::MeshEntityPtr> vertsIn;
  auto& meshInElementsToMesh2Elements = *(m_relationalData->meshInElementsToMesh2Elements);

  if (M_OUTPUT)
    std::cout << "\nEntered classificationPass2" << std::endl;
  for (auto& el2 : els2List)
  {
    elsIn.clear();
    vertsIn.clear();

    if (M_OUTPUT)
      std::cout << "el2 = " << el2->get_id() << std::endl;

    get_verts_in(m_elementMeshData.el1, el2, vertsIn);
    if (M_OUTPUT)
    {
      std::cout << "verts_in:" << std::endl;
      for (auto& v : vertsIn)
        std::cout << v << ", id = " << v->get_id() << std::endl;
    }

    if (vertsIn.size() <= 1 || vertsIn.size() >= 4)
      continue;

    auto idxs = agg.get_group_idxs(vertsIn, 2);
    idxs      = remove_already_classified_groups(idxs);
    if (idxs.size() != 1)
    {
      continue;
    }

    // std::cout << "classifying group " << idxs[0] << " on el2" << std::endl;

    for (auto& elIn : agg.get_group_elements(idxs[0]))
    {
      if (M_OUTPUT)
      {
        std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
        get_downward(elIn, 0, verts.data());
        std::cout << "el_in vert ids = " << verts[0]->get_id() << ", " << verts[1]->get_id() << ", "
                  << verts[2]->get_id() << std::endl;
      }
      bool isElInUnclassified = !meshInElementsToMesh2Elements(elIn, 0, 0); // TODO: should be an assert?

      if (isElInUnclassified)
      {
        if (M_OUTPUT)
        {
          std::cout << "classifying element with centroid " << compute_centroid_3d(elIn) << " on el2" << std::endl;
          mesh::MeshEntityPtr vertsInOutput[mesh::MAX_DOWN];
          get_downward(elIn, 0, vertsInOutput);
          std::cout << "vert_ids = " << vertsInOutput[0]->get_id() << ", " << vertsInOutput[1]->get_id() << ", "
                    << vertsInOutput[2]->get_id() << std::endl;
        }
        // Note: normally we would classify *all* elements in the same mesh agglomeration
        //       group on el2, but pass 1 handles all multi-element groups.  All that is
        //       left to do in this pass is cut off corners, which only have 1 element
        //       in the group.
        meshInElementsToMesh2Elements(elIn, 0, 0) = el2;
        m_agglomerationsClassified[idxs[0]]       = true;
      }
    }
  }
}

void ElementMeshClassifier::classification_pass3(const std::vector<mesh::MeshEntityPtr>& els2List,
                                                 mesh::impl::MeshAgglomerator& agg)
{
  if (M_OUTPUT)
    std::cout << "\nEntered classificationPass3" << std::endl;

  if (agg.get_num_groups() > 2)
    throw std::runtime_error("too many groups for pass 3");

  std::vector<mesh::MeshEntityPtr> elsIn;
  SetType<mesh::MeshEntityPtr> vertsIn;
  auto& meshInElementsToMesh2Elements = *(m_relationalData->meshInElementsToMesh2Elements);
  for (auto& el2 : els2List)
  {
    elsIn.clear();
    vertsIn.clear();

    if (M_OUTPUT)
      std::cout << "el2 = " << el2 << std::endl;

    get_verts_in(m_elementMeshData.el1, el2, vertsIn);

    if (M_OUTPUT)
    {
      std::cout << "verts_in:" << std::endl;
      for (auto& v : vertsIn)
        std::cout << v << ", id = " << v->get_id() << std::endl;
    }

    if (vertsIn.size() <= 1 || vertsIn.size() >= 4)
      continue;

    get_common_elements(vertsIn, 2, elsIn);

    for (auto& elIn : elsIn)
    {
      mesh::MeshEntityPtr elInVerts[mesh::MAX_DOWN];
      get_downward(elIn, 0, elInVerts);
      if (M_OUTPUT)
      {
        std::cout << "el_in has vert ids = " << elInVerts[0]->get_id() << ", " << elInVerts[1]->get_id() << ", "
                  << elInVerts[2]->get_id() << std::endl;
      }
      if (!meshInElementsToMesh2Elements(elIn, 0, 0) && classify_third_vert(el2, vertsIn, elIn))
      {
        // classify all elements in same agglomeration group on mesh2 element
        for (int i = 0; i < 3; ++i)
          vertsIn.insert(elInVerts[i]);

        auto idxs = agg.get_group_idxs(vertsIn, vertsIn.size());
        assert(idxs.size() == 1);
        classify_group(agg, idxs[0], el2);
      }
    }
  }

  // in the unlikely event the above loop is unable to classify an element due to
  // numerical tolerances, pick the best possible element
  for (auto& elementMeshElIn : m_elementMeshData.elementMeshIn->get_elements())
    if (elementMeshElIn)
    {
      mesh::MeshEntityPtr elIn = (*m_elementMeshData.elementMeshEntitiesToMeshInEntities)(elementMeshElIn, 0, 0);
      if (elIn && !(meshInElementsToMesh2Elements(elIn, 0, 0)))
      {
        mesh::MeshEntityPtr el2                   = get_most_inside_element(elIn, els2List);
        meshInElementsToMesh2Elements(elIn, 0, 0) = el2;
      }
    }
}

void ElementMeshClassifier::classify_group(mesh::impl::MeshAgglomerator& agg, int group, mesh::MeshEntityPtr el2)
{
  auto& meshInElementsToMesh2Elements = *(m_relationalData->meshInElementsToMesh2Elements);

  m_agglomerationsClassified[group] = true;
  for (auto& elIn : agg.get_group_elements(group))
    if (!meshInElementsToMesh2Elements(elIn, 0, 0))
    {
      if (M_OUTPUT)
      {
        std::cout << "classifying element with centroid " << compute_centroid_3d(elIn) << " on el2" << std::endl;
        mesh::MeshEntityPtr vertsIn[mesh::MAX_DOWN];
        get_downward(elIn, 0, vertsIn);
        std::cout << "vert_ids = " << vertsIn[0]->get_id() << ", " << vertsIn[1]->get_id() << ", "
                  << vertsIn[2]->get_id() << std::endl;
      }
      meshInElementsToMesh2Elements(elIn, 0, 0) = el2;
    }
}

// returns true if the third vert is on el2
bool ElementMeshClassifier::classify_third_vert(mesh::MeshEntityPtr el2, SetType<mesh::MeshEntityPtr>& vertsIn,
                                                mesh::MeshEntityPtr elIn)
{
  mesh::MeshEntityPtr elInVerts[mesh::MAX_DOWN];
  int nverts = get_downward(elIn, 0, elInVerts);

  mesh::MeshEntityPtr thirdVert = nullptr;
  for (int i = 0; i < nverts; ++i)
    if (vertsIn.count(elInVerts[i]) == 0)
    {
      thirdVert = elInVerts[i];
      break;
    }

  if (!thirdVert)
    return false;

  auto r = m_classifier->classify_reverse(el2, thirdVert->get_point_orig(0));
  return r.type != PointClassification::Exterior;
}

// get all verts on mesh_in that are associated with edges
// of the given element on mesh2
void ElementMeshClassifier::get_verts_in(mesh::MeshEntityPtr /*el1*/, mesh::MeshEntityPtr el2,
                                         SetType<mesh::MeshEntityPtr>& vertsIn)
{
  auto& edges2ToFakeVertsIn    = *(m_relationalData->edges2ToFakeVertsIn);
  auto& fakeVertsToRealVertsIn = m_relationalData->fakeVertsToVertsIn;
  auto& vertsInClassOnMesh1    = *(m_relationalData->vertsInClassOnMesh1);

  for (int i = 0; i < el2->count_down(); ++i)
  {
    mesh::MeshEntityPtr edgeI = el2->get_down(i);
    for (int j = 0; j < edges2ToFakeVertsIn.get_num_comp(edgeI, 0); ++j)
    {
      FakeVert fv                                 = edges2ToFakeVertsIn(edgeI, 0, j).vert;
      mesh::MeshEntityPtr vertIn                  = fakeVertsToRealVertsIn[fv.id];
      if (vertIn)
      {
        const predicates::impl::PointRecord& record = vertsInClassOnMesh1(vertIn, 0, 0);

        if (is_classified_on_element(record))
          vertsIn.insert(vertIn);
      }
    }
  }
}

std::vector<int> ElementMeshClassifier::remove_already_classified_groups(const std::vector<int>& groups)
{
  std::vector<int> filteredGroups;
  for (auto& group : groups)
    if (!m_agglomerationsClassified[group])
      filteredGroups.push_back(group);

  return filteredGroups;
}

// gets elements that have nthres or more vertices from verts_in
void ElementMeshClassifier::get_common_elements(const SetType<mesh::MeshEntityPtr>& vertsIn, const unsigned int nthres,
                                                std::vector<mesh::MeshEntityPtr>& commonEls)
{
  MapType<mesh::MeshEntityPtr, unsigned int> elCounts;
  std::vector<mesh::MeshEntityPtr> upEls;

  for (auto& v : vertsIn)
  {
    get_upward(v, 2, upEls);
    for (auto& el : upEls)
      elCounts[el] += 1;
  }

  commonEls.clear();
  for (auto& p : elCounts)
    if (p.second >= nthres)
      commonEls.push_back(p.first);
}

mesh::MeshEntityPtr ElementMeshClassifier::get_most_inside_element(mesh::MeshEntityPtr elIn,
                                                                   const std::vector<mesh::MeshEntityPtr>& els2)
{
  auto& elementMeshEntitesToMeshInEntities = *m_elementMeshData.elementMeshEntitiesToMeshInEntities;
  auto& meshInElementsToMesh2Elements      = *(m_relationalData->meshInElementsToMesh2Elements);

  double minDist             = std::numeric_limits<double>::max();
  mesh::MeshEntityPtr minEl2 = nullptr;
  utils::Point centroid      = compute_centroid_3d(elIn);
  for (auto& el2 : els2)
  {
    // TODO: this is inefficient, but this function only gets called when
    //       the numerical classification pass fails (ie. there are only 2
    //       agglomeration groups and one of them is a sliver element)
    bool el2AlreadyClassified = false;
    for (auto& elementMeshEl : m_elementMeshData.elementMeshIn->get_elements())
      if (elementMeshEl)
      {
        mesh::MeshEntityPtr meshInElement = elementMeshEntitesToMeshInEntities(elementMeshEl, 0, 0);
        if (meshInElementsToMesh2Elements(meshInElement, 0, 0) == el2)
        {
          el2AlreadyClassified = true;
          break;
        }
      }

    if (!el2AlreadyClassified)
    {
      auto r = m_classifier->classify_reverse(el2, centroid);
      for (int i = 0; i < el2->count_down(); ++i)
      {
        double dist = m_classifier->compute_orthogonal_dist(r, i);
        if (dist < minDist)
        {
          minDist = dist;
          minEl2  = el2;
        }
      }
    }
  }

  return minEl2;
}

// classifies all remaining mesh_in elements on el2
void ElementMeshClassifier::classify_all(mesh::MeshEntityPtr el2)
{
  auto& elementMeshEntitiesToMeshInEntities = *m_elementMeshData.elementMeshEntitiesToMeshInEntities;
  auto& meshInElementsToMesh2Elements       = *(m_relationalData->meshInElementsToMesh2Elements);
  assert(m_relationalData->meshInElementsToMesh2Elements);

  for (auto& elementMeshIn : m_elementMeshData.elementMeshIn->get_elements())
    if (elementMeshIn)
    {
      mesh::MeshEntityPtr meshInElement = elementMeshEntitiesToMeshInEntities(elementMeshIn, 0, 0);
      assert(meshInElement);
      if (!(meshInElementsToMesh2Elements(meshInElement, 0, 0)))
        meshInElementsToMesh2Elements(meshInElement, 0, 0) = el2;
    }
}

// returns true if all verts in the set are on the same mesh1 edge (or
// the vertices that bound it)
bool are_all_verts_on_same_edge(std::set<mesh::MeshEntityPtr, mesh::MeshEntityCompare>& vertsIn,
                                mesh::FieldPtr<predicates::impl::PointRecord> vertsInClassOnMesh1Ptr)
{
  auto& vertsInClassOnMesh1 = *vertsInClassOnMesh1Ptr;

  if (vertsIn.size() == 0)
    return false;
  else if (vertsIn.size() == 1)
  {
    return vertsInClassOnMesh1(*(vertsIn.begin()), 0, 0).type != PointClassification::Interior;
  }

  // do the fast test first: if any point on interior
  for (auto& vertIn : vertsIn)
  {
    const predicates::impl::PointRecord& record = vertsInClassOnMesh1(vertIn, 0, 0);
    assert(record.type != PointClassification::Exterior);
    if (record.type == PointClassification::Interior)
      return false;
  }

  if (vertsIn.size() == 2)
  {
    const predicates::impl::PointRecord& record1 = vertsInClassOnMesh1(*(vertsIn.begin()), 0, 0);
    const predicates::impl::PointRecord& record2 = vertsInClassOnMesh1(*(++vertsIn.begin()), 0, 0);
    if (record1.type == PointClassification::Edge || record2.type == PointClassification::Edge)
    {
      if (record1.type == PointClassification::Edge && record2.type == PointClassification::Edge)
        return get_entity(record1) == get_entity(record2);

      const predicates::impl::PointRecord& recordEdge = record1.type == PointClassification::Edge ? record1 : record2;
      const predicates::impl::PointRecord& recordVert = record1.type == PointClassification::Edge ? record2 : record1;

      mesh::MeshEntityPtr edge = get_entity(recordEdge);
      mesh::MeshEntityPtr vert = get_entity(recordVert);

      return edge->get_down(0) == vert || edge->get_down(1) == vert;
    } else // both are verts
    {
      mesh::MeshEntityPtr vert1 = get_entity(record1);
      mesh::MeshEntityPtr vert2 = get_entity(record2);
      return get_common_edge(vert1, vert2) != nullptr;
    }
  }

  // 3 or more points
  mesh::MeshEntityPtr vertOnEdge = nullptr;
  for (auto& vertIn : vertsIn)
    if (vertsInClassOnMesh1(vertIn, 0, 0).type == PointClassification::Edge)
    {
      vertOnEdge = vertIn;
      break;
    }

  if (!vertOnEdge) // if 3 points and none are on edge (all on vertices),
                   // they can't be on the same edge
    return false;

  mesh::MeshEntityPtr edge = get_entity(vertsInClassOnMesh1(vertOnEdge, 0, 0));
  for (auto& vertIn : vertsIn)
    if (vertIn != vertOnEdge)
    {
      const predicates::impl::PointRecord& record = vertsInClassOnMesh1(vertIn, 0, 0);
      mesh::MeshEntityPtr entity                  = get_entity(record);
      if (record.type == PointClassification::Edge && edge != entity)
        return false;
      else if (record.type == PointClassification::Vert && (edge->get_down(0) != entity && edge->get_down(1) != entity))
        return false;
    }

  return true;
}

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk

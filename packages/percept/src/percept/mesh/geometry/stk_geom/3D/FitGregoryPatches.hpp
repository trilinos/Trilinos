// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef FitGregoryPatches_hpp
#define FitGregoryPatches_hpp

#include <map>
#include <vector>

#include <percept/PerceptMesh.hpp>
#if HAVE_YAML
#include <percept/YamlUtils.hpp>
#endif

/** Computes Gregory patch control points for G1 continuity for a set of
 * surface sets - each surface set consists of surfaces to be considered
 * together so seams between the surfaces are ignored, modulo the angle
 * criterion.  If @param surfaceSets is empty, all surfaces are processed
 * together.
 *
 * Angle criterion: can be specified for each surface set or individual
 * surfaces, or for all surfaces. @param angleMap specifies a map of surface
 * name to angle criterion - if empty, use the specified scalar @param
 * globalAngleCriterion.  If any surface has an empty angleMap entry, use
 * globalAngleCriterion.  If any angle criterion is < 0, it is ignored and
 * averaged normals are used.
 *
 * Note: angles are specified in degrees.
 */
namespace percept {

  struct FGP_QA {
    FGP_QA() : m_activate(false), m_num_divisions(2) {}
    bool m_activate;
    std::string m_file;
    std::vector<std::string> m_filenames;
    int m_num_divisions;
  };


class FitGregoryPatches
{
  static const double m_globalAngleCriterionDefault;
public:
  typedef std::map<std::string,  std::vector<std::string> > SurfaceSets;
  typedef std::map<std::string, double> AngleMap;

  FitGregoryPatches(PerceptMesh& eMesh, SurfaceSets surfaceSets = SurfaceSets(), AngleMap angleMap = AngleMap(), double globalAngleCriterion = m_globalAngleCriterionDefault)
    : m_reverseAll(true), m_eMesh(eMesh), m_surfaceSets(surfaceSets), m_angleMap(angleMap), m_globalAngleCriterion(globalAngleCriterion),
      m_debug(false), m_edgeSeamsPart(0)
  {
  }

  // call this before committing the MetaData to register required fields
  void register_or_set_fields(bool doRegister=true);

  void computeControlPoints(bool doGetNormals = true, bool createEdgeSeamsPart = false);

  /// creates a surface mesh from the surfaces of the PerceptMesh, and optionally fits and refines the mesh
  ///   for display for Q/A purposes
  void createQA(const std::string& QA_file);

  // read surface sets and angle map from YAML file
  void parse(const std::string& file_name);

  static int MaxControlPoints() { return 20; }
  static int NumTriControlPoints() { return 18; }
  static int NumQuadControlPoints() { return 20; }

public:

  // utilities
  static
  std::string
  printForMathematica(PerceptMesh& eMesh, stk::mesh::Entity face, bool convert=true, bool printHeader=true, unsigned lineLength=132);

  // utilities
  std::string
  printForMathematica(stk::mesh::Entity face, bool convert=true, bool printHeader=true, unsigned lineLength=132)
  {
    return printForMathematica(m_eMesh, face, convert, printHeader, lineLength);
  }

  // advanced use

  /// if parts is null, all surfaces are processed
  void
  getNormals(stk::mesh::PartVector* parts = 0);

  void
  faceNormal(stk::mesh::Entity face, double normal[3]);

  void
  getCurrentParts(std::vector<stk::mesh::PartVector>& currentParts);

  void
  findSeams(stk::mesh::PartVector& parts, bool createEdgeSeamsPart = false);

  static bool is_surface_topology(stk::topology topo);
  static bool is_tri_topology(stk::topology topo);
  static bool is_quad_topology(stk::topology topo);

  typedef std::map<stk::mesh::Entity, std::vector<int> > EdgeSeamMap;
  typedef std::pair<stk::mesh::Entity, stk::mesh::Entity > Edge;
  typedef std::set<Edge> EdgeSet;
  typedef std::vector<Edge> EdgeVector;
  typedef std::map<stk::mesh::Entity, EdgeVector > NodeToEdgeMap;
  typedef std::list<stk::mesh::Entity> EdgeList;
  typedef std::vector<stk::mesh::Entity> EntityVector;
  typedef std::set<stk::mesh::Entity> EntitySet;

protected:

  bool in_surface_sets(const std::string& partToTest);

#if HAVE_YAML
  void parse(const YAML::Node& node);
  void emit(const YAML::Node& node);
#endif

  void
  fitCubics(stk::mesh::PartVector& parts);

  void
  extractRibbon(stk::mesh::Entity face, int edge, bool reverse, MDArray& p, MDArray& q);

  void
  putRibbon(stk::mesh::Entity face, int edge, bool reverse, MDArray& p_in);

  void
  fitRibbons(stk::mesh::PartVector& parts);

  /**
   * FIXME - consider moving the following methods dealing with edges to another class
   */

  bool
  isSeam(stk::mesh::Entity face_0, stk::mesh::Entity face_1, double angleCriterion);

  Edge
  create_edge(stk::mesh::Entity face, unsigned edge);

  template<typename NodeType>
  std::pair<NodeType, NodeType>
  create_edge(NodeType n0, NodeType n1)
  {
    typedef std::pair<NodeType, NodeType> LEdge;
    if (n0 < n1)
      return LEdge(n0,n1);
    else
      return LEdge(n1,n0);
  }

  template<typename NodeType, typename Comp>
  std::pair<NodeType, NodeType>
  create_edge(NodeType n0, NodeType n1, const Comp& comp) const
  {
    typedef std::pair<NodeType, NodeType> LEdge;
    if (comp(n0,n1))
      return LEdge(n0,n1);
    else
      return LEdge(n1,n0);
  }

  /**
   * the following methods split edges into loops/contiguous subsets
   *   then further sub-sets them based on finding corners with an
   *   angle criterion; also adds edges to the QA mesh for debugging
   *   and Q/A purposes.
   */
  typedef std::array<double,3> Point;

  double
  edgeAngle(stk::mesh::Entity node, const Edge& e0, const Edge& e1);

  void
  processSeams(stk::mesh::PartVector& parts, bool createEdgeSeamsParts);

  void
  findContiguousEdgeSets(std::vector<EdgeSet>& contiguousEdgeSets, NodeToEdgeMap& nodeToEdgeMap,
                         const EdgeSet& mainEdgeSet);

  void
  findGhostEdges(EdgeSet& mainEdgeSet, NodeToEdgeMap& nodeToEdgeMap);

  void
  findEdgeSetCorners(EntitySet& corners, double angleCriterion, const std::vector<EdgeSet>& contiguousEdgeSets, NodeToEdgeMap& nodeToEdgeMap);

  Point
  getTangent(stk::mesh::Entity n0, stk::mesh::Entity n1);

  void
  findTangentVectors(const Edge& edge, const EntitySet& corners, const NodeToEdgeMap& nodeToEdgeMap);

  void
  findTangentVectors(const EntitySet& corners, const std::vector<EdgeSet>& contiguousEdgeSets, const NodeToEdgeMap& nodeToEdgeMap);

  void
  fitCubicWithTangents(MDArray& c, const MDArray& pi, const MDArray& pj, const Point& ti, const Point& tj);

  void
  sortEdges(NodeToEdgeMap& nodeToEdgeMap);

  void
  addEdgesToMesh(EdgeSet& edgeSet);

  void
  addEdgesToQAMesh(std::vector<EdgeSet>& contiguousEdgeSets);

  /** return -1 if @param face doesn't contain @param node, thus the face/node pair
   *    should use the averaged normal from the node (i.e. it's a "regular" node).
   *  return 0, 1, ... Nedge-1 (where Nedge is number of edges attached to the node)
   *  which is the base edge of the "quadrant" spanned by the two edges going around
   *  the node in edge-sorted order.  This allows determining if the face is on one
   *  side or the other (the simple case where the node has two edges attached),
   *  or which quadrant the face is in thus telling which of the normals to use
   *
   *  high-level interface
   *
   *  Notes:
   *
   *  When finding bi-cubic coefficients, and there is an edge seam,
   *  the node contains the averaged normal, but we want to use a
   *  normal from the side that the face is on, not the averaged one,
   *  so we need to know which side of the edge we are on and we
   *  choose one of the two normals computed there (actually, which "quadrant"
   *  the face is in).
   *
   */
  int
  orient(stk::mesh::Entity face, stk::mesh::Entity node, stk::mesh::Selector& sel, bool debug = false);

  /** low-level interface - for top-level, a recursive function
   *  is called, see orient and orientRecurse
   *
   */
  int
  orientLowLevel(stk::mesh::Entity face, stk::mesh::Entity node, bool debug = false);

  int
  orientRecurse(stk::mesh::Entity face, stk::mesh::Entity node, EntitySet& visited, stk::mesh::Selector& sel, bool debug = false);

  void
  getEdgeNeighborsSharingNode(stk::mesh::Entity face, stk::mesh::Selector& sel, stk::mesh::Entity node, EntitySet& neighbors_in);


  std::string printEdge(const Edge& edge);

public:
  // used for unit testing - this parameter reverses all p,q,r arrays going into and
  // out of the fitting routines to be consistent with the FarinHansford paper
  bool m_reverseAll;

private:

  PerceptMesh& m_eMesh;
  SurfaceSets m_surfaceSets;
  AngleMap m_angleMap;
  double m_globalAngleCriterion;
  bool m_debug;
#if HAVE_YAML
  YAML::Node m_node, m_node1;
#endif

public:
  FGP_QA m_QA;

  EdgeSeamMap m_edgeSeamsMap;
  stk::mesh::Part *m_edgeSeamsPart;
  EdgeSet m_edgeSet; // all edges
  EdgeSet m_edgeSetAll; // all edges
  std::vector<EdgeSet> m_contiguousEdgeSets; // m_edgeSet broken into smooth subsets
  std::vector<EdgeSet> m_contiguousEdgeSetsAll; // m_edgeSet broken into smooth subsets
  typedef std::vector<Point> PointVector;
  std::map<Edge, PointVector > m_tangentVectors;

  NodeToEdgeMap m_nodeToEdgeMap;
  EntitySet m_corners;

  typedef std::map<stk::mesh::Entity, PointVector > NodeNormalsMap;
  NodeNormalsMap m_nodeNormalsMap;
  typedef std::map<stk::mesh::Entity, std::vector<int> > NodeNormalsSetMap;
  NodeNormalsSetMap m_nodeNormalsSetMap;
};

#if HAVE_BOOST_GRAPH
  class BGraphExternal {
  public:
    // external interface
    static void
    findContiguousEdgeSets(PerceptMesh& eMesh,
                           std::vector<FitGregoryPatches::EdgeSet>& contiguousEdgeSets,
                           FitGregoryPatches::NodeToEdgeMap& nodeToEdgeMap,
                           const FitGregoryPatches::EdgeSet& mainEdgeSet);
  };
#endif

}

#endif

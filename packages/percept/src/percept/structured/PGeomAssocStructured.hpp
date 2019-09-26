
#include <memory>
#include <vector>
#include <map>

class RefFace;
class RefEdge;
class RefVertex;
class PGeom;
class PGeomRangeSearch;
class CubitVector;

#include <percept/structured/BlockStructuredGrid.hpp>

enum IJK_DIRECTION
{
    I_DIRECTION = 0,
    J_DIRECTION = 1,
    K_DIRECTION = 2
};

struct CloseGeometry
{
    int geometry_id;
    double distance;

    CloseGeometry(int id, double dist)
    : geometry_id(id), distance(dist)
    {}
    virtual ~CloseGeometry() {};
};

struct NodeAssociationData
{
  std::vector<struct CloseGeometry> close_surfs;
  std::vector<struct CloseGeometry> close_curves;
  std::vector<struct CloseGeometry> close_verts;
  double edge_len_tol;

  void add_surface(int id, double distance)
  {
      close_surfs.push_back(CloseGeometry(id, distance));
  }
  void add_curve(int id, double distance)
  {
      close_curves.push_back(CloseGeometry(id, distance));
  }
  void add_vertex(int id, double distance)
  {
      close_verts.push_back(CloseGeometry(id, distance));
  }
};

struct NodeAssociatedGeom
{
    int dimension;
    std::vector<int> geom_ids;
};

struct EdgeAssociationData
{
  std::vector<int> close_surfs;
  std::vector<int> close_curves;
};

struct FaceAssociationData
{
  std::vector<int> close_surfs;
};

struct BlockAssociationData
{
    std::map<int, NodeAssociationData> node_to_geom_map;

    std::multimap< int, int > surface_to_node_map;
    std::multimap< int, int > curve_to_node_map;
    std::multimap< int, int > vertex_to_node_map;

    std::map<int, EdgeAssociationData> edge_to_geom_map[3];

    std::multimap<int, int> surface_to_edge_map[3];
    std::multimap<int, int> curve_to_edge_map[3];

    std::map<int, FaceAssociationData> face_to_geom_map[3];

    std::multimap<int, int> surface_to_face_map[3];
};

class PGeomAssocStructured
{
public:

  void find_mesh_association(std::shared_ptr<PGeom> pgeom,
                             std::shared_ptr<percept::BlockStructuredGrid> grid,
                             const bool local_tol,
                             const double associate_tolerance,
                             const int print_level = 0);

  void make_refined_mesh_association(std::shared_ptr<percept::BlockStructuredGrid> input_grid,
                                     std::shared_ptr<percept::BlockStructuredGrid> output_grid);


  std::map<uint64_t, NodeAssociatedGeom> & get_refined_node_map(unsigned iblock);

public: // for testing purposes
  void check_cylinder_surface_projection(std::shared_ptr<PGeom> pgeom,
                                         std::shared_ptr<percept::BlockStructuredGrid> grid,
                                         int surface_id,
                                         double radius,
                                         double tolerance) const;

private:

  void find_node_association(std::shared_ptr<PGeom> pgeom,
                             std::shared_ptr<percept::BlockStructuredGrid> grid,
                             bool local_tol,
                             double associate_tolerance);
  void find_edge_association(std::shared_ptr<PGeom> pgeom,
                             std::shared_ptr<percept::BlockStructuredGrid> grid);
  void find_face_association(std::shared_ptr<PGeom> pgeom,
                             std::shared_ptr<percept::BlockStructuredGrid> grid);


  void create_node_association_data(std::shared_ptr<PGeom> pgeom,
                                    std::shared_ptr<PGeomRangeSearch> range_search,
                                    std::shared_ptr<percept::StructuredBlock> sgi,
                                    std::array<unsigned, 3> range_min,
                                    std::array<unsigned, 3> range_max,
                                    bool local_tol,
                                    double associate_tolerance,
                                    unsigned block_num);
  void create_edge_association_data(std::shared_ptr<percept::StructuredBlock> sgi,
                                    std::array<unsigned, 3> range_min,
                                    std::array<unsigned, 3> range_max,
                                    enum IJK_DIRECTION edge_direction,
                                    unsigned block_num);
  void create_face_association_data(std::shared_ptr<percept::StructuredBlock> sgi,
                                    std::array<unsigned, 3> range_min,
                                    std::array<unsigned, 3> range_max,
                                    enum IJK_DIRECTION edge_direction,
                                    unsigned block_num);

  //TODO - these functions (get_*_close_to_nodes)are candidates to be common for both structured and unstructured meshes
  // instead of passing in StructuredBlock, pass in the node coordinate array for the block
  // node_map contains tolerance data for each node, so that needs to also be passed in when the tolerance is locally computed
  void get_surfaces_close_to_nodes(std::shared_ptr<percept::StructuredBlock> sgi,
                                   std::shared_ptr<PGeom> pgeom,
                                   const std::multimap< int, int > &candidate_surf_map,
                                   std::map<int, NodeAssociationData> &node_map,
                                   std::multimap< int, int > &surface_map);
  void get_curves_close_to_nodes(std::shared_ptr<percept::StructuredBlock> sgi,
                                 std::shared_ptr<PGeom> pgeom,
                                 const std::multimap< int, int > &candidate_curve_map,
                                 std::map<int, NodeAssociationData> &node_map,
                                 std::multimap< int, int > &curve_map);
  void get_vertices_close_to_nodes(std::shared_ptr<percept::StructuredBlock> sgi,
                                   std::shared_ptr<PGeom> pgeom,
                                   const std::multimap< int, int > &candidate_vertex_map,
                                   std::map<int, NodeAssociationData> &node_map,
                                   std::multimap< int, int > &vertex_map);

  // TODO - these functions (find_common_*) are candidates to be shared for both structued and unstructured meshes
  static void find_common_curves(std::vector<int> const &elem_nodes,
                          std::map<int, NodeAssociationData> &node_map,
                          std::vector<int> &common_curves);
  static void find_common_surfaces(std::vector<int> const &elem_nodes,
                            std::map<int, NodeAssociationData> &node_map,
                            std::vector<int> &common_surfs);



  static void structured_node_edge_lengths(std::shared_ptr<percept::StructuredBlock> sgi,
                         std::array<unsigned,3> indx,
                         std::vector<double> &edge_lengths);

  static double structured_node_local_tolerance(std::shared_ptr<percept::StructuredBlock> sgi,
                              std::array<unsigned,3> indx,
                              double default_tol);

  static bool find_matching_geometry(int geom_id_to_match, std::vector<CloseGeometry> &close_geometry);





// TODO - these node_coordinates functions and get_face_nodes probably belong in the StructuredBlock class
  static CubitVector node_coordinates(std::shared_ptr<percept::StructuredBlock> sgi,
                                      int node_index);
  static void gather_node_coordinates(std::shared_ptr<percept::StructuredBlock> sgi,
                                      const std::multimap<int,int>::const_iterator lower,
                                      const std::multimap<int,int>::const_iterator upper,
                                      std::vector<CubitVector> &node_positions);
static void get_face_nodes(std::shared_ptr<percept::StructuredBlock> sgi,
                           uint64_t i,
                           uint64_t j,
                           uint64_t k,
                           enum IJK_DIRECTION face_direction,
                           std::vector<int> &face_nodes);



  uint64_t refined_edge_node_id(int start_node,
                                enum IJK_DIRECTION edge_direction,
                                std::shared_ptr<percept::StructuredBlock> sg_in,
                                std::shared_ptr<percept::StructuredBlock> sg_out);
  uint64_t refined_face_node_id(int start_node,
                                enum IJK_DIRECTION face_direction,
                                std::shared_ptr<percept::StructuredBlock> sg_in,
                                std::shared_ptr<percept::StructuredBlock> sg_out);

  void make_refined_edge_mesh_association(std::shared_ptr<percept::BlockStructuredGrid> input_grid,
                                          std::shared_ptr<percept::BlockStructuredGrid> output_grid,
                                          enum IJK_DIRECTION edge_direction);
  void make_refined_face_mesh_association(std::shared_ptr<percept::BlockStructuredGrid> input_grid,
                                          std::shared_ptr<percept::BlockStructuredGrid> output_grid,
                                          enum IJK_DIRECTION face_direction);


  void print_block_geometry_association() const;
  void print_node_association(const std::map<int, NodeAssociationData> &node_to_geom_map) const;
  void print_edge_association_map(const std::map<int, EdgeAssociationData> &edge_map,
                                  const std::string &direction_string) const;
  void print_face_association_map(const std::map<int, FaceAssociationData> &face_map,
                                  const std::string &direction_string) const;
  static void print_ids(const std::string &label,
                        const std::vector<int> &id_list);
  static void print_multimap_ids(const std::string &label,
                                 const std::multimap<int,int>::const_iterator lower,
                                 const std::multimap<int,int>::const_iterator upper);
  static void print_close_geometry_ids(const std::string &label,
                                const std::vector<struct CloseGeometry> &close_geom);
  static void print_geometry_association(const std::string &key_name,
                                         const std::string &value_name,
                                         const std::string &direction_string,
                                         const std::multimap< int, int > &geometry_map);

private:
  static const std::string direction_string[3];

  std::vector<BlockAssociationData> mBlockAssociation;

  std::vector<std::map<uint64_t, NodeAssociatedGeom>> mRefinedNodeToGeomMaps;
};

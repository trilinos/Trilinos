#ifndef MESHGEOMETRY_HPP
#define MESHGEOMETRY_HPP

#include <stk_percept/PerceptMesh.hpp>
#include "GeometryKernel.hpp"

#include <boost/unordered_map.hpp>

#define DEBUG_GEOM_SNAP false


typedef std::vector<double> PointSet;
typedef int GeometryHandle;
using namespace stk;
using namespace mesh;
using namespace percept;

struct GeometryEvaluator
{
    GeometryEvaluator(Part* part) : mMesh(Selector(*part)) {}
    GeometryHandle mGeometry;
    Selector mMesh;
};

class MeshGeometry
{
public:
  typedef std::pair<int, size_t> CacheBucketClassifyValueType;
  typedef boost::unordered_map<const stk::mesh::Bucket *, CacheBucketClassifyValueType > CacheBucketClassifyType;

    MeshGeometry(GeometryKernel* geom, bool cache_bucket_selectors_is_active=false);
    ~MeshGeometry();

    void add_evaluator(GeometryEvaluator* evaluator);
    void add_evaluators(std::vector<GeometryEvaluator*> evaluators);

    // snaps all points in the mesh to their associated geometry
    void snap_points_to_geometry(PerceptMesh* mesh_data);

    // snaps only specified points in the mesh to their associated geometry
    void snap_points_to_geometry(PerceptMesh* mesh_data, std::vector<stk::mesh::Entity *>& nodes);

    // gets normal at a surface (or curve, in which case it returns the curvature vector)
    void normal_at(PerceptMesh* eMesh, stk::mesh::Entity * node, std::vector<double>& normal);

  /**
   * Return 0,1,2,3 if the node or bucket is on a geometry vertex, curve, surface or domain.
   * Return the found evaluators in the curveEvaluators and surfEvaluators.
   */
  int classify_node(const stk::mesh::Entity& node, size_t& curveOrSurfaceEvaluator);
  int classify_bucket(const stk::mesh::Bucket& bucket, size_t& curveOrSurfaceEvaluator);

    const std::vector<GeometryEvaluator*>& getGeomEvaluators();
private:

  int classify_bucket_internal(const stk::mesh::Bucket& bucket, size_t& curveOrSurfaceEvaluator);
  //int classify_bucket_internal(const stk::mesh::Bucket& bucket, std::vector<size_t>& curveEvaluators, std::vector<size_t>& surfEvaluators);

protected:
    std::vector<GeometryEvaluator*> geomEvaluators;
    GeometryKernel* geomKernel;
    CacheBucketClassifyType m_cache_bucket_classify;

public:
  bool m_cache_classify_bucket_is_active;
protected:
    void snap_point_to_geometry(stk::mesh::Entity *node);

private:

    double mDbgNodeCoords[3];

    bool contains_dbg_node( PerceptMesh *mesh_data,
                            Bucket &bucket );
    bool is_dbg_node( double node_coord[3] );

    void snap_nodes( PerceptMesh* mesh_data,
                     Bucket &bucket,
                     size_t evalautor_idx );

    void snap_node( PerceptMesh* mesh_data,
                    Entity &node,
                    size_t evalautor_idx );

    void normal_at( PerceptMesh* mesh_data,
                    Entity &node,
                    size_t evalautor_idx,
                    std::vector<double>& normal);
};

#endif // MESHGEOMETRY_HPP

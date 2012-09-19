#ifndef MESHGEOMETRY_HPP
#define MESHGEOMETRY_HPP

#include <stk_percept/PerceptMesh.hpp>
#include "GeometryKernel.hpp"

#include <boost/unordered_map.hpp>

#define DEBUG_GEOM_SNAP 0


typedef std::vector<double> PointSet;
typedef int GeometryHandle;
using namespace stk;
using namespace mesh;
using namespace percept;

struct GeometryEvaluator
{
    GeometryEvaluator(Part* part) : mMesh(*part), mPart(part) {}
    GeometryHandle mGeometry;
    Selector mMesh;
    Part *mPart;
};

class MeshGeometry
{
public:
  typedef std::pair<int, size_t> CacheBucketClassifyValueType;
  typedef boost::unordered_map<const stk::mesh::Bucket *, CacheBucketClassifyValueType > CacheBucketClassifyType;
  typedef boost::unordered_map<size_t, double> MaxDeltaOnGeometryType;

  MeshGeometry(GeometryKernel* geom, double doCheckMovement=0.0, double doCheckCpuTime=0.0, bool cache_bucket_selectors_is_active=false, bool doPrint=false);
    ~MeshGeometry();

  void print_node_movement_summary();

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
  
  // hold info for which nodes took maximum cpu time
  //struct CpuMaxInfo

private:

  int classify_bucket_internal(const stk::mesh::Bucket& bucket, size_t& curveOrSurfaceEvaluator);
  //int classify_bucket_internal(const stk::mesh::Bucket& bucket, std::vector<size_t>& curveEvaluators, std::vector<size_t>& surfEvaluators);

protected:
    std::vector<GeometryEvaluator*> geomEvaluators;
    GeometryKernel* geomKernel;
    CacheBucketClassifyType m_cache_bucket_classify;

  double m_doCheckMovement;
  double m_checkCPUTime;
  MaxDeltaOnGeometryType m_checkMovementMap;
  MaxDeltaOnGeometryType m_checkCPUTimeMap;

public:
  bool m_cache_classify_bucket_is_active;
  bool m_doPrint;
protected:
    void snap_point_to_geometry(stk::mesh::Entity *node);

private:

    double mDbgNodeCoords[3];
    int m_type;

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

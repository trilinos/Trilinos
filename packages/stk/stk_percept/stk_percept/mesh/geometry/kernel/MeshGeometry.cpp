#include "MeshGeometry.hpp"

MeshGeometry::MeshGeometry(GeometryKernel* geom, bool cache_classify_bucket_is_active)
{
  geomKernel = geom;
  mDbgNodeCoords[0] = -0.00477133907617983;
  mDbgNodeCoords[1] = -0.00477133907617983;
  mDbgNodeCoords[2] =  0.260484055257467;

  m_cache_classify_bucket_is_active = cache_classify_bucket_is_active;
}

MeshGeometry::~MeshGeometry()
{

}

void MeshGeometry::add_evaluator(GeometryEvaluator* evaluator)
{
  geomEvaluators.push_back(evaluator);
}

void MeshGeometry::add_evaluators(std::vector<GeometryEvaluator*> evaluators)
{
  geomEvaluators.insert(geomEvaluators.end(), evaluators.begin(), evaluators.end());
  // std::for_each( evaluators.begin(),
  //                evaluators.end(),
  //                std::mem_fun( &MeshGeometry::add_evaluator ) );
}

const std::vector<GeometryEvaluator*>& MeshGeometry::getGeomEvaluators()
{
  return geomEvaluators;
}

/**
 * Return 0,1,2,3 if the node is on a geometry vertex, curve, surface or domain.
 */
int MeshGeometry::classify_node(const stk::mesh::Entity& node, size_t& curveOrSurfaceEvaluator)
{
  const stk::mesh::Bucket& bucket = node.bucket();
  return classify_bucket(bucket, curveOrSurfaceEvaluator);
}

/**
 * Return 0,1,2,3 if the bucket is on a geometry vertex, curve, surface or domain.
 */
int MeshGeometry::classify_bucket(const stk::mesh::Bucket& bucket, size_t& curveOrSurfaceEvaluator)
{
  int classify_value = 0;
   if (m_cache_classify_bucket_is_active)
     {
       CacheBucketClassifyType::iterator iter = m_cache_bucket_classify.find(&bucket);
       if (iter == m_cache_bucket_classify.end())
         {
           classify_value = classify_bucket_internal(bucket, curveOrSurfaceEvaluator);
           m_cache_bucket_classify[&bucket] = CacheBucketClassifyValueType(classify_value, curveOrSurfaceEvaluator);
         }
       else
         {
           CacheBucketClassifyValueType& val = iter->second;
           classify_value = val.first;
           curveOrSurfaceEvaluator = val.second;
           bool debug=false;
           if (debug)
             {
               size_t t;
               int cv=classify_bucket_internal(bucket, t);
               if (cv != classify_value || t != curveOrSurfaceEvaluator)
                 {
                   std::cout << "cv = " << cv << " classify_value= " << classify_value << " t= " << t << " curveOrSurfaceEvaluator= " << curveOrSurfaceEvaluator << std::endl;
                   exit(123);
                 }
             }
         }
     }
   else
    {
      classify_value = classify_bucket_internal(bucket, curveOrSurfaceEvaluator);
    }
   return classify_value;
}

/**
 * Return 0,1,2,3 if the bucket is on a geometry vertex, curve, surface or domain.
 */
int MeshGeometry::classify_bucket_internal(const stk::mesh::Bucket& bucket, size_t& curveOrSurfaceEvaluator)
{
  // Each bucket contains the set of nodes with unique part intersections.
  // This means that every nodes will be in exactly one bucket.  But, the
  // nodes on curves are also in the part for the adjacent surfaces which
  // means more than one evaluator will be selected for those buckets which
  // are on the boundary of an entity (i.e. buckets representing curves
  // and vertices).  We first create a list of all the evaluators that
  // might be relevant.

  static std::vector<size_t> curveEvaluators(0);
  static std::vector<size_t> surfEvaluators(0);
  curveEvaluators.resize(0);
  surfEvaluators.resize(0);

  size_t s;
  for (s=0; s<geomEvaluators.size(); s++)
    {
      bool selector_has_bucket = geomEvaluators[s]->mMesh(bucket);
      if (selector_has_bucket)
        {
          if (geomKernel->is_curve(s))
            {
              curveEvaluators.push_back(s);
            }
          else if (geomKernel->is_surface(s))
            {
              surfEvaluators.push_back(s);
            }
        }
    }

  curveOrSurfaceEvaluator = 0;
  if ( curveEvaluators.size() > 1 )
    {
      // This is a bucket representing a vertex.  No need to do anything
      // since the node will already be on the vertex, and no new nodes
      // ever get created assigned to vertices during refinement.
      //std::cout << "Vertex node encountered" << std::endl;
      curveOrSurfaceEvaluator = 0;
      return 0;
    }
  if ( curveEvaluators.size() == 1 )
    {
      // This bucket represents a geometric curve.  Snap to it.
      //std::cout << "Snapping to curve" << curveEvaluators[0] << std::endl;
      //snap_nodes( eMesh, bucket, curveEvaluators[0] );
      curveOrSurfaceEvaluator = curveEvaluators[0];
      return 1;
    }

  if ( surfEvaluators.size() == 1 )
    {
      //std::cout << "Snapping to surface" << surfEvaluators[0] << std::endl;
      // This bucket represents a geometric surface.  Snap to it.
      //snap_nodes( eMesh, bucket, surfEvaluators[0] );
      curveOrSurfaceEvaluator = surfEvaluators[0];
      return 2;
    }
  return -1;  // error condition, or it is an interior node
}

void MeshGeometry::snap_points_to_geometry(PerceptMesh* eMesh)
{
  BulkData& bulkData = *eMesh->getBulkData();

  const std::vector<Bucket*> & buckets = bulkData.buckets( stk::mesh::fem::FEMMetaData::NODE_RANK );

  for ( std::vector<Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
  {
    Bucket& bucket = **k;
#if CONTAINS_DEBGU_NODE
    if ( contains_dbg_node( eMesh, bucket ) )
    {
      std::cout << "     DBG Node FOUND" << std::endl;
    }
#endif

    size_t curveOrSurfaceEvaluator;
    int type = classify_bucket(bucket, curveOrSurfaceEvaluator);
    switch (type) {
    case 0:
      // This is a bucket representing a vertex.  No need to do anything
      // since the node will already be on the vertex, and no new nodes
      // ever get created assigned to vertices during refinement.
      //std::cout << "Vertex node encountered" << std::endl;
      break;
    case 1:
      // This bucket represents a geometric curve.  Snap to it.
      //std::cout << "Snapping to curve" << curveEvaluators[0] << std::endl;
      snap_nodes( eMesh, bucket, curveOrSurfaceEvaluator );
      break;
    case 2:
      //std::cout << "Snapping to surface" << surfEvaluators[0] << std::endl;
      // This bucket represents a geometric surface.  Snap to it.
      snap_nodes( eMesh, bucket, curveOrSurfaceEvaluator );
      break;
    case -1:
    default:
      //printf( "ERROR: A bucket found without a geometric evaluator.\n" );
      break;
    }


  }
}

void MeshGeometry::normal_at(PerceptMesh* eMesh, stk::mesh::Entity * node, std::vector<double>& normal)
{
  {
    Bucket& bucket = node->bucket();

    // Each bucket contains the set of nodes with unique part intersections.
    // This means that every nodes will be in exactly one bucket.  But, the
    // nodes on curves are also in the part for the adjacent surfaces which
    // means more than one evaluator will be selected for those buckets which
    // are on the boundary of an entity (i.e. buckets representing curves
    // and vertices).  We first create a list of all the evaluators that
    // might be relevant.

    size_t curveOrSurfaceEvaluator;
    int type = classify_bucket(bucket, curveOrSurfaceEvaluator);
    switch (type) {
    case 0:
      // This is a bucket representing a vertex.  No need to do anything
      // since the node will already be on the vertex, and no new nodes
      // ever get created assigned to vertices during refinement.
      //std::cout << "Vertex node encountered" << std::endl;
      break;
    case 1:
      // This bucket represents a geometric curve.  Snap to it.
      //std::cout << "Snapping to curve" << curveEvaluators[0] << std::endl;
      normal_at( eMesh, *node, curveOrSurfaceEvaluator, normal );
      break;
    case 2:
      //std::cout << "Snapping to surface" << surfEvaluators[0] << std::endl;
      // This bucket represents a geometric surface.  Snap to it.
      normal_at( eMesh, *node, curveOrSurfaceEvaluator, normal );
      break;
    case -1:
    default:
      //printf( "ERROR: A bucket found without a geometric evaluator.\n" );
      break;
    }


  }
}

void MeshGeometry::snap_points_to_geometry(PerceptMesh* eMesh, std::vector<stk::mesh::Entity *>& nodes)
{
  for (unsigned inode=0; inode < nodes.size(); inode++)
  {
    Bucket& bucket = nodes[inode]->bucket();

#if CONTAINS_DEBGU_NODE
    if ( contains_dbg_node( eMesh, bucket ) )
    {
      std::cout << "     DBG Node FOUND" << std::endl;
    }
#endif

    // Each bucket contains the set of nodes with unique part intersections.
    // This means that every nodes will be in exactly one bucket.  But, the
    // nodes on curves are also in the part for the adjacent surfaces which
    // means more than one evaluator will be selected for those buckets which
    // are on the boundary of an entity (i.e. buckets representing curves
    // and vertices).  We first create a list of all the evaluators that
    // might be relevant.

    size_t curveOrSurfaceEvaluator;
    int type = classify_bucket(bucket, curveOrSurfaceEvaluator);
    switch (type) {
    case 0:
      // This is a bucket representing a vertex.  No need to do anything
      // since the node will already be on the vertex, and no new nodes
      // ever get created assigned to vertices during refinement.
      //std::cout << "Vertex node encountered" << std::endl;
      break;
    case 1:
      // This bucket represents a geometric curve.  Snap to it.
      //std::cout << "Snapping to curve" << curveEvaluators[0] << std::endl;
      snap_node( eMesh, *nodes[inode], curveOrSurfaceEvaluator );
      break;
    case 2:
      //std::cout << "Snapping to surface" << surfEvaluators[0] << std::endl;
      // This bucket represents a geometric surface.  Snap to it.
      snap_node( eMesh, *nodes[inode], curveOrSurfaceEvaluator );
      break;
    case -1:
    default:
      //printf( "ERROR: A bucket found without a geometric evaluator.\n" );
      break;
    }


  }
}

void MeshGeometry::snap_node
(
  PerceptMesh *eMesh,
  Entity & node,
  size_t evaluator_idx
)
{
  VectorFieldType* coordField = eMesh->getCoordinatesField();

  /*
  Part* new_nodes_part = eMesh->getNonConstPart("refine_new_nodes_part");
  Selector new_nodes_part_selector;
  if (new_nodes_part) new_nodes_part_selector = Selector(*new_nodes_part);
  */

  {

    double * coord = stk::mesh::field_data( *coordField , node );
    double delta[3] = {coord[0], coord[1], coord[2]};
    bool doPrint = DEBUG_GEOM_SNAP;
    if (doPrint)
    {
      std::string str = geomKernel->get_attribute(evaluator_idx);
      std::cout << "tmp geom snap_points_to_geometry eval name= " << str << " node id= " << node.identifier() 
                << " coords b4= " << coord[0] << " " << coord[1] << " " << coord[2];
    }

    if ( is_dbg_node( coord ) )
    {
      std::cout << "Node in question being projected" << std::endl;
    }

    geomKernel->snap_to(coord, geomEvaluators[evaluator_idx]->mGeometry);

    if (doPrint)
    {
      std::string str = geomKernel->get_attribute(evaluator_idx);
      delta[0] = coord[0] - delta[0];
      delta[1] = coord[1] - delta[1];
      delta[2] = coord[2] - delta[2];
      double dtot = std::sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2]);
      std::cout << " coords af= " << coord[0] << " " << coord[1] << " " << coord[2] 
                << " delta= " << delta[0] << " " << delta[1] << " " << delta[2] 
                << " deltaTot= " << dtot
                << " " << str.substr(6,5)
                << std::endl;
      //if (str.substr(6,5)=="20004") block_20004
      //{
      //  std::cout << "found 20004" << std::endl;
      //}
    }
  }
}

void MeshGeometry::normal_at
(
  PerceptMesh *eMesh,
  Entity & node,
  size_t evaluator_idx,
  std::vector<double>& normal
)
{
  VectorFieldType* coordField = eMesh->getCoordinatesField();

  {

    double * coord = stk::mesh::field_data( *coordField , node );
    bool doPrint = DEBUG_GEOM_SNAP;
    if (doPrint)
    {
      std::string str = geomKernel->get_attribute(evaluator_idx);
      std::cout << "tmp geom snap_points_to_geometry eval name= " << str << " node id= " << node.identifier() 
                << " coords b4= " << coord[0] << " " << coord[1] << " " << coord[2];
    }

    if ( is_dbg_node( coord ) )
    {
      std::cout << "Node in question being projected" << std::endl;
    }

    geomKernel->normal_at(coord, geomEvaluators[evaluator_idx]->mGeometry, normal);

    if (doPrint)
    {
      std::cout << " normal = " << normal[0] << " " << normal[1] << " " << normal[2] 
                << std::endl;
      //if (str.substr(6,5)=="20004") block_20004
      //{
      //  std::cout << "found 20004" << std::endl;
      //}
    }
  }
}

void MeshGeometry::snap_nodes
(
  PerceptMesh *eMesh,
  Bucket &bucket,
  size_t evaluator_idx
)
{
  //VectorFieldType* coordField = eMesh->getCoordinatesField();
  const unsigned num_nodes_in_bucket = bucket.size();

  //std::string str = geomKernel->get_attribute(evaluator_idx);
  for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
  {
    Entity& node = bucket[iNode];

    snap_node(eMesh, node, evaluator_idx);
  }
}

bool MeshGeometry::contains_dbg_node
(
  PerceptMesh *eMesh,
  Bucket &bucket
)
{
  VectorFieldType* coordField = eMesh->getCoordinatesField();
  const unsigned num_nodes_in_bucket = bucket.size();

  for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
  {
    Entity& node = bucket[iNode];

    double * coord = stk::mesh::field_data( *coordField , node );
    if ( is_dbg_node( coord ) )
    {
      return true;
    }
  }
  return false;
}

bool MeshGeometry::is_dbg_node( double node_coord[3] )
{
  double dx = node_coord[0] - mDbgNodeCoords[0];
  double dy = node_coord[1] - mDbgNodeCoords[1];
  double dz = node_coord[2] - mDbgNodeCoords[2];
  double dist = sqrt( (dx*dx) + (dy*dy) + (dz*dz) );
  if ( dist < 0.0001 )
  {
    return true;
  }
  return false;
}


#include "MeshGeometry.hpp"
#include <stk_util/environment/CPUTime.hpp>

MeshGeometry::MeshGeometry(GeometryKernel* geom, double doCheckMovement, double doCheckCPUTime, bool cache_classify_bucket_is_active, bool doPrint) 
  : m_doCheckMovement(doCheckMovement), m_checkCPUTime(doCheckCPUTime), m_cache_classify_bucket_is_active(cache_classify_bucket_is_active), m_doPrint(doPrint),
    m_type(-1)
{
  geomKernel = geom;
  mDbgNodeCoords[0] = -0.00477133907617983;
  mDbgNodeCoords[1] = -0.00477133907617983;
  mDbgNodeCoords[2] =  0.260484055257467;

  mDbgNodeCoords[0] = 183.049;
  mDbgNodeCoords[1] = 174.609;
  mDbgNodeCoords[2] = 27.2434;
  //big cpu = 0.462929 coords_0 = 5.60577 1.48796 19.3677 delta= -0.0963596 -0.0259963 0.00280205 coords_1= 5.50941 1.46197 19.3705 deltaMag= 0.099844 block= TOPOLOGY_SURFACE_BLOCK_20387

  mDbgNodeCoords[0] = 5.60577;
  mDbgNodeCoords[1] = 1.48796;
  mDbgNodeCoords[2] = 19.3677;

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
   m_type = classify_value;
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
#define CONTAINS_DEBUG_NODE 0
#if CONTAINS_DEBUG_NODE
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

#if CONTAINS_DEBUG_NODE
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
    double coord_0[3] = {coord[0], coord[1], coord[2]};
    bool doPrint = m_doPrint; //DEBUG_GEOM_SNAP
    bool doCheckMovement = m_doCheckMovement != 0.0;
    bool doCheckCPUTime = m_checkCPUTime != 0.0;

    if ( is_dbg_node( coord ) )
    {
      std::cout << "Node in question being projected" << std::endl;
      doPrint=true;
    }

    if (doPrint)
    {
      std::string str = geomKernel->get_attribute(evaluator_idx);
      std::cout << "tmp geom snap_points_to_geometry eval name= " << str << " node id= " << node.identifier() 
                << " coords b4= " << coord[0] << " " << coord[1] << " " << coord[2] << " type= " << m_type << std::endl;
    }

    double cpu0 = 0.0, cpu1 = 0.0;
    if (doCheckCPUTime) cpu0 = stk::cpu_time();

    // Look at the neighboring edges and calculate an average edge length.  This will be 
    // used as a tolerance to tell the projection code that if the projected point
    // has moved less than this value we can consider it a valid solution.  This will
    // greatly reduce the number of iterations the projection code has to do.
    double edge_length_ave=0.0;
    // get edge lengths
    stk::mesh::PairIterRelation node_elements = node.relations(eMesh->element_rank());
    for(unsigned ii=0; ii < node_elements.size(); ii++)
    {
      edge_length_ave += eMesh->edge_length_ave(*node_elements[ii].entity());
    }
    edge_length_ave /= ((double)node_elements.size());

    geomKernel->snap_to(coord, geomEvaluators[evaluator_idx]->mGeometry, &edge_length_ave);

    if (doCheckCPUTime) cpu1 = stk::cpu_time() - cpu0;

    delta[0] = coord[0] - delta[0];
    delta[1] = coord[1] - delta[1];
    delta[2] = coord[2] - delta[2];

    if (doPrint || (doCheckCPUTime && cpu1 > m_checkCPUTime))
      {
        if (m_checkCPUTimeMap.find(evaluator_idx) == m_checkCPUTimeMap.end())
          m_checkCPUTimeMap[evaluator_idx] = 0.0;

        m_checkCPUTimeMap[evaluator_idx] = std::max(m_checkCPUTimeMap[evaluator_idx], cpu1);

        bool print_every = true;
        if (print_every)
          {
            std::string str = geomKernel->get_attribute(evaluator_idx);
            double dtot = std::sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2]);

            std::cout << "big cpu = " << cpu1 << " coords_0 = " << coord_0[0] << " " << coord_0[1] << " " << coord_0[2] 
                      << " delta= " << delta[0] << " " << delta[1] << " " << delta[2] 
                      << " coords_1= " << coord[0] << " " << coord[1] << " " << coord[2]
                      << " deltaMag= " << dtot
                      << " block= " << str
                      << std::endl;
          }
      }
    
    if (doPrint || doCheckMovement)
    {
      std::string str = geomKernel->get_attribute(evaluator_idx);
      double dtot = std::sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2]);
      if(doPrint)
        std::cout << " coords af= " << coord[0] << " " << coord[1] << " " << coord[2] 
                  << " delta= " << delta[0] << " " << delta[1] << " " << delta[2] 
                  << " deltaTot= " << dtot
                  << " block= " << str
                  << std::endl;
      if (doCheckMovement)
        {
          if (m_doCheckMovement >= 0.0)
            {
              edge_length_ave = m_doCheckMovement;
            }
          double factor = 1.0;
          if (dtot > edge_length_ave*factor)
            {
              if (m_checkMovementMap.find(evaluator_idx) == m_checkMovementMap.end())
                m_checkMovementMap[evaluator_idx] = 0.0;
              m_checkMovementMap[evaluator_idx] = std::max(m_checkMovementMap[evaluator_idx], dtot);
              const bool print_every = true;
              if (print_every)
                {
                  std::cout << "big delta coords_0 = " << coord_0[0] << " " << coord_0[1] << " " << coord_0[2] 
                            << " delta= " << delta[0] << " " << delta[1] << " " << delta[2] 
                            << " coords_1= " << coord[0] << " " << coord[1] << " " << coord[2]
                            << " deltaMag= " << dtot
                            << " block= " << str
                            << std::endl;
                }

            }
        }

      //if (str.substr(6,5)=="20004") block_20004
      //{
      //  std::cout << "found 20004" << std::endl;
      //}
    }
  }
}

void MeshGeometry::print_node_movement_summary()
{
  //if (m_checkMovementMap.size()) std::cout << "MeshGeometry::print_node_movement_summary, size= " << m_checkMovementMap.size() << std::endl;
  for (MaxDeltaOnGeometryType::iterator iter = m_checkMovementMap.begin(); iter != m_checkMovementMap.end(); ++iter)
    {
      std::cout << "MeshGeometry::print_node_movement_summary, delta= " << iter->second << " on geometry= " << iter->first << " = " << geomKernel->get_attribute(iter->first) << std::endl;
    }
  for (MaxDeltaOnGeometryType::iterator iter = m_checkCPUTimeMap.begin(); iter != m_checkCPUTimeMap.end(); ++iter)
    {
      std::cout << "MeshGeometry::print_node_movement_summary, cpu= " << iter->second << " on geometry= " << iter->first << " = " << geomKernel->get_attribute(iter->first) << std::endl;
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
  if ( dist < 0.001 )
  {
    return true;
  }
  return false;
}


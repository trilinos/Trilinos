// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/mesh/geometry/kernel/MeshGeometry.hpp>
#include <stk_util/environment/CPUTime.hpp>

#include <percept/mesh/geometry/kernel/GeometryKernelGregoryPatch.hpp>
#include <percept/PerceptMesh.hpp>


namespace percept {

GeometryEvaluator::GeometryEvaluator(stk::mesh::Part* part) : mGeometry(-1,INVALID,""), mMesh(*part), mPart(part)  {}
  
MeshGeometry::MeshGeometry(const PerceptMesh& eMesh, GeometryKernel* geom, double doCheckMovement, double doCheckCPUTime, bool cache_classify_bucket_is_active, bool doPrint)
  : m_eMesh(eMesh), m_doCheckMovement(doCheckMovement), m_checkCPUTime(doCheckCPUTime), m_cache_classify_bucket_is_active(cache_classify_bucket_is_active), m_doPrint(doPrint),
    m_type(-1), m_doCurvesOnly(false)
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

  mDbgNodeCoords[0] = 10.5;
  mDbgNodeCoords[1] = -1.5;
  mDbgNodeCoords[2] = 1.5;

  mDbgNodeCoords[0] = 9.;
  mDbgNodeCoords[1] = 0.;
  mDbgNodeCoords[2] = 1.5;

  // mDbgNodeCoords[0] = 0.0968663;
  // mDbgNodeCoords[1] = 0.00277812;
  // mDbgNodeCoords[2] = -0.136322;
   // mDbgNodeCoords[0] = 0.106339;
   // mDbgNodeCoords[1] = 0.00598256;
   // mDbgNodeCoords[2] = -0.136442 ;
}

MeshGeometry::~MeshGeometry()
{
  for (unsigned i = 0; i < geomEvaluators.size(); i++)
    {
      if (geomEvaluators[i])
        delete geomEvaluators[i];
      geomEvaluators[i] = 0;
    }
}

void MeshGeometry::add_evaluator(GeometryEvaluator* evaluator)
{
  geomEvaluators.push_back(evaluator);
}

void MeshGeometry::add_evaluators(std::vector<GeometryEvaluator*> evaluators)
{
  geomEvaluators.insert(geomEvaluators.end(), evaluators.begin(), evaluators.end());
}

const std::vector<GeometryEvaluator*>& MeshGeometry::getGeomEvaluators()
{
  return geomEvaluators;
}

/**
 * Return 0,1,2,3 if the node is on a geometry vertex, curve, surface or domain.
 */
int MeshGeometry::classify_node(const stk::mesh::Entity node, GeometryHandle & curveOrSurfaceEvaluator/*size_t& curveOrSurfaceEvaluator*/)
{
  const stk::mesh::Bucket& bucket = m_eMesh.bucket(node);
  return classify_bucket(bucket, curveOrSurfaceEvaluator);
}

/**
 * Return 0,1,2,3 if the bucket is on a geometry vertex, curve, surface or domain.
 */
int MeshGeometry::classify_bucket(const stk::mesh::Bucket& bucket, GeometryHandle & curveOrSurfaceEvaluator/*size_t& curveOrSurfaceEvaluator*/)
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
int MeshGeometry::classify_bucket_internal(const stk::mesh::Bucket& bucket, GeometryHandle & curveOrSurfaceEvaluator/*size_t& curveOrSurfaceEvaluator*/)
{
  // Each bucket contains the set of nodes with unique part intersections.
  // This means that every nodes will be in exactly one bucket.  But, the
  // nodes on curves are also in the part for the adjacent surfaces which
  // means more than one evaluator will be selected for those buckets which
  // are on the boundary of an entity (i.e. buckets representing curves
  // and vertices).  We first create a list of all the evaluators that
  // might be relevant.

  static std::vector<GeometryHandle> curveEvaluators(0);
  static std::vector<GeometryHandle> surfEvaluators(0);
  curveEvaluators.resize(0);
  surfEvaluators.resize(0);

  size_t s;
  for (s=0; s<geomEvaluators.size(); s++)
    {
      bool selector_has_bucket = geomEvaluators[s]->mMesh(bucket);
      if (selector_has_bucket)
        {
          if (geomKernel->is_curve(geomEvaluators[s]->mGeometry))
            {
              curveEvaluators.push_back(geomEvaluators[s]->mGeometry);
            }
          else if (geomKernel->is_surface(geomEvaluators[s]->mGeometry))
            {
              surfEvaluators.push_back(geomEvaluators[s]->mGeometry);
            }
          // bcarnes: what if neither is true? can that happen? should we throw an error?
        }
    }

  curveOrSurfaceEvaluator.m_id = 0;
  if ( curveEvaluators.size() > 1 )
    {
      // This is a bucket representing a vertex.  No need to do anything
      // since the node will already be on the vertex, and no new nodes
      // ever get created assigned to vertices during refinement.
      curveOrSurfaceEvaluator.m_id = 0;
      return 0;
    }
  if ( curveEvaluators.size() == 1 )
    {
      // This bucket represents a geometric curve.  Snap to it.
      curveOrSurfaceEvaluator = curveEvaluators[0];
      return 1;
    }

  if ( surfEvaluators.size() == 1 )
    {
      // This bucket represents a geometric surface.  Snap to it.
      curveOrSurfaceEvaluator = surfEvaluators[0];
      return 2;
    }
  return -1;  // error condition, or it is an interior node
}

void MeshGeometry::classify_bucket_all(const stk::mesh::Bucket& bucket, std::set<size_t>& curveEvaluators, std::set<size_t>& surfEvaluators)
{
  curveEvaluators.clear();
  surfEvaluators.clear();

  size_t s;
  for (s=0; s<geomEvaluators.size(); s++)
    {
      bool selector_has_bucket = geomEvaluators[s]->mMesh(bucket);
      if (selector_has_bucket)
        {
          if (geomKernel->is_curve(geomEvaluators[s]->mGeometry))
            {
              curveEvaluators.insert(s);
            }
          else if (geomKernel->is_surface(geomEvaluators[s]->mGeometry))
            {
              surfEvaluators.insert(s);
            }
        }
    }
}

void MeshGeometry::pre_process(PerceptMesh* eMesh)
{
  stk::mesh::PartVector partv;
  for (unsigned ii=0; ii < geomEvaluators.size(); ++ii)
    {
      if (geomEvaluators[ii]->mPart->name() == "edgeseams")
        continue;
      partv.push_back(geomEvaluators[ii]->mPart);
    }

  if (geomKernel) geomKernel->pre_process(eMesh, partv);
}


void MeshGeometry::snap_points_to_geometry(PerceptMesh* eMesh)
{
//  stk::mesh::BulkData& bulkData = *eMesh->get_bulk_data();
  pre_process(eMesh);
  stk::mesh::FieldBase* coordField = eMesh->get_coordinates_field();





  	std::map<stk::mesh::Entity, std::list<GeometryEvaluator*>> nodeMap;

	for (size_t iEval = 0; iEval < geomEvaluators.size(); iEval++) {
		stk::mesh::Selector elementSelector = *(geomEvaluators[iEval]->mPart);
		const stk::mesh::BucketVector & bucketsOfElems =
				elementSelector.get_buckets(stk::topology::NODE_RANK);

		for (size_t iBucket = 0; iBucket < bucketsOfElems.size(); iBucket++) {
			stk::mesh::Bucket & bucketOfElems = *(bucketsOfElems[iBucket]);
			const unsigned numElemsInBucket = bucketOfElems.size();

			for (unsigned iNode = 0; iNode < numElemsInBucket; iNode++) {

				stk::mesh::Entity node = bucketOfElems[iNode];



				std::map<stk::mesh::Entity,std::list<GeometryEvaluator*>>::iterator nodeIter;
				nodeIter = nodeMap.find(node);
				if(nodeIter==nodeMap.end()){
					std::list<GeometryEvaluator*> interimList;
					interimList.push_back(geomEvaluators[iEval]);
					nodeMap.insert( std::pair<stk::mesh::Entity,std::list<GeometryEvaluator*>>(node,interimList) );
				}
				else{
					(*nodeIter).second.push_front(geomEvaluators[iEval]);
				}




			}

		}
	}





  size_t s;
  std::vector<GeometryHandle> evaluators, curve_evals, surf_evals;
  for(std::map<stk::mesh::Entity,std::list<GeometryEvaluator*>>::iterator nodeMapIter=nodeMap.begin();nodeMapIter!=nodeMap.end();nodeMapIter++ )
  {




	  std::pair<stk::mesh::Entity,
	  			std::list<GeometryEvaluator*>> nodeToEval;
	  nodeToEval.first=(*nodeMapIter).first;
	  nodeToEval.second=(*nodeMapIter).second;


	  for(std::list<GeometryEvaluator*>::const_iterator evalIter = nodeToEval.second.begin(); evalIter!=nodeToEval.second.end();evalIter++){
		  if((*evalIter)->mGeometry.m_type==percept::GeomEvalType::CURVE)
			 curve_evals.push_back((*evalIter)->mGeometry);
		  else if((*evalIter)->mGeometry.m_type==percept::GeomEvalType::SURFACE)
			  surf_evals.push_back((*evalIter)->mGeometry);
	  }




    // favor lower order evaluators if they exist
    if(surf_evals.size() > 0)
      evaluators = surf_evals;
    if(curve_evals.size() > 0)
      evaluators = curve_evals;

    // bcarnes: can we avoid specific code for one Geometry type?
    if (m_doCurvesOnly && typeid(*geomKernel) == typeid(GeometryKernelGregoryPatch))
      {
        evaluators.clear();
        if (curve_evals.size() > 0 && surf_evals.size() > 0)
          {
            evaluators = surf_evals;
          }
      }

    // If we have more than one evaluator we will project to each and
    // keep the best one.
    int spatialDim = eMesh->get_spatial_dim();
    if(evaluators.size() > 1)
    {
      double * coords = static_cast<double*>(eMesh->field_data(*coordField , nodeToEval.first));
      double orig_pos[3] = {coords[0], coords[1], (spatialDim==3 ? coords[2] : 0) };
      double smallest_dist_sq = std::numeric_limits<double>::max();
      double best_pos[3] = {0,0,0};
//      std::cout<< evaluators.size() << std::endl;
      for(s=0; s<evaluators.size(); s++)
      {
        // Always start with the original position.
        for(int f=0; f < spatialDim; f++)
          coords[f] = orig_pos[f];
        // Do the projection (changes "coords" variable).
        snap_node(eMesh, nodeToEval.first, evaluators[s]);
        // See if this projection is closer to the original point
        // than any of the previous ones.
        double dist_sq = 0.0;
        for(int f=0; f < spatialDim; f++)
          dist_sq += (orig_pos[f] - coords[f]) * (orig_pos[f] - coords[f]);
        if(dist_sq < smallest_dist_sq)
        {
          smallest_dist_sq = dist_sq;
          for(int f=0; f < spatialDim; f++){

        	  best_pos[f] = coords[f];

          }

        }
      }
      // Load the best projection into the actual coordinates.
      for(int f=0; f < spatialDim; f++)
        coords[f] = best_pos[f];
    }
    // If we know we only have one evaluator don't bother with all
    // of the saving/restoring of positions.
    else if(evaluators.size() > 0)
      snap_node(eMesh, nodeToEval.first, evaluators[0]);

    surf_evals.clear();
    curve_evals.clear();
    evaluators.clear();
  }
}

void MeshGeometry::normal_at(PerceptMesh* eMesh, stk::mesh::Entity node, std::vector<double>& normal)
{
  {
    stk::mesh::Bucket& bucket = eMesh->bucket(node);

    // Each bucket contains the set of nodes with unique part intersections.
    // This means that every nodes will be in exactly one bucket.  But, the
    // nodes on curves are also in the part for the adjacent surfaces which
    // means more than one evaluator will be selected for those buckets which
    // are on the boundary of an entity (i.e. buckets representing curves
    // and vertices).  We first create a list of all the evaluators that
    // might be relevant.

    GeometryHandle curveOrSurfaceEvaluator;
    int type = classify_bucket(bucket, curveOrSurfaceEvaluator);
    switch (type) {
    case 0:
      // This is a bucket representing a vertex.  No need to do anything
      // since the node will already be on the vertex, and no new nodes
      // ever get created assigned to vertices during refinement.
      break;
    case 1:
      // This bucket represents a geometric curve.  Snap to it.
      normal_at( eMesh, node, curveOrSurfaceEvaluator, normal );
      break;
    case 2:
      // This bucket represents a geometric surface.  Snap to it.
      normal_at( eMesh, node, curveOrSurfaceEvaluator, normal );
      break;
    case -1:
    default:
      //printf( "ERROR: A bucket found without a geometric evaluator.\n" );
      break;
    }
  }
}

void MeshGeometry::point_at(PerceptMesh* eMesh, stk::mesh::Entity node, std::vector<double>& coords, bool use_node_coords)
{
  {
    stk::mesh::Bucket& bucket = eMesh->bucket(node);

    GeometryHandle curveOrSurfaceEvaluator;
    int type = classify_bucket(bucket, curveOrSurfaceEvaluator);
    switch (type) {
    case 0:
      // This is a bucket representing a vertex.  No need to do anything
      // since the node will already be on the vertex, and no new nodes
      // ever get created assigned to vertices during refinement.
      break;
    case 1:
      // This bucket represents a geometric curve.  Snap to it.
      point_at( eMesh, node, curveOrSurfaceEvaluator, coords, use_node_coords );
      break;
    case 2:
      // This bucket represents a geometric surface.  Snap to it.
      point_at( eMesh, node, curveOrSurfaceEvaluator, coords, use_node_coords );
      break;
    case -1:
    default:
      //printf( "ERROR: A bucket found without a geometric evaluator.\n" );
      break;
    }
  }
}

void MeshGeometry::snap_points_to_geometry(PerceptMesh* eMesh, std::vector<stk::mesh::Entity>& nodes)
{
  for (unsigned inode=0; inode < nodes.size(); inode++)
  {
    stk::mesh::Bucket& bucket = eMesh->bucket(nodes[inode]);

    // Each bucket contains the set of nodes with unique part intersections.
    // This means that every nodes will be in exactly one bucket.  But, the
    // nodes on curves are also in the part for the adjacent surfaces which
    // means more than one evaluator will be selected for those buckets which
    // are on the boundary of an entity (i.e. buckets representing curves
    // and vertices).  We first create a list of all the evaluators that
    // might be relevant.

    GeometryHandle curveOrSurfaceEvaluator;
    int type = classify_bucket(bucket, curveOrSurfaceEvaluator);
    switch (type) {
    case 0:
      // This is a bucket representing a vertex.  No need to do anything
      // since the node will already be on the vertex, and no new nodes
      // ever get created assigned to vertices during refinement.
      break;
    case 1:
      // This bucket represents a geometric curve.  Snap to it.
      snap_node( eMesh, nodes[inode], curveOrSurfaceEvaluator );
      break;
    case 2:
      // This bucket represents a geometric surface.  Snap to it.
      snap_node( eMesh, nodes[inode], curveOrSurfaceEvaluator );
      break;
    case -1:
    default:
      //printf( "ERROR: A bucket found without a geometric evaluator.\n" );
      break;
    }
  }
}

void MeshGeometry::snap_node(PerceptMesh *eMesh, stk::mesh::Entity node, GeometryHandle geomHand /*size_t evaluator_idx*/)
{
  stk::mesh::FieldBase* coordField = eMesh->get_coordinates_field();

  {
    int spatialDim = eMesh->get_spatial_dim();
    double * coord_in = static_cast<double*>(eMesh->field_data( *coordField , node ));
    double coord[3] = {coord_in[0], coord_in[1], (spatialDim==3?coord_in[2]:0)};
    double delta[3] = {coord[0], coord[1], coord[2]};

    // Look at the neighboring edges and calculate an average edge length.  This will be
    // used as a tolerance to tell the projection code that if the projected point
    // has moved less than this value we can consider it a valid solution.  This will
    // greatly reduce the number of iterations the projection code has to do.
    double edge_length_ave=0.0;
    // get edge lengths
    const MyPairIterRelation node_elements(*eMesh, node, eMesh->element_rank() );
    for(unsigned ii=0; ii < node_elements.size(); ii++)
    {
      edge_length_ave += eMesh->edge_length_ave(node_elements[ii].entity());
    }
    edge_length_ave /= ((double)node_elements.size());

    double *crd = &coord[0];
    geomKernel->snap_to(crd, geomHand /*geomEvaluators[evaluator_idx]->mGeometry*/, &edge_length_ave, 0, 0, &node);
    coord_in[0] = coord[0];
    coord_in[1] = coord[1];
    if (spatialDim == 3)
      coord_in[2] = coord[2];

    delta[0] = coord[0] - delta[0];
    delta[1] = coord[1] - delta[1];
    delta[2] = coord[2] - delta[2];
  }
}

void MeshGeometry::print_node_movement_summary()
{
//  for (MaxDeltaOnGeometryType::iterator iter = m_checkMovementMap.begin(); iter != m_checkMovementMap.end(); ++iter)
//    {
//      std::cout << "MeshGeometry::print_node_movement_summary, delta= " << iter->second << " on geometry= " << iter->first << " = " << geomKernel->get_attribute(iter->first) << std::endl;
//    }
//  for (MaxDeltaOnGeometryType::iterator iter = m_checkCPUTimeMap.begin(); iter != m_checkCPUTimeMap.end(); ++iter)
//    {
//      std::cout << "MeshGeometry::print_node_movement_summary, cpu= " << iter->second << " on geometry= " << iter->first << " = " << geomKernel->get_attribute(iter->first) << std::endl;
//    }
}

void MeshGeometry::normal_at(PerceptMesh *eMesh, stk::mesh::Entity node, GeometryHandle & curveOrSurfaceEvaluator/*size_t& curveOrSurfaceEvaluator*//*size_t evalautor_idx*/, std::vector<double>& normal)
{
  stk::mesh::FieldBase* coordField = eMesh->get_coordinates_field();

  {
    double * coord = static_cast<double*>(eMesh->field_data( *coordField , node ));

    geomKernel->normal_at(coord, curveOrSurfaceEvaluator, /*geomEvaluators[evaluator_idx]->mGeometry,*/ normal, (void *)(&node));
  }
}

void MeshGeometry::point_at(PerceptMesh *eMesh, stk::mesh::Entity node, GeometryHandle & curveOrSurfaceEvaluator/*size_t& curveOrSurfaceEvaluator*//*size_t evalautor_idx*/, std::vector<double>& coords, bool use_node_coords)
{
  stk::mesh::FieldBase* coordField = eMesh->get_coordinates_field();

  {
    double * coord_in = static_cast<double*>(eMesh->field_data( *coordField , node ));
    if (!use_node_coords)
      {
        coord_in = &coords[0];
      }
    int spatialDim = eMesh->get_spatial_dim();
    double coord[3] = {coord_in[0], coord_in[1], (spatialDim==3?coord_in[2]:0)};

    double *coord_out = &coords[0];

    // get edge lengths
    double edge_length_ave=0.0;
    const MyPairIterRelation node_elements(*eMesh, node, eMesh->element_rank() );
    for(unsigned ii=0; ii < node_elements.size(); ii++)
      {
        edge_length_ave += eMesh->edge_length_ave(node_elements[ii].entity());
      }
    edge_length_ave /= ((double)node_elements.size());

    double *crd = &coord[0];
    geomKernel->snap_to(crd, curveOrSurfaceEvaluator /*geomEvaluators[evaluator_idx]->mGeometry*/, &edge_length_ave, 0, 0, &node);
    coord_out[0] = coord[0];
    coord_out[1] = coord[1];
    if (spatialDim == 3)
      coord_out[2] = coord[2];
  }
}

//void MeshGeometry::snap_nodes(PerceptMesh *eMesh, stk::mesh::Bucket &bucket, size_t evaluator_idx)
//{
//  const unsigned num_nodes_in_bucket = bucket.size();
//
//  for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
//  {
//    stk::mesh::Entity node = bucket[iNode];
//    snap_node(eMesh, node, evaluator_idx);
//  }
//}

bool MeshGeometry::contains_dbg_node
(
  PerceptMesh *eMesh,
  stk::mesh::Bucket &bucket
)
{
  stk::mesh::FieldBase* coordField = eMesh->get_coordinates_field();
  const unsigned num_nodes_in_bucket = bucket.size();

  for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
  {
    stk::mesh::Entity node = bucket[iNode];

    double * coord = static_cast<double*>(eMesh->field_data( *coordField , node ));
    if ( is_dbg_node( coord ) )
    {
      return true;
    }
  }
  return false;
}

bool MeshGeometry::is_dbg_node( double node_coord[3] , stk::mesh::Entity *node, PerceptMesh *eMesh)
{
  stk::mesh::EntityId nodeIdToFind = 0;
  // e.g. nodeIdToFind = 55239;
  if (node && eMesh && eMesh->id(*node) == nodeIdToFind)
    return true;

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

}

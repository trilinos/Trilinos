// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef MESHGEOMETRY_HPP
#define MESHGEOMETRY_HPP


#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <percept/mesh/geometry/kernel/GeometryKernel.hpp>

#include <unordered_map>

namespace percept {

class PerceptMesh;

typedef std::vector<double> PointSet;

struct GeometryEvaluator
{
  //MAY NEED BETTER DESTRUCTOR
  //MAY NEED COPY CONSTRUCTOR
  ~GeometryEvaluator(){};
  GeometryEvaluator(stk::mesh::Part* part);
  GeometryHandle mGeometry;
  stk::mesh::Selector mMesh;
  stk::mesh::Part *mPart;
  
  //ADD BOOLEAN OPERATORS:      <     AND   ==
};


struct nodeToEvals {

	std::pair<stk::mesh::Entity, std::list<GeometryEvaluator*>> n2e;
	//MAY NEED DESTUCTOR
	nodeToEvals(){};

    STK_FUNCTION
    bool operator==(nodeToEvals entityToEvals) const { return n2e.first.m_value == entityToEvals.n2e.first.m_value; }

//    STK_FUNCTION
//    bool operator==(entity_value_type val) const { return m_value == val; }

    STK_FUNCTION
    bool operator!=(nodeToEvals entityToEvals) const { return n2e.first.m_value != entityToEvals.n2e.first.m_value; }

//    STK_FUNCTION
//    bool operator!=(entity_value_type val) const { return m_value != val; }

    STK_FUNCTION
    bool operator<(nodeToEvals entityToEvals) const { return n2e.first.m_value < entityToEvals.n2e.first.m_value; }

};

struct nodeToEvaluatorsList
{


	nodeToEvaluatorsList() {};
	std::list<std::pair<
		stk::mesh::Entity,
		std::list<GeometryEvaluator*> >> n2e;

	size_t size() {return n2e.size();}

	bool isIn(std::pair<stk::mesh::Entity,
		std::list<GeometryEvaluator*>> toCheck, std::list<std::pair<stk::mesh::Entity,
		std::list<GeometryEvaluator*>>>::iterator & index){
		//returns an iterator containing position of location within n2e
		for (index=n2e.begin();index!=n2e.end();index++/*size_t i=0;i<n2e.size();i++*/)
			if (toCheck.first==(*index).first){
//				n2e.begin();
				return true;
			}
//		n2e.begin();
		return false;
	}
	bool insert(std::pair<stk::mesh::Entity,
		std::list<GeometryEvaluator*>> toCheck)
	{
		std::list<std::pair<stk::mesh::Entity,
			std::list<GeometryEvaluator*>>>::iterator pos;
//		size_t orgSize = 0;
//		size_t newSize = 0;

		if(isIn(toCheck, pos)){//if you find an entity you're equal to ...
//			iterator = intList.begin(); iterator != intList.end(); ++iterator
//			orgSize = n2e[index].mEvalsWithIndex.size();
//			newSize = n2e[index].mEvalsWithIndex.size();
//			std::cout<<"Got through isIn just fine. Starting Iteration through my geomveallist..." << std::endl;
			for(std::list<GeometryEvaluator*>::iterator iterPos = pos->second.begin();iterPos!=pos->second.end();iterPos++
			/*size_t i=0;i<toCheck.mEvalsWithIndex.size();i++*/)
				// ... check to see that your evaluator's part names don't coincide

			{
//				size_t j=0;
//				std::cout<<(*iterPos)->mPart->name() <<std::endl;
				bool toInsert=true;
				std::list<GeometryEvaluator*>::iterator iterToCheck;
//				std::cout<<"the size of toCheck's evaluators is: " <<toCheck.second.size() <<std::endl;
				for( iterToCheck = toCheck.second.begin();iterToCheck!=toCheck.second.end();iterToCheck++
				/*;j<n2e[index].mEvalsWithIndex.size();j++*/){
//					std::cout<<(*iterToCheck)->mPart->name() <<std::endl;
					if ( (*iterToCheck)->mPart->name()==(*iterPos)->mPart->name()
							/*toCheck.mEvalsWithIndex[i]->mPart->name()==n2e[index].mEvalsWithIndex[j]->mPart->name()*/){ //if you find the part, stop looking
						toInsert=false;
						break;}
				} //loop seems fine
				iterToCheck--;
//				std::cout<<"Made it through toCheck loop" <<std::endl;
				if (toInsert /*j==n2e[index].mEvalsWithIndex.size()*/){ //if you didn't find the part, add it to the vector
					//segfaults before it gets here
//					std::cout<<"About to insert geom evaluator with part" << (*iterToCheck)->mPart->name() <<std::endl;
					(*pos).second.push_back(*iterToCheck); //Pretty sure it's messing up here

//					n2e[index].mEvalsWithIndex.push_back(toCheck.mEvalsWithIndex[i]);
//					newSize++;
				}
//				*iterToCheck=NULL;
//				*iterPos=NULL;
//				delete iterToCheck;
//				delete iterPos;
			}



//			if(newSize>orgSize)
//				return true;

			return true;
		}



		n2e.push_back(toCheck);
		return true;
	}//endinsert

};


class GeometryKernel;


class MeshGeometry
{

public:
  const PerceptMesh& m_eMesh;
  //PerceptMesh * m_eMesh_pntr;
  typedef std::pair<int, GeometryHandle> CacheBucketClassifyValueType;
  typedef std::unordered_map<const stk::mesh::Bucket *, CacheBucketClassifyValueType > CacheBucketClassifyType;

  MeshGeometry(const PerceptMesh& eMesh, GeometryKernel* geom, double doCheckMovement=0.0, double doCheckCpuTime=0.0, bool cache_bucket_selectors_is_active=false, bool doPrint=false);
    ~MeshGeometry();

  void print_node_movement_summary();

    void add_evaluator(GeometryEvaluator* evaluator);
    void add_evaluators(std::vector<GeometryEvaluator*> evaluators);



    // snaps all points in the mesh to their associated geometry
    void snap_points_to_geometry(PerceptMesh* mesh_data);

    // snaps only specified points in the mesh to their associated geometry
    void snap_points_to_geometry(PerceptMesh* mesh_data, std::vector<stk::mesh::Entity>& nodes);

    // gets normal at a surface (or curve, in which case it returns the curvature vector)
    void normal_at(PerceptMesh* eMesh, stk::mesh::Entity node, std::vector<double>& normal);

    // find projection of point to the surface -
    //    If @param use_node_coords is false, @param coords is both the input and output coordinates
    //       and the @param node's coordinates are ignored, it is only used for geometry classification
    //    If @param use_node_coords is true, @param coord is only output, the input coordinates
    //       come from the node
    //    In either case, the node's coordinates are not modified
    void point_at(PerceptMesh* eMesh, stk::mesh::Entity node, std::vector<double>& coords, bool use_node_coords=true);

    void pre_process(PerceptMesh *eMesh);

  /**
   * Return 0,1,2,3 if the node or bucket is on a geometry vertex, curve, surface or domain.
   * Return the found evaluators in the curveEvaluators and surfEvaluators.
   */
  int classify_node(const stk::mesh::Entity node, GeometryHandle & curveOrSurfaceEvaluator/*size_t& curveOrSurfaceEvaluator*/);
  int classify_bucket(const stk::mesh::Bucket& bucket, GeometryHandle & curveOrSurfaceEvaluator/*size_t& curveOrSurfaceEvaluator*/);

  // for more fine-grained access by applications - return all classifications of the bucket
  void classify_bucket_all(const stk::mesh::Bucket& bucket, std::set<size_t>& curveEvaluators, std::set<size_t>& surfEvaluators);

  const std::vector<GeometryEvaluator*>& getGeomEvaluators();

  // hold info for which nodes took maximum cpu time
  //struct CpuMaxInfo

private:

  int classify_bucket_internal(const stk::mesh::Bucket& bucket, GeometryHandle & curveOrSurfaceEvaluator/*size_t& curveOrSurfaceEvaluator*/);
  //int classify_bucket_internal(const stk::mesh::Bucket& bucket, std::vector<size_t>& curveEvaluators, std::vector<size_t>& surfEvaluators);

protected:
    std::vector<GeometryEvaluator*> geomEvaluators;
public:
    GeometryKernel* geomKernel;
protected:
    CacheBucketClassifyType m_cache_bucket_classify;

  double m_doCheckMovement;
  double m_checkCPUTime;
//  MaxDeltaOnGeometryType m_checkMovementMap;
//  MaxDeltaOnGeometryType m_checkCPUTimeMap;

public:
  bool m_cache_classify_bucket_is_active;
  bool m_doPrint;
protected:
    void snap_point_to_geometry(stk::mesh::Entity node);

private:

    double mDbgNodeCoords[3];
    int m_type;

    bool contains_dbg_node( PerceptMesh *mesh_data,
                            stk::mesh::Bucket &bucket );
    bool is_dbg_node( double node_coord[3] , stk::mesh::Entity *node=0, PerceptMesh *eMesh=0);

    void snap_nodes( PerceptMesh* mesh_data,
                     stk::mesh::Bucket &bucket,
                     size_t evalautor_idx );

    void snap_node( PerceptMesh* mesh_data,
                    stk::mesh::Entity node,
					GeometryHandle geomHand /*size_t evalautor_idx*/ );

public:

    void normal_at( PerceptMesh* mesh_data,
                    stk::mesh::Entity node,
					GeometryHandle & curveOrSurfaceEvaluator/*size_t& curveOrSurfaceEvaluator*//*size_t evalautor_idx*/,
                    std::vector<double>& normal);


    void point_at( PerceptMesh *eMesh,
                   stk::mesh::Entity node,
				   GeometryHandle & curveOrSurfaceEvaluator/*size_t& curveOrSurfaceEvaluator*//*size_t evaluator_idx*/,
                   std::vector<double>& coords,
                   bool use_node_coords);

public:
  bool m_doCurvesOnly;
};

}
#endif // MESHGEOMETRY_HPP

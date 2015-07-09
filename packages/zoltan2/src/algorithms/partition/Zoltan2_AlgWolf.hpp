// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//                    Michael Wolf (mmwolf@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

//MMW need to specify that this requires Zoltan                  

#ifndef _ZOLTAN2_ALGWOLF_HPP_
#define _ZOLTAN2_ALGWOLF_HPP_

#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_AlgZoltan.hpp>

#include <sstream>
#include <string>
#include <bitset>

/*! \file Zoltan2_AlgWolf.hpp
 *  \brief The algorithm for Wolf partitioning.
 */


void buildPartTree(int level, int leftPart, int splitPart, int rightPart, std::vector<int> &partTree);


namespace Zoltan2
{

// /*! \brief The boolean parameters of interest to the Block algorithm.
//  */
// enum blockParams{
//   block_balanceCount,            /*!< objective = balance_object_count */
//   block_balanceWeight,          /*!< objective = balance_object_weight */
//   block_minTotalWeight,      /*!< objective = mc_minimize_total_weight */
//   block_minMaximumWeight,  /*!< objective = mc_minimize_maximum_weight */
//   block_balanceTotalMaximum, /*!< objective = mc_balance_total_maximum */
//   NUM_BLOCK_PARAMS
// };

////////////////////////////////////////////////////////////////////////////////
/*! Wolf partitioning method.
 *
 *  \param env   library configuration and problem parameters
 *  \param problemComm  the communicator for the problem
 *  \param ids    an Identifier model
 *
 *  Preconditions: The parameters in the environment have been
 *    processed (committed).  No special requirements on the
 *    identifiers.
 *
 */
////////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
class AlgWolf : public Algorithm<Adapter>
{

private:

  typedef typename Adapter::part_t part_t;

  const RCP<const Environment> mEnv;
  const RCP<Comm<int> > mProblemComm;

  const RCP<const GraphModel<typename Adapter::base_adapter_t> > mGraphModel;
  const RCP<const CoordinateModel<typename Adapter::base_adapter_t> > mIds;

  const RCP<const typename Adapter::base_adapter_t> mBaseInputAdapter;


  void getBoundLayerSep(int levelIndx, const std::vector<part_t> &partMap,
			const part_t * parts, 
			std::vector<int> &boundVerts,
			std::vector<std::vector<int> > &boundVertsST,
			const std::set<int> &sepVerts);

public:
  // Constructor
  AlgWolf(const RCP<const Environment> &env_,
	  const RCP<Comm<int> > &problemComm_,
	  const RCP<const GraphModel<typename Adapter::base_adapter_t> > &gModel_,
	  const RCP<const CoordinateModel<typename Adapter::base_adapter_t> > &cModel_,
	  const RCP<const typename Adapter::base_adapter_t> baseInputAdapter_)
    :mEnv(env_), mProblemComm(problemComm_), mGraphModel(gModel_), mIds(cModel_), 
     mBaseInputAdapter(baseInputAdapter_)
  {
#ifndef INCLUDE_ZOLTAN2_EXPERIMENTAL
    Z2_THROW_EXPERIMENTAL("Zoltan2 Wolf is strictly experimental software ")
#endif

#ifndef INCLUDE_ZOLTAN2_EXPERIMENTAL_WOLF
    Z2_THROW_EXPERIMENTAL_WOLF("Zoltan2 Wolf is strictly experimental software ")
#endif

    if(mProblemComm->getSize()!=1)
    {
      Z2_THROW_SERIAL("Zoltan2 Wolf is strictly serial!");
    }

  }

  // Partitioning method
  void partition(const RCP<PartitioningSolution<Adapter> > &solution_);

};
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void AlgWolf<Adapter>::partition(const RCP<PartitioningSolution<Adapter> > &solution_)
{
    // typedef typename Adapter::lno_t lno_t;     // local ids
    // typedef typename Adapter::gno_t gno_t;     // global ids
    // typedef typename Adapter::scalar_t scalar_t;   // scalars

    mEnv->debug(DETAILED_STATUS, std::string("Entering AlgWolf"));

    //////////////////////////////////////////////////////////////////////
    // First, let's partition with RCB using Zoltan.  Eventually, we will change this
    // to use PHG
    //////////////////////////////////////////////////////////////////////

    // Q: can I use solution passed into alg or do I need to create a different one?
    //    For now using the one passed into alg
    // TODO: use new partitioning solution


    AlgZoltan<Adapter> algZoltan(this->mEnv, mProblemComm, this->mBaseInputAdapter);
    algZoltan.partition(solution_);

    size_t numGlobalParts = solution_->getTargetGlobalNumberOfParts();

    const part_t *parts = solution_->getPartListView();
    //////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////
    // Build up tree that represents partitioning subproblems, which will 
    // be used for determining separators at each level
    //   -- for now, this is built up artificially
    //   -- eventually this will be built from PHG output
    //////////////////////////////////////////////////////////////////////
    // change int to something, part_t?
    std::vector<int> partTree;

    buildPartTree( 0, 0, (numGlobalParts-1)/2 + 1, numGlobalParts, partTree);
    unsigned int numSeparators = partTree.size() / 4;
    //////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////
    // Create a map that maps each part number to a new number based on
    // the level of the hiearchy of the separator tree.  This allows us
    // to easily identify the boundary value vertices
    //////////////////////////////////////////////////////////////////////
    int numLevels = partTree[4*(numSeparators-1)]+1;

    std::vector<std::vector<int> > partLevelMap(numLevels,std::vector<int>(numGlobalParts));

    std::vector<int> sepsInLev(numLevels,0);

    for(unsigned int i=0;i<numSeparators;i++)
    {
      int level = partTree[4*i];
      int leftPart = partTree[4*i+1];
      int splitPart = partTree[4*i+2];
      int rightPart = partTree[4*i+3];
      
      for(int part=leftPart; part<splitPart; part++)
      {
        partLevelMap[level][part] = 2*sepsInLev[level];
      }

      for(int part=splitPart; part<rightPart; part++)
      {
        partLevelMap[level][part] = 2*sepsInLev[level]+1;
      }

      sepsInLev[level]++;
    }
    //////////////////////////////////////////////////////////////////////

    // Set of separator vertices.  Used to keep track of what vertices are
    // already in previous calculated separators.  These vertices should be
    // excluded from future separator calculations
    const std::set<int> sepVerts;

    //////////////////////////////////////////////////////////////////////
    // Loop over each cut
    //    1. Build boundary layer between parts
    //    2. Build vertex separator from boundary layer
    //////////////////////////////////////////////////////////////////////
    for(unsigned int level=0;level<numLevels;level++)
    {
      for(unsigned int levIndx=0;levIndx<sepsInLev[level];levIndx++)
      {

	std::vector<int> boundVerts;
	std::vector<std::vector<int> > boundVertsST(2);

        ///////////////////////////////////////////////////////////////
        // Build boundary layer between parts (edge separator)
        ///////////////////////////////////////////////////////////////
        getBoundLayerSep(levIndx, partLevelMap[level], parts, boundVerts,
			 boundVertsST, sepVerts);
        ///////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////
        // Calculate vertex separator from boundary layer
        ///////////////////////////////////////////////////////////////

	//VCOfBoundLayer

        ///////////////////////////////////////////////////////////////


	}
    }
    //////////////////////////////////////////////////////////////////////

    //TODO: calculate vertex separator for each layer, 
    //TODO: using vertex separators, compute new ordering and store in solution
    //TODO: move to ordering directory

    mEnv->debug(DETAILED_STATUS, std::string("Exiting AlgWolf"));
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Create boundary layer of vertices between 2 partitions
////////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void AlgWolf<Adapter>::getBoundLayerSep(int levelIndx, const std::vector<part_t> &partMap,
					const part_t * parts, 
					std::vector<int> &boundVerts,
					std::vector<std::vector<int> > &boundVertsST,
					const std::set<int> &sepVerts)
{
  typedef typename Adapter::lno_t lno_t;         // local ids
  typedef typename Adapter::scalar_t scalar_t;   // scalars
  typedef StridedData<lno_t, scalar_t> input_t;

  int numVerts = mGraphModel->getLocalNumVertices();

  //Teuchos ArrayView
  ArrayView< const lno_t > eIDs;
  ArrayView< const lno_t > vOffsets;
  ArrayView< const lno_t > procIDs;
  ArrayView< input_t > wgts;

  // For some reason getLocalEdgeList seems to be returning empty eIDs
  //size_t numEdges = ( (GraphModel<typename Adapter::base_adapter_t>)  *mGraphModel).getLocalEdgeList(eIDs, vOffsets, wgts);

  size_t numEdges = ( (GraphModel<typename Adapter::base_adapter_t>)  *mGraphModel).getEdgeList(eIDs, procIDs, vOffsets, wgts);

//   size_t Zoltan2::GraphModel< Adapter >::getEdgeList(ArrayView< const gno_t > & edgeIds,
// 						     ArrayView< const int > & procIds,
// 						     ArrayView< const lno_t > & offsets,
// 						     ArrayView< input_t > & wgts 
//						     )

  //  for(int v1=0;v1<numEdges;v1++)
  for(int v1=0;v1<numVerts;v1++)
  {

    part_t vpart1 = partMap[parts[v1]];

    bool correctBL = (vpart1 >= 2*levelIndx && vpart1 < 2*(levelIndx+1) );

    // if vertex is not in the correct range of parts, it cannot be a member of 
    // this boundary layer
    if(!correctBL)
    {
      continue;
    }

    // If this vertex belongs to a previous separator, it cannot belong to this
    // separator
    if(sepVerts.find(v1)!=sepVerts.end())
    {
      continue;
    }

    //Loop over edges connected to v1
    //MMW figure out how to get this from Zoltan2
    for(int j=vOffsets[v1];j<vOffsets[v1+1];j++)
    {

      int v2 = eIDs[j];

      part_t vpart2 = partMap[parts[v2]];

      correctBL = (vpart2 >= 2*levelIndx && vpart2 < 2*(levelIndx+1) );

      // if vertex is not in the correct range of parts, it cannot be a member of 
      // this boundary layer
      if(!correctBL)
      {
        continue;
      }

      // If this vertex belongs to a previous separator, it cannot belong to this
      // separator
      if(sepVerts.find(v2)!=sepVerts.end())
      {
        continue;
      }

      if ( vpart1 !=  vpart2  )
      {
        // Vertex added to set of all boundary vertices
        boundVerts.push_back(v1);

        // Vertex added to 1st set of boundary vertices
	if(vpart1<vpart2)
        {
	  boundVertsST[0].push_back(v1);
	}
        // Vertex added to 2nd set of boundary vertices
	else
	{
	  boundVertsST[1].push_back(v1);
	}
	break;
      }

    }
  }

}
//////////////////////////////////////////////////////////////////////////////

}   // namespace Zoltan2


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void buildPartTree(int level, int leftPart, int splitPart, int rightPart, std::vector<int> &partTree)
{
  // Insert information for this separator
  partTree.push_back(level);
  partTree.push_back(leftPart);
  partTree.push_back(splitPart);
  partTree.push_back(rightPart);

  // Recurse down left side of tree
  if(splitPart-leftPart > 1)
  {
    int newSplit = leftPart+(splitPart-leftPart-1)/2 + 1;
    buildPartTree(level+1,leftPart,newSplit,splitPart,partTree);
  }

  // Recurse down right side of tree
  if(rightPart-splitPart>1)
  {
    int newSplit = splitPart+(rightPart-splitPart-1)/2 + 1;
    buildPartTree(level+1,splitPart,newSplit,rightPart,partTree);
  }
}
////////////////////////////////////////////////////////////////////////////////




#endif

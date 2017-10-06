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
//
// ***********************************************************************
//
// @HEADER

//MMW need to specify that this requires Zoltan                  

#ifndef _ZOLTAN2_ALGND_HPP_
#define _ZOLTAN2_ALGND_HPP_

#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_AlgZoltan.hpp>

#include <Zoltan2_MatcherHelper.hpp>


#include <sstream>
#include <string>
#include <bitset>

/*! \file Zoltan2_AlgND.hpp
 *  \brief The algorithm for ND based ordering
 */


namespace Zoltan2
{

////////////////////////////////////////////////////////////////////////////////
/*! Nested dissection based ordering method.
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
class AlgND : public Algorithm<typename Adapter::base_adapter_t>
//class AlgND : public Algorithm<Adapter>
{

private:

  typedef typename Adapter::part_t part_t;

  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;


  const RCP<const Environment> mEnv;
  const RCP<const Comm<int> > mProblemComm;

  //  const RCP<const GraphModel<Adapter> > mGraphModel;
  const RCP<GraphModel<typename Adapter::base_adapter_t> > mGraphModel;
  //  const RCP<const CoordinateModel<Adapter> > mIds;
  const RCP<CoordinateModel<typename Adapter::base_adapter_t> > mIds;

  //const RCP<const Adapter> mBaseInputAdapter;
  //const RCP<const Adapter> mInputAdapter;
  const RCP<const typename Adapter::base_adapter_t> mBaseInputAdapter;                                                                                                                                 

  void getBoundLayer(int levelIndx, const std::vector<part_t> &partMap,
		     const part_t * parts, 
		     const std::set<int> &excVerts,
		     int &bigraphNumS, int &bigraphNumT, int &bigraphNumE,
		     std::vector<int> &bigraphCRSRowPtr, std::vector<int> &bigraphCRSCols,
	             std::vector<int> &bigraphVMapU, std::vector<int> &bigraphVMapV);

  void buildPartTree(int level, int leftPart, int splitPart, int rightPart, 
                     std::vector<unsigned int> &sepLevels, 
		     std::vector<std::set<part_t> > &sepVerts1, std::vector<std::set<part_t> > &sepVerts2,
                     int &maxLev);


public:
  // Constructor
  AlgND(const RCP<const Environment> &env_,
        const RCP<const Comm<int> > &problemComm_,
	const RCP<GraphModel<typename Adapter::base_adapter_t> > &gModel_,
	const RCP<CoordinateModel<typename Adapter::base_adapter_t> > &cModel_,
	const RCP<const typename Adapter::base_adapter_t> baseInputAdapter_
       )
    :mEnv(env_), mProblemComm(problemComm_), mGraphModel(gModel_), 
     mIds(cModel_), mBaseInputAdapter(baseInputAdapter_)
  {
#ifndef INCLUDE_ZOLTAN2_EXPERIMENTAL
    Z2_THROW_EXPERIMENTAL("Zoltan2 AlgND is strictly experimental software ")
#endif

#ifndef INCLUDE_ZOLTAN2_EXPERIMENTAL_WOLF
    Z2_THROW_EXPERIMENTAL_WOLF("Zoltan2 algND is strictly experimental software ")
#endif

    if(mProblemComm->getSize()!=1)
    {
      Z2_THROW_SERIAL("Zoltan2 AlgND is strictly serial!");
    }

  }

  // Ordering method
  int localOrder(const RCP<LocalOrderingSolution<lno_t> > &solution_);
  int globalOrder(const RCP<GlobalOrderingSolution<gno_t> > &solution_);

};
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
int AlgND<Adapter>::globalOrder(
  const RCP<GlobalOrderingSolution<gno_t> > &solution_)
{
  throw std::logic_error("AlgND does not support global ordering.");
}

////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
int AlgND<Adapter>::localOrder(const RCP<LocalOrderingSolution<lno_t> > &solution_)
{
    typedef typename Adapter::lno_t lno_t;     // local ids
    // typedef typename Adapter::gno_t gno_t;     // global ids
    // typedef typename Adapter::scalar_t scalar_t;   // scalars

    mEnv->debug(DETAILED_STATUS, std::string("Entering AlgND"));

    //////////////////////////////////////////////////////////////////////
    // First, let's partition with RCB using Zoltan.  Eventually, we will change 
    // to give an option to use PHG
    //////////////////////////////////////////////////////////////////////

    RCP<PartitioningSolution<Adapter> > partSoln;
    int nUserWts=0;

    partSoln =
      RCP<PartitioningSolution<Adapter> > (new PartitioningSolution<Adapter>(this->mEnv, mProblemComm, nUserWts));

    AlgZoltan<Adapter> algZoltan(this->mEnv, mProblemComm, this->mBaseInputAdapter);

    algZoltan.partition(partSoln);

    size_t numGlobalParts = partSoln->getTargetGlobalNumberOfParts();

    const part_t *parts = partSoln->getPartListView();

    int numVerts = mGraphModel->getLocalNumVertices();
    for(size_t i=0; i< numVerts; i++)
    {
      std::cout << "part[" << i << "] = " << parts[i] <<std::endl;
    }
    //////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////
    // Build up tree that represents partitioning subproblems, which will 
    // be used for determining separators at each level
    //   -- for now, this is built up artificially
    //   -- eventually this will be obtained from PHG output
    //
    // Each separator i is represented by 1 integer and two sets of part_t's in the 
    // the following 3 structure:  
    //             
    //             sepLevels[i] - level of separator
    //             sepVerts1[i] - 1st set of vertices on one side of separator
    //             sepVerts2[i] - 2nd set of vertices on other side of separator
    // 
    //////////////////////////////////////////////////////////////////////
    std::vector<unsigned int> sepLevels;
    std::vector<std::set<part_t> > sepVerts1;
    std::vector<std::set<part_t> > sepVerts2;

    int maxLevel = 0;

    buildPartTree( 0, 0, (numGlobalParts-1)/2 + 1, numGlobalParts, 
                   sepLevels, sepVerts1, sepVerts2, maxLevel);
    unsigned int numSeparators = sepLevels.size();

    std::cout << "NumSeparators: " << numSeparators << std::endl;
    std::cout << "Max Level: " << maxLevel << std::endl;

    for(unsigned int i=0;i<sepLevels.size(); i++)
    {
      std::cout << sepLevels[i] << " Separator Vs: ";
      for(typename std::set<part_t>::const_iterator iter = sepVerts1[i].begin();
          iter!=sepVerts1[i].end();++iter)
      {
	std::cout << *iter << " ";
      }
      std::cout << " - ";

      for(typename std::set<part_t>::const_iterator iter = sepVerts2[i].begin();
          iter!=sepVerts2[i].end();++iter)
      {
	std::cout << *iter << " ";
      }


      std::cout  << std::endl;
    }
    //////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////
    // Create a map that maps each part number to a new number based on
    // the level of the hiearchy of the separator tree.  This allows us
    // to easily identify the boundary value vertices
    //////////////////////////////////////////////////////////////////////
    int numLevels = maxLevel+1;

    std::vector<std::vector<int> > partLevelMap(numLevels,std::vector<int>(numGlobalParts));

    std::vector<int> sepsInLev(numLevels,0);

    for(unsigned int i=0;i<numSeparators;i++)
    {
      int level = sepLevels[i];

      for(typename std::set<part_t>::const_iterator iter = sepVerts1[i].begin();
          iter!=sepVerts1[i].end();++iter)
      {
	partLevelMap[level][*iter] = 2*sepsInLev[level];
      }


      for(typename std::set<part_t>::const_iterator iter = sepVerts2[i].begin();
          iter!=sepVerts2[i].end();++iter)
      {
	partLevelMap[level][*iter] = 2*sepsInLev[level]+1;
      }





      sepsInLev[level]++;
    }

    std::cout << "partLevelMap[0][0] = " << partLevelMap[0][0] << std::endl; 
    std::cout << "partLevelMap[0][1] = " << partLevelMap[0][1] << std::endl; 

    //////////////////////////////////////////////////////////////////////

    // Set of separator vertices.  Used to keep track of what vertices are
    // already in previous calculated separators.  These vertices should be
    // excluded from future separator calculations
    std::set<lno_t> sepVerts;
    std::vector<std::vector< std::set<lno_t> > > sepVertsByLevel(numLevels);

    //////////////////////////////////////////////////////////////////////
    // Loop over each cut
    //    1. Build boundary layer between parts
    //    2. Build vertex separator from boundary layer
    //////////////////////////////////////////////////////////////////////
    for(unsigned int level=0;level<numLevels;level++)
    {
      sepVertsByLevel[level].resize(sepsInLev[level]);

      for(unsigned int levIndx=0;levIndx<sepsInLev[level];levIndx++)
      {
        ///////////////////////////////////////////////////////////////
        // Build boundary layer between parts (edge separator)
        ///////////////////////////////////////////////////////////////
        int bigraphNumU=0, bigraphNumV=0, bigraphNumE=0;
	std::vector<int> bigraphVMapU; 
        std::vector<int> bigraphVMapV;

	std::vector<int> bigraphCRSRowPtr;
	std::vector<int> bigraphCRSCols;


        getBoundLayer(levIndx, partLevelMap[level], parts, sepVerts,
		      bigraphNumU,bigraphNumV,bigraphNumE,
		      bigraphCRSRowPtr, bigraphCRSCols,
		      bigraphVMapU,bigraphVMapV);

	std::cout << "Bipartite graph: " << bigraphNumU << " " << bigraphNumV << " " 
		  << bigraphNumE << std::endl;

        for (unsigned int i=0;i<bigraphVMapU.size();i++)
	{
	  std::cout << "boundVertU: " << bigraphVMapU[i] << std::endl;
        }

        for (unsigned int i=0;i<bigraphVMapV.size();i++)
	{
	  std::cout << "boundVertV: " << bigraphVMapV[i] << std::endl;
        }



        for (int rownum=0;rownum<bigraphNumU;rownum++)
	{

           for (int eIdx=bigraphCRSRowPtr[rownum];eIdx<bigraphCRSRowPtr[rownum+1];eIdx++)
	   {          
	      std::cout << "bipartite E: " << bigraphVMapU[rownum] << ", " << bigraphVMapV[ bigraphCRSCols[eIdx]]
			<< " ( "  << rownum << "," << bigraphCRSCols[eIdx] << " )" << std::endl;
           }

	}
        ///////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////
        // Calculate bipartite matching from boundary layer
        ///////////////////////////////////////////////////////////////
        if (bigraphNumU > 0)
	{
          assert(bigraphNumV>0);

	  Matcher<lno_t> bpMatch(bigraphCRSRowPtr.data(), bigraphCRSCols.data(), bigraphNumU, bigraphNumV, bigraphNumE);
          bpMatch.match();

	  const std::vector<int> &vertUMatches = bpMatch.getVertexUMatches();
	  const std::vector<int> &vertVMatches = bpMatch.getVertexVMatches();
          ///////////////////////////////////////////////////////////////

          ///////////////////////////////////////////////////////////////
          // Calculate vertex cover (which is vertex separator) from matching
          ///////////////////////////////////////////////////////////////
	  std::vector<int> VC;

          bpMatch.getVCfromMatching(bigraphCRSRowPtr,bigraphCRSCols,vertUMatches,vertVMatches,
	 		  bigraphVMapU,bigraphVMapV,VC);

          for(unsigned int i=0;i<VC.size();i++)
	  {
            sepVerts.insert(VC[i]);

            sepVertsByLevel[level][levIndx].insert(VC[i]);
	    std::cout << "VC: " << VC[i] << std::endl;
	  }        
          ///////////////////////////////////////////////////////////////
	}

        //TODO: Copy data into separator structures?


      }
    }

    //////////////////////////////////////////////////////////////////////
    // Output separators
    //////////////////////////////////////////////////////////////////////
    std::cout << "Separators: " << std::endl;
    for(unsigned int level=0;level<sepVertsByLevel.size();level++)
    {
      sepVertsByLevel[level].resize(sepsInLev[level]);

      for(unsigned int levIndx=0;levIndx<sepVertsByLevel[level].size();levIndx++)
      {
	std::cout << "  Separator " << level << " " << levIndx << ": ";




	typename std::set<lno_t>::const_iterator iterS;
	for (iterS=sepVertsByLevel[level][levIndx].begin();iterS!=sepVertsByLevel[level][levIndx].end();++iterS)
	{
	  std::cout << *iterS << " ";
	}
	std::cout << std::endl;


      }
    }
    //////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////

    // //TODO: calculate vertex separator for each layer, 
    // //TODO: using vertex separators, compute new ordering and store in solution
    // //TODO: move to ordering directory

    mEnv->debug(DETAILED_STATUS, std::string("Exiting AlgND"));
    return 0;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Create boundary layer of vertices between 2 partitions
//
// Could improve the efficiency here by first creating a boundary layer graph
// between all parts
////////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void AlgND<Adapter>::getBoundLayer(int levelIndx, const std::vector<part_t> &partMap,
				   const part_t * parts,
				   const std::set<int> &excVerts,
				   int &bigraphNumS, int &bigraphNumT, int &bigraphNumE,
				   std::vector<int> &bigraphCRSRowPtr, std::vector<int> &bigraphCRSCols,
				   std::vector<int> &bigraphVMapS, std::vector<int> &bigraphVMapT)
{
  typedef typename Adapter::lno_t lno_t;         // local ids
  typedef typename Adapter::offset_t offset_t;   // offset_t
  typedef typename Adapter::scalar_t scalar_t;   // scalars
  typedef StridedData<lno_t, scalar_t> input_t;

  int numVerts = mGraphModel->getLocalNumVertices();

  //Teuchos ArrayView
  ArrayView< const lno_t > eIDs;
  ArrayView< const offset_t > vOffsets;
  ArrayView< input_t > wgts;

  // For some reason getLocalEdgeList seems to be returning empty eIDs
  //size_t numEdges = ( (GraphModel<typename Adapter::base_adapter_t>)  *mGraphModel).getLocalEdgeList(eIDs, vOffsets, wgts);

  //size_t numEdges = ( (GraphModel<typename Adapter::base_adapter_t>)  *mGraphModel).getEdgeList(eIDs, vOffsets, wgts);
  ( (GraphModel<typename Adapter::base_adapter_t>)  *mGraphModel).getEdgeList(eIDs, vOffsets, wgts);


  std::map<int,std::set<int> > bigraphEs;
  std::set<int> vSetS;
  std::set<int> vSetT;

  bigraphNumE=0;

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

    // Ignore vertices that belong to set of vertices to exclude
    if(excVerts.find(v1)!=excVerts.end())
    {
      continue;
    }

    //Loop over edges connected to v1
    //MMW figure out how to get this from Zoltan2
    for(offset_t j=vOffsets[v1];j<vOffsets[v1+1];j++)
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

      // Ignore vertices that belong to set of vertices to exclude
      if(excVerts.find(v2)!=excVerts.end())
      {
        continue;
      }

      if ( vpart1 !=  vpart2  )
      {
        // Vertex added to 1st set of boundary vertices
	if(vpart1<vpart2)
        {
          vSetS.insert(v1);

          // v1, v2          
          if(bigraphEs.find(v1)==bigraphEs.end())
	  {
            bigraphEs[v1] = std::set<int>();
	  }
          bigraphEs[v1].insert(v2);
          bigraphNumE++;

	}
        // Vertex added to 2nd set of boundary vertices
	else
	{
          vSetT.insert(v1);
	}
      }

    }
  }

  /////////////////////////////////////////////////////////////////////////
  // Set size of two vertex sets for bipartite graph
  /////////////////////////////////////////////////////////////////////////
  bigraphNumS = vSetS.size();
  bigraphNumT = vSetT.size();
  /////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  bigraphVMapS.resize(bigraphNumS);

  std::map<int,int> glob2LocTMap;

  unsigned int indx=0;
  for(std::set<int>::const_iterator iter=vSetS.begin(); iter!=vSetS.end(); ++iter)
  {
    bigraphVMapS[indx] = *iter;
    indx++;
  }


  bigraphVMapT.resize(bigraphNumT);
  indx=0;
  for(std::set<int>::const_iterator iter=vSetT.begin();iter!=vSetT.end();++iter)
  {
    bigraphVMapT[indx] = *iter;
    glob2LocTMap[*iter]=indx;
    indx++;
  }
  /////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////
  // Set sizes for bipartite graph data structures
  /////////////////////////////////////////////////////////////////////////
  bigraphCRSRowPtr.resize(bigraphNumS+1);
  bigraphCRSCols.resize(bigraphNumE);
  /////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // Copy bipartite graph edges into CRS format
  /////////////////////////////////////////////////////////////////////////
  bigraphCRSRowPtr[0]=0;

  unsigned int rownum=0;
  unsigned int nzIndx=0;
  std::map<int,std::set<int> >::const_iterator iterM;
  for (iterM=bigraphEs.begin();iterM!=bigraphEs.end();++iterM)
  {
    bigraphCRSRowPtr[rownum+1] = bigraphCRSRowPtr[rownum] + (*iterM).second.size();

    for(std::set<int>::const_iterator iter=(*iterM).second.begin(); iter!=(*iterM).second.end(); ++iter)
    {
      bigraphCRSCols[nzIndx] = glob2LocTMap[(*iter)];

      nzIndx++;
    }
 
    rownum++;
  }
  /////////////////////////////////////////////////////////////////////////

}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void AlgND<Adapter>::
buildPartTree(int level, int leftPart, int splitPart, int rightPart, 
	      std::vector<unsigned int> &sepLevels,
	      std::vector<std::set<part_t> > &sepVerts1, std::vector<std::set<part_t> > &sepVerts2,
              int &maxLev)
{
  // Insert information for this separator
  maxLev=level;

  //////////////////////////////////////////////////////////////////////
  // Insert information for this separator
  //////////////////////////////////////////////////////////////////////
  sepLevels.push_back(level);

  sepVerts1.push_back(std::set<part_t>());
  typename std::vector<std::set<part_t> >::iterator iter = sepVerts1.end();
  iter--;
  for(part_t v= leftPart; v < splitPart; v++)
  {
    (*iter).insert(v);
  }


  sepVerts2.push_back(std::set<part_t>());
  iter = sepVerts2.end();
  iter--;
  for(part_t v= splitPart; v < rightPart; v++)
  {
    (*iter).insert(v);
  }
  //////////////////////////////////////////////////////////////////////


  int tmpMaxLev=maxLev;

  // Recurse down left side of tree
  if(splitPart-leftPart > 1)
  {

    int newSplit = leftPart+(splitPart-leftPart-1)/2 + 1;
    buildPartTree(level+1,leftPart,newSplit,splitPart,sepLevels,sepVerts1,sepVerts2,tmpMaxLev);

    if(tmpMaxLev>maxLev)
    {
      maxLev = tmpMaxLev;
    }

  }

  // Recurse down right side of tree
  if(rightPart-splitPart>1)
  {
    int newSplit = splitPart+(rightPart-splitPart-1)/2 + 1;
    buildPartTree(level+1,splitPart,newSplit,rightPart,sepLevels,sepVerts1,sepVerts2,tmpMaxLev);

    if(tmpMaxLev>maxLev)
    {
      maxLev = tmpMaxLev;
    }

  }
}
//////////////////////////////////////////////////////////////////////////////



}   // namespace Zoltan2





#endif

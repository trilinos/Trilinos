// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  class AlgND : public Algorithm<Adapter>
  {

  private:
    typedef typename Adapter::part_t part_t;

    typedef typename Adapter::lno_t lno_t;
    typedef typename Adapter::gno_t gno_t;

    const RCP<const Environment> mEnv;
    const RCP<const Comm<int>> mProblemComm;

    std::string mPartitionMethod;

    const RCP<const typename Adapter::base_adapter_t> mBaseInputAdapter;
    modelFlag_t graphFlags;

    void getBoundLayer(part_t levelIndx, const std::vector<part_t> &partMap,
                       const part_t *parts, const std::set<lno_t> &excVerts,
                       lno_t &bigraphNumS, lno_t &bigraphNumT, lno_t &bigraphNumE,
                       std::vector<lno_t> &bigraphCRSRowPtr, std::vector<lno_t> &bigraphCRSCols,
                       std::vector<lno_t> &bigraphVMapU, std::vector<lno_t> &bigraphVMapV,
                       const RCP<const GraphModel<Adapter>> &graphModel);

    void buildPartTree(part_t level, std::vector<part_t> &levIndx,
                       part_t startV, const std::vector<part_t> &permPartNums,
                       const std::vector<part_t> &splitRangeBeg,
                       const std::vector<part_t> &splitRangeEnd,
                       const std::vector<part_t> &treeVertParents,
                       std::vector<part_t> &sepLevels, std::vector<std::set<part_t>> &sepParts1,
                       std::vector<std::set<part_t>> &sepParts2, part_t &maxLev,
                       part_t &sepTreeIndx, part_t *sepTreeView,
                       std::vector<std::pair<part_t, part_t>> &treeIndxToSepLev);

    void fillSolutionIperm(const RCP<const GraphModel<Adapter>> &graphModel,
                           const part_t *parts, const std::set<lno_t> &sepVerts,
                           const std::vector<std::vector<std::set<lno_t>>> &sepVertsByLevel,
                           const std::vector<std::pair<part_t, part_t>> &treeIndxToSepLev,
                           lno_t *ipermView, lno_t *sepRangeView);

    void getIdsOfPart(const RCP<const GraphModel<Adapter>> &graphModel, const part_t *parts,
                      part_t targetPart, const std::set<lno_t> &idsToExcl, std::set<lno_t> &outIds);

  public:
    // Constructor
    AlgND(const RCP<const Environment> &env_,
          const RCP<const Comm<int>> &problemComm_,
          const RCP<const typename Adapter::base_adapter_t> &baseInputAdapter_,
          const modelFlag_t &graphFlags_)
        : mEnv(env_), mProblemComm(problemComm_), mPartitionMethod("rcb"),
          mBaseInputAdapter(baseInputAdapter_), graphFlags(graphFlags_)
    {
#ifndef INCLUDE_ZOLTAN2_EXPERIMENTAL
      Z2_THROW_EXPERIMENTAL("Zoltan2 AlgND is strictly experimental software ")
#endif

      if (mProblemComm->getSize() != 1)
      {
        Z2_THROW_SERIAL("Zoltan2 AlgND is strictly serial!");
      }

      const Teuchos::ParameterList &pl = mEnv->getParameters();
      const Teuchos::ParameterEntry *pe;

      pe = pl.getEntryPtr("edge_separator_method");

      if (pe)
      {
        mPartitionMethod = pe->getValue<std::string>(&mPartitionMethod);
      }
    }

    // Ordering method
    int localOrder(const RCP<LocalOrderingSolution<lno_t>> &solution_);
    int globalOrder(const RCP<GlobalOrderingSolution<gno_t>> &solution_);

    /*! \brief Set up validators specific to this algorithm
   */
    static void getValidParameters(ParameterList &pl)
    {

      RCP<Teuchos::StringValidator> es_method_Validator =
          Teuchos::rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("rcb", "phg")));

      pl.set("edge_separator_method", "rcb", "ND ordering - Edge separator method", es_method_Validator);
    }
  };
  ////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  template <typename Adapter>
  int AlgND<Adapter>::globalOrder(
      const RCP<GlobalOrderingSolution<gno_t>> &solution_)
  {
    throw std::logic_error("AlgND does not support global ordering.");
  }

  ////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  template <typename Adapter>
  int AlgND<Adapter>::localOrder(const RCP<LocalOrderingSolution<lno_t>> &solution_)
  {
    mEnv->debug(DETAILED_STATUS, std::string("Entering AlgND"));

    //////////////////////////////////////////////////////////////////////
    // First, let's partition with RCB using Zoltan.  Eventually, we will change
    // to give an option to use PHG
    //////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    // Create parameter list for partitioning environment
    /////////////////////////////////////////////////////////////////
    Teuchos::ParameterList partParams;

    part_t numParts = mEnv->getParameters().template get<part_t>("num_global_parts");

    partParams.set("num_global_parts", numParts);

    // Keeping partitioning tree
    partParams.set("keep_partition_tree", true);

    // Set Zoltan parameter lists
    Teuchos::ParameterList &zparams = partParams.sublist("zoltan_parameters", false);
    zparams.set("LB_METHOD", mPartitionMethod);
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    //Create new environment with parameters for partitioning
    /////////////////////////////////////////////////////////////////
    const RCP<const Environment> partEnv = rcp(new Zoltan2::Environment(partParams, this->mEnv->comm_));
    /////////////////////////////////////////////////////////////////

    int nUserWts = 0;

    RCP<AlgZoltan<Adapter>> algZoltan = rcp(new AlgZoltan<Adapter>(partEnv, mProblemComm, this->mBaseInputAdapter));

    RCP<PartitioningSolution<Adapter>> partSoln;
    partSoln =
        RCP<PartitioningSolution<Adapter>>(new PartitioningSolution<Adapter>(partEnv, mProblemComm, nUserWts, algZoltan));

    algZoltan->partition(partSoln);

    size_t numGlobalParts = partSoln->getTargetGlobalNumberOfParts();

    const part_t *parts = partSoln->getPartListView();

    // size_t numVerts = mGraphModel->getLocalNumVertices();
    // for(size_t i=0; i< numVerts; i++)
    // {
    //   std::cout << "part[" << i << "] = " << parts[i] <<std::endl;
    // }
    //////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////
    // Obtain partitioning tree info from solution
    //////////////////////////////////////////////////////////////////////

    // Need to guarantee partitioning tree is binary
    assert(partSoln->isPartitioningTreeBinary() == true);

    /////////////////////////////////////////////////////////////////
    // Get partitioning tree from solution
    /////////////////////////////////////////////////////////////////
    part_t numTreeVerts = 0;
    std::vector<part_t> permPartNums; // peritab in Scotch

    std::vector<part_t> splitRangeBeg;
    std::vector<part_t> splitRangeEnd;
    std::vector<part_t> treeVertParents;

    partSoln->getPartitionTree(numTreeVerts, permPartNums, splitRangeBeg,
                               splitRangeEnd, treeVertParents);
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    // Grab part numbers that are to be split by the separators
    //
    // Each separator i is represented by 1 integer and two sets of part_t's in the
    // the following 3 structure:
    //
    //             sepLevels[i] - level of separator
    //             sepParts1[i] - 1st set of parts on one side of separator
    //             sepParts2[i] - 2nd set of parts on other side of separator
    //
    /////////////////////////////////////////////////////////////////
    std::vector<part_t> levInds;

    std::vector<part_t> sepLevels;
    std::vector<std::set<part_t>> sepParts1;
    std::vector<std::set<part_t>> sepParts2;

    std::vector<std::pair<part_t, part_t>> treeIndxToSepLev(treeVertParents.size());

    part_t maxLevel = 0;

    //View of separator tree structure of solution
    part_t *sepTreeView = solution_->getSeparatorTreeView();

    part_t sepTreeIndx = treeVertParents.size() - 1;

    sepTreeView[sepTreeIndx] = -1;

    buildPartTree(0, levInds, numTreeVerts, permPartNums, splitRangeBeg, splitRangeEnd, treeVertParents,
                  sepLevels, sepParts1, sepParts2, maxLevel, sepTreeIndx, sepTreeView, treeIndxToSepLev);

    solution_->setHaveSeparatorTree(true);

    // std::cout << "sepTreeView: ";
    // for(part_t i=0; i<treeVertParents.size(); i++)
    // {
    //   std::cout << sepTreeView[i] << " ";
    // }
    // std::cout << std::endl;

    // std::cout << "treeIndxToSepLev:" << std::endl;
    // for(part_t i=0; i<treeVertParents.size(); i++)
    // {
    //   std::cout << treeIndxToSepLev[i].first << " " << treeIndxToSepLev[i].second <<std::endl;
    // }

    part_t numSeparators = sepLevels.size();
    /////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////
    // Create a map that maps each part number to a new number based on
    // the level of the hiearchy of the separator tree.  This allows us
    // to easily identify the boundary value vertices
    //////////////////////////////////////////////////////////////////////
    part_t numLevels = maxLevel + 1;

    std::vector<std::vector<part_t>> partLevelMap(numLevels, std::vector<part_t>(numGlobalParts));

    std::vector<part_t> sepsInLev(numLevels, 0);

    const auto graphModel = rcp(new GraphModel<Adapter>(mBaseInputAdapter, mEnv, mProblemComm, graphFlags));

    for (part_t i = 0; i < numSeparators; i++)
    {
      part_t level = sepLevels[i];

      for (typename std::set<part_t>::const_iterator iter = sepParts1[i].begin();
           iter != sepParts1[i].end(); ++iter)
      {
        partLevelMap[level][*iter] = 2 * sepsInLev[level];
      }

      for (typename std::set<part_t>::const_iterator iter = sepParts2[i].begin();
           iter != sepParts2[i].end(); ++iter)
      {
        partLevelMap[level][*iter] = 2 * sepsInLev[level] + 1;
      }

      sepsInLev[level]++;
    }
    //////////////////////////////////////////////////////////////////////

    // Set of separator vertices.  Used to keep track of what vertices are
    // already in previous calculated separators.  These vertices should be
    // excluded from future separator calculations
    std::set<lno_t> sepVerts;
    std::vector<std::vector<std::set<lno_t>>> sepVertsByLevel(numLevels);

    //////////////////////////////////////////////////////////////////////
    // Loop over each cut
    //    1. Build boundary layer between parts
    //    2. Build vertex separator from boundary layer
    //////////////////////////////////////////////////////////////////////
    for (part_t level = 0; level < numLevels; level++)
    {
      sepVertsByLevel[level].resize(sepsInLev[level]);

      for (part_t levIndx = 0; levIndx < sepsInLev[level]; levIndx++)
      {

        ///////////////////////////////////////////////////////////////
        // Build boundary layer between parts (edge separator)
        ///////////////////////////////////////////////////////////////
        lno_t bigraphNumU = 0, bigraphNumV = 0;
        lno_t bigraphNumE = 0; // Should probably be size_t, but making lno_t for Matcher
        std::vector<lno_t> bigraphVMapU;
        std::vector<lno_t> bigraphVMapV;

        std::vector<lno_t> bigraphCRSRowPtr;
        std::vector<lno_t> bigraphCRSCols;

        getBoundLayer(levIndx, partLevelMap[level], parts, sepVerts,
                      bigraphNumU, bigraphNumV, bigraphNumE,
                      bigraphCRSRowPtr, bigraphCRSCols,
                      bigraphVMapU, bigraphVMapV, graphModel);

        // std::cout << "Bipartite graph: " << bigraphNumU << " " << bigraphNumV << " "
        // 	  << bigraphNumE << std::endl;

        // for (size_t i=0;i<bigraphVMapU.size();i++)
        // {
        //   std::cout << "boundVertU: " << bigraphVMapU[i] << std::endl;
        // }

        // for (size_t i=0;i<bigraphVMapV.size();i++)
        // {
        //   std::cout << "boundVertV: " << bigraphVMapV[i] << std::endl;
        // }

        // for (lno_t rownum=0;rownum<bigraphNumU;rownum++)
        // {

        //    for (lno_t eIdx=bigraphCRSRowPtr[rownum];eIdx<bigraphCRSRowPtr[rownum+1];eIdx++)
        //    {
        //       std::cout << "bipartite E: " << bigraphVMapU[rownum] << ", " << bigraphVMapV[ bigraphCRSCols[eIdx]]
        // 		<< " ( "  << rownum << "," << bigraphCRSCols[eIdx] << " )" << std::endl;
        //    }

        // }
        ///////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////
        // Calculate bipartite matching from boundary layer
        ///////////////////////////////////////////////////////////////
        if (bigraphNumU > 0)
        {
          assert(bigraphNumV > 0);

          Matcher<lno_t> bpMatch(bigraphCRSRowPtr.data(), bigraphCRSCols.data(), bigraphNumU, bigraphNumV, bigraphNumE);
          bpMatch.match();

          const std::vector<lno_t> &vertUMatches = bpMatch.getVertexUMatches();
          const std::vector<lno_t> &vertVMatches = bpMatch.getVertexVMatches();
          ///////////////////////////////////////////////////////////////

          ///////////////////////////////////////////////////////////////
          // Calculate vertex cover (which is vertex separator) from matching
          ///////////////////////////////////////////////////////////////
          std::vector<lno_t> VC;

          bpMatch.getVCfromMatching(bigraphCRSRowPtr, bigraphCRSCols, vertUMatches, vertVMatches,
                                    bigraphVMapU, bigraphVMapV, VC);

          for (size_t i = 0; i < VC.size(); i++)
          {
            sepVerts.insert(VC[i]);

            sepVertsByLevel[level][levIndx].insert(VC[i]);
            //	    std::cout << "VC: " << VC[i] << std::endl;
          }
          ///////////////////////////////////////////////////////////////
        }
      }
    }
    //////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////
    // Fill solution structures: invperm and separatorRange
    //////////////////////////////////////////////////////////////////////
    bool inverse = true;
    lno_t *ipermView = solution_->getPermutationView(inverse);
    lno_t *sepRangeView = solution_->getSeparatorRangeView();

    fillSolutionIperm(graphModel, parts, sepVerts, sepVertsByLevel, treeIndxToSepLev, ipermView, sepRangeView);

    solution_->setHaveInverse(true);
    solution_->setHaveSeparatorRange(true);
    //////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////
    // Output separators
    //////////////////////////////////////////////////////////////////////
    std::cout << "Separators: " << std::endl;

    part_t nLevels = sepVertsByLevel.size();
    for (part_t level = 0; level < nLevels; level++)
    {
      //sepVertsByLevel[level].resize(sepsInLev[level]);
      part_t nSepsOnLev = sepVertsByLevel[level].size();

      for (part_t levIndx = 0; levIndx < nSepsOnLev; levIndx++)
      {
        std::cout << "  Separator " << level << " " << levIndx << ": ";

        typename std::set<lno_t>::const_iterator iterS;
        for (iterS = sepVertsByLevel[level][levIndx].begin(); iterS != sepVertsByLevel[level][levIndx].end(); ++iterS)
        {
          std::cout << *iterS << " ";
        }
        std::cout << std::endl;
      }
    }
    //////////////////////////////////////////////////////////////////////

    // std::cout << "iPerm: ";
    // for(part_t i=0; i<numVerts; i++)
    // {
    //   std::cout << ipermView[i] << " ";
    // }
    // std::cout << std::endl;

    // std::cout << "sepRange: ";
    // for(part_t i=0; i<treeVertParents.size()+1; i++)
    // {
    //   std::cout << sepRangeView[i] << " ";
    // }
    // std::cout << std::endl;

    //////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////

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
  void AlgND<Adapter>::getBoundLayer(part_t levelIndx, const std::vector<part_t> &partMap,
                                     const part_t *parts,
                                     const std::set<lno_t> &excVerts,
                                     lno_t &bigraphNumS, lno_t &bigraphNumT, lno_t &bigraphNumE,
                                     std::vector<lno_t> &bigraphCRSRowPtr, std::vector<lno_t> &bigraphCRSCols,
                                     std::vector<lno_t> &bigraphVMapS, std::vector<lno_t> &bigraphVMapT,
                                     const RCP<const GraphModel<Adapter>> &graphModel)
  {
    typedef typename Adapter::offset_t offset_t; // offset_t
    typedef StridedData<lno_t, typename Adapter::scalar_t> input_t;

    lno_t numVerts = graphModel->getLocalNumVertices();

    //Teuchos ArrayView
    // Original --  ArrayView< const lno_t > eIDs;
    ArrayView<const gno_t> eIDs;
    ArrayView<const offset_t> vOffsets;
    ArrayView<input_t> wgts;

    // MMW:
    // For some reason getLocalEdgeList seems to be returning empty eIDs
    // getEdgeList expects eIDs to be an array of gno_t
    // I wanted eIDs to be lno_t since this ordering is computed on a single node and
    // it seems unnecessary to use the potentially larger gno_t.
    // The problem might be that the partitioning is being calculated on the gno_t.
    // Perhaps a solution would be set gno_t = lno_t in the partitioning.
    // For now, I'll leave this since the edgelist is unlikely to be prohibitively big

    graphModel->getEdgeList(eIDs, vOffsets, wgts);

    // original
    //  ( (GraphModel<typename Adapter::base_adapter_t>)  *mGraphModel).getEdgeList(eIDs, vOffsets, wgts);

    std::map<lno_t, std::set<lno_t>> bigraphEs;
    std::set<lno_t> vSetS;
    std::set<lno_t> vSetT;

    bigraphNumE = 0;

    for (lno_t v1 = 0; v1 < numVerts; v1++)
    {

      part_t vpart1 = partMap[parts[v1]];

      bool correctBL = (vpart1 >= 2 * levelIndx && vpart1 < 2 * (levelIndx + 1));

      // if vertex is not in the correct range of parts, it cannot be a member of
      // this boundary layer
      if (!correctBL)
      {
        continue;
      }

      // Ignore vertices that belong to set of vertices to exclude
      if (excVerts.find(v1) != excVerts.end())
      {
        continue;
      }

      //Loop over edges connected to v1
      //MMW figure out how to get this from Zoltan2
      for (offset_t j = vOffsets[v1]; j < vOffsets[v1 + 1]; j++)
      {

        lno_t v2 = eIDs[j];

        part_t vpart2 = partMap[parts[v2]];

        correctBL = (vpart2 >= 2 * levelIndx && vpart2 < 2 * (levelIndx + 1));

        // if vertex is not in the correct range of parts, it cannot be a member of
        // this boundary layer
        if (!correctBL)
        {
          continue;
        }

        // Ignore vertices that belong to set of vertices to exclude
        if (excVerts.find(v2) != excVerts.end())
        {
          continue;
        }

        if (vpart1 != vpart2)
        {
          // Vertex added to 1st set of boundary vertices
          if (vpart1 < vpart2)
          {
            vSetS.insert(v1);

            // v1, v2
            if (bigraphEs.find(v1) == bigraphEs.end())
            {
              bigraphEs[v1] = std::set<lno_t>();
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

    std::map<lno_t, lno_t> glob2LocTMap;

    lno_t indx = 0;
    for (typename std::set<lno_t>::const_iterator iter = vSetS.begin(); iter != vSetS.end(); ++iter)
    {
      bigraphVMapS[indx] = *iter;
      indx++;
    }

    bigraphVMapT.resize(bigraphNumT);
    indx = 0;
    for (typename std::set<lno_t>::const_iterator iter = vSetT.begin(); iter != vSetT.end(); ++iter)
    {
      bigraphVMapT[indx] = *iter;
      glob2LocTMap[*iter] = indx;
      indx++;
    }
    /////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////
    // Set sizes for bipartite graph data structures
    /////////////////////////////////////////////////////////////////////////
    bigraphCRSRowPtr.resize(bigraphNumS + 1);
    bigraphCRSCols.resize(bigraphNumE);
    /////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////
    // Copy bipartite graph edges into CRS format
    /////////////////////////////////////////////////////////////////////////
    bigraphCRSRowPtr[0] = 0;

    lno_t rownum = 0;
    size_t nzIndx = 0;
    typename std::map<lno_t, std::set<lno_t>>::const_iterator iterM;
    for (iterM = bigraphEs.begin(); iterM != bigraphEs.end(); ++iterM)
    {
      bigraphCRSRowPtr[rownum + 1] = bigraphCRSRowPtr[rownum] + (*iterM).second.size();

      for (typename std::set<lno_t>::const_iterator iter = (*iterM).second.begin(); iter != (*iterM).second.end(); ++iter)
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
      buildPartTree(part_t level, std::vector<part_t> &levIndx,
                    part_t startV,
                    const std::vector<part_t> &permPartNums,
                    const std::vector<part_t> &splitRangeBeg,
                    const std::vector<part_t> &splitRangeEnd,
                    const std::vector<part_t> &treeVertParents,
                    std::vector<part_t> &sepLevels,
                    std::vector<std::set<part_t>> &sepParts1, std::vector<std::set<part_t>> &sepParts2,
                    part_t &maxLev, part_t &gIndx,
                    part_t *sepTreeView, std::vector<std::pair<part_t, part_t>> &treeIndxToSepLev)
  {
    // Insert information for this separator
    maxLev = level;
    part_t tmpMaxLev = maxLev;

    //////////////////////////////////////////////////////////////////////
    // Search for indices that have parent startV
    //////////////////////////////////////////////////////////////////////
    typename std::vector<part_t>::const_iterator iter;
    std::vector<part_t> inds;
    part_t ind = 0;
    for (iter = treeVertParents.begin(); iter != treeVertParents.end(); ++iter)
    {
      if (*iter == startV)
      {
        inds.push_back(ind);
      }
      ind++;
    }
    //////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////
    // If startV has children, this will correspond to a separator.  Construct
    // appropriate data structure and then recurse
    //////////////////////////////////////////////////////////////////////
    assert(inds.size() == 2 || inds.size() == 0);

    // If startV has children
    if (inds.size() == 2)
    {

      if ((part_t)levIndx.size() < level + 1)
      {
        levIndx.push_back(0);
      }
      else
      {
        levIndx[level]++;
      }

      // std::cout << "gIndx " << gIndx << ": separator " << level << " " << levIndx[level] << std::endl;

      treeIndxToSepLev[gIndx].first = level;
      treeIndxToSepLev[gIndx].second = levIndx[level];

      part_t v0 = inds[0];
      part_t v1 = inds[1];

      sepLevels.push_back(level);

      sepParts1.push_back(std::set<part_t>());
      typename std::vector<std::set<part_t>>::iterator setIter = sepParts1.end();
      setIter--; // set iterator to point to new set

      for (part_t k = splitRangeBeg[v0]; k < splitRangeEnd[v0]; k++)
      {
        (*setIter).insert(permPartNums[k]);
      }

      sepParts2.push_back(std::set<part_t>());
      setIter = sepParts2.end();
      setIter--; // set iterator to point to new set

      for (part_t k = splitRangeBeg[v1]; k < splitRangeEnd[v1]; k++)
      {
        (*setIter).insert(permPartNums[k]);
      }

      part_t parentNode = gIndx;
      gIndx--;
      sepTreeView[gIndx] = parentNode;

      // Recursively call function on children
      buildPartTree(level + 1, levIndx, v0,
                    permPartNums, splitRangeBeg, splitRangeEnd, treeVertParents,
                    sepLevels, sepParts1, sepParts2,
                    tmpMaxLev,
                    gIndx, sepTreeView, treeIndxToSepLev);
      if (tmpMaxLev > maxLev)
      {
        maxLev = tmpMaxLev;
      }

      gIndx--;
      sepTreeView[gIndx] = parentNode;

      buildPartTree(level + 1, levIndx, v1,
                    permPartNums, splitRangeBeg, splitRangeEnd, treeVertParents,
                    sepLevels, sepParts1, sepParts2,
                    tmpMaxLev,
                    gIndx, sepTreeView, treeIndxToSepLev);
      if (tmpMaxLev > maxLev)
      {
        maxLev = tmpMaxLev;
      }
    }
    else // no children, so this is not a separator
    {
      // std::cout << "gIndx " << gIndx << " leaf: " << permPartNums[splitRangeBeg[startV]] << std::endl;
      treeIndxToSepLev[gIndx].first = -1;
      treeIndxToSepLev[gIndx].second = permPartNums[splitRangeBeg[startV]];

      maxLev--;
    }
    //////////////////////////////////////////////////////////////////////
  }

  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  template <typename Adapter>
  void AlgND<Adapter>::
      fillSolutionIperm(const RCP<const GraphModel<Adapter>> &graphModel,
                        const part_t *parts, const std::set<lno_t> &sepVerts,
                        const std::vector<std::vector<std::set<lno_t>>> &sepVertsByLevel,
                        const std::vector<std::pair<part_t, part_t>> &treeIndxToSepLev,
                        lno_t *ipermView, lno_t *sepRangeView)
  {
    lno_t permIndx = 0;
    lno_t sepRangeIndx = 0;
    sepRangeView[sepRangeIndx] = 0;

    for (size_t i = 0; i < treeIndxToSepLev.size(); i++)
    {
      part_t lev = treeIndxToSepLev[i].first;
      ////////////////////////////////////////////////////////////////////
      // Leaf node of separator tree
      ////////////////////////////////////////////////////////////////////
      if (lev == -1)
      {
        std::set<lno_t> idSet;
        getIdsOfPart(graphModel, parts, treeIndxToSepLev[i].second, sepVerts, idSet);

        for (typename std::set<lno_t>::const_iterator setIter = idSet.begin(); setIter != idSet.end();
             ++setIter)
        {
          ipermView[permIndx] = *setIter;
          permIndx++;
        }
        sepRangeView[sepRangeIndx + 1] = sepRangeView[sepRangeIndx] + idSet.size();
        sepRangeIndx++;
      }
      ////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////
      // Internal "separator node" of separator tree
      ////////////////////////////////////////////////////////////////////
      else
      {
        const std::set<lno_t> &idSet = sepVertsByLevel[lev][treeIndxToSepLev[i].second];

        for (typename std::set<lno_t>::const_iterator setIter = idSet.begin();
             setIter != idSet.end(); ++setIter)
        {
          ipermView[permIndx] = *setIter;
          permIndx++;
        }
        sepRangeView[sepRangeIndx + 1] = sepRangeView[sepRangeIndx] + idSet.size();
        sepRangeIndx++;
      }
      ////////////////////////////////////////////////////////////////////
    }
  }
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  template <typename Adapter>
  void AlgND<Adapter>::
      getIdsOfPart(const RCP<const GraphModel<Adapter>> &graphModel, const part_t *parts,
                   part_t targetPart, const std::set<lno_t> &idsToExcl, std::set<lno_t> &outIds)
  {
    size_t numVerts = graphModel->getLocalNumVertices();

    for (size_t i = 0; i < numVerts; i++)
    {
      if (parts[i] == targetPart && idsToExcl.find(i) == idsToExcl.end())
      {
        outIds.insert(i);
      }
    }
  }
  //////////////////////////////////////////////////////////////////////////////

} // namespace Zoltan2

#endif

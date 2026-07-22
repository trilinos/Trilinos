// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_ALGSERIALGREEDY_HPP_
#define _ZOLTAN2_ALGSERIALGREEDY_HPP_

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_ColoringSolution.hpp>

////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgSerialGreedy.hpp
//! \brief Serial greedy first-fit graph coloring (local graph only)

namespace Zoltan2{

template <typename Adapter>
class AlgSerialGreedy : public Algorithm<Adapter>
{
  private:
    typedef typename Adapter::lno_t lno_t;
    typedef typename Adapter::gno_t gno_t;
    typedef typename Adapter::offset_t offset_t;
    typedef typename Adapter::scalar_t scalar_t;
    // Class member variables
    const RCP<const typename Adapter::base_adapter_t> adapter_;
    RCP<Teuchos::ParameterList> pl_;
    RCP<Environment> env_;
    RCP<const Teuchos::Comm<int> > comm_;
    modelFlag_t graphFlags_;

  public:
  AlgSerialGreedy(
    const RCP<const typename Adapter::base_adapter_t> &adapter,
    const RCP<Teuchos::ParameterList> &pl,
    const RCP<Environment> &env,
    const RCP<const Teuchos::Comm<int> > &comm,
    modelFlag_t graphFlags
  ) : adapter_(adapter), pl_(pl), env_(env), comm_(comm), graphFlags_(graphFlags)
  {
  }

  // Main entry point for graph coloring.
  void color(
    const RCP<ColoringSolution<Adapter> > &solution
  )
  {
    HELLO;

    // Color local graph. Global coloring is supported in Zoltan (not Zoltan2).
    // Get local graph.
    ArrayView<const gno_t> edgeIds;
    ArrayView<const offset_t> offsets;
    ArrayView<StridedData<lno_t, scalar_t> > wgts; // Not used; needed by getLocalEdgeList

    const auto model = rcp(new GraphModel<typename Adapter::base_adapter_t>(adapter_, env_, comm_, graphFlags_));
    const size_t nVtx = model->getLocalNumVertices(); // Assume (0,nvtx-1)
    model->getEdgeList(edgeIds, offsets, wgts); // Don't need wgts

#if 0
    // Debug
    cout << "Debug: Local graph from getLocalEdgeList" << endl;
    cout << "rank " << comm_->getRank() << ": nVtx= " << nVtx << endl;
    cout << "rank " << comm_->getRank() << ": edgeIds: " << edgeIds << endl;
    cout << "rank " << comm_->getRank() << ": offsets: " << offsets << endl;
#endif

    // Get color array to fill.
    // TODO: Allow user to input an old coloring.
    ArrayRCP<int> colors = solution->getColorsRCP();
    for (size_t i=0; i<nVtx; i++){
      colors[i] = 0;
    }

    // Let colorCrsGraph do the real work.
    env_->timerStart(MACRO_TIMERS, "Coloring algorithm");
    colorCrsGraph(nVtx, edgeIds, offsets, colors);
    env_->timerStop(MACRO_TIMERS, "Coloring algorithm");
    return;
  }

  // Color graph given by two arrays. API may change. Expert users only!
  void colorCrsGraph(
    const size_t nVtx,
    ArrayView<const gno_t> edgeIds,
    ArrayView<const offset_t> offsets,
    ArrayRCP<int> colors
  )
  {
    HELLO;

    // Find max degree, since (max degree)+1 is an upper bound.
    offset_t maxDegree = 0;
    for (size_t i=0; i<nVtx; i++){
      if (offsets[i+1]-offsets[i] > maxDegree)
        maxDegree = offsets[i+1]-offsets[i];
    }

    // Greedy coloring.
    // Use natural order for now.
    // TODO: Support better orderings (e.g., Smallest-Last)
    int maxColor = 0;

    // array of size #colors: forbidden[i]=v means color[v]=i so i is forbidden
    Teuchos::Array<int> forbidden(maxDegree+2, 0);

    // LeastUsed: need array of size #colors
    Teuchos::Array<lno_t> numVerticesWithColor(maxDegree+2, 0);

    // Get colorChoice from parameter list.
    Teuchos::ParameterList &pl = env_->getParametersNonConst();
    std::string colorChoice = pl.get<std::string>("color_choice", "FirstFit");

    for (size_t i=0; i<nVtx; i++){
      //std::cout << "Debug: i= " << i << std::endl;
      lno_t v=i; // TODO: Use ordering here.
      for (offset_t j=offsets[v]; j<offsets[v+1]; j++){
        gno_t nbor = edgeIds[j];
        //std::cout << "Debug: nbor= " << nbor << ", color= " << colors[nbor] << std::endl;
        if (colors[nbor] > 0){
          // Neighbors' colors are forbidden
          forbidden[colors[nbor]] = v;
        }
      }

      // Pick color for v

      // Keep colors[v] if possible, otherwise find valid color.
      if ((colors[v]==0) || ((colors[v]>0) && forbidden[colors[v]] == v)){

        if (colorChoice.compare("FirstFit")){
          // Pick first (smallest) available color > 0
          for (int c=1; c <= maxColor+1; c++){
            if (forbidden[c] != v){
              colors[v] = c;
              break;
            }
          }
        }
        else if (colorChoice.compare("Random")){
          // Pick random available color.
          // Truely random is slow, please consider RandomFast instead.
          int numAvail = 0;
          Teuchos::Array<int> avail(maxColor+1);
          for (int c=1; c < maxColor+1; c++){
            if (forbidden[c] != v){
              avail[numAvail++] = c;
            }
          }
          if (numAvail==0)
            colors[v] = maxColor+1;
          else
            colors[v] = avail[rand()%numAvail];
        }
        else if (colorChoice.compare("RandomFast")){
          // Pick random color, then find first available color after that.
          bool foundColor = false;
          int r = (rand() % maxColor) +1;
          for (int c=r; c <= maxColor; c++){
            if (forbidden[c] != v){
              colors[v] = c;
              foundColor = true;
              break;
            }
          }
          if (!foundColor){ // Look for colors in [1, r)
            for (int c=1; c < r; c++){
              if (forbidden[c] != v){
                colors[v] = c;
                foundColor = true;
                break;
              }
            }
          }
          if (!foundColor) colors[v] = maxColor+1;
        }
        else if (colorChoice.compare("LeastUsed")){
          // Pick least used available color.
          // Simple linear algorithm; could maintain a priority queue but not sure any faster?
          int leastUsedColor = 1;
          for (int c=1; c <= maxColor; c++){
            if (forbidden[c] != v){
              if (numVerticesWithColor[c] < leastUsedColor){
                leastUsedColor = c;
              }
            }
          }
          colors[v] = leastUsedColor;

          // Update color counts
          numVerticesWithColor[colors[v]]++;
        }

        if ((v==0) && colors[v]==0) colors[v]=1; // Corner case for first vertex

        // If we used a new color, increase maxColor.
        if (colors[v] > maxColor){
          maxColor = colors[v];
        }
      }
    }

    return;
  }

  /*! \brief Set up validators specific to this algorithm
  */
  static void getValidParameters(ParameterList & pl)
  {
    RCP<Teuchos::StringValidator> color_choice_Validator = Teuchos::rcp(
      new Teuchos::StringValidator(
      Teuchos::tuple<std::string>(
        "FirstFit", "Random", "RandomFast", "LeastUsed" )));
    pl.set("color_choice", "FirstFit", "selection criterion for coloring",
      color_choice_Validator);
  }
};
}
#endif

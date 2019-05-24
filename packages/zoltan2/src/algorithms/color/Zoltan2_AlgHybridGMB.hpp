#ifndef _ZOLTAN2_ALGHYBRIDGMB_HPP_
#define _ZOLTAN2_ALGHYBRIDGMB_HPP_

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_ColoringSolution.hpp>

//////////////////////////////////////////////
//! \file Zoltan2_AlgHybridGMB.hpp
//! \brief A hybrid version of the framework proposed by Gebremedhin, Manne, and Boman

namespace Zoltan2{

template <typename Adapter>

class AlgHybridGMB : public Algorithm<Adapter>
{
  private:
    typedef typename Adapter::lno_t lno_t;
    typedef typename Adapter::gno_t gno_t;
    typedef typename Adapter::offset_t offset_t;
    typedef typename Adapter::scalar_t scalar_t;

    RCP<GraphModel<typename Adapter::base_adapter_t> > model_;
    RCP<Teuchos::ParameterList> pl_;
    RCP<Environment> env_;
    RCP<const Teuchos::Comm<int> > comm_;

  public:
    AlgHybridGMB(
      const RCP<GraphModel<typename Adapter::base_adapter_t> > &model,
      const RCP<Teuchos::ParameterList> &pl,
      const RCP<Environment> &env,
      const RCP<const Teuchos::Comm<int> > &comm)
    : model_(model), pl_(pl), env_(env), comm_(comm) {
      //This page is intentionally blank
    }

    //Main entry point for graph coloring
    void color( const RCP<ColoringSolution<Adapter> > &solution ) {
      
      //this will color the global graph in a manner similar to Zoltan
      //once everything is implemented, that is.
      ArrayView<const gno_t> edgeIds;
      ArrayView<const offset_t> offsets;
      ArrayView<StridedData<lno_t, scalar_t> > wgts;
      
      const size_t nVtx = model_ -> getLocalNumVertices();
      model_->getEdgeList(edgeIds, offsets, wgts); //wgts is not used

      //Get color array to fill
      ArrayRCP<int> colors = solution->getColorsRCP();
      for(size_t i=0; i<nVtx; i++){
        colors[i] = 0;
      } 

      //create maps for FEMultiVector

      // call actual coloring function
      // THESE ARGUMENTS WILL NEED TO CHANGE,
      // THESE ARE A COPY OF THE EXISTING ALGORITHM CLASS.
      hybridGMB(nVtx, edgeIds, offsets,colors);
    }
    
    void hybridGMB(const size_t nVtx, ArrayView<const gno_t> edgeIds, 
                   ArrayView<const offset_t> offsets, ArrayRCP<int> colors){
      //color the interior vertices (maybe)

      //color boundary vertices using FEMultiVector (definitely)

      //color interior vertices if not colored yet.
      
    }
}

}

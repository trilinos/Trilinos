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

/*! \file Zoltan2_TpetraRowGraphAdapter.hpp
    \brief Defines TpetraRowGraphAdapter class.
*/

#ifndef _ZOLTAN2_TPETRAROWGRAPHADAPTER_HPP_
#define _ZOLTAN2_TPETRAROWGRAPHADAPTER_HPP_

#include <Zoltan2_GraphAdapter.hpp>
#include <Zoltan2_StridedData.hpp>
#include <Zoltan2_PartitioningHelpers.hpp>
#include <Tpetra_RowGraph.hpp>

namespace Zoltan2 {

/*!  \brief Provides access for Zoltan2 to Tpetra::RowGraph data.

    \todo test for memory alloc failure when we resize a vector
    \todo we assume FillComplete has been called.  Should we support
                objects that are not FillCompleted.

    The template parameter is the user's input object:
     \li Tpetra::CrsGraph
     \li Tpetra::RowGraph
     \li Epetra_CrsGraph

    The \c scalar_t type, representing use data such as matrix values, is
    used by Zoltan2 for weights, coordinates, part sizes and
    quality metrics.
    Some User types (like Tpetra::CrsMatrix) have an inherent scalar type,
    and some
    (like Tpetra::CrsGraph) do not.  For such objects, the scalar type is
    set by Zoltan2 to \c float.  If you wish to change it to double, set
    the second template parameter to \c double.
*/

template <typename User, typename UserCoord=User>
  class TpetraRowGraphAdapter : public GraphAdapter<User,UserCoord> {

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::offset_t    offset_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::part_t   part_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef User user_t;
  typedef UserCoord userCoord_t;
#endif

  /*! \brief Destructor
   */
  ~TpetraRowGraphAdapter() { }

  /*! \brief Constructor for graph with no weights or coordinates.
   *  \param ingraph the  Tpetra::RowGraph
   *  \param numVtxWeights  the number of weights per vertex (default = 0)
   *  \param numEdgeWeights the number of weights per edge  (default = 0)
   *
   * Most adapters do not have RCPs in their interface.  This
   * one does because the user is obviously a Trilinos user.
   */

  TpetraRowGraphAdapter(const RCP<const User> &ingraph,
                        int nVtxWeights=0, int nEdgeWeights=0);

  /*! \brief Provide a pointer to weights for the primary entity type.
   *    \param val A pointer to the weights for index \c idx.
   *    \param stride    A stride for the \c val array.  If \stride is
   *             \c k, then val[n * k] is the weight for the
   *             \c n th entity for index \idx.
   *    \param idx A number from 0 to one less than
   *          weight idx specified in the constructor.
   *
   *  The order of the weights should match the order that
   *  entities appear in the input data structure.
   */

  void setWeights(const scalar_t *val, int stride, int idx);

  /*! \brief Provide a pointer to vertex weights.
   *    \param val A pointer to the weights for index \c idx.
   *    \param stride    A stride for the \c val array.  If \stride is
   *             \c k, then val[n * k] is the weight for the
   *             \c n th vertex for index \idx.
   *    \param idx A number from 0 to one less than
   *          number of vertex weights specified in the constructor.
   *
   *  The order of the vertex weights should match the order that
   *  vertices appear in the input data structure.
   *     \code
   *       TheGraph->getRowMap()->getNodeElementList()
   *     \endcode
   */

  void setVertexWeights(const scalar_t *val, int stride, int idx);

  /*! \brief Specify an index for which the weight should be
              the degree of the entity
   *    \param idx Zoltan2 will use the entity's
   *         degree as the entity weight for index \c idx.
   */
  void setWeightIsDegree(int idx);

  /*! \brief Specify an index for which the vertex weight should be
              the degree of the vertex
   *    \param idx Zoltan2 will use the vertex's
   *         degree as the vertex weight for index \c idx.
   */
  void setVertexWeightIsDegree(int idx);

  /*! \brief Provide a pointer to edge weights.
   *    \param val A pointer to the weights for index \c idx.
   *    \param stride    A stride for the \c val array.  If \stride is
   *             \c k, then val[n * k] is the weight for the
   *             \c n th edge for index \idx.
   *    \param dim A number from 0 to one less than the number
   *          of edge weights specified in the constructor.
   *
   *  The order of the edge weights should follow the order that the
   *  the vertices and edges appear in the input data structure.
   *
   *  By vertex:
   *     \code
   *       TheGraph->getRowMap()->getNodeElementList()
   *     \endcode
   *
   *  Then by vertex neighbor:
   *     \code
   *       TheGraph->getLocalRowView(vertexNum, neighborList);
   *     \endcode
   */

  void setEdgeWeights(const scalar_t *val, int stride, int idx);

  ////////////////////////////////////////////////////
  // The Adapter interface.
  ////////////////////////////////////////////////////

  ////////////////////////////////////////////////////
  // The GraphAdapter interface.
  ////////////////////////////////////////////////////

  // TODO:  Assuming rows == objects;
  // TODO:  Need to add option for columns or nonzeros?
  size_t getLocalNumVertices() const { return graph_->getNodeNumRows(); }

  void getVertexIDsView(const gno_t *&ids) const
  {
    ids = NULL;
    if (getLocalNumVertices())
      ids = graph_->getRowMap()->getNodeElementList().getRawPtr();
  }

  size_t getLocalNumEdges() const { return graph_->getNodeNumEntries(); }

  void getEdgesView(const offset_t *&offsets, const gno_t *&adjIds) const
  {
    offsets = offs_.getRawPtr();
    adjIds = (getLocalNumEdges() ? adjids_.getRawPtr() : NULL);
  }

  int getNumWeightsPerVertex() const { return nWeightsPerVertex_;}

  void getVertexWeightsView(const scalar_t *&weights, int &stride,
                            int idx) const
  {
    if(idx<0 || idx >= nWeightsPerVertex_)
    {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid vertex weight index " << idx << std::endl;
      throw std::runtime_error(emsg.str());
    }

    size_t length;
    vertexWeights_[idx].getStridedList(length, weights, stride);
  }

  bool useDegreeAsVertexWeight(int idx) const {return vertexDegreeWeight_[idx];}

  int getNumWeightsPerEdge() const { return nWeightsPerEdge_;}

  void getEdgeWeightsView(const scalar_t *&weights, int &stride, int idx) const
  {
    if(idx<0 || idx >= nWeightsPerEdge_)
    {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid edge weight index " << idx << std::endl;
      throw std::runtime_error(emsg.str());
    }

    size_t length;
    edgeWeights_[idx].getStridedList(length, weights, stride);
  }


  template <typename Adapter>
    void applyPartitioningSolution(const User &in, User *&out,
      const PartitioningSolution<Adapter> &solution) const;

  template <typename Adapter>
    void applyPartitioningSolution(const User &in, RCP<User> &out,
      const PartitioningSolution<Adapter> &solution) const;

private:

  RCP<const User> graph_;

  ArrayRCP<const offset_t> offs_;
  ArrayRCP<const gno_t> adjids_;

  int nWeightsPerVertex_;
  ArrayRCP<StridedData<lno_t, scalar_t> > vertexWeights_;
  ArrayRCP<bool> vertexDegreeWeight_;

  int nWeightsPerEdge_;
  ArrayRCP<StridedData<lno_t, scalar_t> > edgeWeights_;

  int coordinateDim_;
  ArrayRCP<StridedData<lno_t, scalar_t> > coords_;

  RCP<User> doMigration(const User &from, size_t numLocalRows,
                        const gno_t *myNewRows) const;
};

/////////////////////////////////////////////////////////////////
// Definitions
/////////////////////////////////////////////////////////////////

template <typename User, typename UserCoord>
  TpetraRowGraphAdapter<User,UserCoord>::TpetraRowGraphAdapter(
    const RCP<const User> &ingraph, int nVtxWgts, int nEdgeWgts):
      graph_(ingraph), offs_(),
      adjids_(),
      nWeightsPerVertex_(nVtxWgts), vertexWeights_(), vertexDegreeWeight_(),
      nWeightsPerEdge_(nEdgeWgts), edgeWeights_(),
      coordinateDim_(0), coords_()
{
  typedef StridedData<lno_t,scalar_t> input_t;

  size_t nvtx = graph_->getNodeNumRows();
  size_t nedges = graph_->getNodeNumEntries();
  size_t maxnumentries =
         graph_->getNodeMaxNumRowEntries(); // Diff from CrsMatrix

  // Unfortunately we have to copy the offsets and edge Ids
  // because edge Ids are not usually stored in vertex id order.

  size_t n = nvtx + 1;
  offset_t *offs = new offset_t [n];

  if (!offs)
  {
    std::cerr << "Error: " << __FILE__ << ", " << __LINE__<< std::endl;
    std::cerr << n << " objects" << std::endl;
    throw std::bad_alloc();
  }

  gno_t *adjids = NULL;
  if (nedges)
  {
    adjids = new gno_t [nedges];

    if (!adjids)
    {
      std::cerr << "Error: " << __FILE__ << ", " << __LINE__<< std::endl;
      std::cerr << nedges << " objects" << std::endl;
      throw std::bad_alloc();
    }
  }

  ArrayRCP<lno_t> nbors(maxnumentries); // Diff from CrsGraph

  offs[0] = 0;
  for (size_t v=0; v < nvtx; v++){
    graph_->getLocalRowCopy(v, nbors(), nedges);  // Diff from CrsGraph
    offs[v+1] = offs[v] + nedges;
    for (offset_t e=offs[v], i=0; e < offs[v+1]; e++) {
        adjids[e] = graph_->getColMap()->getGlobalElement(nbors[i++]);
    }
  }

  offs_ = arcp(offs, 0, n, true);
  adjids_ = arcp(adjids, 0, nedges, true);

  if (nWeightsPerVertex_ > 0) {
    vertexWeights_ =
          arcp(new input_t[nWeightsPerVertex_], 0, nWeightsPerVertex_, true);
    vertexDegreeWeight_ =
          arcp(new bool[nWeightsPerVertex_], 0, nWeightsPerVertex_, true);
    for (int i=0; i < nWeightsPerVertex_; i++)
      vertexDegreeWeight_[i] = false;
  }

  if (nWeightsPerEdge_ > 0)
    edgeWeights_ = arcp(new input_t[nWeightsPerEdge_], 0, nWeightsPerEdge_, true);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  void TpetraRowGraphAdapter<User,UserCoord>::setWeights(
    const scalar_t *weightVal, int stride, int idx)
{
  if (this->getPrimaryEntityType() == GRAPH_VERTEX)
    setVertexWeights(weightVal, stride, idx);
  else
    setEdgeWeights(weightVal, stride, idx);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  void TpetraRowGraphAdapter<User,UserCoord>::setVertexWeights(
    const scalar_t *weightVal, int stride, int idx)
{
  typedef StridedData<lno_t,scalar_t> input_t;
  if(idx<0 || idx >= nWeightsPerVertex_)
  {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid vertex weight index " << idx << std::endl;
      throw std::runtime_error(emsg.str());
  }

  size_t nvtx = getLocalNumVertices();
  ArrayRCP<const scalar_t> weightV(weightVal, 0, nvtx*stride, false);
  vertexWeights_[idx] = input_t(weightV, stride);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  void TpetraRowGraphAdapter<User,UserCoord>::setWeightIsDegree(
    int idx)
{
  if (this->getPrimaryEntityType() == GRAPH_VERTEX)
    setVertexWeightIsDegree(idx);
  else {
    std::ostringstream emsg;
    emsg << __FILE__ << "," << __LINE__
         << " error:  setWeightIsNumberOfNonZeros is supported only for"
         << " vertices" << std::endl;
    throw std::runtime_error(emsg.str());
  }
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  void TpetraRowGraphAdapter<User,UserCoord>::setVertexWeightIsDegree(
    int idx)
{
  if(idx<0 || idx >= nWeightsPerVertex_)
  {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid vertex weight index " << idx << std::endl;
      throw std::runtime_error(emsg.str());
  }

  vertexDegreeWeight_[idx] = true;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  void TpetraRowGraphAdapter<User,UserCoord>::setEdgeWeights(
    const scalar_t *weightVal, int stride, int idx)
{
  typedef StridedData<lno_t,scalar_t> input_t;

  if(idx<0 || idx >= nWeightsPerEdge_)
  {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid edge weight index " << idx << std::endl;
      throw std::runtime_error(emsg.str());
  }

  size_t nedges = getLocalNumEdges();
  ArrayRCP<const scalar_t> weightV(weightVal, 0, nedges*stride, false);
  edgeWeights_[idx] = input_t(weightV, stride);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  template<typename Adapter>
    void TpetraRowGraphAdapter<User,UserCoord>::applyPartitioningSolution(
      const User &in, User *&out,
      const PartitioningSolution<Adapter> &solution) const
{
  // Get an import list (rows to be received)
  size_t numNewVtx;
  ArrayRCP<gno_t> importList;
  try{
    numNewVtx = Zoltan2::getImportList<Adapter,
                                       TpetraRowGraphAdapter<User,UserCoord> >
                                      (solution, this, importList);
  }
  Z2_FORWARD_EXCEPTIONS;

  // Move the rows, creating a new graph.
  RCP<User> outPtr = doMigration(in, numNewVtx, importList.getRawPtr());
  out = outPtr.get();
  outPtr.release();
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  template<typename Adapter>
    void TpetraRowGraphAdapter<User,UserCoord>::applyPartitioningSolution(
      const User &in, RCP<User> &out,
      const PartitioningSolution<Adapter> &solution) const
{
  // Get an import list (rows to be received)
  size_t numNewVtx;
  ArrayRCP<gno_t> importList;
  try{
    numNewVtx = Zoltan2::getImportList<Adapter,
                                       TpetraRowGraphAdapter<User,UserCoord> >
                                      (solution, this, importList);
  }
  Z2_FORWARD_EXCEPTIONS;

  // Move the rows, creating a new graph.
  out = doMigration(in, numNewVtx, importList.getRawPtr());
}


////////////////////////////////////////////////////////////////////////////
template < typename User, typename UserCoord>
RCP<User> TpetraRowGraphAdapter<User,UserCoord>::doMigration(
  const User &from,
  size_t numLocalRows,
  const gno_t *myNewRows
) const
{
  typedef Tpetra::Map<lno_t, gno_t, node_t> map_t;
  typedef Tpetra::CrsGraph<lno_t, gno_t, node_t> tcrsgraph_t;

  // We cannot create a Tpetra::RowGraph, unless the underlying type is
  // something we know (like Tpetra::CrsGraph).
  // If the underlying type is something different, the user probably doesn't
  // want a Tpetra::CrsGraph back, so we throw an error.

  // Try to cast "from" graph to a TPetra::CrsGraph
  // If that fails we throw an error.
  // We could cast as a ref which will throw std::bad_cast but with ptr
  // approach it might be clearer what's going on here
  const tcrsgraph_t *pCrsGraphSrc = dynamic_cast<const tcrsgraph_t *>(&from);

  if(!pCrsGraphSrc) {
    throw std::logic_error("TpetraRowGraphAdapter cannot migrate data for "
                           "your RowGraph; it can migrate data only for "
                           "Tpetra::CrsGraph.  "
                           "You can inherit from TpetraRowGraphAdapter and "
                           "implement migration for your RowGraph.");
  }

  // source map
  const RCP<const map_t> &smap = from.getRowMap();
  int oldNumElts = smap->getNodeNumElements();
  gno_t numGlobalRows = smap->getGlobalNumElements();
  gno_t base = smap->getMinAllGlobalIndex();

  // target map
  ArrayView<const gno_t> rowList(myNewRows, numLocalRows);
  const RCP<const Teuchos::Comm<int> > &comm = from.getComm();
  RCP<const map_t> tmap = rcp(new map_t(numGlobalRows, rowList, base, comm));

  // importer
  Tpetra::Import<lno_t, gno_t, node_t> importer(smap, tmap);

  // number of entries in my new rows
  typedef Tpetra::Vector<gno_t, lno_t, gno_t, node_t> vector_t;
  vector_t numOld(smap);
  vector_t numNew(tmap);
  for (int lid=0; lid < oldNumElts; lid++){
    numOld.replaceGlobalValue(smap->getGlobalElement(lid),
      from.getNumEntriesInLocalRow(lid));
  }
  numNew.doImport(numOld, importer, Tpetra::INSERT);

  size_t numElts = tmap->getNodeNumElements();
  ArrayRCP<const gno_t> nnz;
  if (numElts > 0)
    nnz = numNew.getData(0);    // hangs if vector len == 0

  ArrayRCP<const size_t> nnz_size_t;

  if (numElts && sizeof(gno_t) != sizeof(size_t)){
    size_t *vals = new size_t [numElts];
    nnz_size_t = arcp(vals, 0, numElts, true);
    for (size_t i=0; i < numElts; i++){
      vals[i] = static_cast<size_t>(nnz[i]);
    }
  }
  else{
    nnz_size_t = arcp_reinterpret_cast<const size_t>(nnz);
  }

  // target graph
  RCP<tcrsgraph_t> G =
    rcp(new tcrsgraph_t(tmap, nnz_size_t(), Tpetra::StaticProfile));

  G->doImport(*pCrsGraphSrc, importer, Tpetra::INSERT);
  G->fillComplete();
  return Teuchos::rcp_dynamic_cast<User>(G);
}

}  //namespace Zoltan2

#endif

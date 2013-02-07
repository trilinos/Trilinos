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

/*! \file Zoltan2_XpetraVectorAdapter.hpp
    \brief Defines the XpetraVectorAdapter class.
*/

#ifndef _ZOLTAN2_XPETRAVECTORADAPTER_HPP_
#define _ZOLTAN2_XPETRAVECTORADAPTER_HPP_

#include <Zoltan2_VectorInput.hpp>
#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_StridedData.hpp>
#include <Zoltan2_Util.hpp>

#include <Xpetra_EpetraVector.hpp>
#include <Xpetra_TpetraVector.hpp>

namespace Zoltan2 {


/*!  \brief Provides access for Zoltan2 to an Xpetra::Vector.

    The template parameter is the user's input data type, which can be:
   \li Epetra_Vector
   \li Tpetra::Vector
   \li Xpetra::Vector


    The \c scalar_t type, representing use data such as matrix values, is
    used by Zoltan2 for weights, coordinates, part sizes and
    quality metrics.
    Some User types (like Tpetra::CrsMatrix) have an inherent scalar type,
    and some
    (like Tpetra::CrsGraph) do not.  For such objects, the scalar type is
    set by Zoltan2 to \c float.  If you wish to change it to double, set
    the second template parameter to \c double.

*/

template <typename User>
  class XpetraVectorAdapter : public VectorAdapter<User> {
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef VectorAdapter<User> base_adapter_t;
  typedef User user_t;

  typedef Xpetra::Vector<
    scalar_t, lno_t, gno_t, node_t> x_vector_t;
  typedef Xpetra::TpetraVector<
    scalar_t, lno_t, gno_t, node_t> xt_vector_t;
  typedef Xpetra::EpetraVector xe_vector_t;
#endif

  /*! \brief Destructor
   */
  ~XpetraVectorAdapter() { }

  /*! \brief Constructor   
   *
   *  \param invector  the user's Xpetra, Tpetra or Epetra Vector object
   *  \param weights  a list of pointers to arrays of weights.
   *      The number of weights per vector element is assumed to be
   *      \c weights.size().
   *  \param weightStrides  a list of strides for the \c weights.
   *     The weight for weight dimension \c n for element \c k should be
   *     found at <tt>weights[n][weightStrides[n] * k]</tt>.
   *     If \c weightStrides.size() is zero, it is assumed all strides are one
   *
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this Adapter.
   */
  XpetraVectorAdapter( const RCP<const User> &invector,
    vector<const scalar_t *> &weights, vector<int> &weightStrides);

  /*! \brief Access to the xpetra-wrapped vector
   */

  const RCP<const x_vector_t> &getVector() const
  {
    return vector_;
  }

  ////////////////////////////////////////////////////
  // The Adapter interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumberOfObjects() const { return getLocalLength();}

  int getNumberOfWeightsPerObject() const { return numWeights_;}

  size_t getObjectWeights(int dim, const scalar_t *&wgt, int &stride) const
  {
    return getVectorWeights(dim, wgt, stride);
  }

  ////////////////////////////////////////////////////
  // The VectorAdapter interface.
  ////////////////////////////////////////////////////

  int getNumberOfVectors() const { return 1; }

  int getNumberOfWeights() const {return numWeights_;}

  size_t getLocalLength() const {return vector_->getLocalLength();}
  
  size_t getGlobalLength() const {return vector_->getGlobalLength();}

  size_t getVector(const gid_t *&Ids, const scalar_t *&elements, 
    int &stride) const;

  size_t getVector(int vectorNumber, const gid_t *&Ids, 
    const scalar_t *&elements, int &stride) const
  {
    env_->localInputAssertion(__FILE__, __LINE__, "invalid vector",
      vectorNumber==0, BASIC_ASSERTION);

    return getVector(Ids, elements, stride);
  }

  size_t getVectorWeights(int dim, const scalar_t *&weights, int &stride) const
  {
    env_->localInputAssertion(__FILE__, __LINE__, "invalid dimension",
      dim >= 0 && dim < numWeights_, BASIC_ASSERTION);

    size_t length;

    weights_[dim].getStridedList(length, weights, stride);

    return length;
  }

  template <typename Adapter>
    size_t applyPartitioningSolution(const User &in, User *&out,
         const PartitioningSolution<Adapter> &solution) const;

private:

  RCP<const User> invector_;
  RCP<const x_vector_t> vector_;
  RCP<const Xpetra::Map<lno_t, gno_t, node_t> > map_;
  RCP<Environment> env_;
  lno_t base_;

  int numWeights_;
  ArrayRCP<StridedData<lno_t, scalar_t> > weights_;
};

////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////
  
template <typename User>
  XpetraVectorAdapter<User>::XpetraVectorAdapter(
    const RCP<const User> &invector, 
    vector<const scalar_t *> &weights, vector<int> &weightStrides):
      invector_(invector), vector_(), map_(),
      env_(rcp(new Environment)), base_(),
      numWeights_(weights.size()), weights_()
{
  typedef StridedData<lno_t, scalar_t> input_t;

  vector_ = XpetraTraits<User>::convertToXpetra(invector);
  map_ = vector_->getMap();
  base_ = map_->getIndexBase();

  size_t length = vector_->getLocalLength();

  if (numWeights_ > 0)
    weights_ = arcp(new input_t [numWeights_], 0, numWeights_, true);

  if (length > 0 && numWeights_ > 0){
    int stride = 1;
    for (int w=0; w < numWeights_; w++){
      if (weightStrides.size())
        stride = weightStrides[w];
      ArrayRCP<const scalar_t> wtArray(weights[w], 0, stride*length, false);
      weights_[w] = input_t(wtArray, stride);
    }
  }
}

template <typename User>
  size_t XpetraVectorAdapter<User>::getVector(const gid_t *&Ids, 
    const scalar_t *&elements, int &stride) const
{
  stride = 1;
  elements = NULL;
  const x_vector_t *vec =  vector_.get();

  if (map_->lib() == Xpetra::UseTpetra){
    const xt_vector_t *tvector = dynamic_cast<const xt_vector_t *>(vec);

    if (tvector->getLocalLength() > 0){
      // getData hangs if vector length is 0
      ArrayRCP<const scalar_t> data = tvector->getData(0);
      elements = data.get();
    }
  }
  else if (map_->lib() == Xpetra::UseEpetra){
    const xe_vector_t *evector = dynamic_cast<const xe_vector_t *>(vec);
      
    if (evector->getLocalLength() > 0){
      // getData hangs if vector length is 0
      ArrayRCP<const double> data = evector->getData(0);

      // Cast so this will compile when scalar_t is not double,
      // a case when this code should never execute.
      elements = reinterpret_cast<const scalar_t *>(data.get());
    }
  }
  else{
    throw logic_error("invalid underlying lib");
  }

  ArrayView<const gid_t> gids = map_->getNodeElementList();
  Ids = gids.getRawPtr();

  return getLocalLength();
}

template <typename User>
  template <typename Adapter>
    size_t XpetraVectorAdapter<User>::applyPartitioningSolution(
      const User &in, User *&out, 
      const PartitioningSolution<Adapter> &solution) const
{ 
  // Get an import list

  size_t len = solution.getLocalNumberOfIds();
  const gid_t *gids = solution.getIdList();
  const partId_t *parts = solution.getPartList();
  ArrayRCP<gid_t> gidList = arcp(const_cast<gid_t *>(gids), 0, len, false);
  ArrayRCP<partId_t> partList = arcp(const_cast<partId_t *>(parts), 0, len, 
    false);

  ArrayRCP<lno_t> dummyIn;
  ArrayRCP<gid_t> importList;
  ArrayRCP<lno_t> dummyOut;
  size_t numNewRows;

  const RCP<const Comm<int> > comm = map_->getComm(); 

  try{
    numNewRows = solution.convertSolutionToImportList(
      0, dummyIn, importList, dummyOut);
  }
  Z2_FORWARD_EXCEPTIONS;

  RCP<const User> inPtr = rcp(&in, false);
  lno_t localNumElts = numNewRows;

  RCP<const User> outPtr = XpetraTraits<User>::doMigration(
   inPtr, localNumElts, importList.get());

  out = const_cast<User *>(outPtr.get());
  outPtr.release();
  return numNewRows;
}

}  //namespace Zoltan2
  
#endif

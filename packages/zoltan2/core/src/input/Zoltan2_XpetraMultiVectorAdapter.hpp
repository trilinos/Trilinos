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

/*! \file Zoltan2_XpetraMultiVectorAdapter.hpp
    \brief Defines the XpetraMultiVectorAdapter
*/

#ifndef _ZOLTAN2_XPETRAMULTIVECTORADAPTER_HPP_
#define _ZOLTAN2_XPETRAMULTIVECTORADAPTER_HPP_

#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_VectorAdapter.hpp>
#include <Zoltan2_StridedData.hpp>
#include <Zoltan2_PartitioningHelpers.hpp>

#if defined(HAVE_ZOLTAN2_EPETRA) && defined(HAVE_XPETRA_EPETRA)
#include <Xpetra_EpetraMultiVector.hpp>
#endif
#include <Xpetra_TpetraMultiVector.hpp>

namespace Zoltan2 {

/*!  \brief An adapter for Xpetra::MultiVector.

    The template parameter is the user's input object:
    \li \c Epetra_MultiVector
    \li \c Tpetra::MultiVector
    \li \c Xpetra::MultiVector

    The \c scalar_t type, representing use data such as vector values, is
    used by Zoltan2 for weights, coordinates, part sizes and
    quality metrics.
    Some User types (like Tpetra::CrsMatrix) have an inherent scalar type,
    and some
    (like Tpetra::CrsGraph) do not.  For such objects, the scalar type is
    set by Zoltan2 to \c float.  If you wish to change it to double, set
    the second template parameter to \c double.
*/

template <typename User>
  class XpetraMultiVectorAdapter : public VectorAdapter<User> {
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::part_t   part_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef User user_t;
  typedef User userCoord_t;

  typedef Xpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> x_mvector_t;
  typedef Xpetra::TpetraMultiVector<
    scalar_t, lno_t, gno_t, node_t> xt_mvector_t;
#endif

  /*! \brief Destructor
   */
  ~XpetraMultiVectorAdapter() { }

  /*! \brief Constructor
   *
   *  \param invector  the user's Xpetra, Tpetra or Epetra MultiVector object
   *  \param weights  a list of pointers to arrays of weights.
   *      The number of weights per multivector element is assumed to be
   *      \c weights.size().
   *  \param weightStrides  a list of strides for the \c weights.
   *     The weight for weight index \c n for multivector element
   *     \c k should be found at <tt>weights[n][weightStrides[n] * k]</tt>.
   *     If \c weightStrides.size() is zero, it is assumed all strides are one.
   *
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this Adapter.
   */

  XpetraMultiVectorAdapter(const RCP<const User> &invector,
    std::vector<const scalar_t *> &weights, std::vector<int> &weightStrides);

  /*! \brief Constructor for case when weights are not being used.
   *
   *  \param invector  the user's Xpetra, Tpetra or Epetra MultiVector object
   */

  XpetraMultiVectorAdapter(const RCP<const User> &invector);


  ////////////////////////////////////////////////////
  // The Adapter interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumIDs() const { return vector_->getLocalLength();}

  void getIDsView(const gno_t *&ids) const
  {
    ids = map_->getNodeElementList().getRawPtr();
  }

  void getIDsKokkosView(
    Kokkos::View<const gno_t *, typename node_t::device_type> &ids) const {
    if (map_->lib() == Xpetra::UseTpetra) {
      const xt_mvector_t *tvector =
        dynamic_cast<const xt_mvector_t *>(vector_.get());
      ids = tvector->getTpetra_MultiVector()->getMap()->getMyGlobalIndices();
    }
    else if (map_->lib() == Xpetra::UseEpetra) {
#if defined(HAVE_ZOLTAN2_EPETRA) && defined(HAVE_XPETRA_EPETRA)
      // this will call getIDsView to get raw ptr and create a view from it
      return BaseAdapter<User>::getIDsKokkosView(ids);
#else
      throw std::logic_error("Epetra requested, but Trilinos is not "
                           "built with Epetra");
#endif
    }
    else {
      throw std::logic_error("getIDsKokkosView called but not on Tpetra or Epetra!");
    }
  }

  int getNumWeightsPerID() const { return numWeights_;}

  void getWeightsView(const scalar_t *&weights, int &stride, int idx) const
  {
    if(idx<0 || idx >= numWeights_)
    {
        std::ostringstream emsg;
        emsg << __FILE__ << ":" << __LINE__
             << "  Invalid weight index " << idx << std::endl;
        throw std::runtime_error(emsg.str()); 
    }

    size_t length;
    weights_[idx].getStridedList(length, weights, stride);
  }

  void getWeightsKokkos2dView(Kokkos::View<scalar_t **,
    typename node_t::device_type> &wgt) const {
    typedef Kokkos::View<scalar_t**, typename node_t::device_type> view_t;
    wgt = view_t("wgts", vector_->getLocalLength(), numWeights_);
    typename view_t::HostMirror host_wgt = Kokkos::create_mirror_view(wgt);
    for(int idx = 0; idx < numWeights_; ++idx) {
      const scalar_t * weights;
      size_t length;
      int stride;
      weights_[idx].getStridedList(length, weights, stride);
      size_t fill_index = 0;
      for(size_t n = 0; n < length; n += stride) {
        host_wgt(fill_index++,idx) = weights[n];
      }
    }
    Kokkos::deep_copy(wgt, host_wgt);
  }

  ////////////////////////////////////////////////////
  // The VectorAdapter interface.
  ////////////////////////////////////////////////////

  int getNumEntriesPerID() const {return vector_->getNumVectors();}

  void getEntriesView(const scalar_t *&elements, int &stride, int idx=0) const;

  void getEntriesKokkosView(
    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    Kokkos::View<scalar_t **, Kokkos::LayoutLeft,
    typename node_t::device_type> & elements) const;

  template <typename Adapter>
    void applyPartitioningSolution(const User &in, User *&out,
         const PartitioningSolution<Adapter> &solution) const;

  template <typename Adapter>
    void applyPartitioningSolution(const User &in, RCP<User> &out,
         const PartitioningSolution<Adapter> &solution) const;

private:

  RCP<const User> invector_;
  RCP<const x_mvector_t> vector_;
  RCP<const Xpetra::Map<lno_t, gno_t, node_t> > map_;

  int numWeights_;
  ArrayRCP<StridedData<lno_t, scalar_t> > weights_;
};

////////////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////////////

template <typename User>
  XpetraMultiVectorAdapter<User>::XpetraMultiVectorAdapter(
    const RCP<const User> &invector,
    std::vector<const scalar_t *> &weights, std::vector<int> &weightStrides):
      invector_(invector), vector_(), map_(),
      numWeights_(weights.size()), weights_(weights.size())
{
  typedef StridedData<lno_t, scalar_t> input_t;

  try {
    RCP<x_mvector_t> tmp =
            XpetraTraits<User>::convertToXpetra(rcp_const_cast<User>(invector));
    vector_ = rcp_const_cast<const x_mvector_t>(tmp);
  }
  Z2_FORWARD_EXCEPTIONS

  map_ = vector_->getMap();

  size_t length = vector_->getLocalLength();

  if (length > 0 && numWeights_ > 0){
    int stride = 1;
    for (int w=0; w < numWeights_; w++){
      if (weightStrides.size())
        stride = weightStrides[w];
      ArrayRCP<const scalar_t> wgtV(weights[w], 0, stride*length, false);
      weights_[w] = input_t(wgtV, stride);
    }
  }
}


////////////////////////////////////////////////////////////////////////////
template <typename User>
  XpetraMultiVectorAdapter<User>::XpetraMultiVectorAdapter(
    const RCP<const User> &invector):
      invector_(invector), vector_(), map_(),
      numWeights_(0), weights_()
{
  try {
    RCP<x_mvector_t> tmp =
            XpetraTraits<User>::convertToXpetra(rcp_const_cast<User>(invector));
    vector_ = rcp_const_cast<const x_mvector_t>(tmp);
  }
  Z2_FORWARD_EXCEPTIONS

  map_ = vector_->getMap();
}

////////////////////////////////////////////////////////////////////////////
template <typename User>
  void XpetraMultiVectorAdapter<User>::getEntriesView(
    const scalar_t *&elements, int &stride, int idx) const
{
  size_t vecsize;
  stride = 1;
  elements = NULL;
  if (map_->lib() == Xpetra::UseTpetra){
    const xt_mvector_t *tvector =
      dynamic_cast<const xt_mvector_t *>(vector_.get());

    vecsize = tvector->getLocalLength();
    if (vecsize > 0){
      ArrayRCP<const scalar_t> data = tvector->getData(idx);
      elements = data.get();
    }
  }
  else if (map_->lib() == Xpetra::UseEpetra){
#if defined(HAVE_ZOLTAN2_EPETRA) && defined(HAVE_XPETRA_EPETRA)
    typedef Xpetra::EpetraMultiVectorT<gno_t,node_t> xe_mvector_t;
    const xe_mvector_t *evector =
      dynamic_cast<const xe_mvector_t *>(vector_.get());

    vecsize = evector->getLocalLength();
    if (vecsize > 0){
      ArrayRCP<const double> data = evector->getData(idx);

      // Cast so this will compile when scalar_t is not double,
      // a case when this code should never execute.
      elements = reinterpret_cast<const scalar_t *>(data.get());
    }
#else
    throw std::logic_error("Epetra requested, but Trilinos is not "
                           "built with Epetra");
#endif
  }
  else{
    throw std::logic_error("invalid underlying lib");
  }
}

////////////////////////////////////////////////////////////////////////////
template <typename User>
  void XpetraMultiVectorAdapter<User>::getEntriesKokkosView(
    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    Kokkos::View<scalar_t **, Kokkos::LayoutLeft, typename node_t::device_type> & elements) const
{
  if (map_->lib() == Xpetra::UseTpetra){
      const xt_mvector_t *tvector =
        dynamic_cast<const xt_mvector_t *>(vector_.get());
    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    Kokkos::View<scalar_t **, Kokkos::LayoutLeft, typename node_t::device_type> view2d =
      tvector->getTpetra_MultiVector()->template getLocalView<node_t>();
    elements = view2d;
  }
  else if (map_->lib() == Xpetra::UseEpetra){
#if defined(HAVE_ZOLTAN2_EPETRA) && defined(HAVE_XPETRA_EPETRA)
    typedef Xpetra::EpetraMultiVectorT<gno_t,node_t> xe_mvector_t;
    const xe_mvector_t *evector =
      dynamic_cast<const xe_mvector_t *>(vector_.get());
    elements =
      Kokkos::View<scalar_t **, Kokkos::LayoutLeft, typename node_t::device_type>
      ("elements", evector->getLocalLength(), evector->getNumVectors());
    if(evector->getLocalLength() > 0) {
      for(size_t idx = 0; idx < evector->getNumVectors(); ++idx) {
        const scalar_t * ptr;
        int stride;
        getEntriesView(ptr, stride, idx);
        for(size_t n = 0; n < evector->getLocalLength(); ++n) {
          elements(n, idx) = ptr[n];
        }
      }
    }
#else
    throw std::logic_error("Epetra requested, but Trilinos is not "
                           "built with Epetra");
#endif
  }
  else {
    throw std::logic_error("getEntriesKokkosView called but not using Tpetra or Epetra!");
  }
}

////////////////////////////////////////////////////////////////////////////
template <typename User>
  template <typename Adapter>
    void XpetraMultiVectorAdapter<User>::applyPartitioningSolution(
      const User &in, User *&out,
      const PartitioningSolution<Adapter> &solution) const
{
  // Get an import list (rows to be received)
  size_t numNewRows;
  ArrayRCP<gno_t> importList;
  try{
    numNewRows = Zoltan2::getImportList<Adapter,
                                        XpetraMultiVectorAdapter<User> >
                                       (solution, this, importList);
  }
  Z2_FORWARD_EXCEPTIONS;

  // Move the rows, creating a new vector.
  RCP<User> outPtr = XpetraTraits<User>::doMigration(in, numNewRows,
                                                     importList.getRawPtr());
  out = outPtr.get();
  outPtr.release();
}

////////////////////////////////////////////////////////////////////////////
template <typename User>
  template <typename Adapter>
    void XpetraMultiVectorAdapter<User>::applyPartitioningSolution(
      const User &in, RCP<User> &out,
      const PartitioningSolution<Adapter> &solution) const
{
  // Get an import list (rows to be received)
  size_t numNewRows;
  ArrayRCP<gno_t> importList;
  try{
    numNewRows = Zoltan2::getImportList<Adapter,
                                        XpetraMultiVectorAdapter<User> >
                                       (solution, this, importList);
  }
  Z2_FORWARD_EXCEPTIONS;

  // Move the rows, creating a new vector.
  out = XpetraTraits<User>::doMigration(in, numNewRows,
                                        importList.getRawPtr());
}

}  //namespace Zoltan2

#endif

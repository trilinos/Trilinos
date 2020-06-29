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

/*! \file Zoltan2_CoordinateModel.hpp
    \brief Defines the CoordinateModel classes.
*/


#ifndef _ZOLTAN2_COORDINATEMODEL_HPP_
#define _ZOLTAN2_COORDINATEMODEL_HPP_

#include <Zoltan2_Model.hpp>
#include <Zoltan2_MeshAdapter.hpp>
#include <Zoltan2_MatrixAdapter.hpp>
#include <Zoltan2_GraphAdapter.hpp>
#include <Zoltan2_IdentifierAdapter.hpp>
#include <Zoltan2_VectorAdapter.hpp>
#include <Zoltan2_StridedData.hpp>

namespace Zoltan2 {

/*!  \brief This class provides geometric coordinates with optional weights
           to the Zoltan2 algorithm.

    The template parameter is an Input Adapter.  Input adapters are
    templated on the basic user input type.
*/
template <typename Adapter>
class CoordinateModel : public Model<Adapter>
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename Adapter::scalar_t    scalar_t;
  typedef typename Adapter::gno_t       gno_t;
  typedef typename Adapter::lno_t       lno_t;
  typedef typename Adapter::node_t      node_t;
  typedef typename Adapter::user_t      user_t;
  typedef typename Adapter::userCoord_t userCoord_t;
  typedef StridedData<lno_t, scalar_t>  input_t;
#endif

  ////////////////////////////////////////////////////
  // Constructors for each Adapter type
  ////////////////////////////////////////////////////

  // VectorAdapter
  CoordinateModel(const RCP<const VectorAdapter<user_t> > &ia,
                  const RCP<const Environment> &env,
                  const RCP<const Comm<int> > &comm,
                  modelFlag_t &flags):
                  numGlobalCoordinates_(), env_(env), comm_(comm),
                  coordinateDim_(), gids_(), 
                  xyz_(), userNumWeights_(0), weights_()
  {
    typedef VectorAdapter<user_t> adapterWithCoords_t;
    sharedConstructor<adapterWithCoords_t>(&(*ia), env, comm, flags);
  }

  // MatrixAdapter
  CoordinateModel(const RCP<const MatrixAdapter<user_t,userCoord_t> > &ia,
                  const RCP<const Environment> &env,
                  const RCP<const Comm<int> > &comm,
                  modelFlag_t &flags) :
                  numGlobalCoordinates_(), env_(env), comm_(comm),
                  coordinateDim_(), gids_(), 
                  xyz_(), userNumWeights_(0), weights_()
  {
    if (!(ia->coordinatesAvailable()))
      throw std::logic_error("No coordinate info provided to MatrixAdapter.");
    else {
      typedef VectorAdapter<userCoord_t> adapterWithCoords_t;
      adapterWithCoords_t *va = ia->getCoordinateInput();
      sharedConstructor<adapterWithCoords_t>(va, env, comm, flags);
    }
  }

  // GraphAdapter
  CoordinateModel(const RCP<const GraphAdapter<user_t,userCoord_t> > &ia,
                  const RCP<const Environment> &env,
                  const RCP<const Comm<int> > &comm,
                  modelFlag_t &flags) :
                  numGlobalCoordinates_(), env_(env), comm_(comm),
                  coordinateDim_(), gids_(), 
                  xyz_(), userNumWeights_(0), weights_()
  {
    if (!(ia->coordinatesAvailable()))
      throw std::logic_error("No coordinate info provided to GraphAdapter.");
    else {
      typedef VectorAdapter<userCoord_t> adapterWithCoords_t;
      adapterWithCoords_t *va = ia->getCoordinateInput();
      sharedConstructor<adapterWithCoords_t>(va, env, comm, flags);
    }
  }

  // MeshAdapter
  CoordinateModel(const RCP<const MeshAdapter<user_t> > &ia,
		  const RCP<const Environment> &env,
		  const RCP<const Comm<int> > &comm,
		  modelFlag_t &flags) :
                  numGlobalCoordinates_(), env_(env), comm_(comm),
                  coordinateDim_(), gids_(), 
                  xyz_(), userNumWeights_(0), weights_()
  {
    typedef MeshAdapter<user_t> adapterWithCoords_t;
    sharedConstructor<adapterWithCoords_t>(&(*ia), env, comm, flags);
  }

  // IdentifierAdapter
  CoordinateModel(const RCP<const IdentifierAdapter<user_t> > &ia,
                  const RCP<const Environment> &env,
                  const RCP<const Comm<int> > &comm,
                  modelFlag_t &flags)
  {
    throw std::logic_error(
      "A coordinate model can not be build from an IdentifierAdapter");
  }

  ////////////////////////////////////////////////////
  // CoordinateModel interface.
  ////////////////////////////////////////////////////

  /*! \brief Returns the dimension of the coordinates.
   */
  int getCoordinateDim() const { return coordinateDim_;}

  /*! \brief Returns the number of coordinates on this process.
   */
  size_t getLocalNumCoordinates() const { return gids_.size();}

  /*! \brief Returns the global number coordinates.
   */
  global_size_t getGlobalNumCoordinates() const {return numGlobalCoordinates_;}

  /*! \brief Returns the number (0 or greater) of weights per coordinate
   */
  int getNumWeightsPerCoordinate() const { return userNumWeights_;}

  /*! \brief Returns the coordinate ids, values and optional weights.

      \param Ids will on return point to the list of the global Ids for
        each coordinate on this process.

      \param xyz on return is a list of getCoordinateDim()
          StridedData objects, each containing the coordinates for
          one dimension. If the coordinate dimension is three, then
          the coordinates for <tt>Ids[k]</tt> are
          <tt>xyz[0][k], xyz[1][k], xyz[2][k]</tt>.

      \param wgts on return is a list of getNumWeightsPerCoordinate()
          StridedData objects, each containing the weights for
          one weight index.  For the index \d,
          the weight for <tt>Ids[k]</tt> is <tt>wgts[d][k]</tt>.

       \return The number of ids in the Ids list.

      Memory for this data is allocated either by the user or the Model.
      The caller gets a view of the data.
   */

  size_t getCoordinates(ArrayView<const gno_t>  &Ids,
                        ArrayView<input_t> &xyz,
                        ArrayView<input_t> &wgts) const
  {
    xyz = xyz_.view(0, coordinateDim_);
    wgts = weights_.view(0, userNumWeights_);

    size_t nCoord = getLocalNumCoordinates();
    Ids =  ArrayView<const gno_t>();

    if (nCoord){
      Ids = Teuchos::arrayView<const gno_t>(
        reinterpret_cast<const gno_t *>(gids_.getRawPtr()), nCoord);
    }

    return nCoord;
  }

  /*! \brief Returns the coordinate ids, values and optional weights.
      \param Ids will on return point to the list of the global Ids for
        each coordinate on this process.
      \param xyz on return is a view of xyz coordinates.
      \param wgts on return is a view of the weights.
      \return The number of ids in the Ids list
   */
  size_t getCoordinatesKokkos(
    Kokkos::View<const gno_t *, typename node_t::device_type> &Ids,
    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    Kokkos::View<scalar_t **,
      Kokkos::LayoutLeft, typename node_t::device_type> &xyz,
    Kokkos::View<scalar_t **, typename node_t::device_type> &wgts) const
  {
    Ids = kokkos_gids_;
    xyz = kokkos_xyz_;
    wgts = kokkos_weights_;
    return getLocalNumCoordinates();
  }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const
  {
    return getLocalNumCoordinates();
  }

  size_t getGlobalNumObjects() const
  {
    return getGlobalNumCoordinates();
  }

private:
  size_t numGlobalCoordinates_;
  const RCP<const Environment> env_;
  const RCP<const Comm<int> > comm_;
  int coordinateDim_;

  // TODO: We now have a Kokkos version and non kokkos version so need to clean
  // this up and perhaps eliminate the non-kokkos version completely.
  // However not all tests are converted to Kokkos so keeping both forms around
  // for now is probably necessary.
  Kokkos::View<const gno_t *, typename node_t::device_type> kokkos_gids_;
  // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
  Kokkos::View<scalar_t **, Kokkos::LayoutLeft, typename node_t::device_type> kokkos_xyz_;
  Kokkos::View<scalar_t **, typename node_t::device_type> kokkos_weights_;

  ArrayRCP<const gno_t> gids_;
  ArrayRCP<input_t> xyz_;
  int userNumWeights_;
  ArrayRCP<input_t> weights_;

  template <typename AdapterWithCoords>
  void sharedConstructor(const AdapterWithCoords *ia,
                         const RCP<const Environment> &env,
                         const RCP<const Comm<int> > &comm,
                         modelFlag_t &flags);

};


////////////////////////////////////////////////////////////////////////////

// sharedConstructor
template <typename Adapter>
template <typename AdapterWithCoords>
void CoordinateModel<Adapter>::sharedConstructor(
    const AdapterWithCoords *ia,
    const RCP<const Environment> &/* env */,
    const RCP<const Comm<int> > &comm,
    modelFlag_t &/* flags */)
{
  size_t nLocalIds = ia->getLocalNumIDs();

  // Get coordinates and weights (if any)

  int tmp[2], gtmp[2];
  tmp[0] = ia->getDimension();
  tmp[1] = ia->getNumWeightsPerID();
  Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_MAX, 2, tmp, gtmp);
  coordinateDim_ = gtmp[0];
  userNumWeights_ = gtmp[1];

  env_->localBugAssertion(__FILE__, __LINE__, "coordinate dimension",
    coordinateDim_ > 0, COMPLEX_ASSERTION);

  input_t *coordArray = new input_t [coordinateDim_];
  input_t *weightArray = NULL;
  if (userNumWeights_)
    weightArray = new input_t [userNumWeights_];

  env_->localMemoryAssertion(__FILE__, __LINE__, userNumWeights_+coordinateDim_,
    coordArray && (!userNumWeights_|| weightArray));


  if (nLocalIds){

    ia->getIDsKokkosView(kokkos_gids_);
    ia->getCoordinatesKokkosView(kokkos_xyz_);
    if(userNumWeights_ > 0) {
      ia->getWeightsKokkosView(kokkos_weights_);
    }

    const gno_t *gids=NULL;

    ia->getIDsView(gids);
    gids_ = arcp(gids, 0, nLocalIds, false);

    for (int dim=0; dim < coordinateDim_; dim++){
      int stride;
      const scalar_t *coords=NULL;
      try{
        ia->getCoordinatesView(coords, stride, dim);
      }
      Z2_FORWARD_EXCEPTIONS;

      ArrayRCP<const scalar_t> cArray(coords, 0, nLocalIds*stride, false);
      coordArray[dim] = input_t(cArray, stride);
    }

    for (int idx=0; idx < userNumWeights_; idx++){
      int stride;
      const scalar_t *weights;
      try{
        ia->getWeightsView(weights, stride, idx);
      }
      Z2_FORWARD_EXCEPTIONS;

      ArrayRCP<const scalar_t> wArray(weights, 0, nLocalIds*stride, false);
      weightArray[idx] = input_t(wArray, stride);
    }
  }

  xyz_ = arcp(coordArray, 0, coordinateDim_);

  if (userNumWeights_)
    weights_ = arcp(weightArray, 0, userNumWeights_);
 
  Teuchos::reduceAll<int, size_t>(*comm, Teuchos::REDUCE_SUM, 1, 
                                  &nLocalIds, &numGlobalCoordinates_);

  env_->memory("After construction of coordinate model");
}

}   // namespace Zoltan2

#endif

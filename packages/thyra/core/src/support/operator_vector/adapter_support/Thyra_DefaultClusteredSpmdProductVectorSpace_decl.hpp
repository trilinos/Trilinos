// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_CLUSTERED_SPMD_PRODUCT_VECTOR_SPACE_DECL_HPP
#define THYRA_DEFAULT_CLUSTERED_SPMD_PRODUCT_VECTOR_SPACE_DECL_HPP

#include "Thyra_VectorSpaceBase_decl.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Thyra_VectorSpaceDefaultBase.hpp"
#include "Teuchos_Comm.hpp"


namespace Thyra {


/** \brief Concrete subclass of <tt>VectorSpaceBase</tt> that takes a
 * collection of individual <tt>VectorSpaceBase</tt> objects distributed over
 * many different processes and creates a single vector space.
 *
 * What this class does is to take different vector space objects distributed
 * over a set of different clusters of processes and then creates a single
 * vector space object.
 *
 * Let <tt>numClusters</tt> be the number of clusters that the processes in
 * the global communicator are partitioned into.  Each cluster of processes is
 * represented by another communicator known in this process as
 * <tt>intraClusterComm</tt>.  The communicator that links the root of each
 * cluster is called <tt>interClusterComm</tt>. All communication is performed
 * with just these two communicators.  There is no overall global communicator
 * that encompasses all of the processes.
 *
 * Within a cluster of processes, the number of constituent vector space
 * objects must be the same. Let <tt>v[0]..v[numBlocks-1]</tt> be the
 * <tt>numBlocks</tt> <tt>VectorBase</tt> objects for constituent vectors in
 * this cluster of processes.  There is no assumption make whatsoever about
 * the natrue of the <tt>VectorSpaceBase</tt> objects or the
 * <tt>[Multi]VectorBase</tt> objects that they create.
 *
 * ToDo: Finish documentation!
 *
 * The default copy constructor and assign operator are allowed and they work
 * correctly and perform shallow copies of the constituent vector spaces!
 *
 * \ingroup Thyra_Op_Vec_adapters_Spmd_concrete_std_grp
 */
template<class Scalar>
class DefaultClusteredSpmdProductVectorSpace
  : public ProductVectorSpaceBase<Scalar>
  , protected VectorSpaceDefaultBase<Scalar>
{
public:

  /** \name Constructors/Intializers/Accessors */
  //@{

  /** \brief . */
  DefaultClusteredSpmdProductVectorSpace();

  /** \brief . */
  DefaultClusteredSpmdProductVectorSpace(
    const Teuchos::RCP<const Teuchos::Comm<Ordinal> >          &intraClusterComm
    ,const int                                                       clusterRootRank
    ,const Teuchos::RCP<const Teuchos::Comm<Ordinal> >         &interClusterComm
    ,const int                                                       numBlocks
    ,const Teuchos::RCP<const VectorSpaceBase<Scalar> >      vecSpaces[]
    );

  /** \brief Initalize.
   *
   * \param  intraClusterComm
   *            [in] The communicator over just this cluster.
   * \param  clusterRootRank
   *            [in] The rank of the root process in <tt>*interClusterComm</tt> for this
   *            cluster of processes.  This is also the process that in included in
   *            <tt>*interClusterComm</tt> (which has a different rank ther obviously). 
   * \param  interClusterComm
   *            [in] Defines the communicator between the root processes of all of the
   *            clusters.  For the root process of each cluster
   *            <tt>*interClusterComm!=Spmd_COMM_NULL</tt>,
   *            otherwise <tt>interClusterComm==Teuchos::null</tt> or 
   *            <tt>*interClusterComm==Spmd_COMM_NULL</tt>.
   * \param  numBlocks
   *            [in] Number of vector space blocks for this cluster of processes.
   * \param  vecSpaces
   *            [in] Array (length <tt>numBlocks</tt>) of the vector spaces
   *            for this cluster of processes.
   */
  void initialize(
    const Teuchos::RCP<const Teuchos::Comm<Ordinal> >          &intraClusterComm
    ,const int                                                       clusterRootRank
    ,const Teuchos::RCP<const Teuchos::Comm<Ordinal> >         &interClusterComm
    ,const int                                                       numBlocks
    ,const Teuchos::RCP<const VectorSpaceBase<Scalar> >      vecSpaces[]
    );

  /** \brief . */
  Teuchos::RCP<const Teuchos::Comm<Ordinal> > intraClusterComm() const;

  /** \brief . */
  int clusterRootRank() const;

  /** \brief . */
  Teuchos::RCP<const Teuchos::Comm<Ordinal> > interClusterComm() const;

  /** \bief The some of the dimensions of the block vector spaces in this
   * cluster.
   */
  int clusterSubDim() const;

  /** \bief The offset of the first element in the first constituent vector in
   * this cluster in the w.r.t. the global vector.
   */
  int clusterOffset() const;
  
  //@}
  
  /** @name Overridden form Teuchos::Describable */
  //@{
  /** \brief . */
  std::string description() const;
  //@}

  /** @name Public overridden from VectorSpaceBase */
  //@{
  /** \brief . */
  Ordinal dim() const;
  /** \brief . */
  bool isCompatible(const VectorSpaceBase<Scalar>& vecSpc) const;
  /** \brief . */
  Teuchos::RCP< const VectorSpaceFactoryBase<Scalar> > smallVecSpcFcty() const;
  /** \brief . */
  Scalar scalarProd( const VectorBase<Scalar>& x, const VectorBase<Scalar>& y ) const;
  /** \brief . */
  void scalarProdsImpl(
    const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
    const ArrayView<Scalar> &scalarProds ) const;
  /** \brief . */
  bool isEuclidean() const;
  /** \brief . */
  bool hasInCoreView(
    const Range1D& rng, const EViewType viewType, const EStrideType strideType
    ) const;
  /** \brief . */
  Teuchos::RCP< const VectorSpaceBase<Scalar> > clone() const;
  //@}

  /** @name Protected overridden from ProductVectorSpaceBase */
  //@{

  /** \brief . */
  int numBlocks() const;
  /** \brief . */
  Teuchos::RCP<const VectorSpaceBase<Scalar> > getBlock(const int k) const; 

  //@}

protected:

  /** @name Protected overridden from VectorSpaceBase */
  //@{

  /** \brief . */
  Teuchos::RCP<VectorBase<Scalar> > createMember() const;
  /** \brief . */
  Teuchos::RCP<MultiVectorBase<Scalar> > createMembers(int numMembers) const;

  //@}

private:

  // //////////////////////////////////////
  // Private types

  typedef std::vector<Teuchos::RCP<const VectorSpaceBase<Scalar> > >  vecSpaces_t;

  // //////////////////////////////////////
  // Private data members

  Teuchos::RCP<const Teuchos::Comm<Ordinal> >  intraClusterComm_;
  int                                                      clusterRootRank_;
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >  interClusterComm_;
  vecSpaces_t vecSpaces_; // size == numBlocks
  bool isEuclidean_;
  Ordinal globalDim_;     // The global dimension of all of the block vectors in
                        // all of the clusters.
  Ordinal clusterSubDim_; // The some of the dimensions of the block vector
                        // spaces in this cluster
  Ordinal clusterOffset_; // The offset of the first element in the first
                        // constituent vector in this cluster in the
                        // w.r.t. the global vector.
  
};


// ///////////////////////////
// Inline defintions


template<class Scalar>
Teuchos::RCP<const Teuchos::Comm<Ordinal> >
DefaultClusteredSpmdProductVectorSpace<Scalar>::intraClusterComm() const
{
  return intraClusterComm_;
}


template<class Scalar>
int DefaultClusteredSpmdProductVectorSpace<Scalar>::clusterRootRank() const
{
  return clusterRootRank_;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::Comm<Ordinal> >
DefaultClusteredSpmdProductVectorSpace<Scalar>::interClusterComm() const
{
  return interClusterComm_;
}


template<class Scalar>
int DefaultClusteredSpmdProductVectorSpace<Scalar>::clusterSubDim() const
{
  return clusterSubDim_;
}


template<class Scalar>
int DefaultClusteredSpmdProductVectorSpace<Scalar>::clusterOffset() const
{
  return clusterOffset_;
}


} // end namespace Thyra


#endif // THYRA_DEFAULT_CLUSTERED_SPMD_PRODUCT_VECTOR_SPACE_DECL_HPP

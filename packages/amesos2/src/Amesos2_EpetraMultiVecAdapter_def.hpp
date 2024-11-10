// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
  \file   Amesos2_EpetraMultiVecAdapter_def.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   

  \brief  Amesos2::MultiVecAdapter specialization for the
          Epetra_MultiVector class.
*/

#ifndef AMESOS2_EPETRA_MULTIVEC_ADAPTER_DEF_HPP
#define AMESOS2_EPETRA_MULTIVEC_ADAPTER_DEF_HPP

#include <Teuchos_as.hpp>

#include <Epetra_SerialComm.h>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#endif
#include <Epetra_LocalMap.h>
#include <Epetra_Import.h>
#include <Epetra_Export.h>

#include "Amesos2_EpetraMultiVecAdapter_decl.hpp"
#include "Amesos2_Util.hpp"

namespace Amesos2 {

MultiVecAdapter<Epetra_MultiVector>::MultiVecAdapter(const MultiVecAdapter<multivec_t>& adapter)
  : mv_(adapter.mv_)
  , mv_map_(adapter.mv_map_)
{ }

MultiVecAdapter<Epetra_MultiVector>::MultiVecAdapter(const Teuchos::RCP<multivec_t>& m)
  : mv_(m)
{
  mv_map_ = Teuchos::rcpFromRef(mv_->Map());
}
  

Teuchos::RCP<Epetra_MultiVector>
MultiVecAdapter<Epetra_MultiVector>::clone() const
{
  Teuchos::RCP<Epetra_MultiVector> Y (new Epetra_MultiVector(*mv_map_, mv_->NumVectors(), false));
  return Y;
}



bool MultiVecAdapter<Epetra_MultiVector>::isLocallyIndexed() const
{
  return !mv_->DistributedGlobal();
}

bool MultiVecAdapter<Epetra_MultiVector>::isGloballyIndexed() const
{
  return mv_->DistributedGlobal();
}


const Teuchos::RCP<const Teuchos::Comm<int> >
MultiVecAdapter<Epetra_MultiVector>::getComm() const
{
  return Util::to_teuchos_comm(Teuchos::rcpFromRef(mv_->Comm()));
}


size_t MultiVecAdapter<Epetra_MultiVector>::getLocalLength() const
{
  return Teuchos::as<size_t>(mv_->MyLength());
}


size_t MultiVecAdapter<Epetra_MultiVector>::getLocalNumVectors() const
{
  return Teuchos::as<size_t>(mv_->NumVectors());
}


MultiVecAdapter<Epetra_MultiVector>::global_size_t
MultiVecAdapter<Epetra_MultiVector>::getGlobalLength() const
{
  return Teuchos::as<global_size_t>(mv_->GlobalLength());
}


size_t MultiVecAdapter<Epetra_MultiVector>::getGlobalNumVectors() const
{
  return Teuchos::as<size_t>(mv_->NumVectors());
}


size_t MultiVecAdapter<Epetra_MultiVector>::getStride() const
{
  return Teuchos::as<size_t>(mv_->Stride());
}


bool MultiVecAdapter<Epetra_MultiVector>::isConstantStride() const
{
  return mv_->ConstantStride();
}


Teuchos::RCP<const Tpetra::Vector<MultiVecAdapter<Epetra_MultiVector>::scalar_t,
                                  MultiVecAdapter<Epetra_MultiVector>::local_ordinal_t,
                                  MultiVecAdapter<Epetra_MultiVector>::global_ordinal_t,
                                  MultiVecAdapter<Epetra_MultiVector>::node_t> >
MultiVecAdapter<Epetra_MultiVector>::getVector( size_t j ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayRCP;
  using Tpetra::MultiVector;

  typedef scalar_t st;
  typedef local_ordinal_t lot;
  typedef global_ordinal_t got;
  typedef node_t nt;
  
  RCP<MultiVector<st,lot,got,nt> > vector = rcp(new MultiVector<st,lot,got,nt>(this->getMap(),1));

  // Copy vector contents into Tpetra multi-vector
  ArrayRCP<st> it = vector->getDataNonConst(0);
  double* vector_data = mv_->operator[](Teuchos::as<int>(j)); // values from j^th vector
  Tpetra::global_size_t size = vector->getGlobalLength();

  for( Tpetra::global_size_t i = 0; i < size; ++i ){
    *it = vector_data[i];
  }

  return vector->getVector(j);
}


// Implementation is essentially the same as getVector()
Teuchos::RCP<Tpetra::Vector<MultiVecAdapter<Epetra_MultiVector>::scalar_t,
                            MultiVecAdapter<Epetra_MultiVector>::local_ordinal_t,
                            MultiVecAdapter<Epetra_MultiVector>::global_ordinal_t,
                            MultiVecAdapter<Epetra_MultiVector>::node_t> >
MultiVecAdapter<Epetra_MultiVector>::getVectorNonConst( size_t j )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayRCP;
  using Tpetra::MultiVector;

  typedef scalar_t st;
  typedef local_ordinal_t lot;
  typedef global_ordinal_t got;
  typedef node_t nt;
  
  RCP<MultiVector<st,lot,got,nt> > vector = rcp(new MultiVector<st,lot,got,nt>(this->getMap(),1));

  // Copy vector contents into Tpetra multi-vector
  ArrayRCP<st> it = vector->getDataNonConst(0);
  double* vector_data = mv_->operator[](Teuchos::as<int>(j)); // values from j^th vector
  Tpetra::global_size_t size = vector->getGlobalLength();

  for( Tpetra::global_size_t i = 0; i < size; ++i ){
    *it = vector_data[i];
  }

  return vector->getVectorNonConst(j);
}


double * MultiVecAdapter<Epetra_MultiVector>::getMVPointer_impl() const
{
  TEUCHOS_TEST_FOR_EXCEPTION( this->getGlobalNumVectors() != 1,
		      std::invalid_argument,
		      "Amesos2_EpetraMultiVectorAdapter: getMVPointer_impl should only be called for case with a single vector and single MPI process" );

  double* vector_data = mv_->operator[](Teuchos::as<int>(0)); // raw pointer to data from 0^th vector
  return vector_data;
}


void MultiVecAdapter<Epetra_MultiVector>::get1dCopy(
  const Teuchos::ArrayView<MultiVecAdapter<Epetra_MultiVector>::scalar_t>& av,
  size_t lda,
  Teuchos::Ptr<
    const Tpetra::Map<MultiVecAdapter<Epetra_MultiVector>::local_ordinal_t,
                      MultiVecAdapter<Epetra_MultiVector>::global_ordinal_t,
                      MultiVecAdapter<Epetra_MultiVector>::node_t> > distribution_map,
                      EDistribution /* distribution */) const
{
  using Teuchos::rcpFromPtr;
  using Teuchos::as;

  const size_t num_vecs = getGlobalNumVectors();
  
#ifdef HAVE_AMESOS2_DEBUG
  const size_t requested_vector_length = distribution_map->getLocalNumElements();
  TEUCHOS_TEST_FOR_EXCEPTION( lda < requested_vector_length,
		      std::invalid_argument,
		      "Given stride is not large enough for local vector length" );
  TEUCHOS_TEST_FOR_EXCEPTION( as<size_t>(av.size()) < (num_vecs-1) * lda + requested_vector_length,
		      std::invalid_argument,
		      "MultiVector storage not large enough given leading dimension "
		      "and number of vectors" );
#endif

  // Optimization for ROOTED and single MPI process
  if ( num_vecs == 1 && mv_->Comm().MyPID() == 0 && mv_->Comm().NumProc() == 1 ) {
	  mv_->ExtractCopy(av.getRawPtr(), lda);
  }
  else {
  Epetra_Map e_dist_map
    = *Util::tpetra_map_to_epetra_map<local_ordinal_t,
                                      global_ordinal_t,
                                      global_size_t,
                                      node_t>(*distribution_map);

  multivec_t redist_mv(e_dist_map, as<int>(num_vecs));
  const Epetra_Import importer(e_dist_map, *mv_map_); // Note, target/source order is reversed in Tpetra
  redist_mv.Import(*mv_, importer, Insert);

  // Finally, do copy
  redist_mv.ExtractCopy(av.getRawPtr(), lda);
  }

}

void MultiVecAdapter<Epetra_MultiVector>::get1dCopy_kokkos_view_host(
  Kokkos::View<scalar_t**, Kokkos::LayoutLeft, Kokkos::HostSpace> & host_view,
  size_t lda,
  Teuchos::Ptr<
    const Tpetra::Map<MultiVecAdapter<Epetra_MultiVector>::local_ordinal_t,
                      MultiVecAdapter<Epetra_MultiVector>::global_ordinal_t,
                      MultiVecAdapter<Epetra_MultiVector>::node_t> > distribution_map,
                      EDistribution /* distribution */) const
{
    using Teuchos::rcpFromPtr;
    using Teuchos::as;

    const size_t num_vecs = getGlobalNumVectors();

  #ifdef HAVE_AMESOS2_DEBUG
    const size_t requested_vector_length = distribution_map->getLocalNumElements();
    TEUCHOS_TEST_FOR_EXCEPTION( lda < requested_vector_length,
            std::invalid_argument,
            "Given stride is not large enough for local vector length" );
  #endif

    // First make a host view
    host_view = Kokkos::View<scalar_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>(
      Kokkos::ViewAllocateWithoutInitializing("get1dCopy_kokkos_view"),
      distribution_map->getLocalNumElements(), num_vecs);

    // Optimization for ROOTED and single MPI process
    if ( num_vecs == 1 && this->mv_->Comm().MyPID() == 0 && this->mv_->Comm().NumProc() == 1 ) {
	    mv_->ExtractCopy(host_view.data(), lda);
    }
    else {
      Epetra_Map e_dist_map
        = *Util::tpetra_map_to_epetra_map<local_ordinal_t,
                                          global_ordinal_t,
                                          global_size_t,
                                          node_t>(*distribution_map);

      multivec_t redist_mv(e_dist_map, as<int>(num_vecs));
      const Epetra_Import importer(e_dist_map, *mv_map_); // Note, target/source order is reversed in Tpetra
      redist_mv.Import(*mv_, importer, Insert);

      // Finally, access data

      // Can we consider direct ptr usage with ExtractView?
      // For now I will just copy - this was discussed as low priority for now.
      redist_mv.ExtractCopy(host_view.data(), lda);
    }
}

Teuchos::ArrayRCP<MultiVecAdapter<Epetra_MultiVector>::scalar_t>
MultiVecAdapter<Epetra_MultiVector>::get1dViewNonConst(bool local)
{
  ((void) local);
  // TEUCHOS_TEST_FOR_EXCEPTION( !this->isConstantStride(),
  //   std::logic_error,
  //   "get1dViewNonConst() : can only get 1d view if stride is constant");

  // if( local ){
  //   TEUCHOS_TEST_FOR_EXCEPTION(
  //     true,
  //     std::logic_error,
  //     "Amesos2::MultiVecAdapter<Epetra_MultiVector> : 1D views not yet supported for local-local Epetra multi-vectors");
    
  //   // localize();
  //   // /* Use the global element list returned by
  //   //  * mv_->getMap()->getLocalElementList() to get a subCopy of mv_ which we
  //   //  * assign to l_l_mv_, then return get1dViewNonConst() of l_l_mv_
  //   //  */
  //   // l_l_mv_ = Teuchos::null;
    
  //   // Teuchos::Array<GlobalOrdinal> nodeElements_go(mv_->Map().NumMyElements());
  //   // mv_->Map().MyGlobalElements(nodeElements_go.getRawPtr());
  //   // Teuchos::Array<size_t> nodeElements_st(nodeElements_go.size());

  //   // // Convert the node element to a list of size_t type objects
  //   // typename Teuchos::ArrayView<const GlobalOrdinal>::iterator it_go;
  //   // Teuchos::Array<size_t>::iterator it_st = nodeElements_st.begin();
  //   // for( it_go = nodeElements_go.begin(); it_go != nodeElements_go.end(); ++it_go ){
  //   //   *(it_st++) = Teuchos::as<size_t>(*it_go);
  //   // }

  //   // l_l_mv_ = l_mv_->subViewNonConst(nodeElements_st);

  //   // return(l_l_mv_->get1dViewNonConst());
  // } else {
  //   scalar_t* values;
  //   int lda;

  //   if( !isLocal() ){
  //     this->localize();
  //     l_mv_->ExtractView(&values, &lda);
  //   } else {
  //     mv_->ExtractView(&values, &lda);
  //   }

  //   TEUCHOS_TEST_FOR_EXCEPTION( lda != Teuchos::as<int>(this->getStride()),
  //     std::logic_error,
  //     "Stride reported during extraction not consistent with what multivector reports");

  //   return Teuchos::arcp(values,0,lda*this->getGlobalNumVectors(),false);
  // }
  return Teuchos::null;
}


void
MultiVecAdapter<Epetra_MultiVector>::put1dData(
  const Teuchos::ArrayView<const MultiVecAdapter<Epetra_MultiVector>::scalar_t>& new_data,
  size_t lda,
  Teuchos::Ptr<
    const Tpetra::Map<MultiVecAdapter<Epetra_MultiVector>::local_ordinal_t,
                      MultiVecAdapter<Epetra_MultiVector>::global_ordinal_t,
                      MultiVecAdapter<Epetra_MultiVector>::node_t> > source_map,
                      EDistribution /* distribution */)
{
  using Teuchos::rcpFromPtr;
  using Teuchos::as;

  const size_t num_vecs  = getGlobalNumVectors();
  // TODO: check that the following const_cast is safe
  double* data_ptr = const_cast<double*>(new_data.getRawPtr());

  // Optimization for ROOTED and single MPI process
  if ( num_vecs == 1 && mv_->Comm().MyPID() == 0 && mv_->Comm().NumProc() == 1 ) {
    // First, functioning impl
    //const multivec_t source_mv(Copy, *mv_map_, data_ptr, as<int>(lda), as<int>(num_vecs));
    //const Epetra_Import importer(*mv_map_, *mv_map_); //trivial - map does not change
	  //mv_->Import(source_mv, importer, Insert);
    // Element-wise copy rather than using importer
    auto vector = mv_->Pointers();
    for ( size_t i = 0; i < lda; ++i ) {
      vector[0][i] = data_ptr[i];
    }
  }
  else {
    const Epetra_BlockMap e_source_map
      = *Util::tpetra_map_to_epetra_map<local_ordinal_t,global_ordinal_t,global_size_t,node_t>(*source_map);
    const multivec_t source_mv(Copy, e_source_map, data_ptr, as<int>(lda), as<int>(num_vecs));
    const Epetra_Import importer(*mv_map_, e_source_map);
  
    mv_->Import(source_mv, importer, Insert);
  }
}

void
MultiVecAdapter<Epetra_MultiVector>::put1dData_kokkos_view_host(
  Kokkos::View<scalar_t**, Kokkos::LayoutLeft, Kokkos::HostSpace> & host_new_data,
  size_t lda,
  Teuchos::Ptr<
    const Tpetra::Map<MultiVecAdapter<Epetra_MultiVector>::local_ordinal_t,
                      MultiVecAdapter<Epetra_MultiVector>::global_ordinal_t,
                      MultiVecAdapter<Epetra_MultiVector>::node_t> > source_map,
                      EDistribution /* distribution */)
{
  using Teuchos::rcpFromPtr;
  using Teuchos::as;

  const size_t num_vecs  = getGlobalNumVectors();

  double* data_ptr = host_new_data.data();

  // Optimization for ROOTED and single MPI process
  if ( num_vecs == 1 && mv_->Comm().MyPID() == 0 && mv_->Comm().NumProc() == 1 ) {
    auto vector = mv_->Pointers();
    for ( size_t i = 0; i < lda; ++i ) {
      vector[0][i] = data_ptr[i];
    }
  }
  else {
    const Epetra_BlockMap e_source_map
      = *Util::tpetra_map_to_epetra_map<local_ordinal_t,global_ordinal_t,global_size_t,node_t>(*source_map);
    const multivec_t source_mv(Copy, e_source_map, data_ptr, as<int>(lda), as<int>(num_vecs));
    const Epetra_Import importer(*mv_map_, e_source_map);

    mv_->Import(source_mv, importer, Insert);
  }
}

std::string MultiVecAdapter<Epetra_MultiVector>::description() const
{
  std::ostringstream oss;
  oss << "Amesos2 adapter wrapping: Epetra_MultiVector";
  return oss.str();
}


void MultiVecAdapter<Epetra_MultiVector>::describe(
  Teuchos::FancyOStream& os,
  const Teuchos::EVerbosityLevel verbLevel) const
{
  // TODO: implement!
  if(verbLevel != Teuchos::VERB_NONE)
    {
      os << "TODO: implement! ";
    }
}


Teuchos::RCP<const Tpetra::Map<MultiVecAdapter<Epetra_MultiVector>::local_ordinal_t,
                               MultiVecAdapter<Epetra_MultiVector>::global_ordinal_t,
                               MultiVecAdapter<Epetra_MultiVector>::node_t> >
MultiVecAdapter<Epetra_MultiVector>::getMap() const
{
  return Util::epetra_map_to_tpetra_map<local_ordinal_t,global_ordinal_t,global_size_t,node_t>(*mv_map_);
}


const char* MultiVecAdapter<Epetra_MultiVector>::name
= "Amesos2 adapter for Epetra_MultiVector";


} // end namespace Amesos2

#endif // AMESOS2_EPETRA_MULTIVEC_ADAPTER_DEF_HPP

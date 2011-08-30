// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
//                  Copyright 2011 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
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


void MultiVecAdapter<Epetra_MultiVector>::get1dCopy(
  const Teuchos::ArrayView<MultiVecAdapter<Epetra_MultiVector>::scalar_t>& av,
  size_t lda,
  Teuchos::Ptr<
    const Tpetra::Map<MultiVecAdapter<Epetra_MultiVector>::local_ordinal_t,
                      MultiVecAdapter<Epetra_MultiVector>::global_ordinal_t,
                      MultiVecAdapter<Epetra_MultiVector>::node_t> > distribution_map ) const
{
  using Teuchos::rcpFromPtr;
  using Teuchos::as;

  const size_t num_vecs = getGlobalNumVectors();
  
#ifdef HAVE_AMESOS2_DEBUG
  const size_t requested_vector_length = distribution_map->getNodeNumElements();
  TEST_FOR_EXCEPTION( lda < requested_vector_length,
		      std::invalid_argument,
		      "Given stride is not large enough for local vector length" );
  TEST_FOR_EXCEPTION( as<size_t>(av.size()) < (num_vecs-1) * lda + requested_vector_length,
		      std::invalid_argument,
		      "MultiVector storage not large enough given leading dimension "
		      "and number of vectors" );
#endif

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


Teuchos::ArrayRCP<MultiVecAdapter<Epetra_MultiVector>::scalar_t>
MultiVecAdapter<Epetra_MultiVector>::get1dViewNonConst(bool local)
{
  // TEST_FOR_EXCEPTION( !this->isConstantStride(),
  //   std::logic_error,
  //   "get1dViewNonConst() : can only get 1d view if stride is constant");

  // if( local ){
  //   TEST_FOR_EXCEPTION(
  //     true,
  //     std::logic_error,
  //     "Amesos2::MultiVecAdapter<Epetra_MultiVector> : 1D views not yet supported for local-local Epetra multi-vectors");
    
  //   // localize();
  //   // /* Use the global element list returned by
  //   //  * mv_->getMap()->getNodeElementList() to get a subCopy of mv_ which we
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

  //   TEST_FOR_EXCEPTION( lda != Teuchos::as<int>(this->getStride()),
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
                      MultiVecAdapter<Epetra_MultiVector>::node_t> > source_map)
{
  using Teuchos::rcpFromPtr;
  using Teuchos::as;

  const size_t num_vecs  = getGlobalNumVectors();
  // TODO: check that the following const_cast is safe
  double* data_ptr = const_cast<double*>(new_data.getRawPtr());
  const Epetra_BlockMap e_source_map
    = *Util::tpetra_map_to_epetra_map<local_ordinal_t,global_ordinal_t,global_size_t,node_t>(*source_map);
  const multivec_t source_mv(Copy, e_source_map, data_ptr, as<int>(lda), as<int>(num_vecs));
  const Epetra_Import importer(*mv_map_, e_source_map);
  
  mv_->Import(source_mv, importer, Insert);
}


std::string MultiVecAdapter<Epetra_MultiVector>::description() const
{
  std::ostringstream oss;
  oss << "Amesos2 adapter wrapping: Epetra_MultiVector";
  return oss.str();
}


void MultiVecAdapter<Epetra_MultiVector>::describe(
  Teuchos::FancyOStream& os,
  const Teuchos::EVerbosityLevel verbLevel=
  Teuchos::Describable::verbLevel_default) const
{
  // TODO: implement!
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

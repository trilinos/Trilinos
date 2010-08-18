// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
//                  Copyright 2010 Sandia Corporation
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

namespace Amesos {

MultiVecAdapter<Epetra_MultiVector>&
MultiVecAdapter<Epetra_MultiVector>::scale(
  const MultiVecAdapter<Epetra_MultiVector>::scalar_type alpha )
{
  mv_->Scale(alpha);

  return *this;
}



MultiVecAdapter<Epetra_MultiVector>&
MultiVecAdapter<Epetra_MultiVector>::update(
  const MultiVecAdapter<Epetra_MultiVector>::scalar_type beta,
  const MultiVecAdapter<Epetra_MultiVector>& B,
  const MultiVecAdapter<Epetra_MultiVector>::scalar_type alpha )
{
  mv_->Update(beta, *(B.mv_), alpha);

  return *this;
}



bool MultiVecAdapter<Epetra_MultiVector>::isLocal() const
{
  return !mv_->DistributedGlobal();
}



// ETB: Following code borrowed from Thyra wrappers
const Teuchos::RCP<const Teuchos::Comm<int> >&
MultiVecAdapter<Epetra_MultiVector>::getComm() const
{
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::set_extra_data;

  RCP<const Epetra_SerialComm>
    serialEpetraComm = rcp_dynamic_cast<const Epetra_SerialComm>(rcp(mv_->Comm().Clone()));
  if( serialEpetraComm.get() ) {
    RCP<const Teuchos::SerialComm<int> >
      serialComm = rcp(new Teuchos::SerialComm<int>());
    set_extra_data( serialEpetraComm, "serialEpetraComm", Teuchos::inOutArg(serialComm) );
    return serialComm;
  }

#ifdef HAVE_MPI
  
  RCP<const Epetra_MpiComm>
    mpiEpetraComm = rcp_dynamic_cast<const Epetra_MpiComm>(rcp(mv_->Comm().Clone()));
  if( mpiEpetraComm.get() ) {
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> >
      rawMpiComm = Teuchos::opaqueWrapper(mpiEpetraComm->Comm());
    set_extra_data( mpiEpetraComm, "mpiEpetraComm", Teuchos::inOutArg(rawMpiComm) );
    RCP<const Teuchos::MpiComm<int> >
      mpiComm = rcp(new Teuchos::MpiComm<int>(rawMpiComm));
    return mpiComm;
  }

#endif // HAVE_MPI
  
  // If you get here then the conversion failed!
  return Teuchos::null;
}


size_t MultiVecAdapter<Epetra_MultiVector>::getLocalLength() const
{
  return Teuchos::as<size_t>(mv_->MyLength());
}


size_t MultiVecAdapter<Epetra_MultiVector>::getLocalNumVectors() const
{
  return Teuchos::as<size_t>(mv_->NumVectors());
}


MultiVecAdapter<Epetra_MultiVector>::global_size_type
MultiVecAdapter<Epetra_MultiVector>::getGlobalLength() const
{
  return Teuchos::as<global_size_type>(mv_->GlobalLength());
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


Teuchos::RCP<const Tpetra::Vector<MultiVecAdapter<Epetra_MultiVector>::scalar_type,
                                  MultiVecAdapter<Epetra_MultiVector>::local_ordinal_type,
                                  MultiVecAdapter<Epetra_MultiVector>::global_ordinal_type,
                                  MultiVecAdapter<Epetra_MultiVector>::node_type> >
MultiVecAdapter<Epetra_MultiVector>::getVector( size_t j ) const
{
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Tpetra::MultiVector;

  typedef scalar_type st;
  typedef local_ordinal_type lot;
  typedef global_ordinal_type got;
  typedef node_type nt;
  
  RCP<MultiVector<st,lot,got,nt> > vector = rcp(new MultiVector<st,lot,got,nt>(this->getMap(),1));

  // Copy vector contents into Tpetra multi-vector
  ArrayRCP<st> it = vector->getDataNonConst(0);
  double* vector_data = mv_->operator[](Teuchos::as<int>(j)); // values from j^th vector
  Tpetra::global_size_t size = vector->getGlobalLength();

  for( Tpetra::global_size_t i = 0; i < size; ++i ){
    *it = vector_data[i];
  }

  return vector->getVector(0);
}


// Implementation is essentially the same as getVector()
Teuchos::RCP<Tpetra::Vector<MultiVecAdapter<Epetra_MultiVector>::scalar_type,
                            MultiVecAdapter<Epetra_MultiVector>::local_ordinal_type,
                            MultiVecAdapter<Epetra_MultiVector>::global_ordinal_type,
                            MultiVecAdapter<Epetra_MultiVector>::node_type> >
MultiVecAdapter<Epetra_MultiVector>::getVectorNonConst( size_t j )
{
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Tpetra::MultiVector;

  typedef scalar_type st;
  typedef local_ordinal_type lot;
  typedef global_ordinal_type got;
  typedef node_type nt;
  
  RCP<MultiVector<st,lot,got,nt> > vector = rcp(new MultiVector<st,lot,got,nt>(this->getMap(),1));

  // Copy vector contents into Tpetra multi-vector
  ArrayRCP<st> it = vector->getDataNonConst(0);
  double* vector_data = mv_->operator[](Teuchos::as<int>(j)); // values from j^th vector
  Tpetra::global_size_t size = vector->getGlobalLength();

  for( Tpetra::global_size_t i = 0; i < size; ++i ){
    *it = vector_data[i];
  }

  return vector->getVectorNonConst(0);
}


void MultiVecAdapter<Epetra_MultiVector>::get1dCopy(
  const Teuchos::ArrayView<MultiVecAdapter<Epetra_MultiVector>::scalar_type>& A,
  size_t lda)
{
  if( !isLocal() ){
    this->localize();

    l_mv_->ExtractCopy(A.getRawPtr(),lda);
  } else {
    mv_->ExtractCopy(A.getRawPtr(),lda);
  }
}


Teuchos::ArrayRCP<MultiVecAdapter<Epetra_MultiVector>::scalar_type>
MultiVecAdapter<Epetra_MultiVector>::get1dViewNonConst(bool local)
{
  TEST_FOR_EXCEPTION( !this->isConstantStride(),
    std::logic_error,
    "get1dViewNonConst() : can only get 1d view if stride is constant");

  if( local ){
    TEST_FOR_EXCEPTION(
      true,
      std::logic_error,
      "Amesos::MultiVecAdapter<Epetra_MultiVector> : 1D views not yet supported for local-local Epetra multi-vectors");
    
    // localize();
    // /* Use the global element list returned by
    //  * mv_->getMap()->getNodeElementList() to get a subCopy of mv_ which we
    //  * assign to l_l_mv_, then return get1dViewNonConst() of l_l_mv_
    //  */
    // l_l_mv_ = Teuchos::null;
    
    // Teuchos::Array<GlobalOrdinal> nodeElements_go(mv_->Map().NumMyElements());
    // mv_->Map().MyGlobalElements(nodeElements_go.getRawPtr());
    // Teuchos::Array<size_t> nodeElements_st(nodeElements_go.size());

    // // Convert the node element to a list of size_t type objects
    // typename Teuchos::ArrayView<const GlobalOrdinal>::iterator it_go;
    // Teuchos::Array<size_t>::iterator it_st = nodeElements_st.begin();
    // for( it_go = nodeElements_go.begin(); it_go != nodeElements_go.end(); ++it_go ){
    //   *(it_st++) = Teuchos::as<size_t>(*it_go);
    // }

    // l_l_mv_ = l_mv_->subViewNonConst(nodeElements_st);

    // return(l_l_mv_->get1dViewNonConst());
  } else {
    scalar_type* values;
    int lda;

    if( !isLocal() ){
      this->localize();
      l_mv_->ExtractView(&values, &lda);
    } else {
      mv_->ExtractView(&values, &lda);
    }

    TEST_FOR_EXCEPTION( lda != Teuchos::as<int>(this->getStride()),
      std::logic_error,
      "Stride reported during extraction not consistent with what multivector reports");

    return Teuchos::arcp(values,0,lda*this->getGlobalNumVectors(),false);
  }
}


void MultiVecAdapter<Epetra_MultiVector>::get2dCopy(
  Teuchos::ArrayView<const Teuchos::ArrayView<MultiVecAdapter<Epetra_MultiVector>::scalar_type> > A)
{
  using Teuchos::ArrayView;
  typedef ArrayView<const ArrayView<scalar_type> >::iterator iterator;
  typedef ArrayView<scalar_type>::iterator inner_iterator;

  TEST_FOR_EXCEPTION( Teuchos::as<size_t>(A.size()) != getGlobalNumVectors(),
    std::length_error,
    "get2dCopy() : The size of A is not equal to the global number of vectors");
  TEST_FOR_EXCEPTION( Teuchos::as<size_t>(A[0].size()) != getGlobalLength(),
    std::length_error,
    "get2dCopy() : The size of arrays in A is not equal to the global vector length");
  
  scalar_type** ptr_ptr = NULL;
  TEST_FOR_EXCEPTION( mv_->ExtractCopy(ptr_ptr) != 0,
    std::runtime_error,
    "get2dCopy() : Error extracting vector values from multi-vector");

  int count = 0, inner_count = 0;
  iterator it = A.begin(), end = A.end();

  for( ; it != end; ++it ){
    inner_iterator inner_it = (*it).begin(), inner_end = (*it).end();
    for( ; inner_it != inner_end; inner_it.operator++() ){
      *inner_it = ptr_ptr[count][inner_count++];
    }
    ++count;
  }
}


Teuchos::ArrayRCP<Teuchos::ArrayRCP<MultiVecAdapter<Epetra_MultiVector>::scalar_type> >
MultiVecAdapter<Epetra_MultiVector>::get2dViewNonConst( bool local )
{
  using Teuchos::ArrayRCP;
  using Teuchos::Array;
  
  TEST_FOR_EXCEPTION( !this->isConstantStride(),
    std::logic_error,
    "get1dViewNonConst() : can only get 1d view if stride is constant");

  if( local ){
    TEST_FOR_EXCEPTION(
      true,
      std::logic_error,
      "Amesos::MultiVecAdapter<Epetra_MultiVector> : 1D views not yet supported for local-local Epetra multi-vectors");
    
    // localize();
    // /* Use the global element list returned by
    //  * mv_->getMap()->getNodeElementList() to get a subCopy of mv_ which we
    //  * assign to l_l_mv_, then return get1dViewNonConst() of l_l_mv_
    //  */
    // l_l_mv_ = Teuchos::null;
    
    // Teuchos::Array<GlobalOrdinal> nodeElements_go(mv_->Map().NumMyElements());
    // mv_->Map().MyGlobalElements(nodeElements_go.getRawPtr());
    // Teuchos::Array<size_t> nodeElements_st(nodeElements_go.size());

    // // Convert the node element to a list of size_t type objects
    // typename Teuchos::ArrayView<const GlobalOrdinal>::iterator it_go;
    // Teuchos::Array<size_t>::iterator it_st = nodeElements_st.begin();
    // for( it_go = nodeElements_go.begin(); it_go != nodeElements_go.end(); ++it_go ){
    //   *(it_st++) = Teuchos::as<size_t>(*it_go);
    // }

    // l_l_mv_ = l_mv_->subViewNonConst(nodeElements_st);

    // return(l_l_mv_->get1dViewNonConst());
  } else {
    scalar_type** vec_values;

    if( !isLocal() ){
      this->localize();
      l_mv_->ExtractView(&vec_values);
    } else {
      mv_->ExtractView(&vec_values);
    }

    size_t numVectors = this->getGlobalNumVectors(), vec_len = this->getGlobalLength();
    Array<ArrayRCP<scalar_type> > ret(numVectors);
    size_t i;

    for( i = 0; i < numVectors; ++i ){
      ret[i] = Teuchos::arcp(vec_values[Teuchos::as<int>(i)], 0, vec_len, false);
    }
    
    return Teuchos::arcpFromArray(ret);
  }
}


void
MultiVecAdapter<Epetra_MultiVector>::globalize()
{
  TEST_FOR_EXCEPTION(
    true,
    std::logic_error,
    "MultiVecAdapter<Epetra_MultiVector> : arbitrary globalize not yet supported");
}


template <typename Value_t>
void
MultiVecAdapter<Epetra_MultiVector>::globalize(const Teuchos::ArrayView<Value_t>& newVals)
{
  if( !isLocal() ){
    if( l_mv_.is_null() ){
      localize();
    }
    Teuchos::ArrayRCP<scalar_type> l_ptr = this->get1dViewNonConst();

    typename Teuchos::ArrayRCP<scalar_type>::iterator it = l_ptr.begin();
    typename Teuchos::ArrayRCP<scalar_type>::iterator end = l_ptr.end();

    typename Teuchos::ArrayView<Value_t>::iterator val_it = newVals.begin();

    // Update local view, doing conversion
    for( ; it != end; ++it ){
      *it = Teuchos::as<scalar_type>(*val_it++);
    }

    mv_->Import(*l_mv_, *exporter_, Insert);
  } else {
    Teuchos::ArrayRCP<scalar_type> ptr = this->get1dViewNonConst();

    typename Teuchos::ArrayRCP<scalar_type>::iterator it = ptr.begin();
    typename Teuchos::ArrayRCP<scalar_type>::iterator end = ptr.end();

    typename Teuchos::ArrayView<Value_t>::iterator val_it = newVals.begin();

    // Update view, doing conversion
    for( ; it != end; ++it ){
      *it = Teuchos::as<scalar_type>(*val_it++);
    }
  }
}


// TODO: adapt to needs
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


void MultiVecAdapter<Epetra_MultiVector>::localize()
{
  using Teuchos::rcp;
  
  if( l_mv_.is_null() ){
    // create local multi-vector, maps, and do import
    *o_map_ = mv_->Map();
    Epetra_Comm* comm = mv_->Comm().Clone();

    size_t numVectors = getGlobalNumVectors();

    l_map_ = rcp(new Epetra_LocalMap(Teuchos::as<int>(numVectors),0,*comm));

    l_mv_ = rcp(new Epetra_MultiVector(*l_map_,Teuchos::as<int>(numVectors)));

    importer_ = rcp(new Epetra_Import(*l_map_, *o_map_));
    exporter_ = rcp(new Epetra_Export(*l_map_, *o_map_));

    l_mv_->Import(*mv_, *importer_, Insert);
  } else {
    // Just update local values
    l_mv_->Import(*mv_, *importer_, Insert);
  }
}


Teuchos::RCP<const Tpetra::Map<MultiVecAdapter<Epetra_MultiVector>::local_ordinal_type,
                               MultiVecAdapter<Epetra_MultiVector>::global_ordinal_type,
                               MultiVecAdapter<Epetra_MultiVector>::node_type> >
MultiVecAdapter<Epetra_MultiVector>::getMap() const
{
  Teuchos::Array<int> local_element_list(mv_->Map().NumMyElements());
  // Get element list with mv_->MyGlobalElements(), then use that list, the
  // number of global elements, and the index base to create a Tpetra::Map
  mv_->Map().MyGlobalElements( local_element_list.getRawPtr() );
  int num_global_elements = mv_->Map().NumGlobalElements();
  int index_base = mv_->Map().IndexBase();

  return( rcp(new Tpetra::Map<local_ordinal_type,global_ordinal_type>(
        num_global_elements,
        local_element_list,
        index_base,
        this->getComm())) );
}


const char* MultiVecAdapter<Epetra_MultiVector>::name
= "Amesos2 adapter for Epetra_MultiVector";


} // end namespace Amesos

#endif // AMESOS2_EPETRA_MULTIVEC_ADAPTER_DEF_HPP

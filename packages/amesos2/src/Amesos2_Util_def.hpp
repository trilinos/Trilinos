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
  \file   Amesos2_Util_def.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Wed May 26 19:37:37 CDT 2010

  \brief  Utility functions for Amesos2
*/

#ifndef AMESOS2_UTIL_DEF_HPP
#define AMESOS2_UTIL_DEF_HPP

#include "Amesos2_config.h"

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>

#include <Teuchos_DefaultSerialComm.hpp>
#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#endif

#include <Teuchos_ScalarTraits.hpp>

#ifdef HAVE_AMESOS2_EPETRA
template <typename LO, typename GO, typename GS, typename Node>
Teuchos::RCP<Tpetra::Map<LO,GO,Node> >
Amesos::Util::epetra_map_to_tpetra_map(const Epetra_BlockMap& map){
  int num_my_elements = map.NumMyElements();
  Teuchos::Array<int> my_global_elements(num_my_elements);
  map.MyGlobalElements(my_global_elements.getRawPtr());

  using Teuchos::as;
  using Teuchos::rcp;
  typedef Tpetra::Map<LO,GO,Node> map_t;
  RCP<map_t> tmap = rcp(new map_t(Teuchos::OrdinalTraits<GS>::invalid(),
				  my_global_elements(),
				  as<GO>(map.IndexBase()),
				  to_teuchos_comm(Teuchos::rcpFromRef(map.Comm()))));
  return tmap;
}

template <typename LO, typename GO, typename GS, typename Node>
Teuchos::RCP<Epetra_Map>
Amesos::Util::tpetra_map_to_epetra_map(const Tpetra::Map<LO,GO,Node>& map){
  using Teuchos::as;
  Teuchos::Array<GO> elements_tmp;
  elements_tmp = map.getNodeElementList();
  int num_my_elements = elements_tmp.size();
  Teuchos::Array<int> my_global_elements(num_my_elements);
  for (int i = 0; i < num_my_elements; ++i){
    my_global_elements[i] = as<int>(elements_tmp[i]);
  }

  using Teuchos::rcp;
  RCP<Epetra_Map> emap = rcp(new Epetra_Map(-1,
					    num_my_elements,
					    my_global_elements.getRawPtr(),
					    as<GO>(map.getIndexBase()),
					    *to_epetra_comm(map.getComm())));
  return emap;
}

const Teuchos::RCP<const Teuchos::Comm<int> >
Amesos::Util::to_teuchos_comm(Teuchos::RCP<const Epetra_Comm> c){
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::set_extra_data;

  Teuchos::RCP<const Epetra_SerialComm>
    serialEpetraComm = rcp_dynamic_cast<const Epetra_SerialComm>(c);
  if( serialEpetraComm.get() ) {
    Teuchos::RCP<const Teuchos::SerialComm<int> >
      serialComm = rcp(new Teuchos::SerialComm<int>());
    set_extra_data( serialEpetraComm, "serialEpetraComm", Teuchos::inOutArg(serialComm) );
    return serialComm;
  }
  
#ifdef HAVE_MPI
  Teuchos::RCP<const Epetra_MpiComm>
    mpiEpetraComm = rcp_dynamic_cast<const Epetra_MpiComm>(c);
  if( mpiEpetraComm.get() ) {
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> >
      rawMpiComm = Teuchos::opaqueWrapper(mpiEpetraComm->Comm());
    set_extra_data( mpiEpetraComm, "mpiEpetraComm", Teuchos::inOutArg(rawMpiComm) );
    Teuchos::RCP<const Teuchos::MpiComm<int> >
      mpiComm = rcp(new Teuchos::MpiComm<int>(rawMpiComm));
    return mpiComm;
  }
#endif

  return(Teuchos::null);
}

const Teuchos::RCP<const Epetra_Comm>
Amesos::Util::to_epetra_comm(Teuchos::RCP<const Teuchos::Comm<int> > c){
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::set_extra_data;

  Teuchos::RCP<const Teuchos::SerialComm<int> >
    serialTeuchosComm = rcp_dynamic_cast<const Teuchos::SerialComm<int> >(c);
  if( serialTeuchosComm.get() ){
    Teuchos::RCP<const Epetra_SerialComm> serialComm = rcp(new Epetra_SerialComm());
    set_extra_data( serialTeuchosComm, "serialTeuchosComm", Teuchos::inOutArg(serialComm) );
    return serialComm;
  }

#ifdef HAVE_MPI
  Teuchos::RCP<const Teuchos::MpiComm<int> >
    mpiTeuchosComm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(c);
  if( mpiTeuchosComm.get() ){
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> >
      rawMpiComm = mpiTeuchosComm->getRawMpiComm();
    set_extra_data( mpiTeuchosComm, "mpiTeuchosComm", Teuchos::inOutArg(rawMpiComm) );
    Teuchos::RCP<const Epetra_MpiComm>
      mpiComm = rcp(new Epetra_MpiComm(*rawMpiComm()));
    return mpiComm;
  }
#endif

  return Teuchos::null;
}

#endif	// HAVE_AMESOS2_EPETRA


template <typename Scalar,
	  typename GlobalOrdinal,
	  typename GlobalSizeT>
void Amesos::Util::transpose(Teuchos::ArrayView<Scalar> vals,
	       Teuchos::ArrayView<GlobalOrdinal> indices,
	       Teuchos::ArrayView<GlobalSizeT> ptr,
	       Teuchos::ArrayView<Scalar> trans_vals,
	       Teuchos::ArrayView<GlobalOrdinal> trans_indices,
	       Teuchos::ArrayView<GlobalSizeT> trans_ptr)
{
  /* We have a compressed-row storage format of this matrix.  We
   * transform this into a compressed-column format using a
   * distribution-counting sort algorithm, which is described by
   * D.E. Knuth in TAOCP Vol 3, 2nd ed pg 78.
   */

#ifdef HAVE_AMESOS2_DEBUG
  typename Teuchos::ArrayView<GlobalOrdinal>::iterator ind_it, ind_begin, ind_end;
  ind_begin = indices.begin();
  ind_end = indices.end();
  size_t min_trans_ptr_size = *std::max_element(ind_begin, ind_end) + 1;
  TEST_FOR_EXCEPTION( trans_ptr.size() < min_trans_ptr_size,
		      std::invalid_argument,
		      "Transpose pointer size not large enough." );
  TEST_FOR_EXCEPTION( trans_vals.size() < vals.size(),
		      std::invalid_argument,
		      "Transpose values array not large enough." );
  TEST_FOR_EXCEPTION( trans_indices.size() < indices.size(),
		      std::invalid_argument,
		      "Transpose indices array not large enough." );
#else
  typename Teuchos::ArrayView<GlobalOrdinal>::iterator ind_it, ind_end;
#endif
  
  // Count the number of entries in each column
  Teuchos::Array<GlobalSizeT> count(trans_ptr.size(), 0);
  ind_end = indices.end();
  for( ind_it = indices.begin(); ind_it != ind_end; ++ind_it ){
    ++(count[(*ind_it) + 1]);
  }

  // Accumulate
  typename Teuchos::Array<GlobalSizeT>::iterator cnt_it, cnt_end;
  cnt_end = count.end();
  for( cnt_it = count.begin() + 1; cnt_it != cnt_end; ++cnt_it ){
    *cnt_it = *cnt_it + *(cnt_it - 1);
  }
  // This becomes the array of column pointers
  trans_ptr.assign(count);

  /* Move the nonzero values into their final place in nzval, based on the
   * counts found previously.
   *
   * This sequence deviates from Knuth's algorithm a bit, following more
   * closely the description presented in Gustavson, Fred G. "Two Fast
   * Algorithms for Sparse Matrices: Multiplication and Permuted
   * Transposition" ACM Trans. Math. Softw. volume 4, number 3, 1978, pages
   * 250--269, http://doi.acm.org/10.1145/355791.355796.
   *
   * The output indices end up in sorted order
   */

  GlobalSizeT size = ptr.size();
  for( GlobalSizeT i = 0; i < size - 1; ++i ){
    GlobalOrdinal u = ptr[i];
    GlobalOrdinal v = ptr[i + 1];
    for( GlobalOrdinal j = u; j < v; ++j ){
      GlobalOrdinal k = count[indices[j]];
      trans_vals[k] = vals[j];
      trans_indices[k] = i;
      ++(count[indices[j]]);
    }
  }
}


template <typename Scalar1, typename Scalar2>
void
Amesos::Util::scale(Teuchos::ArrayView<Scalar1> vals, size_t l,
		    size_t ld, Teuchos::ArrayView<Scalar2> s){
  size_t vals_size = vals.size();
#ifdef HAVE_AMESOS2_DEBUG
  size_t s_size = s.size();
  TEST_FOR_EXCEPTION( s_size < l,
		      std::invalid_argument,
		      "Scale vector must have length at least that of the vector" );
#endif
  size_t i, s_i;
  for( i = 0, s_i = 0; i < vals_size; ++i, ++s_i ){
    if( s_i == l ){
      // bring i to the next multiple of ld
      i += ld - s_i;
      s_i = 0;
    }
    vals[i] *= s[s_i];
  }
}

template <typename Scalar1, typename Scalar2, class BinaryOp>
void
Amesos::Util::scale(Teuchos::ArrayView<Scalar1> vals, size_t l,
		    size_t ld, Teuchos::ArrayView<Scalar2> s,
		    BinaryOp binary_op){
  size_t vals_size = vals.size();
#ifdef HAVE_AMESOS2_DEBUG
  size_t s_size = s.size();
  TEST_FOR_EXCEPTION( s_size < l,
		      std::invalid_argument,
		      "Scale vector must have length at least that of the vector" );
#endif
  size_t i, s_i;
  for( i = 0, s_i = 0; i < vals_size; ++i, ++s_i ){
    if( s_i == l ){
      // bring i to the next multiple of ld
      i += ld - s_i;
      s_i = 0;
    }
    vals[i] = binary_op(vals[i], s[s_i]);
  }
}


/*
 * TODO: Use Matrix and MultiVecAdapters instead of straight-up matrix and
 * vector arguments
 */
template <typename Matrix,
          typename Vector>
void Amesos::Util::computeTrueResidual(
  const Teuchos::RCP<Matrix>& A,
  const Teuchos::RCP<Vector>& X,
  const Teuchos::RCP<Vector>& B,
  const Teuchos::ETransp trans=Teuchos::NO_TRANS,
  const std::string prefix="")
{
  // typename Teuchos::ScalarTraits<typename Matrix::scalar_type>::magnitudeType norm;
  // Tpetra::Vector<typename Matrix::scalar_type> Ax(B.getMap());
  // size_t numVectors = X.getNumVectors();

  // for( int i = 0; i < numVectors; ++i ){
  //   // TODO: Tpetra::Operator::apply only accepts
  //   // Tpetra::MultiVector objects, and indeed, only those that are
  //   // templated on the same types.  We may have to limit the
  //   // objects to Tpetra::RowMatrix and Tpetra::MultiVector
  //   // arguments.
  //   A.apply(*X.getVector(i), Ax, trans);
  //   Ax.update(1.0,*B.getVector(i),-1.0);
  //   norm = Ax.norm2();

  //   if (A.getComm().getRank() == 0){
  //     std::cout << prefix << " : vector "
  //               << i << ", ||b - Ax|| = "
  //               << norm << std::endl;
  //   }
  // }
}


/* We assume that Matrix and Vector are some instance of a
 * Amesos::MatrixAdapter or a Amesos::MultiVecAdapter, or at least implement
 * the required methods
 */
template< typename Matrix,
          typename Vector>
void Amesos::Util::computeVectorNorms(
  const Teuchos::RCP<Matrix> X,
  const Teuchos::RCP<Vector> B,
  std::string prefix="")
{
  typename Matrix::scalar_type normLHS, normRHS;
  size_t numVectors = X->getNumVectors();

  for (int i=0; i<numVectors; ++i){
    normLHS = X->getVector(i)->norm2();
    normRHS = B->getVector(i)->norm2();
    if (X->getMap()->getComm()->getRank() == 0){
      std::cout << prefix << " : vector "
                << ", ||x|| = " << normLHS
                << ", ||b|| = " << normRHS
                << std::endl;
    }
  }
}

template< typename Matrix>
void Amesos::Util::setMaxProcesses(
  const Teuchos::RCP<Matrix>& A,
  int& maxProcesses)
{
  int maxProcs = A->getComm()->getSize();

  switch(maxProcesses){
    case -3:
      maxProcesses = maxProcs; break;
    case -2:
      maxProcesses = (int) sqrt((double)maxProcs); break;
    case -1:			// We should do some testing on this
      // heuristic
      maxProcesses =
        1 + TEUCHOS_MAX(A.getGlobalNumRows() / 10000,
          A.getGlobalNumEntries() / 1000000);
      break;
  }

  if(maxProcesses <= 0) maxProcesses = 1;
  if(maxProcesses > maxProcs) maxProcesses = maxProcs;

  return;
}

/// Prints a line of 80 "-"s on std::cout.
void Amesos::Util::printLine( Teuchos::FancyOStream &out )
{
  out << "----------------------------------------"
      << "----------------------------------------"
      << std::endl;
}


#endif	// #ifndef AMESOS2_UTIL_DEF_HPP

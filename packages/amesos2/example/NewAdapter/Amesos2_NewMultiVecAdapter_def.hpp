// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
  \file   Amesos2_NewMultiVecAdapter_def.hpp
  \author John Joe <jd@sandia.gov>
  \date   

  \brief  Amesos2::MultiVecAdapter specialization for the
          NewMultiVector class.
*/

#ifndef AMESOS2_NEWMULTIVEC_ADAPTER_DEF_HPP
#define AMESOS2_NEWMULTIVEC_ADAPTER_DEF_HPP


namespace Amesos {


MultiVecAdapter<multivec_type>&
MultiVecAdapter<NewMultiVector>::scale( const scalar_type alpha )
{
  // TODO: implement!
}



MultiVecAdapter<multivec_type>&
MultiVecAdapter<NewMultiVector>::update(
  const scalar_type beta,
  const MultiVecAdapter<multivec_type>& B,
  const scalar_type alpha )
{
  // TODO: implement!
}



bool MultiVecAdapter<NewMultiVector>::isLocal() const
{
  // TODO: implement!
}



const Teuchos::RCP<const Teuchos::Comm<int> >&
MultiVecAdapter<NewMultiVector>::getComm() const
{
  // TODO: implement!
}


size_t MultiVecAdapter<NewMultiVector>::getLocalLength() const
{
  // TODO: implement!
}


size_t MultiVecAdapter<NewMultiVector>::getLocalNumVectors() const
{
  // TODO: implement!
}


global_size_type MultiVecAdapter<NewMultiVector>::getGlobalLength() const
{
  // TODO: implement!
}


size_t MultiVecAdapter<NewMultiVector>::getGlobalNumVectors() const
{
  // TODO: implement!
}


size_t MultiVecAdapter<NewMultiVector>::getStride() const
{
  // TODO: implement!
}


bool MultiVecAdapter<NewMultiVector>::isConstantStride() const
{
  // TODO: implement!
}


Teuchos::RCP<const Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >
MultiVecAdapter<NewMultiVector>::getVector( size_t j ) const
{
  // TODO: implement!
}


Teuchos::RCP<Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >
MultiVecAdapter<NewMultiVector>::getVectorNonConst( size_t j )
{
  // TODO: implement!
}


void MultiVecAdapter<NewMultiVector>::get1dCopy(
  const Teuchos::ArrayView<scalar_type>& A,
  size_t lda) const
{
  // TODO: implement!
}


Teuchos::ArrayRCP<scalar_type>
MultiVecAdapter<NewMultiVector>::get1dViewNonConst(bool local)
{
  // TODO: implement!
}


void MultiVecAdapter<NewMultiVector>::get2dCopy(
  Teuchos::ArrayView<const Teuchos::ArrayView<scalar_type> > A) const
{
  // TODO: implement!
}


Teuchos::ArrayRCP<Teuchos::ArrayRCP<scalar_type> >
MultiVecAdapter<NewMultiVector>::get2dViewNonConst( bool local ) const
{
  // TODO: implement!
}


void
MultiVecAdapter<NewMultiVector>::globalize()
{
  // TODO: implement!
}


template <typename Value_t>
void
MultiVecAdapter<NewMultiVector>::globalize(const Teuchos::ArrayView<Value_t>& newVals)
{
  // TODO: implement!
}


// TODO: adapt to needs
std::string MultiVecAdapter<NewMultiVector>::description() const
{
  std::ostringstream oss;
  oss << "Amesos2 adapter wrapping: NewMultiVector";
  return oss.str();
}


void MultiVecAdapter<NewMultiVector>::describe(
  Teuchos::FancyOStream& os,
  const Teuchos::EVerbosityLevel verbLevel=
  Teuchos::Describable::verbLevel_default) const
{
  // TODO: implement!
}


void MultiVecAdapter<NewMultiVector>::localize()
{
  // TODO: implement!
}


const char* MultiVecAdapter<NewMultiVector>::name = "Amesos2 adapter for NewMultiVector";


} // end namespace Amesos

#endif // AMESOS2_NEWMULTIVEC_ADAPTER_DEF_HPP

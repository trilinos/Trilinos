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

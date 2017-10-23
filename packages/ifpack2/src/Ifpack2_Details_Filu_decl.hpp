/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
//@HEADER
*/

/// @file Ifpack2_filu_decl.hpp

#ifndef __IFPACK2_FILU_DECL_HPP__ 
#define __IFPACK2_FILU_DECL_HPP__ 

#include <Ifpack2_Details_FastILU_Base.hpp>

//forward-declare the local preconditioner type
template<typename LocalOrdinal, typename Scalar, typename execution_space>
class FastILUPrec;

namespace Ifpack2
{
namespace Details
{

/// \class Filu
/// \brief The Ifpack2 wrapper to the ILU preconditioner of ShyLU FastILU.
template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
class Filu : public FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>
{
  public:
    typedef FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node> Base;
    typedef typename Base::TRowMatrix TRowMatrix;
    typedef typename Base::ScalarArray ScalarArray;
    typedef FastILUPrec<LocalOrdinal, Scalar, typename Base::execution_space> LocalFILU;

    //! Constructor
    Filu(Teuchos::RCP<const TRowMatrix> mat_);

    //! Get the sweeps (\"nFact\") from localPrec_
    int getSweeps() const;

    //! Get the number of triangular solves (\"nTrisol\") from localPrec_
    int getNTrisol() const;

    //! Verify and print debug info about the internal ILU preconditioner
    void checkLocalILU() const;

    //! Verify and print debug info about the internal IC preconditioner
    void checkLocalIC() const;

  protected:
    mutable Teuchos::RCP<LocalFILU> localPrec_;

    void initLocalPrec();
    //compute() takes A's local values
    void computeLocalPrec();
    void applyLocalPrec(ScalarArray x, ScalarArray y) const;
    std::string getName() const;
};

} //namespace Details
} //namespace Ifpack2

#endif


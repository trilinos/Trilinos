#ifndef INTREPID2_BASISSET_HPP
#define INTREPID2_BASISSET_HPP

// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_BasisSet.hpp
    \brief  Provide a family of basis to construct cell basis and gluing
    inter-element base
    \author Created by Kyungjoo Kim
*/

#if defined( INTREPID_USING_EXPERIMENTAL_HIGH_ORDER )

#include "Intrepid2_FieldContainer.hpp"
#include "Intrepid2_FieldContainer_Kokkos.hpp"
#include "Intrepid2_RealSpaceTools.hpp"

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"
#include "Intrepid2_Basis.hpp"

#include "Shards_CellTopology.hpp"
#include "Shards_BasicTopologies.hpp"

#include "Intrepid2_HGRAD_LINE_Cn_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_Cn_FEM.hpp"
#include "Intrepid2_HGRAD_TET_Cn_FEM.hpp"

#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"

#include <Intrepid2_KokkosRank.hpp>
#include "Kokkos_Core.hpp"

namespace Intrepid2 {

  // Basis super class cannot be instanciated as it includes pure virtual functions.
  template<class Scalar, class ArrayType>
  class Basis_Dummy_FEM: public Basis<Scalar, ArrayType> {
  private:
    virtual void initializeTags();
  public:
    Basis_Dummy_FEM() = default;
    void getValues(ArrayType &       outputValues,
                   const ArrayType & inputPoints,
                   const EOperator   operatorType) const;

    void getValues(ArrayType &       outputValues,
                   const ArrayType & inputPoints,
                   const ArrayType & cellVertices,
                   const EOperator   operatorType) const;

  };

  template<class Scalar, class ArrayType>
  class BasisSet {
  private:
    EFunctionSpace _space;
    Basis_Dummy_FEM<Scalar,ArrayType> _dummy;

  public:
    BasisSet() = default;
    BasisSet(const EFunctionSpace space);
    EFunctionSpace getFunctionSpace() const;

    // these get functions return exceptions and should be over-ridden properly
    virtual const Basis<Scalar,ArrayType>& getCellBasis()          const;
    virtual const Basis<Scalar,ArrayType>& getTriangleBasis()      const;
    virtual const Basis<Scalar,ArrayType>& getQuadrilateralBasis() const;
    virtual const Basis<Scalar,ArrayType>& getLineBasis()          const;
  };

  template<class Scalar, class ArrayType>
  class BasisSet_HGRAD_TRI_Cn_FEM : public BasisSet<Scalar,ArrayType> {
  private:
    Basis_HGRAD_TRI_Cn_FEM<Scalar,ArrayType>  _cellBasis;
    Basis_HGRAD_LINE_Cn_FEM<Scalar,ArrayType> _lineBasis;
  public:
    BasisSet_HGRAD_TRI_Cn_FEM() = delete;
    BasisSet_HGRAD_TRI_Cn_FEM( const int n , const EPointType pointType );
    const Basis<Scalar,ArrayType>& getCellBasis() const;
    const Basis<Scalar,ArrayType>& getLineBasis() const;
  };

  template<class Scalar, class ArrayType>
  class BasisSet_HGRAD_TET_Cn_FEM : public BasisSet<Scalar,ArrayType> {
  private:
    Basis_HGRAD_TET_Cn_FEM<Scalar,ArrayType>  _cellBasis;
    Basis_HGRAD_TRI_Cn_FEM<Scalar,ArrayType>  _trigBasis;
    Basis_HGRAD_LINE_Cn_FEM<Scalar,ArrayType> _lineBasis;
  public:
    BasisSet_HGRAD_TET_Cn_FEM() = delete;
    BasisSet_HGRAD_TET_Cn_FEM( const int n , const EPointType pointType );
    const Basis<Scalar,ArrayType>& getCellBasis()     const;
    const Basis<Scalar,ArrayType>& getTriangleBasis() const;
    const Basis<Scalar,ArrayType>& getLineBasis()     const;
  };

}

// include templated function definitions
#include "Intrepid2_BasisSetDef.hpp"

#endif

#endif

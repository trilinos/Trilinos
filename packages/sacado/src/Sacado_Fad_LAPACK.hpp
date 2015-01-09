// $Id$
// $Source$
// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef SACADO_FAD_LAPACK_HPP
#define SACADO_FAD_LAPACK_HPP

#include "Teuchos_LAPACK.hpp"
#include "Sacado.hpp"
#include "Sacado_CacheFad_DFad.hpp"
#include "Sacado_dummy_arg.hpp"

namespace Sacado {

  namespace Fad {

    template <typename OrdinalType, typename FadType>
    class ArrayTraits {

      typedef typename Sacado::ValueType<FadType>::type ValueType;
      typedef typename Sacado::ScalarType<FadType>::type scalar_type;
      typedef typename Sacado::dummy<ValueType,scalar_type>::type ScalarType;

    public:

      ArrayTraits(bool use_dynamic = true,
                  OrdinalType workspace_size = 0);

      ArrayTraits(const ArrayTraits& a);

      ~ArrayTraits();

      void unpack() const;

      void pack() const;

      void free() const;

      ValueType* allocate_array(OrdinalType size) const;

      void free_array(const ValueType* ptr, OrdinalType size) const;

      bool is_array_contiguous(const FadType* a, OrdinalType n,
                               OrdinalType n_dot) const;

    protected:

      //! Use dynamic memory allocation
      bool use_dynamic;

      //! Size of static workspace
      OrdinalType workspace_size;

      //! Workspace for holding contiguous values/derivatives
      mutable ValueType *workspace;

      //! Pointer to current free entry in workspace
      mutable ValueType *workspace_pointer;

    };

    template <typename T> struct ArrayValueType { typedef T type; };

    //! Fad specializations for Teuchos::LAPACK wrappers
    template <typename OrdinalType, typename FadType>
    class Fad_LAPACK {

      typedef typename Teuchos::ScalarTraits<FadType>::magnitudeType MagnitudeType;
      typedef typename Sacado::ValueType<FadType>::type ValueType;
      typedef typename Sacado::ScalarType<FadType>::type scalar_type;
      typedef typename Sacado::dummy<ValueType,scalar_type>::type ScalarType;
      typedef Teuchos::LAPACK<OrdinalType,FadType> LAPACKType;

    public:
      //! @name Constructor/Destructor.
      //@{

      //! Default constructor.
      Fad_LAPACK(bool use_default_impl = true,
                 bool use_dynamic = true, 
                 OrdinalType static_workspace_size = 0);

      //! Copy constructor.

      Fad_LAPACK(const Fad_LAPACK& x);

      //! Destructor.
      virtual ~Fad_LAPACK();

      //@}

      //! Computes the solution to a real system of linear equations
      void GESV(const OrdinalType n, const OrdinalType nrhs, FadType* A, const OrdinalType lda,
                OrdinalType* IPIV, FadType* B, const OrdinalType ldb, OrdinalType* info) const; 

    protected:

      //! ArrayTraits for packing/unpacking value/derivative arrays
      ArrayTraits<OrdinalType,FadType> arrayTraits;

      //! LAPACK for values
      Teuchos::LAPACK<OrdinalType, ValueType> lapack;

      //! Use custom or default implementation
      bool use_default_impl;

    protected:

      //! Implementation of GESV
      void Fad_GESV() const;

    }; // class FadLAPACK

  } // namespace Fad

} // namespace Sacado

namespace Teuchos {
  // Specialization of Teuchos::LAPACK for Sacado::Fad::DFad
  template <typename OrdinalType, typename ScalarType>
  class LAPACK< OrdinalType, Sacado::Fad::DFad<ScalarType> > :
    public Sacado::Fad::Fad_LAPACK< OrdinalType, Sacado::Fad::DFad<ScalarType> > {};
}

#include "Sacado_Fad_LAPACKImp.hpp"

#endif // SACADO_FAD_LAPACK_HPP

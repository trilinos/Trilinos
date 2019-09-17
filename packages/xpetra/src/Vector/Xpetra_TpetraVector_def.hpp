// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_TPETRAVECTOR_DEF_HPP
#define XPETRA_TPETRAVECTOR_DEF_HPP
#include "Xpetra_TpetraVector_decl.hpp"


namespace Xpetra {


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
TpetraVector(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
                                bool                                                              zeroOut)
    : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, 1, zeroOut)
{
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
TpetraVector(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
             const Teuchos::ArrayView<const Scalar>&                           A)
    : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, A, map->getNodeNumElements(), 1)
{
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
~TpetraVector()
{
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
replaceGlobalValue(GlobalOrdinal globalRow, const Scalar& value)
{
    XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceGlobalValue");
    getTpetra_Vector()->replaceGlobalValue(globalRow, value);
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar& value)
{
    XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoGlobalValue");
    getTpetra_Vector()->sumIntoGlobalValue(globalRow, value);
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
replaceLocalValue(LocalOrdinal myRow, const Scalar& value)
{
    XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceLocalValue");
    getTpetra_Vector()->replaceLocalValue(myRow, value);
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
sumIntoLocalValue(LocalOrdinal myRow, const Scalar& value)
{
    XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoLocalValue");
    getTpetra_Vector()->sumIntoLocalValue(myRow, value);
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
norm1() const
{
    XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm1");
    return getTpetra_Vector()->norm1();
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
norm2() const
{
    XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm2");
    return getTpetra_Vector()->norm2();
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
normInf() const
{
    XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normInf");
    return getTpetra_Vector()->normInf();
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
meanValue() const
{
    XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::meanValue");
    return getTpetra_Vector()->meanValue();
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
description() const
{
    XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::description");
    return getTpetra_Vector()->description();
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
    XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::describe");
    getTpetra_Vector()->describe(out, verbLevel);
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
dot(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& a) const
{
    XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dot");
    return getTpetra_Vector()->dot(*toTpetra(a));
}


#ifdef TPETRA_ENABLE_DEPRECATED_CODE
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TPETRA_DEPRECATED typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
normWeighted(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& weights) const
{
    XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normWeighted");
    return getTpetra_Vector()->normWeighted(*toTpetra(weights));
}
#endif


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
TpetraVector(const Teuchos::RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& vec)
    : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(vec)
{
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getTpetra_Vector() const
{
    return this->TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getTpetra_MultiVector()->getVectorNonConst(0);
}


#ifdef HAVE_XPETRA_KOKKOS_REFACTOR

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type::t_host_um
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getHostLocalView() const
{
    return this->TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getHostLocalView();
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type::t_dev_um
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getDeviceLocalView() const
{
    return this->TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getDeviceLocalView();
}


#if 0
 template<class TargetDeviceType>
    typename Kokkos::Impl::if_c<
      Kokkos::Impl::is_same<
        typename dual_view_type::t_dev_um::execution_space::memory_space,
        typename TargetDeviceType::memory_space>::value,
        typename dual_view_type::t_dev_um,
        typename dual_view_type::t_host_um>::type
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
 TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalView () const {
      return this->TpetraMultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >::template getLocalView<TargetDeviceType>();
    }
#endif  // #if 0

#endif  // HAVE_XPETRA_KOKKOS_REFACTOR



#ifdef HAVE_XPETRA_EPETRA



#if((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) \
    || (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))

// specialization of TpetraVector for GO=int and NO=SerialNode
template<class Scalar>
class TpetraVector<Scalar, int, int, EpetraNode>
    : public virtual Vector<Scalar, int, int, EpetraNode>
    , public TpetraMultiVector<Scalar, int, int, EpetraNode>
{
    typedef int        LocalOrdinal;
    typedef int        GlobalOrdinal;
    typedef EpetraNode Node;

#undef XPETRA_TPETRAMULTIVECTOR_SHORT
#undef XPETRA_TPETRAVECTOR_SHORT
#include "Xpetra_UseShortNames.hpp"
#define XPETRA_TPETRAMULTIVECTOR_SHORT
#define XPETRA_TPETRAVECTOR_SHORT

  public:

    using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dot;                     // overloading, not hiding
    using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::norm1;                   // overloading, not hiding
    using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::norm2;                   // overloading, not hiding
    using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::normInf;                 // overloading, not hiding
    using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::meanValue;               // overloading, not hiding
    using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceGlobalValue;      // overloading, not hiding
    using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sumIntoGlobalValue;      // overloading, not hiding
    using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceLocalValue;       // overloading, not hiding
    using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sumIntoLocalValue;       // overloading, not hiding

    //! @name Constructor/Destructor Methods
    //@{

    //! Sets all vector entries to zero.
    TpetraVector(const Teuchos::RCP<const Map>& map, bool zeroOut = true)
        : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, 1, zeroOut)
    {
        XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                    typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                    "int",
                                    typeid(EpetraNode).name());
    }


    //! Set multi-vector values from an array using Teuchos memory management classes. (copy)
    TpetraVector(const Teuchos::RCP<const Map>& map, const Teuchos::ArrayView<const Scalar>& A)
        : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, A, map->getNodeNumElements(), 1)
    {
        XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                    typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                    "int",
                                    typeid(EpetraNode).name());
    }


    virtual ~TpetraVector() {}


    void replaceGlobalValue(GlobalOrdinal globalRow, const Scalar& value) {}


    //! Adds specified value to existing value at the specified location.
    void sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar& value) {}


    //! Replace current value at the specified location with specified values.
    void replaceLocalValue(LocalOrdinal myRow, const Scalar& value) {}


    //! Adds specified value to existing value at the specified location.
    void sumIntoLocalValue(LocalOrdinal myRow, const Scalar& value) {}


    //@}

    //! @name Mathematical methods
    //@{


    //! Return 1-norm of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm1() const
    {
        return Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero());
    }


    //! Compute 2-norm of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm2() const
    {
        return Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero());
    }


    //! Compute Inf-norm of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType normInf() const
    {
        return Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero());
    }


    //! Compute mean (average) value of this Vector.
    Scalar meanValue() const { return Teuchos::ScalarTraits<Scalar>::zero(); }

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{


    //! Return a simple one-line description of this object.
    std::string description() const { return std::string(""); }


    //! Print the object with some verbosity level to an FancyOStream object.
    void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const {}

    //@}


    //! Computes dot product of this Vector against input Vector x.
    Scalar dot(const Vector& a) const { return Teuchos::ScalarTraits<Scalar>::zero(); }


#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    //! Compute Weighted 2-norm (RMS Norm) of this Vector.
    TPETRA_DEPRECATED typename Teuchos::ScalarTraits<Scalar>::magnitudeType normWeighted(const Vector& weights) const
    {
        return Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero());
    }
#endif      // TPETRA_ENABLE_DEPRECATED_CODE


    //! @name Xpetra specific
    //@{

    //! TpetraMultiVector constructor to wrap a Tpetra::MultiVector object
    TpetraVector(const Teuchos::RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& vec)
        : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(vec)
    {
        XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                    typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                    "int",
                                    typeid(EpetraNode).name());
    }


    //! Get the underlying Tpetra multivector
    RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
    getTpetra_Vector() const
    {
        return Teuchos::null;
    }



#ifdef HAVE_XPETRA_KOKKOS_REFACTOR

    using dual_view_type = typename Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type dual_view_type;

    typename dual_view_type::t_host_um getHostLocalView() const
    {
        typename dual_view_type::t_host_um ret;
        return ret;
    }

    typename dual_view_type::t_dev_um getDeviceLocalView() const
    {
        typename dual_view_type::t_dev_um ret;
        return ret;
    }



    /// \brief Return an unmanaged non-const view of the local data on a specific device.
    /// \tparam TargetDeviceType The Kokkos Device type whose data to return.
    ///
    /// \warning DO NOT USE THIS FUNCTION! There is no reason why you are working directly
    ///          with the Xpetra::TpetraVector object. To write a code which is independent
    ///          from the underlying linear algebra package you should always use the abstract class,
    ///          i.e. Xpetra::Vector!
    ///
    /// \warning Be aware that the view on the vector data is non-persisting, i.e.
    ///          only valid as long as the vector does not run of scope!
    template<class TargetDeviceType>
    typename Kokkos::Impl::if_c<
      Kokkos::Impl::is_same<typename dual_view_type::t_dev_um::execution_space::memory_space, typename TargetDeviceType::memory_space>::value,
      typename dual_view_type::t_dev_um,
      typename dual_view_type::t_host_um>::type
    getLocalView() const
    {
        typename Kokkos::Impl::if_c<
          Kokkos::Impl::is_same<typename dual_view_type::t_dev_um::execution_space::memory_space, typename TargetDeviceType::memory_space>::value,
          typename dual_view_type::t_dev_um,
          typename dual_view_type::t_host_um>::type ret;
        return ret;
    }

#endif  // HAVE_XPETRA_KOKKOS_REFACTOR


    //@}

};      // TpetraVector class (specialization on GO=int, NO=EpetraNode)


#endif    // #if((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT)))
          //    || (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))



#if((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))) \
    || (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))))

// specialization of TpetraVector for GO=int and NO=SerialNode
template<class Scalar>
class TpetraVector<Scalar, int, long long, EpetraNode>
    : public virtual Vector<Scalar, int, long long, EpetraNode>
    , public TpetraMultiVector<Scalar, int, long long, EpetraNode>
{
    typedef int        LocalOrdinal;
    typedef long long  GlobalOrdinal;
    typedef EpetraNode Node;

#undef XPETRA_TPETRAMULTIVECTOR_SHORT
#undef XPETRA_TPETRAVECTOR_SHORT
#include "Xpetra_UseShortNames.hpp"
#define XPETRA_TPETRAMULTIVECTOR_SHORT
#define XPETRA_TPETRAVECTOR_SHORT

  public:
    using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dot;                     // overloading, not hiding
    using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::norm1;                   // overloading, not hiding
    using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::norm2;                   // overloading, not hiding
    using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::normInf;                 // overloading, not hiding
    using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::meanValue;               // overloading, not hiding
    using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceGlobalValue;      // overloading, not hiding
    using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sumIntoGlobalValue;      // overloading, not hiding
    using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceLocalValue;       // overloading, not hiding
    using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sumIntoLocalValue;       // overloading, not hiding

    //! @name Constructor/Destructor Methods
    //@{

    //! Sets all vector entries to zero.
    TpetraVector(const Teuchos::RCP<const Map>& map, bool zeroOut = true)
        : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, 1, zeroOut)
    {
        XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                    typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                    "long long",
                                    typeid(EpetraNode).name());
    }

    //! Set multi-vector values from an array using Teuchos memory management classes. (copy)
    TpetraVector(const Teuchos::RCP<const Map>& map, const Teuchos::ArrayView<const Scalar>& A)
        : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, A, map->getNodeNumElements(), 1)
    {
        XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                    typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                    "long long",
                                    typeid(EpetraNode).name());
    }

    //! Destructor.
    virtual ~TpetraVector() {}

    //@}

    //! @name Post-construction modification routines
    //@{

    //! Replace current value at the specified location with specified value.
    void replaceGlobalValue(GlobalOrdinal globalRow, const Scalar& value) {}

    //! Adds specified value to existing value at the specified location.
    void sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar& value) {}

    //! Replace current value at the specified location with specified values.
    void replaceLocalValue(LocalOrdinal myRow, const Scalar& value) {}

    //! Adds specified value to existing value at the specified location.
    void sumIntoLocalValue(LocalOrdinal myRow, const Scalar& value) {}

    //@}

    //! @name Mathematical methods
    //@{

    //! Return 1-norm of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm1() const
    {
        return Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero());
    }

    //! Compute 2-norm of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm2() const
    {
        return Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero());
    }

    //! Compute Inf-norm of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType normInf() const
    {
        return Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero());
    }

    //! Compute mean (average) value of this Vector.
    Scalar meanValue() const { return Teuchos::ScalarTraits<Scalar>::zero(); }

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const { return std::string(""); }

    //! Print the object with some verbosity level to an FancyOStream object.
    void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const {}

    //@}

    //! Computes dot product of this Vector against input Vector x.
    Scalar dot(const Vector& a) const { return Teuchos::ScalarTraits<Scalar>::zero(); }

    //! Compute Weighted 2-norm (RMS Norm) of this Vector.
#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    TPETRA_DEPRECATED typename Teuchos::ScalarTraits<Scalar>::magnitudeType normWeighted(const Vector& weights) const
    {
        return Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero());
    }
#endif      // TPETRA_ENABLE_DEPRECATED_CODE

    //! @name Xpetra specific
    //@{


    //! TpetraMultiVector constructor to wrap a Tpetra::MultiVector object
    TpetraVector(const Teuchos::RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& vec)
    {
        XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                    typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                    "long long",
                                    typeid(EpetraNode).name());
    }


    //! Get the underlying Tpetra multivector
    RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
    getTpetra_Vector() const
    {
        return Teuchos::null;
    }


#ifdef HAVE_XPETRA_KOKKOS_REFACTOR

    typedef typename Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type dual_view_type;

    typename dual_view_type::t_host_um getHostLocalView() const
    {
        typename dual_view_type::t_host_um ret;
        return ret;
    }

    typename dual_view_type::t_dev_um getDeviceLocalView() const
    {
        typename dual_view_type::t_dev_um ret;
        return ret;
    }

    /// \brief Return an unmanaged non-const view of the local data on a specific device.
    /// \tparam TargetDeviceType The Kokkos Device type whose data to return.
    ///
    /// \warning DO NOT USE THIS FUNCTION! There is no reason why you are working directly
    ///          with the Xpetra::TpetraVector object. To write a code which is independent
    ///          from the underlying linear algebra package you should always use the abstract class,
    ///          i.e. Xpetra::Vector!
    ///
    /// \warning Be aware that the view on the vector data is non-persisting, i.e.
    ///          only valid as long as the vector does not run of scope!
    template<class TargetDeviceType>
    typename Kokkos::Impl::if_c<
      Kokkos::Impl::is_same<typename dual_view_type::t_dev_um::execution_space::memory_space, typename TargetDeviceType::memory_space>::value,
      typename dual_view_type::t_dev_um,
      typename dual_view_type::t_host_um>::type
    getLocalView() const
    {
        typename Kokkos::Impl::if_c<
          Kokkos::Impl::is_same<typename dual_view_type::t_dev_um::execution_space::memory_space, typename TargetDeviceType::memory_space>::value,
          typename dual_view_type::t_dev_um,
          typename dual_view_type::t_host_um>::type ret;
        return ret;
    }
#endif  // HAVE_XPETRA_KOKKOS_REFACTOR

    //@}

};      // TpetraVector class (specialization on GO=long long, NO=EpetraNode)


#endif   // #if((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG)))
         //    || (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))))



#endif      // HAVE_XPETRA_EPETRA
}      // namespace Xpetra

#endif

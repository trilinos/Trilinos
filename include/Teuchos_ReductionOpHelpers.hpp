// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_REDUCTION_OP_HELPERS_HPP
#define TEUCHOS_REDUCTION_OP_HELPERS_HPP

#include "Teuchos_ReductionOp.hpp"
#include "Teuchos_SerializationTraitsHelpers.hpp"
#include "Teuchos_SerializerHelpers.hpp"

namespace Teuchos {

/** \brief Decorator class that uses traits to convert to and from
 * <tt>char[]</tt> to typed buffers for objects that use value semantics and
 * then call a type-specific reduction object.
 *
 * ToDo: Finish Documentation!
 */
template<typename Ordinal, typename T, typename Serializer>
class CharToValueTypeReductionOpImp : public ValueTypeReductionOp<Ordinal,char>
{
public:
  /** \brief . */
  CharToValueTypeReductionOpImp(
    const RCP<const ValueTypeReductionOp<Ordinal,T> >  &reductOp,
    const RCP<const Serializer>& serializer
    );
  /** \brief . */
  void reduce(
    const Ordinal     charCount
    ,const char       charInBuffer[]
    ,char             charInoutBuffer[]
    ) const;
private:
  RCP<const ValueTypeReductionOp<Ordinal,T> >  reductOp_;
  RCP<const Serializer> serializer_;
  // Not defined and not to be called!
  CharToValueTypeReductionOpImp();
  CharToValueTypeReductionOpImp(const CharToValueTypeReductionOpImp&);
  CharToValueTypeReductionOpImp& operator=(const CharToValueTypeReductionOpImp&);
};

/** \brief Decorator class that uses traits to convert to and from
 * <tt>char[]</tt> to typed buffers for objects that use value semantics and
 * then call a type-specific reduction object.
 *
 * ToDo: Finish Documentation!
 */
template<typename Ordinal, typename T,
	 typename Serializer = typename DefaultSerializer<Ordinal,T>::DefaultSerializerType>
class CharToValueTypeReductionOp :
    public CharToValueTypeReductionOpImp<Ordinal,T,Serializer>
{
public:
  typedef CharToValueTypeReductionOpImp<Ordinal,T,Serializer> Base;
  /** \brief . */
  CharToValueTypeReductionOp(
    const RCP<const ValueTypeReductionOp<Ordinal,T> >  &reductOp,
    const RCP<const Serializer>& serializer
    ) : Base(reductOp, serializer) {}
};

/** \brief Decorator class that uses traits to convert to and from
 * <tt>char[]</tt> to typed buffers for objects that use value semantics and
 * then call a type-specific reduction object.
 *
 * Specialization for the default serializer object type with a default
 * argument for the serializer object parameter.
 *
 * ToDo: Finish Documentation!
 */
template<typename Ordinal, typename T>
class CharToValueTypeReductionOp<Ordinal,T,typename DefaultSerializer<Ordinal,T>::DefaultSerializerType> :
    public CharToValueTypeReductionOpImp<Ordinal,T,typename DefaultSerializer<Ordinal,T>::DefaultSerializerType>
{
public:
  typedef DefaultSerializer<Ordinal,T> DS;  // work around for parsing bug in gcc 4.1-4.2
  typedef typename DS::DefaultSerializerType Serializer;
  typedef CharToValueTypeReductionOpImp<Ordinal,T,Serializer> Base;
  /** \brief . */
  CharToValueTypeReductionOp(
    const RCP<const ValueTypeReductionOp<Ordinal,T> >  &reductOp,
    const RCP<const Serializer>& serializer = DS::getDefaultSerializerRCP()
    ) : Base(reductOp, serializer) {}
};

/** \brief Decorator class that uses a strategy object to convert to and from
 * <tt>char[]</tt> to typed buffers for objects that use reference semantics
 * and then call a type-specific reduction object.
 *
 * ToDo: Finish Documentation!
 */
template<typename Ordinal, typename T>
class CharToReferenceTypeReductionOp : public ValueTypeReductionOp<Ordinal,char>
{
public:
  /** \brief . */
  CharToReferenceTypeReductionOp(
    const RCP<const Serializer<Ordinal,T> >                 &serializer
    ,const RCP<const ReferenceTypeReductionOp<Ordinal,T> >  &reductOp
    );
  /** \brief . */
  void reduce(
    const Ordinal     charCount
    ,const char       charInBuffer[]
    ,char             charInoutBuffer[]
    ) const;
private:
  RCP<const Serializer<Ordinal,T> >                serializer_;
  RCP<const ReferenceTypeReductionOp<Ordinal,T> >  reductOp_;
  // Not defined and not to be called!
  CharToReferenceTypeReductionOp();
  CharToReferenceTypeReductionOp(const CharToReferenceTypeReductionOp&);
  CharToReferenceTypeReductionOp& operator=(const CharToReferenceTypeReductionOp&);
};

// /////////////////////////////////////
// Template implementations

//
// CharToValueTypeReductionOpImp
//

template<typename Ordinal, typename T, typename Serializer>
CharToValueTypeReductionOpImp<Ordinal,T,Serializer>::CharToValueTypeReductionOpImp(
  const RCP<const ValueTypeReductionOp<Ordinal,T> >  &reductOp,
  const RCP<const Serializer>& serializer
  )
  :reductOp_(reductOp), serializer_(serializer)
{}

template<typename Ordinal, typename T, typename Serializer>
void CharToValueTypeReductionOpImp<Ordinal,T,Serializer>::reduce(
  const Ordinal     charCount
  ,const char       charInBuffer[]
  ,char             charInoutBuffer[]
  ) const
{
  ConstValueTypeDeserializationBuffer<Ordinal,T,Serializer>
    inBuffer(charCount,charInBuffer,serializer_);
  ValueTypeDeserializationBuffer<Ordinal,T,Serializer>
    inoutBuffer(charCount,charInoutBuffer,serializer_);
  reductOp_->reduce(
    inBuffer.getCount(),inBuffer.getBuffer(),inoutBuffer.getBuffer()
    );
}

//
// CharToReferenceTypeReductionOp
//

template<typename Ordinal, typename T>
CharToReferenceTypeReductionOp<Ordinal,T>::CharToReferenceTypeReductionOp(
  const RCP<const Serializer<Ordinal,T> >                 &serializer
  ,const RCP<const ReferenceTypeReductionOp<Ordinal,T> >  &reductOp
  )
  :serializer_(serializer), reductOp_(reductOp)
{}

template<typename Ordinal, typename T>
void CharToReferenceTypeReductionOp<Ordinal,T>::reduce(
  const Ordinal     charCount
  ,const char       charInBuffer[]
  ,char             charInoutBuffer[]
  ) const
{
  ConstReferenceTypeDeserializationBuffer<Ordinal,T>
    inBuffer(*serializer_,charCount,charInBuffer);
  ReferenceTypeDeserializationBuffer<Ordinal,T>
    inoutBuffer(*serializer_,charCount,charInoutBuffer);
  reductOp_->reduce(
    inBuffer.getCount(),inBuffer.getBuffer(),inoutBuffer.getBuffer()
    );
}

} // namespace Teuchos

#endif // TEUCHOS_REDUCTION_OP_HELPERS_HPP

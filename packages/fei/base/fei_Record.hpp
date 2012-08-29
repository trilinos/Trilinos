/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#ifndef _fei_Record_hpp_
#define _fei_Record_hpp_

#include <fei_macros.hpp>

namespace fei {
class FieldMask;

/** Container for record attributes. A Record is basically the
  FEI representation of a "mesh-object". */
template<typename GlobalIDType>
class Record {
public:
  /** Constructor */
  Record();

  /** copy constructor */
  Record(const Record<GlobalIDType>& src)
    : ID_(src.ID_),
      number_(src.number_),
      fieldMask_(src.fieldMask_),
      offsetIntoEqnNumbers_(src.offsetIntoEqnNumbers_),
      ownerProc_(src.ownerProc_),
      isInLocalSubdomain_(src.isInLocalSubdomain_),
      hasSlaveDof_(src.hasSlaveDof_)
    {}

  /** Destructor */
  virtual ~Record();

  /** Set globally unique identifier. */
  void setID(const GlobalIDType& ID)
  {
    ID_ = ID;
  }

  /** Query globally unique identifier. */
  GlobalIDType getID() const
  {
    return(ID_);
  }

  /** Set globally zero-based number. */
  void setNumber(const GlobalIDType& num)
  {
    number_ = num;
  }

  /** Get globally zero-based number. */
  GlobalIDType getNumber() const
  {
    return(number_);
  }

  /** operator== */
  bool operator==(const Record<GlobalIDType>& rcd) const
  {
    return( ID_ == rcd.ID_ );
  }

  /** operator!= */
  bool operator!=(const Record<GlobalIDType>& rcd) const
  {
    return( ID_ != rcd.ID_ );
  }

  /** operator< */
  bool operator<(const Record<GlobalIDType>& rcd) const
  {
    return( ID_ < rcd.ID_ );
  }

  /** operator> */
  bool operator>(const Record<GlobalIDType>& rcd) const
  {
    return( ID_ > rcd.ID_ );
  }

  /** setOwnerProc */
  void setOwnerProc(int owner)
  {
    ownerProc_ = owner;
  }

  /** getOwnerProc */
  int getOwnerProc() const
  {
    return(ownerProc_);
  }

  /** setFieldMask */
  void setFieldMask(fei::FieldMask* fm)
  {
    fieldMask_ = fm;
  }

  /** getFieldMask */
  fei::FieldMask* getFieldMask()
  {
    return( fieldMask_ );
  }

  /** getFieldMask */
  const fei::FieldMask* getFieldMask() const
  {
    return( fieldMask_ );
  }

  /** Set offset-into-equation-numbers.
   */
  void setOffsetIntoEqnNumbers(int offset)
  {
    offsetIntoEqnNumbers_ = offset;
  }

  /** Return offset-into-equation-numbers.
   */
  int getOffsetIntoEqnNumbers() const
  {
    return(offsetIntoEqnNumbers_);
  }

  bool hasSlaveDof() const
  { return( hasSlaveDof_ ); }

  void hasSlaveDof(bool flag)
  { hasSlaveDof_ = flag; }

  Record<GlobalIDType>& operator=(const Record<GlobalIDType>& src)
  {
    ID_ = src.ID_;
    number_ = src.number_;
    fieldMask_ = src.fieldMask_;
    offsetIntoEqnNumbers_ = src.offsetIntoEqnNumbers_;
    ownerProc_ = src.ownerProc_;
    isInLocalSubdomain_ = src.isInLocalSubdomain_;
    hasSlaveDof_ = src.hasSlaveDof_;
    return *this;
  }

private:

  GlobalIDType ID_;
  GlobalIDType number_;

  fei::FieldMask* fieldMask_;

  int offsetIntoEqnNumbers_;

  int ownerProc_;

public:
  /** ugh, public data member... */
  bool isInLocalSubdomain_;

private:
  bool hasSlaveDof_;
};

/** implementation detail, power-users only */
template<class GlobalIDType>
class Record_Operator {
  public:
    /** destructor */
    virtual ~Record_Operator(){}

    /** operator() */
    virtual void operator()(Record<GlobalIDType>& record) = 0;

};//class Record_Operator

template<class GlobalIDType>
fei::Record<GlobalIDType>::Record()
  : ID_(-1),
    number_(-1),
    fieldMask_(NULL),
    offsetIntoEqnNumbers_(0),
    ownerProc_(-1),
    isInLocalSubdomain_(false),
    hasSlaveDof_(false)
{
}

template<class GlobalIDType>
fei::Record<GlobalIDType>::~Record()
{
}


} //namespace fei

#endif // _fei_Record_hpp_


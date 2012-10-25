/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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


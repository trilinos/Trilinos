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
  class Record {
  public:
    /** Constructor */
    Record();

    Record(const Record& src)
     : isInLocalSubdomain_(src.isInLocalSubdomain_), ID_(src.ID_),
       fieldMask_(NULL), offsetIntoEqnNumbers_(src.offsetIntoEqnNumbers_),
       ownerProc_(src.ownerProc_), hasSlaveDof_(false)
    {}

    /** Destructor */
    virtual ~Record();

    /** Set globally unique identifier. */
    void setID(int ID)
      {
	ID_ = ID;
      }

    /** Query globally unique identifier. */
    int getID() const
      {
	return(ID_);
      }

    /** Set globally zero-based number. */
    void setNumber(int num)
      {
	number_ = num;
      }

    /** Get globally zero-based number. */
    int getNumber() const
      {
	return(number_);
      }

    /** operator== */
    bool operator==(const Record& rcd) const
      {
	return( ID_ == rcd.ID_ );
      }

    /** operator!= */
    bool operator!=(const Record& rcd) const
      {
	return( ID_ != rcd.ID_ );
      }

    /** operator< */
    bool operator<(const Record& rcd) const
      {
	return( ID_ < rcd.ID_ );
      }

    /** operator> */
    bool operator>(const Record& rcd) const
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

    /** Like a copy constructor */
    int deepCopy(const Record& rcd);

    /** ugh, public data member... */
    bool isInLocalSubdomain_;

  private:
    Record& operator=(const Record& src);

    int ID_;
    int number_;

    fei::FieldMask* fieldMask_;

    int offsetIntoEqnNumbers_;

    int ownerProc_;

    bool hasSlaveDof_;
  };

  /** lessthan operator for Record */
  struct record_lessthan {
    /** operator() */
    bool operator()(const Record* lhs,
		    const Record* rhs) const
      {
	return( *lhs < *rhs );
      }
  };

  /** implementation detail, power-users only */
  class Record_Operator {
  public:
    /** destructor */
    virtual ~Record_Operator(){}

    /** operator() */
    virtual void operator()(Record& record) = 0;

  };//class Record_Operator

} //namespace fei

#endif // _fei_Record_hpp_


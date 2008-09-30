/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _snl_fei_Constraint_hpp_
#define _snl_fei_Constraint_hpp_

#include <fei_macros.hpp>
#include <fei_fwd.hpp>
#include <fei_VectorSpace.hpp>
#include <snl_fei_RecordCollection.hpp>

namespace snl_fei {

  /** container for constraint attributes */
  template<class RecordType, class RecordType_COMPARE>
  class Constraint {
  public:
    /** constructor */
    Constraint(int id=0, bool isPenalty=false);

    /** constructor */
    Constraint(int id,
	       int constraintIDType,
	       bool isSlave,
	       bool isPenalty,
	       int numIDs,
	       const int* idTypes,
	       const int* IDs,
	       const int* fieldIDs,
	       int offsetOfSlave,
	       int offsetIntoSlaveField,
	       const double* weights,
	       double rhsValue,
	       fei::VectorSpace* vspace);

    /** destructor */
    virtual ~Constraint();

    /** get constraint identifier */
    int getConstraintID() { return( constraintID_ ); }

    /** set constraint identifier. power-users only */
    void setConstraintID(int id) { constraintID_ = id; }

    /** get the identifier-type that the fei uses to reference constraints */ 
    int getIDType() { return( idType_ ); }

    /** set the identifier-type that the fei uses to reference constraints.
      power-users only, this is a dangerous function with side-effects */ 
    void setIDType(int idType) { idType_ = idType; }

    /** query whether constraint is a penalty constraint */
    bool isPenalty() { return( isPenalty_ ); }

    /** set whether constraint is a penalty constraint. Another dangerous
     function for power-users. */
    void setIsPenalty(bool isPenalty) { isPenalty_ = isPenalty; }

    /** get equation-number of constraint. (only valid if lagrange-multiplier)
    */
    int getEqnNumber() { return( eqnNumber_ ); }

    /** set equation-number of constraint. (only valid if lagrange-multiplier)
    */
    void setEqnNumber(int eqn) { eqnNumber_ = eqn; }

    /** get block-equation number of constraint. (only valid if
     lagrange-multiplier */
    int getBlkEqnNumber() { return( blkEqnNumber_ ); }

    /** set block-equation number of constraint. (only valid if
     lagrange-multiplier */
    void setBlkEqnNumber(int blkEqn) { blkEqnNumber_ = blkEqn; }

    /** not for public use (note to self: review this class, remove the
     many dangerous, ill-defined functions) */
    int allocate();


    /** intended for fei-implementation use only */
    RecordType getSlave() { return( slave_ ); }

    /** intended for fei-implementation use only */
    void setSlave(const RecordType& slv) { slave_ = slv; }

    /** if slave constraint, return constrained field-identifier of
        slaved mesh-object */
    int getSlaveFieldID() { return( slaveField_ ); }

    /** if slave constraint, set constrained field-identifier of
        slaved mesh-object */
    void setSlaveFieldID(int f) { slaveField_ = f; }

    /** get offset of slaved field-component */
    int getOffsetIntoSlaveField() { return( offsetIntoSlaveField_ ); }

    /** set offset of slaved field-component */
    void setOffsetIntoSlaveField(int offset) { offsetIntoSlaveField_ = offset; }


    /** get master mesh-objects */
    feiArray<RecordType,RecordType_COMPARE>* getMasters() { return( masters_ ); }

    /** get identifier-types of master mesh-objects */
    feiArray<int>* getMasterIDTypes() { return( masterIDTypes_ ); }

    /** get field-identifiers of master mesh-objects */
    feiArray<int>* getMasterFieldIDs() { return( masterFields_ ); }

    /** get weight-coefficients of master mesh-objects */
    feiArray<double>* getMasterWeights() { return( masterWeights_ ); }


    /** get right-hand-side value of constraint */
    double getRHSValue() { return( rhsValue_ ); }

    /** set right-hand-side value of constraint */
    void setRHSValue(double rhs) { rhsValue_ = rhs; }
 
    /** operator!= */
    bool operator!=(const Constraint<RecordType,RecordType_COMPARE>& rhs);

    /** query whether connectivity is the same as specified constraint */
    bool structurallySame(const Constraint<RecordType,RecordType_COMPARE>& rhs);

  private:
    Constraint(const Constraint<RecordType,RecordType_COMPARE>& src);
    Constraint<RecordType,RecordType_COMPARE>& operator=(const Constraint<RecordType,RecordType_COMPARE>& src)
      {
	return(*this);
      }

    int constraintID_;
    int idType_;
    bool isPenalty_;

    int eqnNumber_;
    int blkEqnNumber_;

    RecordType slave_;
    int slaveField_;
    int offsetIntoSlaveField_;

    feiArray<RecordType,RecordType_COMPARE>* masters_;
    feiArray<int>* masterIDTypes_;
    feiArray<int>* masterFields_;
    feiArray<double>* masterWeights_;

    double rhsValue_;

  };//class Constraint
} //namespace snl_fei

#include <snl_fei_Constraint.hpp>

#include <feiArray.hpp>

//----------------------------------------------------------------------------
template<class RecordType,class RecordType_COMPARE>
inline snl_fei::Constraint<RecordType,RecordType_COMPARE>::Constraint(int id, bool isPenalty)
  : constraintID_(id),
    idType_(0),
    isPenalty_(isPenalty),
    slaveField_(0),
    offsetIntoSlaveField_(0),
    masters_(NULL),
    masterIDTypes_(NULL),
    masterFields_(NULL),
    masterWeights_(NULL),
    rhsValue_(0.0)
{
}

//----------------------------------------------------------------------------
template<class RecordType,class RecordType_COMPARE>
inline snl_fei::Constraint<RecordType,RecordType_COMPARE>::Constraint(int id,
					    int constraintIDType,
					    bool isSlave,
					    bool isPenalty,
					    int numIDs,
					    const int* idTypes,
					    const int* IDs,
					    const int* fieldIDs,
					    int offsetOfSlave,
					    int offsetIntoSlaveField,
					    const double* weights,
					    double rhsValue,
					    fei::VectorSpace* vspace)
  : constraintID_(id),
    idType_(constraintIDType),
    isPenalty_(isPenalty),
    eqnNumber_(-1),
    blkEqnNumber_(-1), 
    slaveField_(0),
    offsetIntoSlaveField_(offsetIntoSlaveField),
    masters_(NULL),
    masterIDTypes_(NULL),
    masterFields_(NULL),
    masterWeights_(NULL),
    rhsValue_(rhsValue)
{
}

//----------------------------------------------------------------------------
namespace snl_fei {
template<>
inline snl_fei::Constraint<fei::Record*,fei::record_lessthan>::Constraint(int id,
					    int constraintIDType,
					    bool isSlave,
					    bool isPenalty,
					    int numIDs,
					    const int* idTypes,
					    const int* IDs,
					    const int* fieldIDs,
					    int offsetOfSlave,
					    int offsetIntoSlaveField,
					    const double* weights,
					    double rhsValue,
					    fei::VectorSpace* vspace)
  : constraintID_(id),
    idType_(constraintIDType),
    isPenalty_(isPenalty),
    eqnNumber_(-1),
    blkEqnNumber_(-1), 
    slaveField_(0),
    offsetIntoSlaveField_(offsetIntoSlaveField),
    masters_(NULL),
    masterIDTypes_(NULL),
    masterFields_(NULL),
    masterWeights_(NULL),
    rhsValue_(rhsValue)
{
  allocate();

  int weightsOffset = 0;
  for(int i=0; i<numIDs; ++i) {
    snl_fei::RecordCollection* recordCollection = NULL;
    vspace->getRecordCollection(idTypes[i],recordCollection);

    fei::Record* rec = NULL;

    vspace->addDOFs(fieldIDs[i], 1, idTypes[i], 1, &(IDs[i]));
    rec = recordCollection->getRecordWithID(IDs[i]);

    unsigned fieldSize = vspace->getFieldSize(fieldIDs[i]);

    if (isSlave && i == offsetOfSlave) {
      rec->hasSlaveDof(true);
      setSlave(rec);
      setSlaveFieldID(fieldIDs[i]);
      setOffsetIntoSlaveField(offsetIntoSlaveField);
      weightsOffset += fieldSize;
    }
    else {
      getMasters()->append(rec);
      getMasterIDTypes()->append(idTypes[i]);
      getMasterFieldIDs()->append(fieldIDs[i]);

      if (weights != NULL) {
	for(unsigned j=0; j<fieldSize; ++j) {
	  masterWeights_->append(weights[weightsOffset++]);
	}
      }
    }
  }
}

}//namespace snl_fei

//----------------------------------------------------------------------------
template<class RecordType,class RecordType_COMPARE>
inline snl_fei::Constraint<RecordType,RecordType_COMPARE>::Constraint(const Constraint<RecordType,RecordType_COMPARE>& src)
  : constraintID_(-1),
    idType_(0),
    isPenalty_(false),
    eqnNumber_(-1),
    blkEqnNumber_(-1), 
    slaveField_(0),
    offsetIntoSlaveField_(0),
    masters_(NULL),
    masterIDTypes_(NULL),
    masterFields_(NULL),
    masterWeights_(NULL),
    rhsValue_(0.0)
{
}

//----------------------------------------------------------------------------
template<class RecordType,class RecordType_COMPARE>
inline snl_fei::Constraint<RecordType,RecordType_COMPARE>::~Constraint()
{
  delete masters_;
  delete masterIDTypes_;
  delete masterFields_;
  delete masterWeights_;
}

//----------------------------------------------------------------------------
template<class RecordType,class RecordType_COMPARE>
inline int snl_fei::Constraint<RecordType,RecordType_COMPARE>::allocate()
{
  if (masters_ != NULL) delete masters_;
  masters_ = new feiArray<RecordType,RecordType_COMPARE>(0, 5);

  if (masterIDTypes_ != NULL) delete masterIDTypes_;
  masterIDTypes_ = new feiArray<int>(0, 5);

  if (masterFields_ != NULL) delete masterFields_;
  masterFields_ = new feiArray<int>(0, 5);

  if (masterWeights_ != NULL) delete masterWeights_;
  masterWeights_ = new feiArray<double>;

  return(0);
}

//----------------------------------------------------------------------------
template<class RecordType,class RecordType_COMPARE>
inline bool snl_fei::Constraint<RecordType,RecordType_COMPARE>::operator!=(const snl_fei::Constraint<RecordType,RecordType_COMPARE>& rhs)
{
  if (constraintID_ != rhs.constraintID_ ||
      idType_ != rhs.idType_ ||
      isPenalty_ != rhs.isPenalty_ ||
      eqnNumber_ != rhs.eqnNumber_ ||
      blkEqnNumber_ != rhs.blkEqnNumber_ ||
      slaveField_ != rhs.slaveField_ ||
      offsetIntoSlaveField_ != rhs.offsetIntoSlaveField_ ||
      rhsValue_ != rhs.rhsValue_) {
    return( true );
  }

  if (masters_ != NULL) {
    if (rhs.masters_ == NULL) return(true);
    if (*masters_ != *(rhs.masters_)) return(true);
  }

  if (masterIDTypes_ != NULL) {
    if (rhs.masterIDTypes_ == NULL) return(true);
    if (*masterIDTypes_ != *(rhs.masterIDTypes_)) return(true);
  }

  if (masterFields_ != NULL) {
    if (rhs.masterFields_ == NULL) return(true);
    if (*masterFields_ != *(rhs.masterFields_)) return(true);
  }

  if (masterWeights_ != NULL) {
    if (rhs.masterWeights_ == NULL) return(true);
    if (*masterWeights_ != *(rhs.masterWeights_)) return(true);
  }

  return(false);
}

//----------------------------------------------------------------------------
template<class RecordType,class RecordType_COMPARE>
inline bool snl_fei::Constraint<RecordType,RecordType_COMPARE>::structurallySame(const Constraint<RecordType,RecordType_COMPARE>& rhs)
{
  if (constraintID_ != rhs.constraintID_ ||
      idType_ != rhs.idType_ ||
      isPenalty_ != rhs.isPenalty_ ||
      eqnNumber_ != rhs.eqnNumber_ ||
      blkEqnNumber_ != rhs.blkEqnNumber_ ||
      slaveField_ != rhs.slaveField_ ||
      offsetIntoSlaveField_ != rhs.offsetIntoSlaveField_) {
    return( false );
  }

  if (masters_ != NULL) {
    if (rhs.masters_ == NULL) return(false);
    if (*masters_ != *(rhs.masters_)) return(false);
  }

  if (masterIDTypes_ != NULL) {
    if (rhs.masterIDTypes_ == NULL) return(false);
    if (*masterIDTypes_ != *(rhs.masterIDTypes_)) return(false);
  }

  if (masterFields_ != NULL) {
    if (rhs.masterFields_ == NULL) return(false);
    if (*masterFields_ != *(rhs.masterFields_)) return(false);
  }

  return(true);
}

#endif // _snl_fei_Constraint_hpp_

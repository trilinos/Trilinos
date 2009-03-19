#ifndef __PB_EpetraInverseOpWrapper_hpp__
#define __PB_EpetraInverseOpWrapper_hpp__

#include "PB_EpetraOperatorWrapper.hpp"

namespace PB {
namespace Epetra {

class EpetraInverseOpWrapper : public EpetraOperatorWrapper {
public:
    EpetraInverseOpWrapper(const RCP<const MappingStrategy> & forwardMaps) 
       : EpetraOperatorWrapper(forwardMaps) {}

    /** */
    virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    /** */
    virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    /** */
    virtual const Epetra_Map& OperatorDomainMap() const; 

    /** */
    virtual const Epetra_Map& OperatorRangeMap() const;
};

} // end Epetra
} // end PB

#endif

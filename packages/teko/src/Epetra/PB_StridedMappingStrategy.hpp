#ifndef __PB_StridedMappingStrategy_hpp__
#define __PB_StridedMappingStrategy_hpp__

// stl includes
#include <vector>

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"

// Epetra includes
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"

// PB includes
#include "Epetra/PB_EpetraOperatorWrapper.hpp"

namespace PB {
namespace Epetra {

class StridedMappingStrategy : public MappingStrategy {
public:
   // constructors

   // Creates a strided mapping strategy. This class is useful
   // for breaking up nodally ordered matrices (i.e. the unknowns
   // in a FEM problem are ordered [u0,v0,p0,u1,v1,p1,...]). Current
   // implimentation only supports a fixed number of variables
   //
   //    arguments: 
   //       vars - Vector describing the blocking of variables
   //       map  - original Epetra_Map to be broken up
   //       comm - Epetra_Comm object related to the map
   //
   StridedMappingStrategy(const std::vector<int> & vars,const Teuchos::RCP<const Epetra_Map> & map, const Epetra_Comm & comm);

   // Virtual member functions inherited from the
   // PB::Epetra::MappingStrategy parent class.
   /////////////////////////////////////////////////////////

   // Virtual function defined in MappingStrategy.  This copies
   // an Epetra_MultiVector into a Thyra::MultiVectorBase with
   // blocking handled by the strides defined in the constructor.
   //
   //   arguments:
   //      X       - source Epetra_MultiVector
   //      thyra_X - destination Thyra::MultiVectorBase
   //      eow     - Operator that defines the transition...this may
   //                be removed in the future
   //
   virtual void copyEpetraIntoThyra(const Epetra_MultiVector& x, 
                                    const Teuchos::Ptr<Thyra::MultiVectorBase<double> > & thyraVec,
                                    const PB::Epetra::EpetraOperatorWrapper & eow) const;

   // Virtual function defined in MappingStrategy.  This copies
   // an Epetra_MultiVector into a Thyra::MultiVectorBase with
   // blocking handled by the strides defined in the constructor.
   //
   //   arguments:
   //      thyra_Y - source Thyra::MultiVectorBase
   //      Y       - destination Epetra_MultiVector
   //      eow     - Operator that defines the transition...this may
   //                be removed in the future
   //
   virtual void copyThyraIntoEpetra(const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & thyraVec, 
                                    Epetra_MultiVector& v,
                                    const PB::Epetra::EpetraOperatorWrapper & eow) const;

   // Returns the domain and range maps used by this class.
   // This faciliates building an Epetra_Operator around this
   // class with its core functionality being a Thyra::LinearOpBase
   // operator
   //
   //    returns: Range map corresponding to this class
   //
   virtual const Teuchos::RCP<const Epetra_Map> domainMap() const
   { return domainMap_; }

   // Returns the domain and range maps used by this class.
   // This faciliates building an Epetra_Operator around this
   // class with its core functionality being a Thyra::LinearOpBase
   // operator
   //
   //    returns: Range map corresponding to this class
   //
   virtual const Teuchos::RCP<const Epetra_Map> rangeMap() const
   { return rangeMap_; }

   // A function for my sanity
   //
   //    returns: String with description of this class
   //
   virtual std::string toString() const
   { return std::string("StridedMappingStrategy"); }

   // Locally (a concrete insantiation of this this class)
   // useful functions
   /////////////////////////////////////////////////////////

   // this is the core routine that builds the maps
   // and importers/exporters neccessary for all the
   // transfers. Currently it simply calls out to the
   // interlaced epetra functions. (Comment: this
   // routine should probably be private or protected
   // ... it is basically the meat of the constructor)
   //
   //    arguments:
   //       vars    - Vector describing the blocking of variables
   //       baseMap - basic map to use in the transfers
   //       comm    - Epetra_Comm object
   //
   void buildBlockTransferData(const std::vector<int> & vars,const Teuchos::RCP<const Epetra_Map> & baseMap, const Epetra_Comm & comm);

   // Get the individual block maps underlying that
   // make up a strided vector/operator. These are
   // useful if you want to tear something ... i.e.
   // a matrix ... apart.
   //
   //    returns: Return a vector of block maps 
   //             created for this strided operator 
   //
   const std::vector<std::pair<int,Teuchos::RCP<Epetra_Map> > > & getMaps() const
   { return blockMaps_; }

   // Builds a blocked Thyra operator that uses the strided
   // mapping strategy to define sub blocks.
   //
   //    arguments:
   //       mat - Epetra_CrsMatrix with FillComplete called, this
   //             matrix is assumed to be square, with the same
   //             range and domain maps
   //    returns: Blocked Thyra linear operator with sub blocks
   //             defined by this mapping strategy
   //
   const Teuchos::RCP<const Thyra::LinearOpBase<double> > 
   buildBlockedThyraOp(const Teuchos::RCP<const Epetra_CrsMatrix> & mat) const;

protected:
   // member variables
   ///////////////////////////////////////////

   // storage for sanity
   Teuchos::RCP<const Epetra_Map> domainMap_; 
   Teuchos::RCP<const Epetra_Map> rangeMap_;

   // block transfer data
   std::vector<std::pair<int,Teuchos::RCP<Epetra_Map> > > blockMaps_;
   std::vector<Teuchos::RCP<Epetra_Export> > blockExport_;
   std::vector<Teuchos::RCP<Epetra_Import> > blockImport_;
};

} // end namespace Epetra
} // end namespace PB

#endif

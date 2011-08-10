#ifndef __Panzer_SGEpetraLinearObjContainer_hpp__
#define __Panzer_SGEpetraLinearObjContainer_hpp__

#include "Panzer_config.hpp"
#ifdef HAVE_STOKHOS

#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Teuchos_RCP.hpp"
#include "Stokhos_OrthogPolyExpansion.hpp"

#include <vector>

namespace panzer {

/** Linear object container for SG-Epetra objects.
  */
class SGEpetraLinearObjContainer : public LinearObjContainer {
public:
   typedef std::vector<Teuchos::RCP<EpetraLinearObjContainer> > CoeffVector;
   typedef CoeffVector::iterator iterator;
   typedef CoeffVector::const_iterator const_iterator;

   SGEpetraLinearObjContainer(const CoeffVector & coeffs,
                              const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > & basis);

   virtual void initialize();

   CoeffVector::iterator begin() { return coeffs_.begin(); }
   CoeffVector::iterator end() { return coeffs_.end(); }

   CoeffVector::const_iterator begin() const { return coeffs_.begin(); }
   CoeffVector::const_iterator end() const { return coeffs_.end(); }

   Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > getExpansion() const
   { return expansion_; }

private:
   CoeffVector coeffs_;
   Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion_;
};

}

#endif
#endif

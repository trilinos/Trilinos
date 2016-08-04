#ifndef ROL_EQUALITYCONSTRAINT_PARTITIONED_H
#define ROL_EQUALITYCONSTRAINT_PARTITIONED_H

#include "ROL_EqualityConstraint.hpp"
#include "ROL_PartitionedVector.hpp"

namespace ROL {

/** @ingroup func_group
    \class ROL::EqualityConstraint_Partitioned
    \brief Allows composition of equality constraints
*/

template<class Real>
class EqualityConstraint_Partitioned : public EqualityConstraint<Real> {

  typedef EqualityConstraint<Real>  EC;
  typedef Vector<Real>              V;
  typedef PartitionedVector<Real>   PV;
  typedef typename std::vector<Real>::size_type uint;

private:

  const std::vector<Teuchos::RCP<EC> > con_;
  uint dim_;

public:

  EqualityConstraint_Partitioned( const std::vector<Teuchos::RCP<EqualityConstraint<Real> > > &con ) :
    con_(con), dim_(con.size()) {

  }
  
  virtual ~EqualityConstraint_Partitioned() {}

  virtual void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
    // Downcast c to PartitionedVector
    PV& cpv = Teuchos::dyn_cast<PV>(c);

    // Iterate over constraints
    for( uint k=0; k<dim_; ++k ) {
      // Evaluate constraint for each c contribution
      con_[k]->value(*(cpv.get(k)),x,tol);
    }

  }

  virtual void applyJacobian(Vector<Real> &jv, const Vector<Real> &v,
                     const Vector<Real> &x, Real &tol) {

    // Downcast jv to PartitionedVector
    PV& jvpv = Teuchos::dyn_cast<PV>(jv);

    // Iterate over constraints
    for( uint k=0; k<dim_; ++k ) {
      // Evaluate jacobian contributions
      con_[k]->applyJacobian( *(jvpv.get(k)), v, x, tol);
    } 

  }

  virtual void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v,
                            const Vector<Real> &x, Real &tol) {

    const PV& vpv = Teuchos::dyn_cast<const PV>(v);

    Teuchos::RCP<V> temp = ajv.clone();
    ajv.zero();
    
    for( uint k=0; k<dim_; ++k ) {
      con_[k]->applyAdjointJacobian( *temp, *(vpv.get(k)), x, tol);
      ajv.plus(*temp);
    }   
  }

  virtual void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u,
                           const Vector<Real> &v, const Vector<Real> &x,
                           Real &tol) {

    const PV& upv = Teuchos::dyn_cast<const PV>(u);
    
    Teuchos::RCP<V> temp = ahuv.clone();
    ahuv.zero();

    for( uint k=0; k<dim_; ++k ) {
      con_[k]->applyAdjointHessian( *temp, *(upv.get(k)), v, x, tol );
      ahuv.plus( *temp );
    }
  }  

  virtual void applyPreconditioner(Vector<Real> &pv, const Vector<Real> &v,
                                   const Vector<Real> &x, const Vector<Real> &g,
                                   Real &tol) {

    const PV& vpv = Teuchos::dyn_cast<const PV>(v); 
    PV& pvpv = Teuchos::dyn_cast<PV>(pv);

    for( uint k=0; k<dim_; ++k ) {
      con_[k]->applyPreconditioner( *(pvpv.get(k)), *(vpv.get(k)), x, g, tol );
    }

  }

}; // class EqualityConstraint_Partitioned

// Helper methods
template<class Real> 
Teuchos::RCP<EqualityConstraint<Real> > 
CreateEqualityConstraintPartitioned( const Teuchos::RCP<EqualityConstraint<Real> > &con1,
                                     const Teuchos::RCP<EqualityConstraint<Real> > &con2 ) {
  using Teuchos::RCP; using Teuchos::rcp;

  typedef EqualityConstraint_Partitioned<Real> ECP;
  typedef RCP<EqualityConstraint<Real> > RCPEC;
  RCPEC con[] = { con1, con2 };

  return rcp( new ECP( std::vector<RCPEC>( con, con+2 ) ) );
}

template<class Real> 
Teuchos::RCP<EqualityConstraint<Real> > 
CreateEqualityConstraintPartitioned( const Teuchos::RCP<EqualityConstraint<Real> > &con1,
                                     const Teuchos::RCP<EqualityConstraint<Real> > &con2,
                                     const Teuchos::RCP<EqualityConstraint<Real> > &con3 ) {
  using Teuchos::RCP; using Teuchos::rcp;

  typedef EqualityConstraint_Partitioned<Real> ECP;
  typedef RCP<EqualityConstraint<Real> > RCPEC;
  RCPEC con[] = { con1, con2, con3 };

  return rcp( new ECP( std::vector<RCPEC>( con, con+3 ) ) );
}


} // namespace ROL


#endif //  ROL_PARTITIONEQUALITYCONSTRAINT_H

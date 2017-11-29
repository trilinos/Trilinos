
#ifndef ROL_PDEOPT_ENERGY_OBJECTIVE
#define ROL_PDEOPT_ENERGY_OBJECTIVE

#include "ROL_Objective.hpp"
#include "../TOOLS/assembler.hpp"

// Do not instantiate the template in this translation unit.
extern template class Assembler<double>;

template<class Real>
class EnergyObjective : public ROL::Objective<Real> {
private:
  const ROL::SharedPointer<PDE<Real> > pde_;
  ROL::SharedPointer<Assembler<Real> > assembler_;
  bool assembleRHS_, assembleJ1_;

  ROL::SharedPointer<Tpetra::MultiVector<> > cvec_;
  ROL::SharedPointer<Tpetra::MultiVector<> > uvec_;

  ROL::SharedPointer<Tpetra::MultiVector<> > res_;
  ROL::SharedPointer<Tpetra::CrsMatrix<> >   jac_;

  void assemble(void) {
    // Assemble affine term.
    if (assembleRHS_) {
      assembler_->assemblePDEResidual(res_,pde_,uvec_);
    }
    assembleRHS_ = false;
    // Assemble jacobian_1.
    if (assembleJ1_) {
      assembler_->assemblePDEJacobian1(jac_,pde_,uvec_);
    }
    assembleJ1_ = false;
  }

public:
  EnergyObjective(const ROL::SharedPointer<PDE<Real> > &pde,
                  const ROL::SharedPointer<MeshManager<Real> > &meshMgr,
                  const ROL::SharedPointer<const Teuchos::Comm<int> > &comm,
                  Teuchos::ParameterList &parlist,
                  std::ostream &outStream = std::cout)
    : pde_(pde), assembleRHS_(true), assembleJ1_(true) {
    // Construct assembler.
    assembler_ = ROL::makeShared<Assembler<Real>>(pde_->getFields(),meshMgr,comm,parlist,outStream);
    assembler_->setCellNodes(*pde_);
    // Initialize zero vectors.
    cvec_ = assembler_->createResidualVector();
    uvec_ = assembler_->createStateVector();
    uvec_->putScalar(static_cast<Real>(0));
    assemble();
  }

  const ROL::SharedPointer<Assembler<Real> > getAssembler(void) const {
    return assembler_;
  }

  Real value(const ROL::Vector<Real> &u, Real &tol) {
    ROL::SharedPointer<const Tpetra::MultiVector<> > uf = getConstField(u);
    const Real half(0.5), one(1);
    jac_->apply(*uf,*cvec_);
    cvec_->update(one,*res_,half);
    Teuchos::Array<Real> val(1,0);
    cvec_->dot(*uf, val.view(0,1));
    return val[0];
  }

  void gradient(ROL::Vector<Real> &g, const ROL::Vector<Real> &u, Real &tol) {
    ROL::SharedPointer<Tpetra::MultiVector<> >       gf = getField(g);
    ROL::SharedPointer<const Tpetra::MultiVector<> > uf = getConstField(u);
    const Real one(1);
    gf->scale(one,*res_);
    jac_->apply(*uf,*cvec_);
    gf->update(one,*cvec_,one);
  }

  void hessVec(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, Real &tol) {
    ROL::SharedPointer<Tpetra::MultiVector<> >      hvf = getField(hv);
    ROL::SharedPointer<const Tpetra::MultiVector<> > vf = getConstField(v);
    ROL::SharedPointer<const Tpetra::MultiVector<> > uf = getConstField(u);
    jac_->apply(*vf,*hvf);
  }

  void precond(ROL::Vector<Real> &Pv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, Real &tol) {
    Pv.set(v.dual());
  }

private: // Vector accessor functions

  ROL::SharedPointer<const Tpetra::MultiVector<> > getConstField(const ROL::Vector<Real> &x) const {
    ROL::SharedPointer<const Tpetra::MultiVector<> > xp;
    try {
      xp = dynamic_cast<const ROL::TpetraMultiVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      ROL::SharedPointer<const ROL::TpetraMultiVector<Real> > xvec
        = dynamic_cast<const PDE_OptVector<Real>&>(x).getField();
      if (xvec == ROL::nullPointer) {
        xp = ROL::nullPointer;
      }
      else {
        xp = xvec->getVector();
      }
    }
    return xp;
  }

  ROL::SharedPointer<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    ROL::SharedPointer<Tpetra::MultiVector<> > xp;
    try {
      xp = dynamic_cast<ROL::TpetraMultiVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      ROL::SharedPointer<ROL::TpetraMultiVector<Real> > xvec
        = dynamic_cast<PDE_OptVector<Real>&>(x).getField();
      if ( xvec == ROL::nullPointer ) {
        xp = ROL::nullPointer;
      }
      else {
        xp = xvec->getVector();
      }
    }
    return xp;
  }

};

#endif

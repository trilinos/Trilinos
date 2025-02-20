#ifndef PYTRILINOS2_THYRA_CUSTOM_HPP
#define PYTRILINOS2_THYRA_CUSTOM_HPP

#include <Thyra_LinearOpWithSolveBase_decl.hpp>

template<typename T>
void define_solve(T cl) {
  using SCALAR = typename T::type::scalar_type;

  cl.def("solve",[](Teuchos::RCP<Thyra::LinearOpWithSolveBase<SCALAR> > &m,
                    const Thyra::EOpTransp A_trans,
                    const Thyra::MultiVectorBase<SCALAR> &B,
                    const Teuchos::RCP<Thyra::MultiVectorBase<SCALAR> > &X)
  { return m->solve(A_trans, B, X.ptr()); });
}

#endif

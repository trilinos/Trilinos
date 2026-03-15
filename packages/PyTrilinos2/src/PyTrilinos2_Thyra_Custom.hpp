#ifndef PYTRILINOS2_THYRA_CUSTOM_HPP
#define PYTRILINOS2_THYRA_CUSTOM_HPP

#include <Thyra_LinearOpWithSolveBase_decl.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>

template<typename T>
void define_solve(T cl) {
  using SCALAR = typename T::type::scalar_type;

  using scalar_type = Tpetra::Details::DefaultTypes::scalar_type;
  using local_ordinal_type = Tpetra::Details::DefaultTypes::local_ordinal_type;
  using global_ordinal_type = Tpetra::Details::DefaultTypes::global_ordinal_type;
  using node_type = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType;

  cl.def("solve",[](Teuchos::RCP<Thyra::LinearOpWithSolveBase<SCALAR> > &m,
                    const Thyra::EOpTransp A_trans,
                    const Thyra::MultiVectorBase<SCALAR> &B,
                    const Teuchos::RCP<Thyra::MultiVectorBase<SCALAR> > &X)
    { return m->solve(A_trans, B, X.ptr()); });

  cl.def("solve",[](Teuchos::RCP<Thyra::LinearOpWithSolveBase<SCALAR> > &A,
                    const Thyra::EOpTransp A_trans,
                    const Teuchos::RCP<const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type >> &B,
                    const Teuchos::RCP<Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> > &X)
    { return Thyra::solve(*A, A_trans, B, X); });
}

#endif

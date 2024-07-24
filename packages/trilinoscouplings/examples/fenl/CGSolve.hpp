//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_EXAMPLE_CG_SOLVE
#define KOKKOS_EXAMPLE_CG_SOLVE

#include <cmath>
#include <limits>
#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <Kokkos_Timer.hpp>
#include <Kokkos_ArithTraits.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_CrsMatrix.hpp>
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

template <typename Scalar>
inline
Scalar all_reduce( Scalar local , const Teuchos::RCP<const Teuchos::Comm<int> >& comm )
{
  Scalar global = 0 ;
  Teuchos::reduceAll( *comm , Teuchos::REDUCE_SUM , 1 , & local , & global );
  return global ;
}

template <typename Scalar>
inline
Teuchos::Array<Scalar>
all_reduce( Teuchos::Array<Scalar> local ,
            const Teuchos::RCP<const Teuchos::Comm<int> >& comm )
{
  const int sz = local.size();
  Teuchos::Array<Scalar> global(sz);
  Teuchos::reduceAll( *comm , Teuchos::REDUCE_SUM , sz , local.getRawPtr() , global.getRawPtr() );
  return global ;
}

template <typename Scalar>
inline
Scalar all_reduce_max( Scalar local , const Teuchos::RCP<const Teuchos::Comm<int> >& comm )
{
  Scalar global = 0 ;
  Teuchos::reduceAll( *comm , Teuchos::REDUCE_MAX , 1 , & local , & global );
  return global ;
}

struct result_struct {
  double add_time, dot_time, matvec_time, iter_time, prec_setup_time,
    prec_apply_time, total_time;
  int iteration;
  std::vector<int> ensemble_its;
  double norm_res;
  result_struct() :
    add_time(0), dot_time(0), matvec_time(0), iter_time(0), prec_setup_time(0),
    prec_apply_time(0), total_time(0), iteration(0), norm_res(0) {}
  result_struct(double add, double dot, double matvec, double prec_setup,
                double prec_apply, double total, int niter, double res) :
    add_time(add),
    dot_time(dot),
    matvec_time(matvec),
    iter_time(add+dot+matvec),
    prec_setup_time(prec_setup),
    prec_apply_time(prec_apply),
    total_time(total),
    iteration(niter),
    norm_res(res) {};
};

template<class CrsMatrix, class Vector>
result_struct cg_solve(
  Teuchos::RCP<CrsMatrix> A, Teuchos::RCP<Vector> b, Teuchos::RCP<Vector> x,
  int max_iter = 200,
  typename Kokkos::ArithTraits<typename CrsMatrix::scalar_type>::mag_type tolerance =
    Kokkos::ArithTraits<typename CrsMatrix::scalar_type>::epsilon(),
  int print = 0)
{
  Kokkos::Timer total_timer;

  typedef typename CrsMatrix::scalar_type ScalarType;
  typedef Kokkos::ArithTraits<ScalarType> KAT;
  typedef typename KAT::mag_type MagnitudeType;
  typedef typename CrsMatrix::local_ordinal_type LocalOrdinalType;
  Teuchos::RCP<Vector> r,p,Ap;

  // create temporary Vectors
  r = Tpetra::createVector<ScalarType>(A->getRangeMap());
  p = Tpetra::createVector<ScalarType>(A->getRangeMap());
  Ap = Tpetra::createVector<ScalarType>(A->getRangeMap());

  // fill with initial Values (make this a functor call or something)
  x->putScalar(0.0);

  MagnitudeType normr = 0;
  MagnitudeType rtrans = 0;
  MagnitudeType oldrtrans = 0;

  LocalOrdinalType print_freq = max_iter/10;
  if (print_freq>50) print_freq = 50;
  if (print_freq<1)  print_freq = 1;

  double dottime = 0;
  double addtime = 0;
  double matvectime = 0;

  Kokkos::Timer timer;
  p->update(1.0,*x,0.0,*x,0.0);
  addtime += timer.seconds(); timer.reset();

  A->apply(*p, *Ap);
  matvectime += timer.seconds(); timer.reset();

  r->update(1.0,*b,-1.0,*Ap,0.0);
  addtime += timer.seconds(); timer.reset();

  rtrans = r->dot(*r);
  dottime += timer.seconds(); timer.reset();

  normr = std::sqrt(rtrans);

  MagnitudeType normr_orig = normr;

  if (print) {
    std::cout << "Initial Residual = "<< normr << std::endl;
  }

  MagnitudeType brkdown_tol = KAT::epsilon();

  // Count external so that we keep iteration count in the end
  LocalOrdinalType k;
  for(k=1; k <= max_iter && normr / normr_orig > tolerance; ++k) {
    if (k == 1) {
      p->update(1.0,*r,0.0,*r,0.0);
      addtime += timer.seconds(); timer.reset();
    }
    else {
      oldrtrans = rtrans;
      rtrans = r->dot(*r);
      dottime += timer.seconds(); timer.reset();
      MagnitudeType beta = rtrans/oldrtrans;
      p->update(beta,*p,1.0,*r,0.0);
      addtime += timer.seconds(); timer.reset();
    }
    normr = std::sqrt(rtrans);
    if (print && (k%print_freq==0 || k==max_iter)) {
      std::cout << "Iteration = "<<k<<"   Residual = "<<normr<<std::endl;
    }

    MagnitudeType alpha = 0;
    MagnitudeType p_ap_dot = 0;
    A->apply(*p, *Ap);
    matvectime += timer.seconds(); timer.reset();

    p_ap_dot = Ap->dot(*p);
    dottime += timer.seconds(); timer.reset();

   if (p_ap_dot < brkdown_tol) {
      if (p_ap_dot < 0 ) {
        std::cerr << "cg_solve ERROR, numerical breakdown!"<<std::endl;
        return result_struct();
      }
      else brkdown_tol = 0.1 * p_ap_dot;
    }
    alpha = rtrans/p_ap_dot;

    x->update(1.0,*x,alpha,*p,0.0);
    r->update(1.0,*r,-alpha,*Ap,0.0);
    addtime += timer.seconds(); timer.reset();

  }
  rtrans = r->dot(*r);

  normr = std::sqrt(rtrans);

  double totaltime = total_timer.seconds();
  double prec_setup = 0;
  double prec_apply = 0;
  return result_struct(addtime,dottime,matvectime,prec_setup,prec_apply,
                       totaltime,k-1,FENL::scalar_norm(normr));
}




} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXAMPLE_CG_SOLVE */

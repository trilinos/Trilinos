#include <iostream>
#include "Kokkos_Core.hpp"
#include "Kokkos_View.hpp"
#include "KokkosBlas.hpp"
#include "Teuchos_TimeMonitor.hpp"


int main(int narg, char **arg) 
{
  Kokkos::ScopeGuard scope(narg, arg);
  MPI_Init(&narg, &arg);

  if (narg < 4) {
    std::cout << "Usage: a.out nVectors vectorLen nIterations" << std::endl;
  }
  int nvec = std::atoi(arg[1]);
  int len = std::atoi(arg[2]);
  int niter = std::atoi(arg[3]);

  using kv_t = Kokkos::View<double **, Kokkos::LayoutLeft, Kokkos::Serial>;

  kv_t A("A", len, nvec);
  kv_t B("B", len, nvec);
  kv_t C("C", len, nvec);

  double alpha = 2.1;
  double beta = 1.2;
  double gamma = 0.5;

  using time = Teuchos::TimeMonitor;
  for (int it = 0; it < niter; it++) {

    {
      time t(*time::getNewTimer("AXPBY -- KokkosBlas"));
      KokkosBlas::axpby(alpha, A, beta, B);
    }

    {
      time t(*time::getNewTimer("AXPBY -- Epetra-like"));
      for (int k = 0; k < nvec; k++) {
        const double *avec = Kokkos::subview(A, Kokkos::ALL(), k).data();
        double *bvec = Kokkos::subview(B, Kokkos::ALL(), k).data();
        for (int j = 0; j < len; j++) {
          bvec[j] = beta * bvec[j] + alpha * avec[j];
        }
      }
    }

    {
      time t(*time::getNewTimer("AXPBY -- hand-rolled"));
      for (int k = 0; k < nvec; k++) {
        for (int j = 0; j < len; j++) {
          B(j,k) = beta * B(j,k) + alpha * A(j,k);
        }
      }
    }

    {
      time t(*time::getNewTimer("Update -- KokkosBlas"));
      KokkosBlas::update(alpha, A, beta, B, gamma, C);
    }

    {
      time t(*time::getNewTimer("Update -- Epetra-like"));
      for (int k = 0; k < nvec; k++) {
        const double *avec = Kokkos::subview(A, Kokkos::ALL(), k).data();
        const double *bvec = Kokkos::subview(B, Kokkos::ALL(), k).data();
        double *cvec = Kokkos::subview(C, Kokkos::ALL(), k).data();
        for (int j = 0; j < len; j++) {
          cvec[j] = gamma * cvec[j] + beta * bvec[j] + alpha * avec[j];
        }
      }
    }

    {
      time t(*time::getNewTimer("Update -- hand-rolled"));
      for (int k = 0; k < nvec; k++) {
        for (int j = 0; j < len; j++) {
          C(j,k) = gamma * C(j,k) + beta * B(j,k) + alpha * A(j,k);
        }
      }
    }
  }

  Teuchos::TimeMonitor::summarize();
  MPI_Finalize();
  return 0;
}

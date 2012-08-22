#include "AnasaziTpetraAdapter.hpp"

#ifdef HAVE_ANASAZI_TPETRA_TIMERS

#include <Kokkos_SerialNode.hpp>
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<float, Tpetra::MultiVector<float, int, int, Kokkos::SerialNode   > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<float, Tpetra::MultiVector<float, int, int, Kokkos::SerialNode   > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::SerialNode   > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::SerialNode   > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<std::complex<float>, Tpetra::MultiVector<std::complex<float>, int, int, Kokkos::SerialNode   > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<std::complex<float>, Tpetra::MultiVector<std::complex<float>, int, int, Kokkos::SerialNode   > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<std::complex<double>, Tpetra::MultiVector<std::complex<double>, int, int, Kokkos::SerialNode   > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<std::complex<double>, Tpetra::MultiVector<std::complex<double>, int, int, Kokkos::SerialNode   > >::mvTimesMatAddMvTimer_ = Teuchos::null;

#ifdef HAVE_KOKKOSCLASSIC_TBB
#include <Kokkos_TBBNode.hpp>
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<float, Tpetra::MultiVector<float, int, int, Kokkos::TBBNode      > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<float, Tpetra::MultiVector<float, int, int, Kokkos::TBBNode      > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::TBBNode      > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::TBBNode      > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<std::complex<float>, Tpetra::MultiVector<std::complex<float>, int, int, Kokkos::TBBNode      > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<std::complex<float>, Tpetra::MultiVector<std::complex<float>, int, int, Kokkos::TBBNode      > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<std::complex<double>, Tpetra::MultiVector<std::complex<double>, int, int, Kokkos::TBBNode      > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<std::complex<double>, Tpetra::MultiVector<std::complex<double>, int, int, Kokkos::TBBNode      > >::mvTransMvTimer_ = Teuchos::null;
#endif

#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
#include <Kokkos_TPINode.hpp>
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<float, Tpetra::MultiVector<float, int, int, Kokkos::TPINode      > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<float, Tpetra::MultiVector<float, int, int, Kokkos::TPINode      > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::TPINode      > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::TPINode      > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<std::complex<float>, Tpetra::MultiVector<std::complex<float>, int, int, Kokkos::TPINode      > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<std::complex<float>, Tpetra::MultiVector<std::complex<float>, int, int, Kokkos::TPINode      > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<std::complex<double>, Tpetra::MultiVector<std::complex<double>, int, int, Kokkos::TPINode      > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<std::complex<double>, Tpetra::MultiVector<std::complex<double>, int, int, Kokkos::TPINode      > >::mvTimesMatAddMvTimer_ = Teuchos::null;
#endif

#ifdef HAVE_KOKKOSCLASSIC_THRUST
#include <Kokkos_ThrustGPUNode.hpp>
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<float, Tpetra::MultiVector<float, int, int, Kokkos::ThrustGPUNode> >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<float, Tpetra::MultiVector<float, int, int, Kokkos::ThrustGPUNode> >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::ThrustGPUNode> >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::ThrustGPUNode> >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<std::complex<float>, Tpetra::MultiVector<std::complex<float>, int, int, Kokkos::ThrustGPUNode> >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<std::complex<float>, Tpetra::MultiVector<std::complex<float>, int, int, Kokkos::ThrustGPUNode> >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<std::complex<double>, Tpetra::MultiVector<std::complex<double>, int, int, Kokkos::ThrustGPUNode> >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Anasazi::MultiVecTraits<std::complex<double>, Tpetra::MultiVector<std::complex<double>, int, int, Kokkos::ThrustGPUNode> >::mvTimesMatAddMvTimer_ = Teuchos::null;
#endif

#endif

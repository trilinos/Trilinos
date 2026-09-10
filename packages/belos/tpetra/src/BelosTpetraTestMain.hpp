// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOSTPETRATESTMAIN_HPP
#define BELOSTPETRATESTMAIN_HPP

#include "BelosTypes.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_MultiVector.hpp"

// In the test test.cpp we have:
//
// template <class ScalarType, class DM>
// int my_run (Teuchos::CommandLineProcessor& cmdp, int argc, char *argv[]) {
// ....
// }
//
// #define BELOS_MAIN_FUNC my_run
// #define BELOS_DEFAULT_SCALAR double
//
// #include "BelosTpetraTestMain.hpp"
//
// int main(int argc, char* argv[]) {
//   return common_main(argc, argv);
// }
//
// The function common_main defined below parses out the command line argument
// "--denseMatrix" and calls the corresponding instatiation of
// "BELOS_MAIN_FUNC". If BELOS_MAIN_FUNC or BELOS_DEFAULT_SCALAR are not defined
// they will default to the values below.

#ifndef BELOS_MAIN_FUNC
#define BELOS_MAIN_FUNC run
#endif

#ifndef BELOS_DEFAULT_SCALAR
#define BELOS_DEFAULT_SCALAR typename Tpetra::MultiVector<>::scalar_type
#endif

int common_main(int argc, char *argv[]) {
  using Ordinal = int;
  Teuchos::CommandLineProcessor clp(false, false);

  // CAG: The commented out logic allows to select the scalar type via
  //      command line flag. This requires building the solvers for all
  //      enabled scalar types. We should ETI the solvers before we
  //      uncomment this.

  // std::string scalarOption = "double";
  // if (std::is_same_v<BELOS_DEFAULT_SCALAR, double>) {
  //   scalarOption = "double";
  // } else if (std::is_same_v<BELOS_DEFAULT_SCALAR, std::complex<double>>) {
  //   scalarOption = "complex_double";
  // }
  // clp.setOption("scalar", &scalarOption, "Choice of scalar type");

  Belos::DenseMatrixType denseMatrixTypes[2] = {Belos::TeuchosSerialDenseMatrix, Belos::KokkosDualView};
  const char *denseMatrixNames[2] = {"Teuchos", "Kokkos"};

  Belos::DenseMatrixType denseMatrix = Belos::defaultDenseMatrixType;
  clp.setOption("denseMatrix", &denseMatrix, 2, denseMatrixTypes, denseMatrixNames,
                "Choice of dense matrix ( Teuchos | Kokkos )");

  switch (clp.parse(argc, argv, NULL)) {
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
    return EXIT_FAILURE;
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
    break;
  }

  clp.recogniseAllOptions(true);

  if (denseMatrix == Belos::TeuchosSerialDenseMatrix) {
    using DM = Teuchos::SerialDenseMatrix<Ordinal, BELOS_DEFAULT_SCALAR>;
    return BELOS_MAIN_FUNC<BELOS_DEFAULT_SCALAR, DM>(clp, argc, argv);
  } else if (denseMatrix == Belos::KokkosDualView) {
    using MV = Tpetra::MultiVector<BELOS_DEFAULT_SCALAR>;
    using DM = typename MV::wrapped_dual_view_type::DVT;
    return BELOS_MAIN_FUNC<BELOS_DEFAULT_SCALAR, DM>(clp, argc, argv);
  }

//   if (denseMatrix == TeuchosSerialDenseMatrix) {
//     if (scalarOption == "double") {
// #if defined(HAVE_TPETRA_INST_DOUBLE)
//       using DM = Teuchos::SerialDenseMatrix<Ordinal, double>;
//       return BELOS_MAIN_FUNC<double, DM>(clp, argc, argv);
// #else
//       std::cout << "Tpetra has not been instantiated for scalar type double\n";
//       return EXIT_FAILURE;
// #endif
//     } else if (scalarOption == "float") {
// #if defined(HAVE_TPETRA_INST_FLOAT)
//       using DM = Teuchos::SerialDenseMatrix<Ordinal, float>;
//       return BELOS_MAIN_FUNC<float, DM>(clp, argc, argv);
// #else
//       std::cout << "Tpetra has not been instantiated for scalar type float\n";
//       return EXIT_FAILURE;
// #endif
//     } else if (scalarOption == "complex_double") {
// #if defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
//       using DM = Teuchos::SerialDenseMatrix<Ordinal, std::complex<double>>;
//       return BELOS_MAIN_FUNC<std::complex<double>, DM>(clp, argc, argv);
// #else
//       std::cout << "Tpetra has not been instantiated for scalar type complex<double>\n";
//       return EXIT_FAILURE;
// #endif
//     } else if (scalarOption == "complex_float") {
//       // #if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
//       //   using DM = Teuchos::SerialDenseMatrix<Ordinal, std::complex<float>>;
//       //   return BELOS_MAIN_FUNC<std::complex<float>, DM>(clp, argc, argv);
//       // #endif
//     }
//   } else if ((denseMatrix == KokkosDualView)) {
//     if (scalarOption == "double") {
// #if defined(HAVE_TPETRA_INST_DOUBLE)
//       using IST = KokkosKernels::ArithTraits<double>::val_type;
//       using DM = Kokkos::DualView<IST **, Kokkos::LayoutLeft>;
//       return BELOS_MAIN_FUNC<double, DM>(clp, argc, argv);
// #else
//       std::cout << "Tpetra has not been instantiated for scalar type double\n";
//       return EXIT_FAILURE;
// #endif
//     } else if (scalarOption == "float") {
// #if defined(HAVE_TPETRA_INST_FLOAT)
//       using IST = KokkosKernels::ArithTraits<float>::val_type;
//       using DM = Kokkos::DualView<IST **, Kokkos::LayoutLeft>;
//       return BELOS_MAIN_FUNC<float, DM>(clp, argc, argv);
// #else
//       std::cout << "Tpetra has not been instantiated for scalar type float\n";
//       return EXIT_FAILURE;
// #endif
//     } else if (scalarOption == "complex_double") {
// #if defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
//       using IST =
//             KokkosKernels::ArithTraits<std::complex<double>>::val_type;
//       using DM = Kokkos::DualView<IST **, Kokkos::LayoutLeft>;
//       return BELOS_MAIN_FUNC<std::complex<double>, DM>(clp, argc, argv);
// #else
//       std::cout << "Tpetra has not been instantiated for scalar type complex<double>\n";
//       return EXIT_FAILURE;
// #endif
//     }
// // #if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
// //     else if (scalarOption == "complex_float") {
// //   using IST = KokkosKernels::ArithTraits<std::complex<float>>::val_type;
// //   using DM = Kokkos::DualView<IST **, Kokkos::LayoutLeft>;
// //   return BELOS_MAIN_FUNC<std::complex<float>, DM>(clp, argc, argv);
// // }
// // #endif
//   }
  return EXIT_FAILURE;
}

#endif

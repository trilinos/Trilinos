// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef GALERI_XPETRAPARAMETERS_DECL_HPP
#define GALERI_XPETRAPARAMETERS_DECL_HPP

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_VerboseObject.hpp>

#if defined(_MSC_VER) && !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

namespace Galeri {

namespace Xpetra {

// TODO nx/ny/nz == GO or global_size_t ? But what is the best to do?

template <typename GO>
class Parameters : public Teuchos::VerboseObject<Parameters<GO> >, public Teuchos::Describable {
 public:
  Parameters(Teuchos::CommandLineProcessor& clp, GO nx = 16, GO ny = -1, GO nz = -1, const std::string& matrixType = "Laplace1D",
             int keepBCs = 0, double stretchx = 1.0, double stretchy = 1.0, double stretchz = 1.0,
             double Kxx = 1.0, double Kxy = 0.0, double Kyy = 1.0, double dt = 1.0, const std::string& meshType = "tri",
             double h = 1.0, double delta = 0.0,
             int PMLXL = 0, int PMLXR = 0, int PMLYL = 0, int PMLYR = 0, int PMLZL = 0, int PMLZR = 0,
             double omega = 2.0 * M_PI, double shift = 0.5, GO mx = -1, GO my = -1, GO mz = -1, int model = 0,
             double lx = 1., double ly = 1., double lz = 1., double conv = 1., double diff = 1.);

  GO GetNumGlobalElements() const;

  const std::string& GetMatrixType() const;

  void check() const;

  Teuchos::ParameterList& GetParameterList() const;

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const;

  //! Print the object with some verbosity level to an FancyOStream object.
  void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = verbLevel_default) const;

  //@}

 private:
  // See Teuchos BUG 5249: https://software.sandia.gov/bugzilla/show_bug.cgi?id=5249
  mutable GO nx_, ny_, nz_;
  mutable GO mx_, my_, mz_;
  mutable double stretchx_, stretchy_, stretchz_;
  mutable double Kxx_, Kxy_, Kyy_;
  mutable double dt_;
  mutable std::string meshType_;

  std::string matrixType_;

  mutable int keepBCs_;

  mutable double h_;
  mutable double delta_;
  mutable int PMLx_left, PMLx_right;
  mutable int PMLy_left, PMLy_right;
  mutable int PMLz_left, PMLz_right;
  mutable double omega_;
  mutable double shift_;
  mutable int model_;

  mutable double lx_;
  mutable double ly_;
  mutable double lz_;
  mutable double conv_;
  mutable double diff_;

  // There is a major assumption here:
  // As soon as somebody call GetParameterList(), we freeze all other variables into the list,
  // and ignore them. This allows us to make the modification of the list from outside.
  mutable Teuchos::RCP<Teuchos::ParameterList> paramList_;
};

}  // namespace Xpetra
}  // namespace Galeri

#endif

// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef GALERI_XPETRAPARAMETERS_DEF_HPP
#define GALERI_XPETRAPARAMETERS_DEF_HPP

#include "Galeri_XpetraParameters_decl.hpp"

namespace Galeri {

namespace Xpetra {

template <class GO>
Parameters<GO>::Parameters(Teuchos::CommandLineProcessor& clp, GO nx, GO ny, GO nz, const std::string& matrixType,
                           int keepBCs, double stretchx, double stretchy, double stretchz,
                           double Kxx, double Kxy, double Kyy, double dt, const std::string& meshType,
                           double h, double delta,
                           int PMLXL, int PMLXR, int PMLYL, int PMLYR, int PMLZL, int PMLZR,
                           double omega, double shift, GO mx, GO my, GO mz, int model,
                           double lx, double ly, double lz, double conv, double diff)
  : nx_(nx)
  , ny_(ny)
  , nz_(nz)
  , mx_(mx)
  , my_(my)
  , mz_(mz)
  , stretchx_(stretchx)
  , stretchy_(stretchy)
  , stretchz_(stretchz)
  , Kxx_(Kxx)
  , Kxy_(Kxy)
  , Kyy_(Kyy)
  , dt_(dt)
  , meshType_(meshType)
  , matrixType_(matrixType)
  , keepBCs_(keepBCs)
  , h_(h)
  , delta_(delta)
  , PMLx_left(PMLXL)
  , PMLx_right(PMLXR)
  , PMLy_left(PMLYL)
  , PMLy_right(PMLYR)
  , PMLz_left(PMLZL)
  , PMLz_right(PMLZR)
  , omega_(omega)
  , shift_(shift)
  , model_(model)
  , lx_(lx)
  , ly_(ly)
  , lz_(lz)
  , conv_(conv)
  , diff_(diff) {
  clp.setOption("nx", &nx_, "mesh points in x-direction.");
  clp.setOption("ny", &ny_, "mesh points in y-direction.");
  clp.setOption("nz", &nz_, "mesh points in z-direction.");
  clp.setOption("mx", &mx_, "processors in x-direction.");
  clp.setOption("my", &my_, "processors in y-direction.");
  clp.setOption("mz", &mz_, "processors in z-direction.");
  clp.setOption("stretchx", &stretchx_, "stretch mesh in x-direction.");
  clp.setOption("stretchy", &stretchy_, "stretch mesh in y-direction.");
  clp.setOption("stretchz", &stretchz_, "stretch mesh in z-direction.");
  clp.setOption("Kxx", &Kxx_, "diffusion coefficient in xx-direction.");
  clp.setOption("Kxy", &Kxy_, "diffusion coefficient in xy-direction.");
  clp.setOption("Kyy", &Kyy_, "diffusion coefficient in yy-direction.");
  clp.setOption("dt", &dt_, "time step size.");
  clp.setOption("meshType", &meshType_, "meshType.");
  clp.setOption("keepBCs", &keepBCs_, "keep Dirichlet boundary rows in matrix (0=false,1=true)");
  clp.setOption("matrixType", &matrixType_, "matrix type: Laplace1D, Laplace2D, Laplace3D, ...");  // TODO: Star2D, numGlobalElements=...
  clp.setOption("h", &h_, "mesh width for uniform h");
  clp.setOption("delta", &delta_, "maximum PML damping value");
  clp.setOption("PMLx_left", &PMLx_left, "PML grid points in x-direction (left boundary)");
  clp.setOption("PMLx_right", &PMLx_right, "PML grid points in x-direction (right boundary)");
  clp.setOption("PMLy_left", &PMLy_left, "PML grid points in y-direction (left boundary)");
  clp.setOption("PMLy_right", &PMLy_right, "PML grid points in y-direction (right boundary)");
  clp.setOption("PMLz_left", &PMLz_left, "PML grid points in z-direction (left boundary)");
  clp.setOption("PMLz_right", &PMLz_right, "PML grid points in z-direction (right boundary)");
  clp.setOption("omega", &omega_, "angular frequency omega");
  clp.setOption("shift", &shift_, "complex frequency shift");
  clp.setOption("mx", &mx_, "processors in x-direction.");
  clp.setOption("my", &my_, "processors in y-direction.");
  clp.setOption("mz", &mz_, "processors in z-direction.");
  clp.setOption("model", &model_, "velocity model");
  clp.setOption("lx", &lx_, "length in x-direction");
  clp.setOption("ly", &ly_, "length in y-direction");
  clp.setOption("lz", &lz_, "length in z-direction");
  clp.setOption("convection", &conv_, "convection coefficient");
  clp.setOption("diffusion", &diff_, "diffusion coefficient");
}

template <class GO>
GO Parameters<GO>::GetNumGlobalElements() const {
  const Teuchos::ParameterList& pL = GetParameterList();

  const std::string& matrixType = pL.get<std::string>("matrixType");
  GO nx, ny, nz;
  if (pL.isType<int>("nx"))
    nx = Teuchos::as<GO>(pL.get<int>("nx"));
  else
    nx = pL.get<GO>("nx");
  if (pL.isType<int>("ny"))
    ny = Teuchos::as<GO>(pL.get<int>("ny"));
  else
    ny = pL.get<GO>("ny");
  if (pL.isType<int>("nz"))
    nz = Teuchos::as<GO>(pL.get<int>("nz"));
  else
    nz = pL.get<GO>("nz");

  GO numGlobalElements = -1;
  if (matrixType == "Laplace1D" || matrixType == "Helmholtz1D" || matrixType == "Identity")
    numGlobalElements = nx;

  else if (matrixType == "Laplace2D" ||
           matrixType == "Star2D" ||
           matrixType == "BigStar2D" ||
           matrixType == "AnisotropicDiffusion" ||
           matrixType == "Elasticity2D" ||
           matrixType == "Helmholtz2D" ||
           matrixType == "Recirc2D")
    numGlobalElements = nx * ny;

  else if (matrixType == "Laplace3D" ||
           matrixType == "Brick3D" ||
           matrixType == "Scalar3D_27Pt" ||
           matrixType == "HexFEM_LapStiff" ||
           matrixType == "HexFEM_Mass" ||
           matrixType == "BigStar2D" ||
           matrixType == "Elasticity3D" ||
           matrixType == "Helmholtz3D")
    numGlobalElements = nx * ny * nz;

  TEUCHOS_TEST_FOR_EXCEPTION(numGlobalElements < 0, std::runtime_error,
                             "Gallery: numGlobalElements < 0 (did you forget --ny (or --nz) for 2D (3D) problems?)");

  return numGlobalElements;
}

template <class GO>
const std::string& Parameters<GO>::GetMatrixType() const {
  const Teuchos::ParameterList& paramList = GetParameterList();
  return paramList.get<std::string>("matrixType");
}

template <class GO>
void Parameters<GO>::check() const {}

template <class GO>
Teuchos::ParameterList& Parameters<GO>::GetParameterList() const {
  if (!paramList_.is_null())
    return *paramList_;

  paramList_ = rcp(new Teuchos::ParameterList());

  paramList_->set("nx", nx_);
  paramList_->set("ny", ny_);
  paramList_->set("nz", nz_);
  paramList_->set("mx", mx_);
  paramList_->set("my", my_);
  paramList_->set("mz", mz_);
  paramList_->set("model", model_);
  paramList_->set("stretchx", stretchx_);
  paramList_->set("stretchy", stretchy_);
  paramList_->set("stretchz", stretchz_);
  paramList_->set("Kxx", Kxx_);
  paramList_->set("Kxy", Kxy_);
  paramList_->set("Kyy", Kyy_);
  paramList_->set("dt", dt_);
  paramList_->set("meshType", meshType_);
  paramList_->set("keepBCs", static_cast<bool>(keepBCs_));
  paramList_->set("matrixType", matrixType_);
  paramList_->set("h", h_);
  paramList_->set("delta", delta_);
  paramList_->set("PMLx_left", PMLx_left);
  paramList_->set("PMLx_right", PMLx_right);
  paramList_->set("PMLy_left", PMLy_left);
  paramList_->set("PMLy_right", PMLy_right);
  paramList_->set("PMLz_left", PMLz_left);
  paramList_->set("PMLz_right", PMLz_right);
  paramList_->set("omega", omega_);
  paramList_->set("shift", shift_);
  paramList_->set("lx", lx_);
  paramList_->set("ly", ly_);
  paramList_->set("lz", lz_);
  paramList_->set("convection", conv_);
  paramList_->set("diffusion", diff_);

  check();

  return *paramList_;
}

template <class GO>
std::string Parameters<GO>::description() const {
  std::ostringstream out;
  out << Teuchos::Describable::description();
  out << "{type = " << GetMatrixType() << ", size = " << GetNumGlobalElements() << "} ";

  return out.str();
}

template <class GO>
void Parameters<GO>::describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
  int vl = (verbLevel == Teuchos::VERB_DEFAULT) ? Teuchos::VERB_LOW : verbLevel;
  if (vl == Teuchos::VERB_NONE)
    return;

  if (vl == Teuchos::VERB_LOW)
    out << description() << std::endl;
  else
    out << Teuchos::Describable::description() << std::endl;

  if (vl == Teuchos::VERB_MEDIUM || vl == Teuchos::VERB_HIGH || vl == Teuchos::VERB_EXTREME) {
    Teuchos::OSTab tab1(out);

    const Teuchos::ParameterList& paramList = GetParameterList();
    std::string matrixType                  = paramList.get<std::string>("matrixType");
    GO nx, ny, nz;
    if (paramList.isType<int>("nx"))
      nx = Teuchos::as<GO>(paramList.get<int>("nx"));
    else
      nx = paramList.get<GO>("nx");
    if (paramList.isType<int>("ny"))
      ny = Teuchos::as<GO>(paramList.get<int>("ny"));
    else
      ny = paramList.get<GO>("ny");
    if (paramList.isType<int>("nz"))
      nz = Teuchos::as<GO>(paramList.get<int>("nz"));
    else
      nz = paramList.get<GO>("nz");

    out << "Matrix type: " << matrixType << std::endl
        << "Problem size: " << GetNumGlobalElements();

    if (matrixType == "Laplace2D" || matrixType == "AnisotropicDiffusion" || matrixType == "Elasticity2D" || matrixType == "Helmholtz2D")
      out << " (" << nx << "x" << ny << ")";
    else if (matrixType == "Laplace3D" || matrixType == "Elasticity3D" || matrixType == "Helmholtz3D")
      out << " (" << nx << "x" << ny << "x" << nz << ")";

    out << std::endl;

    if (matrixType == "AnisotropicDiffusion") {
      double Kxx = paramList.get<double>("Kxx");
      double Kxy = paramList.get<double>("Kxy");
      double Kyy = paramList.get<double>("Kyy");
      double dt  = paramList.get<double>("dt");
      out << "K  = [[ " << Kxx << ", " << Kxy << " ], [ " << Kxy << ", " << Kyy << " ]]" << std::endl;
      out << "dt = " << dt << std::endl;
    }
  }
}

}  // namespace Xpetra
}  // namespace Galeri

#define GALERI_XPETRAPARAMETERS_INSTANT(GO) \
  template class Galeri::Xpetra::Parameters<GO>;

#endif

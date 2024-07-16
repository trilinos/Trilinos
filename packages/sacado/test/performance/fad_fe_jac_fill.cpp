// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "fe_jac_fill_funcs.hpp"

#include "Fad/fad.h"
#include "TinyFadET/tfad.h"

// A performance test that computes a finite-element-like Jacobian using
// several Fad variants

void FAD::error(const char *msg) {
  std::cout << msg << std::endl;
}

int main(int argc, char* argv[]) {
  int ierr = 0;

  try {
    double t, ta, tr;
    int p = 2;
    int w = p+7;

    // Maximum number of derivative components for SLFad
    const int slfad_max = 130;

    // Set up command line options
    Teuchos::CommandLineProcessor clp;
    clp.setDocString("This program tests the speed of various forward mode AD implementations for a finite-element-like Jacobian fill");
    int num_nodes = 100000;
    int num_eqns = 2;
    int rt = 0;
    clp.setOption("n", &num_nodes, "Number of nodes");
    clp.setOption("p", &num_eqns, "Number of equations");
    clp.setOption("rt", &rt, "Include ADOL-C retaping test");

    // Parse options
    Teuchos::CommandLineProcessor::EParseCommandLineReturn
      parseReturn= clp.parse(argc, argv);
    if(parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
      return 1;

    double mesh_spacing = 1.0 / static_cast<double>(num_nodes - 1);

    std::cout.setf(std::ios::scientific);
    std::cout.precision(p);
    std::cout << "num_nodes =  " << num_nodes
        << ", num_eqns = " << num_eqns << ":  " << std::endl
        << "               " << "   Time   " << "\t"<< "Time/Analytic" << "\t"
        << "Time/(2*p*Residual)" << std::endl;

    ta = 1.0;
    tr = 1.0;

    tr = residual_fill(num_nodes, num_eqns, mesh_spacing);

    ta = analytic_jac_fill(num_nodes, num_eqns, mesh_spacing);
    std::cout << "Analytic:      " << std::setw(w) << ta << "\t" << std::setw(w) << ta/ta << "\t" << std::setw(w) << ta/(2.0*num_eqns*tr) << std::endl;

#ifdef HAVE_ADOLC
#ifndef ADOLC_TAPELESS
    t = adolc_jac_fill(num_nodes, num_eqns, mesh_spacing);
    std::cout << "ADOL-C:        " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;

    if (rt != 0) {
      t = adolc_retape_jac_fill(num_nodes, num_eqns, mesh_spacing);
      std::cout << "ADOL-C(rt):  " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }

#else
    t = adolc_tapeless_jac_fill(num_nodes, num_eqns, mesh_spacing);
    std::cout << "ADOL-C(tl):    " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
#endif
#endif

#ifdef HAVE_ADIC
    t = adic_jac_fill(num_nodes, num_eqns, mesh_spacing);
    std::cout << "ADIC:          " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
#endif

    if (num_eqns*2 == 4) {
      t = fad_jac_fill< FAD::TFad<16,double> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "TFad:          " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }
    else if (num_eqns*2 == 16) {
      t = fad_jac_fill< FAD::TFad<16,double> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "TFad:          " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }
    else if (num_eqns*2 == 32) {
      t = fad_jac_fill< FAD::TFad<32,double> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "TFad:          " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }
    else if (num_eqns*2 == 64) {
      t = fad_jac_fill< FAD::TFad<64,double> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "TFad:          " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }

    t = fad_jac_fill< FAD::Fad<double> >(num_nodes, num_eqns, mesh_spacing);
    std::cout << "Fad:           " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;

    if (num_eqns*2 == 4) {
      t = fad_jac_fill< Sacado::Fad::SFad<double,4> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "SFad:          " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }
    else if (num_eqns*2 == 16) {
      t = fad_jac_fill< Sacado::Fad::SFad<double,16> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "SFad:          " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }
    else if (num_eqns*2 == 32) {
      t = fad_jac_fill< Sacado::Fad::SFad<double,32> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "SFad:          " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }
    else if (num_eqns*2 == 64) {
      t = fad_jac_fill< Sacado::Fad::SFad<double,64> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "SFad:          " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }

    if (num_eqns*2 < slfad_max) {
      t = fad_jac_fill< Sacado::Fad::SLFad<double,slfad_max> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "SLFad:         " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }

    t = fad_jac_fill< Sacado::Fad::DFad<double> >(num_nodes, num_eqns, mesh_spacing);
    std::cout << "DFad:          " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;

    t = fad_jac_fill< Sacado::Fad::SimpleFad<double> >(num_nodes, num_eqns, mesh_spacing);
    std::cout << "SimpleFad:     " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;

    if (num_eqns*2 == 4) {
      t = fad_jac_fill< Sacado::ELRFad::SFad<double,4> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "ELRSFad:       " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }
    else if (num_eqns*2 == 16) {
      t = fad_jac_fill< Sacado::ELRFad::SFad<double,16> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "ELRSFad:       " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }
    else if (num_eqns*2 == 32) {
      t = fad_jac_fill< Sacado::ELRFad::SFad<double,32> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "ELRSFad:       " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }
    else if (num_eqns*2 == 64) {
      t = fad_jac_fill< Sacado::ELRFad::SFad<double,64> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "ELRSFad:       " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }

    if (num_eqns*2 < slfad_max) {
      t = fad_jac_fill< Sacado::ELRFad::SLFad<double,slfad_max> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "ELRSLFad:      " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }

    t = fad_jac_fill< Sacado::ELRFad::DFad<double> >(num_nodes, num_eqns, mesh_spacing);
    std::cout << "ELRDFad:       " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;

    if (num_eqns*2 == 4) {
      t = fad_jac_fill< Sacado::CacheFad::SFad<double,4> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "CacheSFad:     " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }
    else if (num_eqns*2 == 16) {
      t = fad_jac_fill< Sacado::CacheFad::SFad<double,16> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "CacheSFad:     " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }
    else if (num_eqns*2 == 32) {
      t = fad_jac_fill< Sacado::CacheFad::SFad<double,32> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "CacheSFad:     " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }
    else if (num_eqns*2 == 64) {
      t = fad_jac_fill< Sacado::CacheFad::SFad<double,64> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "CacheSFad:     " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }

    if (num_eqns*2 < slfad_max) {
      t = fad_jac_fill< Sacado::CacheFad::SLFad<double,slfad_max> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "CacheSLFad:    " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }

    t = fad_jac_fill< Sacado::CacheFad::DFad<double> >(num_nodes, num_eqns, mesh_spacing);
    std::cout << "CacheFad:      " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;

    if (num_eqns*2 == 4) {
      t = fad_jac_fill< Sacado::ELRCacheFad::SFad<double,4> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "ELRCacheSFad:  " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }
    else if (num_eqns*2 == 16) {
      t = fad_jac_fill< Sacado::ELRCacheFad::SFad<double,16> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "ELRCacheSFad:  " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }
    else if (num_eqns*2 == 32) {
      t = fad_jac_fill< Sacado::ELRCacheFad::SFad<double,32> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "ELRCacheSFad:  " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }
    else if (num_eqns*2 == 64) {
      t = fad_jac_fill< Sacado::ELRCacheFad::SFad<double,64> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "ELRCacheSFad:  " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }

    if (num_eqns*2 < slfad_max) {
      t = fad_jac_fill< Sacado::ELRCacheFad::SLFad<double,slfad_max> >(num_nodes, num_eqns, mesh_spacing);
      std::cout << "ELRCacheSLFad: " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;
    }

    t = fad_jac_fill< Sacado::ELRCacheFad::DFad<double> >(num_nodes, num_eqns, mesh_spacing);
    std::cout << "ELRCacheFad:   " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;

    t = fad_jac_fill< Sacado::Fad::DVFad<double> >(num_nodes, num_eqns, mesh_spacing);
    std::cout << "DVFad:         " << std::setw(w) << t << "\t" << std::setw(w) << t/ta << "\t" << std::setw(w) << t/(2.0*num_eqns*tr) << std::endl;

  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    ierr = 1;
  }
  catch (const char *s) {
    std::cout << s << std::endl;
    ierr = 1;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" << std::endl;
    ierr = 1;
  }

  return ierr;
}

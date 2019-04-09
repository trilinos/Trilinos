// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#if HAVE_OPENNURBS

#include <percept/mesh/geometry/stk_geom/LocalCubicSplineFit.hpp>

  namespace geom {

    void LocalCubicSplineFit::
    fit_internal(int n, Vectors2D& Q, Vectors2D& T)
    {
      const bool debug_print = SplineFit::s_debug_print;
      //std::cout << "LocalCubicSplineFit::fit_internal Q= " << Q << std::endl;

      std::vector<double>& U = m_U;
      Vectors2D& CV = m_CV;
      m_n = n;

      Vectors2D Pk1(n), Pk2(n);
      std::vector<double> u(n+1), alpha(n+1);
      u[0] = 0.0;
      for (int k=0; k <= n-1; ++k)
        {
          // eq (9.50)
          Vector2D Tkp = T[k+1];
          Vector2D Tk = T[k];
          if (m_isCorner[k] && m_isCorner[k+1])
            {
              Tkp = Q[k+1] - Q[k];
              double dTkp = Tkp.Length();
              if (dTkp < 1.e-12) throw std::runtime_error("can't normalize Tkp vector");
              Tkp = Tkp / dTkp;
              Tk = Tkp;
            }
          else if (m_isCorner[k+1])
            {
              Tkp = Q[k+1] - Q[k];
              double dTkp = Tkp.Length();
              if (dTkp < 1.e-12) throw std::runtime_error("can't normalize Tkp vector");
              Tkp = Tkp / dTkp;
              Tk = Tkp;
            }
          else if (m_isCorner[k])
            {
              Tk = Q[k+1] - Q[k];
              double dTk = Tk.Length();
              if (dTk < 1.e-12) throw std::runtime_error("can't normalize Tk vector");
              Tk = Tk / dTk;
              Tkp = Tk;
            }

          double a = 16.0 - (Tkp + Tk).LengthSquared();
          double b = 12.0 * (Q[k+1] - Q[k])*(Tkp + Tk);
          double c = -36.0*(Q[k+1] - Q[k]).LengthSquared();
          double disc = b*b - 4.0*a*c;
          if (disc < 0 || a < 1.e-10)
            throw std::runtime_error("disc < 0 || a <= 0");
          alpha[k] = (-b + std::sqrt(disc))/(2.0*a);
          if (alpha[k] < -1.e-10)
            throw std::runtime_error("alpha < 0");
          DPRINTLN2(k,alpha[k]);
          // eq (9.47)
          // if (m_isCorner[k] || m_isCorner[k+1])
          //   alpha[k] = 1.0/(Q[k+1]-Q[k]).Length();
          Pk1[k] = Q[k] + alpha[k]*Tk/3.0;
          Pk2[k] = Q[k+1] - alpha[k]*Tkp/3.0;

          // eq (9.52)
          u[k+1] = u[k] + 3.0*(Pk2[k] - Pk1[k]).Length();
          DPRINTLN2(k,u[k+1]);
        }
      int order = 4; // degree=3, order = degree+1
      int offset = order-1; // note: P&T convention is order, not order-1 as OpenNURBS seems to require

      // fixup for corners
      int ncorner=0;
      for (int k=1; k <= n-1; k++)
        {
          if (m_isCorner[k])
            {
              ++ncorner;
            }
        }
      DPRINTLN(ncorner);

      int nU = 2*offset+2*(n-1) + ncorner;

      U.resize(nU);
      for (int k=0; k < offset; k++)
        {
          U[k] = 0.0;
          U[nU-1-k] = 1.0;
        }
      int ki=offset;
      for (int k=1; k <= n-1; k++)
        {
          // eq (9.54)
          double uk=u[k]/u[n];
          U[ki++] = uk;
          U[ki++] = uk;
          if (m_isCorner[k])
            U[ki++] = uk;
        }

      int nCV = 2 + 2*n + ncorner;
      m_CV.resize(nCV);
      // eq (9.53)
      CV[0] = Q[0];
      ki=1;
      int inc=0;
      for (int k=0; k <= n-1; k++)
        {
          CV[ki++] = Pk1[k];
          CV[ki++] = Pk2[k];
          if (k < n-1 && m_isCorner[k+1])
            {
              ++inc;
              CV[ki++] = Q[k+1];
            }
        }
      if (debug_print) std::cout << "ki= " << ki << " CV.size= " << CV.size() << " n= " << n <<  " inc= " << inc << " ncorner= " << ncorner <<std::endl;
      CV[ki] = Q[n];
      if (debug_print) std::cout << "\nLocalCubicSplineFit::fit_internal" << std::endl;
      DPRINTLN(Q);
      DPRINTLN(T);
      DPRINTLN(CV);
      DPRINTLN(U);

    }

  }

#endif

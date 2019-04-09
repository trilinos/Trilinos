// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#if HAVE_OPENNURBS

#include <percept/mesh/geometry/stk_geom/BSplineFit.hpp>

  namespace geom {


    double BSplineFit::alpha_five_point(int k, Vectors2D& d, int doff, bool allow_corner)
    {
      double dkp1kp2 = std::fabs(ON_CrossProduct(d[k+1+doff],d[k+2+doff])[2]);
      double dkm1k = std::fabs(ON_CrossProduct(d[k-1+doff],d[k+doff])[2]);
      double den = dkm1k+dkp1kp2;
      if (den < 1.e-8)   // FIXME tol
        {
          if (allow_corner)
            {
              std::cout << "tmp srk FOUND CORNER at k= " << k << std::endl;
              m_isCorner[k] = 1;
            }
          return 1.0;  // implies a corner - for smoothing the corner, return 0.5
        }
      if (m_isCorner[k])
        return 1.0;
      return dkm1k/den;
    }

    ON_Curve * BSplineFit::
    fit(Vectors2D& points_in)
    {
      const bool debug_print = SplineFit::s_debug_print;
      int n_i = points_in.size();
      //int n_stencil = option; // option.n_stencil;
      int n = n_i - 1; // 0...n define the points

      /// using Piegl and Tiller "The NURBS Book" notation and equation numbers
      Vectors2D Q = points_in;
      m_Q = Q;
      if (m_isCorner.size() == 0)
        m_isCorner.assign(n+1,0);

      double arclen = 0.0;
      for (int k=1; k <= n; k++)
        {
          // eq (9.3)
          //u[k] = double(k)/double(n+1);
          // eq (9.4,5)
          arclen += (Q[k] - Q[k-1]).Length();
        }

      // decide if points are closed
      m_isClosed = false;
      double reltol = 1.e-8;
      double tol = reltol*arclen;
      Vector2D diff = Q[0];
      diff -= Q[n];
      if (diff.Length() < tol)
        {
          m_isClosed = true;
        }

      Vectors2D q(n+1), d(n+1), D(n+1), V(n+1), T(n+1);
      std::vector<double> u(n+1), alpha(n+1);
      u[0] = 0.0;
      u[n] = 1.0;

      for (int k=1; k <= n-1; k++)
        {
          // eq (9.4,5)
          u[k] = u[k-1] + (Q[k] - Q[k-1]).Length() / arclen;
          DPRINTLN2(k,u[k]);
        }

      for (int k=1; k <= n; k++)
        {
          // eq (9.28)
          q[k] = Q[k] - Q[k-1];
          d[k] = q[k]/(u[k] - u[k-1]);
        }

      if (m_option == ThreePoint)
        {
          /**
           *        0   1   2   ...  n-1   n   n+1   n+2   n+3
           *
           *                               0    1     2   ....  periodic overlap
           *
           *   ThreePoint:
           *      u[n+1] defined as u[n]+(u[1]-u[0]);
           *      q[n+1] = q[1]
           *      d[n+1] = d[1]
           *      allows definition of alpha[n] and T[n], then
           *      T0 == Tn
           *
           *   FivePoint:
           *      u[n+1] = u[n] + (u[1]-u[0]);
           *      u[n+2] = u[n+1] + (u[2]-u[1]);
           *      u[n+3] = u[n+2] + (u[3]-u[2]);
           *      alpha[n, n+1, n+2] can now be defined, and T[n, n+1, n+2]
           *      then T[0] = T[n], T[1] = T[n+1], T[2] = T[n+2]
           *
           */

          int N = n-1;
          if (m_isPeriodic)
            {
              double unp1 = u[n] + (u[1]-u[0]);
              u.push_back(unp1);
              q.push_back(q[1]);
              d.push_back(d[1]);
              N = n;
            }
          DPRINTLN(u);
          for (int k=1; k <= N; k++)
            {
              // eq (9.30)
              alpha[k] = (u[k] - u[k-1]) / (u[k+1] - u[k-1]);
              if (m_isCorner[k])
                alpha[k]=1.0;
              // eq (9.28,29)
              D[k] = (1.0-alpha[k])*d[k] + alpha[k]*d[k+1];
              V[k] = (1.0-alpha[k])*q[k] + alpha[k]*q[k+1];
              double d0l =  D[k].Length();
              if (d0l < 1.e-12) throw std::runtime_error("can't normalize T vector");
              DPRINTLN2(k,alpha[k]);
              // eq (9.29)
#if 1
              T[k] = D[k] / d0l;
#else
              Vector2D Tkp = d[k+1]/d[k+1].Length();
              Vector2D Tk = d[k]/d[k].Length();
              T[k] = 0.5*(Tkp+Tk);
              if (m_isCorner[k])
                T[k] = Tk;
              T[k] = T[k]/T[k].Length();
#endif
              DPRINTLN2(k,T[k]);
            }
          if (m_isPeriodic)
            {
              T[0] = T[n];
            }
          else
            {
              // endpoints
              // eq (9.32)
              D[0] = 2*d[1] - D[1];
              D[n] = 2*d[n] - D[n-1];
              V[0] = 2*q[1] - V[1];
              V[n] = 2*q[n] - V[n-1];
              double d0l = D[0].Length();
              double dnl = D[n].Length();
              if (d0l < 1.e-12) throw std::runtime_error("can't normalize T vector[0]");
              if (dnl < 1.e-12) throw std::runtime_error("can't normalize T vector[n]");
              T[0] = D[0]/d0l;
              T[n] = D[n]/dnl;

#if 0
              T[0] = (Q[1] - Q[0]);
              T[0] = T[0]/ T[0].Length();
              T[n] = (Q[n] - Q[n-1]);
              T[n] = T[n]/ T[n].Length();
#endif
              DPRINTLN(d[1]);
              DPRINTLN(D[1]);
              DPRINTLN(T[0]);
              DPRINTLN(T[n]);
            }
        }

      //////////////////////////////////////////////////////////////////////////////////
      ///  FivePoint
      //////////////////////////////////////////////////////////////////////////////////
      else if (m_option == FivePoint)
        {
          int doff = 0;
          int N0 = 2;
          int N = n-2;
          alpha.resize(n+1+3);
          if (m_isPeriodic)
            {
              double unp1 = u[n] + (u[1]-u[0]);
              double unp2 = unp1 + (u[2]-u[1]);
              double unp3 = unp2 + (u[3]-u[2]);

              u.push_back(unp1);
              u.push_back(unp2);
              u.push_back(unp3);

              q.push_back(q[1]);
              q.push_back(q[2]);
              q.push_back(q[3]);

              d.push_back(d[1]);
              d.push_back(d[2]);
              d.push_back(d[3]);

              Vector2D v2;
              T.push_back(v2);
              T.push_back(v2);
              T.push_back(v2);

              N = n;
            }
          else
            {
              // -1,...,n+2
              Vectors2D dx0;
              dx0.push_back(Vector2D());
              dx0.insert(dx0.end(), d.begin(), d.end());
              dx0.push_back(Vector2D());
              dx0.push_back(Vector2D());
              DPRINTLN(d);
              d = dx0;
              // endpoints
              // eq (9.33)
              doff = 1;
              Vector2D *dx = &d[1];
              DPRINTLN(d);
              dx[0]   = 2*dx[1]   - dx[2];
              dx[-1]  = 2*dx[0]   - dx[1];
              dx[n+1] = 2*dx[n]   - dx[n-1];
              dx[n+2] = 2*dx[n+1] - dx[n];
              DPRINTLN(d);
            }
          DPRINTLN(u);
          for (int k = N0; k <= N; k++)
            {
              // eq (9.31)
              alpha[k] = alpha_five_point(k, d, doff);
              double akm1 = 1.0-alpha[k];
              // eq (9.29)
              //Vector2D Dk = akm1*d[k] + alpha[k]*d[k+1];
              Vector2D Dk = akm1*d[k+doff] + alpha[k]*d[k+1+doff];
              double dkl =  Dk.Length();
              if (dkl < 1.e-12) throw std::runtime_error("can't normalize T vector");
              DPRINTLN2(k,alpha[k]);
              // eq (9.29)
              T[k] = Dk / dkl;
              DPRINTLN2(k,T[k]);
            }

          if (m_isPeriodic)
            {
              T[0] = T[n];
              T[1] = T[n+1];
              T[2] = T[n+2];
              T.resize(n+1);
            }
          else
            {
              for (int k=0; k < 2; k++)
                {
                  alpha[k] = alpha_five_point(k, d, doff, false);
                  D[k] = (1-alpha[k])*d[k+doff] + alpha[k]*d[k+1+doff];
                  double dkl = D[k].Length();
                  if (dkl < 1.e-12) throw std::runtime_error("can't normalize T vector[0]");
                  T[k] = D[k]/dkl;
                }
              for (int k=n-1; k <= n; k++)
                {
                  alpha[k] = alpha_five_point(k, d, doff, false);
                  D[k] = (1-alpha[k])*d[k+doff] + alpha[k]*d[k+1+doff];
                  double dkl = D[k].Length();
                  if (dkl < 1.e-12) throw std::runtime_error("can't normalize T vector[0]");
                  T[k] = D[k]/dkl;
                }
            }
          DPRINTLN(alpha);
        }
      if (0 && m_isPeriodic && m_isClosed)
        {
          Vector2D Tave = 0.5*(T[0] + T[n]);
          T[0] = Tave;
          T[n] = Tave;
        }
      fit_internal(n, Q, T);
      return create();
    }

    ON_Curve * BSplineFit::create()
    {
      const bool debug_print = SplineFit::s_debug_print;
      const int dimension = 2;
      const int order = 4;

      ON_NurbsCurve* pNurbsCurve =
        ON_NurbsCurve::New(
                           dimension,         // int dimension,
                           false,             // ON_BOOL32 bIsRational,
                           order,             // int order = degree+1
                           (int)m_CV.size()   // int cv_count
                           );
      if (!pNurbsCurve) throw std::runtime_error("failed in ON_NurbsCurve");

      for (int k = 0; k < (int)m_CV.size(); k++)
        {
          Point3D v;
          v.Zero();
          v[0] = m_CV[k][0];
          v[1] = m_CV[k][1];
          ON_BOOL32 success = pNurbsCurve->SetCV(k, v);
          if (!success) {
            PRINTLN(k);
            throw std::runtime_error("failed in SetCV");
          }
        }
      DPRINTLN(m_CV.size());
      DPRINTLN(m_U.size());
      DPRINTLN(pNurbsCurve->KnotCount());

      for (int k = 0; k < (int)m_U.size(); k++)
        {
          ON_BOOL32 success = pNurbsCurve->SetKnot(k, m_U[k]);
          if (!success) {
            PRINTLN(k);
            throw std::runtime_error("failed in SetKnot");
          }
        }
      m_curve = pNurbsCurve;
      return pNurbsCurve;
    }


    void BSplineFit::print()
    {
      std::cout << "BSplineFit::print: getIsPeriodic= " << getIsPeriodic() << " getIsClosed= " << getIsClosed() << std::endl;
      std::cout << "n = " << m_n << " cv count= " << m_CV.size() << " U count = " << m_U.size() << std::endl;
      std::cout << "Q= " << m_Q << std::endl;
      std::cout << "CV= " << m_CV << std::endl;
      std::cout << "U= " << m_U << std::endl;
      if (m_curve)
        {
          PRINTLN(m_curve->IsValid());
          PRINTLN(m_curve->IsRational());
          PRINTLN(m_curve->CVCount());
          PRINTLN(m_curve->CVSize());
          PRINTLN(m_curve->Order());
          PRINTLN(m_curve->KnotCount());
          ON_TextLog log;
          m_curve->Dump(log);
        }
    }

  }

#endif

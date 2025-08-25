// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef BSplineFit_hpp
#define BSplineFit_hpp

#if HAVE_OPENNURBS

#include <percept/mesh/geometry/stk_geom/SplineFit.hpp>

  namespace geom {

    /// using Piegl and Tiller "The NURBS Book, 2nd Edition" notation and equation numbers

    /// base class for "local" spline fit methods creating a B-spline
    /// (note: local methods depend only on a few points around each point
    /// and don't require a tridiagonal solve)
    class BSplineFit : public SplineFit
    {
    public:
      /** Choice for how many points around a given input point to use to calculate tangent vectors
       *
       * FivePoint has the advantage or recognizing straight lines (3 points) and corners with
       * two straight lines eminating from it.  Also, known corners can be identified (in a later
       * version of this code) by the user and enforced with the 5pt scheme.
       *
       * Three point is less susceptible to noisy data.
       */
      enum Option { ThreePoint = 3, FivePoint = 5 };

      BSplineFit(Option option = ThreePoint) : m_option(option), m_n(0), m_isRational(false),
                     m_isPeriodic(false), m_isClosed(false), m_curve(0) {}

      // periodic means it is closed and we want to treat the 0 and n point as a smooth point
      void setIsPeriodic(bool per) { m_isPeriodic=per; }
      bool getIsPeriodic() { return m_isPeriodic; }
      void setIsCorner(const std::vector<int>& isCorner) { m_isCorner = isCorner; }
      //void setIsCorner(int index, bool isCorner) { if (isCorner) m_isCorner[index] = 1; }

      bool getIsClosed() { return m_isClosed; }

      /// create an OpenNURBS curve that fits the given input points
      virtual ON_Curve * fit(Vectors2D& input) override;
      void print();

      const Vectors2D& Q() { return m_Q; }

      const std::vector<double>& U() { return m_U; }
      const Vectors2D& CV() { return m_CV; }

    protected:
      // computes m_CV, m_U
      virtual void fit_internal(int n, Vectors2D& Q, Vectors2D& T) = 0;

      // given m_CV, m_U, gives the ON_NurbsCurve
      ON_Curve * create();

      double alpha_five_point(int k, Vectors2D& d, int doff, bool allow_corner=true);

    protected:
      Option m_option;
      // control vertices, knots
      Vectors2D m_CV;
      Vectors2D m_Q;
      std::vector<int> m_isCorner; // intentional use of int to avoid packed bool
      std::vector<double> m_U;
      int m_n;
      bool m_isRational;
      bool m_isPeriodic;
      bool m_isClosed;
      ON_NurbsCurve *m_curve;
    };
  }

#endif

#endif

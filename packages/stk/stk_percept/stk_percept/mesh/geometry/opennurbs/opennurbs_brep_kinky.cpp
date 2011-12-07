/* $NoKeywords: $ */
/*
//
// Copyright (c) 1993-2007 Robert McNeel & Associates. All rights reserved.
// Rhinoceros is a registered trademark of Robert McNeel & Assoicates.
//
// THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY.
// ALL IMPLIED WARRANTIES OF FITNESS FOR ANY PARTICULAR PURPOSE AND OF
// MERCHANTABILITY ARE HEREBY DISCLAIMED.
//				
// For complete openNURBS copyright information see <http://www.opennurbs.org>.
//
////////////////////////////////////////////////////////////////
*/

#include "opennurbs.h"

bool ON_Brep::SplitKinkyFaces( 
        double kink_tol_radians,
        bool bCompactIfNeeded
        )
{
  bool rc = true;
  // see if splitting is required
  const int ecount = m_E.Count();
  const int fcount = m_F.Count();
  for (int j=0; j<fcount; j++)
  {
    if ( !SplitKinkyFace(j,kink_tol_radians) )
      rc = false;
  }
  if (bCompactIfNeeded && ( fcount != m_F.Count() || ecount != m_E.Count()) )
  {
    Compact();
  }
  return true;
}


bool ON_Brep::SplitKinkyFace( 
  int,   // face_index - formal parameter intentionally ignored in this virtual function
  double // kink_tol_radians - formal parameter intentionally ignored in this virtual function
  )
{
  // works in RHino SDK - not part of free opennurbs
  return false;
}

bool ON_Brep::SplitKinkyEdge( 
  int edge_index, 
  double kink_tol_radians
  )
{
  // Default kink_tol_radians MUST BE ON_PI/180.0.
  //
  // The default kink tol must be kept in sync with the default for 
  // TL_Brep::SplitKinkyFace() and ON_Brep::SplitKinkyFace().
  // See comments in TL_Brep::SplitKinkyFace() for more details.

  bool rc = true;
  if (kink_tol_radians < ON_ZERO_TOLERANCE) 
    kink_tol_radians = ON_ZERO_TOLERANCE;
  else if (kink_tol_radians > ON_PI - ON_ZERO_TOLERANCE) 
    kink_tol_radians = ON_PI - ON_ZERO_TOLERANCE;
  double atol = cos(kink_tol_radians);
  if (edge_index < 0 || edge_index >= m_E.Count()) 
    return false;
  ON_BrepEdge& E = m_E[edge_index];
  if (E.m_c3i < 0) 
    return false;
  ON_SimpleArray<double> split_t(4);
  double t0 = E.Domain()[0];
  int hint = 0;
  ON_Curve* curve = m_C3[E.m_c3i];
  if (!curve) return false;
  int scount = curve->SpanCount();
  while (split_t.Count() < scount){
    double t;
    if (!E.GetNextDiscontinuity(ON::G1_continuous, t0, E.Domain()[1], 
      &t, &hint, NULL, atol)) break;
    split_t.Append(t);
    t0 = t;
  }
  if (split_t.Count() >= scount) 
    return false;

  if (split_t.Count() == 0) 
    return true;//no kinks

  split_t.Reverse();
  for (int i=0; i<split_t.Count(); i++){
    //if split parameter is near start or end, just adjust domain.
    double t0, t1;
    m_E[edge_index].GetDomain(&t0, &t1);
    if (t1 - t0 < 10.0*ON_ZERO_TOLERANCE) continue;

    //6 Dec 2002 Dale Lear:
    //   I added the relative edge_split_s and trm_split_s tests to detect
    //   attempts to trim a nano-gnats-wisker of the end of a trim.

    // set to true if edge should be trimmed instead of split.
    bool bTrimEdgeEnd = false; 

    double edge_split_s = ON_Interval(t0,t1).NormalizedParameterAt(split_t[i]);
    double trim_split_s = 0.5;

    if (split_t[i] - t0 <= ON_ZERO_TOLERANCE || edge_split_s <= ON_SQRT_EPSILON )
    {
      continue;
    }
    if (t1 - split_t[i] <= ON_ZERO_TOLERANCE || edge_split_s >= 1.0-ON_SQRT_EPSILON)
    {
      continue;
    }

    // trim_t[] = corresponding trim parameters
    ON_SimpleArray<double> trim_t(m_E[edge_index].m_ti.Count());


    if ( !bTrimEdgeEnd )
    {
      for (int j=0; j<m_E[edge_index].m_ti.Count(); j++)
      {
        double t;
        if (!GetTrimParameter(m_E[edge_index].m_ti[j], split_t[i], &t)){
          rc = false;
          continue;
        }
        trim_t.Append(t);
        const ON_BrepTrim& trim = m_T[m_E[edge_index].m_ti[j]];
        ON_Interval trim_domain = trim.Domain();
        trim_split_s = trim_domain.NormalizedParameterAt(t);
        if ( trim_split_s <= ON_SQRT_EPSILON || t - trim_domain[0] <= ON_ZERO_TOLERANCE )
        {
          bTrimEdgeEnd = true;
          break;
        }
        if ( trim_split_s >= 1.0-ON_SQRT_EPSILON || trim_domain[1] - t <= ON_ZERO_TOLERANCE )
        {
          bTrimEdgeEnd = true;
          break;
        }
      }
    }

    if ( bTrimEdgeEnd )
    {
      // if we get here, a split parameter we got was too close to
      // the end of the edge or a trim for us to split it.
      if ( edge_split_s <= 0.01 )
      {
        if ( t0 < split_t[i] )
          m_E[edge_index].ON_CurveProxy::Trim(ON_Interval(split_t[i], t1));
      }
      else if ( edge_split_s >= 0.99 )
      {
        if ( split_t[i] < t1 )
          m_E[edge_index].ON_CurveProxy::Trim(ON_Interval(t0, split_t[i]));
      }
      else
      {
        // no decent agreement between trims and edges - continue
        // with other parameters but this is basically the same
        // thing as having SplitEdge() fail.
        rc = false;
      }

      continue;
    }

    if (!SplitEdge(edge_index, split_t[i], trim_t)){
      rc = false;
      continue;
    }
    ON_Curve* new_curve0 = m_E[edge_index].DuplicateCurve();
    if (new_curve0){
      m_E[edge_index].m_c3i = AddEdgeCurve(new_curve0);
      m_E[edge_index].SetProxyCurve(new_curve0);
    }
    ON_Curve* new_curve1 = m_E.Last()->DuplicateCurve();
    if (new_curve1){
      m_E.Last()->m_c3i = AddEdgeCurve(new_curve1);
      m_E.Last()->SetProxyCurve(new_curve1);
    }
  }

  return rc;
}

static ON_Curve* TuneupEdgeOrTrimRealCurve( ON_CurveProxy& curve, bool bEdge )
{
  const ON_Curve* real_curve = curve.ProxyCurve();
  if ( 0 == real_curve )
    return 0;

  ON_Curve* new_real_curve = 0;
  const ON_Interval curve_domain = curve.Domain();
  if (    curve.ProxyCurveIsReversed()
       || real_curve->Domain() != curve_domain 
       || curve.ProxyCurveDomain() != curve_domain 
     )
  {
    // The proxy is a subset or reversed real curve.
    // Make a new real curve that is the exact geometry
    // the edge or trim requires.

    if ( bEdge )
    {
      new_real_curve = curve.DuplicateCurve();
    }
    else
    {
      // Trim curves end up being converted to piecewise
      // bezier nurbs curves when their pbox is set, 
      // so do the conversion here and save time and
      // memory thrashing.
      ON_NurbsCurve* nurbs_curve = curve.NurbsCurve();
      if ( 0 == nurbs_curve )
      {
        // This is a serious error - investigate it.
        // If you need help, ask Dale Lear.
        ON_ERROR("trim.NurbsCurve() returned NULL");
        new_real_curve = curve.DuplicateCurve();
      }
      else
      {
        nurbs_curve->ChangeDimension(2);
        nurbs_curve->MakePiecewiseBezier();
        new_real_curve = nurbs_curve;
      }
    }

    if (0 != new_real_curve)
    {
      new_real_curve->SetDomain(curve_domain);
    }
  }

  return new_real_curve;
}

static 
ON_Curve* TuneupSplitteeHelper( const ON_Curve* curve )
{
  if ( 0 == curve )
    return 0;
  ON_Curve* replacement_curve = 0;
  const ON_PolyCurve* polycurve = ON_PolyCurve::Cast(curve);
  if ( 0 != polycurve )
  {
    // This was being done, but the casts and RemoveNesting() calls
    // were copied-and-pasted in multiple different places.
    const_cast<ON_PolyCurve*>(polycurve)->RemoveNesting();

    // I added this to make spliting much more robust.
    const_cast<ON_PolyCurve*>(polycurve)->SynchronizeSegmentDomains();

    if ( 1 == polycurve->Count() )
    {
      // I added this to get rid of superfluous polycurve wrapping
      ON_Curve* segment0 = polycurve->SegmentCurve(0);
      if ( 0 != segment0 && segment0->Domain() == polycurve->Domain() )
      {
        replacement_curve = segment0->DuplicateCurve();
      }
    }
  }
  return replacement_curve;
}


int ON_Brep::SplitEdgeAtParameters(
  int edge_index,
  int edge_t_count,
  const double* edge_t
  )
{
  // Default kink_tol_radians MUST BE ON_PI/180.0.
  //
  // The default kink tol must be kept in sync with the default for 
  // TL_Brep::SplitKinkyFace() and ON_Brep::SplitKinkyFace().
  // See comments in TL_Brep::SplitKinkyFace() for more details.

  if (0 == edge_t_count)
    return 0;
  if (0 == edge_t)
    return 0;
  if (edge_index < 0 || edge_index >= m_E.Count()) 
    return 0;
  ON_BrepEdge& E = m_E[edge_index];
  if (E.m_c3i < 0) 
    return 0;
  ON_Curve* curve = m_C3[E.m_c3i];
  if (!curve) 
    return 0;

  ON_Interval Edomain;
  if ( !E.GetDomain(&Edomain.m_t[0],&Edomain.m_t[1]) )
    return 0;
  if ( !Edomain.IsIncreasing() )
    return 0;

  // get a list of unique and valid splitting parameters
  ON_SimpleArray<double> split_t(edge_t_count);
  {
    for (int i = 0; i < edge_t_count; i++)
    {
      double e = edge_t[i];
      if ( !ON_IsValid(e) )
      {
        ON_ERROR("Invalid edge_t[] value");
        continue;
      }
      if ( e <= Edomain.m_t[0] )
      {
        ON_ERROR("edge_t[] <= start of edge domain");
        continue;
      }
      if ( e >= Edomain.m_t[1] )
      {
        ON_ERROR("edge_t[] >= end of edge domain");
        continue;
      }
      split_t.Append(e);
    }
    if ( split_t.Count() > 1 )
    {
      // sort split_t[] and remove duplicates
      ON_SortDoubleArray( ON::heap_sort, split_t.Array(), split_t.Count() );
      int count = 1;
      for ( int i = 1; i < split_t.Count(); i++ )
      {
        if ( split_t[i] > split_t[count-1] )
        {
          if ( i > count )
            split_t[count] = split_t[i];
          count++;
        }
      }
      split_t.SetCount(count);
    }
  }

  if (split_t.Count() <= 0) 
    return 0;

  // Reverse split_t[] so the starting segment of the original
  // edge m_E[edge_index] is the one at m_E[edge_index].
  split_t.Reverse();

  ON_Curve* new_curve = TuneupSplitteeHelper(m_E[edge_index].ProxyCurve());
  if ( 0 != new_curve )
  {
    m_E[edge_index].m_c3i = AddEdgeCurve(new_curve);
    m_E[edge_index].SetProxyCurve(new_curve);
    new_curve = 0;
  }

  int eti, ti;
  int successful_split_count = 0;
  for (int i=0; i<split_t.Count(); i++)
  {
    double t0, t1;
    m_E[edge_index].GetDomain(&t0, &t1);
    if (t1 - t0 < 10.0*ON_ZERO_TOLERANCE) 
      break;

    //6 Dec 2002 Dale Lear:
    //   I added the relative edge_split_s and trm_split_s tests to detect
    //   attempts to trim a nano-gnats-wisker off the end of a trim.
    double edge_split_s = ON_Interval(t0,t1).NormalizedParameterAt(split_t[i]);
    double trim_split_s = 0.5;

    if (split_t[i] - t0 <= ON_ZERO_TOLERANCE || edge_split_s <= ON_SQRT_EPSILON )
    {
      // this split is not possible
      continue;
    }
    
    if (t1 - split_t[i] <= ON_ZERO_TOLERANCE || edge_split_s >= 1.0-ON_SQRT_EPSILON)
    {
      // this split is not possible
      continue;
    }

    // trim_t[] = corresponding trim parameters
    ON_SimpleArray<double> trim_t(m_E[edge_index].m_ti.Count());

    for ( eti = 0; eti < m_E[edge_index].m_ti.Count(); eti++)
    {
      ti = m_E[edge_index].m_ti[eti];
      if ( ti < 0 || ti >= m_T.Count() )
        continue;
      ON_BrepTrim& trim = m_T[ti];
      if ( 0 == i )
      {
        // On the first split, make sure the trim curve is up to snuff.
        new_curve = TuneupSplitteeHelper(trim.ProxyCurve());
        if (new_curve)
        {
          trim.m_c2i = AddTrimCurve(new_curve);
          trim.SetProxyCurve(new_curve);
          new_curve = 0;
        }
      }
      double t = ON_UNSET_VALUE;
      if (!GetTrimParameter(ti, split_t[i], &t) || !ON_IsValid(t))
        break;
      trim_t.Append(t);
      const ON_Interval trim_domain = trim.Domain();
      trim_split_s = trim_domain.NormalizedParameterAt(t);
      if ( trim_split_s <= ON_SQRT_EPSILON || t - trim_domain[0] <= ON_ZERO_TOLERANCE )
        break;

      if ( trim_split_s >= 1.0-ON_SQRT_EPSILON || trim_domain[1] - t <= ON_ZERO_TOLERANCE )
        break;
    }

    if ( trim_t.Count() != m_E[edge_index].m_ti.Count() )
      continue;

    if (!SplitEdge(edge_index, split_t[i], trim_t))
    {
      continue;
    }

    // SplitEdge generally adjusts proxy domains instead
    // of trimming the orginal curve. These DuplicateCurve()
    // calls make a new curve whose domain exactly matches
    // the edge.
    for ( int epart = 0; epart < 2; epart++ )
    {
      ON_BrepEdge* newE = (0 == epart) ? &m_E[edge_index] : m_E.Last();
      if ( 0 == newE )
        continue;
      new_curve = TuneupEdgeOrTrimRealCurve(*newE,true);
      if (new_curve)
      {
        newE->m_c3i = AddEdgeCurve(new_curve);
        newE->SetProxyCurve(new_curve);
      }
      for ( eti = 0; eti < newE->m_ti.Count(); eti++ )
      {
        ti = newE->m_ti[eti];
        if ( ti < 0 || ti >= m_T.Count() )
          continue;
        ON_BrepTrim& trim = m_T[ti];
        new_curve = TuneupEdgeOrTrimRealCurve(trim,false);
        if (new_curve)
        {
          trim.m_c2i = AddTrimCurve(new_curve);
          trim.SetProxyCurve(new_curve);
        }
      }
    }

    successful_split_count++;
  }

  return successful_split_count;
}

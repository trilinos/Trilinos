// Emacs will be in -*- Mode: c++ -*-
// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY: Tools for Automatic Differentiaton (order 1)
// RELEASE: 0.1
// USAGE  : You may copy freely these files and use it for
//          teaching or research. These or part of these may
//          not be sold or used for a commercial purpose with-
//          out our consent : fax (33)1 44 27 72 00
//
// AUTHOR : Nicolas Di cesare
// ORG    :
// E-MAIL : Nicolas.Dicesare@ann.jussieu.fr
//
// ORIG-DATE: September 97
// LAST-MOD : 28/07/98
// ********************************************************
// FILE   : tinyfadfunc.h
// ********************************************************
#ifndef _tinyfadfunc_h_
#define _tinyfadfunc_h_

template <int Num, class T> TinyFad<Num,T> exp (const TinyFad<Num,T>& in)
{
  TinyFad<Num,T> tmp(exp(in.val()));

  for (int i=0; i< Num; i++)
      tmp.dx(i) = in.dx(i)*exp(in.val());

  return tmp;
}

template <int Num, class T> TinyFad<Num,T> log (const TinyFad<Num,T>& in)
{
  if ( in.val() <= 0) error("TinyFad log (const TinyFad& in) : zero or negative value");
  TinyFad<Num,T> tmp(log(in.val()));

  for (int i=0; i< Num; i++)
      tmp.dx(i) = in.dx(i) / in.val();

  return tmp;
}

template <int Num, class T> TinyFad<Num,T> sqrt (const TinyFad<Num,T>& in)
{
  if ( in.val() < 0. ) error("TinyFad<Num,T> sqrt (const TinyFad& in) : negative value");
  TinyFad<Num,T> tmp(sqrt(in.val()));

  bool test=true;
  if ( in.val() == T(0.) ){
    for (int i=0; i< Num; i++)
      if ( in.dx(i) != T(0.) ) test = false;

    if ( !test )
      error("TinyFad<Num,T> sqrt (const TinyFad& in) : null value");
  }
  else
    for (int i=0; i< Num; i++)
      tmp.dx(i) = in.dx(i) / sqrt(in.val()) / 2.;

  return tmp;
}

template <int Num, class T> TinyFad<Num,T> sin (const TinyFad<Num,T>& in)
{
  TinyFad<Num,T> tmp( sin(in.val()) );

  for (int i=0; i< Num; i++)
    tmp.dx(i) = in.dx(i) * cos( in.val() );


  return tmp;
}

template <int Num, class T> TinyFad<Num,T> cos (const TinyFad<Num,T>& in)
{
  TinyFad<Num,T> tmp(cos(in.val()));

  for (int i=0; i< Num; i++)
    tmp.dx(i) = - in.dx(i) * sin( in.val() );

  return tmp;
}

template <int Num, class T> TinyFad<Num,T> tan (const TinyFad<Num,T>& in)
{
  if ( in.val() == 0) error("TinyFad tan (const TinyFad& in) undiefined in 0.");
  TinyFad<Num,T> tmp(tan(in.val()));

  for (int i=0; i< Num; i++){
     T cosinus = cos(in.val());
     tmp.dx(i) = in.dx(i) / cosinus / cosinus;
  }

  return tmp;
}

template <int Num, class T> TinyFad<Num,T> pow (const TinyFad<Num,T>& in, double e)
{
  TinyFad<Num,T> tmp(pow(in.val(), e));

  for (int i=0; i< Num; i++){
      tmp.dx(i) = e*in.dx(i)*pow(in.val(), e-1);
  }

  return tmp;
}

template <int Num, class T> TinyFad<Num,T> pow (const TinyFad<Num,T>& un, const TinyFad<Num,T>& deux)
{
  if (un.val() == 0) error("TinyFad pow (const TinyFad& un, const TinyFad& deux), un = 0. ");
  TinyFad<Num,T> tmp(pow(un.val(), deux.val()));

  for (int i=0; i< Num; i++){
      tmp.dx(i) = deux.dx(i) * log(un.val()) * pow(un.val(), deux.val())
	+ deux.val() * un.dx(i) * pow(un.val(), deux.val()-1);
  }

  return tmp;
}

template <int Num, class T> TinyFad<Num,T> pow (const TinyFad<Num,T>& in, const int e)
{
  TinyFad<Num,T> tmp( std::pow((double)in.val(), (double)e) );

  for (int i=0; i< Num; i++){
     tmp.dx(i) = e*in.dx(i)*std::pow((double)in.val(), (double)e-1);
  }

  return tmp;
}

template <int Num, class T> TinyFad<Num,T> abs (const TinyFad<Num,T>& in)
{
  int sign = in.val() > 0? 1:0;

  if (sign) return in;
  else return (-in);
}


#endif

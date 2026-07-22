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
// FILE   : tinyfadbin.h
// ********************************************************
#ifndef _tinyfadbin_h_
#define _tinyfadbin_h_


template <int Num, class L, class R> inline
TinyFad<Num,typename NumericalTraits<L,R>::promote > 
operator +(const TinyFad<Num,L>& un, const TinyFad<Num,R>& deux) {
  
  typedef typename NumericalTraits<L,R>::promote value_type;

  No_Initialization nothing;
  TinyFad<Num,value_type> tmp( nothing );

  for (int i=0; i<Num; ++i)
    tmp.dx(i) = un.dx(i) + deux.dx(i);
  
  tmp.val() = un.val() + deux.val();

  return tmp;
}

template <int Num, class L, class R> inline
TinyFad<Num,typename NumericalTraits<L,R>::promote> 
operator +(const TinyFad<Num,L>& un, const R& deux) {

  typedef typename NumericalTraits<L,R>::promote value_type;

  No_Initialization nothing;
  TinyFad<Num,value_type> tmp( nothing );

  for (int i=0; i<Num; ++i)
    tmp.dx(i) = un.dx(i);

  tmp.val() = un.val() + deux;

  return tmp;
}

template <int Num, class L, class R> inline
TinyFad<Num,typename NumericalTraits<L,R>::promote > 
operator +(const L& un, const TinyFad<Num,R>& deux) {
  return operator +(deux,un);
}

template <int Num, class L, class R> inline
TinyFad<Num,typename NumericalTraits<L,R>::promote > 
operator *(const TinyFad<Num,L>& un, const TinyFad<Num,R>& deux) {

  typedef typename NumericalTraits<L,R>::promote value_type;

  No_Initialization nothing;
  TinyFad<Num,value_type> tmp( nothing );

  for (int i=0; i<Num; ++i)
    tmp.dx(i) = un.dx(i)*deux.val() + un.val() * deux.dx(i);
  
  tmp.val() = un.val() * deux.val();

  return tmp;
}

template <int Num, class L, class R> inline
TinyFad<Num,typename NumericalTraits<L,R>::promote > 
operator *(const TinyFad<Num,L>& un, const R& deux) {

  typedef typename NumericalTraits<L,R>::promote value_type;

  No_Initialization nothing;
  TinyFad<Num,value_type> tmp( nothing );

  for (int i=0; i<Num; ++i)
    tmp.dx(i) = un.dx(i)*deux;
  
  tmp.val() = un.val() * deux;


  return tmp;
}

template <int Num, class L, class R> inline
TinyFad<Num,typename NumericalTraits<L,R>::promote > 
operator *(const L& un, const TinyFad<Num,R>& deux) {

  return operator *(deux,un);
}

//-------------------------------------------------------
template <int Num, class L, class R> inline
TinyFad<Num,typename NumericalTraits<L,R>::promote > 
operator -(const TinyFad<Num,L> & un, const TinyFad<Num,R> & deux) {

  typedef typename NumericalTraits<L,R>::promote value_type;

  No_Initialization nothing;
  TinyFad<Num,value_type> tmp( nothing );

  for (int i=0; i<Num; ++i)
    tmp.dx(i) = un.dx(i) - deux.dx(i);
  
  tmp.val() = un.val() - deux.val();

  return tmp;
}

template <int Num, class L, class R> inline
TinyFad<Num,typename NumericalTraits<L,R>::promote> 
operator -(const L & un, const TinyFad<Num,R> & deux) {

  typedef typename NumericalTraits<L,R>::promote value_type;

  No_Initialization nothing;
  TinyFad<Num,value_type> tmp( nothing );

  for (int i=0; i<Num; ++i)
    tmp.dx(i) -= deux.dx(i);
  
  tmp.val() = un - deux.val();

  return tmp;
}

template <int Num, class L, class R> inline
TinyFad<Num, typename NumericalTraits<L,R>::promote > 
operator -(const TinyFad<Num,L> & un, const R & deux) {

  typedef typename NumericalTraits<L,R>::promote value_type;

  No_Initialization nothing;
  TinyFad<Num,value_type> tmp( nothing );

  for (int i=0; i<Num; ++i)
    tmp.dx(i) = un.dx(i);
  
  tmp.val() = un.val() - deux;

  return tmp;
}

template <int Num, class L, class R> inline
TinyFad<Num, typename NumericalTraits<L,R>::promote > 
operator /(const TinyFad<Num,L> & un, const TinyFad<Num,R> & deux) {

  typedef typename NumericalTraits<L,R>::promote value_type;

  if (deux.val() == 0.) error("TinyFad & TinyFad::operator /(const TinyFad<Num,L> & un, const TinyFad<Num,R> & deux), dividing by 0");

  No_Initialization nothing;
  TinyFad<Num,value_type> tmp( nothing );
  value_type dval = deux.val();

  for (int i=0; i<Num; ++i)
    tmp.dx(i) = ( un.dx(i)* deux.val() - un.val() * deux.dx(i) ) / dval / dval ;
  
  tmp.val() = un.val() / dval;

  return tmp;
}

template <int Num, class L, class R> inline
TinyFad<Num, typename NumericalTraits<L,R>::promote > 
operator /(const L & un, const TinyFad<Num,R> & deux) {

  typedef typename NumericalTraits<L,R>::promote value_type;

  if (deux.val() == 0.) error("TinyFad & TinyFad::operator /(const L & un, const TinyFad<Num,R> & deux), dividing by 0");

  No_Initialization nothing;
  TinyFad<Num,value_type> tmp( nothing );
  value_type dval = deux.val();

  for (int i=0; i<Num; ++i)
    tmp.dx(i) = - un * deux.dx(i)  / dval / dval ;
  
  tmp.val() = un / dval;

  return tmp;
}

template <int Num, class  L, class R> inline
TinyFad<Num, typename NumericalTraits<L,R>::promote > 
operator /(const TinyFad<Num,L> & un, const R & deux) {

  typedef typename NumericalTraits<L,R>::promote value_type;

  if (deux == 0.) error("TinyFad & TinyFad::operator /(const TinyFad<Num,L> & un, const R & deux), dividing by 0");

  No_Initialization nothing;
  TinyFad<Num,value_type> tmp( nothing );

  for (int i=0; i<Num; ++i)
    tmp.dx(i) = un.dx(i)  / deux;
  
  tmp.val() = un.val() / deux;

   return tmp;
}

#endif

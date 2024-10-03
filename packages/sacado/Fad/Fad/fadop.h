// Emacs will be in -*- Mode: c++ -*-
// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// ************ DO NOT REMOVE THIS BANNER ****************
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//            CEMRACS 98 : C++ courses,
//         templates : new C++ techniques
//            for scientific computing
//
//********************************************************
//
//  A short implementation ( not all operators and
//  functions are overloaded ) of 1st order Automatic
//  Differentiation in forward mode (FAD) using
//  EXPRESSION TEMPLATES.
//
//********************************************************
#ifndef _fadop_h_
#define _fadop_h_

#include <algorithm>

//------------------------------- Fad binary operators ------------------------------------------


//------------------------------- Fad addition operators ------------------------------------------
template <class L, class R> class FadBinaryAdd {
public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;

  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

protected:

  const L& left_; const R& right_;

public:
  FadBinaryAdd(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryAdd() {;}


  const value_type val() const {return left_.val() + right_.val();}
  const value_type dx(int i) const {return left_.dx(i) + right_.dx(i);}
  int size() const {
    int lsz = left_.size(), rsz = right_.size();
    return std::max(lsz, rsz);
  }

  bool hasFastAccess() const { return left_.hasFastAccess() && right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i)+right_.fastAccessDx(i);}
};


template <class L> class FadBinaryAdd<L, FadCst<typename L::value_type> > {
public:
  typedef typename L::value_type value_type;
  typedef FadCst<value_type> R;

protected:

  const L& left_; const  R right_;

public:
  FadBinaryAdd(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryAdd() {;}


  const value_type val() const {return left_.val() + right_.val();}
  const value_type dx(int i) const {return left_.dx(i);}
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i);}
};


template <class R> class FadBinaryAdd< FadCst<typename R::value_type>, R> {
public:
  typedef typename R::value_type value_type;
  typedef FadCst<value_type> L;

protected:

  const L left_; const R& right_;

public:
  FadBinaryAdd(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryAdd() {;}


  const value_type val() const {return left_.val() + right_.val();}
  value_type dx(int i) const {return right_.dx(i);}
  int size() const {return right_.size();}

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return right_.fastAccessDx(i);}
};






//------------------------------- Fad substraction operators ------------------------------------------
template <class L, class R> class FadBinaryMinus {
public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

protected:

  const L& left_; const R& right_;

public:
  FadBinaryMinus(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryMinus() {;}


  const value_type val() const {return left_.val() - right_.val();}
  const value_type dx(int i) const {return left_.dx(i) - right_.dx(i);}
  int size() const {
    int lsz = left_.size(), rsz = right_.size();
    return std::max(lsz, rsz);
  }

  bool hasFastAccess() const { return left_.hasFastAccess() && right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i) - right_.fastAccessDx(i);}
};


template <class L> class FadBinaryMinus<L, FadCst<typename L::value_type> > {
public:
  typedef typename L::value_type value_type;
  typedef FadCst<value_type> R;

protected:

  const L& left_; const R right_;

public:
  FadBinaryMinus(const L& left, const R & rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryMinus() {;}


  const value_type val() const {return left_.val() - right_.val();}
  const value_type dx(int i) const {return left_.dx(i);}
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i);}
};


template <class R> class FadBinaryMinus< FadCst<typename R::value_type>, R> {
public:
  typedef typename R::value_type value_type;
  typedef FadCst<value_type> L;

protected:

  const L left_; const R& right_;

public:
  FadBinaryMinus(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryMinus() {;}


  const value_type val() const {return left_.val() - right_.val();}
  const value_type dx(int i) const {return - right_.dx(i);}
  int size() const { return right_.size(); }

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return - right_.fastAccessDx(i);}
};


//------------------------------- Fad multiplication operators ------------------------------------------
template <class L, class R> class FadBinaryMul {
 public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

 protected:

  const L& left_; const R& right_;

 public:
  FadBinaryMul(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryMul() {;}

  const value_type val() const {return left_.val() * right_.val() ;}
  const value_type dx(int i) const {return  left_.dx(i) * right_.val() + right_.dx(i) * left_.val();}
  int size() const {
    int lsz = left_.size(), rsz = right_.size();
    return lsz > rsz ? lsz : rsz;
  }

  bool hasFastAccess() const { return left_.hasFastAccess() && right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i) * right_.val() + right_.fastAccessDx(i) * left_.val();}

};

template <class L> class FadBinaryMul<L, FadCst<typename L::value_type> > {
 public:
  typedef typename L::value_type value_type;
  typedef FadCst<value_type> R;

 protected:

  const L& left_; const R right_;

 public:
  FadBinaryMul(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryMul() {;}

  const value_type val() const {return left_.val() * right_.val() ;}
  const value_type dx(int i) const {return  left_.dx(i) * right_.val();}
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i) * right_.val();}
};

template <class R> class FadBinaryMul< FadCst<typename R::value_type>, R> {
 public:
  typedef typename R::value_type value_type;
  typedef FadCst<value_type> L;

 protected:

  const L left_; const R& right_;

 public:
  FadBinaryMul(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryMul() {;}

  const value_type val() const {return left_.val() * right_.val() ;}
  const value_type dx(int i) const {return  right_.dx(i) * left_.val();}
  int size() const { return right_.size();}

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return right_.fastAccessDx(i) * left_.val();}
};


//------------------------------- Fad division operators ------------------------------------------
template <class L, class R> class FadBinaryDiv {
 public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

 protected:

  const L& left_; const R& right_;

 public:
  FadBinaryDiv(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryDiv() {;}


  const value_type val() const {return left_.val() / right_.val();}
  const value_type dx(int i) const {return  (left_.dx(i) * right_.val() - right_.dx(i) * left_.val() ) / (right_.val() * right_.val()) ;}
  int size() const {
    int lsz = left_.size(), rsz = right_.size();
    return std::max(lsz, rsz);
  }

  bool hasFastAccess() const { return left_.hasFastAccess() && right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return (left_.fastAccessDx(i) * right_.val() - right_.fastAccessDx(i) * left_.val() )
             / (right_.val() * right_.val()) ;}
};


template <class L> class FadBinaryDiv<L, FadCst<typename L::value_type> > {
 public:
  typedef typename L::value_type value_type;
  typedef FadCst<value_type> R;

 protected:

  const L& left_; const R right_;

 public:
  FadBinaryDiv(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryDiv() {;}


  const value_type val() const {return left_.val() / right_.val();}
  const value_type dx(int i) const {return  left_.dx(i) / right_.val();}
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i) / right_.val() ;}
};


template <class R> class FadBinaryDiv< FadCst<typename R::value_type>, R> {
 public:
  typedef typename R::value_type value_type;
  typedef FadCst<value_type> L;

 protected:

  const L left_; const R& right_;

 public:
  FadBinaryDiv(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryDiv() {;}

  const value_type val() const {return left_.val() / right_.val();}
  const value_type dx(int i) const {return  (- right_.dx(i) * left_.val() ) / (right_.val() * right_.val()) ;}
  int size() const { return right_.size();}

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return (- right_.fastAccessDx(i) * left_.val() )
             / (right_.val() * right_.val()) ;}
};


//------------------------------- Fad pow function ------------------------------------------
template <class L, class R> class FadBinaryPow {
public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;

  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

protected:

  const L& left_; const R& right_;

public:
  FadBinaryPow(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryPow() {;}


  const value_type val() const {return std::pow( left_.val(), right_.val() );}
  const value_type dx(int i) const
    {
      return  (right_.dx(i)*std::log(left_.val())+right_.val()*left_.dx(i)/left_.val())
  *std::pow( left_.val(), right_.val() );
    }
  int size() const {
    int lsz = left_.size(), rsz = right_.size();
    return std::max(lsz, rsz);
  }

  bool hasFastAccess() const { return left_.hasFastAccess() && right_.hasFastAccess();}
  value_type fastAccessDx(int i) const
    {
      return  (right_.fastAccessDx(i)*std::log(left_.val())+right_.val()*left_.fastAccessDx(i)/left_.val())
  *std::pow( left_.val(), right_.val() );
    }
};


template <class L> class FadBinaryPow<L, FadCst<typename L::value_type> > {
public:
  typedef typename L::value_type value_type;
  typedef FadCst<value_type> R;

protected:

  const L& left_; const  R right_;

public:
  FadBinaryPow(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryPow() {;}


  const value_type val() const {return std::pow(left_.val(),right_.val()) ;}
  const value_type dx(int i) const
    {
      return  (right_.val()*left_.dx(i)/left_.val())*std::pow( left_.val(), right_.val() );
    }
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const
    {
      return  (right_.val()*left_.fastAccessDx(i)/left_.val())
  *std::pow( left_.val(), right_.val() );
    }
};


template <class R> class FadBinaryPow< FadCst<typename R::value_type>, R> {
public:
  typedef typename R::value_type value_type;
  typedef FadCst<value_type> L;

protected:

  const L left_; const R& right_;

public:
  FadBinaryPow(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryPow() {;}

  const value_type val() const {return std::pow(left_.val(),right_.val());}
  value_type dx(int i) const
    {
      return (right_.dx(i)*std::log(left_.val()))*std::pow( left_.val(), right_.val() );
    }
  int size() const {return right_.size();}

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const
    {
      return  (right_.fastAccessDx(i)*std::log(left_.val()))
  *std::pow( left_.val(), right_.val() );
    }
};

template <class L> class FadBinaryPow< L , int> {
public:
  typedef typename L::value_type value_type;
  typedef FadCst<int> R;

protected:

  const L& left_; const R right_;

public:
  FadBinaryPow(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryPow() {;}


  const value_type val() const {return std::pow(left_.val(),right_.val());}
  value_type dx(int i) const
    {
      return right_.val()*std::pow( left_.val(), right_.val()-1);
    }
  int size() const {return right_.size();}

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const
    {
      return  right_.val() * std::pow( left_.val(), right_.val()-1 );
    }
};


//------------------------------- Fad operators ------------------------------------------
#define FAD_BIN_MACRO(OP,TYPE)                                    \
template<class E1, class E2>                                      \
inline FadExpr< TYPE< FadExpr<E1>, FadExpr<E2> > >                \
OP  (const FadExpr<E1> &v, const FadExpr<E2> &w){                 \
    typedef TYPE<FadExpr<E1>, FadExpr<E2> > expr_t;               \
    return FadExpr<expr_t> (expr_t (v , w ));                     \
}                                                                 \
                                                                  \
template<class E>                                                 \
inline FadExpr<TYPE<FadExpr<E>,FadCst<typename E::value_type> > > \
OP (const FadExpr<E> &e, const typename E::value_type &t){        \
    typedef typename E::value_type A;                             \
    typedef TYPE<FadExpr<E>,FadCst<A> > expr_t;                   \
    return FadExpr<expr_t>(expr_t (e, FadCst<A> (t)));            \
}                                                                 \
                                                                  \
template<typename A>                                              \
inline FadExpr<TYPE<FadCst<A>,Fad<A> > >                          \
OP (const A& a, const Fad<A> &e){                                 \
    typedef TYPE<FadCst<A>,Fad<A> > expr_t;                       \
    return FadExpr<expr_t> (expr_t (FadCst<A>(a), e  ));          \
}                                                                 \
                                                                  \
template<typename A>                                              \
inline FadExpr<TYPE<Fad<A>,FadCst<A> > >                          \
OP (const Fad<A> &e, const A& a){                                 \
    typedef TYPE<Fad<A>,FadCst<A> > expr_t;                       \
    return FadExpr<expr_t>(expr_t (e ,FadCst<A>(a)));             \
}                                                                 \
                                                                  \
template<class E>                                                 \
inline FadExpr<TYPE<FadCst<typename E::value_type>,FadExpr<E> > > \
OP (const typename E::value_type &t, const FadExpr<E> &e){        \
    typedef typename E::value_type A;                             \
    typedef TYPE<FadCst<A>,FadExpr<E> > expr_t;                   \
    return FadExpr<expr_t> (expr_t (FadCst<A> (t),e ));           \
}                                                                 \
                                                                  \
template<class E>                                                 \
inline FadExpr<TYPE<FadExpr<E>,Fad<typename E::value_type> > >    \
OP (const FadExpr<E> &e,const Fad<typename E::value_type>& v){    \
    typedef TYPE<FadExpr<E>,Fad<typename E::value_type> > expr_t; \
    return FadExpr<expr_t>(expr_t (e, v ));                       \
}                                                                 \
                                                                  \
template<typename A>                                              \
inline FadExpr<TYPE<Fad<A>,Fad<A> > >                             \
OP (const Fad<A> &e1,const Fad<A>& e2){                           \
    typedef TYPE<Fad<A>,Fad<A> > expr_t;                          \
    return FadExpr<expr_t>(expr_t (e1 , e2 ));                    \
}                                                                 \
                                                                  \
template<class E>                                                 \
inline FadExpr<TYPE<Fad<typename E::value_type>,FadExpr<E> > >    \
OP (const Fad<typename E::value_type> &v, const FadExpr<E> &e){   \
    typedef TYPE<Fad<typename E::value_type>,FadExpr<E> > expr_t; \
    return FadExpr<expr_t> (expr_t (v , e ));                     \
}

FAD_BIN_MACRO(operator+,FadBinaryAdd)
FAD_BIN_MACRO(operator-,FadBinaryMinus)
FAD_BIN_MACRO(operator*,FadBinaryMul)
FAD_BIN_MACRO(operator/,FadBinaryDiv)

FAD_BIN_MACRO(pow,FadBinaryPow)

#undef FAD_BIN_MACRO


#endif

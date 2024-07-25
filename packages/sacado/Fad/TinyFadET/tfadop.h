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
#ifndef _tfadop_h_
#define _tfadop_h_

//------------------------------- TFad binary operators ------------------------------------------


//------------------------------- TFad addition operators ------------------------------------------
template <class L, class R> class TFadBinaryAdd {
public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;

  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

protected:

  const L& left_; const R& right_;

public:
  TFadBinaryAdd(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~TFadBinaryAdd() {;}


  const value_type val() const {return left_.val() + right_.val();}
  const value_type dx(int i) const {return left_.dx(i) + right_.dx(i);}
  int size() const {
    int lsz = left_.size(), rsz = right_.size();
    return std::max(lsz, rsz);
  }

  bool hasFastAccess() const { return left_.hasFastAccess() && right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i)+right_.fastAccessDx(i);}
};


template <class L> class TFadBinaryAdd<L, TFadCst<typename L::value_type> > {
public:
  typedef typename L::value_type value_type;
  typedef TFadCst<value_type> R;

protected:

  const L& left_; const  R right_;

public:
  TFadBinaryAdd(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~TFadBinaryAdd() {;}


  const value_type val() const {return left_.val() + right_.val();}
  const value_type dx(int i) const {return left_.dx(i);}
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i);}
};


template <class R> class TFadBinaryAdd< TFadCst<typename R::value_type>, R> {
public:
  typedef typename R::value_type value_type;
  typedef TFadCst<value_type> L;

protected:

  const L left_; const R& right_;

public:
  TFadBinaryAdd(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~TFadBinaryAdd() {;}


  const value_type val() const {return left_.val() + right_.val();}
  value_type dx(int i) const {return right_.dx(i);}
  int size() const {return right_.size();}

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return right_.fastAccessDx(i);}
};






//------------------------------- TFad substraction operators ------------------------------------------
template <class L, class R> class TFadBinaryMinus {
public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

protected:

  const L& left_; const R& right_;

public:
  TFadBinaryMinus(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~TFadBinaryMinus() {;}


  const value_type val() const {return left_.val() - right_.val();}
  const value_type dx(int i) const {return left_.dx(i) - right_.dx(i);}
  int size() const {
    int lsz = left_.size(), rsz = right_.size();
    return std::max(lsz, rsz);
  }

  bool hasFastAccess() const { return left_.hasFastAccess() && right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i) - right_.fastAccessDx(i);}
};


template <class L> class TFadBinaryMinus<L, TFadCst<typename L::value_type> > {
public:
  typedef typename L::value_type value_type;
  typedef TFadCst<value_type> R;

protected:

  const L& left_; const R right_;

public:
  TFadBinaryMinus(const L& left, const R & rigth) : left_(left), right_(rigth) {;}
  ~TFadBinaryMinus() {;}


  const value_type val() const {return left_.val() - right_.val();}
  const value_type dx(int i) const {return left_.dx(i);}
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i);}
};


template <class R> class TFadBinaryMinus< TFadCst<typename R::value_type>, R> {
public:
  typedef typename R::value_type value_type;
  typedef TFadCst<value_type> L;

protected:

  const L left_; const R& right_;

public:
  TFadBinaryMinus(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~TFadBinaryMinus() {;}


  const value_type val() const {return left_.val() - right_.val();}
  const value_type dx(int i) const {return - right_.dx(i);}
  int size() const { return right_.size(); }

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return - right_.fastAccessDx(i);}
};


//------------------------------- TFad multiplication operators ------------------------------------------
template <class L, class R> class TFadBinaryMul {
 public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

 protected:

  const L& left_; const R& right_;

 public:
  TFadBinaryMul(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~TFadBinaryMul() {;}

  const value_type val() const {return left_.val() * right_.val() ;}
  const value_type dx(int i) const {return  left_.dx(i) * right_.val() + right_.dx(i) * left_.val();}
  int size() const {
    int lsz = left_.size(), rsz = right_.size();
    return std::max(lsz, rsz);
  }

  bool hasFastAccess() const { return left_.hasFastAccess() && right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i) * right_.val() + right_.fastAccessDx(i) * left_.val();}

};

template <class L> class TFadBinaryMul<L, TFadCst<typename L::value_type> > {
 public:
  typedef typename L::value_type value_type;
  typedef TFadCst<value_type> R;

 protected:

  const L& left_; const R right_;

 public:
  TFadBinaryMul(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~TFadBinaryMul() {;}

  const value_type val() const {return left_.val() * right_.val() ;}
  const value_type dx(int i) const {return  left_.dx(i) * right_.val();}
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i) * right_.val();}
};

template <class R> class TFadBinaryMul< TFadCst<typename R::value_type>, R> {
 public:
  typedef typename R::value_type value_type;
  typedef TFadCst<value_type> L;

 protected:

  const L left_; const R& right_;

 public:
  TFadBinaryMul(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~TFadBinaryMul() {;}

  const value_type val() const {return left_.val() * right_.val() ;}
  const value_type dx(int i) const {return  right_.dx(i) * left_.val();}
  int size() const { return right_.size();}

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return right_.fastAccessDx(i) * left_.val();}
};


//------------------------------- TFad division operators ------------------------------------------
template <class L, class R> class TFadBinaryDiv {
 public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

 protected:

  const L& left_; const R& right_;

 public:
  TFadBinaryDiv(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~TFadBinaryDiv() {;}


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


template <class L> class TFadBinaryDiv<L, TFadCst<typename L::value_type> > {
 public:
  typedef typename L::value_type value_type;
  typedef TFadCst<value_type> R;

 protected:

  const L& left_; const R right_;

 public:
  TFadBinaryDiv(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~TFadBinaryDiv() {;}


  const value_type val() const {return left_.val() / right_.val();}
  const value_type dx(int i) const {return  left_.dx(i) / right_.val();}
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i) / right_.val() ;}
};


template <class R> class TFadBinaryDiv< TFadCst<typename R::value_type>, R> {
 public:
  typedef typename R::value_type value_type;
  typedef TFadCst<value_type> L;

 protected:

  const L left_; const R& right_;

 public:
  TFadBinaryDiv(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~TFadBinaryDiv() {;}

  const value_type val() const {return left_.val() / right_.val();}
  const value_type dx(int i) const {return  (- right_.dx(i) * left_.val() ) / (right_.val() * right_.val()) ;}
  int size() const { return right_.size();}

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return (- right_.fastAccessDx(i) * left_.val() )
             / (right_.val() * right_.val()) ;}
};


//------------------------------- TFad pow function ------------------------------------------
template <class L, class R> class TFadBinaryPow {
public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;

  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

protected:

  const L& left_; const R& right_;

public:
  TFadBinaryPow(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~TFadBinaryPow() {;}


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


template <class L> class TFadBinaryPow<L, TFadCst<typename L::value_type> > {
public:
  typedef typename L::value_type value_type;
  typedef TFadCst<value_type> R;

protected:

  const L& left_; const  R right_;

public:
  TFadBinaryPow(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~TFadBinaryPow() {;}


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


template <class R> class TFadBinaryPow< TFadCst<typename R::value_type>, R> {
public:
  typedef typename R::value_type value_type;
  typedef TFadCst<value_type> L;

protected:

  const L left_; const R& right_;

public:
  TFadBinaryPow(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~TFadBinaryPow() {;}

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

template <class L> class TFadBinaryPow< L , int> {
public:
  typedef typename L::value_type value_type;
  typedef TFadCst<int> R;

protected:

  const L& left_; const R right_;

public:
  TFadBinaryPow(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~TFadBinaryPow() {;}


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


//------------------------------- TFad operators ------------------------------------------
#define FAD_BIN_MACRO(OP,TYPE)                                    \
template<class E1, class E2>                                      \
inline TFadExpr< TYPE< TFadExpr<E1>, TFadExpr<E2> > >                \
OP  (const TFadExpr<E1> &v, const TFadExpr<E2> &w){                 \
    typedef TYPE<TFadExpr<E1>, TFadExpr<E2> > expr_t;               \
    return TFadExpr<expr_t> (expr_t (v , w ));                     \
}                                                                 \
                                                                  \
template<class E>                                                 \
inline TFadExpr<TYPE<TFadExpr<E>,TFadCst<typename E::value_type> > > \
OP (const TFadExpr<E> &e, const typename E::value_type &t){        \
    typedef typename E::value_type A;                             \
    typedef TYPE<TFadExpr<E>,TFadCst<A> > expr_t;                   \
    return TFadExpr<expr_t>(expr_t (e, TFadCst<A> (t)));            \
}                                                                 \
                                                                  \
template<typename A,int Num>                                              \
inline TFadExpr<TYPE<TFadCst<A>,TFad<Num,A> > >                          \
OP (const A& a, const TFad<Num,A> &e){                                 \
    typedef TYPE<TFadCst<A>,TFad<Num,A> > expr_t;                       \
    return TFadExpr<expr_t> (expr_t (TFadCst<A>(a), e  ));          \
}                                                                 \
                                                                  \
template<typename A,int Num>                                              \
inline TFadExpr<TYPE<TFad<Num,A>,TFadCst<A> > >                          \
OP (const TFad<Num,A> &e, const A& a){                                 \
    typedef TYPE<TFad<Num,A>,TFadCst<A> > expr_t;                       \
    return TFadExpr<expr_t>(expr_t (e ,TFadCst<A>(a)));             \
}                                                                 \
                                                                  \
template<class E>                                                 \
inline TFadExpr<TYPE<TFadCst<typename E::value_type>,TFadExpr<E> > > \
OP (const typename E::value_type &t, const TFadExpr<E> &e){        \
    typedef typename E::value_type A;                             \
    typedef TYPE<TFadCst<A>,TFadExpr<E> > expr_t;                   \
    return TFadExpr<expr_t> (expr_t (TFadCst<A> (t),e ));           \
}                                                                 \
                                                                  \
template<class E,int Num>                                                 \
inline TFadExpr<TYPE<TFadExpr<E>,TFad<Num,typename E::value_type> > >    \
OP (const TFadExpr<E> &e,const TFad<Num,typename E::value_type>& v){    \
    typedef TYPE<TFadExpr<E>,TFad<Num,typename E::value_type> > expr_t; \
    return TFadExpr<expr_t>(expr_t (e, v ));                       \
}                                                                 \
                                                                  \
template<typename A,int Num>                                              \
inline TFadExpr<TYPE<TFad<Num,A>,TFad<Num,A> > >                             \
OP (const TFad<Num,A> &e1,const TFad<Num,A>& e2){                           \
    typedef TYPE<TFad<Num,A>,TFad<Num,A> > expr_t;                          \
    return TFadExpr<expr_t>(expr_t (e1 , e2 ));                    \
}                                                                 \
                                                                  \
template<class E, int Num>                                                 \
inline TFadExpr<TYPE<TFad<Num,typename E::value_type>,TFadExpr<E> > >    \
OP (const TFad<Num,typename E::value_type> &v, const TFadExpr<E> &e){   \
    typedef TYPE<TFad<Num,typename E::value_type>,TFadExpr<E> > expr_t; \
    return TFadExpr<expr_t> (expr_t (v , e ));                     \
}

FAD_BIN_MACRO(operator+,TFadBinaryAdd)
FAD_BIN_MACRO(operator-,TFadBinaryMinus)
FAD_BIN_MACRO(operator*,TFadBinaryMul)
FAD_BIN_MACRO(operator/,TFadBinaryDiv)

FAD_BIN_MACRO(pow,TFadBinaryPow)

#undef FAD_BIN_MACRO


#endif

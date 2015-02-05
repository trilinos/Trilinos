
#include "Sacado.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"

using namespace ROL;

// This is a tweak of example_01c which uses a couple of C++11 features
// To build this example, you must add the cmake flag -D ENABLE_CPP11:BOOL=ON 


template<class Real> 
Teuchos::RCP<const std::vector<Real>> 
access_elements(const Vector<Real> &x) {
    return (Teuchos::dyn_cast<StdVector<Real>>
            (const_cast<Vector<Real>&>(x))).getVector();
}
  

template<class Real>
Teuchos::RCP<std::vector<Real>>
access_elements(Vector<Real> &x) {
    return Teuchos::rcp_const_cast<std::vector<Real>> 
           ((Teuchos::dyn_cast<StdVector<Real>> (x)).getVector());   
}



template<class Real>
class Zakharov {

    public:
    Real value(const Vector<Real> &x, Real &tol) {

        auto xp = access_elements(x);

        int n = xp->size();

        Real xdotx = 0;

        Real kdotx = 0;
        Real J = 0;
   
        // Compute dot products 
        for(int i=0; i<n; ++i) {
            xdotx += pow((*xp)[i],2);       // (k,x)
            kdotx += double(i+1)*(*xp)[i];  // (x,x)
        }

        // Sum terms in objective function
        J = xdotx + pow(kdotx,2)/4.0 + pow(kdotx,4)/16.0;
    
       return J;
   }
};



template<class Real>
class ObjectiveEnabler : public Objective<Real> {
    Zakharov<Real> zakharov;
    public:
        ObjectiveEnabler(const Zakharov<Real>& zak) : zakharov(zak) {}
        Real value(const Vector<Real> &x, Real &tol) { return zakharov.value(x,tol); }
};


template<class Real>
class GradientEnabler : public Objective<Real> {

    typedef Sacado::Fad::DFad<Real> FadType; 
    typedef StdVector<FadType> FadVector;

    Teuchos::RCP<Objective<FadType>> obj_;

    public:
    GradientEnabler(Teuchos::RCP<Objective<FadType>> obj) : obj_(obj) {} 

    Real value(const Vector<Real> &x, Real &tol) {

        auto xp = access_elements(x);
 
        int n = xp->size();

        auto x_fad_rcp = Teuchos::rcp(new std::vector<FadType>);

        x_fad_rcp->reserve(n);

        for(int i=0; i<n; ++i) {
            x_fad_rcp->push_back((*xp)[i]);
        }

        FadVector x_fad(x_fad_rcp);

        FadType tol_fad = tol;  
        FadType J_fad = obj_->value(x_fad,tol_fad);
            
        return J_fad.val();
    }

    void gradient(Vector<Real> &g, const Vector<Real> &x, Real &tol) {

        auto xp = access_elements(x);
        auto gp = access_elements(g);
        int  n  = xp->size();
 
        auto x_fad_rcp = Teuchos::rcp( new std::vector<FadType> );

        x_fad_rcp->reserve(n);
   
        for(int i=0; i<n; ++i) {
            x_fad_rcp->push_back(FadType(n,i,(*xp)[i])); 
        }

        FadVector x_fad(x_fad_rcp);

        FadType tol_fad = tol;  
        FadType J_fad = obj_->value(x_fad,tol_fad);

        for(int i=0; i<n; ++i) {
            (*gp)[i] = J_fad.dx(i);
        }
    } 
};




template<class Real>
class HessianEnabler : public Objective<Real> {

    typedef Sacado::Fad::SFad<Real,1> FadType; 
    typedef StdVector<FadType> FadVector;

    Teuchos::RCP<Objective<FadType>> obj_;

    public: 
        HessianEnabler(Teuchos::RCP<Objective<FadType>> obj) : obj_(obj) {}

        Real value(const Vector<Real> &x, Real &tol) {

            auto xp = access_elements(x);
 
            int n = xp->size();

            auto x_fad_rcp = Teuchos::rcp(new std::vector<FadType>);

            x_fad_rcp->reserve(n);
            for(int i=0; i<n; ++i) {
                 x_fad_rcp->push_back((*xp)[i]);
            }

            FadVector x_fad(x_fad_rcp);

            FadType tol_fad = tol;  
            FadType J_fad = obj_->value(x_fad,tol_fad);
            
            return J_fad.val();
        }

        void gradient(Vector<Real> &g, const Vector<Real> &x, Real &tol) {
            
            auto xp = access_elements(x);
            auto gp = access_elements(g);

            int n = xp->size();

            auto x_fad_rcp = Teuchos::rcp(new std::vector<FadType>);
            auto g_fad_rcp = Teuchos::rcp(new std::vector<FadType>);

            x_fad_rcp->reserve(n);
            g_fad_rcp->reserve(n);

            for(int i=0; i<n; ++i) {
                 x_fad_rcp->push_back((*xp)[i]);
                 g_fad_rcp->push_back((*xp)[i]);
            }

            FadVector x_fad(x_fad_rcp);
            FadVector g_fad(g_fad_rcp);
            FadType tol_fad = tol;

            obj_->gradient(g_fad,x_fad,tol_fad);

            for(int i=0; i<n; ++i) {
                (*gp)[i] = (*g_fad_rcp)[i].val();    
            }
        }
 
        void hessVec(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
            auto xp  = access_elements(x);
            auto vp  = access_elements(v);
            auto hvp = access_elements(hv);
            int n = xp->size(); 

            // Create a vector of independent variables
            auto x_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
            x_fad_rcp->reserve(n);
 
            // Allocate for gradient   
            auto g_fad_rcp = Teuchos::rcp( new std::vector<FadType> );

            g_fad_rcp->reserve(n);

            for(int i=0; i<n; ++i) {
                x_fad_rcp->push_back(FadType(1,(*xp)[i]));
                g_fad_rcp->push_back(0);
                (*x_fad_rcp)[i].fastAccessDx(0) = (*vp)[i];
            }

            FadVector x_fad(x_fad_rcp);
            FadVector g_fad(g_fad_rcp);

            FadType tol_fad = tol; 
            obj_->gradient(g_fad,x_fad,tol_fad);

            for(int i=0; i<n; ++i) {
                (*hvp)[i] = (*g_fad_rcp)[i].dx(0); 
            } 
        }
};


template<class Real>
class ObjectiveAD : public Objective<Real> {

    typedef Sacado::Fad::SFad<Real,1> FadType;
    typedef Sacado::Fad::DFad<FadType> FadFadType; 

    private:
        
        Teuchos::RCP<Objective<FadFadType>> obj_;
        Teuchos::RCP<Objective<FadType>>    grad_;
        Teuchos::RCP<Objective<Real>>       hess_;

    public:
        ObjectiveAD( const Zakharov<FadFadType> &zakharov ) :
            obj_(Teuchos::rcp(new ObjectiveEnabler<FadFadType>(zakharov))),
            grad_(Teuchos::rcp(new GradientEnabler<FadType>(obj_))),
            hess_(Teuchos::rcp(new HessianEnabler<Real>(grad_))) {}
        

            Real value(const Vector<Real> &x, Real &tol) { return hess_->value(x,tol); }

            void gradient(Vector<Real> &g, const Vector<Real> &x, Real &tol) {
                hess_->gradient(g,x,tol);
            }           

            void hessVec(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
                hess_->hessVec(hv,v,x,tol); 
            }
};




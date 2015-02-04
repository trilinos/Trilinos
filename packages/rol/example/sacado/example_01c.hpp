#include "Sacado.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"

using namespace ROL;

// Nested class variation on example_01

template<class Real>
class Zakharov : public Objective<Real> {

    public:
    Real value(const Vector<Real> &x, Real &tol) {

        Teuchos::RCP<const std::vector<Real> > xp =
            (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();

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
class GradientEnabler : public Objective<Real> {

    typedef Sacado::Fad::DFad<Real> FadType; 

    Teuchos::RCP<Objective<FadType> > obj_;

    public:
    GradientEnabler(Teuchos::RCP<Objective<FadType> > obj) : obj_(obj) {} 

    Real value(const Vector<Real> &x, Real &tol) {

        Teuchos::RCP<const std::vector<Real> > xp =
            (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
 
        int n = xp->size();

        Teuchos::RCP<std::vector<FadType> > x_fad_rcp = Teuchos::rcp(new std::vector<FadType>);
        x_fad_rcp->reserve(n);
        for(int i=0; i<n; ++i) {
            x_fad_rcp->push_back((*xp)[i]);
        }
        StdVector<FadType> x_fad(x_fad_rcp);

        FadType tol_fad = tol;  
        FadType J_fad = obj_->value(x_fad,tol_fad);
            
        return J_fad.val();
    }

    void gradient(Vector<Real> &g, const Vector<Real> &x, Real &tol) {

        Teuchos::RCP<const std::vector<Real> > xp =
            (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();

        Teuchos::RCP<std::vector<Real> > gp =
            Teuchos::rcp_const_cast<std::vector<Real> > ((Teuchos::dyn_cast<StdVector<Real> > (g)).getVector());
    
        int n = xp->size();
 
        Teuchos::RCP<std::vector<FadType> > x_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
        x_fad_rcp->reserve(n);
   
        for(int i=0; i<n; ++i) {
            x_fad_rcp->push_back(FadType(n,i,(*xp)[i])); 
        }

        StdVector<FadType> x_fad(x_fad_rcp);

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

    Teuchos::RCP<Objective<FadType> > obj_;

    public: 
        HessianEnabler(Teuchos::RCP<Objective<FadType> > obj) : obj_(obj) {}

        Real value(const Vector<Real> &x, Real &tol) {

            Teuchos::RCP<const std::vector<Real> > xp =
                (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
 
            int n = xp->size();

            Teuchos::RCP<std::vector<FadType> > x_fad_rcp = Teuchos::rcp(new std::vector<FadType>);
            x_fad_rcp->reserve(n);
            for(int i=0; i<n; ++i) {
                 x_fad_rcp->push_back((*xp)[i]);
            }
            StdVector<FadType> x_fad(x_fad_rcp);

            FadType tol_fad = tol;  
            FadType J_fad = obj_->value(x_fad,tol_fad);
            
            return J_fad.val();
        }

        void gradient(Vector<Real> &g, const Vector<Real> &x, Real &tol) {
            
            Teuchos::RCP<const std::vector<Real> > xp =
                (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
            
            Teuchos::RCP<std::vector<Real> > gp =
                Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(g)).getVector());

            int n = xp->size();

            Teuchos::RCP<std::vector<FadType> > x_fad_rcp = Teuchos::rcp(new std::vector<FadType>);
            Teuchos::RCP<std::vector<FadType> > g_fad_rcp = Teuchos::rcp(new std::vector<FadType>);

            x_fad_rcp->reserve(n);
            g_fad_rcp->reserve(n);

            for(int i=0; i<n; ++i) {
                 x_fad_rcp->push_back((*xp)[i]);
                 g_fad_rcp->push_back((*xp)[i]);
            }
            StdVector<FadType> x_fad(x_fad_rcp);
            StdVector<FadType> g_fad(g_fad_rcp);
            FadType tol_fad = tol;

            obj_->gradient(g_fad,x_fad,tol_fad);

            for(int i=0; i<n; ++i) {
                (*gp)[i] = (*g_fad_rcp)[i].val();    
            }
        }
     
 
        void hessVec(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
            Teuchos::RCP<const std::vector<Real> > xp =
                (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
            Teuchos::RCP<const std::vector<Real> > vp =
                (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
            Teuchos::RCP<std::vector<Real> > hvp =
                Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());

            int n = xp->size(); 

            // Create a vector of independent variables
            Teuchos::RCP<std::vector<FadType> > x_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
            x_fad_rcp->reserve(n);
 
            // Allocate for gradient   
            Teuchos::RCP<std::vector<FadType> > g_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
            g_fad_rcp->reserve(n);

            for(int i=0; i<n; ++i) {
                x_fad_rcp->push_back(FadType(1,(*xp)[i]));
                g_fad_rcp->push_back(0);
                (*x_fad_rcp)[i].fastAccessDx(0) = (*vp)[i];
            }

            StdVector<FadType> x_fad(x_fad_rcp);
            StdVector<FadType> g_fad(g_fad_rcp);
            FadType tol_fad = tol; 
            obj_->gradient(g_fad,x_fad,tol_fad);

            for(int i=0; i<n; ++i) {
                (*hvp)[i] = (*g_fad_rcp)[i].dx(0); 
            } 
        }
       
};





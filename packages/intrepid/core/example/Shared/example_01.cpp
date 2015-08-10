
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_GlobalMPISession.hpp"
using namespace std;
using namespace Intrepid;
 
 int main(){
	 Teuchos::GlobalMPISession mpiSession();
	typedef ArrayTools art; 
	typedef RealSpaceTools<double> rst;
	    const int c=50, p=90, f=70, d1=7, d2=13;

      FieldContainer<double> in_c_f_p(c, f, p);
      FieldContainer<double> in_c_f_p_d(c, f, p, d1);
      FieldContainer<double> in_c_f_p_d_d(c, f, p, d1, d2);
      FieldContainer<double> data_c_p(c, p);
      FieldContainer<double> datainv_c_p(c, p);
      FieldContainer<double> data_c_1(c, 1);
      FieldContainer<double> datainv_c_1(c, 1);
      FieldContainer<double> out_c_f_p(c, f, p);
      FieldContainer<double> outi_c_f_p(c, f, p);
      FieldContainer<double> out_c_f_p_d(c, f, p, d1);
      FieldContainer<double> outi_c_f_p_d(c, f, p, d1);
      FieldContainer<double> out_c_f_p_d_d(c, f, p, d1, d2);
      FieldContainer<double> outi_c_f_p_d_d(c, f, p, d1, d2);
      double zero = INTREPID_TOL*10000.0;
      
      

      // fill with random numbers
      for (int i=0; i<in_c_f_p.size(); i++) {
        in_c_f_p[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_c_f_p_d.size(); i++) {
        in_c_f_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_c_f_p_d_d.size(); i++) {
        in_c_f_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<data_c_p.size(); i++) {
        data_c_p[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_p[i] = 1.0 / data_c_p[i];
      }
      for (int i=0; i<data_c_1.size(); i++) {
        data_c_1[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_1[i] = 1.0 / data_c_1[i];
      }
   	 art::scalarMultiplyDataField<double>(out_c_f_p, data_c_p, in_c_f_p);

      art::scalarMultiplyDataField<double>(outi_c_f_p, datainv_c_p, out_c_f_p);
      rst::subtract(&outi_c_f_p[0], &in_c_f_p[0], outi_c_f_p.size());
      if (rst::vectorNorm(&outi_c_f_p[0], outi_c_f_p.size(), NORM_ONE) > zero) {
      
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d, data_c_p, in_c_f_p_d);
      art::scalarMultiplyDataField<double>(outi_c_f_p_d, datainv_c_p, out_c_f_p_d);
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
       
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d_d, data_c_p, in_c_f_p_d_d);
      art::scalarMultiplyDataField<double>(outi_c_f_p_d_d, datainv_c_p, out_c_f_p_d_d);
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
       
      }
      art::scalarMultiplyDataField<double>(out_c_f_p_d_d, data_c_1, in_c_f_p_d_d);
      art::scalarMultiplyDataField<double>(outi_c_f_p_d_d, datainv_c_1, out_c_f_p_d_d);
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
     
      }
	
	
	
	

      FieldContainer<double> in_c_p(c, p);
      FieldContainer<double> in_c_p_d(c, p, d1);
      FieldContainer<double> in_c_p_d_d(c, p, d1, d2);
      FieldContainer<double> out_c_p(c, p);
      FieldContainer<double> outi_c_p(c, p);
      FieldContainer<double> out_c_p_d(c, p, d1);
      FieldContainer<double> outi_c_p_d(c, p, d1);
      FieldContainer<double> out_c_p_d_d(c, p, d1, d2);
      FieldContainer<double> outi_c_p_d_d(c, p, d1, d2);
  

      // fill with random numbers
      for (int i=0; i<in_c_p.size(); i++) {
        in_c_p[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_c_p_d.size(); i++) {
        in_c_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_c_p_d_d.size(); i++) {
        in_c_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<data_c_p.size(); i++) {
        data_c_p[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_p[i] = 1.0 / data_c_p[i];
      }
      for (int i=0; i<data_c_1.size(); i++) {
        data_c_1[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_1[i] = 1.0 / data_c_1[i];
      }

      art::scalarMultiplyDataData<double>(out_c_p, data_c_p, in_c_p, true);

      art::scalarMultiplyDataData<double>(outi_c_p, datainv_c_p, out_c_p, true);
      rst::subtract(&outi_c_p[0], &in_c_p[0], outi_c_p.size());
  
      art::scalarMultiplyDataData<double>(out_c_p_d, data_c_p, in_c_p_d, true);
      art::scalarMultiplyDataData<double>(outi_c_p_d, datainv_c_p, out_c_p_d, true);	 
	 
	 
	 }


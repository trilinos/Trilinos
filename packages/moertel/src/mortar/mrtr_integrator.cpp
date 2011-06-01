/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Glen Hansen (Glen.Hansen@inl.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#include "mrtr_integrator.H"
#include "mrtr_node.H"
#include "mrtr_segment.H"
#include "mrtr_interface.H"
#include "mrtr_utils.H"

#include <iostream>
#include <fstream>


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 07/05|
 *----------------------------------------------------------------------*/
MOERTEL::Integrator::Integrator(int ngp, bool oneD, int outlevel) :
oneD_(oneD),
ngp_(ngp),
outputlevel_(outlevel)
{
  if (oneD)
  {
    coords_.resize(ngp_);
    weights_.resize(ngp_);
    switch (ngp_)
    {
      case 1:
        coords_[0]  = 0.0;
        weights_[0] = 2.0; 
      break;
      case 2:
        coords_[0]  = -sqrt((1./3.));
        coords_[1]  = -coords_[0];
        weights_[0] = 1.0; 
        weights_[1] = 1.0; 
      break;
      case 3:
        coords_[0]  = -0.77459667;
        coords_[1]  = 0.0;
        coords_[2]  = 0.77459667;
        weights_[0] = 0.55555555; 
        weights_[1] = 0.88888889; 
        weights_[2] = 0.55555555; 
      break;
      case 4:
        coords_[0]  = -0.86113631;
        coords_[1]  = -0.33998104;
        coords_[2]  = 0.33998104;
        coords_[3]  = 0.86113631;
        weights_[0] = 0.34785485; 
        weights_[1] = 0.65214515; 
        weights_[2] = 0.65214515; 
        weights_[3] = 0.34785485; 
      break;
      case 5:
        coords_[0]  = -0.90617985;
        coords_[1]  = -0.53846931;
        coords_[2]  = 0.0;
        coords_[3]  = 0.53846931;
        coords_[4]  = 0.90617985;
        weights_[0] = 0.23692689; 
        weights_[1] = 0.47862867; 
        weights_[2] = 0.56888889; 
        weights_[3] = 0.47862867; 
        weights_[4] = 0.23692689; 
      break;
      case 6:
        coords_[0]  = -0.93246951;
        coords_[1]  = -0.66120939;
        coords_[2]  = -0.23861918;
        coords_[3]  = 0.23861918;
        coords_[4]  = 0.66120939;
        coords_[5]  = 0.93246951;
        weights_[0] = 0.17132449; 
        weights_[1] = 0.36076157; 
        weights_[2] = 0.46791393; 
        weights_[3] = 0.46791393; 
        weights_[4] = 0.36076157; 
        weights_[5] = 0.17132449; 
      break;
      case 7:
        coords_[0]  = -0.94910791;
        coords_[1]  = -0.74153119;
        coords_[2]  = -0.40584515;
        coords_[3]  = 0.0;
        coords_[4]  = 0.40584515;
        coords_[5]  = 0.74153119;
        coords_[6]  = 0.94910791;
        weights_[0] = 0.12948497; 
        weights_[1] = 0.27970539; 
        weights_[2] = 0.38183005; 
        weights_[3] = 0.41795918; 
        weights_[4] = 0.38183005; 
        weights_[5] = 0.27970539; 
        weights_[6] = 0.12948497; 
      break;
      case 8:
        coords_[0]  = -0.96028986;
        coords_[1]  = -0.79666648;
        coords_[2]  = -0.52553241;
        coords_[3]  = -0.18343464;
        coords_[4]  = 0.18343464;
        coords_[5]  = 0.52553241;
        coords_[6]  = 0.79666648;
        coords_[7]  = 0.96028986;
        weights_[0] = 0.10122854; 
        weights_[1] = 0.22238103; 
        weights_[2] = 0.31370665; 
        weights_[3] = 0.36268378; 
        weights_[4] = 0.36268378; 
        weights_[5] = 0.31370665; 
        weights_[6] = 0.22238103; 
        weights_[7] = 0.10122854; 
      break;
      case 10:
        coords_[0]  = -0.97390653;
        coords_[1]  = -0.86506337;
        coords_[2]  = -0.67940957;
        coords_[3]  = -0.43339539;
        coords_[4]  = -0.14887434;
        coords_[5]  = 0.14887434;
        coords_[6]  = 0.43339539;
        coords_[7]  = 0.67940957;
        coords_[8]  = 0.86506337;
        coords_[9]  = 0.97390653;
        weights_[0] = 0.06667134; 
        weights_[1] = 0.14945135; 
        weights_[2] = 0.21908636; 
        weights_[3] = 0.26926672; 
        weights_[4] = 0.29552422; 
        weights_[5] = 0.29552422; 
        weights_[6] = 0.26926672; 
        weights_[7] = 0.21908636; 
        weights_[8] = 0.14945135; 
        weights_[9] = 0.06667134; 
      break;
      default:

	    std::stringstream oss;
			oss << "***ERR*** MOERTEL::Integrator::Integrator:\n"
             << "***ERR*** given number of gaussian points " << ngp_ << "does not exist\n"
             << "***ERR*** use 1, 2, 3, 4, 5, 6, 7, 8, 10 instead\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
		throw ReportError(oss);

      break;
    }
  } // if (oneD)
  else
  {
    coords_.resize(2*ngp_);
    weights_.resize(ngp_);
    switch (ngp_)
    {
      case 3:
        coords_[0]  = 0.5;
        coords_[1]  = 0.0;
        coords_[2]  = 0.5;
        coords_[3]  = 0.5;
        coords_[4]  = 0.0;
        coords_[5]  = 0.5;
        weights_[0] = 1./3.; 
        weights_[1] = weights_[0]; 
        weights_[2] = weights_[0]; 
      break;
      case 6:
      weights_[0] =  0.05497587 ;
      weights_[1] =  0.05497587 ;
      weights_[2] =  0.05497587 ;
      weights_[3] =  0.1116908  ;
      weights_[4] =  0.1116908  ;
      weights_[5] =  0.1116908  ;

      coords_[0 ] =  0.8168476  ;
      coords_[1 ] =  0.09157621 ;
      coords_[2 ] =  0.09157621 ;
      coords_[3 ] =  0.8168476  ;
      coords_[4 ] =  0.09157621 ;
      coords_[5 ] =  0.09157621 ;
      coords_[6 ] =  0.1081030  ;
      coords_[7 ] =  0.4459485  ;
      coords_[8 ] =  0.4459485  ;
      coords_[9 ] =  0.1081030  ;
      coords_[10] =  0.4459485  ;
      coords_[11] =  0.4459485  ;

      break;
      case 12:
      weights_[0 ]  = 0.05839314 ;
      weights_[1 ]  = 0.05839314 ;
      weights_[2 ]  = 0.05839314 ;
      weights_[3 ]  = 0.02542245 ;
      weights_[4 ]  = 0.02542245 ;
      weights_[5 ]  = 0.02542245 ;
      weights_[6 ]  = 0.04142554 ;
      weights_[7 ]  = 0.04142554 ;
      weights_[8 ]  = 0.04142554 ;
      weights_[9 ]  = 0.04142554 ;
      weights_[10]  = 0.04142554 ;
      weights_[11]  = 0.04142554 ;

      coords_[0 ]   = 0.5014265  ;
      coords_[1 ]   = 0.2492867  ;
      coords_[2 ]   = 0.2492867  ;
      coords_[3 ]   = 0.5014265  ;
      coords_[4 ]   = 0.2492867  ;
      coords_[5 ]   = 0.2492867  ;
      coords_[6 ]   = 0.8738220  ;
      coords_[7 ]   = 0.06308901 ;
      coords_[8 ]   = 0.06308901 ;
      coords_[9 ]   = 0.8738220  ;
      coords_[10]   = 0.06308901 ;
      coords_[11]   = 0.06308901 ;
      coords_[12]   = 0.6365025  ;
      coords_[13]   = 0.05314505 ;
      coords_[14]   = 0.6365025  ;
      coords_[15]   = 0.3103525  ;
      coords_[16]   = 0.05314505 ;
      coords_[17]   = 0.6365025  ;
      coords_[18]   = 0.05314505 ;
      coords_[19]   = 0.3103525  ;
      coords_[20]   = 0.3103525  ;
      coords_[21]   = 0.6365025  ;
      coords_[22]   = 0.3103525  ;
      coords_[23]   = 0.05314505 ;


      break;
      case 13:
      weights_[0 ]  =-0.07478502 ;
      weights_[1 ]  = 0.08780763 ;
      weights_[2 ]  = 0.08780763 ;
      weights_[3 ]  = 0.08780763 ;
      weights_[4 ]  = 0.02667362 ;
      weights_[5 ]  = 0.02667362 ;
      weights_[6 ]  = 0.02667362 ;
      weights_[7 ]  = 0.03855688 ;
      weights_[8 ]  = 0.03855688 ;
      weights_[9 ]  = 0.03855688 ;
      weights_[10]  = 0.03855688 ;
      weights_[11]  = 0.03855688 ;
      weights_[12]  = 0.03855688 ;
                               
      coords_[0 ]   = 0.3333333  ;
      coords_[1 ]   = 0.3333333  ;
      coords_[2 ]   = 0.4793081  ;
      coords_[3 ]   = 0.2603460  ;
      coords_[4 ]   = 0.2603460  ;
      coords_[5 ]   = 0.4793081  ;
      coords_[6 ]   = 0.2603460  ;
      coords_[7 ]   = 0.2603460  ;
      coords_[8 ]   = 0.8697398  ;
      coords_[9 ]   = 0.06513010 ;
      coords_[10]   = 0.06513010 ;
      coords_[11]   = 0.8697398  ;
      coords_[12]   = 0.06513010 ;
      coords_[13]   = 0.06513010 ;
      coords_[14]   = 0.6384442  ;
      coords_[15]   = 0.04869032 ;
      coords_[16]   = 0.6384442  ;
      coords_[17]   = 0.3128655  ;
      coords_[18]   = 0.04869032 ;
      coords_[19]   = 0.6384442  ;
      coords_[20]   = 0.04869032 ;
      coords_[21]   = 0.3128655  ;
      coords_[22]   = 0.3128655  ;
      coords_[23]   = 0.6384442  ;
      coords_[24]   = 0.3128655  ;
      coords_[24]   = 0.04869032 ;
                                        
                                        
      break;                            
      case 16:                          
      weights_[0 ] = 0.07215780 ;
      weights_[1 ] = 0.04754582 ;
      weights_[2 ] = 0.04754582 ;
      weights_[3 ] = 0.04754582 ;
      weights_[4 ] = 0.01622925 ;
      weights_[5 ] = 0.01622925 ;
      weights_[6 ] = 0.01622925 ;
      weights_[7 ] = 0.05160869 ;
      weights_[8 ] = 0.05160869 ;
      weights_[9 ] = 0.05160869 ;
      weights_[10] = 0.01361516 ;
      weights_[11] = 0.01361516 ;
      weights_[12] = 0.01361516 ;
      weights_[13] = 0.01361516 ;
      weights_[14] = 0.01361516 ;
      weights_[15] = 0.01361516 ;     
       
      coords_[0 ]  = 0.3333333  ;
      coords_[1 ]  = 0.3333333  ;
      coords_[2 ]  = 0.08141482 ;
      coords_[3 ]  = 0.4592926  ;
      coords_[4 ]  = 0.4592926  ;
      coords_[5 ]  = 0.08141482 ;
      coords_[6 ]  = 0.4592926  ;
      coords_[7 ]  = 0.4592926  ;
      coords_[8 ]  = 0.8989055  ;
      coords_[9 ]  = 0.05054723 ;
      coords_[10]  = 0.05054723 ;
      coords_[11]  = 0.8989055  ;
      coords_[12]  = 0.05054723 ;
      coords_[13]  = 0.05054723 ;
      coords_[14]  = 0.6588614  ;
      coords_[15]  = 0.1705693  ;
      coords_[16]  = 0.1705693  ;
      coords_[17]  = 0.6588614  ;
      coords_[18]  = 0.1705693  ;
      coords_[19]  = 0.1705693  ;
      coords_[20]  = 0.008394777;
      coords_[21]  = 0.7284924  ;
      coords_[22]  = 0.008394777;
      coords_[23]  = 0.2631128  ;
      coords_[24]  = 0.7284924  ;
      coords_[25]  = 0.008394777;
      coords_[26]  = 0.7284924  ;
      coords_[27]  = 0.2631128  ;
      coords_[28]  = 0.2631128  ;
      coords_[29]  = 0.008394777;
      coords_[30]  = 0.2631128  ;
      coords_[31]  = 0.7284924  ;
                                        
                                        
      break;                            
      case 19:                          
      weights_[0 ] = 0.04856790 ;
      weights_[1 ] = 0.01566735 ;
      weights_[2 ] = 0.01566735 ;
      weights_[3 ] = 0.01566735 ;
      weights_[4 ] = 0.03891377 ;
      weights_[5 ] = 0.03891377 ;
      weights_[6 ] = 0.03891377 ;
      weights_[7 ] = 0.03982387 ;
      weights_[8 ] = 0.03982387 ;
      weights_[9 ] = 0.03982387 ;
      weights_[10] = 0.01278884 ;
      weights_[11] = 0.01278884 ;
      weights_[12] = 0.01278884 ;
      weights_[13] = 0.02164177 ;
      weights_[14] = 0.02164177 ;
      weights_[15] = 0.02164177 ;
      weights_[16] = 0.02164177 ;
      weights_[17] = 0.02164177 ;
      weights_[18] = 0.02164177 ;
                              
      coords_[0 ]  = 0.3333333  ;
      coords_[1 ]  = 0.3333333  ;
      coords_[2 ]  = 0.02063496 ;
      coords_[3 ]  = 0.4896825  ;
      coords_[4 ]  = 0.4896825  ;
      coords_[5 ]  = 0.02063496 ;
      coords_[6 ]  = 0.4896825  ;
      coords_[7 ]  = 0.4896825  ;
      coords_[8 ]  = 0.1258208  ;
      coords_[9 ]  = 0.4370896  ;
      coords_[10]  = 0.4370896  ;
      coords_[11]  = 0.1258208  ;
      coords_[12]  = 0.4370896  ;
      coords_[13]  = 0.4370896  ;
      coords_[14]  = 0.6235929  ;
      coords_[15]  = 0.1882035  ;
      coords_[16]  = 0.1882035  ;
      coords_[17]  = 0.6235929  ;
      coords_[18]  = 0.1882035  ;
      coords_[19]  = 0.1882035  ;
      coords_[20]  = 0.9105410  ;
      coords_[21]  = 0.04472951 ;
      coords_[22]  = 0.04472951 ;
      coords_[23]  = 0.9105410  ;
      coords_[24]  = 0.04472951 ;
      coords_[25]  = 0.04472951 ;
      coords_[26]  = 0.03683841 ;
      coords_[27]  = 0.7411986  ;
      coords_[28]  = 0.03683841 ;
      coords_[29]  = 0.2219630  ;
      coords_[30]  = 0.7411986  ;
      coords_[31]  = 0.03683841 ;
      coords_[32]  = 0.7411986  ;
      coords_[33]  = 0.2219630  ;
      coords_[34]  = 0.2219630  ;
      coords_[35]  = 0.03683841 ;
      coords_[36]  = 0.2219630  ;
      coords_[37]  = 0.7411986  ;
                                        
                                        
      break;
      case 27:
      weights_[0 ] = 0.006829866  ;
      weights_[1 ] = 0.006829866  ;
      weights_[2 ] = 0.006829866  ;
      weights_[3 ] = 0.01809227   ;
      weights_[4 ] = 0.01809227   ;
      weights_[5 ] = 0.01809227   ;
      weights_[6 ] = 0.0004635032 ;
      weights_[7 ] = 0.0004635032 ;
      weights_[8 ] = 0.0004635032 ;
      weights_[9 ] = 0.02966149   ;
      weights_[10] = 0.02966149   ;
      weights_[11] = 0.02966149   ;
      weights_[12] = 0.03857477   ;
      weights_[13] = 0.03857477   ;
      weights_[14] = 0.03857477   ;
      weights_[15] = 0.02616856   ;
      weights_[16] = 0.02616856   ;
      weights_[17] = 0.02616856   ;
      weights_[18] = 0.02616856   ;
      weights_[19] = 0.02616856   ;
      weights_[20] = 0.02616856   ;
      weights_[21] = 0.01035383   ;
      weights_[22] = 0.01035383   ;
      weights_[23] = 0.01035383   ;
      weights_[24] = 0.01035383   ;
      weights_[25] = 0.01035383   ;
      weights_[26] = 0.01035383   ;
                                
      coords_[0 ]  = 0.9352701    ;
      coords_[1 ]  = 0.03236495   ;
      coords_[2 ]  = 0.03236495   ;
      coords_[3 ]  = 0.9352701    ;
      coords_[4 ]  = 0.03236495   ;
      coords_[5 ]  = 0.03236495   ;
      coords_[6 ]  = 0.7612982    ;
      coords_[7 ]  = 0.1193509    ;
      coords_[8 ]  = 0.1193509    ;
      coords_[9 ]  = 0.7612982    ;
      coords_[10]  = 0.1193509    ;
      coords_[11]  = 0.1193509    ;
      coords_[12]  =-0.06922210   ;
      coords_[13]  = 0.5346110    ;
      coords_[14]  = 0.5346110    ;
      coords_[15]  =-0.06922210   ;
      coords_[16]  = 0.5346110    ;
      coords_[17]  = 0.5346110    ;
      coords_[18]  = 0.5933802    ;
      coords_[19]  = 0.2033099    ;
      coords_[20]  = 0.2033099   ;
      coords_[21]  = 0.5933802   ;
      coords_[22]  = 0.2033099   ;
      coords_[23]  = 0.2033099   ;
      coords_[24]  = 0.2020614   ;
      coords_[25]  = 0.3989693   ;
      coords_[26]  = 0.3989693   ;
      coords_[27]  = 0.2020614   ;
      coords_[28]  = 0.3989693   ;
      coords_[29]  = 0.3989693   ;
      coords_[30]  = 0.05017814  ;
      coords_[31]  = 0.5932012   ;
      coords_[32]  = 0.05017814  ;
      coords_[33]  = 0.3566206   ;
      coords_[34]  = 0.5932012   ;
      coords_[35]  = 0.05017814  ;
      coords_[36]  = 0.5932012   ;
      coords_[37]  = 0.3566206   ;
      coords_[38]  = 0.3566206   ;
      coords_[39]  = 0.05017814  ;
      coords_[40]  = 0.3566206   ;
      coords_[41]  = 0.5932012   ;
      coords_[42]  = 0.02102202  ;
      coords_[43]  = 0.8074890   ;
      coords_[44]  = 0.02102202  ;
      coords_[45]  = 0.1714890   ;
      coords_[46]  = 0.8074890   ;
      coords_[47]  = 0.02102202  ;
      coords_[48]  = 0.8074890   ;
      coords_[49]  = 0.1714890   ;
      coords_[50]  = 0.1714890   ;
      coords_[51]  = 0.02102202  ;
      coords_[52]  = 0.1714890   ;
      coords_[53]  = 0.8074890   ;   
                                         
      break;                             

      default:

	    std::stringstream oss;
			oss << "***ERR*** MOERTEL::Integrator::Integrator:\n"
             << "***ERR*** given number of gaussian points " << ngp_ << "does not exist\n"
             << "***ERR*** use 3 6 12 13 16 19 27 instead\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
		throw ReportError(oss);

      break;
    }
  }

}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 07/05|
 *----------------------------------------------------------------------*/
MOERTEL::Integrator::~Integrator()
{
  coords_.clear();
  weights_.clear();
}

/*----------------------------------------------------------------------*
 |  integrate a 1D overlap between 2 segments (public)       mwgee 07/05|
 | This function integrates over an overlap of 2 1D segments            |
 | sxia,sxib are start/end local coordinates of the sseg                |
 | mxia,mxib are start/end local coordinates of the mseg                |
 | This method will introduce a transformation into local coordinate eta|
 | to integrate the overlap of the 2 xi-ranges                          |
 | The transformation is                                                |
 | sxi = 0.5*(1-eta)*sxia + 0.5*(1+eta)*sxib                            |
 | mxi = 0.5*(1+eta)*mxia + 0.5*(1-eta)*mxib                            |
 | Note the sign change as orientation of sseg and mseg is in opposite  |
 | directions                                                           |
 | Output is a Epetra_SerialDenseMatrix holding values of integration   |
 | The calling routine is responsible for destroying the                |
 | Epetra_SerialDenseMatrix object                                      |
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix* MOERTEL::Integrator::Integrate(MOERTEL::Segment& sseg, 
                                                      double sxia, double sxib,
                                                      MOERTEL::Segment& mseg, 
                                                      double mxia, double mxib)
{
  int nrow = sseg.Nnode();
  int ncol = mseg.Nnode();
  Epetra_SerialDenseMatrix* Mdense = new Epetra_SerialDenseMatrix(nrow,ncol);

  for (int gp=0; gp<Ngp(); ++gp)
  {
    double eta = Coordinate(gp);
    double wgt = Weight(gp);
    
    // make the transformation for the slave side
    double sxi = 0.5*(1-eta)*sxia + 0.5*(1+eta)*sxib;
    
    // make the transformation for the master side 
    // (note master side is xi positiv in the other direction
    double mxi = 0.5*(1+eta)*mxia + 0.5*(1-eta)*mxib;
    
    // compute the Jacobian dsxi / deta on the slave side
    double dxideta = -0.5*sxia + 0.5*sxib;
    
    // evaluate the Jacobian dx / dsxi (metric) on the slave side
    double dxdsxi = sseg.Metric(&sxi,NULL,NULL);
    
    // calculate value wgt*dxideta*dxdsxi
    double weight = wgt*dxideta*dxdsxi;
    
    // evaluate function 1 of the slave side (supposed to be the LM function)
    double sval[4];
    sseg.EvaluateFunction(1,&sxi,sval,sseg.Nnode(),NULL);
    
    // evaluate function 0 of the master side (supposed to be the trace function)
    double mval[4];
    mseg.EvaluateFunction(0,&mxi,mval,mseg.Nnode(),NULL);
    
    // loop over nodes of the slave side
    for (int slave=0; slave<sseg.Nnode(); ++slave)
    {
      // loop over nodes of the master side
      for (int master=0; master<mseg.Nnode(); ++master)
      {
        // multiply functions for node slave and node master
        double N1N2 = sval[slave]*mval[master];
        (*Mdense)(slave,master) += (N1N2*weight);
      }
    }
  } // for (int gp=0; gp<Ngp(); ++gp)  


  //cout << *Mdense;

  return Mdense;
}

/*----------------------------------------------------------------------*
 |  assemble the result -Mdense into M (public)              mwgee 08/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Integrator::Assemble(MOERTEL::Interface& inter, 
                                   MOERTEL::Segment& sseg, 
                                   MOERTEL::Segment& mseg, 
                                   Epetra_CrsMatrix& M, 
                                   Epetra_SerialDenseMatrix& Mdense)
{
  MOERTEL::Node** snodes = sseg.Nodes();
  MOERTEL::Node** mnodes = mseg.Nodes();
  
  for (int slave=0; slave<sseg.Nnode(); ++slave)
  {
    // only do slave node rows that belong to this proc
    if (inter.NodePID(snodes[slave]->Id()) != inter.lComm()->MyPID())
      continue;
    
    // we want to add the row Mdense(slave,...) to the rows lmdof[sdof] 
    // get the dofs of slave node snodes[slave];
    int        snlmdof = snodes[slave]->Nlmdof();
    const int* slmdof  = snodes[slave]->LMDof();
    
    // this slave node might not have a projection, then the number
    // of lagrange multipliers snlmdof of it is zero
    // in this case, do nothing
    if (!snlmdof) continue;
    
    // loop nodes on master segment
    for (int master=0; master<mseg.Nnode(); ++master)
    {
      // do not add a zero from (*Mdense)(slave,master)
      double val = -(Mdense(slave,master));
      if (abs(val)<1.e-10) continue;
      
      int mndof = mnodes[master]->Ndof();
      const int* mdof = mnodes[master]->Dof();
      
      if (mndof != snlmdof) {

	    std::stringstream oss;
			oss << "***ERR*** MOERTEL::Integrator::Assemble:\n"
				<< "***ERR*** mismatch in number of Lagrange multipliers and primal degrees of freedom:\n"
				<< "***ERR*** slave node " << snodes[slave]->Id()
				<< " master node " << mnodes[master]->Id() << "\n"
				<< "***ERR*** # Lagrange multipliers " << snlmdof << " # dofs " << mndof << "\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
		throw ReportError(oss);

      }
      
      // loop dofs on slave node and insert a value for each master dof
      for (int i=0; i<snlmdof; ++i)
      {
        int row = slmdof[i];
        int col = mdof[i];
        int err = M.SumIntoGlobalValues(row,1,&val,&col);

        if (err)

          err = M.InsertGlobalValues(row,1,&val,&col);

        if (err<0) {

			std::stringstream oss;
				oss << "***ERR*** MOERTEL::Interface::Integrate_2D_Section:\n"
               << "***ERR*** Epetra_CrsMatrix::InsertGlobalValues returned an error\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
			throw ReportError(oss);

        }
      } // for (int i=0; i<snlmdof; ++i)
    } // for (int master=0; master<mseg.Nnode(); ++master)
  } // for (int slave=0; slave<sseg.Nnode(); ++slave)
  return true;
}

/*----------------------------------------------------------------------*
 |  assemble the result Ddense into D (public)               mwgee 08/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Integrator::Assemble(MOERTEL::Interface& inter, 
                                   MOERTEL::Segment& sseg, 
                                   Epetra_CrsMatrix& D, 
                                   Epetra_SerialDenseMatrix& Ddense)
{
  MOERTEL::Node** snodes = sseg.Nodes();
  
  for (int rownode=0; rownode<sseg.Nnode(); ++rownode)
  {
    // only insert in rows that I own
    if (inter.NodePID(snodes[rownode]->Id()) != inter.lComm()->MyPID())
      continue;
      
    // get row dofs
    int        nlmdof = snodes[rownode]->Nlmdof();
    const int* lmdof  = snodes[rownode]->LMDof();
    
    // this slave node might not have a projection and therefore might not
    // carry lagrange multipliers. In this case, do not insert anything
    if (nlmdof==0) continue;
    
    // loop column nodes
    for (int colnode=0; colnode<sseg.Nnode(); ++colnode)
    {
      // do not add a zero from Ddense
      double val = Ddense(rownode,colnode);
      if (abs(val)<1.e-10) continue;
      
      int ndof = snodes[colnode]->Ndof();
      const int* dof = snodes[colnode]->Dof();
      
      if (nlmdof != ndof) {

		std::stringstream oss;
			oss << "***ERR*** MOERTEL::Interface::Integrate_2D_Section:\n"
				<< "***ERR*** mismatch in number of Lagrange multipliers and primal degrees of freedom:\n"
				<< "***ERR*** slave node " << snodes[rownode]->Id()
				<< "***ERR*** # Lagrange multipliers " << nlmdof << " # dofs " << ndof << "\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
		throw ReportError(oss);

      }
      
      // loop lm dofs and insert a value for each dof
      for (int i=0; i<nlmdof; ++i)
      {
        int row = lmdof[i];
        int col = dof[i];
        int err = D.SumIntoGlobalValues(row,1,&val,&col);

        if (err)

          err = D.InsertGlobalValues(row,1,&val,&col);

        if (err<0) {

			std::stringstream oss;
				oss << "***ERR*** MOERTEL::Interface::Integrate_2D_Section:\n"
               << "***ERR*** Epetra_CrsMatrix::SumIntoGlobalValues returned an error\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
			throw ReportError(oss);

        }
      } // for (int i=0; i<nlmdof; ++i)
    } // for (int colnode=0; colnode<sseg.Nnode(); ++colnode)
  } // for (int rownode=0; rownode<sseg.Nnode(); ++rownode)

  return true;
}


/*----------------------------------------------------------------------*
 |  integrate a 1D slave segment (public)                    mwgee 07/05|
 |  This method integrates 2 functions on the same (slave) segment      |
 |  from given local coordinates sxia to sxib                           |
 | Output is a Epetra_SerialDenseMatrix holding values of integration   |
 | The calling routine is responsible for destroying the                |
 | Epetra_SerialDenseMatrix object                                      |
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix* MOERTEL::Integrator::Integrate(MOERTEL::Segment& sseg, 
                                                      double sxia, double sxib)
{
  int nrow = sseg.Nnode();
  int ncol = nrow;
  Epetra_SerialDenseMatrix* Ddense = new Epetra_SerialDenseMatrix(nrow,ncol);

  for (int gp=0; gp<Ngp(); ++gp)
  {
    double eta = Coordinate(gp);
    double wgt = Weight(gp);
    
    // make the coordinate transformation
    double sxi = 0.5*(1-eta)*sxia + 0.5*(1+eta)*sxib;
    
    // compute the Jacobian dsxi / deta
    double dxideta = -0.5*sxia + 0.5*sxib;
    
    // evaluate the Jacobian dx / dsxi (metric)
    double dxdsxi = sseg.Metric(&sxi,NULL,NULL);
    
    // calculate value wgt*dxideta*dxdsxi
    double weight = wgt*dxideta*dxdsxi;

    // evaluate function 1 (supposed to be the LM function)
    double lmval[4];
    sseg.EvaluateFunction(1,&sxi,lmval,sseg.Nnode(),NULL);
    
    // evaluate function 0 (supposed to be the trace space function)
    double val[4];
    sseg.EvaluateFunction(0,&sxi,val,sseg.Nnode(),NULL);    
    
    // loop over all nodes (lm loop)
    for (int lm=0; lm<sseg.Nnode(); ++lm)
    {
      // loop over all nodes (dof loop)
      for (int dof=0; dof<sseg.Nnode(); ++dof)
      {
        // multiply the 2 functions
        double N1N2 = lmval[lm]*val[dof];
        (*Ddense)(lm,dof) += (N1N2*weight);
      }
    }
  } // for (int gp=0; gp<Ngp(); ++gp)

  //cout << *Ddense;

  return Ddense;
}




/*----------------------------------------------------------------------*
 |                                            (public)       mwgee 08/05|
 | integrate the modification of the master side                        |
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix* MOERTEL::Integrator::Integrate_2D_Mmod(
                                                      MOERTEL::Segment& sseg, 
                                                      double sxia, double sxib,
                                                      MOERTEL::Segment& mseg, 
                                                      double mxia, double mxib)
{
  Epetra_SerialDenseMatrix* Mmod = new Epetra_SerialDenseMatrix(mseg.Nnode(),1);

  for (int gp=0; gp < Ngp(); ++gp)
  {
    double eta = Coordinate(gp);
    double wgt = Weight(gp);
    
    // make the transformation for the slave side
    double sxi = 0.5*(1-eta)*sxia + 0.5*(1+eta)*sxib;
    
    // make the transformation for the master side 
    // (note master side is xi positiv in the other direction
    double mxi = 0.5*(1+eta)*mxia + 0.5*(1-eta)*mxib;
    
    // compute the Jacobian dsxi / deta on the slave side
    double dxideta = -0.5*sxia + 0.5*sxib;
    
    // evaluate the Jacobian dx / dsxi (metric) on the slave side
    double dxdsxi = sseg.Metric(&sxi,NULL,NULL);
    
    // calculate value wgt*dxideta*dxdsxi
    double weight = wgt*dxideta*dxdsxi;
    
    // evaluate function 0 of the slave side (supposed to be the trace function)
    double sval[4];
    sseg.EvaluateFunction(0,&sxi,sval,sseg.Nnode(),NULL);
    
    // make the delta function phi12 = phi1 - phi2
    double val = sval[0] - sval[1];
    
    // evaluate function 0 of the master side (supposed to be the trace function)
    double mval[4];
    mseg.EvaluateFunction(0,&mxi,mval,mseg.Nnode(),NULL);
    
    // loop over nodes of the master side
    for (int master=0; master<mseg.Nnode(); ++master)
    {
      // multiply functions for each master node
      double N1N2 = -0.5 * val * mval[master];
      (*Mmod)(master,0) += (N1N2*weight);
    }
  } // for (int gp=0; gp<integrator.Ngp(); ++gp)  


  //cout << *Mmod;

  return Mmod;
}

#if 0 // old version
/*----------------------------------------------------------------------*
 |  assemble the modification -Mmod into M (public)          mwgee 08/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Integrator::Assemble_2D_Mod(MOERTEL::Interface& inter, 
                                          MOERTEL::Segment& sseg, 
                                          MOERTEL::Segment& mseg, 
                                          Epetra_CrsMatrix& M, 
                                          Epetra_SerialDenseMatrix& Mmod)
{
  MOERTEL::Node** snodes = sseg.Nodes();
  MOERTEL::Node** mnodes = mseg.Nodes();
  
  for (int slave=0; slave<sseg.Nnode(); ++slave)
  {
    // only do slave node rows that belong to this proc
    if (inter.NodePID(snodes[slave]->Id()) != inter.lComm()->MyPID())
      continue;
    
    // we want to add the row Mdense(slave,...) to the rows lmdof[sdof] 
    // and                row Ddense(slave,...) to the rows lmdof[sdof] 
    // get the dofs of slave node snodes[slave];
    int        snlmdof = snodes[slave]->Nlmdof();
    const int* slmdof  = snodes[slave]->LMDof();
    
    // this slave node might not have a projection, then the number
    // of lagrange multipliers snlmdof of it is zero
    // in this case, do nothing
    if (!snlmdof) continue;
    
    // loop lm dofs on node slave
    for (int sdof=0; sdof<snlmdof; ++sdof)
    {
      int row = slmdof[sdof];

      // loop nodes on master segment
      for (int master=0; master<mseg.Nnode(); ++master)
      {
        int mndof = mnodes[master]->Ndof();
        const int* mdofs = mnodes[master]->Dof();
        
        // dofs on node master
        for (int mdof=0; mdof<mndof; ++mdof)
        {
          int col = mdofs[mdof];
          
          // do not add a zero from (*Mdense)(slave,master)
          double val = -(Mmod(slave*snlmdof+sdof,master*mndof+mdof));
          if (abs(val)<1.e-9) continue;
          
          int err = M.SumIntoGlobalValues(row,1,&val,&col);
          if (err)
            err = M.InsertGlobalValues(row,1,&val,&col);

          if (err<0) {

			throw ReportError(
				string("***ERR*** MOERTEL::Interface::Assemble_2D_Mod:\n") +
				string("***ERR*** Epetra_CrsMatrix::SumIntoGlobalValues returned an error\n") +
				string("***ERR*** file/line: ") + string(__FILE__) + "/" + toString(__LINE__) + "\n");

          } // if (err)
        } // for (int mdof=0; mdof<mndof; ++mdof)
      } // for (int master=0; master<mseg.Nnode(); ++master)
    } // for (int sdof=0; sdof<snlmdof; ++sdof)
  } // for (int slave=0; slave<sseg.Nnode(); ++slave)
  return true;
}
#endif


/*----------------------------------------------------------------------*
 |  assemble the modification -Mmod into slave nodes (public)mwgee 08/05|
 |  Note that this is not scalar as the other assemble routines         |
 *----------------------------------------------------------------------*/
bool MOERTEL::Integrator::Assemble_2D_Mod(MOERTEL::Interface& inter, 
                                          MOERTEL::Segment& sseg, 
                                          MOERTEL::Segment& mseg, 
                                          Epetra_SerialDenseMatrix& Mmod)
{
  MOERTEL::Node** snodes = sseg.Nodes();
  MOERTEL::Node** mnodes = mseg.Nodes();
  
  for (int slave=0; slave<sseg.Nnode(); ++slave)
  {
    // only do slave node rows that belong to this proc
    if (inter.NodePID(snodes[slave]->Id()) != inter.lComm()->MyPID())
      continue;
    
    // we want to add the row Mdense(slave,...) to the rows lmdof[sdof] 
    // and                row Ddense(slave,...) to the rows lmdof[sdof] 
    // get the dofs of slave node snodes[slave];
    
    // what we would like to do is
    //int        snlmdof = snodes[slave]->Nlmdof();
    //const int* slmdof  = snodes[slave]->LMDof();

    // what we temporarily do is
    int        snlmdof = snodes[slave]->Ndof();
    // as we don't know the LMdofs yet, we just know that their number is
    // equal to that of the primary dofs
    
    // this slave node might not have a projection, then the number
    // of lagrange multipliers snlmdof of it is zero
    // in this case, do nothing
    //if (!snlmdof) continue;
    // this has to be decided later on
    
    // loop lm dofs on node slave
    for (int sdof=0; sdof<snlmdof; ++sdof)
    {
      // loop nodes on master segment
      for (int master=0; master<mseg.Nnode(); ++master)
      {
        int mndof = mnodes[master]->Ndof();
        const int* mdofs = mnodes[master]->Dof();
        
        // dofs on node master
        for (int mdof=0; mdof<mndof; ++mdof)
        {
          int col = mdofs[mdof];
          
          // do not add a zero from (*Mdense)(slave,master)
          double val = -(Mmod(slave*snlmdof+sdof,master*mndof+mdof));
          if (abs(val)<1.e-9) continue;
          
          // standard assembly for internal nodes
          if (!snodes[slave]->IsOnBoundary())
            snodes[slave]->AddMmodValue(sdof,val,col);
          // if slave node is on boundary
          else
          {
			std::map<int,MOERTEL::Node*>& suppnodes = snodes[slave]->GetSupportedByNode();
            double w = 1./(double)(snodes[slave]->NSupportSet());
			std::map<int,MOERTEL::Node*>::iterator curr;
            for (curr=suppnodes.begin(); curr!=suppnodes.end(); ++curr)
              curr->second->AddMmodValue(sdof,w*val,col);
          }
          
        } // for (int mdof=0; mdof<mndof; ++mdof)
      } // for (int master=0; master<mseg.Nnode(); ++master)
    } // for (int sdof=0; sdof<snlmdof; ++sdof)
  } // for (int slave=0; slave<sseg.Nnode(); ++slave)
  return true;
}



/*----------------------------------------------------------------------*
 |  integrate a 2D triangle/quad overlap segment (public)    mwgee 11/05|
 |  contribution from the master/slave side M/D                         |
 *----------------------------------------------------------------------*/
bool MOERTEL::Integrator::Integrate(Teuchos::RCP<MOERTEL::Segment> actseg,
                                    MOERTEL::Segment& sseg, 
                                    MOERTEL::Segment& mseg, 
                                    Epetra_SerialDenseMatrix** Ddense, 
                                    Epetra_SerialDenseMatrix** Mdense, 
                                    MOERTEL::Overlap& overlap, double eps,
                                    bool exactvalues)
{
//	static int cnt = 0;

  if (oneD_) {

	std::stringstream oss;
		oss << "***ERR*** MOERTEL::Integrator::Integrate:\n"
         << "***ERR*** Integrator was not constructed for 2D integration\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  
  // we integrate the scalar function here
  int nrow = sseg.Nnode();
  int ncol = mseg.Nnode();
  // get the points
  const int np = actseg->Nnode(); // this is an overlap segment - always a triangle
  const int* nodeid = actseg->NodeIds();
  std::vector<MOERTEL::Point*> points;
  overlap.PointView(points,nodeid,np);

  
  // get the area of the overlap segment
  double area = actseg->Area();
  if (area<0.0)
  {
    if (OutLevel()>3)
    cout << "MOERTEL: ***ERR***  MOERTEL::Integrator::Integrate:\n" << "MOERTEL: ***ERR***  overlap segment area is negative: " << area << endl
         << "MOERTEL: ***ERR***  skipping....\n"
         << "MOERTEL: ***ERR***  file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  // get the area of the slave segment
  double sarea = sseg.Area();
  
  if (abs(area/sarea)<eps)
  {
    if (OutLevel()>10)
      cout << "MOERTEL: ***WRN*** Skipping overlap segment with tiny area " << area << endl;
    points.clear();
    return false;  
  }

  // for all integrations other then 3 its dx*dy = 2A dxi*deta
  if (Ngp()!=3) 
    area *= 2.0;

  *Mdense = new Epetra_SerialDenseMatrix(nrow,ncol);
  *Ddense = new Epetra_SerialDenseMatrix(nrow,nrow);

  //========================================================================
  //========================================================================
  // interpolate values at gaussian points with shape functions from actsseg
  std::vector< std::vector<double>* > vals;
  if (!exactvalues)
  {
    // get the function values at the points
    // slave segment function 0 is vals[point][0][0...np-1]
    // slave segment function 1 is vals[point][1][0...np-1]
    // master segment function 0 is vals[point][2][0...np-1]
    vals.resize(np);
    for (int i=0; i<np; ++i)
      vals[i] = points[i]->FunctionValues();

    // loop integration points
    for (int gp=0; gp<Ngp(); ++gp)
    {
      double* eta   = Coordinate(&gp);
      double weight = area*Weight(gp);
      
      
      // evaluate the linear function 0 of the actseg at gaussian point
      double val[20];
      actseg->EvaluateFunction(0,eta,val,actseg->Nnode(),NULL);

      // interpolate shape function 0 from sseg
      // interpolate shape function 1 from sseg
      double val_sfunc0[20];
      double val_sfunc1[20];
      double val_mfunc0[20];
      for (int i=0; i<np; ++i) 
      {
        val_sfunc0[i] = 0.0;
        val_sfunc1[i] = 0.0;
        val_mfunc0[i] = 0.0;
      }
      for (int point=0; point<np; ++point)
      {
        for (int i=0; i<nrow; ++i) // could be either 3 or 4 values
        {
          val_sfunc0[i] += val[point]*vals[point][0][i];
          val_sfunc1[i] += val[point]*vals[point][1][i];
        }
        for (int i=0; i<ncol; ++i) // could be either 3 or 4 values
          val_mfunc0[i] += val[point]*vals[point][2][i];
      }
      
      //cout << val_sfunc0[0] << " " << val_sfunc0[1] << " " << val_sfunc0[2] << " " << val_sfunc0[3] << endl;
      //cout << val_sfunc1[0] << " " << val_sfunc1[1] << " " << val_sfunc1[2] << " " << val_sfunc1[3] << endl;
      //cout << val_mfunc0[0] << " " << val_mfunc0[1] << " " << val_mfunc0[2] << " " << val_mfunc0[3] << endl;

      // loop over all nodes (lm loop)
      for (int lm=0; lm<sseg.Nnode(); ++lm)
      {
        // loop over all nodes (dof loop master)
        for (int dof=0; dof<mseg.Nnode(); ++dof)
        {
          // multiply the 2 functions
          double N1N2 = val_sfunc1[lm]*val_mfunc0[dof];
          (**Mdense)(lm,dof) += (N1N2*weight);
          //cout << "Adding gausspoint M value M(" << lm << "," << dof << ") += " << val_sfunc1[lm] << " * " << val_mfunc0[dof] << " * " << weight << endl;
        }
        
        // loop over all nodes (dof loop slave)
        for (int dof=0; dof<sseg.Nnode(); ++dof)
        {
          // multiply the 2 functions
          double N1N2 = val_sfunc1[lm]*val_sfunc0[dof];
          (**Ddense)(lm,dof) += (N1N2*weight);
          //cout << "Adding gausspoint D value D(" << lm << "," << dof << ") += " << val_sfunc1[lm] << " * " << val_sfunc0[dof] << " * " << weight << endl;
        }
      }
      
    } // for (int gp=0; gp<Ngp(); ++gp)
    //cout << **Ddense;
    //cout << **Mdense;
    
    vals.clear();
    points.clear();
  } // if (!exactvalues)
  //========================================================================
  //========================================================================
  // use exact values at gaussian points
  else
  {
    // compute local coordinates of actseg's nodes in slave element
    double psxi[3][2];
    for (int i=0; i<np; ++i)
    {
      psxi[i][0] = points[i]->Xi()[0];
      psxi[i][1] = points[i]->Xi()[1];
      //cout << "psxi[" << i << "] = " << psxi[i][0] << " / " << psxi[i][1] << endl; 
    }
    // create a node to use for projection
    double x[3];
    double n[3];
    int dof[3]; dof[0] = dof[1] = dof[2] = -1;
	Teuchos::RCP<MOERTEL::Node> gpnode = Teuchos::rcp(new MOERTEL::Node(-2,x,3,dof,false,OutLevel()));
    
    // create a projector to project gaussian points
    MOERTEL::Projector projector(false,OutLevel());
    
    for (int gp=0; gp<Ngp(); ++gp)
    {
      double* eta    = Coordinate(&gp);
      double  weight = area*Weight(gp);
      
      //-------------------------------------------------------------------
      // compute the local coord of the gaussian point in the slave element
      double val[3];
      actseg->EvaluateFunction(0,eta,val,actseg->Nnode(),NULL);

      double sxi[2]; sxi[0] = sxi[1] = 0.0;
      for (int point=0; point<np; ++point)
      {
        sxi[0] += val[point] * psxi[point][0];
        sxi[1] += val[point] * psxi[point][1];
      }
      //cout << "sxi = " << sxi[0] << " / " << sxi[1] << endl;
      
      //-------------------------------------------------------------------
      // get shape function values from function 0 and 1 from slave element here
      double val_sfunc0[20];
      double val_sfunc1[20];
      sseg.EvaluateFunction(0,sxi,val_sfunc0,sseg.Nnode(),NULL);  
      sseg.EvaluateFunction(1,sxi,val_sfunc1,sseg.Nnode(),NULL);  
      
      //-------------------------------------------------------------------
      // compute the global coordinate of the gaussian point and project it
      // to the master element
      x[0] = x[1] = x[2] = 0.0;
      n[0] = n[1] = n[2] = 0.0;
      for (int point=0; point<np; ++point)
        for (int j=0; j<3; ++j)
        {
          x[j] += val[point] * points[point]->Node()->X()[j];
          n[j] += val[point] * points[point]->Node()->N()[j];
        }
      const double length = MOERTEL::length(n,3);
      for (int j=0; j<3; ++j) n[j] /= length;
      //cout << "x = " << x[0] << " / " << x[1] << " / " << x[2] << endl;
      //cout << "n = " << n[0] << " / " << n[1] << " / " << n[2] << endl;
      gpnode->SetX(x);             
      gpnode->SetN(n);
      double mxi[2];
	  double gap;
      bool ok = projector.ProjectNodetoSegment_NodalNormal(*gpnode,mseg,mxi,gap);
      // if we have problems projecting here, we better skip this gauss point
      if (!ok)
      { 
        cout << "MOERTEL: ***WRN***------------------projection failed in integration\n";
		fflush(stdout);
        continue;
      }
      //cout << "mxi = " << mxi[0] << " / " << mxi[1] << endl;
                   
      //-------------------------------------------------------------------
      // get shape function value from mseg
      double val_mfunc0[20];
      mseg.EvaluateFunction(0,mxi,val_mfunc0,mseg.Nnode(),NULL);  

      //-------------------------------------------------------------------
      // loop over all slave nodes (lm loop)
      for (int lm=0; lm<sseg.Nnode(); ++lm)
      {
        // loop over all nodes (dof loop master)
        for (int dof=0; dof<mseg.Nnode(); ++dof){
          (**Mdense)(lm,dof) += (weight * val_sfunc1[lm] * val_mfunc0[dof]);
        
        // loop over all nodes (dof loop slave)
        for (int dof=0; dof<sseg.Nnode(); ++dof){

          (**Ddense)(lm,dof) += (weight * val_sfunc1[lm] * val_sfunc0[dof]); 

      }
    } // for (int gp=0; gp<Ngp(); ++gp)
    
    // tidy up
    gpnode = Teuchos::null;
    
  } // if (exactvalues)



  return true;
}


/*----------------------------------------------------------------------*
 |  assemble contributions Ddense into slave nodes (public)  mwgee 11/05|
 |  the scalar integrated functions are assembled here                  |
 *----------------------------------------------------------------------*/
bool MOERTEL::Integrator::Assemble(MOERTEL::Interface& inter,MOERTEL::Segment& sseg,
                                   Epetra_SerialDenseMatrix& Ddense)
{
  // get nodes
  const int nnode       = sseg.Nnode();
  MOERTEL::Node** snode = sseg.Nodes();

#if 0
  // set a row of Ddense in each snode
  for (int row=0; row<nnode; ++row)
  {
    // assemble only to my own nodes
    if (inter.NodePID(snode[row]->Id()) != inter.lComm()->MyPID())
      continue;
    
    // row node is internal node
    for (int col=0; col<nnode; ++col)
      snode[row]->AddDValue(Ddense(row,col),snode[col]->Id());
  } // for (int row=0; row<nnode; ++row)
#endif


#if 1
  // set a row of Ddense in each snode
  for (int row=0; row<nnode; ++row)
  {
    // assemble only to my own nodes
    if (inter.NodePID(snode[row]->Id()) != inter.lComm()->MyPID())
      continue;
    
    // row node is internal node
    if (!snode[row]->IsOnBoundary())
    {
      for (int col=0; col<nnode; ++col)
      { 
        // row/col are both internal
        if (!snode[col]->IsOnBoundary())
          snode[row]->AddDValue(Ddense(row,col),snode[col]->Id());
        // row is internal node, col is boundary node
        // As those entries would ruin the diagonal structure of D they are
        // simply assembled into M (which is not diagonal anyway)
        else
          snode[row]->AddMValue(Ddense(row,col),snode[col]->Id());
      }
    }
    // row node is a boundary node
    else
    {
	  std::map<int,MOERTEL::Node*>& suppnodes = snode[row]->GetSupportedByNode();
      double w = 1./(double)(snode[row]->NSupportSet());
	  std::map<int,MOERTEL::Node*>::iterator curr;
      for (curr=suppnodes.begin(); curr!=suppnodes.end(); ++curr)
        for (int col=0; col<nnode; ++col)
        { 
          // col node is internal, assemble into D
          if (!snode[col]->IsOnBoundary())
            curr->second->AddDValue((w*Ddense(row,col)),snode[col]->Id());
          // col node is boundary, assemble into M
          else
            curr->second->AddMValue((w*Ddense(row,col)),snode[col]->Id());
        }
    }
  } // for (int row=0; row<nnode; ++row)
#endif

  return true;
}

/*----------------------------------------------------------------------*
 |  assemble contributions Mdense into slave nodes (public)  mwgee 11/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Integrator::Assemble(MOERTEL::Interface& inter,
                                   MOERTEL::Segment& sseg,
                                   MOERTEL::Segment& mseg,
                                   Epetra_SerialDenseMatrix& Mdense)
{
  // get nodes
  const int nsnode    = sseg.Nnode();
  const int nmnode    = mseg.Nnode();
  MOERTEL::Node** snode  = sseg.Nodes();
  MOERTEL::Node** mnode  = mseg.Nodes();

#if 0 // the orig version
  // set a row of Mdense in each snode
  for (int row=0; row<nsnode; ++row)
  {
    // assemble only to my own nodes
    if (inter.NodePID(snode[row]->Id()) != inter.lComm()->MyPID()) 
      continue;
      
    // standard assembly for internal nodes
    for (int col=0; col<nmnode; ++col)
      // note the sign change here!!!!
      snode[row]->AddMValue(-Mdense(row,col),mnode[col]->Id());
  }
#endif

#if 1 // the orig version
  // set a row of Mdense in each snode
  for (int row=0; row<nsnode; ++row)
  {
    // assemble only to my own nodes
    if (inter.NodePID(snode[row]->Id()) != inter.lComm()->MyPID()) 
      continue;
      
    // standard assembly for internal nodes
    if (!snode[row]->IsOnBoundary())
      for (int col=0; col<nmnode; ++col)
        // note the sign change here!!!!
        snode[row]->AddMValue((-Mdense(row,col)),mnode[col]->Id());
    else
    {
	  std::map<int,MOERTEL::Node*>& suppnodes = snode[row]->GetSupportedByNode();
      // note the sign change here!!!!
      double w = -1./(double)(snode[row]->NSupportSet());
	  std::map<int,MOERTEL::Node*>::iterator curr;
      for (curr=suppnodes.begin(); curr!=suppnodes.end(); ++curr)
        for (int col=0; col<nmnode; ++col)
          curr->second->AddMValue((w*Mdense(row,col)),mnode[col]->Id());
    }
  }
#endif
  return true;
}

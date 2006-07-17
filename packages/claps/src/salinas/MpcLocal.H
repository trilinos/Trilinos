//@HEADER
// ************************************************************************
//
//         Claps: A Collection of Domain Decomposition Preconditioners
//                and Solvers
//         Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef Salinas_MpcLocal_H
#define Salinas_MpcLocal_H

/* MpcLocal - a class to hold structures for multipoint constraints.

   The Mpc is an equation of the form,

        Sum  C_i * U_i = G
         i

   where the U_i are the displacements of degree of freedom i, and C_i is
   a real coefficient. G is a constant (usually zero).


   While MPCs may be specified in another coordinate system, they are stored
   and passed in the basic coordinate system.


   MpcLocal contains only the subdomain contribution to the MPC. It
   contains the node ID in local coordinates. It contains two entries
   for the length of the MPC. NumEntries is the number of degrees of
   freedom which contact this subdomain. NumEntrieGlobal is the number
   of entries in the global system.

//    In the event that a MPC touches a degree of freedom on a boundary,
//    then the coefficients of that degree of freedom will be specified
//    identically in all the relevant subdomains. The coefficient will be
//    the same value as the global coefficient. For example, if the
//    global coefficient is 1.0 for a dof on a boundary, then each
//    subdomain which contains that node will have a coefficient of
//    1.0.
  CORRECTION
   As implemented, the above is incorrect.
   In the event that a MPC touches a degree of freedom on a boundary,
   then the coefficients of that degree of freedom will be spread
   across processors so that the sum of the coefficients equals the
   total. Thus, if a coefficient in the serial case were 1, but the
   MPC is a boundary node for 2 subdomains, the MpcLocal coefficients
   on the two subdomains will be 0.5.

   Note, all the subdomains in the solution will have the same number
   of MPC equations. However, if a particular MPC equation does not
   apply to a given subdomain, there will be no data associated with
   it (i.e. the NumEntries will be zero). Only the local contributions
   to a multipoint constraint will be included in the MPC list for a
   subdomain.

   The value of the constant, G, will be identical on all subdomains
   which contact the constraint. That value is the value of the global
   constant (i.e. it will not be divided among the subdomains). If
   the constraint has no entries in a subdomain, the value of G on
   that subdomain will be undefined.

   IMPORTANT - INTERFACE DIVISION
   The coefficients on an interface are divided by the cardinality
   in a manner analogous to the division of the force vector. Thus,
   to reconstruct the global coefficients, a  sum of all appropriate
   terms on all processors must be made. This can be overridden.


   Consider a very simple example. In this example, the element on the
   left is in subdomain 0, the element on the right is in subdomain 1.
   Global Node numbers are shown here.

    2---------4--------6            Constraints:
    |         |        |                x(1)=x(2)  and
    |         |        |                x(1)=x(3)
    1---------3--------5

    After  decomposition into subdomains, and with local numbering

       subdomain 0                subdomain-1

       1-----3                    1-----3
       |     |                    |     |
       |     |                    |     |
       0-----2                    0-----2

first MPC (x components of global nodes 1 and 2 are equal)
  NumEntries=2                   NumEntries=0        
  NumEntrieGlobal=2              NumEntrieGlobal=2   
  LocalId={0 1}                  LocalId=NULL       
  NodalDof={0 0}                 NodalDof=NULL     
  Coef={1.0 -1.0}                Coef=NULL    
  G=0                            G=0      

2nd MPC  (x components of global nodes 1 and 3 are equal)         
  NumEntries=2                   NumEntries=1        
  NumEntrieGlobal=2              NumEntrieGlobal=2   
  LocalId={0 2}                  LocalId={0}       
  NodalDof={0 0}                 NodalDof={0}
  Coef={1.0 -0.5}                Coef={-0.5}    
  G=0                            G=0      



*/
//#include <code_types.h>
// This file is included from feti, clop, etc. It must not have any includes.
// for convenience, we define "Real" as a double, but this is internal data only.
#if defined(DOUBLE128)
#define Real long double
#else
#define Real double
#endif


class MpcLocal {
private:
  unsigned int mNumEntries;       // number of entries in the sum (locally)
  unsigned int mNumEntriesGlobal; // number of entries globally.
  int* mLocalToGlobalEntry;  // CU FIX only required for mpc stiffness scaling
  int* mLocalId;                  // the local node ID of U_i
  short* mNodalDof;               // the node based degree of freedom (0:5)
  Real* mCoef;                  // C_i
  Real  mG;

public:
  MpcLocal() {                       // default constructor
    mNumEntries = 0;
    mNumEntriesGlobal = 0;
    mLocalToGlobalEntry = (int *)0; // CU_FIX
    mLocalId   = (int *)0;
    mCoef       = (Real *)0;
    mNodalDof   = (short *)0;
    mG          = 0.0;
  }
  MpcLocal( const MpcLocal& old );        // copy constructor
  MpcLocal& operator=(const MpcLocal& rhs);
    

  ~MpcLocal() {                      // destructor
    if (mNumEntries){
      delete[] mLocalId; mLocalId = 0;
      delete[] mCoef; mCoef = 0;
      delete[] mNodalDof; mNodalDof = 0;
      delete[] mLocalToGlobalEntry; mLocalToGlobalEntry = 0; // CU FIX
    }
  }

  // methods to Set data
  void NumEntries(int n);             // resizes the data - sets to zero
  void NumEntriesGlobal(int n)        // set only the size of global equation
  { mNumEntriesGlobal=n; }
  void LocalId( unsigned int m, int id) { // sets local ID (0:n-1) of entry "
    if ( m<mNumEntries )
      mLocalId[m]=id;
  }
  void Coef( unsigned int m, Real c ) {// sets coefficent of entry "m"
    if ( m<mNumEntries )
      mCoef[m]=c;
  }
  void NodalDof( unsigned int m, short dof ) {// sets nodal dof (0:5) of entry "m
    if ( m<mNumEntries )
      mNodalDof[m]=dof;
  }
  void G( Real g){ mG=g;}            // sets residual of constraint equation
  void LocalToGlobalEntry( unsigned int m, int global_m ) {  // CU FIX: for mpc stiffness scaling
    if ( m<mNumEntries && (mLocalToGlobalEntry))
      mLocalToGlobalEntry[m] = global_m;
  }
  int LocalToGlobalEntry(unsigned int m) const  // CU FIX
  {
     if ( (m<mNumEntries) && (mLocalToGlobalEntry) )
      return mLocalToGlobalEntry[m];
    else
      return 0;
  }

  // methods to retrieve data
  int NumEntries() const               // returns number of local entries
                           // all variables are dimensioned to this value
    { return mNumEntries; }
  int NumEntriesGlobal() const         // returns number of global entries.
  { return mNumEntriesGlobal; }
  int LocalId( unsigned int m) const   // returns local ID (0:n-1) of entry "m"
  {
    if ( m<mNumEntries )
      return mLocalId[m];
    else
      return 0;
  }
  Real Coef( unsigned int m ) const  // returns coefficent of entry "m"
  {
    if ( m<mNumEntries )
      return mCoef[m];
    return 0.0;
  }
  short NodalDof( unsigned int m )const // returns nodal dof (0:5) of entry "m"
  {
    if ( m<mNumEntries )
      return mNodalDof[m];
    return 0;
  }  
  Real G() const { return mG; }      // returns residual of constraint equation


};   // end of class definition

#undef Real
/**********************************************************************
 * class definition created by Garth Reese - Jan 19,1999

 $Log$
 Revision 1.3  2006/07/17 22:38:22  jmwille
 Adding missing copyright notices to a few files.

 Revision 1.2  2006/04/19 14:07:15  crdohrm
 Added and modified files so that building claps doesn't require Salinas include files

 Revision 1.15  2005/12/15 22:38:37  cwilson
 Remove cf specific ifdefs.  The class is now the same for all solvers.

 Revision 1.14  2004/12/14 16:28:07  gmreese
 minor change for 128 bit compatibility

 Revision 1.13  2004/12/13 19:42:36  gmreese
 corrected an error. MpcLocal may not include other files.

 Revision 1.12  2004/05/13 17:33:22  gmreese
 corrected comments for MPC Local. These specified (previously incorrectly)
 what happens when the MPC is on a subdomain boundary.

 Revision 1.11  2004/04/30 17:23:53  gmreese
 added some comments apropo the division of coefficients by the
 cardinality.

 Revision 1.10  2004/01/23 20:16:52  gmreese
 fixed the boo-boo in changing an interface file.

 Revision 1.9  2004/01/23 17:39:33  gmreese
 small H file changes from CU Feti integration.

 Revision 1.8  2003/02/20 21:54:59  mkbhard
 Changed all include preprocessor statements to use angle brackets < >, and
 not double quotes. NOTE: Please use angle brackets in the future to avoid
 problems that can occur when using SNtools.

 Revision 1.7  2002/10/03 16:08:55  gmreese
 added QA information to all files

 Revision 1.6  2002/03/01 14:47:38  gmreese
 partial port of serial salinas to 128 bit reals.

 Revision 1.5  2001/07/19 00:32:56  mkbhard
 *** empty log message ***

 Revision 1.4  1999/09/20 22:11:55  gmreese
 changes for local MPC definitions to require redundant specifications when
 the MPC lies on a boundary between subdomains.

 Revision 1.3  1999/06/30 21:13:24  gmreese
 changes to define MPCs in local framework for feti

 Revision 1.2  1999/06/14 16:29:31  gmreese
 First cut changes in MPC processing for FETI.

 Revision 1.1  1999/02/19 19:39:50  gmreese
 added a class for passing MPCs to feti.


*/
#endif
// $Date$
// $Revision$

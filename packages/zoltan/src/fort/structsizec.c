/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * Zoltan is distributed under the GNU Lesser General Public License 2.1.    * 
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

void c_structsize(char* a1, char* a2)
{
   printf("number of bytes in the structure is %d\n",a2-a1);
}

void c_structsize_(char* a1, char* a2)
{
   printf("number of bytes in the structure is %d\n",a2-a1);
}

void c_structsize__(char* a1, char* a2)
{
   printf("number of bytes in the structure is %d\n",a2-a1);
}

void C_STRUCTSIZE(char* a1, char* a2)
{
   printf("number of bytes in the structure is %d\n",a2-a1);
}

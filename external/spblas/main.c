main ()
{
   char          bool1, letter1; 
   int           numint1, numint2;
   float         numfloat1;
   double        numdoub1;
   short         numshor1;
   extern        void forts_ ();
   forts_(&bool1,&letter1,&numint1,&numint2,&numfloat1,
          &numdoub1,&numshor1, 1);
   printf(" %s %c %d %d %3.1f %.0f %d\n",
          bool1?"TRUE":"FALSE",letter1,numint1,
          numint2, numfloat1, numdoub1, numshor1);
}

       subroutine forts ( bool1, letter1, numint1, 
     & numint2, numfloat1, numdoub1, numshor1)
       logical*1           bool1
       character           letter1
       integer             numint1, numint2
       double precision    numdoub1
       real                numfloat1
       integer*2           numshor1
         
       bool1 = .true.
       letter1 = "v"
       numint1 = 11
       numint2 = -44
       numdoub1 = 902
       numfloat1 = 39.6
       numshor1 = 299
       return
       end

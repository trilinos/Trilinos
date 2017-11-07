function [] = reg2Comp(name,nProcs)

   for myRank=0:nProcs-1
     if exist(sprintf('myData_%d',myRank),'file'),
       fprintf('reg2Comp: myData_%d already exists\n',i);
       keyboard;
     end;

     fp = fopen(sprintf('myData_%d',myRank),'w');
     if fp == -1,
       fprintf('reg2Comp: cannot open myData_%d\n',myRank);
       keyboard;
     end;
     fprintf(fp,'Regional2Composite\n%s\n',name);
     fclose(fp);
   end




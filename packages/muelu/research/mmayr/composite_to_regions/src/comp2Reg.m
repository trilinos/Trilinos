function [] = comp2Reg(name,nProcs)

   for myRank=0:nProcs-1
     if exist(sprintf('myData_%d',myRank),'file'),
       fprintf('comp2Reg: myData_%d already exists\n',i);
       keyboard;
     end;

     fp = fopen(sprintf('myData_%d',myRank),'w');
     if fp == -1,
       fprintf('comp2Reg: cannot open myData_%d\n',myRank);
       keyboard;
     end;
     fprintf(fp,'Composite2Regional\n%s\n',name);
     fclose(fp);
   end




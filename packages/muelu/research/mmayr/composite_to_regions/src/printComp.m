function [] = printComp(name,nProcs)

   for myRank=0:nProcs-1
     if exist(sprintf('myData_%d',myRank),'file'),
       fprintf('printComp: myData_%d already exists\n',i);
       keyboard;
     end;

     fp = fopen(sprintf('myData_%d',myRank),'w');
     if fp == -1,
       fprintf('printComp: cannot open myData_%d\n',myRank);
       keyboard;
     end;
     fprintf(fp,'PrintComposite\n%s\n',name);
     fclose(fp);
   end




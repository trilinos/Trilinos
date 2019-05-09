function initCompVec(v,name,allMyNodes)

   for myRank=0:length(allMyNodes)-1,
     if exist(sprintf('myData_%d',myRank),'file'),
       fprintf('initCompVec: myData_%d already exists\n',i);
       keyboard;
     end;

     fp = fopen(sprintf('myData_%d',myRank),'w');
     if fp == -1,
       fprintf('initCompVec: cannot open myData_%d\n',myRank);
       keyboard;
     end;
     fprintf(fp,'InitCompositeVec\n%s\n',name);
     myGIDs = getCompositeIDs(allMyNodes{myRank+1},myRank); %myRank starts at 0

     for j=1:length(myGIDs)
       fprintf(fp,'%d\n',v(myGIDs(j)+1));  % myGIDs starts at 0
     end;
     fclose(fp);
   end

      

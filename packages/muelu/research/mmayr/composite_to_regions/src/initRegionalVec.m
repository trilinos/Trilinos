function initRegionalVec(v,name,allMyNodes,allMyRegions)

fprintf('inside initRegional\n');
   for myRank=0:length(allMyNodes)-1,
     if exist(sprintf('myData_%d',myRank),'file'),
       fprintf('initRegionalVec: myData_%d already exists\n',i);
       keyboard;
     end;

     myNodes  = allMyNodes{myRank+1};
     myRegions= allMyRegions{myRank+1};
     myRegions= myRegions.myRegions;

     for j=1:length(myRegions)
        fp = fopen(sprintf('myData_%d',myRank),'w');
        if fp == -1,
          fprintf('initRegionalVec: cannot open myData_%d\n',myRank);
          keyboard;
        end;
        fprintf(fp,'InitRegionalVec\n%s\n',name);
        curRegion = myRegions(j);
        fprintf(fp,'%d\n',curRegion);
        for k=1:length(myNodes)
          z = find(myNodes(k).gRegions == curRegion);
          if ~isempty(z) && (myNodes(k).procIds(z)==myRank),
            fprintf(fp,'%d\n',v(myNodes(k).ID+1));   % myNodes(k).ID starts at 0
          end;
        end;
        fclose(fp);

        while exist(sprintf('myData_%d',myRank),'file') ~= 0,
          pause(.5);
        end;
     end;
   end

      

function [] = commRegAssignments(allMyNodes,myRank,maxRegPerGID,A)

   myNodes = allMyNodes{myRank+1};
   if exist(sprintf('myData_%d',myRank),'file'),
      fprintf('commRegAssignments: myData_%d already exists\n',myRank);
      keyboard;
   end;
   fp = fopen(sprintf('myData_%d',myRank),'w');
   if fp == -1,
      fprintf('commRegAssignments: cannot open myData_%d\n',myRank);
      keyboard;
   end;
   myGIDs = getCompositeIDs(myNodes,myRank);
   [aaa,cols,ccc] = find(A(myGIDs+1,:));
   cols = unique(cols);
   nCols = length(cols);
   colPrinted = -ones(max(cols),1);
   colPrinted(cols) = 0;

fprintf('communicating %d\n',myRank);
   fprintf(fp,'CommRegAssignments\n');
if 1 == 0,
   for i=1:length(myNodes),
      fprintf(fp,'%d ',myNodes(i).ID);
      colPrinted(myNodes(i).ID+1) = 1;
      temp = -ones(maxRegPerGID,1);
      temp(1:length(myNodes(i).gRegions)) = myNodes(i).gRegions;
      fprintf(fp,'%d ',temp);
      fprintf(fp,'\n');
   end;
end

count = 1;
for i=1:length(myGIDs)
   fprintf(fp,'%d ',myGIDs(i));
   colPrinted(myGIDs(i)+1) = 1;
    while myNodes(count).ID ~= myGIDs(i), count = count+1; end;
    temp = -ones(maxRegPerGID,1);
    temp(1:length(myNodes(count).gRegions)) = myNodes(count).gRegions;
    fprintf(fp,'%d ',temp);
    fprintf(fp,'\n');
end;



   % find anything that should be in column map that hasn't already been printed
   % we will do this in an inefficient way by looping over all procs.
   notdone = find(colPrinted == 0);
   for ii=1:length(notdone),
     i = notdone(ii);
     flag = 0;
     for j=1:length(allMyNodes),
       if j~= myRank+1,
          for k=1:length(allMyNodes{j}),
             if allMyNodes{j}(k).ID == i-1,
                fprintf(fp,'%d ',i-1);
                temp = -ones(maxRegPerGID,1);
                temp(1:length(allMyNodes{j}(k).gRegions)) = allMyNodes{j}(k).gRegions;
                fprintf(fp,'%d ',temp);
                fprintf(fp,'\n');
                flag = 1;
                break;
             end;
          end;
       end 
       if flag == 1, break; end;
     end
   end 

   fclose(fp);


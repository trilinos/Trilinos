function [] = send(str,nProcs)

% we can skip the wait step by setting nProcs to -nProcs

if nProcs >= 0, waitForRmDataFiles(nProcs);
else            nProcs = -nProcs; end


for myRank=0:nProcs-1,
   fp=fopen(sprintf('myData_%d',myRank),'w');
   fprintf(fp,'%s\n',str);
   fclose(fp); 
end;

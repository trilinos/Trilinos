function [] = waitForRmDataFiles(nProcs)

   
   pause on;

   flag  = 1;
   while flag,
      flag = 0;
      for i=0:nProcs-1,
         if exist(sprintf('myData_%d',i),'file') ~= 0, flag = 1; break; end;
      end
      if flag, pause(3); end;
   end
         
   pause off;

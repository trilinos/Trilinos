%
% Low level file that takes the output from GrepTimings and produces
% matlab commands that generate pictures.
%
Nproc = length(Proc);

xvals = zeros(2*Nproc,1);
yvals = zeros(2*Nproc,1);
xvals(1:Nproc) = Proc(1:Nproc);
xvals(Nproc+1:Nproc*2) = Proc(Nproc:-1:1);
newsum = zeros(1,Nproc);


clf
hold on;
count = 0;
str = sprintf('AvMaxTotal%d',count);
while (  exist(str) ) 
   oldsum = newsum;
   str = sprintf('newsum = oldsum + %s;',str); eval(str);
   yvals(1:Nproc) = oldsum;
   yvals(Nproc+1:Nproc*2) = newsum(Nproc:-1:1);
   fill(xvals,yvals,[.2 .0 .0]);
   
   oldsum    = newsum;
   str = sprintf('NoCommsum = newsum - AvMaxComm%d;',count); eval(str);
   yvals(1:Nproc) = oldsum;
   yvals(Nproc+1:Nproc*2) = NoCommsum(Nproc:-1:1);
   fill(xvals,yvals,[.9 .9 .9]);

   count = count + 1;
   str = sprintf('AvMaxTotal%d',count);
end;
count = 0;
str = sprintf('RMax%d',count);
while (  exist(str) ) 
   oldsum = newsum;
   str = sprintf('newsum = oldsum + %s;',str); eval(str);
   str = sprintf('newsum = newsum + PMax%d;',count+1); eval(str);
   yvals(1:Nproc) = oldsum;
   yvals(Nproc+1:Nproc*2) = newsum(Nproc:-1:1);
   patch(xvals,yvals,[.9 .0 .0]);

   oldsum    = newsum;
   str = sprintf('NoCommsum = newsum - RMaxComm%d;',count); eval(str);
   str = sprintf('NoCommsum = NoCommsum - PMaxComm%d;',count+1); eval(str);
   yvals(1:Nproc) = oldsum;
   yvals(Nproc+1:Nproc*2) = NoCommsum(Nproc:-1:1);
   fill(xvals,yvals,[.9 .9 .9]);

   count = count + 1;
   str = sprintf('RMax%d',count);
end;


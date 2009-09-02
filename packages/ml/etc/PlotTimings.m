%
% Low level file that takes the output from GrepTimings and produces
% matlab commands that generate pictures. The basic idea is that we will
% draw a bunch of polygons. The x-axis of our plot will correspond to
% processors. The y-axis will corresponding to time. The lower part of
% a polygon is drawn by starting with small processors and going up to large 
% processors. The data on the lower part of the polygon corresponds to the
% top of the previous polygon. Note: If this is the first polygon, then
% the bottom of the polygon is just a bunch of zeros.  The top part of the 
% polygon is drawn by starting with large processors and going down.  The
% top part of the polygon corresponds to the current data.
% 
%  Note: Users may wish to change the first few variables below.

cutoff = .05;                   % determines data to be grouped into
                                % the other category. Basically, we sum
                                % the relevant cases and put a particular
                                % set of data into 'Other' if this field over 
                                % all runs is less than cutoff*sum
                                %
Cases = {'AvMax' 'RMax' 'PMax'};% Lists cases to be plotted
CompColor = lines;              % color scheme for computation times
CommColor = spring;             % color scheme for communication times
%
Nruns = length(Proc);
xvals = zeros(2*Nruns,1);
yvals = zeros(2*Nruns,1);
xvals(1:Nruns) = Proc(1:Nruns);
xvals(Nruns+1:Nruns*2) = Proc(Nruns:-1:1);
CurrData = zeros(1,Nruns);
TheSums  = zeros(1,Nruns);
Other    = zeros(1,Nruns);
LegendCount = 0;
LegendStrings = {};
ScaleData = Proc;
%ScaleData = ones(1,Nruns);     % MIGHT ONE TO UNCOMMENT THIS SO THAT
                                % OUTPUT IS NOT SCALED BY NUMBER OF PROCS

% Sum up all the fields so that we can use these sums in the cutoff scheme.
for jj=1:size(Cases,2)
   for count=0:10
      str = sprintf('%s%d',char(Cases(jj)),count);
      if ( exist(str) )
         PrevData = TheSums;
         % add data corresponding to the file given by 'str'
         % to the previous data.
         tmpstr = sprintf('TmpA = ScaleData .* %s;',str); eval(tmpstr);
         if (length(TmpA) < Nruns) 
            fprintf('Warning: Not all runs have all levels.\n');
            TmpB = zeros(1,Nruns);
            TmpB(1:length(TmpA)) = TmpA;
            TmpA = TmpB;
         end;
         TheSums = PrevData + TmpA;
      end;
   end;
end;

clf
hold on;
for jj=1:size(Cases,2)
   for count=0:10
      flag = 0;
      str = sprintf('%s%d',char(Cases(jj)),count);
      if ( exist(str) )
         

         PrevData = CurrData;
         % add data corresponding to the file given by 'str'
         % to the previous data and plot the polygon. This data
         % corresponds to a matrix vector product with communication.
         tmpstr = sprintf('TmpA = ScaleData .* %s;',str); eval(tmpstr);
         if (length(TmpA) < Nruns) 
            fprintf('Warning: Not all runs have all levels.\n');
            TmpB = zeros(1,Nruns);
            TmpB(1:length(TmpA)) = TmpA;
            TmpA = TmpB;
         end;
         if ( max( TmpA ./ TheSums) > cutoff)
            flag = 1;
            CurrData = PrevData + TmpA;
            yvals(1:Nruns) = PrevData;
            yvals(Nruns+1:Nruns*2) = CurrData(Nruns:-1:1);
            LegendCount = LegendCount + 1;
            if (Nruns == 1)
               t(1) = xvals(1)*.9; t(2) = xvals(2); t(3) = t(2); t(4)=t(1);
               q(1) = yvals(1); q(2) = q(1); q(3) = yvals(2); q(4) = q(3);
               fill(t,q, CompColor(LegendCount,:) );
            else
               fill(xvals,yvals, CompColor(LegendCount,:) );
            end
            xyz=sprintf('LegendStrings{LegendCount} = ''%s'';',str); eval(xyz);
         else
            Other = Other + TmpA;
         end;
   
         % draw a polygon on top of the previous one where the time
         % for communication is carved out???
         PrevData    = CurrData;
         str = sprintf('%sComm%d',char(Cases(jj)),count);
         tmpstr = sprintf('TmpA = ScaleData .* %s;',str); eval(tmpstr);
         if (length(TmpA) < Nruns) 
            fprintf('Warning: Not all runs have all levels.\n');
            TmpB = zeros(1,Nruns);
            TmpB(1:length(TmpA)) = TmpA;
            TmpA = TmpB;
         end;
         if ( flag == 1) 
            NoCommsum = CurrData - TmpA;
            yvals(1:Nruns) = PrevData;
            yvals(Nruns+1:Nruns*2) = NoCommsum(Nruns:-1:1);
            LegendCount = LegendCount + 1;
            if (Nruns == 1)
               t(1) = xvals(1)*.9; t(2) = xvals(2); t(3) = t(2); t(4)=t(1);
               q(1) = yvals(1); q(2) = q(1); q(3) = yvals(2); q(4) = q(3);
               fill(t,q, CommColor(LegendCount,:) );
             else
               fill(xvals,yvals, CommColor(65-LegendCount,:) );
             end;
             xyz=sprintf('LegendStrings{LegendCount} = ''%s'';',str); eval(xyz);
         end;
      end;
   end;
end;
if ( norm(Other) ~= 0.0) 
   CurrData = PrevData + Other;
   yvals(1:Nruns) = PrevData;
   yvals(Nruns+1:Nruns*2) = CurrData(Nruns:-1:1);
   LegendCount = LegendCount + 1;
   if (Nruns == 1)
      t(1) = xvals(1)*.9; t(2) = xvals(2); t(3) = t(2); t(4)=t(1);
      q(1) = yvals(1); q(2) = q(1); q(3) = yvals(2); q(4) = q(3);
      fill(t,q, CompColor(LegendCount,:) );
   else
      fill(xvals,yvals, CompColor(LegendCount,:) );
   end
   xyz=sprintf('LegendStrings{LegendCount} = ''Other'';',str); eval(xyz);
end;
legend(LegendStrings);
xlabel('processors','FontSize',14);
ylabel('time','FontSize',14);
title('Profile Times for Repeated Matvecs','FontSize',14);
print -deps apicture
%print -djpeg95 apicture

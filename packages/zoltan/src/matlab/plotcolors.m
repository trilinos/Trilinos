%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This matlab function returns a unique color/symbol pair for plotting.
% The total number of color/symbol pairs needed and the index of this
% color/symbol are taken as input arguments.
%
% Written by Michael Wolf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [colorV,symb] = plotcolors(cindex,numcolors)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% For 4 or fewer colors
  if numcolors <= 4
    switch cindex
      case 1
        colorV = [1 0 0];
      case 2
        colorV = [0 1 0];
      case 3
        colorV = [0 0 1];
      otherwise
        colorV = [0 1 1];
    end
    symb = '.';
  %% For colors
  else
    %      cmap = lines(numcolors);
    %      cmap = hsv(numcolors);

    %% ten distinct colors, add more later
    diffColors=8;
    cmap(1,:) = [ 1 0 0];       % red
    cmap(2,:) = [ 0 1 0];       % green
    cmap(3,:) = [ 0 0 1];       % blue
    cmap(4,:) = [ 0 1 1];       % cyan
    cmap(5,:) = [ 1 0 1];       % magenta
    cmap(6,:) = [ 0 0 0];       % black
    cmap(7,:) = [ 1 0.67 0];    % orange
    cmap(8,:) = [ 0.6 0.6 0.6]; % gray
%    cmap(9,:) = [0.75 0.75 0]; % dark yellow
%    cmap(10,:) = [ 0 0.4 0];     % dark green

    colorV = cmap(mod(cindex-1,diffColors)+1,:);


     if floor((cindex-1)/diffColors) == 0
       symb = '.';
     else
       symb = 'o';
     end
  

   

  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end

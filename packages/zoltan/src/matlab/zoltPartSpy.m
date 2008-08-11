%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This matlab function produces spy-like plots in which the nonzeros are 
% colored based on the partition that owns them.
%
% The function takes 3 arguments: the filename of a matrix in matrix market 
% format, the name of the zoltan output file containg partition info, and
% the type of partitioning ('1DRow', '1DColumn', '1.5DRow', '1.5DColumn'). 
%
% Written by Michael Wolf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partSpyZoltan(matFilename,partFilename,partType)
  if nargin ~= 3
     error('Wrong number of input arguments. \n %s', ...
           'Usage: partSpyZoltan mtxFile partitionFile PartitionType')
  end

  
  % Parse zdrive output file, Erik Boman's code
  fp = fopen(partFilename, 'r');
  if (~fp)
    error('Could not open partition file\n');
  end
  
  % Skip all lines before 'GID'
  word = 'abc';
  while (~strncmp(word, 'GID', 3))
    % Read only first word, ignore rest of line
    word = fscanf(fp,'%s',1);
    fscanf(fp,'%*[^\n]',1);
  end
  
  % Read the partition numbers; file has 4 fields per line
  [P, num] = fscanf(fp, '%d', [4,inf]);
  % First two rows in P (columns in output file) give the partition vector
  part = zeros(size(P(1,:)));
  part(P(1,:)) = P(2,:);

  
  
  
   if strcmp(partType, '1DColumn')
     partSpyZolt1DC(matFilename,part);
   elseif strcmp(partType, '1DRow')
     partSpyZolt1DR(matFilename,part);
   elseif strcmp(partType, '1.5DColumn')
     partSpyZolt1_5DC(matFilename,part);
   elseif strcmp(partType, '1.5DRow')
     partSpyZolt1_5DR(matFilename,part);
   else
     error('Unsupported partitioning scheme ''%s''',partType);
   end
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for zoltan with stripped partition information
% 1d column partitioning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partSpyZolt1DC(matFilename,part)

   A = mmread(matFilename);

   numParts=max( part ) +1;
     
   numParts

   for i=1:numParts
     [colors(i,:),symbols(i)]=plotcolors(i,numParts); 
     x{i}=[];
     y{i}=[];
   end

   [i,j,val] = find(A);
   n=length(i);

   for cnt=1:n
      x{part(j(cnt))+1}(end+1) = j(cnt);
      y{part(j(cnt))+1}(end+1) = i(cnt);
   end

   figure;
   hold on;
   for cnt=1:numParts
       plot(x{cnt},y{cnt}, symbols(cnt), 'Color', colors(cnt,:), ...
           'MarkerSize', 6);
   end
      
   view([0 0 -1])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for zoltan with stripped partition information
% 1d row partitioning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partSpyZolt1DR(matFilename,part)

   A = mmread(matFilename);

   numParts=max( part ) +1;
   
   numParts

   for i=1:numParts
     [colors(i,:),symbols(i)]=plotcolors(i,numParts);     
     x{i}=[];
     y{i}=[];
   end

   [i,j,val] = find(A);
   n=length(i);


   for cnt=1:n
      x{part(i(cnt))+1}(end+1) = j(cnt);
      y{part(i(cnt))+1}(end+1) = i(cnt);     
   end

   figure;
   hold on;
   for cnt=1:numParts
       plot(x{cnt},y{cnt}, symbols(cnt), 'Color', colors(cnt,:), ...
           'MarkerSize', 6);
   end  
   
   view([0 0 -1])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for zoltan with stripped partition information
% 1.5d col partitioning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partSpyZolt1_5DC(matFilename,part)

   A = mmread(matFilename);

   numParts=max( part ) +1;
   
   numParts

   for i=1:numParts
     [colors(i,:),symbols(i)]=plotcolors(i,numParts);   
     x{i}=[];
     y{i}=[];
   end

   [i,j,val] = find(A);
   n=length(i);

   for cnt=1:n
      if i(cnt) >= j(cnt)
        x{part(j(cnt))+1}(end+1) = j(cnt);
        y{part(j(cnt))+1}(end+1) = i(cnt);     
      else
        x{part(i(cnt))+1}(end+1) = j(cnt);
        y{part(i(cnt))+1}(end+1) = i(cnt);     
      end
   end
   
   figure;
   hold on;
   for cnt=1:numParts
       plot(x{cnt},y{cnt}, symbols(cnt), 'Color', colors(cnt,:), ...
           'MarkerSize', 6);
   end  
   
   view([0 0 -1])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for zoltan with stripped partition information
% 1.5d row partitioning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partSpyZolt1_5DR(matFilename,part)

   A = mmread(matFilename);

   numParts=max( part ) +1;
   
   numParts

   for i=1:numParts
     [colors(i,:),symbols(i)]=plotcolors(i,numParts);  
     x{i}=[];
     y{i}=[];
   end

   [i,j,val] = find(A);
   n=length(i);

   for cnt=1:n
      if i(cnt) >= j(cnt)
        x{part(i(cnt))+1}(end+1) = j(cnt);
        y{part(i(cnt))+1}(end+1) = i(cnt);     
      else
        x{part(j(cnt))+1}(end+1) = j(cnt);
        y{part(j(cnt))+1}(end+1) = i(cnt);     
      end
   end
   
   figure;
   hold on;
   for cnt=1:numParts
       plot(x{cnt},y{cnt}, symbols(cnt), 'Color', colors(cnt,:), ...
           'MarkerSize', 6);
   end  
   
   view([0 0 -1])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute coarse level operator for caseThree.
%
% There are 5 regions.
%

%% Define region-wise quantities
% prolongator for region 1
regP1 = [1 0;
         1 0;
         0 1;
         0 1];
       
% prolongator for region 2
regP2 = [1 0 0;
         1 0 0;
         0 1 0;
         0 1 0;
         0 1 0;
         0 0 1;
         0 0 1];

% prolongator for region 3       
regP3 = [1 0 0 0;
         1 0 0 0;
         0 1 0 0;
         0 1 0 0;
         0 1 0 0;
         0 0 1 0;
         0 0 1 0;
         0 0 1 0;
         0 0 0 1;
         0 0 0 1];
       
% prolongator for region 4
regP4 = [1 0 0;
         1 0 0;
         0 1 0;
         0 1 0;
         0 1 0;
         0 0 1;
         0 0 1];       
       
% prolongator for region 5
regP5 = [1 0;
         1 0;
         0 1;
         0 1];

% matrix for region 1
regA1 = full(oneDimensionalLaplace(4));

% matrix for region 2
regA2 = full(oneDimensionalLaplace(7));
regA2(1,2) = -1;

% matrix for region 3
regA3 = full(oneDimensionalLaplace(10));
regA3(1,2) = -1;

% matrix for region 4
regA4 = full(oneDimensionalLaplace(7));
regA4(1,2) = -1;

% matrix for region 5
regA5 = full(oneDimensionalLaplace(4));
regA5(1,2) = -1;

%% Build global, but regional matrices
regP = [regP1 zeros(size(regP1,1),size(regP2,2)) zeros(size(regP1,1),size(regP3,2)) zeros(size(regP1,1),size(regP4,2)) zeros(size(regP1,1),size(regP5,2));
        zeros(size(regP2,1),size(regP1,2)) regP2 zeros(size(regP2,1),size(regP3,2)) zeros(size(regP2,1),size(regP4,2)) zeros(size(regP2,1),size(regP5,2));
        zeros(size(regP3,1),size(regP1,2)) zeros(size(regP3,1),size(regP2,2)) regP3 zeros(size(regP3,1),size(regP4,2)) zeros(size(regP3,1),size(regP5,2));
        zeros(size(regP4,1),size(regP1,2)) zeros(size(regP4,1),size(regP2,2)) zeros(size(regP4,1),size(regP3,2)) regP4 zeros(size(regP4,1),size(regP5,2));
        zeros(size(regP5,1),size(regP1,2)) zeros(size(regP5,1),size(regP2,2)) zeros(size(regP5,1),size(regP3,2)) zeros(size(regP5,1),size(regP4,2)) regP5];
       
regA = [regA1 zeros(size(regA1,1),size(regA2,2)) zeros(size(regA1,1),size(regA3,2)) zeros(size(regA1,1),size(regA4,2)) zeros(size(regA1,1),size(regA5,2));
        zeros(size(regA2,1),size(regA1,2)) regA2 zeros(size(regA2,1),size(regA3,2)) zeros(size(regA2,1),size(regA4,2)) zeros(size(regA2,1),size(regA5,2));
        zeros(size(regA3,1),size(regA1,2)) zeros(size(regA3,1),size(regA2,2)) regA3 zeros(size(regA3,1),size(regA4,2)) zeros(size(regA3,1),size(regA5,2));
        zeros(size(regA4,1),size(regA1,2)) zeros(size(regA4,1),size(regA2,2)) zeros(size(regA4,1),size(regA3,2)) regA4 zeros(size(regA4,1),size(regA5,2));
        zeros(size(regA5,1),size(regA1,2)) zeros(size(regA5,1),size(regA2,2)) zeros(size(regA5,1),size(regA3,2)) zeros(size(regA5,1),size(regA4,2)) regA5];
      
regAP = regA * regP;
regRAP = regP' * regAP;

 
% %% Construct fine level residual
% 
% fRegRes1 = [0
%             1.3300
%            -3.3600
%            -0.6000
%            -0.7200
%             1.6650
%             0.0350];
% fRegRes2 = [0.0350
%            -0.0250
%             3.0050
%             2.0000
%             2.0000
%             2.0000
%             2.0000
%             2.0000
%           -25.1350
%            -3.6250];
% fRegRes3 = [-3.6250
%            -4.3450
%            20.0900
%             0.0350
%            -1.3650
%             1.6850
%            -0.6800
%             0.6800
%            -0.6793
%             0.3403];
% 
% % scale interface values
% fRegRes1(end) = fRegRes1(end) / 2;
% fRegRes2(1) = fRegRes2(1) / 2;
% fRegRes2(end) = fRegRes2(end) / 2;
% fRegRes3(1) = fRegRes3(1) / 2;
%           
% fRegRes = [fRegRes1; fRegRes2; fRegRes3];
% 
% %% Perform restriction
% 
% cRegB = regP' * fRegRes;



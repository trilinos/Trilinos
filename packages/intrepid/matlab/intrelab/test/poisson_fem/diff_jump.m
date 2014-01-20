function k = diff_jump(usr_par)

%%%%%%%%%%% get the material properties

% x and y coordinates of the circular region with different material properties
center    = [-0.25 -0.25];

% radius of the circular region with different material properties
radius    = 0.3;

% difussion constant value for each material region
mat_prop  = [6e0 2e0];

%%%%%%%%%%% get diffusion spatial distribution
layer1 = ( ( usr_par.mesh.p(:,1) - center(1) ).^2 + ...
    ( usr_par.mesh.p(:,2) - center(2) ).^2 ) <= radius^2;
layer2 = ~layer1;
k = mat_prop(1)*layer1 + mat_prop(2)*layer2;

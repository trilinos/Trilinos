function [Ux, Uy, P, adj, nodes] = parseData(name)

adj   = load('cell_to_node_quad.txt') + 1; %% load node adjacency, increment by 1 for 1-based indexing
nodes = load('nodes.txt');                 %% load node coordinates
N     = 3*length(nodes);                   %% determine number of nodes

%% Read in state
state  = importdata([name,'.txt']);
[M,nt] = size(state);
map    = importdata(['map_',name,'.txt']);
map    = map(1:2:end)+1;
[tmp, perm] = sort(map);
state  = state(perm,:); %% we need to permute the state according to parallel maps

Ux = state(1:3:N,:);
Uy = state(2:3:N,:);
P  = state(3:3:N,:);

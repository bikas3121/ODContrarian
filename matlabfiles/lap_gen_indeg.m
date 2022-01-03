function [L] = lap_gen_indeg (N,c,s,t)
% s = [1 1 1 1  2  3 3  4  5 5  6 6  7 7  7  8 8  9 9  10 10  11 11  11 12 12 12 12  13  14 14  15 15];
% t=  [6 3 5 7  4  2 1  1  1 6  5 11 5 10 12 6 9  2 6  13 7   8  14   6  6  7  8 15   9  12 6   12  6];
% N= 15;
%c =[1 -1 1 1 1 -1 1 1 1 -1 1 1 -1 1 1  ]; 
%% Adjacency Matrix
adj = zeros(N,N);
for i = 1:length(s)
   adj(t(i),s(i)) =1;
end
%% Degree Matrix
deg = zeros(N,N);
dg = 0;
for i = 1:N
    for j  = 1:N
       dg = dg+ adj(i,j);
    end
    deg(i,i) = dg;
    dg = 0;
end

%% The Laplacian is given by the following
%In our case we have the conformist and contrarian nodes so, we identify
%%them with the matrix c. That is, if c(i) = 1, then i node is conformist 
% and if c(i) = -1 then the node i is contrarian. 
%c = [1 -1 1 1 -1];
adj1 = zeros(N,N);
for i = 1:length(s)
   adj1(t(i),s(i)) =c(t(i))*c(s(i))*1;
end

% The Laplacain is the following
L = deg-adj1;
%end

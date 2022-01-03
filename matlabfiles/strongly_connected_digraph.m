%% Generate the permutation of the nodes.
perm_edges = []; 
NumberOfNodes = 10;
for i =1:NumberOfNodes
    for j = 1:NumberOfNodes
        perm_edges = [perm_edges;i,j;];
    end
end
 perm_edges(diff(perm_edges,[],2)==0, :) = []; %remove self-loops
 

 %% Generate the edges
 len_perm = length(perm_edges); 
 % To make some node popular we chose the nodes in random and make the
 % number of edges to and from these nodes more than the average. 
 num_pop_node = 3;
 pop_node = sort(randperm(NumberOfNodes,num_pop_node),'ascend'); % generate random popular nodes and sort 
 
 
 
%% This generates the edges for the graph while considering the popular edges. 
percent_of_edges = 0.2*NumberOfNodes;
pos_edges = randi([1,percent_of_edges],NumberOfNodes,1)';      %number of possible edges for each node.
num_edge_pop_nodes = NumberOfNodes*0.3  ;          % this is the number of edges assigned to the popular node.
% This assigns the number of edges for the popular nodes. 
for i = 1: NumberOfNodes
    for j = 1:length(pop_node)
             if i == pop_node(j)
                 pos_edges(i) = num_edge_pop_nodes;
             end
    end
end

% pos_edges gives the number of possible edges for each node. 

%% Now Generate the graph based on the above construction. 
% The number of maximum possible node for each node is given by the array
% pos_edges.
set_edges= [];
for i = 1:length(pos_edges)
   a = pos_edges(i);
   p = (NumberOfNodes-1)*i;
   if i == 1
        pos = randperm(NumberOfNodes-1,a);
        for j = 1:a
            edg = perm_edges(pos(j),:);
            set_edges = [set_edges;edg];
        end      
   else
        pos = randperm(NumberOfNodes-1,a)+((NumberOfNodes-1)*(i-1));
        for j = 1:a
            edg = perm_edges(pos(j),:);
            set_edges = [set_edges;edg];
        end  
   end
end

%% Plotting of the Graph
 s = set_edges(:,1);
 t = set_edges(:,2);
G = digraph(s,t,'omitselfloops');
% figure
% plot(G,'Layout','force')
bins = conncomp(G);
n_bin = sum(bins);

%% Make the graph strongly connected. 
N = NumberOfNodes;
for i = 1:N
    
    




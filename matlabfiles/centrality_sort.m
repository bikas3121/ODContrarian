%% centrality sort. 
function centrality_sort (v1)
%v1 = [0.0749    0.0543    0.2235   -0.0956    0.1111   -0.0517    0.1085   -0.1602    0.0568    0.0633];
%I = [3     8     7     5     9     6     4    10     2     1];
[centr_sort,ind] = sort(abs(v1),'descend');
figure
 bar(centr_sort)
 set(gca,'xticklabel',ind)
 xlabel('Agent Index')
 ylabel ('Centrality')
end
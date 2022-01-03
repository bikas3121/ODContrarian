%  function GraphPlot(N,s,t,I,v1,c)
N =10;
% load('c_10.mat','c')
c = [1     1    1     1    -1     1    -1     1     1     -1];
s = [1     1     2     2     3     3     4     4     5     5     6     6     7     7     8     8     9  9    10    10];
 t = [2     4     9     6    10     8     3    10     2     9     5     4     6     9     7     6     4  1     1     6];
 %weights = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];
 %I = [3     8     7     5     9     6     4    10     2     1];
 %v1 = [ 0.0749    0.0543    0.2235    0.0956   -0.1111    0.0517   -0.1085    0.1602    0.0568   -0.0633];
G = digraph(s,t);
p=plot(G);
% C = abs(v1)*100;
C = 4;
p.MarkerSize=2*C;
 s1 = [];
 t1 = [];
for i = 1:N
    if c(i) == 1
        s1 = [s1 i];
    else
        t1 = [t1 i];
    end
end
labelnode(p,[1 2 ],{'1' '2'})
highlight (p,s1,'NodeColor','r');
highlight (p,t1,'NodeColor','c');
highlight(p,G,'LineWidth',1.2)

% cont_conf = [];
% for i = 1:length(s)
%     a1 = s(i);
%     b1  = c(a1);
%     a2  = t(i);
%     b2 = c(a2);
%         if b1 == b2
%             cont_conf = [cont_conf; "+"];
%         else
%             cont_conf = [cont_conf; "-"];
%         end
% end
labelText = {'ABC' 'DEF' 'GHI'};
% labeledge(p,s,t,cont_conf')        
% end



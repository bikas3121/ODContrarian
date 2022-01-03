function [gam, gamma, I] = sorting_function(N,c,d,x_0,v)
%x_0 = [0.1667   -1.0000    0.4444    0.7222    1.0000   -0.1667];
gam = [];
% gan = 0;
for i = 1:N
    gan = abs(v(i))*abs((c(i)*d)-x_0(i));
    gam = [gam; gan];
end
   [gamma,I] = sort(gam,'descend') ;
end
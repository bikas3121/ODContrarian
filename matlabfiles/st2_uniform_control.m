function [J_unif] = st2_uniform_control(N,u_b,B,I1,P,c,x_0,d,N_conf,N_cont,v1)
%% Optimal Distribution of the Control
cont_u = [];
for i = 1:N
    B = i*u_b;
    u1 = B/N;
    cont_u1 = [];
    for j = 1:N
        cont_u1(j) = c(j)*u1;
    end
       cont_u = [cont_u;cont_u1];
end

%%
cas_c = zeros(N,N);
for i = 1:N
        cas_c(i,I1(i)) = cont_u(I1(i));
end

idx_nonzerolines = sum(abs(cas_c),2)>0 ;
cas_c = cas_c(idx_nonzerolines,:) ;

cas_c1 = zeros(length(cas_c(:,1)),length(cas_c(1,:)));
for i = 1:length(cas_c(:,1))-1
    cas_c1(1,:) = cas_c(1,:);
    a1 = cas_c1(i,:);
    a2 = cas_c(i+1,:);
    cas_c1(i+1,:) = a1+a2;
end

cas_c1 = [zeros(1,N); cas_c1];

%% Impulse 
len = length(cas_c1(:,1));
x_t = zeros(len,N);
for i = 1:len
    for j = 1:N
        cont_u1 = cas_c1(i,:);
       x_t(i,j)= (1-abs(cont_u1(j)))*x_0(j) + cont_u1(j)*d; 
    end
end
x_t_optimal = x_t(len,:);

%% Cost plot
J_cas = zeros(1,len);
b1= (N_conf-N_cont)/N;
for i = 1:len
 J_cas(i) = abs(b1*(v1*x_t(i,:)')-d);
end
J_unif = J_cas;
end
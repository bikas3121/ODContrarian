load('data_examples/s_300_4.mat','s')
load('data_examples/t_300_4.mat','t')
load('data/c_300.mat','c')
%c= [1 1 -1 1 -1 1 -1 1 1 1];
%load('data/x_0_10.mat','x_0')
N = 300 ;%number of nodes
d = 1; %desired opinion
N_conf = 0; %Number of conforming agents.
N_cont = 0; %Number of contrasting agents.
for i = 1:length(c)
    if c(i) == 1
        N_conf = N_conf +1;
    else
        N_cont = N_cont +1;
    end
end
u_b = 0.7; %upperbound on the budget

B = 210; %Total budget available
P = floor(B/u_b);

 %Initial opinion of the agents
%x_0 = [0.9    0.7    0.5    0.3   -0.3    0.1    -0.4    0.1   -0.5   0.6];
x_0 = zeros(1,N);
for i = 1:N
    x_0(i) = -1 + (2*i)/N;
end

%x_0 =  -1 + (1+1)*rand(1,N);
%% Laplacian of the in-degree graph 
% [L_in] = lap_gen_indeg(N,c,s,t);
% [u1, v1] = normalize_eigenvector(L_in);
% [gam,gamma1, I1] = sorting_function(N,c,d,x_0,v1);

[L_in] = lap_gen_indeg(N,c,s,t);
[u1, v1] = normalize_eigenvector(L_in);
[gam, gamma1, I1] = sorting_function(N,c,d,x_0,v1);


%% Graph Construction
%  GraphPlot(N,s,t,I1,v1,c)
% centrality_sort (v1)
% gamma_plot (gamma1,I1)


%% Optimal Control- Strategy 1
u_ran = [];
for i = 1:P
    u_ran = [u_ran u_b];
end
l = (B-sum(u_ran))/N;
if l >= 0.001
    u_ran = [u_ran (B-sum(u_ran))];
end
for i = 1:length(u_ran)
    if u_ran(i) == 0
        u_ran(i) = [];
    end
end

cont_u = zeros(1,N);
for i = 1:length(u_ran+1)
    cont_u(I1(i)) = u_ran(i);
end
cont_u = cont_u.*c;


cas_c = zeros(N,N);
for i = 1:N
    cas_c (i,I1(i)) = cont_u(I1(i));
end

idx_nonzerolines = sum(abs(cas_c),2)>0 ;
cas_c = cas_c(idx_nonzerolines,:) ;

%% cascade control
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
x_t1 = zeros(len,N);
for i = 1:len
    for j = 1:N
        cont_u1 = cas_c1(i,:);
       x_t1(i,j)= (1-abs(cont_u1(j)))*x_0(j) + cont_u1(j)*d; 
    end
end
% x_t_final = x_t1(len,:);
 ln1 = length(x_t1(:,1));
%x_t_in = x_t1(ln1,:);

%% Cost Calculate
b1= (N_conf-N_cont)/N;
J_initial = abs(b1*(v1*x_0')/N-d);


%% Cost plot
J_optimal = zeros(1,ln1);
for i = 1:len
 J_optimal(i) = abs(b1*(v1*x_t1(i,:)')-d);
end
%t1 = linspace(1,len,len);
% figure
% plot(t1,J_optimal)
J_optimal_value = J_optimal(len);
J_optimal(len) = [];


%% Uniformly distributed control- Strategy 2
cas_c2 = [];

for i = 1:N
    B1 = i*0.7;
    u1 = B1/N;
    cont_u2 = ones(1,N);
    cont_u2 =c.*u1.*cont_u2;
    cas_c2 = [cas_c2;cont_u2];
end
 cas_c2 = [zeros(1,N); cas_c2];   


%% Impulse 
len = length(cas_c2(:,1));
x_t1 = zeros(len,N);
for i = 1:len
    for j = 1:N
        cont_u2 = cas_c2(i,:);
       x_t2(i,j)= (1-abs(cont_u2(j)))*x_0(j) + cont_u2(j)*d; 
    end
end
% x_t_final = x_t1(len,:);
 ln2 = length(x_t2(:,1));
%x_t_in = x_t2(ln1,:);

%% Cost Calculate
b1= (N_conf-N_cont)/N;
J_initial = abs(b1*(v1*x_0')/N-d);


%% Cost plot
J_unif = zeros(1,ln2);
for i = 1:len
 J_unif(i) = abs(b1*(v1*x_t2(i,:)')-d);
end
%t1 = linspace(1,len,len);
% figure
% plot(t1,J_unif)
J_uinf_value = J_unif(len);
J_unif(len) = [];
%% Positive Control- Strategy 3
cas_c3 = [];

for i = 1:N
    B2 = i*0.7;
    u2 = B2/N;
    cont_u3 = ones(1,N);
    cont_u3 =u2.*cont_u3;
    cas_c3 = [cas_c3;cont_u3];
end
    cas_c3 = [zeros(1,N); cas_c3];

%% Impulse 
len = length(cas_c3(:,1));
x_t = zeros(len,N);
for i = 1:len
    for j = 1:N
        cont_u3 = cas_c3(i,:);
       x_t3(i,j)= (1-abs(cont_u3(j)))*x_0(j) + cont_u3(j)*d; 
    end
end

% x_t_final = x_t1(len,:);
 ln3 = length(x_t3(:,1));
%x_t_in = x_t3(ln1,:);

%% Cost Calculate
b1= (N_conf-N_cont)/N;
J_initial = abs(b1*(v1*x_0')/N-d);


%% Cost plot
J_positive = zeros(1,ln3);
for i = 1:len
 J_positive(i) = abs(b1*(v1*x_t3(i,:)')-d);
end
t1 = linspace(1,B,N);
% figure
% plot(t1,J_positive);
J_positive_value = J_positive(len);
J_positive(len) = [];
figure
plot(t1,J_optimal,t1,J_unif, t1,J_positive);
legend('Strategy 1 (Optimal allocation)','Strategy 2 (Uniform Allocation)','Strategy 3 (Positive Uniform Allocation)')
xlabel('Budget')
ylabel('Cost')
grid on 






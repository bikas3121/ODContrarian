load('data/s_10.mat','s')
load('data/t_10.mat','t')
load('data/c_10.mat','c')
%c= [-1 1 -1 1 -1 1 -1 1 -1 1];

N = 10 ;%number of nodes
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

B = 4; %Total budget available
P = floor(B/u_b);
% 
%% Initial Opinion of Agents
%load('data/x_0_10.mat','x_0')
%x_0 =  -1 + (1+1)*rand(1,N);
%x_0 =[ 0.1181    0.7082   -0.3042   -0.1079   -0.8915   -0.6458    0.3256   -0.3383    0.7970   -0.7637];
x_0 = [0.8    0.7    0.5    0.3   -0.3    0    -0.4   0.2   -0.5   0.6];
%x_0 = [1  0.8  0.6  0.4  -0.2  0  -0.4  0.2  -0.6    1];

% x_0 = zeros(1,N);
% for i = 1:N
%     x_0(i) = -1 + (2*i)/N;
% end
%% Laplacian of the in-degree graph 
[L_in] = lap_gen_indeg(N,c,s,t);
[u1, v1] = normalize_eigenvector(L_in);
[gam, gamma1, I1] = sorting_function(N,c,d,x_0,v1);

%% Graph Construction
 figure
 GraphPlot(N,s,t,I1,v1,c)
centrality_sort (v1)
% gamma_plot (gamma1,I1)

%% Allocation of the Budget
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

%% Uniform Budget Allocation
%cont_u = c.*((B/N)*ones(1,N));

%% Positive Budget Allocation
%cont_u = ((B/N)*ones(1,N));

%% Budget Input
x_t_in = zeros(1,N);
for i = 1:N
    x_t_in(i) = (1-abs(cont_u(i)))*x_0(i) + cont_u(i)*d; 
end
b1= (N_conf-N_cont)/N;
J_initial = abs(b1*(v1*x_0')-d);
J_optimal1 = abs(b1*(v1*x_t_in')-d); 
%% Inital Opinion Plot
x_0_1 = [];
for i = 1:10+1
    x_0_1(i,:) = x_0;
end
x_0_1 = [x_0_1;x_t_in];

% figure
% plot(t3,x_0_1)

 
% %% Agent Dynamics
[t1,x] = ode45(@(t,x) -L_in*x, [0 5], x_t_in);
x = [x_0_1;x];
x1 = x(length(x(:,1)),:);
J_optimal2 = abs(sum(x1)/N-d);
t3 = [];
for i = 1:N
    a = 1/10;
    t3(i) = -a*i;
end
t3 = flip(t3);
t3 = [t3 0 0];
t2 = [t3'; t1];
x_conf = [];
x_cont = [];
for i = 1:N
    if c(i) ==1
        x_conf = [x_conf x(:,i)];
    else 
        x_cont = [x_cont x(:,i)];
    end
end
figure

p1 = plot(t2,x_conf,'-o','MarkerIndices',11,'Color','r','linewidth',1.2);
hold on
p2 = plot(t2,x_cont,'-o','MarkerIndices',11,'Color','c','linewidth',1.2);
xlabel ('Time')
ylabel('Opinions')




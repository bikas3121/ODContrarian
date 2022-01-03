function gamma_plot (gamma1,I)
%gamma1 = [0.1037    0.0710    0.0668    0.0539    0.0404    0.0334    0.0144    0.0067         0         0];
%I = [1     2     6     3     8     5     9     7     4    10];
figure
bar(gamma1)
set(gca,'xticklabel',I)
xlabel('Agent Index')
ylabel('Gamma')
end
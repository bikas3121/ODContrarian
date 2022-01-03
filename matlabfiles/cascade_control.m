%% cascade control
function [cas_c1] = cascade_control(cas_c)
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
end
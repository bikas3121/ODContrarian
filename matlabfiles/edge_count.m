load('data/s_300.mat','s')
load('data/t_300.mat','t')
N =300;
ab = [];
for i = 1 :N
    p =0;
    for j = 1:length(s)
        if i == s(j)
            p=p+1;
        end
    end
    ab = [ab; p];
end
ind = [];
edg = [];
for i = 1:length(ab)
    if ab(i) >= 30
        ind = [ind; i];
        edg = [edg; ab(i)];
    end
end

len = length(ind);
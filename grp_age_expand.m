function contact_age = age_expander(mixed, ind)
mixed(isnan(mixed)) = 0;
lower = 1;
upper = 7;
age_bracket_mean = [2, 11, 21, 35, 55, 73, 91, 2, 11, 21, 35, 55, 73, 91, 2, 11, 21, 35, 55, 73, 91, 2, 11, 21, 35, 55, 73, 91, 2, 11, 21, 35, 55, 73, 91];

mixed_dim = size(mixed);
comb_age = cell(1,mixed_dim(1));
for i = 1:mixed_dim(1)
    if sum(mixed(i,:) ~= 0) == 0
        comb_age(1,i) = ind(:,i);
    else
        for j = 1:mixed_dim(2)
            ones(1,mixed(i,j))*age_bracket_mean(j);
        end
    end
end
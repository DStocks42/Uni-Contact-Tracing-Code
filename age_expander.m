function comb_age = age_expander(mixed, ind)
mixed(isnan(mixed)) = 0;
age_bracket_mean = [2, 11, 21, 35, 55, 73, 91, 2, 11, 21, 35, 55, 73, 91, 2, 11, 21, 35, 55, 73, 91, 2, 11, 21, 35, 55, 73, 91, 2, 11, 21, 35, 55, 73, 91];

mixed_dim = size(mixed);
comb_age = cell(1,mixed_dim(1));
for i = 1:mixed_dim(1)
    if sum(mixed(i,:) ~= 0) == 0
        comb_age{1,i} = ind(i,~isnan(ind(i,:)));
    else
        grp_ages = [];
        for j = 1:mixed_dim(2)
            grp_ages = [grp_ages, ones(1,mixed(i,j))*age_bracket_mean(j)];
        end
%         disp(size(ind(i,:)))
%         disp(size(grp_ages))
        comb_age{1,i} = [ind(i,~isnan(ind(i,:))),grp_ages];
    end
end
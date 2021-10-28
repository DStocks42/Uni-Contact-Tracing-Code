function grp_size = grp_size_builder(mix)
mix(isnan(mix)) = 0;
dim_mix = size(mix);
grp_size = zeros(1,dim_mix(1)*5);
for i = 1:dim_mix(1)
    lower = 1;
    upper = 7;
    c = 0;
    while upper <= 35
      grp_size(i+c) = sum(mix(i,lower:upper));
      lower = lower + 7;
      upper = upper + 7;
      c = c+1;
    end
end
end
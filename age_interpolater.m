function mix = age_interpolater(mix)
dim_mix = size(mix);
for i = 1:dim_mix(1)
    reported = ~isnan(mix(i,:));
    if sum(reported) == 0
        continue
    else
        for j = 1:sum(reported)
            if mix(i,j) == 1
                mix(i,j) = randi([0,4]);
            elseif mix(i,j) == 2
                mix(i,j) = randi([5,17]);
            elseif mix(i,j) == 3
                mix(i,j) = randi([18,24]);
            elseif mix(i,j) == 4
                mix(i,j) = randi([25,44]);
            elseif mix(i,j) == 5
                mix(i,j) = randi([45,64]);
            elseif mix(i,j) == 6
                mix(i,j) = randi([65,80]);
            elseif mix(i,j) == 7
                mix(i,j) = randi([81,100]);
            end
        end
    end
end

end
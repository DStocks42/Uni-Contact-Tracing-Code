function out_data = unique_data(org_data)
    org_data_dim = size(org_data);
    out_data = zeros(length(unique(org_data(:,1))), org_data_dim(2));
    seen = [];
    c = 1;
    for i =1:length(org_data(:,1))
        if sum(ismember(seen, org_data(i,1))) > 0
            continue
        else
            out_data(c,:) = org_data(i,:);
            seen = [seen, org_data(i,1)];
            c = c+1;
        end
    end
end
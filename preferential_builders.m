function network_mat = preferential_builders(network_mat, desired_clustering)
    if size(network_mat,1) == 2
        network_mat = network_mat;
    else
        n = size(network_mat(1:end-1,1:end-1),1);
        complete_edges = n*(n-1)/2;
        num_edges = sum(sum(network_mat(1:end-1,1:end-1)));
        clustering = num_edges/(2*complete_edges); %to on the bottom beacuse of double counting
        alts = 1:length(network_mat)-1;
        while clustering < desired_clustering
            prob_first = sum(network_mat(alts, :), 2);
            first_selection_vec = zeros(1,sum(prob_first));
            lower = 1;
            for i = 1:length(prob_first)
                upper = lower + prob_first(i)-1;
                first_selection_vec(lower:upper) = i;
                lower = upper+1;
            end
            first_selection_idx = randi(length(first_selection_vec),1);
            first_node = first_selection_vec(first_selection_idx);

            second_selection_vec = first_selection_vec(first_selection_vec ~= first_node);
            second_selection_idx = randi(length(second_selection_vec),1);
            second_node = second_selection_vec(second_selection_idx);

            network_mat(first_node, second_node) = 1;
            network_mat(second_node, first_node) = 1;

            clustering = sum(sum(network_mat(1:end-1,1:end-1)))/(2*complete_edges);
        end
    end
end
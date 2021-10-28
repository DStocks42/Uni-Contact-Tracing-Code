function resulting_isolates = isolating_contacts(positive_result, n_piv, traced_nodes)
    resulting_isolates = [];
    isolate = [];
    piv_contacts = find(positive_result);
    for i = 1:n_piv
        for j = 1:length(piv_contacts)
            isolate_j = unique(cell2mat(traced_nodes(piv_contacts(j), :)));
            new_isolates = isolate_j(~ismember(isolate_j, isolate));
            isolate = [isolate, new_isolates];
        end
        resulting_isolates = [resulting_isolates, isolate];
        piv_contacts = isolate;
    end
end
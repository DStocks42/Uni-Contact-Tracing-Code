clf;
staff_data_mix = readmatrix('CONQ_staff_edit.csv');
staff_data = unique_data(staff_data_mix);

sin_bin_staff = [133, 187, 346, 439]; %incomplete data entries 
staff_data(sin_bin_staff, :) = [];

prob_asym = 0.589;
prob_age = [0.121, 0.114, 0.132, 0.131, 0.128, 0.134, 0.106, 0.134];
prob_age_asmy = [0.191, 0.147, 0.162, 0.147, 0.074, 0.103, 0.074, 0.074];
prob_asmy_age = (prob_age_asmy.*prob_asym)./prob_age;

%r0 1 = 1.6
%r0 2 = 3.1
%r0 3 = 5.6
%r0 4 = 7.9

alpha_trans = [1.6e-4, 3.1e-4, 5.6e-4, 7.9e-4];

sim_period = 30; %the simulation period
num_res = 1; %the number of realisations

filename_inf = {'total_num_infected_staff_R01_TT1_all', 'total_num_infected_staff_R02_TT1_all', 'total_num_infected_staff_R03_TT1_all', 'total_num_infected_staff_R04_TT1_all'};
filename_untraced = {'total_num_untraced_staff_R01_TT1_all', 'total_num_untraced_staff_R02_TT1_all', 'total_num_untraced_staff_R03_TT1_all', 'total_num_untraced_staff_R04_TT1_all'};
for a = 1:length(alpha_trans)
    num_infected_total = zeros(sim_period,length(staff_data));
    num_infected_total(1,:) = 1;
    num_untraced_staff = zeros(length(staff_data), num_res);
    num_isolating  = zeros(sim_period,length(staff_data));
    true_c_num = zeros(1, length(staff_data));
    ego_trans = zeros(length(staff_data), num_res);

    deg_num_staff = zeros(1,length(staff_data));
    tot_alt_frq_staff = zeros(1,length(staff_data));
    tot_alt_dur_staff = zeros(1,length(staff_data));
    time_ego_isolate_staff = zeros(length(staff_data), num_res);
    num_first_time_total_staff = zeros(1, num_res);
    ego_infected_period_staff = zeros(1,length(staff_data));
    for s = 1:num_res
        isolated_nodes = [];
        for i = 1:length(staff_data)

            node_frq_ind = staff_data(i,28:17:521);
            node_dur_ind = staff_data(i,30:17:523);
            node_frq_grp = staff_data(i,547:22:635);
            node_dur_grp = staff_data(i,545:22:633);


            node_frq_ind(isnan(node_frq_ind)) = [];
            node_dur_ind(isnan(node_dur_ind)) = [];
            node_dur_grp(isnan(node_dur_grp)) = [];
            node_frq_grp(isnan(node_frq_grp)) = [];

            node_frq = [node_frq_ind, node_frq_grp];
            node_dur = [node_dur_ind, node_dur_grp];

            node_frq_dummy = node_frq;
            node_frq_dummy(node_frq_dummy == 1) = randi([8,14]);
            node_frq_dummy(node_frq_dummy == 2) = randi([4,6]);
            node_frq_dummy(node_frq_dummy == 3) = 2;
            node_frq_dummy(node_frq_dummy == 4) = 14;

            node_dur(node_dur == 1) = unifrnd(0,9);
            node_dur(node_dur == 2) = unifrnd(10,60);
            node_dur(node_dur == 3) = unifrnd(60,240);
            node_dur(node_dur == 4) = unifrnd(240,480);

            num_nodes = length(node_frq) + 1;

            if num_nodes == 1
                isolated_nodes = [isolated_nodes, i];
                continue
            end
%             if length(node_frq) ~= num_nodes-1 || length(node_dur) ~= num_nodes-1
%                continue
%             end
            network_mat = zeros(num_nodes, num_nodes);
            network_mat(:,end) = 1;
            network_mat(end,:) = 1;
            network_mat(end,end) = 0;

            network_mat = preferential_builders(network_mat, 0);

            mat_f = [node_frq_dummy'; 0];
            dummy_f = network_mat.*mat_f;
            dummy_f(end,:) = mat_f;
            f_mat = (dummy_f + dummy_f')/2;

            mat_d = [node_dur'; 0];
            dummy_d = network_mat.*mat_d;
            dummy_d(end,:) = mat_d;
            d_mat = (dummy_d + dummy_d')/2;

            alt_age_node_grp_cols = [];
            lower_bracket = 526;
            upper_bracket = 532;
            while upper_bracket <= 620
                alt_age_node_grp_cols = [alt_age_node_grp_cols, lower_bracket:upper_bracket];
                lower_bracket = lower_bracket + 22;
                upper_bracket = lower_bracket + 6;
            end

            alt_age_node_grp_mixed = staff_data(i,alt_age_node_grp_cols);
            grp_size_node = grp_size_builder(alt_age_node_grp_mixed); %we work out each group size by the age allocations
            grp_size_node(grp_size_node == 0) = [];
            mat_s = ones(num_nodes,1);
            mat_s(end-length(grp_size_node):end-1) = grp_size_node;
            group_idx = mat_s > 1;

            alt_age_node_ind = staff_data(i,14:17:507);
            alt_age_node_ind = age_interpolater(alt_age_node_ind);
            alt_age = age_expander(alt_age_node_grp_mixed, alt_age_node_ind);
            alt_age = cell2mat(alt_age);

            node_age_staff = [alt_age, staff_data(i,2)];

            if length(node_age_staff) < num_nodes
                continue
            end

            node_age_staff(node_age_staff <= 9) = 1;
            node_age_staff(node_age_staff >= 10 & node_age_staff <= 19) = 2;
            node_age_staff(node_age_staff >= 20 & node_age_staff <= 29) = 3;
            node_age_staff(node_age_staff >= 30 & node_age_staff <= 39) = 4;
            node_age_staff(node_age_staff >= 40 & node_age_staff <= 49) = 5;
            node_age_staff(node_age_staff >= 50 & node_age_staff <= 59) = 6;
            node_age_staff(node_age_staff >= 60 & node_age_staff <= 69) = 7;
            node_age_staff(node_age_staff >= 70) = 8;

            infected_period = cell(1,num_nodes);
            time_infected = cell(1,num_nodes);
            asym_yes_no = cell(1,num_nodes);
            inc_period = cell(1,num_nodes);
            asym_idx = [];
            for h = 1:num_nodes
                inc_period{1,h} = ceil(wblrnd(7.163, 3.0379, 1, mat_s(h)));
                infected_period{1,h} = cell2mat(inc_period(1,h)) + 3 + 10; %incubation + test + isolation period
                if h ~= num_nodes
                    time_infected{1,h} = zeros(1, mat_s(h));
                else
                    time_infected{1,h} = 1;
                end
                asym_yes_no{1,h} = unifrnd(0,1,1,mat_s(h)) <= prob_asmy_age(node_age_staff(h));
                asym_yes_no{1,end} = 1; %1 is asymptomatic
                if sum(cell2mat(asym_yes_no(1,h))) > 0
                    asym_idx = [asym_idx, h];
                end

            end
            num_recovered = zeros(1,num_nodes);
            positive_result = false(1, num_nodes);
            stop_isolating = zeros(1, num_nodes);
            isolating = false(1,num_nodes);

            first_time_c = node_frq == 5;
            untraceable = first_time_c;

            alt_size = mat_s(1:end-1);
            num_first_time_perday = sum(alt_size(first_time_c));
            num_first_time_c = num_first_time_perday*(sim_period-1); %add this in for new occuring contacts
            true_c_num(i) = num_first_time_c + num_nodes;

            num_first_time_inf_perday = zeros(1,sim_period);
            num_untracable_perday = zeros(1,sim_period);

            traced_nodes = cell(num_nodes, 5);

            infected = false(1,num_nodes);
            infected(end) = 1;
            contagious = false(1,num_nodes);
            num_infected = zeros(1,num_nodes);
            num_infected(end) = 1;
            ego_infected_period_staff(i) = cell2mat(infected_period(end));

            for j = 2:sim_period
                infected_idx = find(infected);
                contagious_idx = find(contagious);
                for h = 1:length(infected_idx)       
                    asym_dummy = cell2mat(asym_yes_no(1,infected_idx(h)));
                    inc_dummy = cell2mat(inc_period(1,infected_idx(h)));
                    time_infected_dummy = cell2mat(time_infected(1,infected_idx(h)));
                    time_infected_dummy(1:num_infected(infected_idx(h))) = time_infected_dummy(1:num_infected(infected_idx(h))) + 1;
                    if sum(time_infected_dummy > 4) > 0
                        contagious(infected_idx(h)) = 1;
                    end
                    if sum(asym_dummy(1:num_infected(infected_idx(h))) == 0 & time_infected_dummy(1:num_infected(infected_idx(h))) == (inc_dummy(1:num_infected(infected_idx(h))) + 3)) > 0
                        positive_result(infected_idx(h)) = 1;
                        stop_isolating(infected_idx(h)) = j + 10;
                    end

                    infected_period_dummy = cell2mat(infected_period(1,infected_idx(h)));
                    if sum(infected_period_dummy(1:num_infected(infected_idx(h))) - time_infected_dummy(1:num_infected(infected_idx(h)))) == 0
                        infected(infected_idx(h)) = 0;
                        contagious(infected_idx(h)) = 0;
                        time_infect_dummy(:) = 0;
                        num_infected(infected_idx(h)) = 0;
                        num_recovered(infected_idx(h)) = mat_s(infected_idx(h));
                        recovered = true(1,mat_s(infected_idx(h)));
                    else
                        recovered = infected_period_dummy - time_infected_dummy == 0;
                        time_infected_dummy(recovered) = 0;
                        num_infected(infected_idx(h)) = num_infected(infected_idx(h)) - sum(recovered);
                        num_recovered(infected_idx(h)) = num_recovered(infected_idx(h)) + sum(recovered);
                    end
                    if isolating(infected_idx(h)) == 0
                        num_untraced_staff(i,s) = num_untraced_staff(i,s) + sum(recovered);
                    end
                    time_infected{1,infected_idx(h)} = time_infected_dummy; 
                end
                all_isolates = isolating_contacts(positive_result, 1, traced_nodes);
                stop_isolating(all_isolates) = j + 14;
                positive_result = false(1,num_nodes); %reset 

                infected_delay = false(1,length(network_mat));
                contagious_d = d_mat(contagious, :);
                mat_dim = size(f_mat);
                prob_meet = f_mat/14;
%                 isolation by setting their probability of meeting to zero for
%                 everyone else to them and from them to everyone else (their column and row respectively)
                isolating = stop_isolating > j;
                if isolating(end) == 1
                    time_ego_isolate_staff(i,s) = j;
                    break
                end
                isolating(untraceable) = 0;
                prob_meet(isolating,:) = 0;
                prob_meet(:,isolating) = 0;

                non_symetric_meet = unifrnd(0,1,mat_dim(1), mat_dim(2)) < prob_meet;
                meet_yes_no = triu(non_symetric_meet,1) + triu(non_symetric_meet,1)'; %ensures the matrix is symetric
                for x = 1:length(contagious_idx)
                    time_infected_dummy2 = cell2mat(time_infected(1,contagious_idx(x)));
                    met_nodes = find(meet_yes_no(contagious_idx(x), :));
                    if isempty(met_nodes)
                        met_nodes = [];
                    end
                    if max(time_infected_dummy2) > 5
                        replace = mod(max(time_infected_dummy2), 5) + 1;
                        traced_nodes{contagious_idx(x), replace} = met_nodes;
                    else
                        traced_nodes{contagious_idx(x), max(time_infected_dummy2)} = met_nodes;
                    end
                end
                if sum(any(meet_yes_no(contagious, :), 1) & group_idx' == 1) > 0 %did an infected node have a group interaction?     
                    contagious_groups = contagious == 1 & group_idx' == 1;
                    if sum(contagious_groups) > 1 % this double loops eliminates the repeated group interactions
                        contagious_groups_idx = find(contagious_groups);
                        for x = 1:length(contagious_groups_idx)
                            other_groups_idx = contagious_groups_idx(contagious_groups_idx ~= x);
                            for y = 1:length(other_groups_idx)
                                repeated_group_interactions = meet_yes_no(contagious_groups_idx(x), :)' == meet_yes_no(:, other_groups_idx(y));
                                meet_yes_no(repeated_group_interactions, other_groups_idx(y)) = 0;
                            end
                        end
                    end
                    met_groups = any(meet_yes_no(contagious, :), 1) & group_idx' == 1;
                    met_groups_idx = find(met_groups);
                    num_groups_met = length(met_groups_idx);
                    resultant_trans = cell(1,length(network_mat));
                    for k = 1:num_groups_met
                        group_size = mat_s(met_groups_idx(k));
                        num_infected_group = num_infected(met_groups_idx(k));
                        infected_parties = any(meet_yes_no(:, met_groups_idx(k)), 2)' & contagious == 1;
                        infected_party_idx = find(infected_parties);        
                        for h = 1:length(infected_party_idx)
                            party_size = mat_s(infected_party_idx(h));
                            num_infected_party = num_infected(infected_party_idx(h));
                            duration_gathering = d_mat(infected_party_idx(h), met_groups_idx(k));
                            gathering_size = group_size + party_size;
                            time_per_person = duration_gathering/(gathering_size - 1); % minus one because people cannot interact with themsleves
                            prob_trans_uniform = 1 - exp(-time_per_person*alpha_trans(a));

                            num_susceptible_group = group_size - num_infected_group;
                            num_susceptible_party = party_size - num_infected_party;

                            trans_party2group_invidual = unifrnd(0,1,num_infected_party, num_susceptible_group) <= prob_trans_uniform;
                            if infected_party_idx(h) == num_nodes && j <= ego_infected_period_staff(i) %check if ego and if still first infection
                                ego_trans(i,s) = ego_trans(i,s) + sum(trans_party2group_invidual);
                            end
                            trans_group2group_invidual = unifrnd(0,1,num_infected_group, num_susceptible_group) <= prob_trans_uniform;
                            resultant_trans{1,met_groups_idx(k)} = [cell2mat(resultant_trans(1,met_groups_idx(k))); trans_party2group_invidual; trans_group2group_invidual]; % concatination all interactions of the node


                            trans_group2party_invidual = unifrnd(0,1,num_infected_group, num_susceptible_party) <= prob_trans_uniform;
                            trans_party2party_invidual = unifrnd(0,1,num_infected_party, num_susceptible_party) <= prob_trans_uniform;
                            resultant_trans{1,infected_party_idx(h)} = [cell2mat(resultant_trans(1,infected_party_idx(h))); trans_group2party_invidual; trans_party2party_invidual]; % concatination all interactions of the node

                        end
                    end
                    nodes_involved_idx = unique([met_groups_idx, infected_party_idx]);
                    for h = 1:length(nodes_involved_idx)
                        whos_infected = any(cell2mat(resultant_trans(1, nodes_involved_idx(h))), 1);
                        num_infected(nodes_involved_idx(h)) = num_infected(nodes_involved_idx(h)) + sum(whos_infected);
                        if sum(whos_infected) > 0 
                            infected_delay(nodes_involved_idx(h)) = 1;
                        end
                    end
                    meet_yes_no(contagious, met_groups) = 0; %this is eliminating the accounted for group interactions
                end
                contagious_mat_dim = size(contagious_d);

                possible_trans = meet_yes_no;
                infected_asym = infected_idx(ismember(infected_idx, asym_idx));
                for p = infected_asym
                    met_p = meet_yes_no(p,:) == 1;
                    possible_trans(p, met_p) = unifrnd(0,1,1,sum(met_p)) <= 0.65;
                end

                prob_trans = 1 - exp(-contagious_d*alpha_trans(a));
                transmit = meet_yes_no(contagious, :) & (unifrnd(0,1,contagious_mat_dim(1), contagious_mat_dim(2)) < prob_trans) & possible_trans(contagious, :);

                if j < ego_infected_period_staff(i) && j > 4
                    ego_trans(i,s) = ego_trans(i,s) + sum(transmit(end,:) & infected == 0);
                end

                update = infected == 0 & any(transmit, 1);
                infected(update) = 1;
                num_infected(update) = 1;
                infected = infected | infected_delay;

                if num_first_time_perday ~= 0 % resests the the first time contacts to be uninfected after everday, is equivalent to a new person
                    num_first_time_inf_perday(j) = sum(num_infected(first_time_c)); % recording transmissions to first time contacts
                    num_infected(first_time_c) = 0;
                    infected(first_time_c) = 0;
                    num_untraced_staff(i,s) = num_untraced_staff(i,s) + num_first_time_inf_perday(j);
                    num_first_time_total_staff(s) = num_first_time_total_staff(s) + num_first_time_inf_perday(j);
                end
            end
            num_infected_total(j,i) = sum(num_infected) + sum(num_first_time_inf_perday);
        end
    end
    num_infected_total(:,isolated_nodes) = [];
    num_isolating(:,isolated_nodes) = [];
    num_untraced_staff(isolated_nodes,:,:) = [];
    ego_trans(isolated_nodes,:,:) = [];
    time_ego_isolate_staff(isolated_nodes,:) = [];
    ego_infected_period_staff(isolated_nodes) = [];

    total_ego_trans_staff = sum(ego_trans,1);
    total_ego_trans_staff_re = reshape(total_ego_trans_staff, [num_res 1]);
    filename_inf_out = [filename_inf(a), '.csv'];
    writematrix(total_ego_trans_staff_re, cell2mat(filename_inf_out))
    total_num_untraced_staff = sum(num_untraced_staff, 1);
    filename_untraced_out = [filename_untraced(a), '.csv'];
    total_num_untraced_staff_re = reshape(total_num_untraced_staff, [num_res 1]);
    writematrix(total_num_untraced_staff_re, cell2mat(filename_untraced_out))
end
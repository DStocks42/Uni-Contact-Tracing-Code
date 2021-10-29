
%A1
ego_trans_staff_A1_TT0 = readmatrix('total_ego_trans_staff_R01_TT0.csv');
ego_trans_staff_A1_TT0_asym = readmatrix('total_num_infected_staff_R01_TT0_asym.csv');
ego_trans_staff_A1_TT0_all = readmatrix('total_num_infected_staff_R01_TT0_all.csv');

ego_trans_staff_A1_TT1 = readmatrix('total_ego_trans_staff_R01_TT1.csv');
ego_trans_staff_A1_TT1_asym = readmatrix('total_ego_trans_staff_R01_TT1_asym.csv');
ego_trans_staff_A1_TT1_all = readmatrix('total_num_infected_staff_R01_TT1_all.csv');

ego_trans_stu_A1_TT0 = readmatrix('total_ego_trans_stu_R01_TT0.csv');
ego_trans_stu_A1_TT0_asym = readmatrix('total_num_infected_stu_R01_TT0_asym.csv');
ego_trans_stu_A1_TT0_all = readmatrix('total_num_infected_stu_R01_TT0_all.csv');

ego_trans_stu_A1_TT1 = readmatrix('total_ego_trans_stu_R01_TT1.csv');
ego_trans_stu_A1_TT1_asym = readmatrix('total_ego_trans_stu_R01_TT1_asym.csv');
ego_trans_stu_A1_TT1_all = readmatrix('total_num_infected_stu_R01_TT1_all.csv');

%A2
ego_trans_staff_A2_TT0 = readmatrix('total_ego_trans_staff_R02_TT0.csv');
ego_trans_staff_A2_TT0_asym = readmatrix('total_num_infected_staff_R02_TT0_asym.csv');
ego_trans_staff_A2_TT0_all = readmatrix('total_num_infected_staff_R02_TT0_all.csv');

ego_trans_staff_A2_TT1 = readmatrix('total_ego_trans_staff_R02_TT1.csv');
ego_trans_staff_A2_TT1_asym = readmatrix('total_ego_trans_staff_R02_TT1_asym.csv');
ego_trans_staff_A2_TT1_all = readmatrix('total_num_infected_staff_R02_TT1_all.csv');

ego_trans_stu_A2_TT0 = readmatrix('total_ego_trans_stu_R02_TT0.csv');
ego_trans_stu_A2_TT0_asym = readmatrix('total_num_infected_stu_R02_TT0_asym.csv');
ego_trans_stu_A2_TT0_all = readmatrix('total_num_infected_stu_R02_TT0_all.csv');

ego_trans_stu_A2_TT1 = readmatrix('total_ego_trans_stu_R02_TT1.csv');
ego_trans_stu_A2_TT1_asym = readmatrix('total_ego_trans_stu_R02_TT1_asym.csv');
ego_trans_stu_A2_TT1_all = readmatrix('total_num_infected_stu_R02_TT1_all.csv');

%A3
ego_trans_staff_A3_TT0 = readmatrix('total_ego_trans_staff_R03_TT0.csv');
ego_trans_staff_A3_TT0_asym = readmatrix('total_num_infected_staff_R03_TT0_asym.csv');
ego_trans_staff_A3_TT0_all = readmatrix('total_num_infected_staff_R03_TT0_all.csv');

ego_trans_staff_A3_TT1 = readmatrix('total_ego_trans_staff_R03_TT1.csv');
ego_trans_staff_A3_TT1_asym = readmatrix('total_ego_trans_staff_R03_TT1_asym.csv');
ego_trans_staff_A3_TT1_all = readmatrix('total_num_infected_staff_R03_TT1_all.csv');


ego_trans_stu_A3_TT0 = readmatrix('total_ego_trans_stu_R03_TT0.csv');
ego_trans_stu_A3_TT0_asym = readmatrix('total_num_infected_stu_R03_TT0_asym.csv');
ego_trans_stu_A3_TT0_all = readmatrix('total_num_infected_stu_R03_TT0_all.csv');

ego_trans_stu_A3_TT1 = readmatrix('total_ego_trans_stu_R03_TT1.csv');
ego_trans_stu_A3_TT1_asym = readmatrix('total_ego_trans_stu_R03_TT1_asym.csv');
ego_trans_stu_A3_TT1_all = readmatrix('total_num_infected_stu_R03_TT1_all.csv');

%A4
ego_trans_staff_A4_TT0 = readmatrix('total_ego_trans_staff_R04_TT0.csv');
ego_trans_staff_A4_TT0_asym = readmatrix('total_num_infected_staff_R04_TT0_asym.csv');
ego_trans_staff_A4_TT0_all = readmatrix('total_num_infected_staff_R04_TT0_all.csv');

ego_trans_staff_A4_TT1 = readmatrix('total_ego_trans_staff_R04_TT1.csv'); %requires resim
ego_trans_staff_A4_TT1_asym = readmatrix('total_num_infected_staff_R04_TT1_asym.csv');
ego_trans_staff_A4_TT1_all = readmatrix('total_num_infected_staff_R04_TT1_all.csv');


ego_trans_stu_A4_TT0 = readmatrix('total_ego_trans_stu_R04_TT0.csv');
ego_trans_stu_A4_TT0_asym = readmatrix('total_num_infected_stu_R04_TT0_asym.csv');
ego_trans_stu_A4_TT0_all = readmatrix('total_num_infected_stu_R04_TT0_all.csv');

ego_trans_stu_A4_TT1 = readmatrix('total_ego_trans_stu_R04_TT1.csv');
ego_trans_stu_A4_TT1_asym = readmatrix('total_ego_trans_stu_R04_TT1_asym.csv');
ego_trans_stu_A4_TT1_all = readmatrix('total_num_infected_stu_R04_TT1_all.csv');

%combined numbers
A1_TT0 = ego_trans_staff_A1_TT0 + ego_trans_stu_A1_TT0;
A1_TT0_asym = ego_trans_staff_A1_TT0_asym + ego_trans_stu_A1_TT0_asym;
A1_TT0_all = ego_trans_staff_A1_TT0_all + ego_trans_stu_A1_TT0_all;
A1_TT1 = ego_trans_staff_A1_TT1 + ego_trans_stu_A1_TT1;
A1_TT1_asym = ego_trans_staff_A1_TT1_asym + ego_trans_stu_A1_TT1_asym;
A1_TT1_all = ego_trans_staff_A1_TT1_all + ego_trans_stu_A1_TT1_all;

A2_TT0 = ego_trans_staff_A2_TT0 + ego_trans_stu_A2_TT0;
A2_TT0_asym = ego_trans_staff_A2_TT0_asym + ego_trans_stu_A2_TT0_asym;
A2_TT0_all = ego_trans_staff_A2_TT0_all + ego_trans_stu_A2_TT0_all;
A2_TT1 = ego_trans_staff_A2_TT1 + ego_trans_stu_A2_TT1;
A2_TT1_asym = ego_trans_staff_A2_TT1_asym + ego_trans_stu_A2_TT1_asym;
A2_TT1_all = ego_trans_staff_A2_TT1_all + ego_trans_stu_A2_TT1_all;

A3_TT0 = ego_trans_staff_A3_TT0 + ego_trans_stu_A3_TT0;
A3_TT0_asym = ego_trans_staff_A3_TT0_asym + ego_trans_stu_A3_TT0_asym;
A3_TT0_all = ego_trans_staff_A3_TT0_all + ego_trans_stu_A3_TT0_all;
A3_TT1 = ego_trans_staff_A3_TT1 + ego_trans_stu_A3_TT1;
A3_TT1_asym = ego_trans_staff_A3_TT1_asym + ego_trans_stu_A3_TT1_asym;
A3_TT1_all = ego_trans_staff_A3_TT1_all + ego_trans_stu_A3_TT1_all;

A4_TT0 = ego_trans_staff_A4_TT0 + ego_trans_stu_A4_TT0;
A4_TT0_asym = ego_trans_staff_A4_TT0_asym + ego_trans_stu_A4_TT0_asym;
A4_TT0_all = ego_trans_staff_A4_TT0_all + ego_trans_stu_A4_TT0_all;
A4_TT1 = ego_trans_staff_A4_TT1 + ego_trans_stu_A4_TT1;
A4_TT1_asym = ego_trans_staff_A4_TT1_asym + ego_trans_stu_A4_TT1_asym;
A4_TT1_all = ego_trans_staff_A4_TT1_all + ego_trans_stu_A4_TT1_all;

%obataining the mean of each alpha
total_num_networks = 1317;

Reff_A1_TT0 = mean(A1_TT0, 1)/total_num_networks;
Reff_A2_TT0 = mean(A2_TT0, 1)/total_num_networks;
Reff_A3_TT0 = mean(A3_TT0, 1)/total_num_networks;
Reff_A4_TT0 = mean(A4_TT0, 1)/total_num_networks;

Reff_A1_TT1 = mean(A1_TT1, 1)/total_num_networks;
Reff_A2_TT1 = mean(A2_TT1, 1)/total_num_networks;
Reff_A3_TT1 = mean(A3_TT1, 1)/total_num_networks;
Reff_A4_TT1 = mean(A4_TT1, 1)/total_num_networks;

Reff_A1_TT0_asym = mean(A1_TT0_asym, 1)/total_num_networks;
Reff_A2_TT0_asym = mean(A2_TT0_asym, 1)/total_num_networks;
Reff_A3_TT0_asym = mean(A3_TT0_asym, 1)/total_num_networks;
Reff_A4_TT0_asym = mean(A4_TT0_asym, 1)/total_num_networks;

Reff_A1_TT1_asym = mean(A1_TT1_asym, 1)/total_num_networks;
Reff_A2_TT1_asym = mean(A2_TT1_asym, 1)/total_num_networks;
Reff_A3_TT1_asym = mean(A3_TT1_asym, 1)/total_num_networks;
Reff_A4_TT1_asym = mean(A4_TT1_asym, 1)/total_num_networks;

sym_data = [A1_TT0(:,1), A1_TT1(:,1), A2_TT0(:,1), A2_TT1(:,1), A3_TT0(:,1), A3_TT1(:,1), A4_TT0(:,1), A4_TT1(:,1)]/total_num_networks;

asym_data = [A1_TT0_asym(:,1), A1_TT1_asym(:,1), A2_TT0_asym(:,1), A2_TT1_asym(:,1), A3_TT0_asym(:,1), A3_TT1_asym(:,1), A4_TT0_asym(:,1), A4_TT1_asym(:,1)]/total_num_networks;

all_data = [A1_TT0_all, A1_TT1_all, A2_TT0_all, A2_TT1_all, A3_TT0_all, A3_TT1_all, A4_TT0_all, A4_TT1_all]/total_num_networks;

fig = figure('Name','Results');
subplot(1,3,1)
positions =[1 1.2 1.5 1.7 2 2.2 2.5 2.7];
boxplot(sym_data, 'positions', positions)
set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) mean(positions(5:6)) mean(positions(7:8))])
set(gca,'TickLabelInterpreter', 'tex','xticklabel',{'α_1','α_2', 'α_3', 'α_4'})
ylim([0.1,1.1])
title('(a)')


subplot(1,3,2)
positions =[1 1.2 1.5 1.7 2 2.2 2.5 2.7];
boxplot(asym_data, 'positions', positions)
set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) mean(positions(5:6)) mean(positions(7:8))])
set(gca,'TickLabelInterpreter', 'tex','xticklabel',{'α_1','α_2', 'α_3', 'α_4'})
ylim([0.1,1.1])
title('(b)')
% boxplot(asym_data_TT0, 'Labels', {1,2,3,4})
% hold on
% boxplot(asym_data_TT1, 'Labels', {1,2,3,4})
% hold off

subplot(1,3,3)
positions =[1 1.2 1.5 1.7 2 2.2 2.5 2.7];
boxplot(all_data, 'positions', positions)
set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) mean(positions(5:6)) mean(positions(7:8))])
set(gca,'TickLabelInterpreter', 'tex','xticklabel',{'α_1','α_2', 'α_3', 'α_4'})
ylim([0.1,1.1])
title('(c)')

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Average Number of Secondary Cases');
xlabel(han,'Transmission Rate');

savepdf('Results')
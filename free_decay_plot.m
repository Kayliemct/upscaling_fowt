clear all; close all; clc

% Tn_IEA = readtable('IEA-15-240-RWT-UMaineSemi-Tn.xlsm');  
Tn_OC4 = readtable('5MW_OC4Semi_Tn_tower off.xls');


%%
% time = Tn_IEA.Time;
time = Tn_OC4.Time;

%%

figure(1)
plot(time, Tn_OC4.PtfmPitch,'LineWidth',2)
hold on 
axis([0 600 -10 10])
xlabel('Time (s)')
ylabel('Pitch angle (deg)')
set(gca,'FontSize',20)
hold off

%%
fie = gca 
% exportgraphics(fi,'IEA_free_decay.png','Resolution',300)
exportgraphics(fie,'OC4_free_decay.png','Resolution',300)
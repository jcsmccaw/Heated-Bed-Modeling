%% 315 Project Data Processing
clear; clc; close all

% Load Data:
filepath = 'C:\Users\JC\Documents\Purdue Spring 2019\ME 315\Project Docs\Data\';
filenames = {'Case1_no_insulation.xlsx'; 'linear_insulation_case2.xlsx'; 'case3_constant insulation.xlsx'; 'case4_rectangular.xlsx'; 'case5_rectangular_with_edges.xlsx'};
files = cell(length(filenames), 1);

data = cell(5,1);
time = cell(5,1);

for i = [1 2 3 5]
    filename = strcat(filepath, filenames{i});
    temp = xlsread(filename);
    time{i} = temp(:,1); % We only need one time vector for each sample
    data{i} = temp(:,2:2:end); % all the thermocouple data, 1:16
end
i = 4;
filename = strcat(filepath, filenames{i});
temp = xlsread(filename);
time{i} = temp(:,1); % We only need one time vector for each sample
data{i} = temp(:,3:end); % all the thermocouple data, 1:16

%% Plot Thermocouple Locations:
% y = [3 3 3 3 2 2 2 2 1 1 1 1 0 0 0 0];
% x = [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4];
% figure
% hold on
% for i = 1:16
%     x(i) = x(i)-1;
%     couple = num2str(i);
%     thermocouple = strcat('T', couple);
%     jitter = 0.15;
%     plot(x(i), y(i), 'o')
%     text(x(i)+ jitter, y(i)- jitter, thermocouple)
% end
% xlim([-0.25 3.75])
% ylim([-0.5 3.5])
% mn = -0;
% mx = 3;
% plot([mn mn], [mn mx], '-k')
% plot([mn mx], [mx mx], '-k')
% plot([mx mx], [mx mn], '--k')
% plot([mx mn], [mn mn], '--k')
% 
% box on
% set(gca, 'Xtick', [])
% set(gca, 'Ytick', [])
% title('Thermocouple Positions')
% %% Plot Raw Data:
% % Case 1:
% figure
% hold on
% for i = 1:15
%     plot(time{1}, (data{1}(:,i)))
% end
% % 
% % for i = 11:15
% %     plot(time{1}, data{1}(:,i))
% % end
% legend({'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15'})
% title('No Insulation')
% legend('boxoff')
% ylim([0 max(max(data{5}(:,:)))+5])
% xlim([0 time{5}(end)+110])
% set(gca, 'Xtick', [0 50 100 150 200 250 300 350])
% legend('Orientation', 'vertical', 'Location', 'east')
% 
% 
% % Case 2:
% figure
% hold on
% for i = 1:15
%     plot(time{2}, data{2}(:,i))
% end
% title('Linear Insulation Profile')
% legend({'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15'})
% 
% legend('boxoff')
% ylim([0 max(max(data{5}(:,:)))+5])
% xlim([0 time{5}(end)+90])
% set(gca, 'Xtick', [0 50 100 150 200 250 300 350])
% legend('Orientation', 'vertical', 'Location', 'east')
% 
% % Case 3:
% figure
% hold on
% for i = 1:16
%     plot(time{3}, data{3}(:,i))
% end
% title('Constant Insulation Profile')
% 
% legend({'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16'})
% 
% legend('boxoff')
% ylim([0 max(max(data{5}(:,:)))+5])
% xlim([0 time{5}(end)+90])
% set(gca, 'Xtick', [0 50 100 150 200 250 300 350])
% legend('Orientation', 'vertical', 'Location', 'east')
% 
% % Case 4:
% figure
% hold on
% for i = 1:16
%     plot(time{4}, data{4}(:,i))
% end
% title('Stepped Insulation Profile')
% legend({'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16'})
% 
% legend('boxoff')
% ylim([0 max(max(data{5}(:,:)))+5])
% xlim([0 time{5}(end)+90])
% set(gca, 'Xtick',[0 50 100 150 200 250 300 350])
% legend('Orientation', 'vertical', 'Location', 'east')
% 
% % Case 5:
% figure
% hold on
% for i = 1:16
%     plot(time{5}, data{5}(:,i))
% end
% 
% legend({'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16'})
% title('Stepped Insulation Profile with Raised Edges')
% legend('boxoff')
% ylim([0 max(max(data{5}(:,:)))+5])
% xlim([0 time{5}(end)+90])
% set(gca, 'Xtick', [0 50 100 150 200 250 300 350])
% legend('Orientation', 'vertical', 'Location', 'east')

% %% Normalized Data and Plot:
% 
% % Case 1:
% figure
% hold on
% for i = 1:15
%     plot(time{1}, (data{1}(:,i) - data{1}(1,i))/data{1}(1,i))
% end
% % 
% % for i = 11:15
% %     plot(time{1}, data{1}(:,i))
% % end
% legend({'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15'})
% title('No Insulation')
% legend('boxoff')
% ylim([0 2])
% xlim([0 time{5}(end)+110])
% set(gca, 'Xtick', [0 50 100 150 200 250 300 350])
% legend('Orientation', 'vertical', 'Location', 'east')
% 
% 
% % Case 2:
% figure
% hold on
% for i = 1:15
%     plot(time{2}, (data{2}(:,i)- data{2}(1,i))/data{2}(1,i))
% end
% title('Linear Insulation Profile')
% legend({'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15'})
% 
% legend('boxoff')
% ylim([0 2])
% xlim([0 time{5}(end)+90])
% set(gca, 'Xtick', [0 50 100 150 200 250 300 350])
% legend('Orientation', 'vertical', 'Location', 'east')
% 
% % Case 3:
% figure
% hold on
% for i = 1:16
%     plot(time{3}, (data{3}(:,i)- data{3}(1,i))/data{3}(1,i))
% end
% title('Constant Insulation Profile')
% 
% legend({'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16'})
% 
% legend('boxoff')
% ylim([0 2])
% xlim([0 time{5}(end)+90])
% set(gca, 'Xtick', [0 50 100 150 200 250 300 350])
% legend('Orientation', 'vertical', 'Location', 'east')
% 
% % Case 4:
% figure
% hold on
% for i = 1:16
%     plot(time{4}, (data{4}(:,i)- data{4}(1,i))/data{4}(1,i))
% end
% title('Stepped Insulation Profile')
% legend({'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16'})
% 
% legend('boxoff')
% ylim([0 2])
% xlim([0 time{5}(end)+90])
% set(gca, 'Xtick',[0 50 100 150 200 250 300 350])
% legend('Orientation', 'vertical', 'Location', 'east')
% 
% % Case 5:
% figure
% hold on
% for i = 1:16
%     plot(time{5}, (data{5}(:,i)- data{5}(1,i))/data{5}(1,i))
% end
% 
% legend({'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16'})
% title('Stepped Insulation Profile with Raised Edges')
% legend('boxoff')
% ylim([0 2])
% xlim([0 time{5}(end)+90])
% set(gca, 'Xtick', [0 50 100 150 200 250 300 350])
% legend('Orientation', 'vertical', 'Location', 'east')
% 
% 
% %% Calculate Gradient from center to edges for each case to determine Insulation Viability:
% % To get a rough estimate, we want the difference between T16 and T1 as an
% % approximate gradient from the center to the corners. 
% 
% % For cases 1 and 2, we can use the finite difference equation to
% % approximate T16:
%         % To do this we need to compute the volumetric heat generation
%         % term, qdot. 
%         
% k = 237; % W/mK;
% dx = 0.02667; % m 
% dy = 0.02167; % m
% dz = 1.8/1000; % m
% 
% kx = k*dx/dy;
% ky = k*dy/dx;
% 
% % Creating a system of 2 equations and 2 unknowns to get qdot and h from
% % the no insulation case at thermocouple (2,2) (T6) and (3,3) (T11):
% % steady state is well-approximated at t=220sec. 
% data{1}(220,10) = mean([data{1}(220,11), data{1}(220,9), data{1}(220,6), data{1}(220,14)]);
% mat = [dx*dy*(data{1}(220,7) - data{1}(1,7))/dz, -dx*dy; dx*dy*(data{1}(220,10)/dz - data{1}(1,10)), -dx*dy];
% rix = [kx*(data{1}(220,11) + data{1}(220,3) - 2*data{1}(220,7)) + ky*(data{1}(220,8) + data{1}(220,6) - 2*data{1}(220,7)); kx*(data{1}(220,14) + data{1}(220,6) - 2*data{1}(220,10)) + ky*(data{1}(220,11) + data{1}(220,9) - 2*data{1}(220,10))];
% % mat[h;q]= trix
% solved = -inv(mat)*rix;
% h = solved(1);
% qdot = solved(2);
% 
% 
% % Finding T16:
% T16_1 = (1/2).*(data{1}(:,15) + data{1}(:,12)) + qdot*(dx*dy*dz)/4;
% T16_2 = (1/2).*(data{2}(:,15) + data{2}(:,12)) + qdot*(dx*dy*dz)/4;
% %% Plotting: 
% % Reformatting the data:
% vid_data = zeros(4,4,max(size(data{1})));
% vid_data(:,1,:) = [data{1}(:,1)'; data{1}(:,2)'; data{1}(:,3)'; data{1}(:,4)'];
% vid_data(:,2,:) = [data{1}(:,5) data{1}(:,6) data{1}(:,7) data{1}(:,8)]';
% vid_data(:,3,:) = [data{1}(:,9) data{1}(:,10) data{1}(:,11) data{1}(:,12)]';
% vid_data(:,4,:) = [data{1}(:,13) data{1}(:,14) data{1}(:,15) T16_1]';
% % 
% figure
% surf(vid_data(:,:, 1));
% set(gca,'nextplot','replacechildren'); 
% % vid = VideoWriter('transient.avi');
% % open(vid);
% dt = 0.001; % s
% end_time = 80; % sec
% time_steps = ceil(end_time/dt);
% 
% 
% [xx, yy] = meshgrid(1:0.1:10); 
% for i =1:360
%     surf(interp2(vid_data(:,:, i), xx, yy), 'EdgeColor', 'interp');
%     view(-37.5+180,30)
%     zlabel('^\circC')
%     xlabel('X, m')
%     ylabel('Y, m')
%     tstring = strcat('t = ', num2str(i));
%     title({'Transient Conduction for Quarter Plate, Model Validation',tstring})
%     tops = max(vid_data(:,:,i));
%     bottoms = min(vid_data(:,:,i));
%     zlim([bottoms(end)-2 (tops(1)+2)])
% %     frame = getframe(gcf);
% %     writeVideo(vid, frame);
% end
% close(vid);
% 
% reconstructed = {T16_1, T16_2};
% diff = cell(5,1);
% 
% for i = 1:5
%      if(i==1)
%          diff{i} = data{i}(75:end, 15) - data{i}(75:end,1);
%      else
%         diff{i} = data{i}(:,15) - data{i}(:,1); 
%      end
%     
% end
% 
% figure
% hold on
% for i = 1:5
%     if(i==1)
%         plot(time{i}(75:end)-75, diff{i})
%     else
%         plot(time{i}, diff{i})
%     end
%     
% end
% legend({'Case 1: No Insulation', 'Case 2: Linear Profile', 'Case 3: Constant Profile', 'Case 4: Stepped Profile', 'Case 5: Stepped with Raised Edges'})


%% Time to Reach set point:
figure
hold on

plot(time{1}(80:end)-75, data{1}(80:end,15))
plot(time{2}, data{2}(:,15))
plot(time{3}+5, data{3}(:,15))
plot(time{4}-2, data{4}(:,15))
plot(time{5}+3, data{5}(:,15))

plot([0 360], [50 50], '--k')
legend({'Case 1: No Insulation', 'Case 2: Linear Profile', 'Case 3: Constant Profile', 'Case 4: Stepped Profile', 'Case 5: Stepped with Raised Edges'})


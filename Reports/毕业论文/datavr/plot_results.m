clear;
close all;
clc;

load('data');
%% 2d temporal
temp_2d = [];
for i_col = 1:3:size(RECORD2d, 2)
    temp_2d = [temp_2d, sortrows(RECORD2d(:, i_col:i_col + 2), [1, 3])];
end
n_sub = size(temp_2d, 2) / 3;

temp_2d_no_sound = zeros(n_sub, length(temp_2d)/6);
temp_2d_collision_sound = zeros(n_sub, length(temp_2d)/6);
temp_2d_all_sound = zeros(n_sub, length(temp_2d)/6);

c_delay = 1;
for i_delay = 1:6:length(temp_2d)
    cur_row = mean(temp_2d(i_delay : i_delay + 1, :));
    temp_2d_no_sound(:, c_delay) = cur_row(2:3:end);
    
    cur_row = mean(temp_2d(i_delay + 2 : i_delay + 3, :));
    temp_2d_collision_sound(:, c_delay) = cur_row(2:3:end);
    
    cur_row = mean(temp_2d(i_delay + 4 : i_delay + 5, :));
    temp_2d_all_sound(:, c_delay) = cur_row(2:3:end);
    
    c_delay = c_delay + 1;
end

% temp_2d_no_sound(1, :) = [];
subplot(1, 3, 1);
boxplot([temp_2d_no_sound(:, 1), ...
    temp_2d_no_sound(:, 2), ...
    temp_2d_no_sound(:, 3), ...
    temp_2d_no_sound(:, 4), ...
    temp_2d_no_sound(:, 5), ...
    temp_2d_no_sound(:, 6), ...
    temp_2d_no_sound(:, 7), ...
    temp_2d_no_sound(:, 8), ...
    temp_2d_no_sound(:, 9)]);
xlabel('Delay (ms)');
ylabel('Rating (%)');
set(gca, 'XTickLabels', {'0','50','100','150','200','250','300','350','400'});
title('no sound');

subplot(1, 3, 2);
boxplot([temp_2d_collision_sound(:, 1), ...
    temp_2d_collision_sound(:, 2), ...
    temp_2d_collision_sound(:, 3), ...
    temp_2d_collision_sound(:, 4), ...
    temp_2d_collision_sound(:, 5), ...
    temp_2d_collision_sound(:, 6), ...
    temp_2d_collision_sound(:, 7), ...
    temp_2d_collision_sound(:, 8), ...
    temp_2d_collision_sound(:, 9)]);
xlabel('Delay (ms)');
ylabel('Rating (%)');
set(gca, 'XTickLabels', {'0','50','100','150','200','250','300','350','400'});
hold on;
set(gca,'ticklength',[0,0]);
title('with collision sound');

subplot(1, 3, 3);
boxplot([temp_2d_all_sound(:, 1), ...
    temp_2d_all_sound(:, 2), ...
    temp_2d_all_sound(:, 3), ...
    temp_2d_all_sound(:, 4), ...
    temp_2d_all_sound(:, 5), ...
    temp_2d_all_sound(:, 6), ...
    temp_2d_all_sound(:, 7), ...
    temp_2d_all_sound(:, 8), ...
    temp_2d_all_sound(:, 9)]);
xlabel('Delay (ms)');
ylabel('Rating (%)');
set(gca, 'XTickLabels', {'0','50','100','150','200','250','300','350','400'});
title('with collision and rolling sound');
%% 3d temporal
temp_3d = [];
for i_col = 1:3:size(RECORD3d, 2)
    temp_3d = [temp_3d, sortrows(RECORD3d(:, i_col:i_col + 2), [1, 3])];
end
n_sub = size(temp_3d, 2) / 3;

temp_3d_no_sound = zeros(n_sub, length(temp_3d)/6);
temp_3d_collision_sound = zeros(n_sub, length(temp_3d)/6);
temp_3d_all_sound = zeros(n_sub, length(temp_3d)/6);

c_delay = 1;
for i_delay = 1:6:length(temp_3d)
    cur_row = mean(temp_3d(i_delay : i_delay + 1, :));
    temp_3d_no_sound(:, c_delay) = cur_row(2:3:end);
    
    cur_row = mean(temp_3d(i_delay + 2 : i_delay + 3, :));
    temp_3d_collision_sound(:, c_delay) = cur_row(2:3:end);
    
    cur_row = mean(temp_3d(i_delay + 4 : i_delay + 5, :));
    temp_3d_all_sound(:, c_delay) = cur_row(2:3:end);
    
    c_delay = c_delay + 1;
end

% temp_3d_no_sound(1, :) = [];
subplot(1, 3, 1);
boxplot([temp_3d_no_sound(:, 1), ...
    temp_3d_no_sound(:, 2), ...
    temp_3d_no_sound(:, 3), ...
    temp_3d_no_sound(:, 4), ...
    temp_3d_no_sound(:, 5), ...
    temp_3d_no_sound(:, 6), ...
    temp_3d_no_sound(:, 7), ...
    temp_3d_no_sound(:, 8), ...
    temp_3d_no_sound(:, 9)]);
xlabel('Delay (ms)');
ylabel('Rating (%)');
set(gca, 'XTickLabels', {'0','50','100','150','200','250','300','350','400'});
title('no sound');

subplot(1, 3, 2);
boxplot([temp_3d_collision_sound(:, 1), ...
    temp_3d_collision_sound(:, 2), ...
    temp_3d_collision_sound(:, 3), ...
    temp_3d_collision_sound(:, 4), ...
    temp_3d_collision_sound(:, 5), ...
    temp_3d_collision_sound(:, 6), ...
    temp_3d_collision_sound(:, 7), ...
    temp_3d_collision_sound(:, 8), ...
    temp_3d_collision_sound(:, 9)]);
xlabel('Delay (ms)');
ylabel('Rating (%)');
set(gca, 'XTickLabels', {'0','50','100','150','200','250','300','350','400'});
title('with collision sound');

subplot(1, 3, 3);
boxplot([temp_3d_all_sound(:, 1), ...
    temp_3d_all_sound(:, 2), ...
    temp_3d_all_sound(:, 3), ...
    temp_3d_all_sound(:, 4), ...
    temp_3d_all_sound(:, 5), ...
    temp_3d_all_sound(:, 6), ...
    temp_3d_all_sound(:, 7), ...
    temp_3d_all_sound(:, 8), ...
    temp_3d_all_sound(:, 9)]);
xlabel('Delay (ms)');
ylabel('Rating (%)');
set(gca, 'XTickLabels', {'0','50','100','150','200','250','300','350','400'});
title('with collision and rolling sound');
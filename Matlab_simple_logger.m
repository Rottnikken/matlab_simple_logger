delete(instrfindall); clear; clc; close;

%% Connect to device
MyCOM = serial('COM6', 'Baudrate', 115200); % Wired FTDI
fopen(MyCOM);
fprintf(MyCOM, 'k');
disp('COM open');
%% Collect and plot data from device
i = 1;
test = 0;
instring = fscanf(MyCOM,'%s'); % Throw away first string.

% Sets of data collected before finishing program
len = 1000;

time = 1:len;
heading = zeros(len,1);

length_r = zeros(len,1);
length_c = zeros(len,1);

m_x_r_plot = zeros(len,1);
m_y_r_plot = zeros(len,1);
m_z_r_plot = zeros(len,1);

m_x_c_plot = zeros(len,1);
m_y_c_plot = zeros(len,1);
m_z_c_plot = zeros(len,1);

figure('units','normalized','outerposition',[0 0 1 1])
graf = plot(time, m_x_r_plot, time,  m_y_r_plot, time, m_z_r_plot);
axis([0 len -1500 1500]);
grid on;

% subplot(1,2,1);
% hold on;
% theta = linspace(0,2*pi,len);
% x = sin(theta);
% y = cos(theta);
% scat_ref = scatter(x,y,'filled');
% scat = scatter(0,0, 'filled');
% title('Uncalibrated');
% axis([-1 1 -1 1]);
% axis square
% grid on;
% xlabel('m_x');
% ylabel('m_y');

% subplot(1,2,2);
% hold on;
% x = sin(theta);
% y = cos(theta);
% scat2_ref = scatter(x,y,'filled');
% scat2 = scatter(0,0,'filled');
% title('Calibrated');
% axis([-1 1 -1 1]);
% axis square
% grid on;
% xlabel('m_x');
% ylabel('m_y');

%legend('mx\_r', 'mx\_c');
%legend('m_x [Gauss]', 'm_y [Gauss]', 'm_z [Gauss]');
%xlabel('Tid [x100 millisek]');
%% 2D-plot of values
while (i<=len)
    instring = fscanf(MyCOM,'%s');
    M_STRING = strsplit(instring, ',');
    M = zeros(1,6);
    for teller = 1:6
        M(teller) = hex2dec(M_STRING(teller));
        if(M(teller)>32768)
            M(teller) = M(teller) - 65536;
        end
    end
    
    % Uncalibrated x, y and z values
    m_x_r = (M(1));
    m_y_r = (M(2));
    m_z_r = (M(3));
    
    % Calibrated x, y and z values
    m_x_c = (M(4))/1000;
    m_y_c = (M(5))/1000;
    m_z_c = (M(6))/1000;
    
    % Append to datasets
    m_x_r_plot(i) = m_x_r;
    m_y_r_plot(i) = m_y_r;
    m_z_r_plot(i) = m_z_r;
    length_r(i) = sqrt((m_x_r)^2 + (m_y_r)^2);
    radius_r = mean(length_r(1:i));
    
    m_x_c_plot(i) = m_x_c;
    m_y_c_plot(i) = m_y_c;
    m_z_c_plot(i) = m_z_c;
    length_c(i) = sqrt((m_x_c)^2 + (m_y_c)^2);
    radius_c = mean(length_c(1:i));
    
    % Create scatterplots of the data
    %     set(scat, 'XData', m_x_r_plot(1:i));
    %     set(scat, 'YData', m_y_r_plot(1:i));
    %     set(scat_ref, 'XData', radius_r.*x);
    %     set(scat_ref, 'YData', radius_r.*y);
    %
    %     set(scat2, 'XData', m_x_c_plot(1:i));
    %     set(scat2, 'YData', m_y_c_plot(1:i));
    %     set(scat2_ref, 'XData', radius_c.*x);
    %     set(scat2_ref, 'YData', radius_c.*y);
    
    %     set(scat, 'ZData', m_z_plot(1:i));
    
    %     set(graf(1), 'YData', m_x_r_plot(1:i));
    %     set(graf(2), 'YData', m_y_r_plot(1:i));
    %     set(graf(3), 'YData', m_z_r_plot(1:i));
    
    % X-AKSE:
    %     set(graf(1), 'XData', time(1:i));
    %     set(graf(2), 'XData', time(1:i));
    %     set(graf(3), 'XData', time(1:i));
    
    drawnow
    i = i + 1;
end
%% Close connection to device
fprintf(MyCOM, 's');
fclose(MyCOM);
disp('COM closed');
%% Compute min-max values, offset and scaling of sensor axes
clear, close, clc;

% Load datasets into matlab.
load('matlab_mxr_dataset.mat');
load('matlab_myr_dataset.mat');
load('matlab_mzr_dataset.mat');

% Sort arrays in descending order.
mxr_sorted = sort(m_x_r_plot);
myr_sorted = sort(m_y_r_plot);
mzr_sorted = sort(m_z_r_plot);

% Extract 50 first and 50 last in sorted arrays (top and bottom 5% of
% measurements).
mxr_bot = mxr_sorted(1:50);
mxr_top = mxr_sorted(950:1000);
myr_bot = myr_sorted(1:50);
myr_top = myr_sorted(950:1000);
mzr_bot = mzr_sorted(1:50);
mzr_top = mzr_sorted(950:1000);

% Find average min-max value for each sensor axis.
mxr_min = mean(mxr_bot);
mxr_max = mean(mxr_top);
myr_min = mean(myr_bot);
myr_max = mean(myr_top);
mzr_min = mean(mzr_bot);
mzr_max = mean(mzr_top);

% Calculate offset for each axis.
mxr_offset = (mxr_min + mxr_max)/2;
myr_offset = (myr_min + myr_max)/2;
mzr_offset = (mzr_min + mzr_max)/2;

% Calculate offset-adjusteded min-max values for each sensor axis.
mxr_min_adj = mxr_min - mxr_offset;
mxr_max_adj = mxr_max - mxr_offset;
myr_min_adj = myr_min - myr_offset;
myr_max_adj = myr_max - myr_offset;
mzr_min_adj = mzr_min - mzr_offset;
mzr_max_adj = mzr_max - mzr_offset;

% Calculate scaling factor so the resulting values lie within [1, -1].
mxr_scale = (mxr_max_adj - mxr_min_adj)/2;
myr_scale = (myr_max_adj - myr_min_adj)/2;
mzr_scale = (mzr_max_adj - mzr_min_adj)/2;

% Plot datasets and calculated min-max values.
plot(1:1000, m_x_r_plot);
axis([0 1000 -700 700]);
hold on;
grid on;
plot([1 1000], [mxr_max mxr_max] , '--', 'LineWidth', 1.5);
plot([1 1000], [mxr_min mxr_min], '--', 'LineWidth', 1.5);
plot([1 1000], [mxr_offset mxr_offset],  '--', 'LineWidth', 1.5);
legend('mxr', 'mxr\_max', 'mxr\_min', 'mxr\_offset');
xlabel('tid [100ms]');
ylabel('Råverdi fra magnetometer');

figure;
plot(1:1000, m_y_r_plot);
axis([0 1000 -700 700]);
hold on;
grid on;
plot([1 1000], [myr_max myr_max] , '--', 'LineWidth', 1.5);
plot([1 1000], [myr_min myr_min], '--', 'LineWidth', 1.5);
plot([1 1000], [myr_offset myr_offset],  '--', 'LineWidth', 1.5);
legend('myr', 'myr\_max', 'myr\_min', 'myr\_offset');
xlabel('tid [100ms]');
ylabel('Råverdi fra magnetometer');

figure;
plot(1:1000, m_z_r_plot);
axis([0 1000 -700 700]);
hold on;
grid on;
plot([1 1000], [mzr_max mzr_max] , '--', 'LineWidth', 1.5);
plot([1 1000], [mzr_min mzr_min], '--', 'LineWidth', 1.5);
plot([1 1000], [mzr_offset mzr_offset],  '--', 'LineWidth', 1.5);
legend('mzr', 'mzr\_max', 'mzr\_min', 'mzr\_offset');
xlabel('tid [100ms]');
ylabel('Råverdi fra magnetometer');
%%

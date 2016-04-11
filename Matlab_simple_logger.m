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
len = 500;

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
% graf = plot(time, m_x_plot, time,  m_y_plot, time, m_z_plot);

subplot(1,2,1);
hold on;
theta = linspace(0,2*pi,len);
x = sin(theta);
y = cos(theta);
scat_ref = scatter(x,y,'filled');
scat = scatter(0,0, 'filled');
title('Uncalibrated');
axis([-1 1 -1 1]);
axis square
grid on;
xlabel('m_x');
ylabel('m_y');

subplot(1,2,2);
hold on;
x = sin(theta);
y = cos(theta);
scat2_ref = scatter(x,y,'filled');
scat2 = scatter(0,0,'filled');
title('Calibrated');
axis([-1 1 -1 1]);
axis square
grid on;
xlabel('m_x');
ylabel('m_y');

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
    m_x_r = (M(1))/1000;
    m_y_r = (M(2))/1000;
    m_z_r = (M(3))/1000;
    
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
    set(scat, 'XData', m_x_r_plot(1:i));
    set(scat, 'YData', m_y_r_plot(1:i));
    set(scat_ref, 'XData', radius_r.*x);
    set(scat_ref, 'YData', radius_r.*y);
    
    set(scat2, 'XData', m_x_c_plot(1:i));
    set(scat2, 'YData', m_y_c_plot(1:i));
    set(scat2_ref, 'XData', radius_c.*x);
    set(scat2_ref, 'YData', radius_c.*y);
    
    %     set(scat, 'ZData', m_z_plot(1:i));
    
    %     set(graf(1), 'YData', m_x_plot(1:i));
    %     set(graf(2), 'YData', m_y_plot(1:i));
    %     set(graf(3), 'YData', m_z_plot(1:i));
    %
    %     % X-AKSE:
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
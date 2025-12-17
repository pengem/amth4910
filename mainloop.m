% main loop!!
clear all; clc; close all;

% how many times to run simulation
loop = 1;

% storing the time that the first TR starts
firstTR = zeros(loop,1);

% main loop
for Nloop = 1:loop
    test1210;
    firstTR(Nloop) = triggertime;
end


% histograms for trigger time
nameFile = ['histogramfirstTRtest1210_' num2str(loop)];
fig = figure;
histogram(firstTR);
xlabel('Time of first TR (s)', 'FontSize', 14);
ylabel('Frequency', 'FontSize', 14);
title('Histogram of time of first TR', 'FontSize', 20);
grid on;

saveas(fig,[nameFile '.png']);

twodvisualization;
threedvisualization;


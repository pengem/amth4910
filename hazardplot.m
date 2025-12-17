% hazard parameters 
Tcrit = 170; 
lambda0 = 1e-4; 
beta_TR = 0.08;

dt = 6.313131313131314e-04; 

% temperature range
T = 100:0.5:300;

lambda = zeros(size(T));
p_step = zeros(size(T));
p_1s   = zeros(size(T));

for idx = 1:length(T)
    if T(idx) <= Tcrit
        lambda(idx) = 0;
    else
        lambda(idx) = lambda0 * exp(beta_TR * (T(idx) - Tcrit));
    end

    % convert hazard into per-timestep probability
    p_step(idx) = 1 - exp(-lambda(idx) * dt);

    % convert hazard into probability over 1 full second
    p_1s(idx) = 1 - exp(-lambda(idx) * 1.0);
end

% plot hazard rate
figure; 
plot(T, lambda, 'LineWidth', 2);
xlabel('Temperature (°C)', 'FontSize', 14);
ylabel('\lambda(T) (1/s)', 'FontSize', 14);
title('Hazard Rate vs Temperature', 'FontSize', 16);
grid on;

% plot probability per timestep
figure;
plot(T, p_step, 'LineWidth', 2);
xlabel('Temperature (°C)', 'FontSize', 14);
ylabel('P(runaway in one timestep)', 'FontSize', 14);
title(['Runaway Probability per Timestep (dt = ' num2str(dt) ' s)'], 'FontSize', 16);
grid on;

% plot probability per 1 sec
nameFilehazard = ['probper1sec'];
fig = figure;
plot(T, p_1s, 'LineWidth', 2);
xlabel('Temperature (°C)', 'FontSize', 14);
ylabel('P(runaway within 1 second)', 'FontSize', 14);
title('Runaway Probability Within 1 Second', 'FontSize', 20);
grid on;

saveas(fig,[nameFilehazard '.png']);

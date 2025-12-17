tic;
% parameters: # batteries, walls, total 
Nbatt_x = 13; % batteries in x
Nbatt_y = 2; % batteries in y

% double wall grid structure lol
% so the total grid in each direction = (2 wall layers per gap) + 2 (one
% wall border each side)
Nx_core = 3*Nbatt_x + 2;  
Ny_core = 3*Nbatt_y + 2;

% add 1-cell-thick air shell around entire structure
Nx = Nx_core + 2;  
Ny = Ny_core + 2;

% this is adding in height. air, battery, concrete/wood 
Nz_base = 1; % one base layer
Nz_batt = 1; % one battery layer
N_midair = 1; % air layer between cover and battery
N_cover = 1; % the actual cover
N_air_layers = 3; 

Nz_core = Nz_base + Nz_batt + N_midair + N_cover + N_air_layers; % core vertical layers
Nz = Nz_core + 2; % +2 for top and bottom air shell

% layer thicknesses/height
dz_base = 1;
dz_batt = 50;
dz_midair = 0.5;
dz_cover = 2;
dz_air = 10;
dz = [dz_air;
      dz_base*ones(Nz_base,1);
      dz_batt*ones(Nz_batt,1);
      dz_midair*ones(N_midair, 1);
      dz_cover*ones(N_cover, 1);
      dz_air*ones(N_air_layers,1);
      dz_air];

% material prop
alpha_b = 97; % thermal diffustivity of battery
alpha_aerogel = 0.1; % aerogel
alpha_h = 90; % heater... this is a guess lol
alpha_air = 19; % air
alpha_base = 99; % wood/concrete base, change depending on which
alpha_metal = 99;
alpha_cover = 0.96;

% width and height for each of the batteries 
dx_b = 270;
dx_w = 1;

dy_b = 135;
dy_w = 1;

% making each square the correct length: see notes for picture
% first going to make all of the squares the same width and height as
% aerogel

% set up pattern: batteries at indices 3,6,9,... in both x and y (core indices)
batt_x_idx_core = 3:3:Nx_core-2; % all internal battery x positions
batt_y_idx_core = 3:3:Ny_core-2; % all internal battery y positions

% shift to full-domain indices (inside air shell)
batt_x_idx = batt_x_idx_core + 1;
batt_y_idx = batt_y_idx_core + 1;

dx = dx_w*ones(Nx,Ny); 
dy = dy_w*ones(Nx,Ny); 

dx(:, batt_y_idx) = dx_b;
dy(batt_x_idx, :) = dy_b;

% setting metal wall length to different width (now just inside air shell)
dy(2:3, :) = 10;
dy(Nx-2:Nx-1, :) = 10;

% setting air on the sides to different length too (inside air shell)
dx(:, 2:3) = 10;
dx(:, Ny-2:Ny-1) = 10;

% set outer air shell cell sizes (one cell thick)
dx(:,1)  = 10;
dx(:,Ny) = 10;
dy(1,:)  = 10;
dy(Nx,:) = 10;

% setting up the diffusivity of the nodes
alpha = alpha_aerogel*ones(Nx,Ny,Nz); % first setting everything to wall diffusitivty

% assign base, battery, and air layers (all indices shifted by +1 in z)
alpha(:,:,2) = alpha_base; % bottom layer
alpha(:,:,Nz_base+Nz_batt+1+1) = alpha_air; % middle air thing
alpha(:,:,Nz_base+Nz_batt+N_midair+1+1) = alpha_cover;
alpha(batt_x_idx, batt_y_idx, Nz_base+1+1) = alpha_b; % battery layer

% metal walls (now at x=2 and x=Nx-1, inside air shell)
alpha(2, :, :) = alpha_metal;
alpha(Nx-1, :, :) = alpha_metal;
% the next two here will overwrite the left and right rows col air on all
% layers
alpha(:, 2:3, :) = alpha_air;
alpha(:, Ny-2:Ny-1, :) = alpha_air;
alpha(:,:,Nz_base+Nz_batt + N_midair + N_cover +1+1:Nz-1) = alpha_air; % air layers above cover (inside top shell)

% outer air shell: one cell thick in all three directions
alpha(1,:,:) = alpha_air;
alpha(Nx,:,:) = alpha_air;
alpha(:,1,:) = alpha_air;
alpha(:,Ny,:) = alpha_air;
alpha(:,:,1) = alpha_air;
alpha(:,:,Nz) = alpha_air;

% would insert multiple modules here btw

% need to make sure dt abides by the CFL. usually courant number is 0.5, 
% but must be less than or equal to 0.5
Co = 0.25; % this is the courant number
dt = Co * min([dx(:); dy(:); dz(:)]).^2 ./ max(alpha(:));

% Tfinal is how many total seconds we want, and Nt will give # steps
Tfinal = 1500;
Nt = ceil(Tfinal/dt);

% initial temperature field (everything starts at room temperature)
U = 20 * ones(Nx, Ny, Nz);  % °C

% heating stuff for heater 
heat_rate = 2; % °C per second (linear heating rate)

% thermal runaway stats
TR_threshold = 170; % trigger temp for thermal runaway
TR_spike = 1000; % spike temp
TR_hold = 15; % seconds to hold at spike temp
TR_hold_steps = round(TR_hold/dt);  

% hazard parameters
Tcrit = 170; % onset temperature for runaway risk
lambda0 = 1e-4; % base hazard rate at T = Tcrit
beta_TR = 0.08; % exponential steepness

% state array for each battery
state = zeros(Nx,Ny,Nz);
hold_counter = zeros(Nx,Ny,Nz); 

% this is just a way to identify where the heater is, and we make a grid
% where it's all false except for where the heater is. 
heater_mask = false(Nx,Ny,Nz);

% this is the middle battery that we're going to put the heater under. 
chosenbatX = batt_x_idx(ceil(length(batt_x_idx)/2));
chosenbatY = batt_y_idx(ceil(length(batt_y_idx)/2));

% this is the aerogel directly below the middle battery, changing it to a
% heater
heater_x = chosenbatX;
heater_y = chosenbatY + 1; % directly below in y, same z layer
heater_z = Nz_base + 1 + 1; % same z-layer as batteries (shifted by +1 in z)
heater_mask(heater_x, heater_y, heater_z) = true;

% give heater same thermal diffusivity as battery
alpha(heater_x, heater_y, heater_z) = alpha_h;

% allocate array to save Nsave amount of screenshots
Nsave = 1000; 
saveEvery = floor(Nt/Nsave); 
Usave = zeros(Nx, Ny, Nz, Nsave);
saveCount = 1;

% initial save
Usave(:,:,:,saveCount) = U;   

% track elapsed time for heater linear ramp 
heater_elapsed_time = 0;

% main loop: only takes into account the interior cells 
for n = 1:Nt
    Uold = U;

    % going through each of the grid points
    for k = 2:Nz-1
        for j = 2:Ny-1
            for i = 2:Nx-1
                % harmonic mean of diffusivity bc resistances add in series
                alpha_L = 2*alpha(i-1,j,k)*alpha(i,j,k)/(alpha(i-1,j,k)+alpha(i,j,k));
                alpha_R = 2*alpha(i+1,j,k)*alpha(i,j,k)/(alpha(i+1,j,k)+alpha(i,j,k));

                % left/right face spacings use dx at the two adjacent cells
                dx_L = 0.5*(dx(i-1,j) + dx(i,j));
                dx_R = 0.5*(dx(i+1,j) + dx(i,j));

                fluxL = alpha_L/dx_L * (Uold(i-1,j,k) - Uold(i,j,k));
                fluxR = alpha_R/dx_R * (Uold(i+1,j,k) - Uold(i,j,k));

                alpha_D = 2*alpha(i,j-1,k)*alpha(i,j,k)/(alpha(i,j-1,k)+alpha(i,j,k));
                alpha_U = 2*alpha(i,j+1,k)*alpha(i,j,k)/(alpha(i,j+1,k)+alpha(i,j,k));

                % down/up face spacings use dy at the two adjacent cells
                dy_D = 0.5*(dy(i,j-1) + dy(i,j));
                dy_U = 0.5*(dy(i,j+1) + dy(i,j));

                fluxD = alpha_D/dy_D * (Uold(i,j-1,k) - Uold(i,j,k));
                fluxU = alpha_U/dy_U * (Uold(i,j+1,k) - Uold(i,j,k));

                alpha_B = 2*alpha(i,j,k-1)*alpha(i,j,k)/(alpha(i,j,k-1)+alpha(i,j,k));
                alpha_T = 2*alpha(i,j,k+1)*alpha(i,j,k)/(alpha(i,j,k+1)+alpha(i,j,k));

                dz_B = 0.5*(dz(k-1)+dz(k));
                dz_T = 0.5*(dz(k+1)+dz(k));

                fluxB = alpha_B/dz_B * (Uold(i,j,k-1) - Uold(i,j,k));
                fluxT = alpha_T/dz_T * (Uold(i,j,k+1) - Uold(i,j,k));

                % finite-volume form
                U(i,j,k) = Uold(i,j,k) + dt * ( ...
                    (fluxL + fluxR) / dx(i,j) + ...
                    (fluxD + fluxU) / dy(i,j) + ...
                    (fluxB + fluxT) / dz(k) );
            end
        end
    end

    % thermal runaway hazard section
    for j = batt_y_idx
        for i = batt_x_idx
    
            % temperature of this battery cell
            Tcell = U(i,j,Nz_base+1+1);
    
            % only consider TR if this battery hasn't ignited yet
            if state(i,j,Nz_base+1+1) == 0 
    
                if Tcell <= Tcrit
                    lambda = 0;
                else
                    lambda = lambda0 * exp(beta_TR * (Tcell - Tcrit));
                end
    
                % convert hazard rate to per-step probability
                p_step = 1 - exp(-lambda * dt);

                % roll the dice
                if rand() < p_step
                    % trigger TR
                    U(i,j,Nz_base+1+1) = TR_spike;
                    state(i,j,Nz_base+1+1) = 1;
                    hold_counter(i,j,Nz_base+1+1) = TR_hold_steps;
                    triggertime = n*dt; 
                end
    
            % if battery is in runaway, hold the spike temperature
            elseif state(i,j,Nz_base+1+1) == 1
                U(i,j,Nz_base+1+1) = TR_spike;
                hold_counter(i,j,Nz_base+1+1) = hold_counter(i,j,Nz_base+1+1) - 1;
    
                % after hold time, go to "post-runaway" state
                if hold_counter(i,j,Nz_base+1+1) <= 0
                    state(i,j,Nz_base+1+1) = 2;
                end
            end
        end
    end

    % heater stuff: controls the steady increase in temperature 
    if heater_mask(heater_x, heater_y, heater_z)
        if state(chosenbatX, chosenbatY, Nz_base+1+1) == 0
            heater_elapsed_time = heater_elapsed_time + dt;
            U(heater_x, heater_y, heater_z) = 20 + heat_rate * heater_elapsed_time; 
        end
    end

    % store snapshots
    if mod(n, saveEvery) == 0 && saveCount < Nsave
        saveCount = saveCount + 1;
        Usave(:,:,:,saveCount) = U;
    end
end

compute_time = toc;
fprintf('compute time = %.3f seconds\n', compute_time);
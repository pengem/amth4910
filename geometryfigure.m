MAT_AIR     = 1;
MAT_AEROGEL = 2;
MAT_BATTERY = 3;
MAT_HEATER  = 4;
MAT_METAL   = 5;
MAT_BASE    = 6;
MAT_COVER   = 7;


% start: everything aerogel
material_map = MAT_AEROGEL * ones(Nx, Ny);

% battery cells
material_map(batt_x_idx, batt_y_idx) = MAT_BATTERY;

% metal walls (inside air shell)
material_map(2,:)    = MAT_METAL;
material_map(Nx-1,:) = MAT_METAL;

% overwrite with inner air columns
material_map(:,2:3)       = MAT_AIR;
material_map(:,Ny-2:Ny-1) = MAT_AIR;

% outer air shell (one cell thick)
material_map(1,:)  = MAT_AIR;
material_map(Nx,:) = MAT_AIR;
material_map(:,1)  = MAT_AIR;
material_map(:,Ny) = MAT_AIR;

% heater (final overwrite)
material_map(heater_x, heater_y) = MAT_HEATER;


material_colors = [
    0.85 0.90 1.00;   % air
    0.90 0.90 0.90;   % aerogel
    0.20 0.60 0.20;   % battery
    0.90 0.10 0.10;   % heater
    0.50 0.50 0.50;   % metal
    0.60 0.40 0.20;   % base (unused)
    0.20 0.20 0.60;   % cover (unused)
];

% plot
figure;

imagesc(rot90(material_map', 2));
axis equal tight;
set(gca,'YDir','normal');

colormap(material_colors);
caxis([1 size(material_colors,1)]);
title('Battery Layer Schematic');

hold on;

% rotated battery outlines
for bx = batt_x_idx
    for by = batt_y_idx
        bx_r = Nx + 1 - bx;
        by_r = Ny + 1 - by;
        rectangle('Position',[bx_r-0.5, by_r-0.5, 1, 1], ...
                  'EdgeColor','k','LineWidth',1.2);
    end
end

% rotated heater outline
hx_r = Nx + 1 - heater_x;
hy_r = Ny + 1 - heater_y;
rectangle('Position',[hx_r-0.5, hy_r-0.5, 1, 1], ...
          'EdgeColor','r','LineWidth',2);

hold off;

cb = colorbar;
cb.Ticks = 1:7;
cb.TickLabels = {'Air','Aerogel','Battery','Heater','Metal','Base','Cover'};

% 3D visualization code!

S = load('Usave1210loop100.mat');
% 'myVariable' is now a field within the structure S.
Usave = S.Usave;

tic;
% 3d
gifFile3d = '3dgraphs.gif';
%exportgraphics(gcf, gifFile3d);

figure;
colormap parula;

% which cells are batteries
batt_mask = false(Nx, Ny, Nz);
batt_mask(batt_x_idx, batt_y_idx, Nz_base+2) = true; % battery layer

% coordinate grids for plotting
[xg, yg, zg] = ndgrid(1:Nx, 1:Ny, 1:Nz);

saveCount = 1000;

for k = 1:saveCount
    clf;
    Uframe = Usave(:,:,:,k);
    t_elapsed = (k-1)*saveEvery*dt;

    % normalize temperature to [0,1] for color mapping
    cmin = 20; cmax = 1000;
    U_norm = (Uframe - cmin) / (cmax - cmin);
    U_norm = min(max(U_norm, 0), 1);

    % get RGB colors from colormap
    cmap = colormap('parula');
    nC = size(cmap,1);
    colors = cmap(ceil(U_norm * (nC-1) + 1), :);

    % flatten
    colors = reshape(colors, Nx*Ny*Nz, 3);
    batt_flat = batt_mask(:);

    % clear axes
    hold on;
    axis equal tight;
    xlabel('X'); 
    ylabel('Y'); 
    zlabel('Z');
    view(45,30);
    title(sprintf('3D Temperature Field (t = %.2f s)', t_elapsed));
    box on; grid on;

    % plot cubes: batteries solid, others transparent
    for idx = 1:numel(Uframe)
        [ix, iy, iz] = ind2sub(size(Uframe), idx);
        color = colors(idx,:);
        alpha_val = batt_flat(idx)*1.0 + (~batt_flat(idx))*0.15;

        % only draw visible cubes
        if alpha_val > 0.05
            % draw a cube centered at (ix,iy,iz)
            plotcube([1 1 1],[ix-0.5, iy-0.5, iz-0.5], alpha_val, color);
        end
    end

    colorbar;
    caxis([cmin cmax]);
    drawnow;
    if k ==1
        exportgraphics(gcf, gifFile3d);
    end
    exportgraphics(gcf, gifFile3d, Append=true);
end
threedvis_time = toc;
fprintf('3d visualization time = %.3f seconds\n', threedvis_time);
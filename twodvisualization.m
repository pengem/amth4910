% 2D visualization code! 
tic;

figure;
gifFile2d = '2dgraphs.gif';
for k = 1:saveCount
    clf;
    t_elapsed = (k-1) * saveEvery * dt;

    % 2d plot from above
    subplot(1,2,1);
    imagesc(Usave(:,:,Nz_base+1+1,k)'); 
    colorbar; 
    caxis([20 1000]);
    axis equal tight;
    title(sprintf('Battery layer, t = %.2f s', t_elapsed));
    hold on;
    for bx = batt_x_idx
        for by = batt_y_idx
            rectangle('Position',[bx-0.5, by-0.5, 1, 1],'EdgeColor','k','LineWidth',1.2);
        end
    end
    rectangle('Position',[heater_x-0.5, heater_y-0.5, 1, 1],'EdgeColor','r','LineWidth',2);
    hold off;

    % vertical cross section
    subplot(1,2,2);
    sliceY = chosenbatY;
    imagesc(squeeze(Usave(:,sliceY,:,k))');
    set(gca,'YDir','normal');
    colorbar; 
    caxis([20 1000]);
    axis equal tight;
    xlabel('X index'); 
    ylabel('Z layer');
    title('Vertical cross-section (x–z plane)');
    drawnow;
    if k ==1
        exportgraphics(gcf, gifFile2d);
    end
    exportgraphics(gcf, gifFile2d, Append=true);

    
end
twodvis_time = toc;
fprintf('2d visualization time = %.3f seconds\n', twodvis_time);

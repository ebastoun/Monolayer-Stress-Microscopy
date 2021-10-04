function Plot_quiver_variables(x_pix,y_pix,variable1,v1_min,v1_max,quiver1x,quiver1y,...
                                           variable2,v2_min,v2_max,quiver2x,quiver2y,...
                                           path,k,n_images,time,title1,title2,plot_variable)
   
         
     
    h1=figure('visible','off');
    h1.Position=[0 0 2000 2000];
    im(1)  = subplot(1,2,1);
    if plot_variable == 3
        color = jet;
    elseif plot_variable == 4
        color = gray;
    end
    imagesc(x_pix(:),y_pix(:),variable1); colormap(im(1),color);axis off ;axis image;
    title({title1},'FontWeight', 'Normal','fontsize',18);
    if plot_variable == 3
        colorbar;
        caxis([v1_min v1_max]);
    end
    
    
    % PIV plot (vectors and arrows)
    hold on
    if plot_variable == 3
        quiver(x_pix(:),y_pix(:),quiver1x/50,quiver1y/50,'AutoScale','off','color', 'w');
    elseif plot_variable == 4
        quiver(x_pix(:),y_pix(:),quiver1x/100,quiver1y/100,'AutoScale','off','color', 'r');
    end
    hold off
    
    
    
    
    im(2)  = subplot(1,2,2); 
     if plot_variable == 3
        color = flipud(jet);
    elseif plot_variable == 4
        color = gray;
    end
    imagesc(x_pix(:),y_pix(:),variable2); colormap(im(2),color);axis image;
    axis off %xlabel('X position [um]'); ylabel('Y position [um]');
    title({title2},'FontWeight', 'Normal','fontsize',18);
    if plot_variable == 3
        colorbar;
        caxis([v2_min v2_max]);
    end
    
    % PIV plot (vectors and arrows)
    hold on
    if plot_variable == 3
        quiver(x_pix(:),y_pix(:),quiver2x/50,quiver2y/50,'AutoScale','off','color', 'w');
    elseif plot_variable == 4
        quiver(x_pix(:),y_pix(:),quiver2x/100,quiver2y/100,'AutoScale','off','color', 'r');
    end
    hold off
    
    
     caption = sprintf('Frame #%d of %d, t = %d h', k, n_images, time(k));
     sgtitle(caption)
     
     name  = sprintf(path,k);
     Image = getframe(gcf);
     imwrite(Image.cdata, name);

     
     
 end
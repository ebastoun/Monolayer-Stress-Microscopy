function Plot_3_variables(variable1,variable2,variable3,path,k,v1_min,v1_max,v2_min,v2_max,v3_min,v3_max,n_images,time,title1,title2,title3,plot_variable)
    
    
     h1=figure('visible','off');
     h1.Position=[100 100 2000 500];
     ax(1)=subplot(1,3,1);
     imagesc(variable1);colormap(ax(1),jet);axis off;axis image; caxis([v1_min v1_max]);
     title({title1},'FontWeight', 'Normal','fontsize',18);         ax(1).Position = [0.15    0.1100    0.1566    0.8150];
     colorbar;%title(colorbar,'[\mum]','FontSize',18);
     %h = colorbar; h.Ruler.TickLabelFormat='%g%%'
     
     if plot_variable == 1
         
         color = flipud(jet);
         
     else 
         color = jet;
         
     end
     ax(2)=subplot(1,3,2);
     imagesc((variable2));colormap(ax(2),color);axis off;axis image; caxis([v2_min v2_max]);
     title({title2},'FontWeight', 'Normal','fontsize',18);           ax(2).Position = [0.4    0.1100    0.1566    0.8150];
     colorbar; %title(colorbar,'[\mum]','FontSize',18);

     ax(3)=subplot(1,3,3);
     imagesc(variable3); colormap(ax(3),jet);axis off;axis image; caxis([v3_min v3_max]);
     title({title3},'FontWeight', 'Normal','fontsize',18);  ax(3).Position = [0.66    0.1100    0.1566    0.8150];
     colorbar;%title(colorbar,'[\mum]','FontSize',18);
     
     caption = sprintf('Frame #%d of %d, t = %d h', k, n_images, time(k));
     sgtitle(caption)
     
     name  = sprintf(path,k);
     Image = getframe(gcf);
     imwrite(Image.cdata, name);

     
     
 end
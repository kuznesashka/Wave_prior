function click = PlotMri(MRI, pos, h)

    k = str2num(get(gca, 'tag'));
    click = get(gca,'CurrentPoint');
        
    if k == 1
        clx = round(click(1,1)); % sagittal
        cly = pos(2);        % coronal
        clz = round(click(1,2)); % axial
        
        subplot(2,3,1)
        imagesc(flipud(squeeze(MRI.Cube(:,cly,:))'))
        set(gca,'tag',num2str(1))
        hold on
        plot(clx,clz,'r+', 'MarkerSize', 30);
        
        subplot(2,3,2)
        imagesc(flipud(squeeze(MRI.Cube(256-clx,:,:))'))
        set(gca,'tag',num2str(2))
        hold on
        plot(pos(2),clz,'r+', 'MarkerSize', 30);
        
        subplot(2,3,3)
        imagesc(flipud(squeeze(MRI.Cube(:,:,256-clz))'))
        set(gca,'tag',num2str(3))
        hold on
        plot(clx,256-pos(2),'r+', 'MarkerSize', 30);
        pos([1,3]) = [clx, clz];
    else
        if k == 2
            clx = pos(1); % sagittal
            cly = round(click(1,1));        % coronal
            clz = round(click(1,2)); % axial
            
            subplot(2,3,2)
            imagesc(flipud(squeeze(MRI.Cube(256-pos(1),:,:))'))
            set(gca,'tag',num2str(2))
            hold on
            plot(cly,clz,'r+', 'MarkerSize', 30);
            
            subplot(2,3,1)
            imagesc(flipud(squeeze(MRI.Cube(:,cly,:))'))
            set(gca,'tag',num2str(1))
            hold on
            plot(256-pos(1),clz,'r+', 'MarkerSize', 30);
            
            subplot(2,3,3)
            imagesc(flipud(squeeze(MRI.Cube(:,:,256-clz))'))
            set(gca,'tag',num2str(3))
            hold on
            plot(pos(1),256-cly,'r+', 'MarkerSize', 30);
            pos([2,3]) = [cly, clz];
            
        else
            clx = round(click(1,1)); % sagittal
            cly = round(click(1,2));        % coronal
            clz = pos(3); % axial
            
            subplot(2,3,3)
            imagesc(flipud(squeeze(MRI.Cube(:,:,256-clz))'))
            set(gca,'tag',num2str(3))
            hold on
            plot(clx,cly,'r+', 'MarkerSize', 30);
            
            subplot(2,3,1)
            imagesc(flipud(squeeze(MRI.Cube(:,256-cly,:))'))
            set(gca,'tag',num2str(1))
            hold on
            plot(clx,clz,'r+', 'MarkerSize', 30);
            
            subplot(2,3,2)
            imagesc(flipud(squeeze(MRI.Cube(clx,:,:))'))
            set(gca,'tag',num2str(2))
            hold on
            plot(256-cly,clz,'r+', 'MarkerSize', 30);
            pos([1,2]) = [clx, 256-cly];
        end
    end

    setappdata(h, 'pos', pos)
 
end
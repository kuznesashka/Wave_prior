 function scrolling_plot2(x, num_val, events, events2)
 
 % x -- [1 x T] values on the x axis
 % Y -- [ncomps x T] matrix of values on the y axis
 % num_val -- number of values shown at one time, step between the two
 % figures is 0.01*T
 % events -- [1 x T] indicator of the event
 % nw -- number of components shown on the one plot, ncomps::nw
 
 
    % Create the first figure
    f = figure('Visible','off');
    % pick the first nw components
    subplot(2,1,1)
    plot(x(1:(1+num_val)), events(1:(1+num_val)),'LineWidth', 1.5)
    xlim([x(1) x(1+num_val)])
    ylim([0 1])
    subplot(2,1,2)
    plot(x(1:(1+num_val)), events2(1:(1+num_val)),'LineWidth', 1.5)
    xlim([x(1) x(1+num_val)])
    ylim([0 1])
    
        
    % Create slide slide the time axis
    sld = uicontrol('Style', 'slider',...
        'Min',1,'Max',size(x,2)-num_val-mod((size(x,2)-num_val), 100),...
        'Value',1,...
        'Position', [400 5 120 20],...
        'SliderStep', [0.01 1], ... % step = 0.01*T
        'Callback', @slidedata);
    
      
    % Make figure visble after adding all components
    f.Visible = 'on';
    
    function slidedata(source,event)
        val = int32(source.Value);
        subplot(2,1,1)
        plot(x(val:(val+num_val)), events(val:(val+num_val)), 'LineWidth', 1.5)
        xlim([x(val) x(val+num_val)])
        ylim([0 1])
        subplot(2,1,2)
        plot(x(val:(val+num_val)), events2(val:(val+num_val)), 'LineWidth', 1.5)
        xlim([x(val) x(val+num_val)])
        ylim([0 1])
     end

   
end
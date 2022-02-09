function  plot_time(y,fs,title_label)
    t = linspace(0, length(y)/fs, length(y));
    figure();
    plot(t,y);
    title(title_label);
    xlabel('Time'); 
    ylabel('Value'); 
    grid;
end
function  plot_frequency(y,fs,title_label)
    f = linspace(-fs/2,fs/2,length(y));
    figure()
    plot(f,y);
    title(title_label);
    xlabel('Frequency (HZ)'); 
    ylabel('Amplitude'); 
    grid;
end
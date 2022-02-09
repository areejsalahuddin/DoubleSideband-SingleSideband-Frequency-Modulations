function filter = generate_filter(signal_length,fs)
    filter = ones(signal_length,1);
    f = linspace(-fs/2,fs/2,signal_length);
    for i = 1: signal_length
        if abs(f(i))>4000
            filter(i)=0;
        end
    end
end
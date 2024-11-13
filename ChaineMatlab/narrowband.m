function noise_add = narrowband(SNR, cf, n, bw)
    low_f = cf/2^(bw/2);
    high_f= cf*2^(bw/2);
    fs = 4096;										%sampling rate 
    x = randn(n,1);
    y = bandpass(x, [low_f high_f], fs);						%creating bandpass filter 
    
    SNR = 10;										%fix the power of noise 
    signal_power = rms(y)^2;
    noise_power_ratio = 10^(SNR/10);
    noise_power = signal_power*noise_power_ratio;					%noise power 
    noise_add = sqrt(noise_power)*randn(n, 1);
    noisy_signal = y + noise_add;							%adding white noise to bp filter
    noise_add = 1.5*noisy_signal;
end

function signal_SNR = SNR_estimation(Signal,frequency,time_shift,CALen,sat_number)
useful_Signal = Signal(time_shift:time_shift+CALen);
CACode = GenCACode(1, sat_number);
P_noise_avg = mean(abs((useful_Signal-CACode)).^2);   
P_signal_avg = mean(abs((CACode)).^2)
signal_SNR = 10*log10(P_signal_avg/P_noise_avg);
clc;
clear;
close all;
%строит вероятность правильного обнаружения

%1 - когерентный; 2 - некогерентный с корнем 20 накоплений; 3 - некогерентный без корня;
%4 некогерентный с корнем 10 накоплений
Method =4;
switch Method
    case 1
        number_of_pieces = 2;
        NumCA2Search = 10;
        Threshold = 4.5; % False alarm = 10^-7 
    case 2
        number_of_pieces = 1;
        NumCA2Search = 20;
        Threshold = 1.64; %3.383713061132416e-07
    case 3
        number_of_pieces = 1;
        NumCA2Search = 20;
        Threshold = 2.48;
    case 4
        number_of_pieces = 1;
        NumCA2Search = 10;
        
       Threshold = 1.64; %3.383713061132416e-07 1.93
end
number_of_experiment = 1000;
SNR = -38:1:-27;

CALen = 2046;
CentralFreqs = -6000 : 1000 : 6000;
NumCFreqs = length(CentralFreqs);


dt=10^-3/CALen;


dt=10^-3/CALen;
CACode = GenCACode(1, 1);
CACode(CACode==0)=-1;
CACode(2,:) = CACode;
CACode = reshape(CACode,1,2*size(CACode,2));
n = 0:number_of_pieces*CALen*NumCA2Search+CALen-2;
number_false_i = zeros(1,length(Threshold));

number_right_i = zeros(1,length(SNR));

for i = 1:number_of_experiment
    
    Noise = [];
    Noise = [1, 1j]*randn(2, number_of_pieces*CALen*NumCA2Search+CALen-1);
    iii = 1;

    for snr = SNR
        correlation=[];
        Signal_shift = [];
        coherent_sum=[];
        norm_coherent_sum=[];
        Signal = [];
        Signal = repmat(CACode,1,number_of_pieces*NumCA2Search);
        Signal = [Signal,CACode(1:CALen-1)]; % добавляем в конец, чтобы посчитать conv valid до конца
        S = mean(Signal.^2);
        N = mean (abs(Noise).^2);
        K = sqrt((10^(snr/10)*N)/S);
        Signal = K.*Signal;
        S_new =  mean(Signal.^2);
%         SNR_estimated = 10*log10(S_new/N);
        Signal = Signal + Noise;
        
        
        for k = 1:NumCFreqs
            Signal_shift(k,:) =  Signal.*exp(-n*1j*2*pi*CentralFreqs(k)*dt);
            correlation(k,:) = conv(Signal_shift(k,:),fliplr(CACode),'valid');
            for ii = 0:number_of_pieces-1
                if NumCA2Search~=1
                    switch Method
                        case 1
                            summa(k,ii*CALen+(1:CALen)) = abs(sum(reshape(correlation(k,ii*CALen*NumCA2Search+(1:CALen*NumCA2Search)),CALen,NumCA2Search)'));
                        case 2
                            summa(k,ii*CALen+(1:CALen)) = sum(abs(reshape(correlation(k,ii*CALen*NumCA2Search+(1:CALen*NumCA2Search)),CALen,NumCA2Search)'));
                        case 3
                            ImValues = [];
                            SqrImValue = [];
                            ImValues = reshape(correlation(k,ii*CALen*NumCA2Search+(1:CALen*NumCA2Search)),CALen,NumCA2Search);
                            SqrImValue = ImValues.*conj(ImValues);
                            summa(k,ii*CALen+(1:CALen)) = sum(SqrImValue');
                        case 4
                            summa(k,ii*CALen+(1:CALen)) = sum(abs(reshape(correlation(k,ii*CALen*NumCA2Search+(1:CALen*NumCA2Search)),CALen,NumCA2Search)'));
                    end
                else
                    summa(k,:) = abs(correlation(k,:));
                end
            end
        end

        for ii = 0:number_of_pieces-1
            norm_summa(:,ii*CALen+(1:CALen)) = summa(:,ii*CALen+(1:CALen))/mean(mean(summa(:,ii*CALen+(1:CALen))));
        end
        
    max_norm_summa = max(max(norm_summa));
    S_new2 =  mean(max_norm_summa.^2);

    if  max_norm_summa >  Threshold
        [frequency_point,sample_point] = find(norm_summa == max_norm_summa);
        sample_shift = mod(sample_point,CALen)-1;
        if (CentralFreqs(frequency_point) == 0)&&(sample_shift==0)
            number_right_i(iii) = number_right_i(iii) + 1;
        end
    end
        
        iii = iii+1;
    end

    disp (100*i/number_of_experiment);

end
right_recieving_probability = n /(number_of_experiment);

figure(1);
plot(SNR,right_recieving_probability,'-o');
grid on;
xlabel ('SNR');
ylabel ('Вероятность правильного обнаружения');



% right_recieving_probability3= right_recieving_probability;
% save('C:\Users\1\Desktop\СНС\Matlab\Figures\right_curve3.mat','right_recieving_probability3');
% load('C:\Users\1\Desktop\СНС\Matlab\Figures\right_curve2.mat');
% 
% figure(2);
% plot(SNR,right_recieving_probability1,'-o');
% hold on;
% plot(SNR,right_recieving_probability2,'-o');
% hold on;
% plot(SNR,right_recieving_probability3,'-o');
% hold on;
% plot(SNR,right_recieving_probability4,'-o');
% grid on;
% xlabel ('SNR');
% ylabel ('Вероятность правильного обнаружения');
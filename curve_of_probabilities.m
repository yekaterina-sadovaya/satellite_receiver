%строит вероятность ложной тревоги
clc;
clear;
close all;
%1 - когерентный; 2 - некогерентный с корнем 20 накоплений; 3 - некогерентный без корня;
%4 некогерентный с корнем 10 накоплений
Method = 4;

switch Method
    case 1
        number_of_pieces = 2; 
        NumCA2Search = 10;
    case 2
        number_of_pieces = 1;
        NumCA2Search = 20;
    case 3
        number_of_pieces = 1;
        NumCA2Search = 20;
    case 4
        number_of_pieces = 1;
        NumCA2Search = 10;
end
number_of_experiment = 100;
CALen = 2046;

CentralFreqs = -6000 : 1000 : 6000;
NumCFreqs = length(CentralFreqs);
Treshold = 1:0.01:7;

dt=10^-3/CALen;
CACode = GenCACode(1, 1);
CACode(CACode==0)=-1;
CACode(2,:) = CACode;
CACode = reshape(CACode,1,2*size(CACode,2));
n = 0:number_of_pieces*CALen*NumCA2Search+CALen-2;
number_false_i = zeros(1,length(Treshold));

for i = 1:number_of_experiment
    correlation=[];
    Signal_shift = [];
    coherent_sum=[];
    norm_coherent_sum=[];
    N = [];
    N = [1, 1j]*randn(2, number_of_pieces*CALen*NumCA2Search+CALen-1);
    
    for k = 1:NumCFreqs
        Signal_shift(k,:) = N.*exp(-n*1j*2*pi*CentralFreqs(k)*dt);
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
    iii = 1;
    for Tr = Treshold
        number_false_i(iii) = number_false_i(iii) + length(find(norm_summa>Tr));
        iii = iii+1;
    end
end
false_alarm_probability = number_false_i./(number_of_experiment*length(CentralFreqs)*CALen);
figure(1);
semilogy(Treshold,false_alarm_probability,'-o');
grid on;
xlabel ('Порог');
ylabel ('Вероятность ложного обнаружения');

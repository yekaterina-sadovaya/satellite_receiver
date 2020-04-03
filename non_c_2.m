clc;
clear all;
close all;

load('C:\Users\1\Desktop\11_semester\GPS\ямя\MATLAB\Results_Rate2_Coh\Rate2.mat');
Res.File.Name =  'C:\Users\1\Desktop\ямя\Matlab\Signals\30_08_2018__19_38_33_x02_1ch_16b_15pos_90000ms.dat';

CALen = 2046;

SatNum = 1;

    CACode = GenCACode(Res.Search.SatNums(SatNum),1);
    CACode(CACode==0)=-1;
    CACode = repmat( CACode,Res.File.R,1);
    CACode = reshape(CACode,1,[]);
    dt=10^-3/CALen;
    CACode = CACode.*exp((0:CALen-1)*1j*2*pi*Res.Search.FreqShifts(SatNum)*dt);
    
    Shift = Res.Search.SamplesShifts(SatNum);
    
    Bn = 10;
    W = 4*Bn;
    Z = 0;
    T = 1e-3;
    
    NCA = 20000;
    
    CAShifts = zeros(1, NCA);
    NCOVals =  zeros(1, NCA);
    Prompts =  zeros(1, NCA);
    Discrs  =  zeros(1, NCA);
    
    for k = 1:NCA
        [Signal, ~] = ReadSignalFromFile(Res.File, Shift-1, CALen+2);
            
        EPL = abs(conv(Signal, fliplr(CACode), 'valid'));
        
        Prompts(k) = EPL(2);
        Discrs(k) = 0.5*(EPL(1) - EPL(3))/(EPL(1) + EPL(3));
        if k == 1
            NCOVals(k) = Discrs(k)*W*T;
        else
            NCOVals(k) = NCOVals(k-1) + Discrs(k)*W*T;
        end
        
        Shift = Shift + CALen;
        if NCOVals(k) > 0.5
            Shift = Shift - 1;
            NCOVals(k) = NCOVals(k) - 1;
        elseif NCOVals(k) < -0.5
            Shift = Shift + 1;
            NCOVals(k) = NCOVals(k) + 1;
        end
    end
    
    figure();
    plot(Prompts);
    figure();
    plot(NCOVals);
       figure();
    plot(Discrs); 
% %     figure();
%     plot(Track.CorVals{i});
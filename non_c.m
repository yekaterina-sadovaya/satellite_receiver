clc;
clear all;
close all;
load('C:\Users\yulya\Desktop\СНС\Matlab\Results_Rate2_Coh\Rate2.mat');
Res.File.Name =  'C:\Users\yulya\Desktop\СНС\Matlab\Signals\30_08_2018__19_38_33_x02_1ch_16b_15pos_90000ms.dat';
Track = struct('a',[]);
 k = 1;
 i = 1;
  CALen = 2046;
 dt=10^-3/CALen;

NumCA2NextSync =1;
   shift_cor_vals = [-1;1]
    CACode = GenCACode(Res.Search.SatNums(k),1);
    CACode(find(CACode==0))=-1;
    CACode = repmat( CACode,Res.File.R,1);
    %     CACode = CACode*(1+1j);
    CACode = reshape(CACode,1,Res.File.R*size(CACode,2));
    CACode = CACode.*exp((0:CALen-1)*1j*2*pi*Res.Search.FreqShifts(k)*dt);
    read_group_number = 0;
    adjustment_nacop = 0;
    while (read_group_number+1) <= 2000
            Signal_check_sum = [];
            %У меня записывается шифт для каждой грцппы из 100 периодов или
            %столько сколько укажем, считывается по 100 и сразу
            %декодируется. В промежутке шифт не находится
            Track.SamplesShifts{i}(read_group_number+1) = Res.Search.SamplesShifts(k) + read_group_number*CALen + adjustment_nacop;
            [Signal, ~] = ReadSignalFromFile(Res.File, Track.SamplesShifts{i}(read_group_number+1), CALen);
            Track.CorVals{i}(read_group_number+1) = abs(CACode*Signal');
            ii=1;
for i2 = 1:2
            sh = shift_cor_vals(i2);
                Signal_check=[];
                
                [Signal_check(ii,:), ~] = ReadSignalFromFile(Res.File, Track.SamplesShifts{i}(read_group_number+1)+sh,CALen);
                promt(ii)=abs(CACode*Signal_check(ii,:)');
                ii = ii+1;
            end
if 1/2*(promt(1)-promt(2))/(promt(1)+promt(2))>0.5
    Signal = [];
    adjustment_nacop = adjustment_nacop + 1;
         Track.CorVals{i}(read_group_number+1) = promt(2);
elseif 1/2*(promt(1)-promt(2))/(promt(1)+promt(2))<0.5
    Signal = [];
            adjustment_nacop = adjustment_nacop - 1;
                 Track.CorVals{i}(read_group_number+1) = promt(1);
end
            cum_curve(read_group_number+1) = adjustment_nacop;


    Track.SamplesShifts{i}(read_group_number+1) = Res.Search.SamplesShifts(k) + adjustment_nacop;


  
       

               read_group_number = read_group_number+1;
    end
    figure();
    plot(cum_curve);
    figure();
    plot(Track.CorVals{i});
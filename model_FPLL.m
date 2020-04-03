clc;
clear all;
close all;
DataSize = 10^4;
% DataSize = 100;
NumCaperBit = 20;
Data = randi([0,1],1,DataSize);
Data(Data==0)=-1;
CAs = reshape(repmat(Data,NumCaperBit,1),1,NumCaperBit*DataSize);


dt = 1*10^-3;
L  = length(CAs);
NumberSummCAs = 20;
% % % for Phi0 = [-pi:0.5*pi:pi]
% % %     for f_0 = [-5*10^3:1*10^3:5*10^3]
% % %         for df =  [-1:0.5:1]
for Phi0 =0
    for f_0 = 3*10^3
        for df =  0
            f  = [];
            PCAs = [];
            Phi = [];
            dPhiValues = [];
            ddPhiValues = [];
            f = f_0 + df*(0:L-1)*dt;
            Phi = Phi0 + 2*pi*dt.*cumsum(f);
            PCAs = CAs.*exp(1j*Phi);
            %трэкинг после интегратора когерентного
            numberGroup = 1;
            EstimatedPhi = 0;
            EstimatedF = 0;
            EstimatedFs(1) = 0;
            EstimatedPhi0 = 0;
            F = 0;
            Phi2 = 0;
            Phi1  = 0;
            PCAValue1 = 0;
            FinalPhi = 0;
            Phi1_nacop = 0;
            while numberGroup*NumberSummCAs<=length(PCAs)
% % %                 PhiValues ((numberGroup-1)*NumberSummCAs+(1:NumberSummCAs)) = PCAs((numberGroup-1)*NumberSummCAs+(1:NumberSummCAs)).*exp(-1j*(EstimatedPhi0+2*pi*dt*F));
                PhiValues ((numberGroup-1)*NumberSummCAs+(1:NumberSummCAs)) = PCAs((numberGroup-1)*NumberSummCAs+(1:NumberSummCAs)).*exp(1j*(- Phi1_nacop));
                PCAValue(numberGroup) = sum(PhiValues((numberGroup-1)*NumberSummCAs+(1:NumberSummCAs)));
                
                %оценка частоты
                
                PCAValue2 = PCAValue(numberGroup);
                %вектор частот всегда будет на один меньше, чем кол-во
                %групп СА  кодов посчитанных, потому что нужно смотреть
                %разницу между соседними, и поэтому брать 2 значения
                %всегда, поэтому будет на1 меньше
                % % %                 EstimatedF  = real(PCAValue1) * imag(PCAValue2) - real(PCAValue2) * imag(PCAValue1) ;
%                 diff = real(PCAValue1) * imag(PCAValue2) - real(PCAValue2) * imag(PCAValue1);
                diff = atan2(imag(PCAValue2),real(PCAValue2)) - atan2(imag(PCAValue1),real(PCAValue1));
%                 if abs(mod(diff,pi)) > pi/2
%                     diff = -(pi - mod(diff,pi));
%                 else
%                     diff  =  mod(diff,pi);
%                 end
                EstimatedF = (diff)/(2*pi*dt);

                Phi1 = 2*pi* EstimatedF*dt;
                if numberGroup ~=1
                Phi1_nacop = Phi1_nacop +  diff;
                end
                EstimatedFs (numberGroup)  = EstimatedF;
                
                % оценка фазы
                EstimatedPhi = atan(imag(PCAValue(numberGroup))/real(PCAValue(numberGroup)));
% % %                 if abs(mod(EstimatedPhi,pi)) > pi/2
% % %                     EstimatedPhi = pi - mod(EstimatedPhi,pi);
% % %                 else
% % %                     EstimatedPhi  = - mod(EstimatedPhi,pi);
% % %                 end
                EstimatedPhis (numberGroup) = EstimatedPhi; 
                Phi2 = EstimatedPhi;
                EstimatedPhi0 = EstimatedPhi0+EstimatedPhi;
                
%                 сравнение
                if numberGroup ~=1
                if Phi1 > Phi2
                    FinalPhi = FinalPhi + Phi1;
                else
                    FinalPhi = FinalPhi + Phi1;
                end
                end
               
                PCAValue1 =  PCAValue2;

                
                numberGroup = numberGroup +1;
            end
            %делаем 2 раза деференцирование, первый раз чтобы
            %линейное отклонение устранить, второе чтобы квадратичное,
            %которе на большом масштабе времени
            dPhiValues = conj(PhiValues(1:end-1)).*PhiValues(2:end);
            ddPhiValues = conj(dPhiValues(1:end-1)).*dPhiValues(2:end);
            figure();
            subplot(3,1,1);
            plot(angle(CAs)/pi,'.','color','blue');
            hold on;
            plot(angle(PCAs)/pi,'.','color','red');
            hold on;
            plot(angle(PhiValues)/pi,'.','color','green');
            subplot(3,1,2);
            plot(angle(dPhiValues)/pi,'.');
            subplot(3,1,3);
            plot(angle(ddPhiValues)/pi,'.');
            figure();
            plot(EstimatedPhis,'.')
            figure();
            plot(EstimatedFs,'.')
%             figure();
%             plot(Phi1_nacop,'.')
            a=1;
        end
    end
end



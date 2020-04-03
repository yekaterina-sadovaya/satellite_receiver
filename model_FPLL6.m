clc;
clear all;
close all;
DataSize = 10^3;
% DataSize = 100;
NumCaperBit = 20;
 Data = randi([0,1],1,DataSize);
%  Data = ones(1,DataSize);
Data(Data==0)=-1;
CAs = reshape(repmat(Data,NumCaperBit,1),1,NumCaperBit*DataSize);


dt = 1*10^-3;
L  = length(CAs);
NumberSummCAs = 1;
% % % for Phi0 = [-pi:0.5*pi:pi]
% % %     for f_0 = [-5*10^3:1*10^3:5*10^3]
% % %         for df =  [-1:0.5:1]
for Phi0 =pi/4
    for f_0 = 9
        for df = 0.1
            f  = [];
            PCAs = [];
            Phi = [];
            dPhiValues = [];
            ddPhiValues = [];
            f = f_0 + df*(0:L-1)*dt;
            Phi = Phi0 + 2*pi*dt.*(cumsum(f));
            PCAs = CAs.*exp(1j*Phi);
            %трэкинг после интегратора когерентного
            numberGroup = 1;
            EstimatedPhi = 0;
            EstimatedF = 0;
            EstimatedFs(1) = 0;
            EstimatedPhi0 = 0;
            OutNCO = 0;
            PCAValue1 = 0;
            InFilterFreq=0;
            Acc = 0;
            VelocAcc = 0;
            AccelAcc = 0;
            while numberGroup<=length(PCAs)
                % % %                 PhiValues ((numberGroup-1)*NumberSummCAs+(1:NumberSummCAs)) = PCAs((numberGroup-1)*NumberSummCAs+(1:NumberSummCAs)).*exp(-1j*(EstimatedPhi0+2*pi*dt*F));
                %                 PhiValues ((numberGroup-1)*NumberSummCAs+(1:NumberSummCAs)) = PCAs((numberGroup-1)*NumberSummCAs+(1:NumberSummCAs)).*exp(1j*(-OutNCO));
                %                 PCAValue(numberGroup) = sum(PhiValues((numberGroup-1)*NumberSummCAs+(1:NumberSummCAs)));
                PCAValue = [];
                InFilterPhase = 0;
                LoopFilterFPLL = ClassFilter;
                Order = [2,3];
                Bn = [10,10];
                
                LoopFilterFPLL.PrepareFilter(Order, Bn, dt, VelocAcc, AccelAcc);
                
                
                PCAValue = PCAs(numberGroup).*exp(1j*(-OutNCO));
                PCAs_correct (numberGroup)= PCAValue;
                PCAValues(1) = PCAValue;
                if length(PCAValues)>1
                    PCAValue1 = PCAValues(2);
                    PCAValue2 = PCAValues(1);
                      InFilterFreq = (atan(imag(PCAValue2)/real(PCAValue2)) - atan(imag(PCAValue1)/real(PCAValue1)))/dt;
%                    InFilterFreq = atan2((real(PCAValue1) * real(PCAValue2) + imag(PCAValue2) * imag(PCAValue1)),(real(PCAValue1) * imag(PCAValue2) - real(PCAValue2) * imag(PCAValue1)))/dt;
                    EstimatedF = InFilterFreq/(2*pi);
                    EstimatedFs (numberGroup)  = EstimatedF;
                    
                    EstimatedPhi = atan(imag(PCAValue2)/real(PCAValue2));
                    EstimatedPhis (numberGroup) = EstimatedPhi;
                    InFilterPhase = EstimatedPhi;
                    
                    [Output(numberGroup), VelocAcc, AccelAcc] = LoopFilterFPLL.Step(InFilterFreq,InFilterPhase);
                    OutNCO =  0.5*dt*Output(numberGroup)+Acc;
                    Acc = Acc + dt*Output(numberGroup);
                    Veloc(numberGroup) = VelocAcc;
                    Accel(numberGroup) = AccelAcc;
                end
                PCAValues(2) = PCAValue(1);
                
                numberGroup = numberGroup +1;
            end
            %делаем 2 раза деференцирование, первый раз чтобы
            %линейное отклонение устранить, второе чтобы квадратичное,
            %которе на большом масштабе времени

            figure();
            subplot(3,1,1);
            plot(angle(CAs)/pi,'.','color','blue');
            subplot(3,1,2);
            plot(angle(PCAs)/pi,'.','color','red');
            subplot(3,1,3);
            plot(angle(PCAs_correct)/pi,'.','color','green');

            figure();
            plot(EstimatedPhis,'.')
            figure();
            plot(EstimatedFs,'.','color','red');
            figure ();
            subplot(2,1,1);
            plot(Accel/(2*pi),'.');
            title('AccelAcc');
            subplot(2,1,2);
            plot(Veloc/(2*pi),'.');
            title('VelocAcc');

%             figure();
%             plot(Phi1_nacop,'.')
        end
    end
end



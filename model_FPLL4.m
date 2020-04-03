clc;
clear all;
close all;
DataSize = 10^3;
% DataSize = 100;
NumCaperBit = 20;
% Data = randi([0,1],1,DataSize);
Data = ones(1,DataSize);
Data(Data==0)=-1;
CAs = reshape(repmat(Data,NumCaperBit,1),1,NumCaperBit*DataSize);


dt = 1*10^-3;
L  = length(CAs);
NumberSummCAs = 1;
% % % for Phi0 = [-pi:0.5*pi:pi]
% % %     for f_0 = [-5*10^3:1*10^3:5*10^3]
% % %         for df =  [-1:0.5:1]
for Phi0 =pi/4
    for f_0 =125
        for df = 0
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
numCA =2;
            InFilterFreq=0;
Acc = 0;
            while numberGroup*numCA<=length(PCAs)
                PCAValues = [];
                PCAValues = PCAs(numCA*(numberGroup-1)+(1:numCA)).*exp(1j*(-OutNCO));

                PCAs_correct (numCA*(numberGroup-1)+(1:numCA))= PCAValues;
                PCAValue1 =  PCAValues(1);
                PCAValue2 = PCAValues(2);
                PCAValue = PCAValue1+PCAValue2;
                
                
%                 PCAValue = [];
%                 PCAValue = PCAs(numberGroup).*exp(1j*(-OutNCO));
%                 PCAs_correct (numberGroup)= PCAValue;
%                 PCAValue1 =  PCAValue2;
%                 PCAValue2 = PCAValue;


%     InFilterFreq = (atan2(imag(PCAValue2),real(PCAValue2)) - atan2(imag(PCAValue1),real(PCAValue1)));
% if numberGroup   ~= 1
         InFilterFreq = atan2(real(PCAValue1) * real(PCAValue2) + imag(PCAValue2) * imag(PCAValue1),real(PCAValue1) * imag(PCAValue2) - real(PCAValue2) * imag(PCAValue1))/dt;

                EstimatedF = InFilterFreq/(dt);
                EstimatedFs (numberGroup)  = EstimatedF;      

                % оценка фазы
                EstimatedPhi = atan(imag(PCAValue)/real(PCAValue));

                EstimatedPhis (numberGroup) = EstimatedPhi; 

                InFilterPhase = EstimatedPhi;
                LoopFilterFPLL = ClassFilter;
                Order = [1,2];
                Bn = [5,5];
                
                LoopFilterFPLL.PrepareFilter(Order, Bn, dt);

                
                [Output(numberGroup), VelocAcc(numberGroup), AccelAcc(numberGroup)] = LoopFilterFPLL.Step(InFilterFreq, InFilterPhase);
                
                OutNCO =  0.5*dt*Output(numberGroup)+Acc;
                Acc = Acc + dt*Output(numberGroup);
%                                 OutNCO =  dt*Output(numberGroup)+OutNCO;
% end

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
            subplot(3,1,1);
            plot(EstimatedFs,'.');
            subplot(3,1,2);
            plot(AccelAcc/(2*pi),'.');
            subplot(3,1,3);
            plot(VelocAcc/(2*pi),'.');

%             figure();
%             plot(Phi1_nacop,'.')
        end
    end
end



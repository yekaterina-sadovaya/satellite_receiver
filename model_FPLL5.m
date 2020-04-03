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
    for f_0 =50
        for df = 0
            f  = [];
            PCAs = [];
            Phi = [];
            dPhiValues = [];
            ddPhiValues = [];
            f = f_0 + df*(0:L-1)*dt;
            Phi = Phi0 + 2*pi*dt.*(cumsum(f));
            PCAs = CAs.*exp(1j*Phi);
            %������� ����� ����������� ������������
            numberGroup = 1;
            EstimatedPhi = 0;
            EstimatedF = 0;
            EstimatedFs(1) = 0;
            EstimatedPhi0 = 0;
            OutNCO = 0;
            PCAValue1 = 0;
            
            FLAG_FLL = 0;
            InFilterFreq=0;
            Acc = 0;
            while numberGroup<=length(PCAs)
                % % %                 PhiValues ((numberGroup-1)*NumberSummCAs+(1:NumberSummCAs)) = PCAs((numberGroup-1)*NumberSummCAs+(1:NumberSummCAs)).*exp(-1j*(EstimatedPhi0+2*pi*dt*F));
                %                 PhiValues ((numberGroup-1)*NumberSummCAs+(1:NumberSummCAs)) = PCAs((numberGroup-1)*NumberSummCAs+(1:NumberSummCAs)).*exp(1j*(-OutNCO));
                %                 PCAValue(numberGroup) = sum(PhiValues((numberGroup-1)*NumberSummCAs+(1:NumberSummCAs)));
                PCAValue = [];
                InFilterPhase = 0;
                LoopFilterFPLL = ClassFilter;
                Order = [1,2];
                Bn = [5,5];
                
                LoopFilterFPLL.PrepareFilter(Order, Bn, dt);
                
                
                PCAValue = PCAs(numberGroup).*exp(1j*(-OutNCO));
                PCAs_correct (numberGroup)= PCAValue;
                PCAValues(1) = PCAValue;
                if length(PCAValues)>1
                    PCAValue1 = PCAValues(2);
                    PCAValue2 = PCAValues(1);
                       InFilterFreq = (atan2(imag(PCAValue2),real(PCAValue2)) - atan2(imag(PCAValue1),real(PCAValue1)))/dt;
%                   InFilterFreq = atan2((real(PCAValue1) * real(PCAValue2) + imag(PCAValue2) * imag(PCAValue1)),(real(PCAValue1) * imag(PCAValue2) - real(PCAValue2) * imag(PCAValue1)));
                    EstimatedF = InFilterFreq/(2*pi);
                    EstimatedFs (numberGroup)  = EstimatedF;
                    
                    EstimatedPhi = atan(imag(PCAValue2)/real(PCAValue2));
                    EstimatedPhis (numberGroup) = EstimatedPhi;
                    InFilterPhase = EstimatedPhi;
                    
                    [Output(numberGroup), VelocAcc(numberGroup), AccelAcc(numberGroup)] = LoopFilterFPLL.Step(InFilterFreq, InFilterPhase);
                    OutNCO =  0.5*dt*Output(numberGroup)+Acc;
                    Acc = Acc + dt*Output(numberGroup);
                    PCAValues = [];

                else
                    InFilterFreq=0;
                    EstimatedPhi = atan(imag(PCAValues(1))/real(PCAValues(1)));
                    EstimatedPhis (numberGroup) = EstimatedPhi;
                    InFilterPhase = EstimatedPhi;
                                     if   FLAG_FLL == 1
                    [Output(numberGroup), VelocAcc(numberGroup), AccelAcc(numberGroup)] = LoopFilterFPLL.Step(InFilterFreq, InFilterPhase);
                    OutNCO =  0.5*dt*Output(numberGroup)+Acc;
                    Acc = Acc + dt*Output(numberGroup);
                                     end
                               PCAValues(2) = PCAValues(1);
                               FLAG_FLL = 1;

                end
                
                
                numberGroup = numberGroup +1;
            end
            %������ 2 ���� ����������������, ������ ��� �����
            %�������� ���������� ���������, ������ ����� ������������,
            %������ �� ������� �������� �������

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



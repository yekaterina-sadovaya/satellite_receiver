function Res = P20_CohTrackSatsAndBitSync(inRes, Params)
%2. 20 ���������� ����������� �������� �������
% ������� �������������� �������� ��������� � ������� �������������
%
% ������� ����������
%   inRes - ��������� � ������������ ������, ����������� � Main;
%
% �������� ����������
%   Res - ���������, ������� ���������� �� inRes ����������� ������ ����,
%       �������� �������� ���� ���� � ����.

% �������������� �����������
Res = inRes;

%% ������������� ����������
Track = struct( ...
    'SamplesShifts', {cell(Res.Search.NumSats, 1)}, ...
    'CorVals',       {cell(Res.Search.NumSats, 1)} ...
    );
% ������ ������ cell-�������� SamplesShifts � CorVals �������� ��������
%   1xN, ��� N - ���������� �������� CA-���� ���������������� ��������,
%   ��������� � �����-������ (N ����� ���� ������ ��� ������
%   ���������).
% ������ ������� ������� SamplesShifts{k} - ���������� ��������,
%   ������� ���� ���������� � �����-������ �� ������ ����������������
%   ������� CA-����.
% ������ ������� ������� CorVals{k} - ����������� �������� ����������
%   ����� �������, ���������� ��������������� ������ CA-����, � �������
%   ��������.

BitSync = struct( ...
    'CAShifts', zeros(Res.Search.NumSats, 1), ...
    'Cors', zeros(Res.Search.NumSats, 19) ...
    );
% ������ ������� ������� CAShifts - ���������� �������� CA-����,
%   ������� ���� ���������� �� ������ ����.
% ������ ������ ������� Cors - ����������, �� ������� �������� �������
%   ������������ ������� �������������.

%% ��������� ����������
% ���������� �������� CA-���� ����� ��������� ��������������� ��
% ������� (NumCA2NextSync >= 1, NumCA2NextSync = 1 - ������������� ���
% ������� CA-����)
NumCA2NextSync = Params.P20_NonCohTrackSatsAndBitSync.NumCA2NextSync;

% �������� ���������� �������������� �������� CA-����, ������������ ���
% ������������� �� �������
HalfNumCA4Sync = Params.P20_NonCohTrackSatsAndBitSync.HalfNumCA4Sync;

% ���������� ����������� �������� ��������/������ ������������� ��
% �������
HalfCorLen = Params.P20_NonCohTrackSatsAndBitSync.HalfCorLen;

% ������, � ������� ������������ ����������� ����� ������������
% CA-�����
NumCA2Disp = Params.P20_NonCohTrackSatsAndBitSync.NumCA2Disp;

% ������������ ����� �������������� CA-����� (inf - �� ����� �����!)
MaxNumCA2Process = Params.P20_NonCohTrackSatsAndBitSync.MaxNumCA2Process;

% ���������� ���, ������������ ��� ������� �������������
NBits4Sync = Params.P20_NonCohTrackSatsAndBitSync.NBits4Sync;

%% ���������� ����������
Track.NumCA2NextSync   = 1;
Track.HalfNumCA4Sync   = HalfNumCA4Sync;
Track.HalfCorLen       = HalfCorLen;
Track.MaxNumCA2Process = MaxNumCA2Process;

BitSync.NBits4Sync     = NBits4Sync;

%% ���ר� ����������
% ����� CA-���� � ������ ������� �������������
CALen = 1023 * Res.File.R;

% ���������� �������� CA-����, ������������ �� ���� ���
CAPerBit = 20;

%% �������� ����� ������� - �������
[~,File] = ReadSignalFromFile(Res.File, 0, 0);
all_period_number = floor(File.SamplesLen/CALen);
dt=10^-3/CALen;
% shift_cor_vals = -HalfCorLen:HalfCorLen;

% ������ ���������
% fprintf('%s ������� ���������\n', datestr(now));
NumCASum = 20;
for SatNum  = 1:Res.Search.NumSats
%  for SatNum  = 6
    LoopFilterFPLL = ClassFilter;
    Order = 2;
    Bn =1;
    VelocAcc = 0;
    AccelAcc = 0;

    CACode = GenCACode(Res.Search.SatNums(SatNum),1);
    CACode(CACode==0)=-1;
    CACode = repmat( CACode,Res.File.R,1);
    CACode = reshape(CACode,1,[]);
    dt=10^-3/CALen;
    CACode = CACode.*exp(-(0:CALen-1)*1j*2*pi*Res.Search.FreqShifts(SatNum)*dt);
    
    Shift = Res.Search.SamplesShifts(SatNum);

      Shift_hard =0; %��� �������
      Shift_soft = 0; %��� �������
      Shift_soft_cum1 = 0; %��� �������
      Shift_hard_cum1 = 0; %��� �������
      Shift_soft_cum2 = 0; %��� �������
      Shift_hard_cum2 = 0; %��� �������
      mod_shifts_soft = []; %��� �������
      mod_shifts_hard = []; %��� �������
      Shifts = 0;
    Z = 0;
    T = NumCASum*1e-3;
    Prompts = [];
    Discrs = 0;
    FilteredDiscrs = 0;
    NCOVals = 0;
    L_nacop = 0 ;
    E_nacop = 0;

    k2 = 1;
    for k = 1:all_period_number
        EPL = [];
        [Signal, ~] = ReadSignalFromFile(Res.File, Shift-1, CALen+2); % ��������� EPL 
        
        % �� ������ 3 �������� --> ��� �� ��������� �� ���� E,P � L
        EPL = conv(Signal, fliplr(CACode), 'valid');

        if length(EPL)==3
            Track.CorVals{SatNum}(k) =  EPL(2); 
 
            EPL = abs(EPL);

            
            Prompts(k) = EPL(2);
            % ����������� �� 20
            E_nacop = E_nacop + EPL(1);
            L_nacop = L_nacop + EPL(3);
            if mod(k,NumCASum)==0
                % ������ �� �������������
                k2 = k2+1;
                Discrs(k2) = 0.5*(E_nacop - L_nacop)/(E_nacop + L_nacop);
                LoopFilterFPLL.PrepareFilter(Order, Bn, T, VelocAcc, AccelAcc);
                [FilteredDiscrs(k2), VelocAcc, AccelAcc] = LoopFilterFPLL.Step(Discrs(k2));
                Veloc(k2) = VelocAcc;
                Accel(k2) = AccelAcc;
                FilteredDiscrs(k2)=FilteredDiscrs(k2)/NumCASum; % ��������� �� ���������� ����������, ����� ���� ������
                L_nacop = 0; % dump 
                E_nacop = 0;
            end
            

            if k == 1
                % Shift - �������� ����
                % Shift-NCOVals(k)  -- soft
                NCOVals(k) = FilteredDiscrs(k2)*T;
                 Track.SamplesShifts{SatNum}(k) = Shift-NCOVals(k);
                        Shifts(k) = Shift; %��� �������
%                  Track.SamplesShifts{SatNum}(k) =  Shift; �� �� �����
%                  ����������
            else
                NCOVals(k) = NCOVals(k-1) + FilteredDiscrs(k2)*T;
                Track.SamplesShifts{SatNum}(k) =  Shift-NCOVals(k);  %!!!!!!
                       Shifts(k) = Shift; %��� �������
            end
            
            
            Shift = Shift + CALen;
            if NCOVals(k) > 0.5 % 0.5 �������� ��� E ��� L �������� ����� -- ��� ����� ������� ������������
                Shift = Shift - 1; %!!!!
                NCOVals(k) = NCOVals(k) - 1; %!!!!%!!!!%!!!!%!!!!
                Shift_hard = Shift_hard -1;  %��� �������
            elseif NCOVals(k) < -0.5
                Shift = Shift + 1; %!!!!
                NCOVals(k) = NCOVals(k) + 1;  %!!!!%!!!!%!!!!%!!!!
                Shift_hard = Shift_hard + 1; %��� �������
            end

        Shift_soft = Shift_soft + NCOVals(k); %��� �������
        
        Shift_soft_cum1 (1,k) = Shift_soft; %��� �������
        Shift_hard_cum1 (1,k) = Shift_hard; %��� �������

        end
        if k ==2 % �� ������ ����� ������ � ������ ����� ��������
            Track.CorVals{SatNum}(k) =  EPL(2); 
        end

    end

        mod_shifts_soft = mod(diff(Track.SamplesShifts{SatNum}),CALen) ;
        mod_shifts_soft(mod_shifts_soft> CALen/2) = mod_shifts_soft(mod_shifts_soft> CALen/2) -CALen ;  
        Shift_soft_cum2 = cumsum(mod_shifts_soft);
    
        mod_shifts_hard = mod(diff(Shifts),CALen) ;
        mod_shifts_hard(mod_shifts_hard> CALen/2) = mod_shifts_hard(mod_shifts_hard> CALen/2) -CALen ;  
        Shift_hard_cum2= cumsum(mod_shifts_hard);
    
    figure();
    plot(Prompts);
    figure();
    plot(abs(Track.CorVals{SatNum}));
    figure();
    plot(NCOVals);
    

    
    dCA_all{SatNum} = Track.CorVals{SatNum}(1:end-1).*conj(Track.CorVals{SatNum}(2:end));
    for nn = 0:NBits4Sync
        dCA{SatNum}(nn+1,:) = Track.CorVals{SatNum}(nn*20+(2:20)).*conj(Track.CorVals{SatNum}(nn*20+(1:20-1)));
    end
    sum_dCA{SatNum} = abs(sum(dCA{SatNum}));
    BitSync.Cors(SatNum,:)=sum_dCA{SatNum};
    [~,BitSync.CAShifts(SatNum)] = min( sum_dCA{SatNum} );
    
    
    figure ();
    plot(abs(Track.CorVals{SatNum}),'.');
    title(num2str(Res.Search.SatNums(SatNum)));
    
    
    figure ();
    title(num2str(Res.Search.SatNums(SatNum)));
    subplot(3,1,1);
    plot(sum_dCA{SatNum});
    subplot(3,1,2);
    plot(angle(Track.CorVals{SatNum})/pi,'.');
    
%     figure();
%     plot(Shift_soft_cum1,'color','green');
%     hold on;
%     plot(Shift_hard_cum1,'color','red');
%     title(num2str(Res.Search.SatNums(SatNum)));
%     grid on;
figure();
plot(Shift_soft_cum2,'color','green');
hold on;
plot(Shift_hard_cum2,'color','red');
title(num2str(Res.Search.SatNums(SatNum)));
grid on;
end
% ������� ����� ���� � ������������ � Res
Res.Track = Track;
Res.BitSync = BitSync;
end






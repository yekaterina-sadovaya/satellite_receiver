function Res = P20_NonCohTrackSatsAndBitSync(inRes, Params)
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
Track.NumCA2NextSync   = NumCA2NextSync;
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
Shifts_EPL = -HalfCorLen:HalfCorLen; % EPL -- early; promt; late 

for SatNum  = 1:Res.Search.NumSats
    dNCO =0; % NCO - numerical control oscillator - ��� ��������
    Time_shifts_nacop = 0;
    CACode = GenCACode(Res.Search.SatNums(SatNum),1);
    CACode(CACode==0)=-1;
    CACode = repmat(CACode,Res.File.R,1);
    CACode = reshape(CACode,1,[]);
    dt=10^-3/CALen;
    
    Shift = Res.Search.SamplesShifts(SatNum);
    dNCO = exp(-1j*2*pi*Res.Search.FreqShifts(SatNum)*dt);
    Init_NCO = 1;
    for k = 1:all_period_number
        Signal = [];
        [Signal, ~] = ReadSignalFromFile(Res.File, Shift, CALen);
        Signal = Signal.*Init_NCO.*dNCO.^(0:CALen-1);
        if ~isnan(Signal)
            
            Track.CorVals{SatNum}(k) =  Signal*CACode';
            Track.SamplesShifts{SatNum}(k) =  Shift;
            
            Time_shifts(k) = 0;

            if mod(k,NumCA2NextSync +1)==0 % ��� �������� ������ - 100; ������� ���� � ������ �����
                
                EPL_sum = [];
                ii = 1;
                
                for sh = Shifts_EPL
                    EPL = [];
                    [EPL, ~] = ReadSignalFromFile(Res.File, Shift-HalfNumCA4Sync*CALen+sh, 2*HalfNumCA4Sync*CALen);
                    EPL = EPL.*dNCO.^(0:2*HalfNumCA4Sync*CALen-1).*Init_NCO*dNCO^sh;
                    EPL_sum(ii) = sum(abs(CACode*reshape(EPL,CALen,2*HalfNumCA4Sync)));
                    ii = ii+1;
                end
 
                [~,max_index]=max(EPL_sum) ;
                Shift = Shift + Shifts_EPL(max_index);
                % Time_shifts(k) - ��� ��������� ���� �� ������ k-�� ����,
                % Time_shifts_nacop(k) - ��� ������ ����� ���� (��� � ���� ������������ ������ ������)
                %  ��� ���� ������ ������ ���������� ������������ ������ (� �������)
                Time_shifts(k) = Shifts_EPL(max_index);
            end
            
            if Time_shifts(k)~=0
                %��� ������������ ������� ������������!
                % Init �����, ����� �������� �� time shift; ����� ��� ���,
                % �� � ���� �� �������
                Init_NCO = Init_NCO*dNCO^Time_shifts(k);
            end
            
            if k~=1
                Time_shifts_nacop(k) = Time_shifts_nacop(k-1) + Time_shifts(k);
            end
            Shift = Shift + CALen;
            
            if mod(k,NumCA2Disp) == 0
                disp(strcat('����������: ',num2str(k),' �������� CA ����'));
            end
        end
    end
    
    %% �������� ����� ������� - ������� �������������
    %dCA_all{SatNum} = Track.CorVals{SatNum}(1:end-1).*conj(Track.CorVals{SatNum}(2:end));
    for nn = 0:NBits4Sync
        % ���� ������� 0/pi
        % � ����� ���� ������� �� nn, ������� ����� ���������� 
        dCA{SatNum}(nn+1,:) = Track.CorVals{SatNum}(nn*20+(2:20)).*conj(Track.CorVals{SatNum}(nn*20+(1:20-1)));
    end
    sum_dCA{SatNum} = abs(sum(dCA{SatNum})); % ���������� ����������, �� �� ���� � ��� ����������� ����� � ��� ���� - ����� �������
                                             % 1*exp(0) = 1; 1*exp(pi) = -1
    BitSync.Cors(SatNum,:)=sum_dCA{SatNum};
    [~,BitSync.CAShifts(SatNum)] = min( sum_dCA{SatNum} );
    
    
    figure ();
    plot(abs(Track.CorVals{SatNum}),'.');
    title(num2str(Res.Search.SatNums(SatNum)));
    
    
    figure ();
    subplot(3,1,1);
    plot(sum_dCA{SatNum});
    title(num2str(Res.Search.SatNums(SatNum)));
    subplot(3,1,2);
    plot(angle(Track.CorVals{SatNum})/pi,'.');
    subplot(3,1,3);
    plot(Time_shifts_nacop); 
    
end
% ������� ����� ���� � ������������ � Res
Res.Track = Track;
Res.BitSync = BitSync;
end




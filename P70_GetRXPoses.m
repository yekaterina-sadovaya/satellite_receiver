function Res = P70_GetRXPoses(inRes, Params)
%������� ��������� ������� � numsubframe ����� �� ��������
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
% Positioning = struct(  ...
%   'RXPoses', cell(N, M), ...
%   'CAStep', CAStep, ...
%   'isCommonRxTime', isCommonRxTime ...
% );

% ���������� ����� N cell-������� RXPoses ��������� � �����������
%   ��������� ��������� � ����������� ���������� TOW. ����������
%   �������� ������������ ����������� ���������� ��������� ��
%   ������������ ������ �������� � ������� �� ��������� CAStep,
%   ������������� ����.
% CAStep - ��� � �������� CA-���� ����� ��������� ������������
%   ���������.
% isCommonRxTime - ��������, ������������ ��� ��������� ���������
%   ���������: � ���� ����� �������� ��� � ������.

%% ��������� ����������
% ��� � �������� CA-���� ����� ��������� ������������ ���������. �����
% � �������� 6000 �������� CA-����, �������, ��������, CAStep = 1000
% ������� � ���������� 6 ��������� �� ���� �������.
CAStep = Params.P70_GetRXPoses.CAStep;
CALen = 1023 * Res.File.R;
% ������� ���������� ���������.
% isCommonRxTime = 1 - ���������� ��������� ����������� � ����������
%   ������  ������� ��������, ��������������� ������ ������
%   ������� GPS
% isCommonRxTime = 0 - ���������� ��������� ����������� � ������
%   ������� ������� ��������, ��������������� ���������� �����
%   ������� GPS
isCommonRxTime = Params.P70_GetRXPoses.isCommonRxTime;

% ���������� ������ ���������, ����������� ��� ���������� ���������:
% 'all' - ��� ��������;
% 'firstX' - ������ � ���������, �������� 'first5';
% [1, 2, 5, 7] - ���������� ������.
SatNums2Pos = Params.P70_GetRXPoses.SatNums2Pos;


Positioning = cell(size(Res.Ephemeris,1), 6000/CAStep);

% struct(  ...
%   'RXPoses', cell(size(Res.Ephemeris,1), 6000/CAStep), ...
%   'CAStep', CAStep, ...
%   'isCommonRxTime', isCommonRxTime ...
% );
%% ���ר� ����������
% �������� ������������� �������
dt = 1/Res.File.Fs;

% ��������� ���������� ������ ���������
if ischar(SatNums2Pos)
    if strcmp(SatNums2Pos, 'all')
        CurSatNums2Pos = 1:Res.Search.NumSats;
    else
        Buf = str2double(SatNums2Pos(6:end));
        CurSatNums2Pos = 1:Buf;
    end
else
    CurSatNums2Pos = SatNums2Pos;
end

%% ���ר� ���������
%���������� ��������� ����, ������� �� ����� �� ���� ������ ���������,
%���� �� ������ �������, ������ ��� � ������ �������� �� ������� ���, ���
%HOW ��������� �� ���� ���� ��������� ������ ������

kk = 1;
NumSubframe = 1;
for k = CurSatNums2Pos
    CAShifs(kk) = Res.Ephemeris{NumSubframe,k}.CANum;
    kk = kk+1;
end

[~,min_index]=min(CAShifs);
NumSubframe = 1;

n=1;
for N = 1:size(Res.Ephemeris,1) % N - ���� �� ���������� ���������
    for M = 1:6000/CAStep % 6000/CAStep ��������� ������� ��� �������� ������� �� ������� 
        
        SampleNum = Res.Track.SamplesShifts{min_index}(min(CAShifs))+1+(n-1)*CALen*CAStep; % Res.Track.SamplesShifts{min_index}(min(CAShifs)) ����� ���������� �����
        
        kk=1;
        Es = [];
        
        for k = CurSatNums2Pos
            
            RefTOW(kk) = Res.Ephemeris{N,k}.TOW;
            RefCANum(kk) = Res.Ephemeris{N,k}.CANum;
            
            TCA = 10^-3;
            [~,NumberCA(kk)] = min(abs(Res.Track.SamplesShifts{k}+1-SampleNum)); % ���������� ����� ��
            SampleNums(kk) = Res.Track.SamplesShifts{k}(NumberCA(kk)) ; % ���� � ������� �� ����� ��
            d(kk) = Res.Track.SamplesShifts{k}(NumberCA(kk)) +1 - SampleNum;
            inGPSTime(kk) = GettGPS (RefCANum(kk), RefTOW(kk), dt,TCA,NumberCA(kk),d(kk)); %��� �� t_c
            
            Es{kk} = Res.Ephemeris{N,k};
            kk = kk+1;
        end
        
        % ������� ������ �� ��� � �������� (��������) �.�. �� �� �����
        % ������ �� ������ �� ��������. ������ ��� � ��� ������ � ��� ��
        % ��� ���������������� -- ����� ������ �� �� ������������� ��
        % ������� ��������. � ������� ����� ���������� ���������
        [min_value_ref,min_index_ref] = min(RefCANum);
        % min_value_ref - �a���� �������� ��������� � ����������� ��. 
        % ������ ����� ������ ������� �������� ������ �� � ������������
        diff_CA =  NumberCA(min_index_ref) - min_value_ref ;
        kk=1;
        % inTimeShifts - ��� T(i)
        for k = CurSatNums2Pos
            inTimeShifts(kk) = (Res.Track.SamplesShifts{k}(RefCANum(kk)+diff_CA) +1- SampleNum).*dt;
            kk = kk+1;
        end
        
        UPos = P71_GetOneRXPos(Es, inGPSTime, inTimeShifts, ...
            SampleNums, Params);
        Positioning{N,M} = struct(  ...
            'RXPoses', struct(), ...
            'CAStep', CAStep, ...
            'isCommonRxTime', isCommonRxTime ...
            );
        
        Positioning{N,M}.RXPoses = UPos;
        
        RXPoses{N,M} = UPos;
    
        n=n+1;   
    end
end

Res.Positioning = Positioning;

P76_ExportResults(RXPoses, Params);

end

function inGPSTime = GettGPS(RefCANum, RefTOW, dt,TCA,NumberCA,d)

% ������� ���������� tGPS ��� ������� ������� SampleNum

% ������� ����������
%   SampleNum - ����� ������� ������, ��� �������� ���� ��������� �����
%       GPS;
%   SamplesNums - ������ ������ �������� CA-����� �������� ��������;
%   RefCANum - ����� CA-����, ������� �������� ������ � ��������, � �������
%       ��������� �������� RefTOW;
%   RefTOW - �������� RefTOW;
%   dt - �������� ������������� ������.
%
% �������� ����������
%   tGPS - ����� GPS � ������ SampleNum.

%  t = (TOW-1) + ��� ������ �
% TOW ������ �� ����� �������� � ������� �� ����������
% ������� ����� �������� ��� ����� TOW; ������� �� ������ �������� (�� 1-��
% ��) �� ������ ������� �������; ������ ����� ���������� ����� �����
% ������� �������� ������� ���� � �������� ������

% ���������

inGPSTime = (RefTOW-1)*6 + (NumberCA - RefCANum)*TCA - d*dt;

end
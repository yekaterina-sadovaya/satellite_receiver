function Res = P60_GatherSatsEphemeris(inRes, Params) %#ok<INUSD>

% ������� ����� ������������� ���������� ��� ���������, � ������� ����
% ������� ���� �� ���� �������� TOW_Count_Message
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
% Ephemeris = cell(N, Res.Search.NumSats);

% ���������� ����� cell-������� Ephemeris ��������� � �����������
% ��������� ��������� � ���������� ��������� TOW (�����������,
% ����������� ������ �� ��������, � ������� SatsData.isSat2Use = 1).
% ���������� Ephemeris �������� ���������, ���������� �������� ����
% ���������� �������, ������� � �������� ���������, � ����� ����������
% ����� ��������, ����� ������� CA ��������, ������������ �������� TOW,
% ��� �������� ����� ��� ����������. 
%?????????????���� ������ ���������� ������� ��
% �������, �� ������� cell-������� ������ ���� ������.

%% ��������� ����������

%% ���ר� ����������
% ����� ���� ����� ��������, ���������� ���������� Ephemeris
ENames = { ...
    ... % ��� ���� �� ��������� � ������������� ����������
    'SFNum', ... % ���������� ����� �������� ��������,
    ... % ���������������� ������� ������ (��������) Ephemeris
    'CANum', ... % ����� CA-���� ��������, � �������� ����������
    ... % ������� � ���������� ������� SFNum
    'TOW', ... % �������� TOW, ������������ � �������� � ����������
    ... % ������� SFNum. ��� �������� ���������� ��� ����
    ... % ��������� ����� ������ Ephemeris
    ...
    ... % ���� � ������������� �����������
    'WeekNumber', ...
    'CodesOnL2', ...
    'URA', ...
    'URA_in_meters', ...
    'SV_Health_Summary', ...
    'SV_Health', ...
    'IODC', ...
    'Flag_L2_P_Code', ...
    'T_GD', ...
    't_oc', ...
    'a_f2', ...
    'a_f1', ...
    'a_f0', ...
    'IODE', ...
    'C_rs', ...
    'Delta_n', ...
    'M_0', ...
    'C_uc', ...
    'e', ...
    'C_us', ...
    'sqrtA', ...
    't_oe', ...
    'Fit_Interval_Flag', ...
    'AODO', ...
    'C_ic', ...
    'Omega_0', ...
    'C_is', ...
    'i_0', ...
    'C_rc', ...
    'omega', ...
    'DOmega', ...
    'IODE', ...
    'IDOT', ...
    };

%% �������� ����� ������� - ���� �� ��������� ���������
% ������ ���������

% ��������� ���������� ������ ���������, ��� ������ �� ����� ��������
% �������� ���������
i = 1;
iiii = 1;
iiiii = 1;
for k = 1:Res.Search.NumSats
    if Res.SatsData.isSat2Use
        for kk = 1:size(Res.SatsData.HOW{k},2)
            Ephemeris{kk,i}.(ENames{1}) = kk; % ���������� ����� ��������� (�� �� ��� �������������� �����)
            Ephemeris{kk,i}.(ENames{2}) = Res.BitSync.CAShifts(k) + Res.SubFrames.BitShift(k)*20 + (kk-1)*20*300+1; % ���������� ����� �� ���� - �� ����� � �� �������� ���� �� ������� ��������
            Ephemeris{kk,i}.(ENames{3}) = Res.SatsData.HOW{k}(kk).TOW;
            
            if Res.SatsData.HOW{k}(kk).Subframe_ID ==  1
                for ii = kk:kk+4
                    if ii<=size(Res.SatsData.HOW{k},2)
                        for iii = 4:16
                            Ephemeris{ii,i}.(ENames{iii}) = Res.SatsData.SF1{k}(kk).(ENames{iii});
                        end
                    end
                end
            end
            if Res.SatsData.HOW{k}(kk).Subframe_ID ==  2
                for ii = kk-1:kk+3
                    if ii<=size(Res.SatsData.HOW{k},2)
                        for iii = 17:27
                            Ephemeris{ii,i}.(ENames{iii}) = Res.SatsData.SF2{k}(kk).(ENames{iii});
                        end
                    end
                end
            end
            if Res.SatsData.HOW{k}(kk).Subframe_ID ==  3
                for ii = kk-2:kk+2
                    if ii<=size(Res.SatsData.HOW{k},2)
                        for iii = 28:36
                            Ephemeris{ii,i}.(ENames{iii}) = Res.SatsData.SF3{k}(kk).(ENames{iii});
                        end
                    end
                end
            end 
            if Res.SatsData.HOW{k}(kk).Subframe_ID ==  4 & kk-3<0
               fill_gaps_4SF(iiii,1) =   k;
               fill_gaps_4SF(iiii,2) =   kk;
               iiii = iiii +1;
            end
            if Res.SatsData.HOW{k}(kk).Subframe_ID ==  5 & kk-4<0
               fill_gaps_5SF(iiiii,1) =   k;
               fill_gaps_5SF(iiiii,2) =   kk;
               iiiii = iiiii+1;
            end
        end     
    end
    

                
    % ��� ������� �������� ��������� ���������� ����� ��������, � �������
    % ����������� ������ �������� TOW ����� � ���������� ����������
    
    % ��������� ���������
    
    % ������ ��������� ��� ������� �������� ������� �������� ����������
    % ���������
    
    % ������� ����� ���� � ������������ � Res
    
    % ������ ���������
    i = i+1;
end
    %��� ������� ����� �������� IODC � IODE � ��� SF4 �������
                for i = 1:size(fill_gaps_5SF,1)
                        for iii = 4:16
                            Ephemeris{fill_gaps_5SF(i,2),fill_gaps_5SF(i,1)}.(ENames{iii}) = Res.SatsData.SF1{fill_gaps_5SF(i,1)}(fill_gaps_5SF(i,2)+1).(ENames{iii});
                        end
                        for iii = 17:27
                            Ephemeris{fill_gaps_5SF(i,2),fill_gaps_5SF(i,1)}.(ENames{iii}) = Res.SatsData.SF2{fill_gaps_5SF(i,1)}(fill_gaps_5SF(i,2)+2).(ENames{iii});
                        end

                        for iii = 28:36
                            Ephemeris{fill_gaps_5SF(i,2),fill_gaps_5SF(i,1)}.(ENames{iii}) = Res.SatsData.SF3{fill_gaps_5SF(i,1)}(fill_gaps_5SF(i,2)+3).(ENames{iii});
                        end
                end
Res.Ephemeris = Ephemeris;

end



function E = MakeEmptyE(ENames)
    % �������� ��� ����

    % ��������� � ����, �� ����������� � ������������� ������, ������������
    % ���������, ����� ���� isGathered �������� ������� (��. CheckAndAddE)
        E.SFNum = -1;
        E.CANum = -1;
        E.TOW   = -1;
end

function [outE, isNew] = CheckAndAddE(inE, SFNum, SFData, ENames)

% � ����������� �� ������ �������� �� ���������� �������� ���� IODC, ����
% IODE, ��������� � InE � ��� �� ��������� � SFData, ����� ����������
% �������� IODC � IODE � inE

% ���� ����� ��������� ����� E, �� ������� ���, � ��������� ������
% ��������� E �� �����
        
% ������� ������ ���� outE ���������� �� SFData
end

function isGathered = CheckE(E, ENames)
% ��������, �������� �� ������ ����
end
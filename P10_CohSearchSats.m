function Res = P10_CohSearchSats(inRes, Params)
    
% ������� ������������ ������ ��������� � �����-������
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
    Search = struct( ...
        'NumSats',       [], ... % ������, ���������� ��������� ���������
        'SatNums',       [], ... % ������ 1�NumSats � �������� ���������
            ... % ���������
        'SamplesShifts', [], ... % ������ 1�NumSats, ������ ������� -
            ... % ���������� ��������, ������� ����� ���������� � �����-
            ... % ������ �� ������ ������� ������� CA-���� ����������������
            ... % ��������
        'FreqShifts',    [], ... % ������ 1�NumSats �� ���������� ���������
            ... % ������� ��������� ��������� � ��
        'CorVals',       [], ... % ������ 1�NumSats ������������ ��������
            ... % ����� �������������� ������� ������������� �� �������
            ... % ��������, �� ������� ���� ������� ��������
        'AllCorVals',    zeros(1, 32) ... % ������ ������������ ��������
            ... % ���� �������������� �������
    );

%% ��������� ����������
    % ���������� ��������, ����������� ��� �����������.
    % ��� ������������ ����������� 1 <= NumCA2Search <= 10
        NumCA2Search = Params.P10_CohSearchSats.NumCA2Search;

    % ������ ����������� ������ ������������� ����������, ��
        CentralFreqs = Params.P10_CohSearchSats.CentralFreqs;

    % ����� �����������
        SearchThreshold = Params.P10_CohSearchSats.SearchThreshold;

%% ���������� ����������
    Search.NumCA2Search    = NumCA2Search;
    Search.CentralFreqs    = CentralFreqs;
    Search.SearchThreshold = SearchThreshold;
    
    %% ���ר� ����������
    % ���������� ��������������� ��������� ����������
    NumCFreqs = length(CentralFreqs);
    %���� ������� ������������� ������, ����� ��� ��������, ������ ��� � ����� �������� � �������
    % Res.File.R = 4;
    % ����� CA-���� � ������ ������� �������������
    CALen = 1023 * Res.File.R;
    %������� ��� ������� ��� �������������� ������� ����� ������� �������
    number_of_pieces = 2;
    %���� ������������� � �������� �������� ������� ���� ������ ��
    FlagOneCA = 1;
%% �������� ����� �������

Search.NumSats = 0;



dt=10^-3/CALen;
if FlagOneCA == 1
    [Signal,~] = ReadSignalFromFile(Res.File, 0, number_of_pieces*NumCA2Search*CALen+CALen-1);
    n = 0:number_of_pieces*CALen*NumCA2Search+CALen-2;
else
    [Signal,~] = ReadSignalFromFile(Res.File, 0, number_of_pieces*NumCA2Search*CALen+NumCA2Search*CALen-1);
    n = 0:number_of_pieces*CALen*NumCA2Search+NumCA2Search*CALen-2;
end

for i = 1:32
    correlation = [];
    summa = [];
    norm_summa = [];
    Signal_shift = [];
    CACode_one=[];
    CACode = [];
    CACode_one = GenCACode(i, 1);
    CACode_one(CACode_one==0)=-1;
    for r = 1:Res.File.R
        CACode(r,:) = CACode_one;
    end
        
    CACode = reshape(CACode,1,Res.File.R*size(CACode,2));
    
    if FlagOneCA == 0
        CACode = repmat(CACode,1,NumCA2Search);
    end
    
    for k = 1:NumCFreqs
        Signal_shift(k,:) = Signal.*exp(-n*1j*2*pi*CentralFreqs(k)*dt);
        correlation(k,:) = conv(Signal_shift(k,:),fliplr(CACode),'valid');
        for ii = 0:number_of_pieces-1
            if NumCA2Search~=1
                summa(k,ii*CALen+(1:CALen)) = abs(sum(reshape(correlation(k,ii*CALen*NumCA2Search+(1:CALen*NumCA2Search)),CALen,NumCA2Search)'));
            else
                summa(k,:) = abs(correlation(k,:));
            end
        end
    end
    for ii = 0:number_of_pieces-1
        norm_summa(:,ii*CALen+(1:CALen)) = summa(:,ii*CALen+(1:CALen))/mean(mean(summa(:,ii*CALen+(1:CALen))));
    end
    max_norm_summa = max(max(norm_summa));
    
    if  max_norm_summa > SearchThreshold
        [frequency_point,sample_point] = find(norm_summa == max_norm_summa);
        matrix(Search.NumSats+1,1) = max_norm_summa;
        matrix(Search.NumSats+1,2) = mod(sample_point,CALen)-1;
        matrix(Search.NumSats+1,3) = CentralFreqs(frequency_point);
        matrix(Search.NumSats+1,4) = i;
        Search.NumSats = Search.NumSats + 1;
    end
    Search.AllCorVals(i) =  max_norm_summa;
    
        figure();
        surf(1:size(norm_summa,2),CentralFreqs,norm_summa);
end

if size(matrix,1)>=0
matrix = sortrows(matrix, -1);
Search.FreqShifts = matrix(:,3)';
Search.SamplesShifts = matrix(:,2)';
Search.SatNums = matrix(:,4)';
Search.CorVals =  matrix(:,1)';
end

Res.Search  = Search;

figure ();
plot(Search.AllCorVals,'*-');
hold on;
plot([1,length(Search.AllCorVals)],[SearchThreshold,SearchThreshold])
grid on;
% 
% Res2 = Res;
% load('C:\Users\1\Desktop\11_semester\GPS\���\MATLAB\Results_Rate2_Coh\Rate2.mat');
% Res.Search.CorVals
% Res2.Search.CorVals
% Res.Search.SamplesShifts
% Res2.Search.SamplesShifts


function Res = P10_NonCohSearchSats(inRes, Params)
% 1) ���������� ��� ����� ���� ����� �� ������ � ������ ��������
% ������� �������������� ������ ��������� � �����-������
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
    'NumSats',       [0], ... % ������, ���������� ��������� ���������
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
NumCA2Search = Params.P10_NonCohSearchSats.NumCA2Search;
%         NumCA2Search = 1;
% ������ ����������� ������ ������������� ����������, ��
CentralFreqs = Params.P10_NonCohSearchSats.CentralFreqs;

% ����� �����������
SearchThreshold = Params.P10_NonCohSearchSats.SearchThreshold;

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
%���� ������������� ����� ��� ������� ������, ��������� ��� �����
FlagSqrt = 1;

%% �������� ����� �������
Search.NumSats = 0;

[Signal,~] = ReadSignalFromFile(Res.File, 0, NumCA2Search*CALen+CALen-1);

dt=10^-3/CALen;
n = 0:CALen-1;

%���� �� ���������� ���������
for i = 1:32
    correlation_result = [];
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

    
    for k = 1:NumCFreqs
         CACode_shift = [];
        CACode_shift = reshape(CACode,1,Res.File.R*size(CACode,2)).*exp(-n*1j*2*pi*CentralFreqs(k)*dt);
        %����� ���� �������� ������, �� � �������� �� ���
        Signal_shift(k,:) = Signal;
        correlation_result(k,:) = conv(Signal_shift(k,:),fliplr(CACode_shift),'valid');
            if NumCA2Search~=1
                if FlagSqrt == 1
                    summa(k,1:CALen) = sum(abs(reshape(correlation_result(k,1:CALen*NumCA2Search),CALen,NumCA2Search)'));
                else
                    ImValues = reshape(correlation_result(k,1:CALen*NumCA2Search),CALen,NumCA2Search);
                    SqrImValue = ImValues.*conj(ImValues); % ��� ������ ����� ���������
                    summa(k,1:CALen) = sum(SqrImValue');
                end
            else
                if FlagSqrt == 1
                    summa(k,:) = abs(correlation_result(k,:));
                else
                    ImValues = correlation_result(k,:);
                    SqrImValue = ImValues.*conj(ImValues);
                    summa(k,:) = sum(SqrImValue');
                end
            end
    end

    norm_summa = summa/mean(mean(summa));

    max_norm_summa = max(max(norm_summa));
 
    %��� ��������� � ������� � ����� �� ������������ �������� ��������
    %�������
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
%��� ��������
if size(matrix,1)>0
    matrix = sortrows(matrix, -1); % same as sortrows(matrix, 1, 'descend')
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

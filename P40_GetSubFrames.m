function Res = P40_GetSubFrames(inRes, Params)
%
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
Res.SubFrames = [];
    SubFrames = struct( ...
        'isSubFrameSync', zeros(Res.Search.NumSats, 1), ... 
        'BitSeqNum',      zeros(Res.Search.NumSats, 1), ...
        'BitShift',       zeros(Res.Search.NumSats, 1), ...
        'Words',          {cell(Res.Search.NumSats, 1)} ...
    );
    % ������ ������� ������� isSubFrameSync - ���� ���������� �����������
    %   �������������.
    % ������ ������� ������� BitSeqNum - ����� �������� ������, � �������
    %   ������� ��������� ������������� � ������� ��������.
    % ������ ������� ������� BitShift - ���������� ���, ������� ����
    %   ���������� �� ������ �������� ������ �� ������ ������� ��������.
    % ������ ������ cell-������� Words - cell-������ (Nx10), ��� N -
    %   ���������� ������������ ���������, ������ ������ - ������ 1�24 ���
    %   ��������������� �����, ���� CRC �������, � ������ ������, ���� CRC
    %   �� �������.

%% ��������� ����������
Bits = Res.Demod.Bits;
SatNums = Res.Search.SatNums;
NumSats = Res.Search.NumSats;
isDraw = Params.Main.isDraw;

%% ���ר� ����������

%% �������� ����� ������� - ���� �� ��������� ���������
for k = 1:Res.Search.NumSats
    [SubFrames.isSubFrameSync(k), SubFrames.BitSeqNum(k), SubFrames.BitShift(k)] = SubFrameSync(Bits{k}, isDraw, ...
        SatNums(k),NumSats);
    if SubFrames.isSubFrameSync(k)
        SubFrames.Words{k} = CheckFrames(Bits{k}(SubFrames.BitSeqNum(k),SubFrames.BitShift(k)-1:end));
    end
end
Res.SubFrames = SubFrames;
end

function Words_reshaped = CheckFrames(Shifted_Bits)
%
% �� �������� ������ ���������� ��� ��������� �����, � ������ �����
% ����������� CRC ������� �����, ���� CRC �������, �� �����������
% �������������� �����, � ��������� ������ ����������� ������ ������

%������ � ���  ���� ����� � �� ����� ���������� �����, ������ �� ����
%������� ��������� �� ��������� ���� ������������������� 
for i = 0:floor(length(Shifted_Bits)/30)-1 % �������� ������� ������ ���� � ��� ������ �� ��� ������� ������������������
    CRCisOk2=[];
    DWord = [];
    if i*30+32 <= length(Shifted_Bits) %������ ��� �������� �� ������� �� �� �� ������� ������������������, ����� ������ �� ����
        [CRCisOk2, DWord] = CheckCRC(Shifted_Bits(i*30+(1:32)));% ��������� �������� (CRC) ����� ��� ������ i
    end
    if CRCisOk2 %���� ��� �������, ���������� � �������
        Words{i+1} = DWord;
    else
        Words{i+1} = [];
    end
end
        Words_reshaped = reshape(Words(1:floor(size(Words,2)/10)*10),10,floor(size(Words,2)/10))'; %��� ���� ���� ������ ������ ��� �� ����� �������,�� ���� �����, ����� ��������
%        Words_reshaped(floor(size(Words,2)/10)+1,:) = Words(floor(size(Words,2)/10)*10+1:end);
end

function [isOk, BitSeqNum, BitShift] = SubFrameSync(Bits, isDraw, ...
   SatNum,NumSats)
%
% ������� ����������� �������������
%
% isOk      - ����, �����������, ������� ������������� ��� ���, ������ ���
%   ������ ���� ������� ������ ���� ���!
% BitSeqNum - ����� ������� ������������������, ��� ������� �������
%   �������������. �.�. ������������������, � ������� ���� ������ ��������.
% BitShift  - ���������� ���, ������� ����� ���������� � �������
%   ������������������ �� ������ ��������.

preambule = [1,0,0,0,1,0,1,1];
i=1;
preambule(preambule==0)=-1; % �������� 0 �� -1 ������ ��� ����� ���� ������ ������ ������������������ 

          Bits_new = Bits;
          Bits_new(Bits_new==0)=-1;
          for ii = 1:4
              subframe_search(ii,:)= abs(conv(Bits_new(ii,1:307),fliplr(preambule),'valid'));   
          end
          [seq_shift,bits_shift]=find(subframe_search==8);
          % bits_shift - �������  ����� ��������, ������� ����� 8, �����
          % ������� ��, ������� ������� ����� TLM (������ ��������). ������
          % ��������� TLM - ��� ����� ����������� �������������. �������
          % ������ ���� �� ��������� bits_shift
          for iii = 1:length(bits_shift)
              if bits_shift(iii)-2>0 %��� ������� �����, ������ ��� ���  �������� �������� 
                  %(CRC) ��� ����� ����� ��� �����(30 ���) � ��� 2 ���� ����������� �����, ����� �������� ������ ����������������
                [CRCisOk,~] = CheckCRC(Bits(seq_shift(iii),bits_shift(iii)+(-2:29)));
              else
                  %�� ���� ���, �������� �� ������ �������, �� �� �� �����
                  %����� -2, �� � ����� ��� ������� ������ ���  �����,
                  %����� �� ���������� ����� 300 ���, ��� �� � ����� �����
                  %300 ���������
                [CRCisOk,~] = CheckCRC(Bits(seq_shift(iii),bits_shift(iii)+300+(-2:29))); 
              end
          if CRCisOk == 1 
                BitShift = bits_shift(iii)-1;
                BitSeqNum = seq_shift(iii);
                isOk = 1;
          end
          end
          if isDraw ~=0
            figure ();
            plot(1:length(subframe_search(seq_shift,:)),subframe_search(seq_shift,:));
            title(num2str(SatNum));
          end
                       
end

function [isOk, DWord] = CheckCRC(EWord)
% ������� ������������ �������� CRC ��� ������ ����� ��������������
% ���������
% �� �����:
% EWord - ����� (������) � ����� ������ ����������� ����� � ������, �.�.
% ����� 32 ����.

% �� ������: 
%   isOk - 1, ���� CRC ��������, 0 � ��������� ������.
%   DWord - �������������� ����� (������), �.�. ����� 24 ����.
if EWord(2) == 1
    DWord = [mod(EWord(3:26)+EWord(2),2)];
else 
    DWord = [EWord(3:26)];
end
D25 = mod(EWord(1) + sum(DWord([1 2 3 5 6 10 11 12 13 14 17 18 20 23] )),2);
D26 = mod(EWord(2) + sum(DWord([2 3 4 6 7 11 12 13 14 15 18 19 21 24] )),2);
D27 = mod(EWord(1) + sum(DWord([1 3 4 5 7 8 12 13 14 15 16 19 20 22] )),2);
D28 = mod(EWord(2) + sum(DWord([2 4 5 6 8 9 13 14 15 16 17 20 21 23] )),2);
D29 = mod(EWord(2) + sum(DWord([1 3 5 6 7 9 10 14 15 16 17 18 21 22 24] )),2);
D30 = mod(EWord(1) + sum(DWord([3 5 6 8 9 10 11 13 15 19 22 23 24] )),2);
computed_CRC = [D25,D26,D27,D28,D29,D30];
if isequal(computed_CRC,EWord(27:32))
    isOk = 1;
else
    isOk = 0;
    DWord = [];
end
end


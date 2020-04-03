function Res = P40_GetSubFrames(inRes, Params)
%
%
% Входные переменные
%   inRes - структура с результатами модели, объявленная в Main;
%
% Выходные переменные
%   Res - структура, которая отличается от inRes добавлением нового поля,
%       описание которого дано ниже в коде.

% Пересохранение результатов
    Res = inRes;

%% ИНИЦИАЛИЗАЦИЯ РЕЗУЛЬТАТА
Res.SubFrames = [];
    SubFrames = struct( ...
        'isSubFrameSync', zeros(Res.Search.NumSats, 1), ... 
        'BitSeqNum',      zeros(Res.Search.NumSats, 1), ...
        'BitShift',       zeros(Res.Search.NumSats, 1), ...
        'Words',          {cell(Res.Search.NumSats, 1)} ...
    );
    % Каждый элемент массива isSubFrameSync - флаг успешности подкадровой
    %   синхронизации.
    % Каждый элемент массива BitSeqNum - номер битового потока, в котором
    %   удалось выполнить синхронизацию с началом подкадра.
    % Каждый элемент массива BitShift - количество бит, которые надо
    %   пропустить от начала битового потока до начала первого подкадра.
    % Каждая ячейка cell-массива Words - cell-массив (Nx10), где N -
    %   количество обработанных подкадров, каждая ячейка - массив 1х24 бит
    %   декодированного слова, если CRC сошлось, и пустой массив, если CRC
    %   не сошлось.

%% УСТАНОВКА ПАРАМЕТРОВ
Bits = Res.Demod.Bits;
SatNums = Res.Search.SatNums;
NumSats = Res.Search.NumSats;
isDraw = Params.Main.isDraw;

%% РАСЧЁТ ПАРАМЕТРОВ

%% ОСНОВНАЯ ЧАСТЬ ФУНКЦИИ - ЦИКЛ ПО НАЙДЕННЫМ СПУТНИКАМ
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
% Из битового потока выделяются все возможные кадры, в каждом кадре
% проверяется CRC каждого слова, если CRC сошлось, то сохраняется
% декодированное слово, в противном случае сохраняется пустой массив

%Теперь у нас  есть сдвиг и мы можем определить слова, подали на вход
%функции сдвинутую на найденный шифт последоваательность 
for i = 0:floor(length(Shifted_Bits)/30)-1 % примерно оценили колько слов у нас влазит во всю битовую последовательность
    CRCisOk2=[];
    DWord = [];
    if i*30+32 <= length(Shifted_Bits) %каждый раз цениваем не вылазим ли мы за пределы последовательности, чтобы ошибок не было
        [CRCisOk2, DWord] = CheckCRC(Shifted_Bits(i*30+(1:32)));% проверяем четность (CRC) слова под буквой i
    end
    if CRCisOk2 %если оно сошлось, записываем в матрицу
        Words{i+1} = DWord;
    else
        Words{i+1} = [];
    end
end
        Words_reshaped = reshape(Words(1:floor(size(Words,2)/10)*10),10,floor(size(Words,2)/10))'; %тут пару слов теряем потому что не целый подкадр,но если нужно, можно дописать
%        Words_reshaped(floor(size(Words,2)/10)+1,:) = Words(floor(size(Words,2)/10)*10+1:end);
end

function [isOk, BitSeqNum, BitShift] = SubFrameSync(Bits, isDraw, ...
   SatNum,NumSats)
%
% Функция подкадровой синхронизации
%
% isOk      - флаг, указывающий, найдена синхронизация или нет, причём она
%   должна быть найдена только один раз!
% BitSeqNum - номер битовой последовательности, для которой найдена
%   синхронизация. т.е. последовательности, с которой надо дальше работать.
% BitShift  - количество бит, которые нужно пропустить в битовой
%   последовательности до начала подкадра.

preambule = [1,0,0,0,1,0,1,1];
i=1;
preambule(preambule==0)=-1; % заменяем 0 на -1 потому что могут быть разных знаков последовательности 

          Bits_new = Bits;
          Bits_new(Bits_new==0)=-1;
          for ii = 1:4
              subframe_search(ii,:)= abs(conv(Bits_new(ii,1:307),fliplr(preambule),'valid'));   
          end
          [seq_shift,bits_shift]=find(subframe_search==8);
          % bits_shift - нашлось  много значений, которые равны 8, нужно
          % выбрать те, которые реально слова TLM (начало подкадра). Первый
          % найденный TLM - это будет подкадровая синхронизация. Поэтому
          % делаем цикл по найденным bits_shift
          for iii = 1:length(bits_shift)
              if bits_shift(iii)-2>0 %это условие нужно, потому что для  проверки четности 
                  %(CRC) нам нужно взять все слово(30 бит) и еще 2 бита предыдущего слова, чтобы избежать ошибки неопределенности
                [CRCisOk,~] = CheckCRC(Bits(seq_shift(iii),bits_shift(iii)+(-2:29)));
              else
                  %но если пик, например на первой позиции, то мы не можем
                  %взять -2, но а вдруг это реально нужный пик  слова,
                  %тогда он повориться через 300 бит, вот мы и берем через
                  %300 проверяем
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
% Функиця осуществляет проверку CRC для одного слова навигационного
% сообщения
% На входе:
% EWord - слово (строка) с двумя битами предыдущего слова в начале, т.е.
% всего 32 бита.

% На выходе: 
%   isOk - 1, если CRC сходится, 0 в противном случае.
%   DWord - декодированное слово (строка), т.е. всего 24 бита.
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


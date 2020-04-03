%этот тоже строит тело неопределенности, но тут еще я проверяла
%корреляционные свойства, а это не надо
clc;
clear all;
close all;
%%
%Формирование СА кода, я формировала не с помощью написанной функции, чтобы
%посмотреть корреляционные свойства последовательностей Голда и М
SigNum = 1;
NumCycles = 1;
delay = [5;6;7;8;17;18;139;140;141;251;252;254;255;256;257;258;469;470;471;472;473;474;509;512;513;514;515;516;859;860;861;862;863;950;947;948;950;67;103;91;19;679;225;625;946;638;161;1001;554;280;710;709;775;864;558;220;397;55;898;759;367;299;1018];
 
G1 = M_seq_gen([0,0,1,0,0,0,0,0,0,1],[1,1,1,1,1,1,1,1,1,1]);
G2 = M_seq_gen([0,1,1,0,0,1,0,1,1,1],[1,1,1,1,1,1,1,1,1,1]);
G2_i = circshift(G2,delay(SigNum));
CACode_1period = mod((G1+G2_i),2);
CACode=[];
for i = 1:NumCycles
    CACode = [CACode,CACode_1period];
end


G1(G1==0)=-1;
G2(G2==0)=-1;
CACode(CACode==0)=-1;

% cor1 = xcorr(G1,G2);
% cor2 = xcorr(G1,CACode);
% 
% cor3 = xcorr(G1,G1);
% cor4 = xcorr(G2,G2);
% cor5 = xcorr(CACode,CACode);

%%
%Сдвиг по Доплеру vs. корреляция для одного периода CA кода  
%Тут получается на один чип по одному отсчету, а в сигнале будет 2 отсчета
% dt=1*10^-3/1023;
% 
% df=-5*10^3:100:5*10^3;
% n = 0:1022;
% 
% CACode_shift=CACode.*exp(1j*2*pi*df'.*(n.*ones(length(df),1))*dt);
% 
% for i = 1:length(df)
%     matrix_cor(i,:)=xcorr(CACode_shift(i,:),CACode);
% end
% 
% figure (1);
% surf((-1022:1022),df,abs(matrix_cor));
% grid on;
% figure (2);
% plot(-1022:1022,abs(matrix_cor(51,:)));
% grid on;
% figure (3);
% plot(df,abs(matrix_cor(:,1023)));
% grid on;
% figure (6);
% plot(df,10*log10(abs(matrix_cor(:,1023))/max(abs(matrix_cor(:,1023)))));
% grid on;
%%
%Сравнение M-sequence и Gold по боковым лепесткам
% figure(4);
% plot(cor1,'color','red');
% hold on;
% plot(cor2,'color','blue');
% 
% figure(5);
% plot(cor5,'color','blue');
% hold on;
% plot(cor4,'color','magenta');
% hold on;
% plot(cor3,'color','yellow');
%%
%Сдвиг по Доплеру vs. корреляция для трех периодов CA кода
dt=1*10^-3/1023;
% %Для одинаковых СА кодов
% CACode_3 = [CACode,CACode,CACode];


% %Для разных СА кодов
% CACode2 = GenCACode(2, 1);
% CACode3 = GenCACode(3, 1);
% CACode2(find(CACode2==0))=-1;
% CACode3(find(CACode3==0))=-1;
% CACode_3 = [CACode2,CACode,CACode3];

%Для одинаковых СА кодов с разными знаками (наихудший случай)
CACode_3 = [CACode,CACode,-CACode];

%Удаляем 1 чип
CACode_3(1) = [];
CACode_3(length(CACode_3)) = [];
df=-5*10^3:100:5*10^3;
n = 0:3*1023-3;

CACode_shift=CACode_3.*exp(1j*2*pi*df'.*(n.*ones(length(df),1))*dt);

for i = 1:length(df)
    matrix_cor(i,:)=conv(CACode_shift(i,:),fliplr(CACode),'valid');
end

figure (7);
surf((-1023+1:1023-1),df,abs(matrix_cor));
grid on;
figure (8);
plot(-1023+1:1023-1,abs(matrix_cor(51,:)));
grid on;
figure (9);
plot(df,abs(matrix_cor(:,1023)));
grid on;
figure (10);
plot(df,10*log10(abs(matrix_cor(:,1023))/max(abs(matrix_cor(:,1023)))));
grid on;

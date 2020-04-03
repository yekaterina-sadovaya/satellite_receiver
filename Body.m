clc;
clear all;
close all;
%—троим тело неопределенности по пунктам из тетради
%1. 1)
LenCA = 1023;

FreqStep1 = 100;

dt=1*10^-3/LenCA;
SigNum1 = 1;
CACode1_3 = GenCACode(SigNum1, 3);
CACode1 = GenCACode(SigNum1, 1);
CACode1_3(CACode1_3==0)=-1;
CACode1(CACode1==0)=-1;
CACode1_3(1)= [];
CACode1_3(length(CACode1_3)) = [];
df1=-5*10^3:FreqStep1:5*10^3;
%¬ычитаем R1*2 убрать должны с боков по чипу, а так как если отсчеты могут
%удваиватьс€, то убираем 2 умножить на коэффициент передискретизации
n = 0:3*LenCA-3;
%это все равно что мы сделаем цикл по t потом по f и заполним матрицу
%столбик f умножаетс€ на каждый столбик из матрицы 
CACode_shift1=CACode1_3.*exp(1j*2*pi*df1'.*repmat(n,length(df1),1)*dt);

for i = 1:length(df1)
    matrix_cor1(i,:)=conv(CACode_shift1(i,:),fliplr(CACode1),'valid');
end

figure ();
surf((-LenCA+1:LenCA-1),df1,abs(matrix_cor1));
grid on;
%2)
figure ();
plot(-LenCA+1:LenCA-1,abs(matrix_cor1(51,:)));
grid on;
%3)
figure ();
plot(df1,abs(matrix_cor1(:,1023)));
grid on;
%4)
figure ();
plot(df1,10*log10(abs(matrix_cor1(:,1023))/max(abs(matrix_cor1(:,1023)))));
grid on;


%%

% %2. ≈сли нужно запустить пункт 2, то нужно закоментировать пункт 3 и
% наоборот
% 
FreqStep2 = 100;

dt=1*10^-3/LenCA;
SigNum2 = 2;
SigNum3 = 3;
CACode2_2 = GenCACode(SigNum2, 2);
CACode2_3 = [CACode2_2,GenCACode(SigNum3, 1)];
CACode2 = GenCACode(SigNum2, 1);
CACode2_3(CACode2_3==0)=-1;
CACode2(CACode2==0)=-1;
CACode2_3(1)= [];
CACode2_3(length(CACode2_3)) = [];
df2=-5*10^3:FreqStep2:5*10^3;
% ¬ычитаем R1*2 убрать должны с боков по чипу, а так как если отсчеты могут
% удваиватьс€, то убираем 2 умножить на коэффициент передискретизации
n = 0:3*LenCA-3;
% это все равно что мы сделаем цикл по t потом по f и заполним матрицу
% столбик f умножаетс€ на каждый столбик из матрицы 
CACode_shift2=CACode2_3.*exp(1j*2*pi*df2'.*repmat(n,length(df2),1)*dt);

for i = 1:length(df2)
    matrix_cor2(i,:)=conv(CACode_shift2(i,:),fliplr(CACode2),'valid');
end

figure ();
surf((-LenCA+1:LenCA-1),df2,abs(matrix_cor2));
grid on;
figure ();
plot(-LenCA+1:LenCA-1,abs(matrix_cor2(51,:)));
grid on;
figure ();
plot(df2,abs(matrix_cor2(:,1023)));
grid on;
figure ();
plot(df2,10*log10(abs(matrix_cor2(:,1023))/max(abs(matrix_cor2(:,1023)))));
grid on;

%%
% 3.
% 
% FreqStep2 = 100;
% 
% dt=1*10^-3/LenCA;
% SigNum2 = 2;
% % такой же как до этого, только одинаковый —ј код генерируем, надо же просто
% % знак помен€ть
% SigNum3 = 2;
% CACode2_2 = GenCACode(SigNum2, 2);
% CACode2_3 = [CACode2_2,-GenCACode(SigNum3, 1)];
% CACode2 = GenCACode(SigNum2, 1);
% CACode2_3(CACode2_3==0)=-1;
% CACode2(CACode2==0)=-1;
% CACode2_3(1)= [];
% CACode2_3(length(CACode2_3)) = [];
% df2=-5*10^3:FreqStep2:5*10^3;
% % ¬ычитаем R1*2 убрать должны с боков по чипу, а так как если отсчеты могут
% % удваиватьс€, то убираем 2 умножить на коэффициент передискретизации
% n = 0:3*LenCA-3;
% % это все равно что мы сделаем цикл по t потом по f и заполним матрицу
% % столбик f умножаетс€ на каждый столбик из матрицы 
% CACode_shift2=CACode2_3.*exp(1j*2*pi*df2'.*repmat(n,length(df2),1)*dt);
% 
% for i = 1:length(df2)
%     matrix_cor2(i,:)=conv(CACode_shift2(i,:),fliplr(CACode2),'valid');
% end
% 
% % figure ();
% % surf((-LenCA+1:LenCA-1),df2,abs(matrix_cor2));
% % grid on;
% figure ();
% plot(-LenCA+1:LenCA-1,abs(matrix_cor2(51,:)));
% grid on;
% % figure ();
% % plot(df2,abs(matrix_cor2(:,1023)));
% % grid on;
% % figure ();
% % plot(df2,10*log10(abs(matrix_cor2(:,1023))/max(abs(matrix_cor2(:,1023)))));
% % grid on;
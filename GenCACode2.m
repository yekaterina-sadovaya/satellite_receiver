function CACode = GenCACode2(SigNum, NumCycles)
%второй способ
%коды от 1 до 37
all_tap = [2,6;3,7;4,8;5,9;1,9;2,10;1,8;2,9;3,10;2,3;3,4;5,6;6,7;7,8;8,9;9,10;1,4;2,5;3,6;4,7;5,8;6,9;1,3;4,6;5,7;6,8;7,9;8,10;1,6;2,7;3,8;4,9;5,10;4,10;1,7;2,8;4,10];
tap = all_tap(SigNum,:);
G1 = M_seq_gen([0,0,1,0,0,0,0,0,0,1],[1,1,1,1,1,1,1,1,1,1]);
G2 = M_seq_gen2([0,1,1,0,0,1,0,1,1,1],[1,1,1,1,1,1,1,1,1,1],tap);

CACode_1period = mod((G1+G2),2);
CACode=[];
for i = 1:NumCycles
    CACode = [CACode,CACode_1period];
end
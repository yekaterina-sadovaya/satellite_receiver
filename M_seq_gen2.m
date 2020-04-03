function [M_sequence] = M_seq_gen2(polinom,phase,tap)

new_phase=phase;
for i=1:2^length(polinom)-1
    a=mod(new_phase*polinom',2);
    M_sequence(i) = mod(sum(new_phase(tap)),2);
    new_phase = circshift(new_phase,1);
    new_phase(1) = a;
end

function [M_sequence] = M_seq_gen(polinom,phase)

new_phase=phase;
for i=1:2^length(polinom)-1
    a=mod(new_phase*polinom',2);
    M_sequence(i) = new_phase(length(new_phase));
    new_phase = circshift(new_phase,1);
    new_phase(1) = a;
end

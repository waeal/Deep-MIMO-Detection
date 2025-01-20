function [QPSK_symbols] = QPSK_mapper(bitseq)
QPSK_table = [1 1i -1i -1]/sqrt(2);
for i=1:length(bitseq)/2
temp = bitseq(2*(i-1)+1)*2 +bitseq(2*(i-1)+2);
QPSK_symbols(i) =QPSK_table(temp+1);
end
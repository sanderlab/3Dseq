%%% function writing a fasta file of reverse complement sequences.
%%% Copyright (c) 2019 Frank J. Poelwijk, poelwijk@gmail.com. 

function FuncQnDrevcompl(filename)

tic;
revseqsDNA=fastaread(filename);
toc;
fwdseqsDNA=revseqsDNA;

for i=1:size(revseqsDNA,1)
    fwdseqsDNA(i).Sequence=seqrcomplement(revseqsDNA(i).Sequence);
end

fastawrite([filename(1:end-7),'rrc.fas'],fwdseqsDNA);
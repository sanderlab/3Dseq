%%% function aligning amino acid sequences with strong anchoring of
%%% N- and C-terminal parts. Copyright (c) 2019 Frank J. Poelwijk, 
%%% poelwijk@gmail.com. 

function SeqsMultiAlign = multialignfp(AlignSeqs,GapOpen,ExtendGap)

Nseq=size(AlignSeqs,2);
for i=1:Nseq
    AlignSeqs(i).Sequence=['WWWWWWWWWW',AlignSeqs(i).Sequence,'WWWWWWWWWW'];
end

SeqsMultiAlign=multialign(AlignSeqs,'GapOpen',GapOpen,'ExtendGap',ExtendGap);

for i=1:Nseq
    SeqsMultiAlign(i).Sequence=SeqsMultiAlign(i).Sequence(11:end-10);
end
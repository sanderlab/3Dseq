%%% function returning unique fasta sequences, and the number of times 
%%% a particular sequence occurs in the input fasta variable. 
%%% Copyright (c) 2019 Frank J. Poelwijk, poelwijk@gmail.com. 

function [fastaunique, multiplseq] = FuncFastaUnique(fastaraw)

tempcell=cell(size(fastaraw,1),1);
for i=1:size(fastaraw,1)
    tempcell{i}=fastaraw(i).Sequence;
end

[~,ia,ic] = unique(tempcell,'stable');
fastaunique=fastaraw(ia);

multiplseq=histcounts(ic,(1:max(ic)+1)-0.5)';

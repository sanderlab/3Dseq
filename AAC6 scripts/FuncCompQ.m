%%% function calculating the compound Q score for an Illumina MiSeq Qscore
%%% sequence. Copyright (c) 2019 Frank J. Poelwijk, poelwijk@gmail.com. 

function CompQ = FuncCompQ(numQscores)

numQscores(isnan(numQscores))=Inf;

CompQ = -10*log10(1-prod(1-10.^(-numQscores/10)));
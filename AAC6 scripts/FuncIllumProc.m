%%% function processing the Qscore filtered MiSeq sequence reads. Applies
%%% length and compound Qscore filtering, contaminant sequence filtering
%%% and performs mutation assessment.
%%% Copyright (c) 2019 Frank J. Poelwijk, poelwijk@gmail.com.
%%% 
%%% dependencies: FuncCompQ.m, FuncFastaUnique.m, multialignfp.m, 
%%% FuncSeqCompareEP.m, and text files containing wt amino acid 
%%% sequences, here PSEAB.txt and DEVLI.txt.

function FuncIllumProc(basename,QCompCutoff,writefastaflag,writedataflag,writecondflag)
%% first perform length and compound Qscore filtering; write the results as fasta
LSel=447; % length of wt DNA sequence
maxlinewidth=510;
MinReadLength=410;

fastaname=[basename,'_SeqsQ15.fas'];
Qscorename = [basename,'_QscoresQ15.txt'];
fastafilteredname=[basename,'_SeqsQM15C',num2str(QCompCutoff),'.fas'];

%%% parameters for conditioned alignment
T=1.15;
PSEABindelsAAallow=0;
SeqLengthAAallow=148;
PSEABsubsAAmin=2;

refseqnames=['PSEAB';'DEVLI'];
refseqdir='./RefSeq/';
refseqnum=size(refseqnames,1);

%%% read-in part etc
disp('Analyzing Qscores');
tic;
DNAseqs=fastaread(fastaname);
num_lines=size(DNAseqs,1);

quality=zeros(num_lines,maxlinewidth);
fid=fopen(Qscorename,'r');
for i=1:num_lines
    dummy=fgets(fid);
    dummy(dummy==9|dummy==10|dummy==13)=[];
    if length(dummy)<=maxlinewidth
        dummy=[dummy 33*ones(1,maxlinewidth-length(dummy))];
    else
        dummy=dummy(1:maxlinewidth);
    end
    quality(i,:)=dummy;
    if mod(i,5000)==0
        disp(['Analyzed ',num2str(i),' reads.']);
    end
end
fclose(fid);
qualitynum=quality-33;
qualitynum(qualitynum==0)=NaN;

CompQscores=zeros(size(qualitynum,1),1);
ReadLengths=zeros(size(qualitynum,1),1);
MinQscores=zeros(size(qualitynum,1),1);
for i=1:size(qualitynum,1)
    if sum(isnan(qualitynum(i,:)))>0
        ReadLengths(i)=find(isnan(qualitynum(i,:)),1,'first')-1;
    else
        ReadLengths(i)=size(qualitynum,2);
    end
    CompQscores(i)=FuncCompQ(qualitynum(i,:));
    MinQscores(i)=min(qualitynum(i,:));
end

KeepSeqs=false(num_lines,1);
KeepSeqs(CompQscores>=QCompCutoff & ReadLengths>=MinReadLength)=true;

DNAseqsFiltered=DNAseqs(KeepSeqs);
if writefastaflag==1
    fastawrite(fastafilteredname,DNAseqsFiltered);
end

QMinCutoffGrid=[1:40];
QCompCutoffGrid=[1:15];
NumPassed=zeros(numel(QMinCutoffGrid),numel(QCompCutoffGrid));
NumAtQs=zeros(numel(QMinCutoffGrid),numel(QCompCutoffGrid));
for i=1:numel(QMinCutoffGrid)
    for j=1:numel(QCompCutoffGrid)
        NumPassed(i,j)=sum(MinQscores>=QMinCutoffGrid(i) & CompQscores>=QCompCutoffGrid(j) & ReadLengths>=MinReadLength);
        NumAtQs(i,j)=sum(MinQscores==QMinCutoffGrid(i) & CompQscores>=QCompCutoffGrid(j) & CompQscores<QCompCutoffGrid(j)+1 & ReadLengths>=MinReadLength);
    end
end

if writedataflag==1
    save([fastafilteredname(1:end-3),'mat']);
end

if writefastaflag==1 %%%%% onlyFL!
    fastawrite([fastafilteredname(1:end-4),'_onlyFL.fas'],DNAseqsFiltered(ReadLengths(KeepSeqs)==LSel));
end
toc;

%% then perform mutation assessment

filename=[fastafilteredname(1:end-4),'_onlyFL.fas'];
runname=filename(1:end-4);

reffilename = [refseqnames(1,:),'.txt'];
fidref=fopen([refseqdir,reffilename],'r');
refseq=fgets(fidref);
fclose(fidref);

altreffilename = [refseqnames(2,:),'.txt'];
fidaltref=fopen([refseqdir,altreffilename],'r');
altrefseq=fgets(fidaltref);
fclose(fidaltref);

startpos=1; % MIND: changing startpos changes which part of the refseq will 
% be aligned, but the position numbering (and AA translation) starts from 
% the beginning of the refseq! Also for now startpos and stoppos need to
% be in frame with coding due to their usage in FuncSeqCompareEP!
stoppos=length(refseq); % or give number.

disp(['Analyzing conditioned alignment ',runname]);
%disp(['Input alignment is ',filename]);

%%%%%%%%%%%%%% data reading and initial filtering  %%%%%%%%%%%%
disp('Reading the alignment');
tic;
seqsDNA=fastaread(filename);
toc;

disp('Determining unique sequences');
tic;
for i=1:size(seqsDNA,1)
    seqsDNA(i).Header=['s',num2str(i)];
end

keepseqs=false(size(seqsDNA,1),1);
for i=1:size(seqsDNA,1)
    if length(seqsDNA(i).Sequence)==LSel
        keepseqs(i)=1;
    end
end
seqsDNA=seqsDNA(keepseqs);

Nseq=size(seqsDNA,1);

[seqsDNAunique, multiplseq] = FuncFastaUnique(seqsDNA);

Nsequnique=size(seqsDNAunique,1);

seqsAAunique=seqsDNAunique;
for i=1:Nsequnique
    seqsAAunique(i).Sequence=nt2aa(seqsDNAunique(i).Sequence, 'AlternativeStartCodons', false);
end

stoplocs=zeros(size(seqsDNAunique,1),1);
for i=1:Nsequnique
    stoplocs(i)=find(seqsAAunique(i).Sequence=='*',1,'first');
    seqsAAunique(i).Sequence(stoplocs(i):end)=[];
end
toc;


%%%%%%%%%%%%%%%%%%%%%%%% mutation assignment %%%%%%%%%%%%%%%%%

disp('Mutation assessment from DNA level');
tic;
%indelsL=false(Nsequnique,1);
for i=1:Nsequnique
	[MUTS,MUTpos,MUTcode,MUTcodeAA,MUTAAname,classMUT,score]=FuncSeqCompareEP(seqsDNAunique(i).Sequence, refseq, startpos, stoppos);
	datareport(i).mutname=seqsDNAunique(i).Header;
    datareport(i).data=cell(numel(MUTpos),1);
    datareport(i).nummut=numel(MUTpos);
    datareport(i).algnscore=score;
    datareport(i).mutpos=MUTpos;
    datareport(i).mutcode=MUTcode; % MUTcode 1-4 is subst to A,C,G,T; 5-8 is ins A,C,G,T, 9 is deletion
 	datareport(i).mutcodeAA=MUTcodeAA; % MUTcodeAA 0 is synon; 1 is non-synon; 2 is stop; 3 is ins/del/frameshift
    datareport(i).mutAAname=MUTAAname; % MUTAAname is amino acid 1-20; 0 is synonymous, 21 is stop, 22 is ins/del/frameshift
    datareport(i).nonfuncflag=0;
    if sum(datareport(i).mutcodeAA>1)>0
        datareport(i).nonfuncflag=1;
    end
    for j=1:numel(MUTpos)
     		datareport(i).data{j,1}=[MUTS(j,1),num2str(MUTpos(j)),MUTS(j,2),'   ',classMUT{j,1}];
    end
	[~,altMUTpos,~,altMUTcodeAA,~,~,~]=FuncSeqCompareEP(seqsDNAunique(i).Sequence, altrefseq, 1, length(altrefseq));
    datareport(i).DEVLImutposAA=unique(ceil(altMUTpos(altMUTcodeAA==1)/3));
    datareport(i).DEVLImuts=numel(datareport(i).DEVLImutposAA);
end

synmut=zeros(size(datareport,2),1);
nonsynmut=zeros(size(datareport,2),1);
mutstretchlength=zeros(size(datareport,2),1);
fracinstretch=zeros(size(datareport,2),1);
DEVLImuts=zeros(size(datareport,2),1);

for i=1:size(datareport,2)
    synmut(i)=sum(datareport(i).mutcodeAA==0);
    nonsynmut(i)=sum(datareport(i).mutcodeAA==1);
    if sum(datareport(i).mutcodeAA>1)>0
        mutstretchlength=1000;
        fracinstretch(i)=1;
    else
        if nonsynmut(i)==0
            mutstretchlength=0;
            fracinstretch(i)=0;
        else
            nonsynmutpos=unique(ceil(datareport(i).mutpos(datareport(i).mutcodeAA==1)/3));
            startstopstretch=diff([0;diff([nonsynmutpos;Inf])==1]);
            startstopnums=-nonsynmutpos.*startstopstretch;
            startstopnums(startstopnums==0)=[];
            if isempty(startstopnums)
                mutstretchlength=0;
            else
                numstretches=numel(startstopnums)/2;
                mutstretchlength=sum(startstopnums((1:numstretches)*2)+startstopnums((1:numstretches)*2-1)+1);
            end
            fracinstretch(i)=mutstretchlength/numel(nonsynmutpos); % changed stretchlength in coding mutations divided by all mutations? 
        end
    end
    DEVLImuts(i)=datareport(i).DEVLImuts; % because datareport isn't saved if it is huge, I need to save these guys    
end
toc;

disp('Mutation assessment AA level');
tic;
refseqAA=nt2aa(refseq);
refseqAA(find(refseqAA=='*',1,'first'):end)=[];
altrefseqAA=nt2aa(altrefseq);
altrefseqAA(find(altrefseqAA=='*',1,'first'):end)=[];

[seqsAAuniqueU, multiplseqAA] = FuncFastaUnique(seqsAAunique); % make algn unique at the AA level
NsequniqueAA = size(seqsAAuniqueU,1);

PSEABmutsAA=zeros(NsequniqueAA,1);
PSEABindelsAA=zeros(NsequniqueAA,1);
PSEABsubsAA=zeros(NsequniqueAA,1);
DEVLImutsAA=zeros(NsequniqueAA,1);
DEVLIindelsAA=zeros(NsequniqueAA,1);
DEVLIsubsAA=zeros(NsequniqueAA,1);
fracinstretchAA=zeros(NsequniqueAA,1);
SeqLengthAA=zeros(NsequniqueAA,1);

AlignSeqs(2).Sequence=refseqAA;
AlignSeqs(3).Sequence=altrefseqAA;
AlignSeqs(1).Header='testseq';
AlignSeqs(2).Header='WT_PSEAB';
AlignSeqs(3).Header='WT_DEVLI';
smoothwidth=4;
smoothmode=5;
smoothendtype=1;
DEVLIstretchAA=zeros(NsequniqueAA,1);

for i=1:NsequniqueAA
    SeqLengthAA(i)=length(seqsAAuniqueU(i).Sequence);
    [~, algn, ~]=nwalign(refseqAA,seqsAAuniqueU(i).Sequence,'GapOpen',15,'ExtendGap',3);
    PSEABmutsAA(i)=numel(algn(2,:))-sum(algn(2,:)=='|');
    PSEABindelsAA(i)=sum(algn(1,:)=='-')+sum(algn(3,:)=='-');
    PSEABsubsAA(i)=PSEABmutsAA(i)-PSEABindelsAA(i);
    mutlocsAA=algn(2,:)~='|';
    startstopnums=-(1:numel(mutlocsAA)).*diff([mutlocsAA,0]);
    startstopnums(startstopnums==0)=[];
        if isempty(startstopnums)
            mutstretchlength=0;
        else
            numstretches=numel(startstopnums)/2;
            mutstretchlength=startstopnums((1:numstretches)*2)+startstopnums((1:numstretches)*2-1);
            mutstretchlength(mutstretchlength==1)=0;
            mutstretchlength=sum(mutstretchlength);
        end
    fracinstretchAA(i)=mutstretchlength/PSEABmutsAA(i); 
    [~, algn, ~]=nwalign(altrefseqAA,seqsAAuniqueU(i).Sequence,'GapOpen',15,'ExtendGap',3);
    DEVLImutsAA(i)=numel(algn(2,:))-sum(algn(2,:)=='|');
    DEVLIindelsAA(i)=sum(algn(1,:)=='-')+sum(algn(3,:)=='-');
    DEVLIsubsAA(i)=DEVLImutsAA(i)-DEVLIindelsAA(i);
    
    AlignSeqs(1).Sequence=seqsAAuniqueU(i).Sequence;
    SeqsMultiAlign=multialignfp(AlignSeqs,15,3);

    dummyA=SeqsMultiAlign(1).Sequence==SeqsMultiAlign(2).Sequence;
    dummyB=SeqsMultiAlign(1).Sequence==SeqsMultiAlign(3).Sequence;

    dummyASmooth=fastsmooth(dummyA,smoothwidth,smoothmode,smoothendtype);
    dummyBSmooth=fastsmooth(dummyB,smoothwidth,smoothmode,smoothendtype);

    DEVLIstretchAA(i)=max((dummyBSmooth+0.1)./(dummyASmooth+0.1));
end
toc;

NoverS=nonsynmut./synmut;
NoverS(nonsynmut==0)=0;
aveN=mean(nonsynmut);
aveNoverS=mean(NoverS(~isinf(NoverS)));

maxbars=max(max(synmut),max(nonsynmut))+2;


if writecondflag==1
    WriteCond=PSEABindelsAA==PSEABindelsAAallow&SeqLengthAA==SeqLengthAAallow&DEVLIstretchAA<=T&PSEABsubsAA>PSEABsubsAAmin;
    
    seqAAwt.Header='WT_PSEAB';
    seqAAwt.Sequence=refseqAA;

    WriteSeqs=[seqAAwt; seqsAAuniqueU(WriteCond)];
    CondAlgnName=[fastafilteredname(1:end-4),'_',num2str(PSEABsubsAAmin),'_M15C',num2str(QCompCutoff),'_St',num2str(T*100),'.fas'];
    fastawrite(CondAlgnName,WriteSeqs);
end

save([runname,'oneFun']);
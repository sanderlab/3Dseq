%%% function loading the computed pdb files from a directory, and apply 
%%% the contact filtering, producing a filtered EC list and plotting the
%%% maps. Copyright (c) 2019 Frank J. Poelwijk, poelwijk@gmail.com.
%%% dependencies & inputs: ContactListsPSEAB.mat, pdbdists3.m, 
%%% FuncPlotContactmapCov.m 

clear 

L=134; % this is for AAC6 with dimer tail cut off
proteinname='aac6';
load('ContactListsPSEAB.mat');

pdbdirectory='./prediction_results/';
ECdatasetname='CouplingScores.csv';
maskedfilesuffix='masked';
pdblist='toprankingPDBs.mat';

writemaskeddataflag=0; % 1 for writing the data to .csv file
load_pdblist_flag=0; % 1 for supplying a set of pdbs, anything else to use all pdbs in the directory
compounding_flag=0; % 0 for scaled mask; 1 for hard cutoff with k_overlay; 2 for scaled mask with factor of 0.5 for 0 predicted contacts
k_overlay=10;

if load_pdblist_flag==1
    load([pdbdirectory,pdblist]);
    datafiles=ranking_filelist;
else
    datafiles=dir([pdbdirectory,'*.pdb']);
end

contacts_ALL=zeros(numel(datafiles),L,L);
for i=1:numel(datafiles)
    if load_pdblist_flag==1
        pred_struct=pdbread([pdbdirectory,char(datafiles{i,1})]);
    else
        pred_struct=pdbread([pdbdirectory,datafiles(i).name]);
    end
    for j=1:numel(pred_struct.Model.Atom)
        pred_struct.Model.Atom(j).element=pred_struct.Model.Atom(j).AtomName(1);
        pred_struct.Model.Atom(j).chainID='A';
    end
    distdata=pdbdists3(pred_struct,'A');
    contacts_ALL(i,:,:)=distdata.contacts;
end
startpos=pred_struct.Model.Atom(1).resSeq;
stoppos=pred_struct.Model.Atom(end).resSeq;

OVERLAY_contacts=squeeze(sum(contacts_ALL,1));

figure()
spy(OVERLAY_contacts);
title(['all compound contacts for protein ',proteinname,' resi ',num2str(startpos),'-',num2str(stoppos)]);

if compounding_flag==1
    figure()
    spy(OVERLAY_contacts>=k_overlay);
    title(['overlay map with k-overlay = ',num2str(k_overlay), ' for protein ',proteinname,' resi ',num2str(startpos),'-',num2str(stoppos)]); 
end

dataset=readtable([pdbdirectory,ECdatasetname]);

ShowSecStruct=0;
ShowActiveSiteRes=0;

startCS=strfind(ECdatasetname,'CouplingScores');
if ~isempty(startCS)
    fignamepart=[ECdatasetname(1:startCS-1),ECdatasetname(startCS+14:end-4)];
else
    fignamepart=datasetname;
end
fignamepart=[fignamepart,'_',maskedfilesuffix];
fignamepart(fignamepart=='_')='-';

datasetmasked=dataset;
datasetmasked.cn(datasetmasked.i==1)=-1;
datasetmasked.probability(datasetmasked.i==1)=0;
datasetmasked = sortrows(datasetmasked,[1 3],'ascend');

if compounding_flag==0
    MASK = OVERLAY_contacts; % so this will scale the ECs by how many structures have this pair of positions as a contact; if 0 structures, then multiply by 0; but might want to do by 0.5 or so.
elseif compounding_flag==1
    MASK = OVERLAY_contacts>=k_overlay;
elseif compounding_flag==2
    MASK = OVERLAY_contacts;
    MASK(MASK==0)=0.5;
end

for p=1:size(datasetmasked,1)
    if MASK(datasetmasked.i(p)-startpos+1,datasetmasked.j(p)-startpos+1)==0
        datasetmasked.cn(p)=-1;
        datasetmasked.probability(p)=0;
    end
end
datasetmasked=sortrows(datasetmasked,'cn','descend');

if writemaskeddataflag==1
    writetable(datasetmasked,[pdbdirectory,ECdatasetname(1:end-4),maskedfilesuffix,'.csv']);
end

EC(:,1)=datasetmasked.i;
EC(:,2)=datasetmasked.j;
EC(:,3)=datasetmasked.cn;

figure()
suptitle(fignamepart);
subplot(1,5,1)
FuncPlotContactmapCov(EC(:,1:2),MonomerContactList,L,'dimercontacts',DimerContactList,'LCutoff',floor(L/2),'PDistCutoff',5);
title([num2str(floor(L/2)),' ECs'],'Interpreter', 'none');
if ShowSecStruct==1
    FuncAddSecStruct(SecStruct);
end
if ShowActiveSiteRes==1
    FuncAddActiveSiteOverlay(ActiveSiteRes);
end

subplot(1,5,2)
FuncPlotContactmapCov(EC(:,1:2),MonomerContactList,L,'dimercontacts',DimerContactList,'LCutoff',L,'PDistCutoff',5);
title([num2str(L),' ECs'],'Interpreter', 'none');
if ShowSecStruct==1
    FuncAddSecStruct(SecStruct);
end
if ShowActiveSiteRes==1
    FuncAddActiveSiteOverlay(ActiveSiteRes);
end

EC(:,4)=zeros(size(EC,1),1);
for i=1:size(MonomerContactList,1)
    EC(ismember(EC(:,1:2),MonomerContactList(i,:),'rows'),4)=1;
end

if exist('DimerContactList','var')
    for i=1:size(DimerContactList,1)
        EC(ismember(EC(:,1:2),DimerContactList(i,:),'rows'),4)=2;
    end
end

PDistCutoff=[0 3 5 10 15 20];
numTopECs=3*L;
fracCorrectECs=zeros(numTopECs,numel(PDistCutoff));
for i=1:numTopECs
    for j=1:numel(PDistCutoff)
        ECtemp=EC(abs(EC(:,1)-EC(:,2))>=PDistCutoff(j),:);
        fracCorrectECs(i,j)=sum(ECtemp(1:i,4)>0)/i;
    end
end

subplot(1,5,3)
plot(1:numTopECs,fracCorrectECs(1:numTopECs,:),'LineWidth',3);
hold on
line([L/2 L/2],[0 1],'Color','red','LineStyle','--');
line([L L],[0 1],'Color','red','LineStyle','--');
hold off
title('TP as f^{ion} of #topECs and |i-j| cutoff');
xlabel('number of top ECs');
ylabel('fraction of correct ECs');
legend('0', '3', '5', '10', '15', '20');
axis square

subplot(1,5,4)
scatter(abs(EC(:,1)-EC(:,2)),EC(:,3),4,[0.5 0.5 0.5],'filled');
hold on
scatter(abs(EC(EC(:,4)==1,1)-EC(EC(:,4)==1,2)),EC(EC(:,4)==1,3),[],[0.7 0 0]);
scatter(abs(EC(EC(:,4)==2,1)-EC(EC(:,4)==2,2)),EC(EC(:,4)==2,3),[],[0 0.65 0]);
hold off
title('ECs as f^{ion} of |i-j|');
xlabel('|i-j|');
ylabel('EC value');
axis square

subplot(1,5,5)
FuncPlotContactmapCov(EC(:,1:2),MonomerContactList,L,'dimercontacts',DimerContactList,'LCutoff',3*L,'PDistCutoff',5);
title([num2str(3*L),' ECs'],'Interpreter', 'none');
if ShowSecStruct==1
    FuncAddSecStruct(SecStruct);
end
if ShowActiveSiteRes==1
    FuncAddActiveSiteOverlay(ActiveSiteRes);
end

disp(fignamepart);
if compounding_flag==1
    disp(['k_overlay is ',num2str(k_overlay)]);
end
disp(['TP (|i-j|=>5) for L/2 is ',num2str(fracCorrectECs(floor(L/2),3))]);
disp(['TP (|i-j|=>5) for L is ',num2str(fracCorrectECs(L,3))]);
disp(['TP (|i-j|=>5) for 2*L is ',num2str(fracCorrectECs(2*L,3))]);
disp(['TP (|i-j|=>5) for 3*L is ',num2str(fracCorrectECs(3*L,3))]);
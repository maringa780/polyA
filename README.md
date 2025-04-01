# Comparing the microbiome abundance landscape between WGS and RNA-seq platforms using the same samples
This repository contains Matlab codes for reproducing work cited below. 
#### Citation
Githaka, J.M. (2025) ‘Misuse’ of RNA-seq data in microbiome studies: A cautionary tale of poly(A). *mLife*. https://doi.org/ 

# Code
#### Step 1: Setup
Download all files in this repository.
Files *TableS8_T2T_KrakenUniq_BIO_Fullset.csv* and *TableS9_metadata_KrakenUniq_BIO_Fullset.csv* are part of Supplementary Tables in Sepich-Poore, G.D., McDonald, D., Kopylova, E. *et al.* Robustness of cancer microbiome signals over a broad range of methodological variation. *Oncogene* 43, 1127–1148 (2024).
#### Step 2: Read files
Read TCGA taxonomy data Table S8: T2T-KrakenUniq-MicrobialDB filtered genera abundances in TCGA (full set).
``` MATLAB
TCGAtaxa=readcell('TableS8_T2T_KrakenUniq_BIO_Fullset.csv')';
```
Read TCGA metadata Table S9: Metadata for T2T-KrakenUniq-MicrobialDB filtered genera in TCGA (full set).
``` MATLAB
TCGAmetadata=readcell('TableS9_metadata_KrakenUniq_BIO_Fullset.csv');
```
#### Step 3: Analysis & plots
Eliminate samples sequenced only once in both WGS and RNAseq (Note, not patient IDs!)
``` MATLAB
TCGAmetadata=TCGAmetadata(:,[1,5,17,19,21,24,29,31]);

n=tabulate(TCGAmetadata(2:end,ismember(TCGAmetadata(1,:),'tcga_sample_id')));% find number of time sample appears
TCGAmetadata(ismember(TCGAmetadata(:,ismember(TCGAmetadata(1,:),'tcga_sample_id')),n(cell2mat(n(:,2))<2,1)),:)=[];% must appear at least twice

idxWGS=contains(TCGAmetadata(:,ismember(TCGAmetadata(1,:),'experimental_strategy')),{'WGS'}); %index  WGS samples
idxRNA=contains(TCGAmetadata(:,ismember(TCGAmetadata(1,:),'experimental_strategy')),{'RNA-Seq'}); % index RNA samples

WGS.main=TCGAmetadata(idxWGS,:);
RNA.main=TCGAmetadata(idxRNA,:);
```

Find samples double sequenced with WGS (WGS1 vs WGS2)
``` MATLAB
n=tabulate(cellfun(@(x,y) strjoin({x,y}),WGS.main(:,4),WGS.main(:,7),'UniformOutput',false));
WGS.wgsVwgs=sortrows(WGS.main(ismember(cellfun(@(x,y) strjoin({x,y}),WGS.main(:,4),WGS.main(:,7),'UniformOutput',false), n(cell2mat(n(:,2))==2,1)),:),4);

[~,locb]=ismember(WGS.wgsVwgs(:,1),TCGAtaxa(1,:)); locb(locb==0)=[];% find the pairs from metadata file
WGS.wgsVwgsData=TCGAtaxa(:,[1;locb]);% align with WGS.twice

a=cell2mat(WGS.wgsVwgsData(2:end,2:2:end));
b=cell2mat(WGS.wgsVwgsData(2:end,3:2:end));
WGS.wgsVwgsRho=corrpercolumn(a,b);% get correlations
```
Find samples double sequenced with WGS & RNA (WGS vs RNA)
``` MATLAB
[lia,locb]=ismember(cellfun(@(x,y) strjoin({x,y}),WGS.main(:,4),WGS.main(:,7),'UniformOutput',false),...
                    cellfun(@(x,y) strjoin({x,y}),RNA.main(:,4),RNA.main(:,7),'UniformOutput',false));
locb(locb==0)=[];
WGS.wgsVrna=sortrows([WGS.main(lia,:);RNA.main(locb,:)],4);

[~,locb]=ismember(WGS.wgsVrna(:,1),TCGAtaxa(1,:)); locb(locb==0)=[];% find the pairs from metadata file
WGS.wgsVrnaData=TCGAtaxa(:,[1;locb]);% align with WGS.twice

a=cell2mat(WGS.wgsVrnaData(2:end,2:2:end));
b=cell2mat(WGS.wgsVrnaData(2:end,3:2:end));
WGS.wgsVrnaRho=corrpercolumn(a,b);% get correlations
```

Plot correlation boxplots
``` MATLAB
Rho=nan(max([numel(WGS.wgsVrnaRho),numel(WGS.wgsVwgsRho)]),2);
Rho(1:numel(WGS.wgsVwgsRho),1)=WGS.wgsVwgsRho;
Rho(1:numel(WGS.wgsVrnaRho),2)=WGS.wgsVrnaRho;

figure; boxplot(Rho,'Notch','on'); box off
pvalue=ranksum(WGS.wgsVwgsRho,WGS.wgsVrnaRho);

title(['pval = ' num2str(pvalue)])
ylabel('correlation coefficient')
H=gca; H.XTickLabel={'WGS1 vs WGS2','RNA vs WGS'};%H.XTickLabelRotation=90;
```
PCA plots and Heatmaps
``` MATLAB
clear ALLcombine
WGS1=WGS.wgsVwgsData(:,[1,2:2:end]); WGS1(1,2:end)={'WGS1'};
WGS2=WGS.wgsVwgsData(:,[1,3:2:end]); WGS2(1,2:end)={'WGS2'};

WGSrna=WGS.wgsVrnaData(:,[1,2:2:end]); WGSrna(1,2:end)={'WGSrna'};
RNAwgs=WGS.wgsVrnaData(:,[1,3:2:end]); RNAwgs(1,2:end)={'RNAwgs'};


[lia,locb]=ismember(WGS1(:,1),WGSrna(:,1)); 
ALLcombine=[WGS1(lia,:),WGS2(lia,2:end),WGSrna(locb,2:end),RNAwgs(locb,2:end)];
ALLcombine(2:end,2:end)=(num2cell(((cell2mat(ALLcombine(2:end,2:end))').*100)./(sum(cell2mat(ALLcombine(2:end,2:end))))'))';%normalize

D=cell2mat(ALLcombine(2:end,2:end));
S=ALLcombine(1,2:end);
S2=[WGS.wgsVwgs(1:2:end,6)',WGS.wgsVwgs(2:2:end,6)',WGS.wgsVrna(1:2:end,6)',WGS.wgsVrna(2:2:end,6)']; % primary site

Sunique=[S',[WGS.wgsVwgs(1:2:end,4);WGS.wgsVwgs(2:2:end,4);WGS.wgsVrna(1:2:end,4);WGS.wgsVrna(2:2:end,4)]];
Sunique=strrep(cellfun(@strjoin,cellfun(@(x,y){x,y},Sunique(:,1),Sunique(:,2),'UniformOutput',false),'UniformOutput',false),' ',':-')';

B=strrep(cellfun(@(x) x(end),cellfun(@(x) strsplit(x,'|s__'),ALLcombine(2:end,1),'UniformOutput',false)),'_',' ');

idx=contains(S,'WGS1') | contains(S,'WGS2');
```
PCA WGS1 vs WGS2
``` MATLAB
[~, zscores, pcvars] = pca(sqrt(D(:,idx))');
pcvars=pcvars./sum(pcvars) * 100;pcvars(1)
colors = distinguishable_colors(numel(unique(S(idx))));
figure('Position',[1000 1062 363 276]);
gscatter(zscores(:,1),zscores(:,2),S(idx)',colors,'.'),box off%legend off
xlabel(['PC1 (' num2str(round(pcvars(1),2)) '%)']);
ylabel(['PC2 (' num2str(round(pcvars(2),2)) '%)']);
```
PCA WGS vs RNA
``` MATLAB
[~, zscores, pcvars] = pca(sqrt(D(:,~idx))');
pcvars=pcvars./sum(pcvars) * 100;pcvars(1)
colors = distinguishable_colors(numel(unique(S(~idx))));
figure('Position',[1000 1062 363 276]);
gscatter(zscores(:,1),zscores(:,2),S(~idx)',colors,'.'),box off%legend off
xlabel(['PC1 (' num2str(round(pcvars(1),2)) '%)']);
ylabel(['PC2 (' num2str(round(pcvars(2),2)) '%)']);
```
Heatmap for top 15 genera
``` MATLAB
[~,idx]=sort(sum(logical(D(:,contains(S,'WGSrna') & contains(S2,'Brain'))),2)*(100/sum(contains(S,'WGSrna') & contains(S2,'Brain')))-...
    sum(logical(D(:,contains(S,'RNAwgs') & contains(S2,'Brain'))),2)*(100/sum(contains(S,'RNAwgs') & contains(S2,'Brain'))),'descend');
idx2=(contains(S,'RNAwgs') | contains(S,'WGSrna')) & contains(S2,'Brain');
n=15;

cgo = clustergram(D(idx(1:n),idx2)','RowPDist','euclidean','ColumnPDist','euclidean','Linkage','ward',...
     'ColumnLabels',B(idx(1:n)),'RowLabels',Sunique(idx2),'Colormap',redbluecmap(11));%,'ImputeFun', @knnimpute);%,'Symmetric',false, 'ColumnLabels',G

d=D(idx(1:n),idx2)';
b=B(idx(1:n));
s=Sunique(idx2);
[Lia1,Locb1]=ismember(cgo.ColumnLabels,b); Locb1(Locb1==0)=[]; disp(num2str([numel(Locb1),sum(Lia1)]))
[Lia2,Locb2]=ismember(cgo.RowLabels,s); Locb2(Locb2==0)=[]; disp(num2str([numel(Locb2),sum(Lia2)]))
cgo.RowLabels=s(1,Locb2);

%IMPORTANT, Add a color bar to the clustergram by clicking the Insert Colorbar button on the toolbar.

d=d(Locb2,Locb1);

cgo2=plot(cgo); %  colorbar
close(findall(groot,'Type','figure','Tag','Clustergram'))
labels=num2cell(d); labels(d>0)={''};labels(d==0)={'o'};
labelheatmap(cgo,labels)
hfig=gcf; hfig.Position=[1000 42 1222 1314]; % hfig=gcf; hfig.Position=[1000 42 366 1314];
```

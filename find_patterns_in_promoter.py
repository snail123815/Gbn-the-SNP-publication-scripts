import lzma, pickle, os
from tempfile import NamedTemporaryFile
from pybedtools import BedTool
from Bio import SeqIO
import pandas as pd


pattern = 'GATCWT'
promoterRegion = (-350, 50)

print('Load pattern')
patternPosBedLzma = f'../../SCO1839_SGR5654/ChIP-Seq_BGI/{pattern}.bed.lzma'
patternPosBedTemp = NamedTemporaryFile()
with lzma.open(patternPosBedLzma, 'rb') as handle:
    patternPosBedTemp.writelines(handle.readlines())
tempPatPosBed = patternPosBedTemp.name

print(f'Loading genome...')
genomeM145pkl = '.genomeM145.pkl.xz'
try:
    with lzma.open(genomeM145pkl, 'rb') as fh:
        genomeM145 = pickle.load(fh)
except FileNotFoundError:
    genomeM145 = SeqIO.read("../../../Resources/Genomes_Info/Streptomyces_coelicolor/M145.gb", 'genbank')
    with lzma.open(genomeM145pkl, 'wb') as fh:
        pickle.dump(genomeM145, fh)
        
print(f'Getting promoter region.')
promoterBedTemp = NamedTemporaryFile()
genomeName = 'NC_003888.3'
with open(promoterBedTemp.name, 'w') as fh:
    for feat in genomeM145.features:
        if feat.type == "gene":
            strand = feat.location.strand
            name = feat.qualifiers['locus_tag'][0]
            if strand == 1:
                pStart = feat.location.start + promoterRegion[0]
                pEnd = feat.location.start + promoterRegion[1]
            else:
                pStart = feat.location.end - promoterRegion[1]
                pEnd = feat.location.end - promoterRegion[0]
            line = '\t'.join(
                [genomeName, str(pStart), str(pEnd), name, '1', ("+" if strand == 1 else "-")]
            ) + '\n'
            fh.write(line)

print(f'Creating intersection')
promoterWithPattern = BedTool(promoterBedTemp.name).intersect(
    BedTool(patternPosBedTemp.name),
    wb=True
)
nPatEachGene = {}
for i, l in enumerate(promoterWithPattern):
    gene = l.fields[3]
    try:
        nPatEachGene[gene] += 1
    except KeyError:
        nPatEachGene[gene] = 1
nPatEachGene = [i for i in nPatEachGene.items()]
nPatEachGene.sort(key=lambda x: x[0])
nPatEachGene.sort(key=lambda x: x[1], reverse=True)
print(f'Found {len(nPatEachGene)} genes with {pattern} in their promoter region {promoterRegion}')

minPats = 2
print(f'Keep only genes with more than {minPats - 1} {pattern}.')
nPatEachGene = [g for g in nPatEachGene if g[1] >= minPats]
print(f'There are {len(nPatEachGene)} genes with more than {minPats-1} {pattern} in their promoter region {promoterRegion}')

targetGenes = [g[0] for g in nPatEachGene]

minPats = 3
print(f'Keep only genes with more than {minPats - 1} {pattern}.')
nPatEachGene = [g for g in nPatEachGene if g[1] >= minPats]
print(f'There are {len(nPatEachGene)} genes with more than {minPats-1} {pattern} in their promoter region {promoterRegion}')

minPats = 4
print(f'Keep only genes with more than {minPats - 1} {pattern}.')
nPatEachGene = [g for g in nPatEachGene if g[1] >= minPats]
print(f'There are {len(nPatEachGene)} genes with more than {minPats-1} {pattern} in their promoter region {promoterRegion}')

print(f'Read RNA-Seq result data')
from ChIP_Expression import *
deseq2_res = '../../../../GitProjects/Gbn-the-SNP-publication-scripts/Transcriptomics/dataTables/transformed/deseq2_comparisonResult_lfcShrink_22d17e_922b2a.tsv'
#deseq2_res = pd.read_csv(deseq2_res, sep='\t', header=0, index_col=0)
#targetTest = ['mu_wt24', 'mu_wt45']
#deseq2_res_target = deseq2_res.loc[[(l in targetTest) for l in deseq2_res.Label], :].loc[targetGenes, :]
#print(deseq2_res_target)

geneCompDict_shrink = readCompRes(deseq2_res)
tpmMeanDf = pd.read_csv('../../../../GitProjects/Gbn-the-SNP-publication-scripts/Transcriptomics/dataTables/mean_417110.tsv', sep='\t', index_col=0)
newtpmMeanDf = tpmMeanDf.loc[~tpmMeanDf.index.str.contains('SCOs02'),:]
#plotDistDiff(targetGenes, newtpmMeanDf, 'Dgbn_24', figsize=(5,5))
#plotDistDiff(targetGenes, newtpmMeanDf, 'Dgbn_45', figsize=(5,5))
#plotDistDiff(targetGenes, newtpmMeanDf, 'WT_24', figsize=(5,5))
#plotDistDiff(targetGenes, newtpmMeanDf, 'WT_45', figsize=(5,5))
#plotCompRes(targetGenes, geneCompDict_shrink,'mu_wt24','log2FC', title='', figsize=(5,5))
#plotCompRes(targetGenes, geneCompDict_shrink,'mu_wt45','log2FC', title='', figsize=(5,5))


#print(f'About the arms')
leftArmGenes = [f'SCO{str(i).zfill(4)}' for i in range(1, 1777)]
rightArmGenes = [f'SCO{str(i).zfill(4)}' for i in range(5762, 7847)]
#plotCompRes(leftArmGenes, geneCompDict_shrink, 'mu_wt24', 'log2FC', title='', figsize=(5,5))
#plotCompRes(rightArmGenes, geneCompDict_shrink, 'mu_wt24', 'log2FC', title='', figsize=(5,5))
#plotCompRes(leftArmGenes+rightArmGenes, geneCompDict_shrink, 'mu_wt24', 'log2FC', title='', figsize=(5,5))
#plotCompRes(leftArmGenes, geneCompDict_shrink, 'mu_wt45', 'log2FC', title='', figsize=(5,5))
#plotCompRes(rightArmGenes, geneCompDict_shrink, 'mu_wt45', 'log2FC', title='', figsize=(5,5))
#plotCompRes(leftArmGenes+rightArmGenes, geneCompDict_shrink, 'mu_wt45', 'log2FC', title='', figsize=(5,5))

print(f'Generate volcano plot for arm and core genes')

from oleveler import plotVolcano
import numpy as np

vstDf = pd.read_csv('../../../../GitProjects/Gbn-the-SNP-publication-scripts/Transcriptomics/dataTables/transformed/vst_2f3e41.tsv',
        sep='\t', index_col=0, header=0)
plotVolcano(geneCompDict_shrink['mu_wt24'],
        vstDf.loc[:, ['D24_1', 'D24_2', 'D24_3', 'C24_1', 'C24_2', 'C24_3']].mean(axis=1),
        highlights={'right arm':[(1, 0.576, 0.141, 0.5),rightArmGenes],
                    'left arm':[(0.149, 0.811, 0.450, 0.5),leftArmGenes]},
        xmax=6, ymax=80,
        figsize=(4,4),
        title='24h diff highlight arms',)
plotVolcano(geneCompDict_shrink['mu_wt45'],
        vstDf.loc[:, ['D45_1', 'D45_2', 'D45_3', 'C45_1', 'C45_2', 'C45_3']].mean(axis=1),
        highlights={'right arm':[(1, 0.576, 0.141, 0.5),rightArmGenes],
                    'left arm':[(0.149, 0.811, 0.450, 0.5),leftArmGenes]},
        xmax=6, ymax=160,
        figsize=(4,4),
        title='45h diff highlight arms',)

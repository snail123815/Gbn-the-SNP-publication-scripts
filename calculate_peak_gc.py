from plot_wholeGenome_Fold_Enrichment import *
import matplotlib.pyplot as plt
import pandas as pd
from pybedtools import BedTool
from tempfile import NamedTemporaryFile
import subprocess
import lzma, pickle
import os

# binning data
bin = 300
step = 1
lowGCthresh = 70
overlapPeak = 0.8
maxPeakWidth = 400


print(f'Loading genome...')
genomeM145pkl = '.genomeM145.pkl.xz'
try:
    with lzma.open(genomeM145pkl, 'rb') as fh:
        genomeM145 = pickle.load(fh)
except FileNotFoundError:
    genomeM145 = SeqIO.read("../../../Resources/Genomes_Info/Streptomyces_coelicolor/M145.gb", 'genbank')
    with lzma.open(genomeM145pkl, 'wb') as fh:
        pickle.dump(genomeM145, fh)
print(f'Binning {int(bin)}, step {int(step)}')

gcbinFile = f'.gcbins_{bin}_{step}.pkl.xz'
if os.path.isfile(gcbinFile):
    with lzma.open(gcbinFile, 'rb') as fh:
        gcContentData = pickle.load(fh)
else:
    gcContentData = calculateGCprecent(genomeM145, bin=bin, step=step)
    with lzma.open(gcbinFile, 'wb') as fh:
        pickle.dump(gcContentData, fh)
print(f'{len(gcContentData)} bins generated.')
print('GC content calculation complete.')

lowGCregions = [r for r in gcContentData if r[1] <= lowGCthresh]

print(f'{len(lowGCregions)} bins ({len(lowGCregions)/len(gcContentData):.1%}) are with low GC (<= {lowGCthresh}%).')

def gcContentToBed(gcContentData, chr='chr'):
    # gcContentData should be [(loc, gc), (loc, gc), ...] gc is 0.0-100.0
    # return:
    # chr    start    end    name    GC
    bedFile = NamedTemporaryFile()
    with open(bedFile.name, 'w') as fh:
        for r in gcContentData:
            start = int(r[0] - bin/2)
            end = int(r[0] + bin/2)
            name = f'{start}_{end}'
            fh.write(f'{chr}\t{start}\t{end}\t{name}\t{r[1]}\n')
    return bedFile


peakComBed = '../../SCO1839_SGR5654/ChIP-Seq_BGI/common_peaks_2402.bed'

with open(peakComBed, 'r') as fh:
    chr = fh.readline().split('\t')[0]

lowGCbedTemp = gcContentToBed(lowGCregions, chr=chr)

lowGCbed = BedTool(lowGCbedTemp.name)
peakComBed = BedTool(peakComBed)

print(f'Keep only narrow peaks with width less than {maxPeakWidth:.0f} (bin/minimum overlap = {bin}/{overlapPeak})')
peakcomBedFilteredBed_fh = NamedTemporaryFile()
BedTool(p for p in peakComBed if p.length < maxPeakWidth).saveas(peakcomBedFilteredBed_fh.name)
# The "generator" type of BedTool will be consumed after each use.
peakcomBedFilteredBed = BedTool(peakcomBedFilteredBed_fh.name)

nFiltered = len(peakcomBedFilteredBed)
print(f'{nFiltered}/{len(peakComBed)} peaks left ({nFiltered/len(peakComBed):.1%})')
lowGCpeakCom = peakcomBedFilteredBed.intersect(lowGCbed, wa=True, u=True, F=overlapPeak)
print(f'Found {len(lowGCpeakCom)} common peaks ({len(lowGCpeakCom)/nFiltered:.1%}) in low GC (<= {lowGCthresh}%) regions.')
peakComlowGC = lowGCbed.intersect(peakcomBedFilteredBed, wo=True, F=overlapPeak)
lowGCbedTemp.close()
fn = f'peak_common_in_lowGC_{lowGCthresh}_{bin}_{step}.tsv'
peakComlowGC.saveas(fn)
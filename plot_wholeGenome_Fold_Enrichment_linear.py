# adapted from /Users/durand.dc/Documents/works/Project_Sporulation/Sporulation_Scripts/ChIP_Seq/plotCoverage_gene_range.py

import lzma
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import bz2
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from math import ceil
from datetime import datetime
from Bio.SeqUtils import GC


def readData(compressedDataPath, rawDataPaths=None):
    """Return fold enrichment array by [[postion, fe, fe], ...]
    rawDataPaths should be a list of pathes"""
    if os.path.isfile(compressedDataPath):
        with bz2.open(compressedDataPath, 'rb') as data:
            feArray = pickle.load(data)
    elif rawDataPaths is None:
        raise FileNotFoundError(f'{compressedDataPath} not found')
    else:
        def parse(file):
            resArray = np.zeros((8667507), dtype=np.float32)
            data = np.genfromtxt(file, dtype=['i4', 'i4', 'f4'],
                                 delimiter='\t', usecols=[1, 2, 3])
            for start, end, fe in data:
                resArray[start:end] = fe
            return resArray
        print(f"Don't find any pre-loaded data at path\n{compressedDataPath}")
        print('Creating data array...')
        feArray = np.zeros((8667507, len(rawDataPaths) + 1), dtype=np.float32)
        feArray[:, 0] = np.arange(1, 8667508)
        print('Now loading from raw...')
        for file, idx in zip(rawDataPaths, range(1, len(rawDataPaths) + 1)):
            print(file)
            feArray[:, idx] = parse(file)
        print(f'Writing to file {compressedDataPath}')
        with bz2.open(compressedDataPath, 'wb') as data:
            pickle.dump(feArray, data)
    return feArray


def binData(dArr, bin=1e4, step=1e3, reverse=False):
    binned = []
    lengthData = len(dArr)
    steps = ceil(lengthData / step)
    for i in range(steps):
        middle = int(((2 * i) * step + bin) / 2)
        start = middle - int(bin / 2)
        end = middle + int(bin / 2)
        if end > lengthData:
            end = lengthData
        mean = float(dArr[start:end, 1].mean())
        mean = (-mean if reverse else mean)
        binned.append((middle, mean))
    return binned


def calculateGCprecent(seq, bin, step):
    seq = seq.seq
    lengthSeq = len(seq)
    binnedGC = []
    steps = ceil(lengthSeq / step)
    for i in range(steps):
        middle = int(((2 * i) * step + bin) / 2)
        start = middle - int(bin / 2)
        end = middle + int(bin / 2)
        if end > lengthSeq:
            end = lengthSeq
        binnedGC.append((middle, GC(seq[start:end])))
    return binnedGC


if __name__ == "__main__":
    # to use reletive path in this script
    dir_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(dir_path)

    print('Loading data...')

    compressedDataPath = 'ChIP_resultOnly/foldEnrichmentDf.pickle.bz2'
    outputDir = "Plots/"
    #
    #
    # dfLoc = pd.DataFrame.from_dict(dict(zip(
    #     dfAnnotation.index, dfAnnotation.Location.str.split('.').values))).transpose().drop(1, axis=1)
    # dfLoc.columns = ['start', 'end']
    # dfLoc['start'] = dfLoc['start'].str.extract('(\d+)', expand=False)
    # dfLoc['end'] = dfLoc['end'].str.extract('(\d+)', expand=False)
    # dfLoc = dfLoc.apply(pd.to_numeric)
    #
    #
    feArray = readData(compressedDataPath)
    print('Load complete.')
    print('Load genome')
    genomeM145pkl = '.genomeM145.pkl.xz'
    try:
        with lzma.open(genomeM145pkl, 'rb') as fh:
            genomeM145 = pickle.load(fh)
    except FileNotFoundError:
        genomeM145 = SeqIO.read("M145.gb", 'genbank')
        with lzma.open(genomeM145pkl, 'wb') as fh:
            pickle.dump(genomeM145, fh)
    print('Genome loaded.')

    dArr25h = feArray[:, (0, 1)]
    dArr48h = feArray[:, (0, 2)]
    # binning data
    # bin = 2e4
    # step = 1e3
    drawSvg = 1
    bin = 5000 # needs to be changed when plotting zoom-in regions
    step = 500

    if drawSvg:
        fmt = 'svg'
        #bin = bin*50
        #step = step*50 
    else:
        fmt = 'png'
    dBinned25h = binData(dArr25h, bin=bin, step=step, reverse=True)
    dBinned48h = binData(dArr48h, bin=bin, step=step)
    print('Binning complete.')

    # initialization of the plot
    m145Diagram = GenomeDiagram.Diagram('M145 genome')
    # gene track for annotation
    # graph track for plot enrichment
    gdTrackGenes = m145Diagram.new_track(1, name='genes', height=4, scale=False)
    gdTrackGC = m145Diagram.new_track(2, name='GC content', height=35,
                                      scale_color=colors.Color(0., 0., 0., 1),
                                      scale_largetick_interval=5e5, scale_smalltick_interval=5e4,
                                      scale_largetick_labels=False,
                                      scale_largeticks=0., scale_smallticks=0.,
                                      scale=True, axis_labels=False)
    gdTrackEnrichment = m145Diagram.new_track(3, name='Fold enrichment', height=22,
                                              scale_color=colors.Color(0.3765, 0.2196, 0.5765, 1),
                                              scale=True,
                                              scale_largetick_interval=5e6, scale_smalltick_interval=5e5,
                                              scale_largetick_labels=False,
                                              scale_largeticks=0., scale_smallticks=0,
                                              axis_labels=False)
    gdTrackGenome = m145Diagram.new_track(4, name='genome',
                                          scale_format='SInt',
                                          height=8,
                                          #scale_largetick_interval=5e5, scale_smalltick_interval=5e4,
                                          scale_largetick_interval=1e6, scale_smalltick_interval=1e5,
                                          scale_fontsize=30, scale_fontangle=0,
                                          scale_largeticks=1, scale_smallticks=0.5,
                                          axis_labels=False,
                                          )
    gdTrackOuterAnno = m145Diagram.new_track(5, name='features', height=4, scale=False)

    gdTGFeatures = gdTrackGenes.new_set(type='feature')
    gdTGCprecent = gdTrackGC.new_set(type='graph')
    gdTE25h = gdTrackEnrichment.new_set(type='graph')
    gdTE48h = gdTrackEnrichment.new_set(type='graph')
    gdTGenomeFeatures = gdTrackGenome.new_set(type='graph')
    gdTGOAFeatures = gdTrackOuterAnno.new_set(type='feature')

    # Fill genome track
    gdTGenomeFeatures.new_graph([(i,0) for i in np.arange(0, 8667507, step)] + [(8667507, 0)], center=0)
    # The genome has to be filled. Else error will occur when slicing the graph
    # oriC: 4270778 - 4272748 inclusive
    oriC = SeqFeature(location=FeatureLocation(4270777, 4272748), type='ori', id='oriC')
    gdTGOAFeatures.add_feature(oriC, color=colors.Color(1,0,0,1), label=True, name='oriC', label_size=25)
    
    # GC content
    gcContentData = calculateGCprecent(genomeM145, bin=bin, step=step)
    print('GC content calculation complete.')
    gdTGCprecent.new_graph(gcContentData, style='bar',
                           color=colors.Color(0.5, 0.5, 0.5, 1),
                           altcolor=colors.Color(0.3, 0.3, 0.3, 1),
                           center=70)

    print('Plot initialization complete.')

    # Add data in the plot
    gdTE25h.new_graph(dBinned25h, style='bar',
                      altcolor=colors.Color(0.2784, 0.4196, 0.8510, 1),
                      center=0)
    gdTE48h.new_graph(dBinned48h, style='bar',
                      color=colors.Color(0.8039, 0.04667, 0.8196, 1),
                      center=0)


    plotGenes = ['SCO1839']
    for feature in genomeM145.features:
        if feature.type != 'gene':
            continue
        if 'locus_tag' not in feature.qualifiers:
            continue
        geneId = feature.qualifiers['locus_tag'][0]
        #if geneId in plotGenes:
        #    gdTGFeatures.add_feature(feature, color=colors.blue,
        #            sigil="ARROW", label=True,
        #            name=geneId, label_size=25)
        #gdTGFeatures.add_feature(feature, color=colors.blue,
        #        sigil="ARROW", label=True,
        #        name=geneId, label_size=25)
    # lowGC region
    #lowGCtable = pd.read_excel('peak_common_in_lowGC_65_800_5.xlsx', sheet_name='4peaks_lowgc')
    #for i, r in lowGCtable.iterrows():
    #    start = int(r['zone start'])
    #    end = int(r['zone end'])
    #    feat = SeqFeature(location=FeatureLocation(start-1, end))
    #    gdTGFeatures.add_feature(feat, color=colors.Color(0,0.5,0.5, 1), name=f'lowGC_{i}')
    # highlighted region
    for i, r in enumerate([(3895278,3913450),(5098264,5166208),(5788306,5818220),(7561246, 7590427)]):
        start = int(r[0])
        end = int(r[1])
        feat = SeqFeature(location=FeatureLocation(start-1, end))
        gdTGFeatures.add_feature(feat, color=colors.Color(0,0.5,0.5, 1), name=f'lowGC_{i}')

    print('Plot data filling complete.')
    outputFigure = os.path.join(outputDir,
                                f"genomeMap_foldEnrichment_{datetime.now().strftime('%Y.%m.%d-%H.%M.%S')}.{fmt}")
    m145Diagram.draw(format='linear', circular=False,
                     pagesize=(140 * cm, 50 * cm),
                     x=0.01, y=0.01,
                     fragments=1,
                     start=0, end=8667507)
                     #pagesize=(30 * cm, 20 * cm), # for plot zoom-in regions
                     #start=3895278, end=3913450)
                     #start=5098264, end=5166208)
                     #start=5788306, end=5818220)
                     #start=7561246, end=7590427)
    m145Diagram.write(outputFigure, fmt)
    # plt.show()
    print(f'Finished.\n{outputFigure}')

import os
import logging
import sys
import json
from tempfile import NamedTemporaryFile
from time import sleep  # may be used in the notebook
from hashlib import md5
from collections import OrderedDict

import pybedtools
import numpy as np
import pandas as pd
import seaborn as sns

from Bio import SeqIO
from BCBio import GFF
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
from matplotlib.legend import Legend
from scipy.stats.stats import _validate_distribution
from scipy.stats import fisher_exact, mannwhitneyu
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression as PLS

from jupyterthemes import jtplot


logger = logging.getLogger()
logging.getLogger('matplotlib.font_manager').disabled = True
if len(logger.handlers) == 0:
    logFormatter = logging.Formatter(
        '%(asctime)s [%(levelname)s] - %(message)s'
    )
    logfhandler = logging.FileHandler('ChIP_Expression.log')
    logfhandler.setLevel(logging.INFO)
    logfhandler.setFormatter(logFormatter)
    logshandler = logging.StreamHandler(stream=sys.stdout)
    logshandler.setLevel(logging.DEBUG)
    logshandler.setFormatter(logFormatter)
    logger.addHandler(logfhandler)
    logger.addHandler(logshandler)
    logger.setLevel(logging.DEBUG)

def setDarkModePlotting(forceDark=False, forceWhite=False):
    if 'current_jupyterthemes_style_dark' not in globals():
        global current_jupyterthemes_style_dark
        current_jupyterthemes_style_dark = False
    if (current_jupyterthemes_style_dark and not forceDark) or forceWhite:
        logger.warning('Change plot style to default')
        jtplot.reset()
        current_jupyterthemes_style_dark = False
    else:
        logger.warning('Change plot style to dark (monokai)')
        jtplot.style(theme='monokai', context='notebook', ticks=True, grid=False)
        current_jupyterthemes_style_dark = True


def displayImage(src, background=None, **kwargs):
    from IPython.display import Image, SVG, display
    from IPython.core.display import HTML
    if src[-4:].lower() == '.svg':
        htmlstr = '<div class="jp-RenderedSVG jp-OutputArea-output " data-mime-type="image/svg+xml">'
        svgstr = SVG(data=src)._repr_svg_()
        if not isinstance(background, type(None)):
            bgstr = f'style="background-color:{background}"'
            svgstr = '<svg ' + bgstr + svgstr[4:]
        htmlstr += svgstr + '</div>'
        display(HTML(htmlstr))
    else:
        img = Image(data=src)
        fmt = img.format
        if fmt == 'png':
            htmlstr = '<div class="jp-RenderedImage jp-OutputArea-output ">'
            htmlstr += f'\n<img src="data:image/png;base64,{img._repr_png_()}"></div>'
            if not isinstance(background, type(None)):
                htmlstr = htmlstr[:-7] + f'style="background-color:{background}"' + htmlstr[-7:]
            display(HTML(htmlstr))
        else:
            logger.warning('background for none PNG image will be ignored')
            display(img)


def calHash(*args) -> str:
    """Produce 6 digit string with unlimited number of arguments passed in
    Designed in mind that all types of data can be calculated
    resulting the same hash across platform.
    Should be safe with nested dict but no guarantee
    """
    def orderDict(di):
        try:
            d = di.copy()
        except:
            d = di
        if isinstance(d, dict):
            od = OrderedDict(sorted(d.items()))
            for k, v in od.items():
                v = orderDict(v)
                od[str(k)] = v
            d = od
        elif isinstance(d, (list, tuple, set)):
            d = [orderDict(el) for el in d]
        else:
            d = str(d).strip()
        return d

    def hashDict(di):
        od = orderDict(di)
        ha = md5(
            json.dumps(
                od,
                sort_keys=True,
                ensure_ascii=True,
                default=str
            ).encode()
        ).digest()
        return ha

    hastr = ''.encode()
    for arg in args:
        if isinstance(arg, str):
            hastr += arg.encode()
        elif isinstance(arg, set):
            hastr += str(sorted(list())).encode()
        elif isinstance(arg, dict):
            hastr += hashDict(arg)
        elif isinstance(arg, pd.core.frame.DataFrame) or isinstance(arg, pd.core.series.Series):
            hastr += md5(arg.to_json().encode()).digest()
        elif isinstance(arg, PCA):
            hastr += arg.components_.tobytes()
        elif isinstance(arg, PLS):
            hastr += arg.x_loadings_.tobytes()
        else:
            hastr += str(arg).encode()
    return md5(hastr).hexdigest()[:6]


def safeLog(x):
    x = x[x>0]
    return np.log10(x)

def minusSafeLog(x):
    return -safeLog(x)

def getChangedGenes(compFile, comps=[],
                    filterTypes=['log2FC', 'adj.pvalue'],
                    thresholds=[(-np.inf,-1,1,np.inf), (-np.inf, 0.05)]):
    """Usage:
    filteredDf = getChangedGenes(
        'Transcriptomics/dataTables/transformed/deseq2_comparisonResult.tsv',
        comps=['mu_wt24', 'mu_wt45'])

    NOTE the returned data frame will contain duplicated id whin the given comparison name (comps) is more than one.

    - [x] Read comparison file from expression analysis (proteomics or transcriptomics)
    - [x] Return a **DataFrame** of genes that either up, down, changed; with significance or not.
    - [x] Able to set significance thresholds: lfc, adj.pvalue etc.
    - [x] test

    Args:
        compFile (str): Path to the comparison data tsv file
        comps (list, optional): comparison name as lsited in "Label" column, will use all if empty. Defaults to [].
        filterTypes (list, optional): column names used to filter. Defaults to ['log2FC', 'adj.pvalue'].
        thresholds (list, optional): list of tuples, each tuple has to be even in number. Will use the first one as left
        boundary; next one as right boundary; etc. The number of tuples needs to match filterTypes. Defaults to [(-np.inf,-1,1,np.inf), (-np.inf, 0.05)].

    Returns:
        filteredDF: filtered DataFrame
    """
    #- [ ] getChangedGenes(compFile, ...)
        #- [ ] Read comparison file from expression analysis (proteomics or transcriptomics)
        #- [ ] Return a **DataFrame** of genes that either up, down, changed; with significance or not.
        #- [ ] Able to set significance thresholds: lfc, adj.pvalue etc.
    assert len(filterTypes) == len(thresholds), f'Please give one threshold (or one tuple) fo reach filter type'
    allCompDf = pd.read_csv(compFile, sep='\t', index_col=0)
    allComps = list(set(allCompDf['Label'].to_list()))
    allComps.sort()

    # if comps not set, then use all data
    if len(comps) == 0:
        comps = allComps
    assert all(c in allComps for c in comps), f'Comparison names should be in {allComps}, not {comps}'
    acceptedFilters = list(allCompDf.columns)
    acceptedFilters.remove('Label')
    assert all(t in acceptedFilters for t in filterTypes), f'filterTypes needs to be subset of {acceptedFilters}'
    
    # filter based on comp names
    filteredDf = allCompDf[allCompDf['Label'].isin(comps)]
    for t,r in zip(filterTypes, thresholds):
        assert len(r)%2 == 0 ,\
            f'Need tuple contain even number of elements for low and high bound. Use np.inf and -np.inf if needed.'
        tfDf = []
        for i in range(len(r)):
            if i%2 == 1: continue
            tfDf.append(filteredDf[filteredDf[t].between(r[i], r[i+1])])
        filteredDf = pd.concat(tfDf, axis=0)
    
    return filteredDf


def parsePeakFile(peakFile):
    """read peak .xls file which is an output of MACS2. The file should actually be a tsv file.

    Args:
        peakFile (str): path to peak file

    Returns:
        peakDf: with columns ['name', 'start', 'end', 'length', 'abs_summit', 'pileup', 'fold_enrichment', '-log10(qvalue)']
    """
    peakDf = pd.read_csv(peakFile, sep='\t', comment='#', skip_blank_lines=True, skipinitialspace=True)
    peakDf = peakDf.loc[:, ['name', 'start', 'end', 'length', 'abs_summit', 'pileup', 'fold_enrichment', '-log10(qvalue)']]
    peakDf = peakDf.set_index('name', drop=True)
    return peakDf



def parseBedFile(bedFile):
    """
    parse bed file generated by MACS2
    simplePeakDf.columns = ['chromosome', 'start', 'end', 'length', 'score', 'strand'] 
    # note there is no summit info

    simplePeakDf = parseBedFile(r"/mnt/c/Users/duc/Documents/tmp_works/ChIP-48h_peaks.bed")

    Args:
        bedFile (str): path to bed file

    Returns:
        simplePeakDf
    """
    cols = ['chromosome', 'start', 'end', 'length', 'score', 'strand'] # note there is no summit info

    lines = []
    n = 3 # to track max columns of input bedFile, minimum is 3
    with open(bedFile, 'r') as b:
        for i,l in enumerate(b): # instead of using pd.read, this is easier to determine the header.
            if l.startswith('track name'):
                continue
            if l.startswith('#'):
                continue
            es = l.strip().split('\t')
            if len(es) < 5:
                continue
            chrom = es[0]
            start = int(es[1])
            end = int(es[2])
            try:
                name = es[3]
                n = 4
            except:
                name = f"name{i}"
            try:
                score = es[4]
                n = 5
            except:
                score = '0'
            try:
                strand = es[5]
                n = 6
            except:
                strand = '+'
            row = [chrom, start, end, end-start, score, strand]
            #print(n)
            #print(es)
            #print(cols)
            #print(row)
            if len(es) > n:
                row.extend(es[n:])
            if len(row) > len(cols):
                cols.extend([f'_{i+1}' for i in range(len(row)-len(cols))])
            #print(cols)
            #print(row)
            #raise ValueError
            lines.append(pd.DataFrame([row], index=[name], columns=cols))
    bedDf = pd.concat(lines, axis=0)
                            
    return bedDf


def annotationToBed(annoFile, relToRef=None, ref='start',
                    bedType='gene', idQualifier='locus_tag',
                    withSynonym=None, bedFile=None):
    """read genbank file or gff file, write info as bed file.
    dependency: parseBedFile()

    Args:
        annoFile (str): path to the file to read info from
        relToRef (tuple/None, optional): control the location of the result bed file.
                                         Should be tuple (start, end), relative to the start
                                         site of the feature. Bothe start and end can be None, which
                                         will refer to the start and end of each feature automatically.
                                         NOTE the location is always relative to the strand of the
                                         selected feature.
                                         Defaults to None. Then the result will be the feature itself
        ref (str, 'start'/'end'): Control the behaviour of calculating position when relToRef != None.
                                  Will select either 'start' of the target feature or 'end' of it.
        bedType (str, optional): feature type of target. Defaults to 'gene'.
        idQualifier (str, optional): qualifier name to parse as ID. Defaults to 'locus_tag'.
        withSynonym (str/None, optional): for M145 gbk file, it can be 'gene_synonym'. Defaults to None. 
                                          NOTE the gene name will be polluted!
        bedFile (str/""/None, optional): path to output bed file. Defaults to None. Then
                                      the bed file will not be written. If bedFile == "", then
                                      the bed file will appear in the same directory of the 
                                      input file.
    """
    filePath, genomeName = os.path.split(annoFile)
    trackName = os.path.splitext(genomeName)[0]
    if isinstance(bedFile, type(None)):
        bedFileIO = NamedTemporaryFile()
        bedFile = bedFileIO.name
    if bedFile == "":
        bedFile = f'''{trackName}{
                   "" if isinstance(relToRef, type(None)) else
                        f"_ss_{relToRef[0]}_{relToRef[1]}"
                    }{f"_{withSynonym}" if not isinstance(withSynonym, type(None)) else
                        "_SCO"}.bed'''
        bedFile = os.path.join(filePath, bedFile)
    ext = os.path.splitext(annoFile)[1].lower()
    if ext in ['.gbk','.gb']:
        fullRec = SeqIO.parse(annoFile, 'genbank')
    elif ext in ['.gff','.gtf']:
        fullRec = GFF.parse(annoFile)
    else:
        raise ValueError(f'File extension {annoFile} not recognisable')
    bedLines = []
    listID = []
    dictDuplicats = {}
    chroms = []
    for myseq in fullRec:
        chrom = myseq.name
        chroms.append(chrom)
        for feature in myseq.features:
            if feature.type == bedType:
                id = feature.qualifiers[idQualifier][0]
                if not isinstance(withSynonym, type(None)):
                    if withSynonym in feature.qualifiers:
                        synonyms = feature.qualifiers[withSynonym][0].split('; ')
                        for synonym in synonyms: # only use the first short name
                            if len(synonym) <= 6:
                                id = f'{id}_{synonym}'
                                break

                start = feature.location.start.position
                end = feature.location.end.position
                strand = ('+' if feature.strand == 1 else '-')
                score = '0'
                if id in listID:
                    if id in dictDuplicats:
                        dictDuplicats[id] += 1
                    else:
                        dictDuplicats[id] = 0
                    id = f'{id}_cp{dictDuplicats[id]}'
                else:
                    listID.append(id)
                if not isinstance(relToRef, type(None)):
                    if isinstance(relToRef[0], type(None)):
                        relToRef = (0, relToRef[1])
                    if isinstance(relToRef[1], type(None)):
                        relToRef = (relToRef[0], len(feature))
                    if ref == 'start':
                        position = (start if strand == "+" else end)
                    elif ref == 'end':
                        position = (end if strand == "+" else start)
                    else:
                        raise ValueError(f'ref needs to be "start" or "end", not {ref}')
                    start = position + (relToRef[0] if strand == '+' else -relToRef[1])
                    end = position + (relToRef[1] if strand == '+' else -relToRef[0])
                bedLines.append(f'{chrom}\t{start}\t{end}\t{id}\t{score}\t{strand}')

    bedLines.sort(key=lambda l: int(l.split('\t')[1]))
    bedLines.insert(0, f'track name={trackName}')
    with open(bedFile, 'w') as f:
        for line in bedLines:
            f.write(f'{line}\n')
    bedDf = parseBedFile(bedFile)
    return bedDf


def writeBed(bedDf, fn, trackName='track', chromName=None):
    """reverse parseBedFile

    Args:
        bedDf (pd.DataFrame): return value of parseBedFile()
        fn (str): path of bed file to write
        trackName (str, optional): will write into the header of bed file. Defaults to 'track'.
        chromName (str/None, optional): will be used in all lines if set, else will try to parse.                
            # Will always parse and return chroms from file, but:
            # if chromName is set, the output file will follow.

    Returns:
        chroms: list of unique chromosomes. Could be empty list if not found. Even if you set chrom 
                parameter.
        fn: path of bed file that has been written
    """
    with open(fn, 'w') as f:
        f.write(f'track name={trackName}\n')
        cols = [c.lower() for c in bedDf.columns]
        for c in ['start', 'end', 'chromosome', 'score', 'strand', 'chrom']:
            try:
                cols.remove(c)
            except:
                pass

        parseChrom = False
        chIdx = [c for c in bedDf.columns if 'chrom' in c.lower()] # parse chromosome name
        if len(chIdx) > 0:
            chIdx = chIdx[0]
            parseChrom = True
            chroms = [] # will store parsed chromosome names.
        else:
            chroms = ['None']
            
        for id, row in bedDf.iterrows():
            start = str(int(row['start']))
            end = str(int(row['end']))
            if parseChrom:
                chrom = row[chIdx]
                if chrom not in chroms:
                    chroms.append(chrom)
            else:
                chrom = chroms[0]
            try:
                score = str(row['score'])
            except:
                score = '0'
            try:
                strand = row['strand']
            except:
                strand = '+'
            el = [chrom, start, end, id, score, strand]
            # if chromName is set, the output file will follow, but return values stays the same
            if not isinstance(chromName, type(None)):
                el[0] = chromName
            if len(cols) > 0:
                el.extend(str(e) for e in row[cols])
            for s, e in zip(["chrom", "start", "end", "id", "score", "strand"], el):
                if isinstance(e, type(None)):
                    print(s, e)
                    raise ValueError
            l = '\t'.join(el) + '\n'
            f.write(l)
    return chroms, fn 

def getIntersection(featDf, peakDf, intersect=True, overlap=0.5):
    """
    Calculate intersection using bedtools.
    peakDf = parsePeakFile("/mnt/c/Users/duc/Documents/tmp_works/ChIP-25h_peaks.xls")
    featDf = annotationToBed("/mnt/c/Users/duc/Documents/tmp_works/M145.gb", relToRef=(-250,50))
    tFeatDf, tPeakDf, tFeatOlDf = getIntersection(featDf, peakDf, intersect=True, overlap=0.5)

    Args:
        featDf (pd.DataFrame): return value of annotationToBed()/parseBedFile()
        peakDf (pd.DataFrame): return value of parsePeakFile. columns = [start, end, length, score, strand] 
        intersect (bool, optional): whether to calculate intersect. Set to False will give NONE-intersection
                                    features. Defaults to True.
        overlap (float, optional): (0,1] overlap precentage of a PEAK to feature. Defaults to 0.5.

    Returns:
        tFeatDf, tPeakDf, tFeatOlDf:
            tFeatDf is the (none-)overlapping features, subset of genome features
            tPeakDf is the (none-)overlapping peaks, subset of peakDf
            tFeatOlDf is the overlapping features, including info of which feature overlaps with
                      which peak, and the length of overlap.
    """

    peakBed = NamedTemporaryFile()
    featBed = NamedTemporaryFile()
    # write peakBed:
    chromsf, _ = writeBed(featDf, featBed.name)
    chromsp, _ = writeBed(peakDf, peakBed.name, chromName=chromsf[0])
    if len(chromsp) > 1: # multiple chroms, will only work if
        # chromsf and p have common elements. Else raise error
        if any(p in chromsf for p in chromsp): # at least one common chromos
            if chromsf != chromsp:
                logger.warning(f'[intersection] not all chroms match eachother, \
                    feature: {chromsf}, peak: {chromsp}. Will only compare matched chromosome(s)')
            writeBed(peakDf, peakBed.name) # write again with its own chromosome names
        else:
            raise ValueError(f'peakDf contains multiple chromosomes: "{chromsp}" \
                    which have no intersection from annotation file: "{chromsf}".\
                    this will not make a intersection. Please check your input.')
    elif len(chromsf) > 1: # when len(chromsp) == 1
        logger.warning(f'[intersection] not all chroms match eachother, \
            feature: {chromsf}, peak: {chromsp}. Will only compare features on {chromsf[0]}.')
    

    peak = pybedtools.BedTool(peakBed.name)
    feat = pybedtools.BedTool(featBed.name)

    pfTargetFeat = feat.intersect(peak, v=not intersect, F=overlap, wa=intersect)
    pfTargetPeak = peak.intersect(feat, v=not intersect, f=overlap, wa=intersect)
    pfTargetFeatOl = feat.intersect(peak, F=overlap, wo=True)
    # https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html

    #with open(pfTargetFeat.fn, 'r') as f:
    #    for i, l in enumerate(f):
    #        print(l)
    #        if i == 5 :
    #            break
    ## the output of intersect will have an extra column of length. This is not standard
    ## bed format that I decide to not cover in parseBedFile(). That function wil calculate
    ## length and added to the resulted dataframe. As a result, there will be
    ## duplicated length column in the returned dataframe.
    tFeatDf = parseBedFile(pfTargetFeat.fn)
    tPeakDf = parseBedFile(pfTargetPeak.fn)
    tFeatOlDf = parseBedFile(pfTargetFeatOl.fn)

    # remove duplicates
    tFeatDf = tFeatDf[~tFeatDf.index.duplicated()]
    tPeakDf = tPeakDf[~tPeakDf.index.duplicated()]
    tFeatOlDf = tFeatOlDf[~tFeatOlDf.index.duplicated()]

    peakBed.close()
    featBed.close()
    pybedtools.cleanup()
    #print(featDf.columns)
    #print(peakDf.columns)

    return tFeatDf, tPeakDf, tFeatOlDf


def mwuTest(expList, ctrList, dataSeries):
    """[summary]
    - Mann-Whitney U test, to report
        - A measure of the central tendencies of the two groups (means or medians; since the Mann–Whitney U test is an ordinal test, medians are usually recommended)
        - The value of U (perhaps with some measure of effect size, such as common language effect size or rank-biserial correlation).
        - The sample sizes
        - The significance level

    Args:
        expList (-): [description]
        dataSeries ([type]): [description]

    Returns:
        [type]: [description]
    """
    expData = dataSeries[expList].dropna()
    ctrData = dataSeries[ctrList].dropna()
    n1, n2 = len(expList), len(ctrList)
    
    u1, p = mannwhitneyu(expData, ctrData, alternative='two-sided')
    u2 = n1*n2 - u1
    u1g, pg = mannwhitneyu(expData, ctrData, alternative='greater')
    u2g = n1*n2 - u1g
    u1l, pl = mannwhitneyu(expData, ctrData, alternative='less')
    u2l = n1*n2 - u1l
    if pg <= 0.05 and pl > 0.05:
        u1, u2, p, alt = (u1g, u2g, pg, 'greater')
    elif pl <= 0.05 and pg > 0.05:
        u1, u2, p, alt = (u1l, u2l, pl, 'less')
    else:
        alt = 'two-sided'
             
    return u1, u2, p, alt, n1, n2

def plotDistDiff(targetList, dataDf, prop, controlList=None, transform=None,
                 figsize=(6,5), log_scale=True, xlims=None, title=''):
                 # Transform needs to be a function.
    # TODO add option to add scatter plot
    # TODO add option to select orientation of the plot

    # Data preparation
    # Remove empty data
    ha = calHash(targetList, dataDf, prop, controlList, transform,
                 figsize, log_scale, xlims, title)
    emptyIdx = dataDf[dataDf[prop].isna() | (~np.isfinite(dataDf[prop]))].index
    dataDf = dataDf.drop(emptyIdx, axis=0)
    # Transform to dataframe in order to add a "hue" column for seaborn to plot
    if isinstance(dataDf, pd.core.series.Series):
        logger.info(f'Convert input Series to DataFrame, {prop} ignored, using {dataDf.name}.')
        prop = dataDf.name
        dataDf = pd.DataFrame(dataDf)
    # Remove items in input list that are not in dataDf
    # targetList
    try:
        assert all(i in dataDf.index for i in targetList)
    except:
        newT = [i for i in targetList if i in dataDf.index]
        excluded = [i for i in targetList if i not in dataDf.index]
        logger.info(f'Items in {excluded} are not found in data table (or is infinite).')
        targetList = newT
        assert len(targetList) > 0
    # controlList
    if isinstance(controlList, type(None)):
        controlList = [i for i in dataDf.index if i not in targetList]
    else:
        try:
            assert all(i in dataDf.index for i in controlList)
        except:
            newC = [i for i in controlList if i in dataDf.index]
            excluded = [i for i in controlList if i not in dataDf.index]
            logger.info(f'Items in {excluded} are not found in data table (or is infinite).')
            controlList = newC
            assert len(controlList) > 0
    # Check property    
    assert prop in dataDf.columns, f'{prop} not in column of peakDf {dataDf.columns}'
    logger.info(f'{dataDf.columns.to_list()} avaliable, using "{prop}"')
    # Transform data
    if not isinstance(transform, type(None)):
        trans = transform(dataDf[prop])
        dataDf = dataDf.loc[trans.index, :]
        dataDf.loc[:,prop] = trans

    # Mann Whitney U test
    u1, u2, p, alt, n1, n2 = mwuTest(targetList, controlList, dataDf[prop])

    # Calculate mean and median
    tMean   = dataDf.loc[targetList, prop].mean()
    tMedian = dataDf.loc[targetList, prop].median()
    cMean   = dataDf.loc[controlList, prop].mean()
    cMedian = dataDf.loc[controlList, prop].median()
    
    # add data column as "hue" for plotting
    dataName = 'data'
    histData = dataDf.copy()
    histData.loc[:,dataName] = ['exp' if i in targetList else 'ctr' for i in histData.index]

    # Init a figure with two subplots
    pname = f'{title}{"_" if title!="" else ""}{prop} distribution_{ha}'
    plt.close(pname)
    fig = plt.figure(
        constrained_layout=True,
        figsize=figsize, num=pname)
    figSpec = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=[1, 4])
    boxPlot = fig.add_subplot(figSpec[0])
    histPlot = fig.add_subplot(figSpec[1], sharex=boxPlot)
    histPlotT = histPlot.twinx()

    # PLOTTING
    # Boxplot
    sns.boxplot(data=histData, x=prop, y=dataName, order=['exp', 'ctr'], ax=boxPlot, showfliers=True, 
                flierprops=dict(markersize=2, alpha=0.5), width=0.5)
    if log_scale:
        boxPlot.set_xscale('log')
        histData = histData.loc[histData[prop]!=0, :] # else histplot will complain
    # Remove x ticks from box plot
    boxPlot.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    boxPlot.set_xlabel('')
    # Set x lims, this will apply to both subplots (sharex)
    if not isinstance(xlims, type(None)):
        boxPlot.set_xlim(xlims)
    # Histplot (distribution plot)
    #sns.histplot(histData, x=prop, hue=dataName, hue_order=['exp', 'ctr'], ax=disPlot, kde=True)
    # Plot each element seperately to "twinx" axes.
    cPalette = dict(exp=sns.color_palette()[0], ctr=sns.color_palette()[1])
    sns.histplot(histData, x=prop, hue=dataName, hue_order=['ctr'], ax=histPlot,
                 palette=cPalette, kde=True, legend=False)
    sns.histplot(histData, x=prop, hue=dataName, hue_order=['exp'], ax=histPlotT,
                 palette=cPalette, kde=True, legend=False)
    # Move 'exp' label to the left while keep the 'exp' ploted on top of 'ctr'
    histPlot.yaxis.tick_right()
    histPlot.yaxis.set_label_position('right')
    histPlotT.yaxis.tick_left()
    histPlotT.yaxis.set_label_position('left')
    histPlot.set_ylabel('Count Ctr')
    histPlotT.set_ylabel('Count Exp')
    
    # Make legend because sns cannot generate legend when the plots are seperate.
    #legends = [
    #    Patch(facecolor=cPalette['exp'],
    #          ec='k'),
    #    Patch(facecolor=cPalette['ctr'],
    #          ec='k')
    #]
    #histPlot.add_artist(Legend(histPlot, handles=legends, labels=['exp', 'ctr']))
    # Show mean and median on plot as legend
    mmLegs = [
        Patch(
            facecolor=cPalette['exp'], edgecolor='k',
        ),
        Patch(
            facecolor=cPalette['ctr'], edgecolor='k',
        ),
    ]
    tMeanStr = f'{tMean:.1f}'
    tMedianStr = f'{tMedian:.1f}'
    cMeanStr = f'{cMean:.1f}'
    cMedianStr = f'{cMedian:.1f}'
    if any([x == '0.0' for x in [tMeanStr, tMedianStr, cMeanStr, cMedianStr]]):
        tMeanStr = f'{tMean:.3f}'
        tMedianStr = f'{tMedian:.3f}'
        cMeanStr = f'{cMean:.3f}'
        cMedianStr = f'{cMedian:.3f}'

    mmLabels = [f'Exp mean {tMeanStr}\nExp median {tMedianStr}',
                f'Ctr mean {cMeanStr}\nCtr median {cMedianStr}']
    histPlot.add_artist(Legend(histPlot, mmLegs, mmLabels))

    # Show p value on plot
    pstr = f'MWU exp(U={u1},n={n1})/ctr(U={u2},n={n2})\nexp {alt} $p$-value = '
    logger.info(pstr + str(p))
    if 0.01 < p <= 0.05:
        sig = "*"
        pstr += f'{p:.4f}'
    elif 0.001 < p <= 0.01:
        sig = '**'
        pstr += f'{p:.5f}'
    elif 0.0001 < p <= 0.001:
        sig = '***'
        pstr += f'{p:.6f}'
    elif p <= 0.0001:
        sig = '****'
        pstr += f'{p:.3e}'
    else:
        sig = 'No significance'
        pstr += f'{p:.4f}'
    boxPlot.text(0,1, pstr, figure=fig, ha='left', va='bottom', transform=boxPlot.transAxes, fontsize=8)

    fig.suptitle(pname)
    # calculate hash and save figure
    outputPath = 'Plots'
    os.makedirs(outputPath, exist_ok=True)
    figFile = os.path.join(outputPath, pname+'.svg')
    if os.path.isfile(figFile):
        logger.info(f'Plot exists: {figFile}')
    else:
        logger.info(f'Save plot at {figFile}')
        fig.savefig(figFile)
    

    plt.show()


def readCompRes(compResultFile):

    #     newCompDf = pd.Dataframe()
    #     newCompTablePath = 'dataTables/transformed/msstats_proposed_comparisonResult_readable.tsv'
    """
    allComparisons - dict of {'comp1':DF, 'comp2':DF,...}

    DF.columns = ['log2FC', 'SE', 'pvalue', 'adj.pvalue', 'ImputationPercentage']
    for MSstats result
    DF.columns = ['baseMean', 'log2FC', 'SE', 'pvalue', 'adj.pvalue']
    for DESeq2 result
    """
    compResultDf = pd.read_csv(compResultFile, sep='\t', index_col=0)

    if 'ImputationPercentage' in compResultDf.columns:
        # MSstats output
        # "Label"	"log2FC"	"SE"	"Tvalue"	"DF"	"pvalue"	"adj.pvalue"	"issue"	"MissingPercentage"	"ImputationPercentage"
        targetCols = ['log2FC', 'SE', 'pvalue', 'adj.pvalue', 'ImputationPercentage']
    elif 'baseMean' in compResultDf.columns:
        # DESeq2 result
        # "Label" is added while combining all data
        # "baseMean"	"log2FoldChange"	"lfcSE"	"pvalue"	"padj"
        # for shrinked data, note the column "lfcSE" is actually "posterior SD"
        # OR
        # "baseMean"	"log2FoldChange"	"lfcSE"	"stat"	"pvalue"	"padj"
        # for non-shrinked data
        targetCols = ['log2FC', 'SE', 'pvalue', 'adj.pvalue', 'baseMean']
    else:
        raise ValueError(f'"baseMean" and "Imputation Percentage" not found in \n{compResultDf.columns}')

    allCompResults = {}
    comparisons = list(compResultDf.Label.unique())
    comparisons.sort()
    for c in comparisons:
        compData = compResultDf[compResultDf.Label == c]
        compData = compData.loc[:, targetCols]
        allCompResults[c] = compData

    return allCompResults



def plotCompRes(geneList, compDict, comp, target, controlList=None, title='', **kwargs):
    logger.info(f'{list(compDict.keys())} avaliable in dict, using "{comp}".')
    compDf = compDict[comp]
    assert target in compDf.columns
    plotDistDiff(geneList, compDf, target, controlList=controlList, title=comp+title, **kwargs)


def countDiff(threshA, threshB, sample, dfA, dfB, featIdx, printOnly=False):
    protList = []
    if threshA[1] >= max(dfA[sample]):
        pl = dfA[(dfA[sample]>=threshA[0])&(dfA[sample]<=threshA[1])].index.to_list()
    else:    
        pl = dfA[(dfA[sample]>=threshA[0])&(dfA[sample]<threshA[1])].index.to_list()
    for p in pl:
        ps = p.split(';')
        if len(ps) > 1:
            continue
        protList.append(p)
    if threshB[1] >= max(dfB[sample]):
        geneList = dfB[(dfB[sample]>=threshB[0])&(dfB[sample]<=threshB[1])].index.to_list()
    else:
        geneList = dfB[(dfB[sample]>=threshB[0])&(dfB[sample]<threshB[1])].index.to_list()
    idList = featIdx

    intProt = idList.intersection(set(protList))
    intGene = idList.intersection(set(geneList))
    intPG = intProt.intersection(intGene)
    difp = idList.symmetric_difference(set(protList))
    difg = idList.symmetric_difference(set(geneList))
    difPG = difp.intersection(difg)

    averageLenInt = max(((len(intGene)+len(intProt))/2), 1)
    averageLenDif = max(((len(difp)+len(difg))/2), 1)
    comRateInt = len(intPG)/averageLenInt
    comRateDif = len(difPG)/averageLenDif

    if printOnly:
        textA = [f"{x:.1e}" for x in threshA]
        textA = f'[{textA[0]}, {textA[1]})'
        textB = [f"{x:.1e}" for x in threshB]
        textB = f'[{textB[0]}, {textB[1]})'
        print(f'{len(intProt)} proteins between {textA} and bound by Gbn.')
        print(f'{len(intGene)} genes between {textB} and bound by Gbn.')
        print(f'{len(intPG)} are common ({comRateInt:0.0%}).')
        print(f'As a back ground:')
        print(f'{len(difp)} proteins between {textA} and not bound by Gbn.')
        print(f'{len(difg)} genes between {textB} and not bound by Gbn.')
        print(f'{len(difPG)} are common for different ids ({comRateDif:0.0%}).')
    else:
        countDict = dict(
            commonIdsRate=comRateInt,
            diffIdsRate=comRateDif,
            commonA=intProt,
            commonB=intGene,
            commonAB=intPG,
            diffA=difp,
            diffB=difg,
            diffAB=difPG
        )
        return countDict


def calBin(data, n=15, log=False):
    # should be all positive or zero data
    data = np.asarray(data)
    data = data[~np.isnan(data)]
    hasZero = (min(data) == 0)
    if hasZero:
        data[data==0] = min(data[data>0])/2

    bins = []
    if log:
        t = lambda x: np.log10(x)
        r = lambda x: np.power(10, x)
    else:
        t = lambda x: x
        r = lambda x: x
    binEdges = np.linspace(t(min(data)), t(max(data)), n+1)
    for i in range(n):
        bins.append((r(binEdges[i]), r(binEdges[i+1])))
    if hasZero:
        bins[0] = (0, bins[0][1])
    if np.isnan(bins[0][1]):
        print(data)
        raise ValueError('data not correct?')
    #print(len(bins))
    return bins

def plotDiffCounts(dfA, dfB, sample, targets,
                   n=25, figsize=(6,5),
                   title='', alabel='a', blabel='b',
                   printN=None, bins=None,
                   ):
    if isinstance(bins, type(None)):
        binsA = calBin(dfA[sample], n=n, log=True)
        binsB = calBin(dfB[sample], n=n, log=True)
    else:
        binsA, binsB = bins
    #binEdgesA = list(set(np.asarray(binsA).flatten()))
    #binEdgesB = list(set(np.asarray(binsB).flatten()))
    #binEdgesA.sort()
    #binEdgesB.sort()

    plotData = pd.DataFrame()

    for i, (pb, gb) in enumerate(zip(binsA, binsB)):
        c = countDiff(pb, gb, sample, dfA, dfB, targets)
        plotData[i] = pd.Series([c['commonIdsRate'], c['diffIdsRate']],
                            index = ['pCommon', 'pDiff'])
    if isinstance(printN, type(None)):
        printN = int(n/2)
    print(f'{printN}th(outof {n}) bin statment:')
    countDiff(binsA[printN-1], binsB[printN-1], sample, dfA, dfB, targets, printOnly=True)
    
    plotData = plotData.T
    
    fname = f'plotDiffCounts_{sample}{title}'
    plt.close(fname)
    fig = plt.figure(
        constrained_layout=True,
        figsize=figsize, num=fname)
    figSpec = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=[1, 4])
    ax = fig.add_subplot(figSpec[1])

    ax.plot(plotData.index, plotData.pDiff, c='C1', label='common')
    ax.fill_between(plotData.index, plotData.pDiff, color='C1',
                    linewidth=0, alpha=0.3)

    tax = ax.twiny()
    tax.plot(plotData.index, plotData.pCommon,c='C0',  label='difference')
    tax.fill_between(plotData.index, plotData.pCommon, color='C0',
                     linewidth=0, alpha=0.3)

    legends = [
        Line2D([],[], c='C0', label='common'),
        Line2D([],[], c='C1', label='difference')
    ]
    ax.legend(handles=legends)
    tax.xaxis.tick_bottom()
    tax.xaxis.set_label_position('bottom')
    ax.set_xlabel('Percentage')

    # change y ticks to %
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))

    # change x ticks to corresponding bin value
    ax.set_xticks(plotData.index)
    tax.set_xticks(plotData.index)
    ax.set_xticklabels(
        [f'{(a+b)/2:.1e}' for a,b in binsA],
        rotation=30,
        fontsize='small'
        )
    ax.set_xlabel(alabel)
    tax.spines.bottom.set_position(('outward', 55))
    tax.set_xticklabels(
        [f'{(a+b)/2:.1e}' for a,b in binsB],
        rotation=30,

        fontsize='small'
    )
    tax.set_xlabel(blabel)

    # plot distribution
    distAx = fig.add_subplot(figSpec[0])
    tdistAx = distAx.twinx()

    distData = pd.DataFrame()
    for i, (ta, tb) in enumerate(zip(binsA, binsB)):
        if i < len(binsA)-1:
            na = sum((ta[0]<=dfA[sample])&(dfA[sample]<ta[1]))
            nb = sum((tb[0]<=dfB[sample])&(dfB[sample]<tb[1]))
            distData[i] = pd.Series([na, nb], index=['A','B'])
        else:
            na = sum((ta[0]<=dfA[sample])&(dfA[sample]<=ta[1]))
            nb = sum((tb[0]<=dfB[sample])&(dfB[sample]<=tb[1]))
            distData[i] = pd.Series([na, nb], index=['A','B'])
    #print(distData.shape)
    distData = distData.T

    distAx.bar(plotData.index, distData['A'], color='C2', alpha=0.3, width=1, align='edge')
    tdistAx.bar(plotData.index, distData['B'], color='C3', alpha=0.3, width=1, align='edge')

    distAx.set_yticks([])
    tdistAx.set_yticks([])

    distAx.set_xlim(ax.get_xlim())
    legendsDist=[
        Patch(facecolor='C2', alpha=0.3,
              label=f'{alabel} dist'),
        Patch(facecolor='C3', alpha=0.3,
              label=f'{blabel} dist'),
    ]
    distAx.legend(handles=legendsDist)

    distAx.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    distAx.set_xlabel('')
    fig.suptitle(fname)

    plt.show()

    return binsA, binsB
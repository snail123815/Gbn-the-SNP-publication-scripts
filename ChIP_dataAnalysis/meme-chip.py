# samtools convert sam to bam
import subprocess
import os

fileDict = {}
path = '.'

for file in os.listdir(path):
    if file.endswith('_peaks.fa'):
        filePath = os.path.join(path, file)
        source = os.path.splitext(file)[0]
        fileDict[source] = filePath


memedb = 'collectf.meme'


def runMemeChip(inputFile, db, outputDir, mod='anr'):
#    if os.path.isdir(outputDir):
#        print(f'{outputDir} already exists, skip')
#        return
    args = ['meme-chip',
            '-oc', outputDir,
            '-db', db,
            '-meme-mod', mod,  # [oops|zoops|anr]
            inputFile,
            '-meme-p', '4'
            ]
    print('Running...')
    print(' '.join(args))
    with subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) as p:
        for line in p.stdout:
            print(line.decode('utf-8'), end='')
    print('*' * 80)


output = 'MEME_CHIP_Motif'
os.makedirs(output, exist_ok=True)
mod='anr'
for source, file in fileDict.items():
    inputFile = os.path.join(path, file)
    outputDir = os.path.join(
        output, f'{os.path.splitext(file)[0]}_{mod}')
    runMemeChip(inputFile, memedb, outputDir, mod = mod)

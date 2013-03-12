import numpy.random;
import math;
import sys;

'''
if(len(sys.argv) != 2):
    print "No input file given! Quitting..."
    exit(1);
else:
    FileName = sys.argv[1];
'''

numXPixels = 0;
numYPixels = 0;

for linetext in open("InputScript.txt"):
    line = linetext.split('\t');
    if line[0] == "CCDNumXPixels":
        numXPixels = int(line[1]);
    if line[0] == "CCDNumYPixels":
        numYPixels = int(line[1]);


if (numXPixels*numYPixels) > 9000000: #9 million
    print "Too many pixels to manage. Disable this check by editing this script!"
    exit(1);


CCDArray = []

for x in range(numXPixels):
    CCDArray.append([]);
    for y in range(numYPixels):
        CCDArray[x].append(float(0))



#'DiffractResultsPostPinhole.txt'

Data = [line.split('\t') for line in open("DiffractResultsPostPinhole.txt")];

nPileupEvents = 0;

for line in Data:
    xPixel = int(line[0])
    yPixel = int(line[1])
    if (0 <= xPixel < numXPixels) and (0 <= yPixel < numYPixels):
        if( CCDArray[xPixel][yPixel] > 0):
            nPileupEvents += 1
        CCDArray[xPixel][yPixel] = CCDArray[xPixel][yPixel] + float( line[2] )


Data = [line.split('\t') for line in open("FluoResultsPostPinhole.txt")];
for line in Data:
    xPixel = int(line[0])
    yPixel = int(line[1])
    if (0 <= xPixel < numXPixels) and (0 <= yPixel < numYPixels):
        if( CCDArray[xPixel][yPixel] > 0):
            nPileupEvents += 1
        CCDArray[xPixel][yPixel] = CCDArray[xPixel][yPixel] + 1.71 #Ta M shell fluorescence energy

numpy.random.seed(42)

'''
    
    delta(E) = 2.355*w*sqrt( r^2 + F*E/w  )
    
    w = eh pair energy, 3.68ev
    r = readout noise (approx 3?)
    F = fano factor 0.117
    E = xray energy in ev
    
    remove 2.355 for sigma instead of FWHM
    
'''

w = 3.65; #in eV not keV
r = 3;
F = 0.117;


Outfile = open('PixelatedCCDSpectrum.txt', 'w')
    
for x in range(numXPixels):
    CCDArray.append([]);
    for y in range(numYPixels):
        if( CCDArray[x][y] > 0.0):
            Energy = CCDArray[x][y]*1000; #energy in ev
            Stdev = w*math.sqrt( r*r + (F*Energy)/w );
            RandomE = numpy.random.normal( Energy, Stdev);
            Outfile.write( str(x) + '\t' + str(y) + '\t' + str(RandomE/1000.0) + '\n');
#Outfile.write( str(x) + '\t' + str(y) + '\t' + str(CCDArray[x][y]) + '\n')


Outfile.close();

print "Number of Pileup events: ", nPileupEvents;

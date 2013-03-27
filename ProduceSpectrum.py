import numpy.random;
import math;
import sys;

if(len(sys.argv) != 2):
	print "No input file given! Quitting..."
	exit(1);
else:
	FileName = sys.argv[1];


numXPixels = 0;
numYPixels = 0;

for linetext in open("InputScript.txt"):
	line = linetext.split('\t');
	if line[0] == "CCDNumXPixels":
		numXPixels = int(line[1]);
	if line[0] == "CCDNumYPixels":
		numYPixels = int(line[1]);



QEData = [line.split() for line in open("QuantEfficiency.txt")];

for lineIndex in range(len(QEData)):
	for itemIndex in range(len(QEData[lineIndex])):
		QEData[lineIndex][itemIndex] = float(QEData[lineIndex][itemIndex])
	
def FindQuantEfficiency( Energy ):
		
	if Energy <= QEData[0][0]:
		return QEData[0][1]
	
	if Energy >= QEData[len(QEData)-1][0]:
		return QEData[len(QEData)-1][1]
	
	Index = 0;
	
	for lineIndex in range(len(QEData)):
		if QEData[lineIndex][0] <= Energy:
			Index = lineIndex;
			
	if Index != (len(QEData)-1):
		binWidth = QEData[Index+1][0] - QEData[Index][0]
		weight = (Energy-QEData[Index][0])/binWidth
		return QEData[Index][1] + ((QEData[Index+1][1] - QEData[Index][1])*weight);
	else:
		exit(1)

'''
for i in range(801):
	E = 2.0+float(i)*(6.0/500.0);
	print str(E) + "\t" + str(FindQuantEfficiency(E))
	
exit(0)
'''


#'DiffractResultsPostPinhole.txt'

#DiffractData = [line.split('\t') for line in open(FileName)];

#print "Input Photons: ", len(DiffractData)

numpy.random.seed(42)

'''

delta(E) = 2.355*w*sqrt( r^2 + F*E/w  )

w = eh pair energy, 3.68ev
r = readout noise (approx 3?)
F = fano factor 0.117
E = xray energy in ev

remove 2.355 for sigma instead of FWHM

'''

Outfile = open('CCDSpectrum.txt', 'w')

w = 3.65; #in eV not keV
r = 3;
F = 0.117;

nOutputPhotons = 0;
nEfficiencyFailed = 0
for line in open(FileName):
	XRay = line.split('\t')
	if (0 <= int(XRay[0]) < numXPixels) and (0 <= int(XRay[1]) < numYPixels):
		Energy = float(XRay[2])*1000; #energy in ev
		if(Energy < 3000.0):
			continue; #throw away secondary fluorescence
		'''
		if FindQuantEfficiency(Energy) > numpy.random.uniform(0.0,1.0):
			nEfficiencyFailed = nEfficiencyFailed + 1
			continue;
		'''
		Stdev = w*math.sqrt( r*r + (F*Energy)/w );
		RandomE = numpy.random.normal( Energy, Stdev);
		Outfile.write( XRay[0] + '\t' + XRay[1] + '\t' + str(RandomE/1000.0) + '\n');
		nOutputPhotons = nOutputPhotons+1;
		
print "Number lost to QE: ", nEfficiencyFailed
print "Output Photons (without secondary fluorescence & CCD Misses & QE Test): ", nOutputPhotons

Outfile.close();
	

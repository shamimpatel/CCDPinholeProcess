import sys;
import math;
import numpy;
import copy;

bWriteToFile = False;


#1st arg radius
#2nd arg distance

PinholeRadius = 0;
PinholeDistance = 0;

if(len(sys.argv) == 3):
	bWriteToFile = True;
	PinholeRadius = float(sys.argv[1])
	PinholeDistance = float(sys.argv[2])

def PickFirstElement( array ):
	return array[0]


RangeData = [line.split('\t') for line in open("DiffractionPeakLimits.txt")];

Ranges = [] 

for line in RangeData:
	Ranges.append([float(line[1]),float(line[2]),line[0],[],[]]) #[LowE, HighE, Name, xcoords, Evalue]

Ranges.sort( key=PickFirstElement );

OriginalRanges = copy.deepcopy(Ranges)


for i in range(len(Ranges)):
	
	if(i == (len(Ranges) - 1) ):
		lowEBound = (OriginalRanges[i-1][1] + OriginalRanges[i][0])/2.0;
		Ranges[i][0] = lowEBound;
		Ranges[i][1] = 12;
		break;
	
	if(i == 0):
		Ranges[i][0] = 2;
		HighEBound = (OriginalRanges[i][1] + OriginalRanges[i+1][0])/2.0
		Ranges[i][1] = HighEBound
		continue; #this sorts the lower energy bound for the lowest energy peak

	LowEBound = (OriginalRanges[i][0] + OriginalRanges[i-1][1])/2.0
	Ranges[i][0] = LowEBound
	HighEBound = (OriginalRanges[i][1] + OriginalRanges[i+1][0])/2.0
	Ranges[i][1] = HighEBound
	
	

	
#CCDData = [line.split('\t') for line in open("DiffractResultsPostPinhole.txt")];
CCDData = [line.split('\t') for line in open("CCDSpectrumDiffract.txt")];

w = 0.00365; #in keV
r = 3;
F = 0.117;


for Range in Ranges:
#Range = Ranges[2];
	for line in CCDData:
		if( Range[0] <= float(line[2]) <= Range[1] ): #if Energy is between min/max energy ranges for this peak
			Range[3].append(int(line[0])) #add it's xpixel coordinate and energy to the stack
			Range[4].append(float(line[2])) #Ener
	if(len(Range[3]) > 1 ):
		mean = numpy.mean(Range[4])
		width = numpy.std(Range[4])
		spectralwidth = w*math.sqrt( r*r + (F*mean)/w )
		print Range[2], "Mean:", mean, "StDev:", width, "Expected Width:", spectralwidth, "Ratio: ", width/spectralwidth;
	else:
		print Range[2], "Less than 2 points"




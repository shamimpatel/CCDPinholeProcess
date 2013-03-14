from scipy.stats import linregress
import numpy;
import sys;

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
	Ranges.append([float(line[1]),float(line[2]),line[0],[],[]]) #[LowE, HighE, Name, xcoords, ycoords]

Ranges.sort( key=PickFirstElement );

#for R in Ranges:
#	print R;

CCDData = [line.split('\t') for line in open("DiffractResultsPostPinhole.txt")];

Correlations = []


#for Range in Ranges:
Range = Ranges[2];
AvgCorrelation = 0;

for line in CCDData:
	if( Range[0] <= float(line[2]) <= Range[1] ):
		Range[3].append(int(line[0]))
		Range[4].append(float(line[2]))
	if(len(Range[3]) > 2 ):
		slope, intercept, r, p, stderr = linregress(Range[3],Range[4])
		#print Range[2] + ":\t" + str(slope) + "\t" + str(r)
		#Correlations.append(r)
		AvgCorrelation = r;
	else:
		AvgCorrelation = 0;
		#Correlations.append(0);


#AvgCorrelation = Correlations[0];  #numpy.mean( Correlations );

if( bWriteToFile == False ):
	print "Average Correlation:\t" + str(AvgCorrelation)
	exit(0)



Outfile = open('CorrelationData.txt', 'a+')
Outfile.write( str(PinholeRadius) + "\t" + str(PinholeDistance) + "\t" + str(AvgCorrelation) + "\t" + str(len(Range[3])) + "\n")
Outfile.close();
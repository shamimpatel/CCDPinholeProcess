'''
Created on Mar 26, 2013

@author: shamim
'''
import matplotlib.pyplot as plt

print "start"

DiffractEnergies = []
DiffractData = [line.split('\t') for line in open("DiffractResultsPostPinhole_normal.txt")];
for line in DiffractData:
	if (0 <= int(line[0]) < 2048) and (0 <= int(line[1]) < 2048):
		DiffractEnergies.append( float( line[2] ) );
plt.hist( DiffractEnergies, bins=400, histtype='step', label="No Blur, Normal")


DiffractEnergies2 = []
DiffractData = [line.split('\t') for line in open("DiffractResultsPostPinhole.txt")];
for line in DiffractData:
	if (0 <= int(line[0]) < 2048) and (0 <= int(line[1]) < 2048):
		DiffractEnergies2.append( float( line[2] ) );
plt.hist( DiffractEnergies2, bins=400, histtype='step', label="No Blur, Compressed")




BlurredEnergies = []
BlurredData = [line.split('\t') for line in open("CCDSpectrum_normal.txt")];
for line in BlurredData:
	BlurredEnergies.append( float( line[2] ) )
plt.hist( BlurredEnergies, bins=400, histtype='step', label="WithCCDBlur, Normal")	
	
BlurredEnergies2 = []
BlurredData = [line.split('\t') for line in open("CCDSpectrum.txt")];
for line in BlurredData:
	BlurredEnergies2.append( float( line[2] ) )
plt.hist( BlurredEnergies2, bins=400, histtype='step', label="WithCCDBlur, Compressed")


plt.xlabel( "Energy (keV)")
plt.ylabel( "Num Photons")


plt.legend();

plt.show()

print "Done!"
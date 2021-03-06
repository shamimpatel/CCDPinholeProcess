import os;

minDistance = 0.00001;
maxDistance = 4;
DistanceStep = 0.01

minRadius = 0.0
maxRadius = 0.02
RadiusStep = 0.0005


ProgressCounter  = 0;

Outfile = open('CorrelationData.txt', 'w')
Outfile.close();


#1st arg radius
#2nd arg distance
for DistanceCounter in range( int(( maxDistance - minDistance )/DistanceStep) ):
    for RadiusCounter in range( int(( maxRadius - minRadius )/RadiusStep) ):
        Distance = minDistance + DistanceStep*DistanceCounter;
        Radius = minRadius + RadiusStep*RadiusCounter;
        os.system( "./CCDPinholeProcessMain " + str(Radius) + " " + str(Distance) + " >/dev/null" )
        os.system( "python2.7 CalculateCorrelations.py " + str(Radius) + " " + str(Distance) + " >/dev/null" )
    Progress = 100.0 * float(DistanceCounter)/int(( maxDistance - minDistance )/DistanceStep)
    print str(Progress) + "%"
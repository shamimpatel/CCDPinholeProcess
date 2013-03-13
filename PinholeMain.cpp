#include <iostream>
#include "Vector.h"
#include "AbsorbCoeffData.h"
#include "XRay.h"

#include "CCD.h"
#include "Pinhole.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

typedef boost::mt11213b base_generator_type;

#include "FileReading.h"

using namespace std;

int main(int argc, char *argv[])
{

    ifstream datafile("InputScript.txt");
    if(datafile.is_open() == false)
    {
        cout << "Error: Failed to open InputScript.txt" << endl;
        exit(1);
    }
    std::map<std::string,std::string> InputData;
    AddToMapFromFile(datafile, InputData);
    datafile.close();
    
   
    CCD CCDCamera = GenerateCCDFromInputScript("InputScript.txt");
    
    Vector PinholeOrigin(1.39895,0,0.3497);
    Vector PinholeNormal(40,0,10);
    
    VectorFromMap("CCDNormal", InputData, PinholeNormal);
    
    PinholeNormal = PinholeNormal.Normalized();
    
    double PinholeDistance;
    DoubleFromMap("TempPinholeDistance", InputData, PinholeDistance);
        
    double PinholeRadius = 0.5;    
    DoubleFromMap("PinholeRadius", InputData, PinholeRadius);
    
    
    //optional override so that pinhole parameters can be modified without changing the input script
    //TODO: use switches/flags instead of just reading in order
    if(argc == 3)
    {        
        std::stringstream S1(argv[1]);
        S1 >> PinholeRadius;
        std::stringstream S2(argv[2]);
        S2 >> PinholeDistance;        
    }   
    
    
    PinholeOrigin = PinholeNormal*PinholeDistance;
    PinholeOrigin.Print();
        
    PinholePlane Pinhole( PinholeOrigin, PinholeNormal, PinholeRadius);

    base_generator_type generator(48u);
    boost::uniform_real<> uni_dist(0,1);
    boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);


    AbsorbCoeffData MuData( 1.0f, 15.0f, 2000);
    MuData.LoadData("FeAbsorbCoeff.txt");

    ofstream FluoResults2( "FluoResultsPostPinhole.txt" );

    ifstream FluoFile("AdvFluoResults.txt");
    string dataline;

    Vector FilterNormal = PinholeNormal;
    FilterNormal = FilterNormal.Normalized();

    double FilterThickness = 25000.0; //2.5 micron in A
    DoubleFromMap("FilterThickness", InputData, FilterThickness);

    float Energy = 1.710; //M shell fluorescence energy
    float AbsorbCoeff = MuData.GetAbsorbCoeffDataPoint(EnergyToWavelength(Energy));

    double RayLength;
    Vector IntersectPoint;
    int XPixel,YPixel;
    double XIntersect,YIntersect;

    while(getline(FluoFile, dataline, '\n'))
    {

        Vector Source(0,0,0);
        Vector Direction(0,0,0);

        stringstream linestream(dataline);
        linestream >> Source.x >> Source.y >> Direction.x >> Direction.y >> Direction.z;

        double CosTheta = fabs(Direction.Dot(FilterNormal)); //fabs for anti-parallel.

        double PathLength;

        if(CosTheta <= 0.001) //handle divide by zero.
        {
            PathLength = 0.0; //just say it passes through.
            //continue; //Just throw it away instead?? Depends on geometry
        }
        else
        {
            PathLength = FilterThickness/CosTheta;
        }

        double ProbTransmit = exp( -1.0f * AbsorbCoeff * PathLength);

        //Don't combine these. It very slightly boosts the fluorescence. This first one *must* fail (ie not transmit) in order to test for fluorescence.
        if(uni() < ProbTransmit)
        {
            /*double CCDIntersectX = Source.x + ((50.0-Source.z)/(Direction.z))*Direction.x;
            double CCDIntersectY = Source.y + ((50.0-Source.z)/(Direction.z))*Direction.y;*/

            if(CCDCamera.GetPixelRayCCDIntersect(Source,Direction,
                                                 RayLength,IntersectPoint,
                                                 XIntersect,YIntersect,
                                                 XPixel,YPixel) == true)
            {
                if(Pinhole.TestRayPinholeIntersect(Source, Direction, RayLength, IntersectPoint))
                {                
                    FluoResults2 << XPixel  << "\t" << YPixel << endl;
                }
            }
        }
        else if(uni() < 0.0063*0.5) //fluorescence yield approx 0.0063 (http://www.nist.gov/data/PDFfiles/jpcrd473.pdf)
        {
            //double CCDIntersectX = Source.x + ((50.0-Source.z)/(Direction.z))*Direction.x;
            //double CCDIntersectY = Source.y + ((50.0-Source.z)/(Direction.z))*Direction.y;
            //cout << "secondary" << endl;
            //FluoResults2 << CCDIntersectX  << "\t" << CCDIntersectY << endl;
            if(CCDCamera.GetPixelRayCCDIntersect(Source,Direction,
                                                 RayLength,IntersectPoint,
                                                 XIntersect,YIntersect,
                                                 XPixel,YPixel) == true)
            {
                //FluoResults2 << XPixel  << "\t" << YPixel << endl;
                if(Pinhole.TestRayPinholeIntersect(Source, Direction, RayLength, IntersectPoint))
                {
                    FluoResults2 << XPixel  << "\t" << YPixel << endl;
                }
            }
        }
    }

    FluoResults2.close();

    ifstream DiffractFile("AdvDiffractResults.txt");

    ofstream DiffractResults2( "DiffractResultsPostPinhole.txt" );

    while(getline(DiffractFile, dataline, '\n'))
    {
        Vector Source(0,0,0);
        Vector Direction(0,0,0);

        stringstream linestream(dataline);
        linestream >> Source.x >> Source.y >> Direction.x >> Direction.y >> Direction.z >> Energy;

        AbsorbCoeff = MuData.GetAbsorbCoeffDataPoint(EnergyToWavelength(Energy));

        double CosTheta = fabs(Direction.Dot(FilterNormal)); //fabs for anti-parallel.

        double PathLength;

        if(CosTheta <= 0.001) //handle divide by zero.
        {
            PathLength = 0.0; //just say it passes through.
            //continue; //Just throw it away instead?? Depends on geometry
        }
        else
        {
            PathLength = FilterThickness/CosTheta;
        }

        double ProbTransmit = exp( -1.0f * AbsorbCoeff * PathLength);

        if(uni() < ProbTransmit)
        {
            if(CCDCamera.GetPixelRayCCDIntersect(Source,Direction,
                                                 RayLength,IntersectPoint,
                                                 XIntersect,YIntersect,
                                                 XPixel,YPixel) == true)
            {
                if(Pinhole.TestRayPinholeIntersect(Source, Direction, RayLength, IntersectPoint))
                {
                    DiffractResults2    << XPixel  << "\t" << YPixel << "\t"
                                        << Energy << "\t" << Source.x << endl;
                }                
            }
        }
        else if(uni() < 0.0063*0.5) //fluorescence yield
        {
            if(CCDCamera.GetPixelRayCCDIntersect(Source,Direction,
                                                 RayLength,IntersectPoint,
                                                 XIntersect,YIntersect,
                                                 XPixel,YPixel) == true)
            {
                if(Pinhole.TestRayPinholeIntersect(Source, Direction, RayLength, IntersectPoint))
                {
                    DiffractResults2 << XPixel  << "\t" << YPixel << "\t" << 0.7 << endl;
                }
                //DiffractResults2 << XPixel  << "\t" << YPixel << "\t" << 0.7 << endl;
            }
        }
    }

    DiffractResults2.close();
    
    ofstream CCDBounds( "CCDBounds.txt" );
    
    int numXPixels, numYPixels;
    
    IntFromMap("CCDNumXPixels", InputData, numXPixels);
    IntFromMap("CCDNumYPixels", InputData, numYPixels);
    
    CCDBounds << "0\t0" << endl;
    CCDBounds << "0\t"  << numYPixels << endl;
    CCDBounds << numXPixels << "\t" << numYPixels << endl;
    CCDBounds << numXPixels << "\t0"    << endl;
    CCDBounds << "0\t0" << endl; //need to complete the loop
    
    CCDBounds.close();

    
    cout << "Done!" << endl;
}

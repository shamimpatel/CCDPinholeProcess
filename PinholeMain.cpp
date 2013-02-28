#include <iostream>
#include "Vector.h"
#include "AbsorbCoeffData.h"
#include "XRay.h"

#include "CCD.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

typedef boost::mt11213b base_generator_type;


using namespace std;

int main()
{
    Vector InputCCDOrigin(110,-1,50); //origin of CCD
    Vector InputCCDNormal(0,0,1); //direction that CCD points in.
    double InputCCDAngle = 0;

    CCD CCDCamera(InputCCDOrigin, InputCCDNormal, InputCCDAngle,
                  0.05, 0.05,
                  0.0,0.0,
                  0.0,0.0);

    base_generator_type generator(48u);
    boost::uniform_real<> uni_dist(0,1);
    boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);


    AbsorbCoeffData TaMuData( 1.0f, 15.0f, 2000);
    TaMuData.LoadData("FeAbsorbCoeff.txt");

    ofstream FluoResults2( "FluoResultsPostPinhole.txt" );

    ifstream FluoFile("AdvFluoResults.txt");
    string dataline;

    Vector FilterNormal( 112.0, 0, 50.0);
    FilterNormal = FilterNormal.Normalized();

    float FilterThickness = 25000.0; //2.5 micron in A

    float Energy = 1.710;
    float AbsorbCoeff = TaMuData.GetAbsorbCoeffDataPoint(EnergyToWavelength(Energy));

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
                //FluoResults2 << XIntersect  << "\t" << YIntersect << endl;
                FluoResults2 << XPixel  << "\t" << YPixel << endl;
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
                FluoResults2 << XPixel  << "\t" << YPixel << endl;
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

        AbsorbCoeff = TaMuData.GetAbsorbCoeffDataPoint(EnergyToWavelength(Energy));

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
            /*double CCDIntersectX = Source.x + ((50.0-Source.z)/(Direction.z))*Direction.x;
            double CCDIntersectY = Source.y + ((50.0-Source.z)/(Direction.z))*Direction.y;
            DiffractResults2 << CCDIntersectX  << "\t" << CCDIntersectY << "\t" << Energy << endl;*/
            if(CCDCamera.GetPixelRayCCDIntersect(Source,Direction,
                                                 RayLength,IntersectPoint,
                                                 XIntersect,YIntersect,
                                                 XPixel,YPixel) == true)
            {
                //DiffractResults2 << XIntersect  << "\t" << YIntersect << "\t" << Energy << endl;
                DiffractResults2 << XPixel  << "\t" << YPixel << "\t" << Energy << endl;
            }
        }
        else if(uni() < 0.0063*0.5) //fluorescence yield
        {
            /*double CCDIntersectX = Source.x + ((50.0-Source.z)/(Direction.z))*Direction.x;
            double CCDIntersectY = Source.y + ((50.0-Source.z)/(Direction.z))*Direction.y;
            DiffractResults2 << CCDIntersectX  << "\t" << CCDIntersectY << "\t" << 0.7 << endl;*/
            if(CCDCamera.GetPixelRayCCDIntersect(Source,Direction,
                                                 RayLength,IntersectPoint,
                                                 XIntersect,YIntersect,
                                                 XPixel,YPixel) == true)
            {
                //DiffractResults2 << XIntersect  << "\t" << YIntersect << "\t" << 0.7 << endl;
                DiffractResults2 << XPixel  << "\t" << YPixel << "\t" << 0.7 << endl;
            }
        }
    }

    DiffractResults2.close();

}

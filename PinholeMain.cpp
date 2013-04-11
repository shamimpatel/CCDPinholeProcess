#include <iostream>
#include "Vector.h"
#include "AbsorbCoeffData.h"
#include "XRay.h"

#include "CCD.h"
#include "Pinhole.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>



#include "FileReading.h"
#include "FluorescenceData.h"

using namespace std;


#ifndef __generator_type
	typedef boost::mt19937 base_generator_type;
	#define __generator_type
#endif


base_generator_type generator(48u);
boost::uniform_real<> uni_dist(0,1);
boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);

class PinholeFileParser
{
private:
	
	struct _DataLine
	{
		Vector Source;
		Vector Direction;
		double Energy;
	};
	
	
	std::vector< _DataLine >	Data;
	
	std::vector< _DataLine >::iterator Iterator; //position of next dataline to return
	
	std::ifstream DataFile;
	bool bCached;
public:
	
	PinholeFileParser( const char* FileName, bool bCached)
	{
		//if not cached just keep a file handle open
		//otherwise load the entire thing into memory (check the file isn't too big)
		this->bCached = bCached;
		
		DataFile.open(FileName, ios::in);
		if(DataFile.is_open() == false)
		{
			cout << "Error: Failed to open input file: " << FileName << endl;
			exit(1);
		}
		
		
		
		if(bCached)
		{
			std::string dataline;
			while( getline(DataFile, dataline, '\n') )
			{
				_DataLine Line;
				stringstream linestream(dataline);
				linestream >> Line.Source.x >> Line.Source.y >> Line.Direction.x >> Line.Direction.y >> Line.Direction.z >> Line.Energy;
				Data.push_back(Line);
			}
			
			
			Iterator = Data.begin();
			
			DataFile .close();
		}
		
	}
	
	~PinholeFileParser()
	{
		if(!bCached)
		{
			DataFile.close();
		}
	}
	
	
	bool GetLine( Vector &Source, Vector &Direction, double &Energy)
	{
		
		if( bCached )
		{
			if(Iterator != Data.end() )
			{
				Source = Iterator->Source;
				Direction = Iterator->Direction;
				Energy = Iterator->Energy;
				++Iterator; //be ready for next line
				return true;
			}
			else
			{
				return false;
			}
		}
		else
		{
			std::string dataline;
			
			if( getline(DataFile, dataline, '\n') ) //if line available fill values and return true
			{
				stringstream linestream(dataline);
				linestream >> Source.x >> Source.y >> Direction.x >> Direction.y >> Direction.z >> Energy;
				return true;
			}
			else
			{
				return false;
			}
		}
		
	}
};

void ProcessFile( PinholeFileParser *Parser, const char* OutputFileName, AbsorbCoeffData* FilterAbsorbData, CCD* CCDCamera,
				 double FilterThickness, PinholePlane* Pinhole, FluorescenceData* FilterFluoData)
{
	/*ifstream datafile;
	datafile.open(InputFileName, ios::in);
	
	if(datafile.is_open() == false)
	{
		cout << "Error: Failed to open input file: " << InputFileName << endl;
		exit(1);
	}*/
	
	
	//PinholeFileParser Parser( InputFileName, bCached);
	
	
	ofstream ResultsFile;
	ResultsFile.open(OutputFileName);
	
	if(ResultsFile.is_open() == false)
	{
		cout << "Error: Failed to open output file: " << OutputFileName << endl;
		exit(1);
	}
	
	std::string dataline;
	
	Vector FilterNormal = CCDCamera->CCDNormal;
	
	double RayLength;
    Vector IntersectPoint;
    int XPixel,YPixel;
    double XIntersect,YIntersect;
	
	Vector Source, Direction;
	double Energy;
	
	//while(getline(datafile, dataline, '\n'))
	while ( Parser->GetLine( Source, Direction, Energy) )
    {
		/*
        Vector Source(0,0,0);
        Vector Direction(0,0,0);
		double Energy;
		
        stringstream linestream(dataline);
        linestream >> Source.x >> Source.y >> Direction.x >> Direction.y >> Direction.z >> Energy;
		*/
		
		
		
		double AbsorbCoeff = FilterAbsorbData->GetAbsorbCoeffDataPoint(EnergyToWavelength(Energy));
		
        double CosTheta = fabs(Direction.Dot(FilterNormal)); //fabs for anti-parallel.
		
        double PathLength; //effective thickness of crystal given angle of xray to filter
		
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
        
        if(FilterThickness < 0.0)
		{
			ProbTransmit = 2.0; //guarantee transmission
		}
		
        //Don't combine these. It very slightly boosts the fluorescence. This first one *must* fail (ie not transmit) in order to test for fluorescence.
        if(uni() < ProbTransmit)
        {
			if(Pinhole->TestRayPinholeIntersect(Source, Direction, RayLength, IntersectPoint))
			{
				if(CCDCamera->GetPixelRayCCDIntersect(Source,Direction,
													  RayLength,IntersectPoint,
													  XIntersect,YIntersect,
													  XPixel,YPixel) == true)
				{
					
					ResultsFile << XPixel  << "\t" << YPixel << "\t" << Energy << endl;
				}
			}
        }
		else
		{
			
			double FluoEnergy = FilterFluoData->PickFluorescenceEnergy(Energy);
			if(FluoEnergy < 0.0) //energy will be negative if the xray emission does not pass the fluorescence yield test
			{
				continue;
			}
			
			
			float AbsorptionLength = (-1.0f / AbsorbCoeff) * log(uni()); //how far xray travelled into material
			
			int Attempts = 0;
			
			while(AbsorptionLength > PathLength)
			{
				//need to guarantee that the photon is absorbed within the crystal.
				//This is a stupid way to do this... If the absorption probablility is very low there's a reasonable chance of hanging here for a very long time.
				AbsorptionLength = (-1.0f / AbsorbCoeff) * log(uni());
				
				if( Attempts >= 50) //check for the loop hanging here....
				{
					AbsorptionLength = PathLength; //just say it got absorbed at the end.
					cout << "Warning: If this message pops up a lot then the secondary fluorescence count is likely to be too high" << endl;
					break;
				}
				
				Attempts++;				
			}
			
			double RemainingLength = PathLength - AbsorptionLength;
			
			float AbsorbCoeffFluo = FilterAbsorbData->GetAbsorbCoeffDataPoint(EnergyToWavelength(FluoEnergy));
			float ProbFluoTransmit = exp( -1.0f * AbsorbCoeffFluo * RemainingLength  );
			
			if(uni() > ProbFluoTransmit) //if absorbed, continue to next photon
			{
				//this checks if the photon would actually be transmitted (pass back out through the material)
				continue;
			}
			
			if(Pinhole->TestRayPinholeIntersect(Source, Direction, RayLength, IntersectPoint))
			{
				if(CCDCamera->GetPixelRayCCDIntersect(Source,Direction,
													  RayLength,IntersectPoint,
													  XIntersect,YIntersect,
													  XPixel,YPixel) == true)
				{
					
                    ResultsFile << XPixel  << "\t" << YPixel << "\t" << FluoEnergy << endl;
                }
            }
		}
    }
	
    //datafile.close();
	ResultsFile.close();	
}





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
    
    /*
    //optional override so that pinhole parameters can be modified without changing the input script
    //TODO: use switches/flags instead of just reading in order
    if(argc == 3)
    {        
        std::stringstream S1(argv[1]);
        S1 >> PinholeRadius;
        std::stringstream S2(argv[2]);
        S2 >> PinholeDistance;        
    }   
    */
    
    PinholeOrigin = PinholeNormal*PinholeDistance;
    PinholeOrigin.Print();
        
    PinholePlane Pinhole( PinholeOrigin, PinholeNormal, PinholeRadius);



    
    double MinE = 3.0;//4.23;
    double MaxE = 9.0;//4.28;
        
    DoubleFromMap("MinEnergy", InputData, MinE);
    DoubleFromMap("MaxEnergy", InputData, MaxE);
    
    double MinWavelength = EnergyToWavelength(MaxE); MinWavelength -= MinWavelength*0.02;
    double MaxWavelength = EnergyToWavelength(MinE);
    
    if( MaxWavelength < EnergyToWavelength(0.4) )
    {
        MaxWavelength = EnergyToWavelength(0.4);
    }
    
    MaxWavelength += MaxWavelength*0.02;
    
	std::string FilterMuDataFilename;
	
	StringFromMap("FilterAbsorbData", InputData, FilterMuDataFilename);
	
    AbsorbCoeffData MuData( MinWavelength, MaxWavelength, 1, 5000, FilterMuDataFilename.c_str());
    
    Vector FilterNormal = PinholeNormal;
    FilterNormal = FilterNormal.Normalized();

    double FilterThickness = 25000.0; //2.5 micron in A
    DoubleFromMap("FilterThickness", InputData, FilterThickness);


	std::string EmissionLineData;
	std::string ShellProbData;
	
	StringFromMap("FilterEmissionLineData", InputData, EmissionLineData);
	StringFromMap("FilterShellProbabilityData", InputData, ShellProbData);
	
	FluorescenceData FilterFluoData( 0.4, MaxE, 5000, EmissionLineData.c_str(), ShellProbData.c_str(), &uni );
	
	PinholeFileParser FluoParser("AdvFluoResults.txt", false );
	PinholeFileParser DiffractParser("AdvDiffractResults.txt", false );
	
	ProcessFile( &FluoParser, "FluoResultsPostPinhole.txt", &MuData, &CCDCamera,
				FilterThickness, &Pinhole, &FilterFluoData);
	
	ProcessFile( &DiffractParser, "DiffractResultsPostPinhole.txt", &MuData, &CCDCamera,
				FilterThickness, &Pinhole, &FilterFluoData);
	
	
   
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




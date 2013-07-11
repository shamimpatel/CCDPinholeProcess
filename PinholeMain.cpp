#include <iostream>
#include "Vector.h"
#include "AbsorbCoeffData.h"
#include "XRay.h"

#include "CCD.h"
#include "Pinhole.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
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
boost::normal_distribution<> normal_dist(0,1);
boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);
boost::variate_generator<base_generator_type&, boost::normal_distribution<> > normal(generator, normal_dist);


class PinholeFileParser
{
private:
	
	
	struct Ray
	{
		Vector Source, Direction;
		int XPixel;
		double Energy;
		double BlurEnergy;
		bool PassPinhole;
		int PlaneId;
		bool bIsNoise;
	};
	
	
	std::vector< Ray > XRays;
	AbsorbCoeffDataEnergy* FilterAbsorbData;
public:
	
	PinholeFileParser( const char* FileName, CCD* CCDCamera, AbsorbCoeffDataEnergy* FilterAbsorbData, double FilterThickness)
	{
		ifstream DataFile;
		DataFile.open(FileName, ios::in);
		if(DataFile.is_open() == false)
		{
			cout << "Error: Failed to open input file: " << FileName << endl;
			exit(1);
		}
		
		std::string dataline;
		while( getline(DataFile, dataline, '\n') )
		{
			Vector Source, Direction;
			double Energy;
			stringstream linestream(dataline);
			linestream >> Source.x >> Source.y >> Direction.x >> Direction.y >> Direction.z >> Energy;

			
			double RayLength, XIntersect, YIntersect;
			int XPixel, YPixel;
			Vector IntersectPoint;
			
			/*double AbsorbCoeff = FilterAbsorbData->GetAbsorbCoeffDataPoint(EnergyToWavelength(Energy));
			
			double CosTheta = fabs(Direction.Dot(CCDCamera->CCDNormal)); //fabs for anti-parallel.
			
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
			
			
			if( uni() > ProbTransmit )
			{
				continue;
			}*/ 			
			
			
			if(CCDCamera->GetPixelRayCCDIntersect(Source,Direction,
												  RayLength,IntersectPoint,
												  XIntersect,YIntersect,
												  XPixel,YPixel) == true)
			{
				if( XPixel >= 0 && XPixel < CCDCamera->nXPixels &&
					YPixel >= 0 && YPixel < CCDCamera->nYPixels)
				{
					Ray XRay;
					
					XRay.Source = Source;
					XRay.Direction = Direction;
					XRay.XPixel = XPixel;
					XRay.Energy = Energy;
					XRay.bIsNoise = false;
					XRays.push_back(XRay);
				}
			}
			
			
		}
				
		
		DataFile.close();
		
		
		
		//std::vector< Ray > NoiseXRays;
		
		//double NumSignalXRays = 60;
		//double Noise = 0.05;
		int NumNoiseXRays = 800;//= (Noise*NumSignalXRays)/(1.0-Noise);
		
		for( int i = 0; i < NumNoiseXRays; i++)
		{
			Ray Ray;
			
			/*Vector Source, Direction;
			 int XPixel;
			 double Energy;
			 double BlurEnergy;
			 bool PassPinhole;
			 int PlaneId;*/
			
			
			Ray.Energy = 5.564 + (uni()*(6.0565-5.564));
			
			/*const static double w = 0.00365; //average energy to produce eh pair in silicon
			const static double r = 3; //readout noise
			const static double F = 0.117; //Fano factor
			double StDev = w*sqrt( r*r + (F*Ray.Energy)/w );
			
			Ray.BlurEnergy = normal()*StDev + Ray.Energy;*/
			
			Ray.XPixel = uni()*2047.0;
			Ray.PassPinhole = true;
			Ray.PlaneId = 0;
			Ray.bIsNoise = true;
			Ray.PassPinhole = true;
			//Ray.yPixel = uni()*2047.0;
			
			XRays.push_back(Ray);
		}
		
	}

	struct BraggLimits
	{
		std::string PlaneName;
		double MinE;
		double MaxE;
	};
	
	std::vector< std::tuple<double,double,int> >	CalculateCorrelations( PinholePlane* Pinhole, std::vector<BraggLimits> Limits)
	{
	
		std::vector< unsigned long int > SumX, NumXRays;
		std::vector< double > SumY, SumXDiffSquare, SumYDiffSquare, SumYBlurDiffSquare, SumXYDiff, SumXYBlurDiff, SumYBlur;
		std::vector< std::tuple<double,double,int> > Correlations;
		
		
		std::vector< double > XMean, YMean, YBlurMean;
		
		for( int i = 0; i < Limits.size(); i++)
		{
			NumXRays.push_back(0);
			SumX.push_back(0); 
			SumY.push_back(0.0); SumXDiffSquare.push_back(0.0); SumYDiffSquare.push_back(0.0);
			SumXYDiff.push_back(0.0);
			
			//Correlations.push_back( std::make_tuple< double, double >(0.0,0.0,0) );
			
			XMean.push_back(0);
			YMean.push_back(0.0);
			
			
			SumYBlur.push_back(0.0);
			SumYBlurDiffSquare.push_back(0.0);
			SumXYBlurDiff.push_back(0.0);
			YBlurMean.push_back(0.0);
		}
		
		
		
		double RayLength;
		Vector IntersectPoint;
		
		//calculate correlation using a two pass algorithm as the single pass version can sometimes be unstable due to a very high number of XRays. 		
		for( auto Ray = XRays.begin(); Ray != XRays.end(); ++Ray)
		{
			if( Ray->bIsNoise || Pinhole->TestRayPinholeIntersect(Ray->Source, Ray->Direction, RayLength, IntersectPoint))
			{
				Ray->PassPinhole = true;
				
				int PlaneId = 0;
				for( PlaneId = 0; PlaneId < (Limits.size()); PlaneId++)
				{
					if( Ray->Energy >= Limits[PlaneId].MinE && Ray->Energy <= Limits[PlaneId].MaxE)
					{
						break;
					}
				}				
				
				Ray->PlaneId = PlaneId;
				
				
				
				int X = Ray->XPixel;
				double Y = Ray->Energy;
				
				const static double w = 0.00365; //average energy to produce eh pair in silicon
				const static double r = 3; //readout noise
				const static double F = 0.117; //Fano factor
				double StDev = w*sqrt( r*r + (F*Y)/w );
				
				Ray->BlurEnergy = normal()*StDev + Y;
				SumYBlur[PlaneId] += Ray->BlurEnergy;
				
				SumX[PlaneId] += X;
				SumY[PlaneId] += Y;
				NumXRays[PlaneId] += 1;
			}
			else
			{
				Ray->PassPinhole = false;
			}
		}
				
		
		
		for( int i = 0; i < Limits.size(); i++)
		{
			XMean[i] = double(SumX[i])/double(NumXRays[i]);
			YMean[i] = double(SumY[i])/double(NumXRays[i]);
			YBlurMean[i] = double(SumYBlur[i])/double(NumXRays[i]);
		}
		
		for( auto Ray = XRays.begin(); Ray != XRays.end(); ++Ray)
		{
			if( Ray->PassPinhole == false)
			{
				continue;
			}
			
			int PlaneId = Ray->PlaneId;
			double XDiff = double(Ray->XPixel) - XMean[PlaneId];
			double YDiff = double(Ray->Energy) - YMean[PlaneId];
			double YBlurDiff = double(Ray->BlurEnergy) - YBlurMean[PlaneId];
			
			SumXYDiff[PlaneId] += XDiff*YDiff;
			SumXYBlurDiff[PlaneId] += XDiff*YBlurDiff;

			
			SumXDiffSquare[PlaneId] += XDiff*XDiff;
			SumYDiffSquare[PlaneId] += YDiff*YDiff;
			SumYBlurDiffSquare[PlaneId] += YBlurDiff*YBlurDiff;
		}
		
		for( int i = 0; i < Limits.size(); i++)
		{			
			double r = SumXYDiff[i]/sqrt( SumXDiffSquare[i] * SumYDiffSquare[i] );
			double r2 = SumXYBlurDiff[i]/sqrt( SumXDiffSquare[i] * SumYBlurDiffSquare[i] );
			
			//Correlations[i][0] = r;
			//Correlations[i][1] = r2;
			
			
			/*if( NumXRays[i] < 10)
			{
				Correlations[i].first = 0.0;
				Correlations[i].second = 0.0;
			}*/
			
			Correlations.push_back(std::make_tuple(r,r2,NumXRays[i]));			
		}
			
			
		
		return Correlations;
	}

};

void ProcessFileNormal( const char* InputFileName, const char* OutputFileName, AbsorbCoeffDataEnergy* FilterAbsorbData, CCD* CCDCamera,
				 double FilterThickness, PinholePlane* Pinhole, FluorescenceData* FilterFluoData)
{	
	ifstream DataFile;
	DataFile.open(InputFileName, ios::in);
	if(DataFile.is_open() == false)
	{
		cout << "Error: Failed to open input file: " << InputFileName << endl;
		exit(1);
	}
	
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
	
	while ( getline(DataFile, dataline, '\n') )
    {		
		stringstream linestream(dataline);
		linestream >> Source.x >> Source.y >> Direction.x >> Direction.y >> Direction.z >> Energy;

		double AbsorbCoeff = FilterAbsorbData->GetAbsorbCoeffDataPointEnergy(Energy);
	
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
					
					ResultsFile << XPixel  << "\t" << YPixel << "\t" << Energy << "\t" << Source.x << endl;
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
						
			float AbsorbCoeffFluo = FilterAbsorbData->GetAbsorbCoeffDataPointEnergy(FluoEnergy);
			
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
                    ResultsFile << XPixel  << "\t" << YPixel << "\t" << FluoEnergy << "\t" << Source.x << endl;
                }
            }
		}
    }

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
    
    /*double MinWavelength = EnergyToWavelength(MaxE); MinWavelength -= MinWavelength*0.02;
    double MaxWavelength = EnergyToWavelength(MinE);
    
    if( MaxWavelength < EnergyToWavelength(0.1) )
    {
        MaxWavelength = EnergyToWavelength(0.1);
    }
    
    MaxWavelength += MaxWavelength*0.02;*/
    
	std::string FilterMuDataFilename;
	
	StringFromMap("FilterAbsorbData_Energy", InputData, FilterMuDataFilename);
	
    AbsorbCoeffDataEnergy MuData( 0.1, MaxE, 10000, FilterMuDataFilename.c_str());
    
    Vector FilterNormal = PinholeNormal;
    FilterNormal = FilterNormal.Normalized();

    double FilterThickness = 25000.0; //2.5 micron in A
    DoubleFromMap("FilterThickness", InputData, FilterThickness);


	std::string EmissionLineData;
	std::string ShellProbData;
	
	StringFromMap("FilterEmissionLineData", InputData, EmissionLineData);
	StringFromMap("FilterShellProbabilityData", InputData, ShellProbData);
	
	FluorescenceData FilterFluoData( 0.2, MaxE, 5000, EmissionLineData.c_str(), ShellProbData.c_str(), &uni );
	
	
	bool bIterations = false;
	
	std::map<std::string,std::string> IterationData;
	
	if (argc == 2)
	{
		std::string IterationFileName( argv[1] );
		
		ifstream IterationDataFile(IterationFileName.c_str());
		if(IterationDataFile.is_open() == false)
		{
			cout << "Error: Failed to open iteration file: " << IterationFileName << endl;
			exit(1);
		}
		
		AddToMapFromFile(IterationDataFile, IterationData);
		
		IterationDataFile.close();
		
		bIterations = true;
	}
	
	
	if(bIterations)
	{		
		//PinholeFileParser FluoParser("AdvFluoResults.txt", &CCDCamera, &MuData, FilterThickness );
		PinholeFileParser DiffractParser("AdvDiffractResults.txt", &CCDCamera, &MuData, FilterThickness );
		
		ifstream BraggLimitsFile( "DiffractionPeakLimits_112.txt", ios::in);
		
		if(BraggLimitsFile.is_open() == false)
		{
			cout << "Error: Failed to open Bragg limits file" << endl;
			exit(1);
		}
		
		
		std::string dataline;	
		std::vector< PinholeFileParser::BraggLimits > Limits;
		
		cout << "Bragg energy ranges for each plane" << endl;
		
		while( getline(BraggLimitsFile,dataline) )
		{
			stringstream stream(dataline);
			
			PinholeFileParser::BraggLimits Limit;
			
			stream >> Limit.PlaneName >> Limit.MinE >> Limit.MaxE;				
			Limits.push_back(Limit);
			
			
			cout << Limit.PlaneName << ": " << Limit.MinE << " --> " << Limit.MaxE << endl;
			
		}
		
		double MinNoiseE, MaxNoiseE;		
		DoubleFromMap( "MinNoiseE", IterationData, MinNoiseE);
		DoubleFromMap( "MaxNoiseE", IterationData, MaxNoiseE);

		std::vector< PinholeFileParser::BraggLimits > NoiseLimits;
		PinholeFileParser::BraggLimits NoiseLimit;
		
		
		
		NoiseLimit.PlaneName = "Noise";
		NoiseLimit.MinE = MinNoiseE;
		NoiseLimit.MaxE = MaxNoiseE;
		NoiseLimits.push_back(NoiseLimit);
		
		
		double MinRadius = 0.0;
		double MaxRadius = 2.0;
		double RadiusStep = 0.1;
		DoubleFromMap("MinRadius", IterationData, MinRadius);
		DoubleFromMap("MaxRadius", IterationData, MaxRadius);
		DoubleFromMap("RadiusStep", IterationData, RadiusStep);
				
		
		double MinDistance = 0.0;
		double MaxDistance = 2.0;
		double DistanceStep = 0.1;
		DoubleFromMap("MinDistance", IterationData, MinDistance);
		DoubleFromMap("MaxDistance", IterationData, MaxDistance);
		DoubleFromMap("DistanceStep", IterationData, DistanceStep);
		
		ofstream IterationOut( "IterationResults.txt" );
		
		
		IterationOut << "Radius\tDistance\t";
		
		for( int PlaneId = 0; PlaneId < int(Limits.size()); PlaneId++ )
		{
			IterationOut << Limits[PlaneId].PlaneName << "\t\t\t";
		}
		
		IterationOut << "Noise";
		
		IterationOut << endl;
		
		double ProgressCounter = 0.0;
		
		for( double Radius = MinRadius; Radius <= MaxRadius; Radius += RadiusStep)
		{
			for( double Distance = MinDistance; Distance <= MaxDistance; Distance += DistanceStep)
			{
				
				PinholeOrigin = PinholeNormal*Distance;
				
				PinholePlane VariablePinhole( PinholeOrigin, PinholeNormal, Radius);
				
				
				IterationOut << Radius << "\t" << Distance << "\t";
				
				std::vector< std::tuple<double,double,int> > Correlations = DiffractParser.CalculateCorrelations(&VariablePinhole, Limits );
				//std::vector< std::pair<double,double> > NoiseCorrelations = FluoParser.CalculateCorrelations(&VariablePinhole, NoiseLimits );
				
				for( int PlaneId = 0; PlaneId < int(Limits.size()); PlaneId++ )
				{
					IterationOut << std::get<0>(Correlations[PlaneId]) << "\t" <<
									std::get<1>(Correlations[PlaneId]) << "\t" <<
									std::get<2>(Correlations[PlaneId]) << "\t";
				}
				
				
				//IterationOut << NoiseCorrelations[0].first << "\t" << NoiseCorrelations[0].second;
				
				IterationOut << endl;
			}
			
			IterationOut << endl; //for pm3d
			
			double Progress = (Radius-MinRadius)/(MaxRadius - MinRadius);
			
			if( Progress >= ProgressCounter)
			{
				printf( "Progress: %.0f%%\n",Progress*100.0);
				ProgressCounter += 0.01;
			}
			
		}
		
		IterationOut.close();

	}
	else
	{
		ProcessFileNormal( "AdvFluoResults.txt", "FluoResultsPostPinhole.txt", &MuData, &CCDCamera,
						  FilterThickness, &Pinhole, &FilterFluoData);
		
		ProcessFileNormal( "AdvDiffractResults.txt", "DiffractResultsPostPinhole.txt", &MuData, &CCDCamera,
						  FilterThickness, &Pinhole, &FilterFluoData);
	}
	

	
	
	//generate some numbers for gnuplot to draw the edges of the CCD
   
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

//single pass. Unstable to large values of X or Y
double CalcPearsonCorrelation( std::vector<int> x, std::vector<double> y)
{
	if( x.size() != y.size() )
	{
		cout << "Error: Arrays not equal size" << endl;
		exit(1);
	}
	
	double SumXY, SumY, SumYSquare;
	int SumX,SumXSquare;
	
	SumXY = SumY = SumX = SumXSquare = SumYSquare = 0;
	
	for( int i = 0; i < int(x.size()); i++)
	{
		int X = x[i];
		double Y = y[i];
		SumX += X;
		SumXSquare += X*X;
		
		SumY += Y;
		SumYSquare += Y*Y;
		
		SumXY += double(X)*Y;
	}
	
	double n = double(x.size());
	double r = (n*SumXY - SumX*SumY)/ ( sqrt(n*double(SumXSquare) - double(SumX)*double(SumX))*sqrt(n*SumYSquare - SumY*SumY ) );
	
	return r;
}



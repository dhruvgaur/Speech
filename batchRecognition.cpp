// HMM.cpp : Defines the entry point for the console application.
//Commenting all the cout
#include "stdafx.h"
#include<iostream>
#include <sstream>
#include<conio.h>
#include<fstream>
#include<string>
#include<vector>
#include<stdlib.h>
#include<cstdlib>
#include<math.h>
#include<limits.h>
#include<cstdlib>
#include <iomanip>

//Declare variable for number of states : N
int N;

//Declare variable for number of different observations possible : M
int M;

//Declare variable for the size of observation sequence : T
int T=0;

//Declare inital a,b and pi 
double** a;
double** b;
double* pi;

//Declare new a,b and pi to be used afterwards
double** aNew;
double** bNew;
double* piNew;

//Declare observation sequence array to store observation sequence
int *observation;

//Declare variables for Forward Procedure
double** alpha;
double** alphaNew;

//Declare variables for Backward Procedure
double** beta;

//Declare variables for B-W Procedure
double*** Xi;
double** gamma;
double* denominator;

//Declare variables for Viterbi Prodecure
double** delta;
double** psi;

//Declaring miscellaneous variables
long double sum=0,sumNew=0,maxV=0,maxTempV=0,minForward=0;
bool flag=true;int indexForward=0;


int *checkB;

std::ofstream logfile("FullLog.txt",std::ios_base::app);
std::ofstream output("RequiredOutput.txt",std::ios_base::app);

void record(){
	std::system("Recording_Module.exe 3 input.wav input.txt");
}

//Speech Detection and Cepstral coefficient generation

std::vector<signed long long> DCShift(std::vector<signed long long> input , long double *max){
	signed long long sum=0;
	for(int i=0;i<input.size();i++){
		sum += input[i];
		if(input[i]>*max) *max=input[i];		
		
	}
	
	sum=(sum/input.size());
	
	for(int i=0;i<input.size();i++)
		input[i] = input[i]-sum;

	return input;
}

std::vector<signed long long> Normalization(std::vector<signed long long> input, long double factor){

	long double temp;
	for(int i=0;i<input.size();i++){
		
		temp=((long double)input[i])*factor;
		input[i] = ceil(temp);

	}
		
	return input;
}

std::vector<signed long long> hamming(std::vector<signed long long> input){

	long double temp=0;

	for(int i=0;i<input.size();i++){
		
		temp=0.54-0.46*cos((2*3.14*i)/319);
		
		temp = temp * float(input[i]);
		input[i] = (long long) floor(temp + 0.5);
		

	}

	return input;
	
}

std::vector<signed long long> autocorrelationCoeff(std::vector<signed long long> input){


	signed long long rvalue=0;
	std::vector <signed long long> autoCoeff;
	for(int i=0;i<13;i++){
		rvalue=0;
		for(int j=0;j<320-i;j++)
			rvalue +=(input[j]*input[j+i]);
					
		autoCoeff.push_back(rvalue);
	}

	return autoCoeff;
}

std::vector<long double> durbinsCoeff(std::vector<long double> values){

	
	long double E=values.at(0),temp;	
	
	const int p=values.size()-1;
	//Constant defining arrray size of coefficients
	std::vector<long double> result;
	
	//Initializing array of a
	int i,j;

	long double** a = new long double*[p];
	for(int i = 0; i < p; ++i)
    a[i] = new long double[p];

	for(i=0;i<p;i++)
		for(j=0;j<p;j++)
				a[i][j]=0;

	
	a[0][0]=values.at(1)/E;

	E*=(1.0-pow(a[0][0],2));

	
	for(i=1;i<p;i++){


	
		temp=0;
			
		for( j=0;j<i;j++)
			temp+=a[j][i-1]*values.at(i-j);
		
		a[i][i]=(values.at(i+1)-temp)/E;

		
		for(j=i-1;j>=0;j--)
		{a[j][i]=a[j][i-1]-a[i][i]*a[i-j-1][i-1];

		}
			
		E*=(1-pow(a[i][i],2));

	
	}
	//std::cout<<std::endl;
	/*
	for(i=0;i<p;i++)
		{	for	(j=0;j<p;j++)
				std::cout<<a[i][j]<<"::";
			std::cout<<std::endl;		
	}
	*/
	for(i=0;i<p;i++)
		result.push_back(a[i][p-1]);
	
	return result;
}//End of durbinsCoeff

std::vector<long double> cepstral(std::vector<long double> valuesTemp){

	
	std::vector<long double> cepstralCoeffTemp;
	long double temp=0.0;
	cepstralCoeffTemp.push_back(valuesTemp[0]);
	for(int i=1;i<valuesTemp.size();i++){
		temp=0.0;

		
		for(int j=0;j<i;j++){
				
				temp += ((long double(j+1)/long double(i+1))*cepstralCoeffTemp[j]*valuesTemp[i-j-1]);

			}

	
		temp += valuesTemp[i];

		cepstralCoeffTemp.push_back(temp);

	}
	/*
	for(int i=0;i<cepstralCoeffTemp.size();i++)	
		std::cout<<std::endl<<cepstralCoeffTemp[i];
*/
	

	return cepstralCoeffTemp;
}

void createCepstral(){

//All constants are based on the specific sample used for the experiment
	
//Reading configuration files and initializing variables
//record(); //Record input from the user

std::ifstream inputConfig("RecordConfig.txt");
std::string line;

std::getline(inputConfig, line);
const int eThreshold =atoi(line.c_str()); 

std::getline(inputConfig, line);
const std::string inputPath= line;

std::getline(inputConfig, line);
const std::string outputPath= line;

std::getline(inputConfig, line);
const std::string cepstralOutput= line;

//Reading the input text file

std::ifstream infile(inputPath);
std::vector <signed long long> inputVector; 

if (infile) {        
    int value;
   	while ( infile >> value ) {
        inputVector.push_back(value);
    }
}
else
	std::cout << "Unable to open file"; 

infile.close();

//Declaration of variables and vectors for main()

int size=inputVector.size();
int i,j=0,k=0,startSegment=0,endSegment=0,ZCRSumAvg=0,check=0;
std::vector <signed long long> segmentEnergyVector;
std::vector <signed long long> segmentZCRVector;
signed long long energy=0, ZCR=0;

std::vector <signed long long> analysisVector;
std::vector <signed long long> autoCoeff;
std::vector<long double> durbinCoeff;
std::vector<long double> cepstralCoeff;

std::vector<std::vector<long double>> Durbin;
std::vector<std::vector<long double>> Cepstral;

int startSample;
int endSample;
std::ofstream outfile(outputPath);
std::ofstream outfileCepstral(cepstralOutput,std::ios_base::app);
std::string word;

//Calculating Energy and ZCR for the segments of size 300 samples in the file. Ignoring first 1000 samples 

for(i=1000;i<size;i++){

	if(j<300){		
		energy+=(inputVector[i]*inputVector[i]);

		if(inputVector[i]<0)
			if(inputVector[i-1]>=0)
				ZCR++;
		if(inputVector[i]>0)
			if(inputVector[i-1]<=0)
				ZCR++;	
		j++;
	}
	else{
		segmentEnergyVector.push_back(energy);		
		segmentZCRVector.push_back(ZCR);
		energy=0;ZCR=0;j=0;

	}
	
}

for(i=0,j=0;i<segmentEnergyVector.size(),j<segmentZCRVector.size();i++,j++)
	std::cout<<i+1<<" : "<<segmentEnergyVector[i]<<"------"<<segmentZCRVector[j]<<std::endl;


//Calculating start segment of the speech, based on average three consequtive segments above energy 
//threshold and below ZCR threshold, which shows the start of the speech

for(i=1;i<segmentEnergyVector.size()-1;i++)
{

	signed long long avgE=(segmentEnergyVector[i-1]+segmentEnergyVector[i]+segmentEnergyVector[i+1])/3;
	
	
	if(avgE>eThreshold)
		{startSegment=i;break;}
}

//Calculating end segment of the speech, based on average three consequtive segments below energy 
//threshold and above ZCR threshold, which shows the end of the speech

for(i=segmentEnergyVector.size()-2;i>1;i--){

	signed long long avgE=(segmentEnergyVector[i-1]+segmentEnergyVector[i]+segmentEnergyVector[i+1])/3;

	if(avgE>eThreshold)
		{endSegment=i;check=1;break;}

}

//Condition to check if speech is till the end of the file or not
if (!check)endSegment=segmentEnergyVector.size()-2;

//Calculating start and end segment
startSample=(startSegment*300)+1000;
endSample=(endSegment*300)+1000;


if (outfile) {   

	//Output the Average Energy, Average ZCR, Start sample point and End sample point
outfile<<"------------Speech start point and end point-----------"<<std::endl;
outfile<<std::endl<<"Start Segment point is : "<<startSegment;
outfile<<std::endl<<"End Segment point is : "<<endSegment;
outfile<<std::endl<<"Start Sample point is : "<<startSample;
outfile<<std::endl<<"End Sample point is : "<<endSample;
}
else
	std::cout << "Unable to open file"; 

outfile.close(); 

//getch();

long double factor=inputVector[startSample];

//Performing DC Shift
inputVector=DCShift(inputVector,&factor);

std::cout<<std::endl<<"Factor="<<factor;
factor=30000/factor;

std::cout<<std::endl<<"Factor="<<factor;

//Perform Normalization
inputVector=Normalization(inputVector,factor);

endSample = startSample+320;

while(endSample <= ((endSegment*320)+1000) ){

	analysisVector.clear();
	autoCoeff.clear();
	durbinCoeff.clear();
	cepstralCoeff.clear();

	//Storing the frame  of 320 samples from speech into a new analysis vector
	for(int i=startSample;i<endSample;i++)
	analysisVector.push_back(inputVector[i]);

	//Applying Hamming Window over the sample
	analysisVector=hamming(analysisVector);
		
	//Calculating auto-correlation Coefficients R's
	autoCoeff=autocorrelationCoeff(analysisVector);

	//Converting signed long long values to long double values

	for(int i=0;i<autoCoeff.size();i++)
		durbinCoeff.push_back(autoCoeff[i]);
	
	//Calculating Durbin's Coefficients
	durbinCoeff=durbinsCoeff(durbinCoeff);

	Durbin.push_back(durbinCoeff);

	//Calculating Cepstral Coefficients
	cepstralCoeff=cepstral(durbinCoeff);

	Cepstral.push_back(cepstralCoeff);

	startSample+=80;
	endSample=startSample+320;

	if (outfileCepstral) {   

	//Output values of Cepstral's Coeff to the output file
	for(int i=0;i<cepstralCoeff.size();i++)
		if(i<cepstralCoeff.size()-1)
			outfileCepstral<<cepstralCoeff[i]<<",";
		else
			outfileCepstral<<cepstralCoeff[i];

	outfileCepstral<<std::endl;

	}
	else
		std::cout << "Unable to open file"; 

}

std::cout<<"Total number of vectors="<<Cepstral.size();


outfileCepstral.close(); 

std::cout<<std::endl<<std::endl<<std::endl<<std::endl<<std::endl<<"Required output is stored in the output file";


}

//HMM part below

void Forward(){

	sum=0;

	/*
	//Printing what has been read
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				std::cout<<a[i][j]<<" : ";
				
			}
			std::cout<<std::endl;						
		}

		std::cout<<std::endl<<"b:"<<std::endl<<std::endl;		

		for(int i=0;i<N;i++){						
			for(int j=0;j<M;j++){
				std::cout<<b[i][j]<<" : ";
				
			}			
			std::cout<<std::endl;					
		}
		std::cout<<std::endl<<"pi:"<<std::endl<<std::endl;			

		for(int i=0;i<N;i++){
			std::cout<<pi[i]<<" : ";	
		}

	*/

	for(int i=0;i<N;i++)alpha[i][0]=pi[i]*b[i][observation[0]];

		for(int j=1;j<T;j++)
			for(int i=0;i<N;i++){
				double temp=0;
				for(int k=0;k<N;k++)
					temp += alpha[k][j-1]*a[k][i];
			
				alpha[i][j]=temp*b[i][observation[j]];
			
			}
		/*		
		std::cout<<std::endl<<"Alpha :"<<std::endl<<std::endl;
		//logfile<<std::endl<<"Alpha :"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){
			for(int j=0;j<T;j++){
				std::cout<<alpha[i][j]<<" : ";
				//logfile<<alpha[i][j]<<" : ";
			}
	std::cout<<std::endl;
		//logfile<<std::endl;
		}
	*/
					
		for(int i=0;i<N;i++)
			sum += alpha[i][T-1];

		std::cout<<std::endl<<"Probability of partial observation sequence given lamda [ P{O/lamda} ]= "<<sum;
		
		//logfile<<std::endl<<"Probability of partial observation sequence given lamda [ P{O/lamda} ]=  "<<sum;
}


//VectorQuantization part creating observation sequence below

double Tokhura(std::vector<long double> test, std::vector<long double> reference){
	
	double distance=0;const double array[]={1,3,7,13,19,22,25,33,42,50,56,61,};
	for(int i=0;i<test.size();i++){
		
		distance+=array[i]*pow((test[i]-reference[i]),2);
	
	}
	//std::cout<<"Distance="<<distance;
	return distance;
}

void VectorQuantization(){

	std::cout<<"Getting in";

	std::vector<std::vector<long double>> cepstralUniverse;
	std::vector<std::vector<long double>> codeBook;
	std::vector<int> observationSequence;
	std::vector<long double> temp;
	long double distortion;
	long double temp3;

	char line[1000];
	char* token;

	//Read codebook and store it in vector of vectors	

	FILE *fpcodebook=fopen("..//CodeBook.txt","r");
	int i;
	if(fpcodebook != NULL){
	
		while(fgets(line, sizeof line, fpcodebook) != NULL)
		{		
				temp.clear();
				token=strtok(line,":");
				while(token != NULL){

					temp.push_back(std::stod(token));

					token=strtok(NULL,":");
				}
				codeBook.push_back(temp);	
		
			
		}
	
	}else std::cout << "Unable to open CodeBook file"; 
	
	//close(fpcodebook);


	//Read cepstral in vector of cepstrals
	//Reading the input file from which Cepstral coefficients are read and stored in a universe called : 'cepstralUniverse'.

	FILE *fp=fopen("Cepstral.csv","r");
	
	if(fp != NULL){
	
		while(fgets(line, sizeof line, fp) != NULL)
		{		
				temp.clear();
				token=strtok(line,",");
				while(token != NULL){

					temp.push_back(std::stod(token));

					token=strtok(NULL,",");
				}
				cepstralUniverse.push_back(temp);	
		
			
		}
	
	}else std::cout << "Unable to open Cepstral file"; 

	//close(fpcodebook);


	//Write observation sequence in comma seperated value

	std::ofstream outfile("Observation.txt");

	for(int i=0;i<cepstralUniverse.size();i++){

		double min=1000;int indexMin=0;double dist;

		//std::cout<<"For i="<<i<<std::endl;

		for(int j=0;j<codeBook.size();j++){				

				dist=Tokhura(cepstralUniverse[i],codeBook[j]);

				//std::cout<<dist<<std::endl;

				if(dist<min){min=dist;indexMin=j;}
				
			}	

		observationSequence.push_back(indexMin);

	}

	//Write the observation sequence (indexMin) in the csv file
	
	for(int i=0;i<observationSequence.size();i++)
		if(i != observationSequence.size()-1)
			outfile<<observationSequence[i]<<",";
		else
			outfile<<observationSequence[i];
		

}


void Initialize(){

	char line[1000];
	char* token;

	//Defining initial lamda

	FILE *fp=fopen("..//ThirdGenerationLamda//Zero.txt","r");

	if(fp != NULL){

		fgets(line, sizeof line, fp);
		token=strtok(line,",");
		N=atoi(token);

		fgets(line, sizeof line, fp);
		token=strtok(line,",");
		M=atoi(token);

		//Declaring inital a,b and pi 
		a = new double*[N];
		b = new double*[N];
		pi =new double[N];

		//Declaring new a,b and pi to be used afterwards
		aNew = new double*[N];
		bNew = new double*[N];
		piNew =new double[N];

		for(int i=0;i<N;i++){
			a[i] = new double[N];
			b[i] = new double[M];

			aNew[i] = new double[N];
			bNew[i] = new double[M];
		}
	}
		fclose(fp);
	

	///////////////////////Definition of observation sequence and reading file////////////

	FILE *fpObservation=fopen("Observation.txt","r");

	if(fpObservation != NULL){
		
		fgets(line, sizeof line, fpObservation);
		token=strtok(line,",");

		while(token != NULL){
				T++;
				token=strtok(NULL,",");
			}	
	
	}else std::cout<<"Unable to open observation.txt";
	
	fclose(fpObservation);

	observation = new int[T];
	T=0;
	FILE *fpOb=fopen("Observation.txt","r");
	if(fpOb != NULL){

		fgets(line, sizeof line, fpOb);
		token=strtok(line,",");

		while(token != NULL){
				observation[T++]=atoi(token);
				token=strtok(NULL,",");
			}
	
	}else std::cout<<"Unable to open observation.txt";

	
	
	fclose(fpOb);


	for(int i=0;i<T;i++){
		std::cout<<observation[i];
		//observation[i] -= 1;
	}


	/////////////////////////Forward Procedure Definition/////////////////////////////////
	alpha = new double*[N];
	alphaNew = new double*[N];

	for(int i=0;i<N;i++){
		alpha[i]= new double[T];
		alphaNew[i]= new double[T];
	}
		
	
	
}

int _tmain(int argc, _TCHAR* argv[])
{	
	//When training 

	//createCepstral();

	//VectorQuantization();

	//Training();

	//When recognizing
	
	createCepstral();

	VectorQuantization();

	Initialize();

	//Readlamda and calclulate forward

	/////////////////////////////////////////For Zero////////////////////////////////////

	std::ifstream input0("..//ThirdGenerationLamda//Zero.txt");if(!input0)std::cout<<"Lamda file for Zero not read";
	char *token1;
	std::string line1;
	std::getline(input0, line1);
	std::getline(input0, line1);

	//Read a

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input0, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			a[i][j++]=check1;		
			token1=strtok(NULL,",");
		}		
	}

	//Read b

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input0, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			b[i][j++]=check1;		
			token1=strtok(NULL,",");
		}
		
	}

	//Read pi
		//std::string line1;
		std::getline(input0, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			pi[j++]=check1;		
			token1=strtok(NULL,",");
		}

		Forward();

		if(sum > minForward){minForward=sum; indexForward=0;}

		std::cout<<std::endl<<"****"<<sum;


/////////////////////////////////////////For One////////////////////////////////////

	std::ifstream input1("..//ThirdGenerationLamda//One.txt");if(!input0)std::cout<<"Lamda file for One not read";

	std::getline(input1, line1);
	std::getline(input1, line1);

	//Read a

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input1, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			a[i][j++]=check1;		
			token1=strtok(NULL,",");
		}		
	}

	//Read b

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input1, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			b[i][j++]=check1;		
			token1=strtok(NULL,",");
		}
		
	}

	//Read pi
		//std::string line1;
		std::getline(input1, line1);
		c1 = strdup(line1.c_str());

		check1=0;
		token1 = strtok(c1, ",");

		j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			pi[j++]=check1;		
			token1=strtok(NULL,",");
		}

		Forward();

		if(sum > minForward){minForward=sum; indexForward=1;}

		std::cout<<std::endl<<"****"<<sum<<"**"<<indexForward;



/////////////////////////////////////////For Two////////////////////////////////////

	std::ifstream input2("..//ThirdGenerationLamda//Two.txt");if(!input0)std::cout<<"Lamda file for Two not read";

	std::getline(input2, line1);
	std::getline(input2, line1);

	//Read a

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input2, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			a[i][j++]=check1;		
			token1=strtok(NULL,",");
		}		
	}

	//Read b

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input2, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			b[i][j++]=check1;		
			token1=strtok(NULL,",");
		}
		
	}

	//Read pi
		//std::string line1;
		std::getline(input2, line1);
		c1 = strdup(line1.c_str());

		check1=0;
		token1 = strtok(c1, ",");

		j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			pi[j++]=check1;		
			token1=strtok(NULL,",");
		}

		Forward();

		if(sum > minForward){minForward=sum; indexForward=2;}

		std::cout<<std::endl<<"****"<<sum<<"**"<<indexForward;


/////////////////////////////////////////For Three////////////////////////////////////

	std::ifstream input3("..//ThirdGenerationLamda//Three.txt");if(!input0)std::cout<<"Lamda file for Three not read";

	std::getline(input3, line1);
	std::getline(input3, line1);

	//Read a

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input3, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			a[i][j++]=check1;		
			token1=strtok(NULL,",");
		}		
	}

	//Read b

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input3, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			b[i][j++]=check1;		
			token1=strtok(NULL,",");
		}
		
	}

	//Read pi
		//std::string line1;
		std::getline(input3, line1);
		c1 = strdup(line1.c_str());

		check1=0;
		token1 = strtok(c1, ",");

		j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			pi[j++]=check1;		
			token1=strtok(NULL,",");
		}

		Forward();

		if(sum > minForward){minForward=sum; indexForward=3;}

		std::cout<<std::endl<<"****"<<sum<<"**"<<indexForward;

/////////////////////////////////////////For Four////////////////////////////////////

	std::ifstream input4("..//ThirdGenerationLamda//Four.txt");if(!input0)std::cout<<"Lamda file for Four not read";

	std::getline(input4, line1);
	std::getline(input4, line1);

	//Read a

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input4, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			a[i][j++]=check1;		
			token1=strtok(NULL,",");
		}		
	}

	//Read b

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input4, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			b[i][j++]=check1;		
			token1=strtok(NULL,",");
		}
		
	}

	//Read pi
		//std::string line1;
		std::getline(input4, line1);
		c1 = strdup(line1.c_str());

		check1=0;
		token1 = strtok(c1, ",");

		j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			pi[j++]=check1;		
			token1=strtok(NULL,",");
		}

		Forward();

		if(sum > minForward){minForward=sum; indexForward=4;}

		std::cout<<std::endl<<"****"<<sum<<"**"<<indexForward;

		/////////////////////////////////////////For Five////////////////////////////////////

	std::ifstream input5("..//ThirdGenerationLamda//Five.txt");if(!input0)std::cout<<"Lamda file for Five not read";

	std::getline(input5, line1);
	std::getline(input5, line1);

	//Read a

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input5, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			a[i][j++]=check1;		
			token1=strtok(NULL,",");
		}		
	}

	//Read b

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input5, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			b[i][j++]=check1;		
			token1=strtok(NULL,",");
		}
		
	}

	//Read pi
		//std::string line1;
		std::getline(input5, line1);
		c1 = strdup(line1.c_str());

		check1=0;
		token1 = strtok(c1, ",");

		j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			pi[j++]=check1;		
			token1=strtok(NULL,",");
		}

		Forward();

		if(sum > minForward){minForward=sum; indexForward=5;}

		std::cout<<std::endl<<"****"<<sum<<"**"<<indexForward;

		/////////////////////////////////////////For Six////////////////////////////////////

	std::ifstream input6("..//ThirdGenerationLamda//Six.txt");if(!input0)std::cout<<"Lamda file for Six not read";

	std::getline(input6, line1);
	std::getline(input6, line1);

	//Read a

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input6, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			a[i][j++]=check1;		
			token1=strtok(NULL,",");
		}		
	}

	//Read b

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input6, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			b[i][j++]=check1;		
			token1=strtok(NULL,",");
		}
		
	}

	//Read pi
		//std::string line1;
		std::getline(input6, line1);
		c1 = strdup(line1.c_str());

		check1=0;
		token1 = strtok(c1, ",");

		j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			pi[j++]=check1;		
			token1=strtok(NULL,",");
		}

		Forward();

		if(sum > minForward){minForward=sum; indexForward=6;}

		std::cout<<std::endl<<"****"<<sum<<"**"<<indexForward;

		/////////////////////////////////////////For Seven////////////////////////////////////

	std::ifstream input7("..//ThirdGenerationLamda//Seven.txt");if(!input0)std::cout<<"Lamda file for Seven not read";

	std::getline(input7, line1);
	std::getline(input7, line1);

	//Read a

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input7, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			a[i][j++]=check1;		
			token1=strtok(NULL,",");
		}		
	}

	//Read b

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input7, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			b[i][j++]=check1;		
			token1=strtok(NULL,",");
		}
		
	}

	//Read pi
		//std::string line1;
		std::getline(input7, line1);
		c1 = strdup(line1.c_str());

		check1=0;
		token1 = strtok(c1, ",");

		j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			pi[j++]=check1;		
			token1=strtok(NULL,",");
		}

		Forward();

		if(sum > minForward){minForward=sum; indexForward=7;}

		std::cout<<std::endl<<"****"<<sum<<"**"<<indexForward;

		/////////////////////////////////////////For Eight////////////////////////////////////

	std::ifstream input8("..//ThirdGenerationLamda//Eight.txt");if(!input0)std::cout<<"Lamda file for Eight not read";

	std::getline(input8, line1);
	std::getline(input8, line1);

	//Read a

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input8, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			a[i][j++]=check1;		
			token1=strtok(NULL,",");
		}		
	}

	//Read b

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input8, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			b[i][j++]=check1;		
			token1=strtok(NULL,",");
		}
		
	}

	//Read pi
		//std::string line1;
		std::getline(input8, line1);
		c1 = strdup(line1.c_str());

		check1=0;
		token1 = strtok(c1, ",");

		j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			pi[j++]=check1;		
			token1=strtok(NULL,",");
		}

		Forward();

		if(sum > minForward){minForward=sum; indexForward=8;}

		std::cout<<std::endl<<"****"<<sum<<"**"<<indexForward;


		/////////////////////////////////////////For Nine////////////////////////////////////

	std::ifstream input9("..//ThirdGenerationLamda//Nine.txt");if(!input0)std::cout<<"Lamda file for Nine not read";

	std::getline(input9, line1);
	std::getline(input9, line1);

	//Read a

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input9, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			a[i][j++]=check1;		
			token1=strtok(NULL,",");
		}		
	}

	//Read b

	for(int i=0;i<N;i++){

		//std::string line1;
		std::getline(input9, line1);
		char *c1 = strdup(line1.c_str());

		double check1=0;
		token1 = strtok(c1, ",");

		int j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			b[i][j++]=check1;		
			token1=strtok(NULL,",");
		}
		
	}

	//Read pi
		//std::string line1;
		std::getline(input9, line1);
		c1 = strdup(line1.c_str());

		check1=0;
		token1 = strtok(c1, ",");

		j=0;
		while(token1 != NULL){

			std::string someString1(token1);
			std::istringstream os1 (someString1);
			os1 >> check1;
			pi[j++]=check1;		
			token1=strtok(NULL,",");
		}

		Forward();

		if(sum > minForward){minForward=sum; indexForward=9;}

		std::cout<<std::endl<<"****"<<sum<<"**"<<indexForward;




		std::cout<<std::endl<<std::endl<<"----------------------------------------------"<<std::endl<<std::endl;

		std::cout<<"Spoken word is :" << indexForward;
		

	getch();
	return 0;
}


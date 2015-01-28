// HMM.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include<iostream>
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
long double sum=0,sumNew=0,maxV=0,maxTempV=0;	
bool flag=true;

int *checkB;

std::ofstream logfile("Input//HMMTraining//Output//FullLog.txt",std::ios_base::app);
std::ofstream output("Input//HMMTraining//Output//RequiredOutput.txt",std::ios_base::app);

void record(){
	std::system("Recording_Module.exe 3 Input\\SpeechDetection\\Input\\input_no_3.wav Input\\SpeechDetection\\Input\\input_no_3.txt");
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
	std::cout<<std::endl;

	for(i=0;i<p;i++)
		{	for	(j=0;j<p;j++)
				std::cout<<a[i][j]<<"::";
			std::cout<<std::endl;		
	}

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

	for(int i=0;i<cepstralCoeffTemp.size();i++)	
		std::cout<<std::endl<<cepstralCoeffTemp[i];

	
	return cepstralCoeffTemp;
}

void createCepstral(){

//All constants are based on the specific sample used for the experiment
	
//Reading configuration files and initializing variables
//record(); //Record input from the user

std::ifstream inputConfig("Input//SpeechDetection//RecordConfig.txt");
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

getch();

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

void Initialize(){

	logfile<<"-----------------Log file for HMM-----------------"<<std::endl;
	logfile<<"Created by : Dhruv Gaur"<<std::endl<<"Roll # 144101068"<<std::endl;
	logfile<<"Description : It is a dump of all the structures which are used during HMM calculations while solving all three problems."<<std::endl<<std::endl<<std::endl<<std::endl;

	output<<"-----------------Output file for HMM ( as required for the assignment) -----------------"<<std::endl;
	output<<"Created by : Dhruv Gaur"<<std::endl<<"Roll # 144101068"<<std::endl;
	output<<"Contains : "<<std::endl;
	output<<"(1) The re-estimated model (matrices a, b and pi)"<<std::endl;
	output<<"(2) The probability of an utterance being from the model."<<std::endl;
	output<<"(3) The estimated state sequence."<<std::endl<<std::endl<<std::endl<<std::endl;


	///////////////////////Definition of initial lamda and reading file//////////////////
	char line[1000];
	char* token;
	
	FILE *fp=fopen("Input//HMMTraining//Input//Initialization.txt","r");

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

		//Reading a from the file
		for(int i=0;i<N;i++){
			int j=0;
			fgets(line, sizeof line, fp);
			token=strtok(line,",");			

			while(token != NULL){
					a[i][j++]=std::stod(token);
					token=strtok(NULL,",");
				}	
		}

		//Reading b from the file
		for(int i=0;i<N;i++){
			int j=0;
			fgets(line, sizeof line, fp);
			token=strtok(line,",");			

			while(token != NULL){
					b[i][j++]=std::stod(token);
					token=strtok(NULL,",");
				}	
		}
		
		//Reading pi from the file

		fgets(line, sizeof line, fp);
		token=strtok(line,",");int j=0;
		while(token != NULL){
			pi[j++]=std::stod(token);
			token=strtok(NULL,",");
		}

	}else std::cout << "Unable to open Initialization.txt file"; 

	fclose(fp);
	
	///////////////////////Definition of observation sequence and reading file////////////

	FILE *fpObservation=fopen("Input//ObservationSequence//Observation.txt","r");

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
	FILE *fpOb=fopen("Input//ObservationSequence//Observation.txt","r");
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
		observation[i] -= 1;
	}

	
	
	//Definition of array to check if there are any 0's in B
	checkB= new int[M];

	/////////////////////////Forward Procedure Definition/////////////////////////////////
	alpha = new double*[N];
	alphaNew = new double*[N];

	for(int i=0;i<N;i++){
		alpha[i]= new double[T];
		alphaNew[i]= new double[T];
	}
		

	/////////////////////////Backward Procedure Definition/////////////////////////////////
	beta = new double*[N];
	for(int i=0;i<N;i++)
		beta[i]= new double[T];	

	/////////////////////////////////B-W Definition////////////////////////////////////////
	//Defining Xi: index of T goes only upto T-1
	Xi = new double**[N];
	for(int i=0;i<N;i++)
		Xi[i] = new double*[N];

	for(int i=0;i<N;i++)
		for(int j=0;j<N;j++)
			Xi[i][j]=new double[T-1];

	//Declaring gamma
	gamma=new double*[N];
	
	for(int i=0;i<N;i++)
		gamma[i]=new double[T];

	//Declaring denominator
	denominator = new double[T];

	/////////////////////////Viterbi Procedure Definition//////////////////////////////////
	//Declaring delta and psi of Viterbi algorithm
	delta = new double*[N];
	psi = new double*[N];
	for(int i=0;i<N;i++){
		delta[i]= new double[T];
		psi[i] = new double[T];
	}

}

void Forward(){
	

	for(int i=0;i<N;i++)alpha[i][0]=pi[i]*b[i][observation[0]];

		for(int j=1;j<T;j++)
			for(int i=0;i<N;i++){
				double temp=0;
				for(int k=0;k<N;k++)
					temp += alpha[k][j-1]*a[k][i];
			
				alpha[i][j]=temp*b[i][observation[j]];
			
			}
					
		//std::cout<<std::endl<<"Alpha :"<<std::endl<<std::endl;
		logfile<<std::endl<<"Alpha :"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){
			for(int j=0;j<T;j++){
				//std::cout<<alpha[i][j]<<" : ";
				logfile<<alpha[i][j]<<" : ";
			}
		//std::cout<<std::endl;
		logfile<<std::endl;
		}
	
					
		for(int i=0;i<N;i++)
			sum += alpha[i][T-1];

		//std::cout<<std::endl<<"Probability of partial observation sequence given lamda [ P{O/lamda} ]= "<<sum;
		logfile<<std::endl<<"Probability of partial observation sequence given lamda [ P{O/lamda} ]=  "<<sum;
}

void Backward(){
	//std::cout<<std::endl<<std::endl;
	logfile<<std::endl<<std::endl;
	

		for(int i=0;i<N;i++) beta[i][T-1]=1;

		for(int j=T-2;j>=0;j--)
			for(int i=0;i<N;i++){
				double temp=0;
				for(int k=0;k<N;k++)
					temp += a[i][k]*b[k][observation[j+1]]*beta[k][j+1];

				beta[i][j]=temp;
			}

			//std::cout<<std::endl<<"Beta :"<<std::endl<<std::endl;
			logfile<<std::endl<<"Beta :"<<std::endl<<std::endl;


		for(int i=0;i<N;i++){
			for(int j=0;j<T;j++){

				//std::cout<<beta[i][j]<<":";
				logfile<<beta[i][j]<<":";

			
			}
				
				//std::cout<<std::endl;
				logfile<<std::endl;
		}		
}

void Viterbi(){	
		//Initialization step for delta and psi

		for(int i=0;i<N;i++){
		delta[i][0]=pi[i]*b[i][observation[0]];
		psi[i][0]=0;
		}

		//Recurrsive step for delta and psi in Viterbi algorithm

		for(int j=1;j<T;j++)
			for(int i=0;i<N;i++){
				double max=delta[0][j-1]*a[0][i];int index=0;
				for(int k=1;k<N;k++)
					if(max < (delta[k][j-1] * a[k][i])){
					
						max=delta[k][j-1] * a[k][i];
						index=k;

					}

				delta[i][j]=max*b[i][observation[j]];

				psi[i][j]=index;
			
			}

		//Displaying delta after the recurrsive step

	//	std::cout<<std::endl<<"Delta:"<<std::endl<<std::endl;
		logfile<<std::endl<<"Delta:"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){
			for(int j=0;j<T;j++){
				//std::cout<<delta[i][j]<<" : ";
				logfile<<delta[i][j]<<" : ";			
			}
				
		//std::cout<<std::endl;
		logfile<<std::endl;
		}

		//Displaying psi after the recurrsive step

		//std::cout<<std::endl<<"Psi:"<<std::endl<<std::endl;
		logfile<<std::endl<<"Psi:"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){
			for(int j=0;j<T;j++){
				//std::cout<<psi[i][j]<<" : ";

				logfile<<psi[i][j]<<" : ";
			}
				
		//std::cout<<std::endl;
		logfile<<std::endl;
		}

		//Termination step for Viterbi algorithm

		long double max=delta[0][T-1];

		int index=0;double temp;
		for(int i=1;i<N;i++)
			if(max < delta[i][T-1]){
			
				max=delta[i][T-1];
				index=i;
			}
			

		std::cout<<std::endl<<"Probability of determined state sequence(P*) from first Viterbi= "<<max;
		logfile<<std::endl<<"Probability of determined state sequence(P*)from first Viterbi= "<<max;

		//Determining the state sequence in reverse order and storing them in correct order in array : sequence

		int *sequence = new int[T];
	
		int count=T-1;

		sequence[count--]=index;

		for(int j=T-2;j>=0;j--){

			temp=psi[index][j+1];
			sequence[count--]=temp;
			index=temp;		
		}

		//Displaying array sequence

		//std::cout<<std::endl<<std::endl<<"Displaying the state sequence for lamda using Viterbi (Qt*):"<<std::endl;
		logfile<<std::endl<<std::endl<<"Displaying the state sequence for lamda using Viterbi (Qt*):"<<std::endl;

		for(int i=0;i<T;i++){
			std::cout<<sequence[i]<<" : ";
			logfile<<sequence[i]<<" : ";		
		}
		
		maxV=max;
}

void Baum_Welch(){

	//std::cout<<std::endl<<std::endl;
	logfile<<std::endl<<std::endl;

		//Calculation the denominator beforehand for all the values of T

		for(int k=0;k<T-1;k++){
		
			denominator[k]=0;

			//Book variable : varable here(i:i,j:j,t:k)

			for(int i=0;i<N;i++)
				for(int j=0;j<N;j++)
					denominator[k] += (alpha[i][k]*a[i][j]*b[j][observation[k+1]]*beta[j][k+1]);
	
		}

		//std::cout<<std::endl<<"Printing denominator for each T :"<<std::endl<<std::endl;
		logfile<<std::endl<<"Printing denominator for each T :"<<std::endl<<std::endl;

		for(int i=0;i<T-1;i++){
			//std::cout<<denominator[i]<<" : ";
			logfile<<denominator[i]<<" : ";
		}
		

		//B-M procedure main part: Calculating Xi 

		for(int i=0;i<N;i++)
			for(int j=0;j<N;j++){

				for(int k=0;k<T-1;k++){

					//(Book variable : variable here)(i:i,j:j,t:k)
					double numerator = alpha[i][k]*a[i][j]*beta[j][k+1]*b[j][observation[k+1]];
					if(denominator[k] != 0)
						Xi[i][j][k]=numerator/denominator[k];
					else
						Xi[i][j][k]=0;
							
				}
						
			}
					
		//B-M procedure main part: Calculating gamma 	
						
			for(int i=0;i<N;i++)
				for(int j=0;j<T-1;j++)
				{double gammaTemp=0;

					for(int k=0;k<N;k++)
						gammaTemp += Xi[i][k][j];

					gamma[i][j]=gammaTemp;
							
							
				}					
				

			//std::cout<<std::endl<<"Xi:"<<std::endl;
			logfile<<std::endl<<"Xi:"<<std::endl;

			for(int i=0;i<N;i++){

				//std::cout<<std::endl<<"layer #"<<i<<std::endl<<std::endl;
				logfile<<std::endl<<"layer #"<<i<<std::endl<<std::endl;

				for(int j=0;j<N;j++){

					for(int k=0;k<T-1;k++){
						//std::cout<<Xi[i][j][k]<<" : "; 
						logfile<<Xi[i][j][k]<<" : ";
					
					}
						

					//std::cout<<std::endl;
					logfile<<std::endl;
							
				}

			//std::cout<<std::endl;	
			logfile<<std::endl;
						
			}

			//std::cout<<std::endl<<"Gamma:"<<std::endl<<std::endl;
			logfile<<std::endl<<"Gamma:"<<std::endl<<std::endl;
			
			for(int i=0;i<N;i++){
			
				for(int j=0;j<T-1;j++)	{
					//std::cout<<gamma[i][j]<<" : "; 
					logfile<<gamma[i][j]<<" : ";
				
				}									
							
					
				//std::cout<<std::endl;
				logfile<<std::endl;
			
			}
}

void EstimateNewLamda(){

	output<<"(1) The re-estimated model (matrices a, b and pi):"<<std::endl<<std::endl;

	for(int i=0;i<N;i++)
			piNew[i]=gamma[i][0];

		for(int i=0;i<N;i++)
			for(int j=0;j<N;j++){
				//(Book variable : variable here)(i:i,j:j,t:k)
				double aNumerator=0,aDenominator=0;

				for(int k=0;k<T-1;k++){

					aNumerator += Xi[i][j][k];
					aDenominator += gamma[i][k];
													
				}
					
				aNew[i][j]=aNumerator/aDenominator;				

			}

			for(int i=0;i<N;i++)
				for(int j=0;j<M;j++){
					//(Book variable : variable here)(j:i,k:j,t:k)
					double bNumerator=0,bDenominator=0;
					
					for(int k=0;k<T-1;k++){
					
						if(observation[k]==j)
							bNumerator += gamma[i][k];

						bDenominator += gamma[i][k];
					}

					bNew[i][j]=bNumerator/bDenominator;


				}		

		//Check for 0 values in b and replace them after splitting 
				
		for(int i=0;i<N;i++){

			for(int i=0;i<M;i++) checkB[i]=0;
			double max=bNew[i][0];int index=0;

			for(int j=0;j<M;j++){

				if(bNew[i][j] > max){max=bNew[i][j];index=j;}

				if(bNew[i][j]==0) checkB[j]=1;

			}

			for(int j=0;j<M;j++)
				if(checkB[j]==1){					
					//Split from largest
					bNew[i][j] = bNew[i][index]*(1.000000e-100);
					bNew[i][index] -= bNew[i][j];

				}			
		
		}

			
		//Printing new a,b and pi

		//std::cout<<std::endl<<"aNew:"<<std::endl<<std::endl;
		logfile<<std::endl<<"aNew:"<<std::endl<<std::endl;
		output<<std::endl<<"aNew:"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){

			for(int j=0;j<N;j++){
				//std::cout<<aNew[i][j]<<" : ";
				logfile<<aNew[i][j]<<" : ";
				output<<aNew[i][j]<<" : ";
			
			}
				

			//std::cout<<std::endl;
			logfile<<std::endl;
			output<<std::endl;
					
		}

		//std::cout<<std::endl<<"bNew:"<<std::endl<<std::endl;
		logfile<<std::endl<<"bNew:"<<std::endl<<std::endl;
		output<<std::endl<<"bNew:"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){
						
			for(int j=0;j<M;j++){
				//std::cout<<bNew[i][j]<<" : ";
				logfile<<bNew[i][j]<<" : ";
				output<<bNew[i][j]<<" : ";
			
			}
				

			//std::cout<<std::endl;
			logfile<<std::endl;
			output<<std::endl;
					
		}


		//std::cout<<std::endl<<"piNew:"<<std::endl<<std::endl;	
		logfile<<std::endl<<"piNew:"<<std::endl<<std::endl;
		output<<std::endl<<"piNew:"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){
			//std::cout<<piNew[i]<<" : ";
			logfile<<piNew[i]<<" : ";
			output<<piNew[i]<<" : ";
		
		}
			
}

void SubsequentForward(){
	output<<std::endl<<std::endl<<std::endl<<"(2) The probability of an utterance being from the model :"<<std::endl<<std::endl;

	std::cout<<std::endl<<observation[0];

	for(int i=0;i<N;i++)alphaNew[i][0]=piNew[i]*bNew[i][observation[0]];

		for(int j=1;j<T;j++)
			for(int i=0;i<N;i++){
				double temp=0;
				for(int k=0;k<N;k++)
					temp += alphaNew[k][j-1]*aNew[k][i];
			
				alphaNew[i][j]=temp*bNew[i][observation[j]];
			
			}
					
		//std::cout<<std::endl<<std::endl<<"AlphaNew :"<<std::endl<<std::endl;
		logfile<<std::endl<<std::endl<<"AlphaNew :"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){
			for(int j=0;j<T;j++){
				//std::cout<<alphaNew[i][j]<<" : ";
				logfile<<alphaNew[i][j]<<" : ";
			
			}
				
		//std::cout<<std::endl;
		logfile<<std::endl;
		}
	
					
		for(int i=0;i<N;i++)
			sumNew += alphaNew[i][T-1];

		//std::cout<<std::endl<<"Probability of partial observation sequence given new lamda [ P{O/lamda} ]= "<<sumNew;
		logfile<<std::endl<<"Probability of partial observation sequence given new lamda [ P{O/lamda} ]="<<sumNew;
		output<<std::endl<<"Probability of partial observation sequence given new lamda [ P{O/lamda} ]="<<sumNew;
}

void SubsequentViterbi(){
	
	output<<std::endl<<std::endl<<std::endl<<"(3) The estimated state sequence :"<<std::endl<<std::endl;


	//Initialization step for delta and psi

		for(int i=0;i<N;i++){
		delta[i][0]=piNew[i]*bNew[i][observation[0]];
		psi[i][0]=0;
		}

		//Recurrsive step for delta and psi in Viterbi algorithm

		for(int j=1;j<T;j++)
			for(int i=0;i<N;i++){
				double maxTemp=delta[0][j-1]*aNew[0][i];int indexTemp=0;
				for(int k=1;k<N;k++)
					if(maxTemp < (delta[k][j-1] * aNew[k][i])){
					
						maxTemp=delta[k][j-1] * aNew[k][i];
						indexTemp=k;

					}

				delta[i][j]=maxTemp*bNew[i][observation[j]];

				psi[i][j]=indexTemp;
			
			}

		//Displaying delta after the recurrsive step

		//std::cout<<std::endl<<std::endl<<"Delta:"<<std::endl<<std::endl;
		logfile<<std::endl<<std::endl<<"Delta:"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){
			for(int j=0;j<T;j++){

				//std::cout<<delta[i][j]<<" : ";
				logfile<<delta[i][j]<<" : ";
			
			}
				
		//std::cout<<std::endl;
		logfile<<std::endl;
		}

		//Displaying psi after the recurrsive step

		//std::cout<<std::endl<<"Psi:"<<std::endl<<std::endl;
		logfile<<std::endl<<"Psi:"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){
			for(int j=0;j<T;j++){

				//std::cout<<psi[i][j]<<" : ";
				logfile<<psi[i][j]<<" : ";
			
			}
				
		//std::cout<<std::endl;
		logfile<<std::endl;
		}

		//Termination step for Viterbi algorithm

		long double maxTemp=delta[0][T-1];
		int indexTemp=0;double tempTemp;
		for(int i=1;i<N;i++)
			if(maxTemp < delta[i][T-1]){
			
				maxTemp=delta[i][T-1];
				indexTemp=i;
			}
			

		//std::cout<<std::endl<<"Probability of determined state sequence(P*)= "<<maxTemp;
		logfile<<std::endl<<"Probability of determined state sequence(P*)= "<<maxTemp;
		output<<std::endl<<"Probability of determined state sequence(P*)= "<<maxTemp;


		//Determining the state sequence in reverse order and storing them in correct order in array : sequenceTemp

		int *sequenceTemp = new int[T];
	
		int countTemp=T-1;

		sequenceTemp[countTemp--]=indexTemp;

		for(int j=T-2;j>=0;j--){

			tempTemp=psi[indexTemp][j+1];
			sequenceTemp[countTemp--]=tempTemp;
			indexTemp=tempTemp;		
		}

		//Displaying array sequence

		//std::cout<<std::endl<<std::endl<<"Displaying the state sequence for lamda using Viterbi (Qt*):"<<std::endl;
		logfile<<std::endl<<std::endl<<"Displaying the state sequence for lamda using Viterbi (Qt*):"<<std::endl;
		output<<std::endl<<std::endl<<"Displaying the state sequence for lamda using Viterbi (Qt*):"<<std::endl;

		for(int i=0;i<T;i++){
			//std::cout<<sequenceTemp[i]<<" : ";
			logfile<<sequenceTemp[i]<<" : ";
			output<<sequenceTemp[i]<<" : ";
		
		}

		maxTempV=maxTemp;
			
}

void Training(){

	std::cout<< std::scientific;
	std::cout<<std::setprecision(6);

	logfile<< std::scientific;
	logfile<<std::setprecision(6);

	output<< std::scientific;
	output<<std::setprecision(6);
	
	//Read data from Initialization.txt file and initialize variables
	Initialize();
	
	int count=0;

	//Here goes the iterations to stabilize lamda	

do{
		
		//Forward procedure : alpha calculation
		Forward();

		//Backward procedure : beta calculation
		Backward();	

		//Viterbi procedure : delta and psi calculation
		Viterbi();

		//B-W procedure : Xi and gamma calculation
		Baum_Welch();

		//Estimating new a,b,and pi
		EstimateNewLamda();

		//Again calculating Forward Procedure for new lamda
		SubsequentForward();

		//Again calculating Viterbi on new lamda
		SubsequentViterbi();
		
		std::cout<<"maxTempV="<<maxTempV<<"maxV="<<maxV;

		if((sumNew==sum) || (maxTempV < maxV) )
			flag=false;
		else{
		
			//Write values of following into the file : a,b,pi,alpha,beta,delta,psi,gamma,Xi,observation/lamda,observation/lamdaNew,p* from Viterbi and qT* from Viterbi
					
			for(int i=0;i<N;i++)
				for(int j=0;j<N;j++)
					a[i][j]=aNew[i][j];

			for(int i=0;i<N;i++)
				for(int j=0;j<M;j++)
					b[i][j]=bNew[i][j];

			for(int i=0;i<N;i++)
				pi[i]=piNew[i];				
				
		}
						
	//std::cout<<std::endl<<"Reaching the end of while iteration"<<std::endl;	
	std::cout<<std::endl<<"Count="<<count++;
}while((flag) && (count < 100));

//write values of following into the file : a,b,pi,alpha,beta,delta,psi,gamma,Xi,observation/lamda,observation/lamdaNew,p* from Viterbi and qT* from Viterbi

		

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

	FILE *fpcodebook=fopen("Input//ObservationSequence//CodeBook.txt","r");
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

	FILE *fp=fopen("Input//SpeechDetection//Output//Cepstral.csv","r");
	
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

	std::ofstream outfile("Input\\ObservationSequence\\Observation.txt");

	for(int i=0;i<cepstralUniverse.size();i++){

		double min=1000;int indexMin=0;double dist;

		std::cout<<"For i="<<i<<std::endl;

		for(int j=0;j<codeBook.size();j++){				

				dist=Tokhura(cepstralUniverse[i],codeBook[j]);

				std::cout<<dist<<std::endl;

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

void Recognition(){

//1. Read all the lamda's stored in Recognition Folderwd

//2. Read the observation sequence from observation sequence folder

//3. Calculate forward procedure for observation sequence and all lamdas and recognize the word

}

int _tmain(int argc, _TCHAR* argv[])
{	
	//When training 

	//createCepstral();

	//VectorQuantization();

	Training();

	//When recognizing

	/*createCepstral();

	VectorQuantization();

	Recognition();*/
	
	/*
	std::cout<< std::scientific;
	std::cout<<std::setprecision(6);

	logfile<< std::scientific;
	logfile<<std::setprecision(6);

	output<< std::scientific;
	output<<std::setprecision(6);
	
	//Read data from Initialization.txt file and initialize variables
	Initialize();
	
	int count=0;

	//Here goes the iterations to stabilize lamda	

do{
		
		//Forward procedure : alpha calculation
		Forward();

		//Backward procedure : beta calculation
		Backward();	

		//Viterbi procedure : delta and psi calculation
		Viterbi();

		//B-W procedure : Xi and gamma calculation
		Baum_Welch();

		//Estimating new a,b,and pi
		EstimateNewLamda();

		//Again calculating Forward Procedure for new lamda
		SubsequentForward();

		//Again calculating Viterbi on new lamda
		SubsequentViterbi();
		
							
		if((sumNew==sum) || (maxTempV > maxV))
			flag=false;
		else{
		
			//Write values of following into the file : a,b,pi,alpha,beta,delta,psi,gamma,Xi,observation/lamda,observation/lamdaNew,p* from Viterbi and qT* from Viterbi
					
			for(int i=0;i<N;i++)
				for(int j=0;j<N;j++)
					a[i][j]=aNew[i][j];

			for(int i=0;i<N;i++)
				for(int j=0;j<M;j++)
					b[i][j]=bNew[i][j];

			for(int i=0;i<N;i++)
				pi[i]=piNew[i];				
				
		}
						
	//std::cout<<std::endl<<"Reaching the end of while iteration"<<std::endl;	
	std::cout<<std::endl<<count++;
}while((flag) && (count < 500));

//write values of following into the file : a,b,pi,alpha,beta,delta,psi,gamma,Xi,observation/lamda,observation/lamdaNew,p* from Viterbi and qT* form Viterbi
*/
		

	getch();
	return 0;
}


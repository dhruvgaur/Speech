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
long double sum=0,sumNew=0;			
bool flag=true;

std::ofstream logfile("Output//FullLog.txt",std::ios_base::app);
std::ofstream output("Output//RequiredOutput.txt",std::ios_base::app);

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

	FILE *fp=fopen("Input//Initialization.txt","r");

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
	std::ifstream is("Input//Observation.txt");
	char c;
	while (is.get(c)) {
		if((c-'0' >= 0) && (c-'0' <= 9))T++;
	}	
		
	is.close();
	observation = new int[T];
	T=0;
	std::ifstream inp("Input//Observation.txt");
	while (inp.get(c)) {
		if((c-'0' >= 0) && (c-'0' <= 9))	observation[T++]=c-'0';
	
	}// loop getting single characters
		

	for(int i=0;i<T;i++){
		std::cout<<observation[i];
		//observation[i] -= 1;
	}
	inp.close();

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
					
		std::cout<<std::endl<<"Alpha :"<<std::endl<<std::endl;
		logfile<<std::endl<<"Alpha :"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){
			for(int j=0;j<T;j++){
				std::cout<<alpha[i][j]<<" : ";
				logfile<<alpha[i][j]<<" : ";
			}
		std::cout<<std::endl;
		logfile<<std::endl;
		}
	
					
		for(int i=0;i<N;i++)
			sum += alpha[i][T-1];

		std::cout<<std::endl<<"Probability of partial observation sequence given lamda [ P{O/lamda} ]= "<<sum;
		logfile<<std::endl<<"Probability of partial observation sequence given lamda [ P{O/lamda} ]=  "<<sum;
}

void Backward(){
	std::cout<<std::endl<<std::endl;
	logfile<<std::endl<<std::endl;
	

		for(int i=0;i<N;i++) beta[i][T-1]=1;

		for(int j=T-2;j>=0;j--)
			for(int i=0;i<N;i++){
				double temp=0;
				for(int k=0;k<N;k++)
					temp += a[i][k]*b[k][observation[j+1]]*beta[k][j+1];

				beta[i][j]=temp;
			}

			std::cout<<std::endl<<"Beta :"<<std::endl<<std::endl;
			logfile<<std::endl<<"Beta :"<<std::endl<<std::endl;


		for(int i=0;i<N;i++){
			for(int j=0;j<T;j++){

				std::cout<<beta[i][j]<<":";
				logfile<<beta[i][j]<<":";

			
			}
				
				std::cout<<std::endl;
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

		std::cout<<std::endl<<"Delta:"<<std::endl<<std::endl;
		logfile<<std::endl<<"Delta:"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){
			for(int j=0;j<T;j++){
				std::cout<<delta[i][j]<<" : ";
				logfile<<delta[i][j]<<" : ";			
			}
				
		std::cout<<std::endl;
		logfile<<std::endl;
		}

		//Displaying psi after the recurrsive step

		std::cout<<std::endl<<"Psi:"<<std::endl<<std::endl;
		logfile<<std::endl<<"Psi:"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){
			for(int j=0;j<T;j++){
				std::cout<<psi[i][j]<<" : ";

				logfile<<psi[i][j]<<" : ";
			}
				
		std::cout<<std::endl;
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
			

		std::cout<<std::endl<<"Probability of determined state sequence(P*)= "<<max;
		logfile<<std::endl<<"Probability of determined state sequence(P*)= "<<max;

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

		std::cout<<std::endl<<std::endl<<"Displaying the state sequence for lamda using Viterbi (Qt*):"<<std::endl;
		logfile<<std::endl<<std::endl<<"Displaying the state sequence for lamda using Viterbi (Qt*):"<<std::endl;

		for(int i=0;i<T;i++){
			std::cout<<sequence[i]<<" : ";
			logfile<<sequence[i]<<" : ";		
		}
			
}

void Baum_Welch(){

	std::cout<<std::endl<<std::endl;
	logfile<<std::endl<<std::endl;

		//Calculation the denominator beforehand for all the values of T

		for(int k=0;k<T-1;k++){
		
			denominator[k]=0;

			//Book variable : varable here(i:i,j:j,t:k)

			for(int i=0;i<N;i++)
				for(int j=0;j<N;j++)
					denominator[k] += (alpha[i][k]*a[i][j]*b[j][observation[k+1]]*beta[j][k+1]);
	
		}

		std::cout<<std::endl<<"Printing denominator for each T :"<<std::endl<<std::endl;
		logfile<<std::endl<<"Printing denominator for each T :"<<std::endl<<std::endl;

		for(int i=0;i<T-1;i++){
			std::cout<<denominator[i]<<" : ";
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
				

			std::cout<<std::endl<<"Xi:"<<std::endl;
			logfile<<std::endl<<"Xi:"<<std::endl;

			for(int i=0;i<N;i++){

				std::cout<<std::endl<<"layer #"<<i<<std::endl<<std::endl;
				logfile<<std::endl<<"layer #"<<i<<std::endl<<std::endl;

				for(int j=0;j<N;j++){

					for(int k=0;k<T-1;k++){
						std::cout<<Xi[i][j][k]<<" : "; 
						logfile<<Xi[i][j][k]<<" : ";
					
					}
						

					std::cout<<std::endl;
					logfile<<std::endl;
							
				}

			std::cout<<std::endl;	
			logfile<<std::endl;
						
			}

			std::cout<<std::endl<<"Gamma:"<<std::endl<<std::endl;
			logfile<<std::endl<<"Gamma:"<<std::endl<<std::endl;
			
			for(int i=0;i<N;i++){
			
				for(int j=0;j<T-1;j++)	{
					std::cout<<gamma[i][j]<<" : "; 
					logfile<<gamma[i][j]<<" : ";
				
				}									
							
					
				std::cout<<std::endl;
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

		//Printing new a,b and pi

		std::cout<<std::endl<<"aNew:"<<std::endl<<std::endl;
		logfile<<std::endl<<"aNew:"<<std::endl<<std::endl;
		output<<std::endl<<"aNew:"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){

			for(int j=0;j<N;j++){
				std::cout<<aNew[i][j]<<" : ";
				logfile<<aNew[i][j]<<" : ";
				output<<aNew[i][j]<<" : ";
			
			}
				

			std::cout<<std::endl;
			logfile<<std::endl;
			output<<std::endl;
					
		}

		std::cout<<std::endl<<"bNew:"<<std::endl<<std::endl;
		logfile<<std::endl<<"bNew:"<<std::endl<<std::endl;
		output<<std::endl<<"bNew:"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){
						
			for(int j=0;j<M;j++){
				std::cout<<bNew[i][j]<<" : ";
				logfile<<bNew[i][j]<<" : ";
				output<<bNew[i][j]<<" : ";
			
			}
				

			std::cout<<std::endl;
			logfile<<std::endl;
			output<<std::endl;
					
		}


		std::cout<<std::endl<<"piNew:"<<std::endl<<std::endl;	
		logfile<<std::endl<<"piNew:"<<std::endl<<std::endl;
		output<<std::endl<<"piNew:"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){
			std::cout<<piNew[i]<<" : ";
			logfile<<piNew[i]<<" : ";
			output<<piNew[i]<<" : ";
		
		}
			
}

void SubsequentForward(){
	output<<std::endl<<std::endl<<std::endl<<"(2) The probability of an utterance being from the model :"<<std::endl<<std::endl;

	for(int i=0;i<N;i++)alphaNew[i][0]=piNew[i]*bNew[i][observation[0]];

		for(int j=1;j<T;j++)
			for(int i=0;i<N;i++){
				double temp=0;
				for(int k=0;k<N;k++)
					temp += alphaNew[k][j-1]*aNew[k][i];
			
				alphaNew[i][j]=temp*bNew[i][observation[j]];
			
			}
					
		std::cout<<std::endl<<std::endl<<"AlphaNew :"<<std::endl<<std::endl;
		logfile<<std::endl<<std::endl<<"AlphaNew :"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){
			for(int j=0;j<T;j++){
				std::cout<<alphaNew[i][j]<<" : ";
				logfile<<alphaNew[i][j]<<" : ";
			
			}
				
		std::cout<<std::endl;
		logfile<<std::endl;
		}
	
					
		for(int i=0;i<N;i++)
			sumNew += alphaNew[i][T-1];

		std::cout<<std::endl<<"Probability of partial observation sequence given new lamda [ P{O/lamda} ]= "<<sumNew;
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

		std::cout<<std::endl<<std::endl<<"Delta:"<<std::endl<<std::endl;
		logfile<<std::endl<<std::endl<<"Delta:"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){
			for(int j=0;j<T;j++){

				std::cout<<delta[i][j]<<" : ";
				logfile<<delta[i][j]<<" : ";
			
			}
				
		std::cout<<std::endl;
		logfile<<std::endl;
		}

		//Displaying psi after the recurrsive step

		std::cout<<std::endl<<"Psi:"<<std::endl<<std::endl;
		logfile<<std::endl<<"Psi:"<<std::endl<<std::endl;

		for(int i=0;i<N;i++){
			for(int j=0;j<T;j++){

				std::cout<<psi[i][j]<<" : ";
				logfile<<psi[i][j]<<" : ";
			
			}
				
		std::cout<<std::endl;
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
			

		std::cout<<std::endl<<"Probability of determined state sequence(P*)= "<<maxTemp;
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

		std::cout<<std::endl<<std::endl<<"Displaying the state sequence for lamda using Viterbi (Qt*):"<<std::endl;
		logfile<<std::endl<<std::endl<<"Displaying the state sequence for lamda using Viterbi (Qt*):"<<std::endl;
		output<<std::endl<<std::endl<<"Displaying the state sequence for lamda using Viterbi (Qt*):"<<std::endl;

		for(int i=0;i<T;i++){
			std::cout<<sequenceTemp[i]<<" : ";
			logfile<<sequenceTemp[i]<<" : ";
			output<<sequenceTemp[i]<<" : ";
		
		}
			
}

int _tmain(int argc, _TCHAR* argv[])
{	
	std::cout<< std::scientific;
	std::cout<<std::setprecision(4);

	logfile<< std::scientific;
	logfile<<std::setprecision(4);

	output<< std::scientific;
	output<<std::setprecision(4);
	
	//Read data forom Initialization.txt file and initialize variables
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
		
							
		if(sumNew==sum)
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
						
	std::cout<<std::endl<<"Reaching the end of while iteration"<<std::endl;			
}while((flag) && (count<0));

//write values of following into the file : a,b,pi,alpha,beta,delta,psi,gamma,Xi,observation/lamda,observation/lamdaNew,p* from Viterbi and qT* form Viterbi

			

	getch();
	return 0;
}


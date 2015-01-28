// KMeans.cpp : Defines the entry point for the console application.

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

//Defines the size of codebook to generate from 
#define sizeofcodebook 32

//Define epsilon value for LBG
#define epsilon 0.03

//Define the pseudoepsilon value to split the largest codeBook vector in case of empty cell if any. However such scenarios never occured.
#define pseudoEpsilon 0.01


//functionMerge : Function to add two vectors and return the added vector

std::vector<long double> functionMerge(int number_of_dimensions,std::vector<long double> a, std::vector<long double> b){


	for(int j=0;j<number_of_dimensions;j++)
			a[j]+=b[j];
		
	return a;

}

//Tokhura : Function to calculate Tokura distance between two vectors

double Tokhura(std::vector<long double> test, std::vector<long double> reference){
	
	double distance=0;const double array[]={1,3,7,13,19,22,25,33,42,50,56,61,};
	for(int i=0;i<test.size();i++){
		
		distance+=array[i]*pow((test[i]-reference[i]),2);
	
	}
	//std::cout<<"Distance="<<distance;
	return distance;
}

//KMeans : Function to carry out KMeans iterations over a universe and a codeBook  

std::vector<std::vector<long double>> KMeans(std::vector<std::vector<long double>> *cepstralUniverse, std::vector<std::vector<long double>> codeBook, long double *distortion){
	
	std::vector<long double> temp1,temp2;

	std::vector<std::vector<long double>> clusteredUniverse;

	*distortion=0;

	std::vector<std::vector<long double>> clusteredUniverseXXXX;
	clusteredUniverseXXXX.clear();
	clusteredUniverseXXXX.resize(codeBook.size());
	
	std::ofstream logfile("Output//LBGLog.txt",std::ios_base::app);
	
	double *dist = new double[codeBook.size()];	
	
	int *counter = new int[codeBook.size()];
	std::vector<long double> dimensions=codeBook[0];
	int number_of_dimensions=dimensions.size();
	dimensions.clear();

	double *temp = new double[number_of_dimensions];	

	for(int j=0;j<codeBook.size();j++) counter[j]=0;
		
	//For every vector in the universe map it to the closest vector in the codeBook also keep updating distortion for every vector and count of vectors being mapped
	//to each cell in the codeBook.

	for(int i=0;i<cepstralUniverse->size();i++){

			for(int j=0;j<codeBook.size();j++) dist[j]=1000;

			double min;int indexMin;
			
			min=dist[0];indexMin=0;
			

			for(int j=0;j<codeBook.size();j++){

				dist[j]=Tokhura(cepstralUniverse->at(i),codeBook[j]);
				if(dist[j]<min){min=dist[j];indexMin=j;}

			}
			
			if(counter[indexMin]==0)
				clusteredUniverseXXXX[indexMin]=cepstralUniverse->at(i);
			else
				clusteredUniverseXXXX[indexMin]=functionMerge(number_of_dimensions,clusteredUniverseXXXX[indexMin],cepstralUniverse->at(i));

			
			counter[indexMin]++;
			*distortion += dist[indexMin];
			
	}	

	//Check for the empty cell situation and fix any by splitting the largest size bin.
	//Also normalize the sum of all vectors mapped to each bin by dividing by the total number of vectors getting mapped.

	int max;

	max=0;
	
	//Normalize each codeBook vector by dividing it with the size of vectors mapped to each bin

	for(int i=0;i<codeBook.size();i++){
		if(counter[i] != 0){

			temp1.clear();
			temp1=clusteredUniverseXXXX[i];

			for(int j=0;j<number_of_dimensions;j++){

				temp1[j] /= counter[i];
			}
			clusteredUniverseXXXX[i]=temp1;
		
		
		}
	
	}

	//Check for the empty cell situation and fix any by splitting the largest size bin.

	for(int i=0;i<codeBook.size();i++){
	
		if(counter[i] ==0){
			temp1.clear();
			temp2.clear();
			//Finding max sized codeBook entry
			for(int i=0;i<codeBook.size();i++)
				if(counter[i]>counter[max]) max=i;

			temp1=clusteredUniverseXXXX[max];
			temp2=temp1;
			
			//Splitting the largest sized codeBook entry
			for(int j=0;j<number_of_dimensions;j++){
			
				temp1[j] = temp1[j] * (1-pseudoEpsilon);
				temp2[j]=temp2[j] * (pseudoEpsilon);
			}
			clusteredUniverseXXXX[max]=temp1;

			clusteredUniverseXXXX[i]=temp2;
		}
		
	}

	for(int i=0;i<codeBook.size();i++){
		logfile<<"Code Book index #"<<i+1<<" : "<<counter[i]<<std::endl;
		std::cout<<"Code Book index #"<<i+1<<" : "<<counter[i]<<std::endl;
		
	}

	*distortion=(*distortion)/cepstralUniverse->size();

	return clusteredUniverseXXXX;
}

int _tmain(int argc, _TCHAR* argv[])
{
	   
std::vector<std::vector<long double>> cepstralUniverse;
std::vector<std::vector<long double>> codeBook;
std::vector<std::vector<long double>> clusteredUniverse;
std::vector<long double> temp;
long double distortion;
long double temp3;

char line[1000];
char* token;

//Reading the input file from which Cepstral coefficients are read and stored in a universe called : 'cepstralUniverse'.

FILE *fp=fopen("Input//Universe.csv","r");
int i;
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
	
}else std::cout << "Unable to open file"; 

std::ofstream logfile("Output//LBGLog.txt",std::ios_base::app);
std::ofstream output("Output//Output.txt",std::ios_base::app);

int r;
temp=cepstralUniverse[0];

//Generating mean values of all the data


for(int i=1;i<cepstralUniverse.size();i++){
	
	for(int j=0;j<temp.size();j++){
		temp[j] += cepstralUniverse[i][j];
	}
}

for(int j=0;j<temp.size();j++)
		temp[j] = temp[j]/cepstralUniverse.size();

codeBook.push_back(temp);
temp.clear();


if (logfile) {   

	logfile<<"-----------------Log file for KMeans Clustering-----------------"<<std::endl;
	logfile<<"Created by : Dhruv Gaur"<<std::endl<<"Roll # 144101068"<<std::endl;
	logfile<<"Description : Contains 'Code Book Vectors',' Distortion' and 'Number of Vectors in each bin of codebook' after every iteration"<<std::endl<<std::endl<<std::endl<<std::endl;

	std::cout<<"-----------------Log file for KMeans Clustering-----------------"<<std::endl;
	std::cout<<"Created by : Dhruv Gaur"<<std::endl<<"Roll # 144101068"<<std::endl;
	std::cout<<"Description : Contains 'Code Book Vectors',' Distortion' and 'Number of Vectors in each bin of codebook' after every iteration"<<std::endl<<std::endl<<std::endl<<std::endl;

	}
	else
		std::cout << "Unable to open file"; 

if (output) {   

	output<<"-----------------Output file for KMeans Clustering-----------------"<<std::endl;
	output<<"Created by : Dhruv Gaur"<<std::endl<<"Roll # 144101068"<<std::endl;
	output<<"Description : Contains final settled 'Code Book' and 'Distortion' after every iteration"<<std::endl<<std::endl<<std::endl<<std::endl;
	}
	else
		std::cout << "Unable to open file";


//Iterating for LBG algorithm until desired size of codebook is not reached.


int outeriteration=1;
do{
	logfile<<std::endl<<std::endl<<"----------------------Splitting Codebook and making it a size of #"<<pow(2.0,outeriteration)<<"----------------------"<<std::endl<<std::endl;
	std::cout<<std::endl<<std::endl<<"----------------------Splitting Codebook and making it a size of #"<<pow(2.0,outeriteration)<<"----------------------"<<std::endl<<std::endl;
	output<<std::endl<<std::endl<<"----------------------Splitting Codebook and making it a size of #"<<pow(2.0,outeriteration++)<<"----------------------"<<std::endl<<std::endl;

	

	//Splitting the codebook

	int size=codeBook.size();
	for(int i=0;i<size;i++){
		temp.clear();
		temp=codeBook[i];
		for(int j=0;j<temp.size();j++){
			codeBook[i][j] *= (1+epsilon);
			temp[j] *= (1-epsilon);
		}
		
		codeBook.push_back(temp);
	}

	//Iterating for inner KMeans algorithm

	int iteration=1;
	distortion=0;
	do{
		temp3=distortion;

		logfile<<"Number of vectors from Universe mapping to each Code Book index during iteration #"<<iteration<<" : "<<std::endl<<std::endl;
		std::cout<<"Number of vectors from Universe mapping to each Code Book index during iteration #"<<iteration<<" : "<<std::endl<<std::endl;

		clusteredUniverse=KMeans(&cepstralUniverse,codeBook,&distortion);

		logfile<<std::endl<<"Code Book after iteration #"<<iteration<<" : "<<std::endl<<std::endl;
		std::cout<<std::endl<<"Code Book after iteration #"<<iteration<<" : "<<std::endl<<std::endl;

		//Printing codebook vectors

		for(int i=0;i<clusteredUniverse.size();i++){
			temp.clear();
			temp=clusteredUniverse[i];
			
				for(int j=0;j<temp.size();j++)	
				if(j != temp.size()-1){
					if(j==0){
						logfile<<"Entry #"<<i+1<<" : "<<temp[j]<<":";
						std::cout<<"Entry #"<<i+1<<" : "<<temp[j]<<":";
					}
						
					else{
						logfile<<temp[j]<<":";
						std::cout<<temp[j]<<":";
					}
						
				}				
				else
					logfile<<temp[j];

			logfile<<std::endl;
					
		
		}

		logfile<<std::endl;
	
		logfile<<"Distortion after iteration #"<<iteration<<" : "<<distortion<<std::endl<<std::endl;
		std::cout<<"Distortion after iteration #"<<iteration<<" : "<<distortion<<std::endl<<std::endl;
		if(temp3==0){

			logfile<<"Change in distortion after iteration #"<<iteration<<" : "<<"This is first iteration"<<std::endl<<std::endl;
		}
		else{

			logfile<<"Change in distortion after iteration #"<<iteration<<" : "<<((temp3-distortion)*100)/temp3<<"%"<<std::endl<<std::endl;
		
		}
	
		output<<"Distortion after iteration #"<<iteration++<<" : "<<distortion<<std::endl;


		codeBook.clear();
		codeBook=clusteredUniverse;

		std::cout<<std::endl<<clusteredUniverse.size()<<"::"<<distortion<<std::endl;


		clusteredUniverse.clear();


	}while((temp3-distortion) != 0);

} while(codeBook.size() < sizeofcodebook);
output<<std::endl<<"Code Book after final settlement is:"<<std::endl<<std::endl;

for(int i=0;i<codeBook.size();i++){
		temp.clear();
		temp=codeBook[i];
			
			for(int j=0;j<temp.size();j++)	
			if(j != temp.size()-1){
				if(j==0)
					output<<"Entry #"<<i+1<<" : "<<temp[j]<<":";
				else
					output<<temp[j]<<":";
			}				
			else
				output<<temp[j];

		output<<std::endl;
					
		
	}
	
std::cout<<"Reaching end of execution";
getch();
return 0;
}//End of main()


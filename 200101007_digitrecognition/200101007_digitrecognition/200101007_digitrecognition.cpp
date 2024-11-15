// 200101007_digitrecognition.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <numeric>

#define ld long double
#define vvld vector<vector<ld>>
#define vld vector<ld>
#define vint vector<int>
#define pi 3.142857142857

using namespace std;

string uni = "Universe.csv";
int framesize = 320;
int frames = 50;
int normalizemax = 5000;
ld w[12] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
int p = 12;
int M = 8;
int N = 5;
int T = 50;
ld const threshold = 1e-30;
int environment_known = 0;
int is_live_testing = 0;


vint readData(string &filename){
	vint data;

    // Open the input text file
    ifstream inputFile(filename);

    // Check if the file is open
    if (!inputFile.is_open()) {
        cerr << "Failed to open the file." << endl;
		return data;
    }

    int value;
    while (inputFile >> value) {
        data.push_back(value);
    }
    inputFile.close();

	return data;
}

//Function to read the values from the csv file
vector<vector<ld>> readcsv(string fname){
	vector<vector<ld>> content;
	vector<ld> row;
	string line, word;
 
	fstream file (fname, ios::in);
	if(file.is_open())
	{
		while(getline(file, line))
		{
			row.clear();
 
			stringstream str(line);
 
			while(getline(str, word, ',')){
				ld num = stold(word);
				row.push_back(num);
			}
			content.push_back(row);
		}
	}
	return content;
}

vector<ld> read1dcsv(string fname){
	vector<ld> content;
 
	fstream file (fname, ios::in);
	string line;
    while (getline(file, line)) {
        istringstream iss(line);
        string token;
        while (getline(iss, token, ',')) {
            long double value;
            stringstream converter(token);
            converter >> value;
            content.push_back(value);
        }
    }
	return content;
}

void printcsv(vld &data, string filename){
	ofstream csvFile(filename, ios::app);
	for (int i = 0; i < data.size(); ++i) {
        csvFile << data[i];
        if (i < data.size() - 1) {
            csvFile << ",";
        }
    }
	csvFile << endl;
	csvFile.close();
}

vld dcshift(vint &data){
	ld dcshift = accumulate(data.begin(), data.end(), 0.0)/data.size();
	vld dcdata(frames*framesize,0);
	for(int i = 0; i<data.size(); i++){
		dcdata[i] = data[i]-dcshift;
	}
	return dcdata;
}

void normalize(vld &data){
	ld mx = 0;
	for(int i = 0; i<data.size(); i++){
		mx = max(data[i],mx);
	}
	for(int i = 0; i<data.size(); i++){
		data[i] = data[i]*normalizemax/mx;
	}
}

vld correction(vint &data){
	vld correcteddata = dcshift(data);
	normalize(correcteddata);
	return correcteddata;
}

//function to apply raised sine window on the cepstral coefficients
vector<ld> applyrsw(vector<ld> ci){

	for(int i = 1; i<=p; i++){
		ci[i]*=1+(p/2)*sin(pi*i/p);
	}

	return ci;
}

//Function for calculation of Ri's
vector<ld> Rcalc(vector<ld> input){

	vector<ld> R(p+1,0.0);

	for(int i = 0; i<=p; i++){
		for(int j = 0; j+i<input.size(); j++){
			R[i]+=input[j]*input[j+i];
		}
		R[i]/=input.size();
	}

	return R;
}

//Function for calculation of the ai's
vector<ld> acalc(vector<ld> R){

	vector<ld> Ei(p+1,0.0);
	vector<ld> ki(p+1,0.0);
	vector<vector<ld>> alpha(p+1,vector<ld>(p+1,0.0));
	vector<ld> ai(p+1,0.0);

	Ei[0]=R[0];
	
	for(int i = 1; i<=p; i++){
		if(i==1){
			ki[i]=R[i]/Ei[i-1];
		}
		else{
			ld alphasum = 0.0;
			for(int j = 1; j<i; j++){
				alphasum+=alpha[i-1][j]*R[i-j];
			}
			ki[i] = (R[i]-alphasum)/Ei[i-1];
		}
		alpha[i][i]=ki[i];
		for(int j = 1; j<i; j++){
			alpha[i][j]=alpha[i-1][j]-ki[i]*alpha[i-1][i-j];
		}
		Ei[i]=(1-ki[i]*ki[i])*Ei[i-1];
	}

	for(int i = 1; i<=p; i++){
		ai[i]=alpha[p][i];
	}

	return ai;
}

//Function for calculation of the cepstral coefficients ci's
vector<ld> cicalc(vector<ld> ai){

	vector<ld> ci(p+1,0.0);

	for(int i = 1; i<=p; i++){
		if(i==1){
			ci[i]=ai[i];
		}
		else{
			ld csum = 0.0;
			for(int j = 1; j<i; j++){
				csum+=j*ci[j]*ai[i-j]/i;
			}
			ci[i]=ai[i]+csum;
		}
	}

	return ci;
}

void calc_universe(){
	for(long long digit = 0; digit<10; digit++){
		for(long long i = 1; i<=20; i++){
			string filename = "English/txt/" + to_string(digit) + "/200101007_E_" + to_string(digit) + "_" + to_string(i) + ".txt";
			vector<int> data = readData(filename);
			vector<ld> correcteddata = correction(data);
			vector<vector<ld>> input(frames,vector<ld>(framesize,0.0));
			for(int j = 0; j<correcteddata.size(); j++){
				input[j/framesize][j%framesize]=correcteddata[j];
			}
	
			for(int j = 0; j<frames; j++){
				vector<ld> R = Rcalc(input[j]);
				vector<ld> ai = acalc(R);
				vector<ld> ci = cicalc(ai);
				vector<ld> cirsw = applyrsw(ci);
				cirsw.erase(cirsw.begin());
				printcsv(cirsw,uni);
			}
		}
	}
}

//Function to calculate the Tokhura's Distance for a given pair of test and reference vector
ld tokhura(vector<ld> test, vector<ld> ref){
	ld dist = 0.0;
	for(int i = 0; i<p; i++){
		dist+=w[i]*(test[i]-ref[i])*(test[i]-ref[i]);
	}
	return dist;
}

//Function to assign every vector in the universe to the nearest region
//Output is a 2D vector with k rows containing the index of vectors in the respective region
vector<vector<int>> assignment(vector<vector<ld>> &universe, vector<vector<ld>> &codebook, int k){
	vector<vector<int>> cbassign(k,vector<int>());
	for(int i = 0; i<universe.size(); i++){
		ld dist = INT_MAX;
		int index = 0;
		for(int j = 0; j<k; j++){
			ld temp = tokhura(universe[i],codebook[j]);
			if(temp < dist) {
				dist = temp;
				index = j;
			}
		}
		cbassign[index].push_back(i);
	}
	return cbassign;
}

//Calculate the total distortion of every region and averaged
ld distortion(vector<vector<ld>> &universe, vector<vector<ld>> &codebook, vector<vector<int>> &cbassign){
	ld dist = 0.0;
	for(int i = 0; i<cbassign.size(); i++){
		for(int j = 0; j<cbassign[i].size(); j++){
			dist+=tokhura(universe[cbassign[i][j]], codebook[i])/universe.size();
		}
	}
	return dist;
}

//Calculate the new centroids based on the assignment of the vectors
vector<vector<ld>> centroid(vector<vector<ld>> &universe, vector<vector<int>> &cbassign, int k){
	vector<vector<ld>> codebook(k,vector<ld>(p,0.0));
	for(int i = 0; i<k; i++){
		for(int j = 0; j<p; j++){
			for(int l = 0; l<cbassign[i].size(); l++){
				codebook[i][j]+=universe[cbassign[i][l]][j]/cbassign[i].size();
			}
		}
	}
	return codebook;
}

vector<vector<ld>> kmeans(vector<vector<ld>> &universe, vector<vector<ld>> &codebook, int k){
	int m = 0;
	vector<vector<int>> cbassign = assignment(universe,codebook, k); //Assigning the region to the vectors in the universe based on the codebook
	
	ld dist = distortion(universe, codebook, cbassign); //Calculate the avg distortion for the codebook

	ld prev = INT_MAX; //Distortion for the previous iteration initialized to a large value
	ld delta = 0.00001; //Threshold

	codebook = centroid(universe, cbassign, k); //Updating the codebook with the new centroid
	m++;
	
	//Repeating these steps until the change in distortion is less than threshold
	while(abs(dist-prev)>delta){
		cbassign.clear();
		cbassign = assignment(universe,codebook, k);

		prev = dist;
		dist = distortion(universe, codebook, cbassign);

		codebook.clear();
		codebook = centroid(universe, cbassign, k);

		m++;
	}
	return codebook;
}

void create_codebook(){
	vector<vector<ld>> universe = readcsv(uni);
	int k = 1; //k is the current codebook size
	vector<vector<ld>> codebook;
	vector<vector<int>> cbassign(k,vector<int>(universe.size(),0)); //Assigning all the elements to the first and only region of the codebook
	codebook = centroid(universe, cbassign, k);// Initializing the codebook to the centroid of all the elements in the universe
	ld ep = 0.03;

	while(k<M){
		//Creating a temporary codebook storing the yi +- epsilon values
		vector<vector<ld>> tempcodebook;
		for(int i = 0; i<k; i++){
			vector<ld> temp;
			for(int j = 0; j<p; j++){
				temp.push_back(codebook[i][j]+ep);
			}
			tempcodebook.push_back(temp);
			temp.clear();
			for(int j = 0; j<p; j++){
				temp.push_back(codebook[i][j]-ep);
			}
			tempcodebook.push_back(temp);
		}
		k*=2;
		codebook.clear();
		//Applying kmeans to find the new and final codebook for the current value of k
		codebook = kmeans(universe, tempcodebook, k);
	}
	for(int i = 0; i<codebook.size(); i++){
		printcsv(codebook[i],"codebook.csv");
	}
}

int calc_obs_seq(vld &ci, vvld &codebook){
	
	ld dist = INT_MAX;
	int index = 0;
	for(int j = 0; j<M; j++){
		ld temp = tokhura(ci,codebook[j]);
		if(temp < dist) {
			dist = temp;
			index = j;
		}
	}
	return index+1;
}

void create_initial_lambda(){

			vld pie;
			pie.push_back(0);
			pie.push_back(1.0);
			for(int i = 1; i<N; i++){
				pie.push_back(0.0);
			}
			printcsv(pie,"pie.csv");
			vvld a(N+1,vld(N+1,0.0));
			for(int i = 1; i<N; i++){
				a[i][i]=0.8;
				a[i][i+1]=0.2;
			}
			a[N][N]=1.0;
			for(int i = 0; i<N+1; i++){
				printcsv(a[i],"a.csv");
			}
			vvld b(N+1,vld(M+1,0.125));
			for(int i = 0; i<=N; i++){
				b[i][0]=0;
				printcsv(b[i],"b.csv");
			}
}

void import_lambda(vvld &a, vvld&b, vld &pie){
	pie = read1dcsv("pie.csv");
	a = readcsv("a.csv");
	b = readcsv("b.csv");
}

vvld forward_procedure(vvld &A, vvld &B, vld &Pi, vint &O){
	int i , j , t;
	long double sum ;
	vvld Alpha(T+1,vld(N+1,0));
	int index = O[1];
	

	for(i=1;i<=N;i++){
		Alpha[1][i] = Pi[i]*B[i][index];
	}
	
	for (t = 1; t < T; t++){
		for (j = 1; j <= N; j++){
			sum = 0;
			for (i = 1; i <= N; i++){
				sum += Alpha[t][i] * A[i][j];
			}
			Alpha[t + 1][j] = sum * B[j][O[t + 1]];
		}
	}
	
	
	return Alpha;
}

//Calculation of Beta variable.
vvld backward_procedure(vvld &A, vvld &B, vld &Pi, vint &O){
	int i , j , t;
	long double sum;
	vvld Beta(T+1,vld(N+1,0));
	int index = 0;
	for(i=1;i<=N;i++){
		Beta[T][i] = 1.0;
	}
	for(t=T-1;t>=1;t--){
		index = O[t+1];
		for(i=1;i<=N;i++){
			sum = 0;
			for(j=1;j<=N;j++){
				sum = sum + B[j][index]*A[i][j]*Beta[t+1][j];
			}
			Beta[t][i]=sum;
		}
	}
	return Beta;
}

//viterbi algorithm
ld viterbi(vvld &A, vvld &B, vld &Pi, vint &O){
	//initialization
	vvld Delta(T+1,vld(N+1,0));
	vvld Psi(T+1,vld(N+1,0));
    for(int i=1; i<=N; i++){
        Delta[1][i] = Pi[i] * B[i][O[1]];
        Psi[1][i] = 0;
    }

	//induction
	for(int j=1; j<=N; j++){
		for(int t=2; t<=T; t++){
            long double max = 0, ti = 0;
            int ind = 0;
            
            for(int i=1; i<=N; i++){
                ti = Delta[t-1][i] * A[i][j];
                if(ti > max){
					max = ti;
					ind = i;
				}
            }

            Delta[t][j] = max * B[j][O[t]];
			Psi[t][j] = ind;
        }
    }

	vint Q(T+1,0);
	ld pstar;

	//termination
    long double max = 0;
    for(int i=1; i<=N; i++){
        if(Delta[T][i] > max) {
			max = Delta[T][i];
			Q[T] = i;
		}

        pstar = max;
    }

	//backtracking
    for(int t = T-1; t>0; t--){
        Q[t] = Psi[t+1][Q[t+1]];
    }
	return pstar;
}

vector<vvld> calculate_xi(vvld &Alpha, vvld &Beta, vvld &A, vvld &B, vld &Pi, vint &O){
	int i , j , t , index;
	vld summation(T+1,0);
	vector<vvld> Xi(T+1,vvld(N+1,vld(N+1,0)));
	for(t=1;t<T;t++){
		summation[t] = 0;
		for(i=1;i<=N;i++){
			for(j=1;j<=N;j++){
				summation[t] += Alpha[t][i]*A[i][j]*B[j][O[t+1]]*Beta[t+1][j];
			}
		}
		
		for(i=1;i<=N;i++){
			long double x;
			for(j=1;j<=N;j++){
				x = Alpha[t][i]*A[i][j]*B[j][O[t+1]]*Beta[t+1][j];
				Xi[t][i][j]= x/summation[t];
			}
		}
	}
	return Xi;
}

//calculating alpha and beta values
vvld calculate_gamma(vvld &Alpha, vvld &Beta){
	vvld Gamma(T+1,vld(N+1,0));
	for(int t=1;t<=T;t++){
		for(int j=1;j<=N;j++){
			long double summation=0;
			for(int k=1;k<=N;k++){
				summation += Alpha[t][k] * Beta[t][k];
			}
			Gamma[t][j]=(Alpha[t][j] * Beta[t][j])/summation;
		}
	}
	return Gamma;
}

void reevaluate_model_parameters(vector<vvld> &Xi, vvld &Gamma, vvld &Abar, vvld &Bbar, vld &Pibar, vint &O){
	int i, j, k, t;
	long double sum1=0 , sum2 =0;
	//Re-evaluating Pi
	for(i=1;i<=N;i++){
		Pibar[i] = Gamma[1][i];
	}
	
	for(i = 1; i<=N; i++){
		for(j = 1; j <= N; j++){
			long double t1 = 0, t2 = 0;
			for(t = 1; t <= T-1; t++){
				t1 += Xi[t][i][j];
				t2 += Gamma[t][i];
			}
			Abar[i][j] = t1/t2;
		}
	}
	//Re-evaluating B
	for(j=1;j<=N;j++){
		int count=0;
		long double max=0;
		int ind_j=0, ind_k=0;
		
		for(k=1;k<=M;k++){
			sum1 =0 , sum2 =0;
			for(t=1;t<T;t++){
				sum1 = sum1 + Gamma[t][j];
				if(O[t]==k){
					sum2 = sum2 + Gamma[t][j];				
				}
			}
			Bbar[j][k] = sum2/sum1;
			
			//finding max
			if(Bbar[j][k]>max){
				max=Bbar[j][k];
				ind_j = j;
				ind_k = k;
			}
			
			//updating new bij with threshold value if it is zero
			if(Bbar[j][k] == 0){
				Bbar[j][k]=threshold;
				count++;
			}
		}
		Bbar[ind_j][ind_k] = max - count*threshold;
	}
}

void export_lambda(vvld &A, vvld &B, vld &Pi, long long digit, long long i){
	string fname = "output/lambda/" + to_string(digit) + "/" + to_string(i) + "_";
	printcsv(Pi,fname+"pie.csv");
	for(int i = 0; i<N+1; i++){
		printcsv(A[i],fname+"a.csv");
	}
	for(int i = 0; i<=N; i++){
		printcsv(B[i],fname+"b.csv");
	}
}

void print_final(){
	vvld a_final(N+1,vld(N+1,0));
	vvld b_final(N+1,vld(M+1,0));
	vld pi_final(N+1,0);
	for(long long digit = 0; digit<10; digit++){
		for(long long i = 1; i<=20; i++){
			string fname = "output/lambda/" + to_string(digit) + "/" + to_string(i) + "_";
			vvld a = readcsv(fname+"a.csv");
			vvld b = readcsv(fname+"b.csv");
			vld pie = read1dcsv(fname+"pie.csv");
			for(int j = 1; j<=N; j++){
				for(int l = 1; l<=N; l++){
					a_final[j][l] += a[j][l]/20;
				}
			}
			
			for(int j = 1; j<=N; j++){
				for(int l = 1; l<=M; l++){
					b_final[j][l] += b[j][l]/20;
				}
			}
			for(int j = 1; j<N; j++){
				pi_final[j] += pie[j]/20;
			}
		}
		
		string fname = "output/lambda/" + to_string(digit) + "/";
		printcsv(pi_final,fname+"pie_final.csv");
		for(int j = 0; j<N+1; j++){
			printcsv(a_final[j],fname+"a_final.csv");
		}
		for(int j = 0; j<=N; j++){
			printcsv(b_final[j],fname+"b_final.csv");
		}

		for(int j = 1; j<=N; j++){
				for(int l = 1; l<=N; l++){
					a_final[j][l] = 0;
				}
			}
			
			for(int j = 1; j<=N; j++){
				for(int l = 1; l<=M; l++){
					b_final[j][l] = 0;
				}
			}
			for(int j = 1; j<N; j++){
				pi_final[j] = 0;
			}
	}
}

vld calc_ci(vld &input){
	vector<ld> R = Rcalc(input);
	vector<ld> ai = acalc(R);
	vector<ld> ci = cicalc(ai);
	vector<ld> cirsw = applyrsw(ci);
	cirsw.erase(cirsw.begin());
	return cirsw;
}


void train_model(){
	vvld codebook = readcsv("codebook.csv");
	for(long long digit = 0; digit<10; digit++){
		for(long long i = 1; i<=20; i++){
			string filename = "English/txt/" + to_string(digit) + "/200101007_E_" + to_string(digit) + "_" + to_string(i) + ".txt";
			vector<int> data = readData(filename);
			vector<ld> correcteddata = correction(data);
			vector<vector<ld>> input(frames,vector<ld>(framesize,0.0));
			for(int j = 0; j<correcteddata.size(); j++){
				input[j/framesize][j%framesize]=correcteddata[j];
			}
			vint obs_seq;
			obs_seq.push_back(0);
			for(int j = 0; j<frames; j++){
				vld ci = calc_ci(input[j]);
				obs_seq.push_back(calc_obs_seq(ci,codebook));
			}
			vvld a, b;
			vld pie;
			import_lambda(a,b,pie);
			ld pstar = 0, prev_pstar = -1;
			int iteration = 1;
			
			while(pstar > prev_pstar && iteration < 200){
				iteration++;
				prev_pstar = pstar;
				
				vvld Alpha = forward_procedure(a,b,pie,obs_seq);
				ld P_O_given_Model = 0;
				
				for(int j=1;j<=N;j++){
					P_O_given_Model = P_O_given_Model + Alpha[T][j];
				}
				
				vvld beta = backward_procedure(a,b,pie,obs_seq);
				pstar = viterbi(a,b,pie,obs_seq);
				
				vector<vvld> Xi = calculate_xi(Alpha, beta, a, b, pie, obs_seq);
				vvld gamma = calculate_gamma(Alpha, beta);
				reevaluate_model_parameters(Xi,gamma, a, b, pie, obs_seq);
			}
			cout << digit << "-" << i << endl;
			export_lambda(a,b,pie,digit,i);
		}
	}
	print_final();
}

void import_final_lambda(vvld &a, vvld&b, vld &pie, long long x){
	string fname = "output/lambda/" + to_string(x) + "/";
	pie = read1dcsv(fname+"pie_final.csv");
	a = readcsv(fname+"a_final.csv");
	b = readcsv(fname+"b_final.csv");
}

void testing(){
	vvld codebook = readcsv("codebook.csv");
	int total_correct = 0;
	for(long long digit = 0; digit<10; digit++){
		int correct = 0;
		for(long long i = 21; i<=30; i++){
			
			string filename = "English/txt/" + to_string(digit) + "/200101007_E_" + to_string(digit) + "_" + to_string(i) + ".txt";
			vector<int> data = readData(filename);
			vector<ld> correcteddata = correction(data);
			vector<vector<ld>> input(frames,vector<ld>(framesize,0.0));
			for(int j = 0; j<correcteddata.size(); j++){
				input[j/framesize][j%framesize]=correcteddata[j];
			}
			vint obs_seq;
			obs_seq.push_back(0);
			for(int j = 0; j<frames; j++){
				vld ci = calc_ci(input[j]);
				obs_seq.push_back(calc_obs_seq(ci,codebook));
			}
			ld prob = 0;
			ld pred = 0;
			for(long long x = 0; x<10; x++){
				vvld a, b;
				vld pie;
				import_final_lambda(a,b,pie,x);
				
				vvld Alpha = forward_procedure(a,b,pie,obs_seq);
				ld P_O_given_Model = 0;
				
				for(int j=1;j<=N;j++){
					P_O_given_Model = P_O_given_Model + Alpha[T][j];
				}
				
				if(P_O_given_Model>prob){
					prob = P_O_given_Model;
					pred = x;
				}
			}
			if(pred == digit) {
				total_correct++;
				correct++;
			}
			cout << "Actual Value : " << digit << " Predicted Value : " << pred << endl;  
		}
		cout << endl << "Accuracy for " << digit << " = " << correct*10 << "%" << endl << endl;
	}
	cout << "Accuracy of the Model = " << total_correct << "%" << endl << endl;
}

void test_prediction(string filename){
	vvld codebook = readcsv("codebook.csv");
	vector<int> data = readData(filename);
	vector<ld> correcteddata = correction(data);
	vector<vector<ld>> input(frames,vector<ld>(framesize,0.0));
	for(int j = 0; j<correcteddata.size(); j++){
		input[j/framesize][j%framesize]=correcteddata[j];
	}
	vint obs_seq;
	obs_seq.push_back(0);
	for(int j = 0; j<frames; j++){
		vector<ld> R = Rcalc(input[j]);
		vector<ld> ai = acalc(R);
		vector<ld> ci = cicalc(ai);
		vector<ld> cirsw = applyrsw(ci);
		cirsw.erase(cirsw.begin());
		obs_seq.push_back(calc_obs_seq(cirsw,codebook));
	}
	ld prob = 0;
	ld pred = 0;
	for(long long x = 0; x<10; x++){
		vvld a, b;
		vld pie;
		import_final_lambda(a,b,pie,x);
		vvld Alpha = forward_procedure(a,b,pie,obs_seq);
		ld P_O_given_Model = 0;
				
		for(int j=1;j<=N;j++){
			P_O_given_Model = P_O_given_Model + Alpha[T][j];
		}
				
		if(P_O_given_Model>prob){
			prob = P_O_given_Model;
			pred = x;
		}
	}
	cout << "Prediction : " << pred << endl;  
}

void live_testing(){
	system("Recording_Module\\Recording_Module.exe 1 input.wav input_file.txt");	
	test_prediction("input_file.txt");
}

int _tmain(int argc, _TCHAR* argv[])
{
	char choice;
	
	while(1){
		cout<<"Press 1. for automated test on test files\nPress 2. for live testing\nEnter your choice: "; cin>>choice;
		switch(choice){
			case '1':
				{
					testing();
					break;
				}
			case '2':
				{
					is_live_testing = 1;
					live_testing();
					is_live_testing = 0;
					break;
				}
			case '0':
				{
					return 0;
				}
		}
	}

	//train_model();
	
	return 0;
}



// YesNoDetection.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include<iostream>
#include<string>
#include<fstream>
#include<vector>
#include<cmath>
#include<unordered_map>

#define WINDOW_SIZE 320
#define MAX_AMP 10000
#define P 12
#define Q 12
#define THRESH_FACTOR 20
#define IGNORE_SAMPLES 6300				//you can observe that after 6300 samples sample start recording noise
#define SAMPLE_FILE "input_file.txt"
#define SILENCE_SAMPLE 1600
#define COMPRESSED_FILE "compress.txt"

#define PRE_PATH_INPUT "Input_Vowels/170101022_"
#define PRE_PATH_OUTPUT "Output_Vowels/170101022_Out_"
#define PRE_PATH_COFF "Coefficient_Vowels/170101022_Coff_"

#define PRE_PATH_TEST "Test_Vowels/170101022_test_"
#define PRE_PATH_TEST_OUTPUT "Test_Output_Vowels/170101022_Out_test_"
#define PRE_PATH_TEST_COFF "Test_Coefficient_Vowels/170101022_Coff_test_"

#define AVERAGE_CEPS "Average_Ceps/average_ceps_"
#define CODE_TABLE "Code_Table.txt"
#define COUNT_STEADY_WINDOW 5

#define RESULT_FILE "Result.txt"

long double Tokhura_Weights[] = {1.0,3.0,7.0,13.0,19.0,22.0,25.0,33.0,42.0,50.0,56.0,61.0};
using namespace std;

const long double PI = 3.141592653589793238;


double ShiftandNormalize(double sampleValue,double dcShift,double normalFactor)
{
	return (sampleValue-dcShift)*normalFactor;
}

double max(double s1,double s2)
{
	return (s1>s2?s1:s2);
}

double CalculateDCShift(char c,int y, string filePath)
{
	string str = to_string(long long(y));
	string inFileName = filePath;
	inFileName = inFileName + c + "_" + str + ".txt";
	ifstream infile;							
	infile.open(inFileName);

	double dcShift= 0;

	if(infile.is_open())
	{
		int countIgnore =0;
		while(!infile.eof() && countIgnore<IGNORE_SAMPLES)
		{
			int sampleValue;
			infile >> sampleValue;
			countIgnore++;
		}

		int countSilence = 0; 

		while(countSilence < SILENCE_SAMPLE && !infile.eof())
		{
			int sampleValue;
			infile >> sampleValue;

			dcShift += sampleValue/double(SILENCE_SAMPLE); 
			countSilence++;
		}
	}

	infile.close();
	return dcShift;
}

pair<double,double> CalculateMaxAmpandSilenceSTE(double dcShift,char c,int y,string filePath)
{
	string str = to_string(long long(y));
	string inFileName = filePath;
	inFileName = inFileName + c + "_" + str + ".txt";
	cout << inFileName << endl;
	ifstream infile;							
	infile.open(inFileName);

	/*ifstream infile;								
	infile.open(SAMPLE_FILE);*/
	
	double steSilence=0;
	double maxAmp = 0;

	if(infile.is_open())
	{
		int countIgnore =0;
		while(!infile.eof() && countIgnore<IGNORE_SAMPLES)
		{
			int sampleValue;
			infile >> sampleValue;
			countIgnore++;
		}

		int countSilence = 0; 
		while(countSilence < SILENCE_SAMPLE && !infile.eof())
		{
			double sampleValue;
			infile >> sampleValue;

			sampleValue = sampleValue-dcShift;

			steSilence += (sampleValue*sampleValue)/double(SILENCE_SAMPLE);
			maxAmp = max(maxAmp,abs(sampleValue));
			countSilence++;
		}

		while(!infile.eof())
		{
			double sampleValue;
			infile >> sampleValue;

			sampleValue = sampleValue-dcShift;

			maxAmp = max(maxAmp,abs(sampleValue));
		}
	}
	infile.close();
	return make_pair(steSilence,maxAmp);
}
void CalculateAs(vector<double> &R,vector<double> &alpha)
{
	alpha.resize(P+1);
	fill(alpha.begin(),alpha.end(),0);

	vector<double> E;
	E.resize(P+1);

	vector<double> k;
	k.resize(P+1);

	vector<double> temp_alp(P+1,0);

	E[0] = R[0];
	for(int i=1;i<=P;i++)
	{
		if(i==1)
		{
			k[1]=R[1]/R[0];
		}
		else
		{
			double sum=0;
			for(int j=1;j<=i-1;j++)
			{
				temp_alp[j] = alpha[j];
				sum += temp_alp[j]*R[i-j];
			}
			k[i]=(R[i]-sum)/E[i-1];
		}

		alpha[i] = k[i];
		for(int j=1;j<=i-1;j++)
		{
			alpha[j] = temp_alp[j]-k[i]*temp_alp[i-j];
		}
		E[i]=(1-k[i]*k[i])*E[i-1];
	}
	return;
}
void CalculateRs(vector<double> window,vector<double> &R)
{
	for(int i=0;i<=P;i++)
	{
		double temp =0;
		for(int j=0;j+i<int(window.size());j++)
		{
			temp += (window[j]*window[i+j])/(window.size());
		}
		R.push_back(temp);
	}
	return;
}
void cepstralcoefficients(vector<double> R,vector<double> alpha,vector<long double> &ceps)
{
	for(int i=0;i<alpha.size();i++)
	{
		alpha[i] = -alpha[i];
	}
	ceps.resize(Q+1);
	unsigned int m,j;
	
	ceps[0]= logl(R[0]); //initial cepstral coefficent computation from energy
	for(m=1;m<=P;m++)
	{
		long double sum=0;
		for(j=1;j<=m-1;j++) //calculate the sum from older cepstral coefficents to compute new cepstral coefficients
		{
			sum += (j/long double(m))*(ceps[j]*alpha[m-j]);
		}
		ceps[m]=alpha[m]+sum; //new cepstral coefficients
	}
	if(m>P) // for our assignment this never get executed as we assume q=p
	{
		for(m=P+1;m<=Q;m++)
		{
			long double sum=0;
			for(j=m-P;j<=m-1;j++)
			{
				sum+=(j/long double(m))*(ceps[j]*alpha[m-j]);
			}
			ceps[m]=sum;
		}
	}
	//cout << "\nCepstral Coefficient values\n";//write the resuls to file
	//for(int i=0;i<=Q;i++)
	//{
	//	cout << ceps[i] << " ";
	//}
	//cout << endl;
	return;
}
vector<long double> RaisedSineWeights()
{
	vector<long double> weights;
	for(int i=0;i<=Q;i++)
	{
		long double temp = 1 + (Q/2)*sin((PI*i)/Q);
		weights.push_back(temp);
	}
	return weights;
}
pair<vector<double>,vector<double>> CalculateandFindCoefficients(vector<double> &wordWindow)
{
	vector<double> R;
	CalculateRs(wordWindow,R);

	vector<double> alpha;	
	if(R[0]!=0)
	{
		CalculateAs(R,alpha);
	}

	return make_pair(R,alpha);
}
int findSteadyIndex(vector<vector<double>> word)
{
	int windowIndex = word.size()/2;
	double maxAmp = 0;
	for(int i=0;i<word.size();i++)
	{
		for(int j=0;j<word[i].size();j++)
		{
			if(abs(word[i][j])>maxAmp)
			{
				maxAmp = abs(word[i][j]);
				windowIndex = i;
			}
		}
	}
	return windowIndex;
}
int _tmain(int argc, _TCHAR* argv[])
{
	char vowels[] = {'a','e','i','o','u'};
	//char vowels[] = {'a'};
	int vowelCount = sizeof(vowels)/sizeof(vowels[0]);

	//for(int x=0;x<vowelCount;x++)
	//{
	//	vector<vector<long double>> cepsTableForEachVowel[COUNT_STEADY_WINDOW];
	//	for(int y=1;y<=10;y++)
	//	{
	//		double dcShift = CalculateDCShift(vowels[x],y,PRE_PATH_INPUT);

	//		pair<double,double> p = CalculateMaxAmpandSilenceSTE(dcShift,vowels[x],y,PRE_PATH_INPUT);
	//
	//		double steSilence = p.first;
	//		double thresholdSound = steSilence*THRESH_FACTOR;

	//		double maxAmp = p.second;
	//		double normalFactor = MAX_AMP/maxAmp;

	//		//cout << steSilence << " " << thresholdSound << " " << maxAmp << " " << normalFactor << endl;
	//		string str = to_string(long long(y));
	//		string inFileName = PRE_PATH_INPUT;

	//		string suffix = "";
	//		suffix = suffix + vowels[x] + "_" + str + ".txt";

	//		inFileName = inFileName + suffix;
	//		//cout << inFileName << endl;
	//		ifstream infile;							
	//		infile.open(inFileName);
	//		/*ifstream infile;
	//		infile.open(SAMPLE_FILE);*/

	//		string outFileName = PRE_PATH_OUTPUT;
	//		outFileName = outFileName + suffix;
	//		ofstream outfile;
	//		outfile.open(outFileName);

	//		/*ofstream outfile;
	//		outfile.open("outputMultipleWords.txt");*/


	//		string compressFileName = PRE_PATH_COFF;
	//		compressFileName = compressFileName + suffix;
	//		ofstream compressFile;
	//		compressFile.open(compressFileName);

	//		/*ofstream compressFile;
	//		compressFile.open(COMPRESSED_FILE);*/

	//		if(infile.is_open())
	//		{
	//			int countIgnore =0;
	//			while(!infile.eof() && countIgnore<IGNORE_SAMPLES)
	//			{
	//				int sampleValue;
	//				infile >> sampleValue;
	//				countIgnore++;
	//			}

	//			vector<vector<double>> word;							//word variable store the sample of word
	//			int wordStarted=0;
	//			int wordCount = 0;
	//			int windowNumber = 0;
	//			while(!infile.eof())
	//			{
	//				double meanSqWindow=0;
	//				vector<double> windowSamples;

	//				while(windowSamples.size()<WINDOW_SIZE && !infile.eof())
	//				{
	//					double sampleValue;
	//					infile >> sampleValue;

	//					sampleValue = ShiftandNormalize(sampleValue,dcShift,normalFactor);

	//					windowSamples.push_back(sampleValue);
	//					meanSqWindow += (sampleValue*sampleValue)/double(WINDOW_SIZE);
	//				}

	//				if(windowSamples.size() < WINDOW_SIZE)
	//				{
	//					if(wordStarted==1)
	//					{
	//						if(word.size()>10)
	//							{
	//								int steadyIndex = findSteadyIndex(word);
	//								compressFile << wordCount << "th word windows are " << word.size() << " and Coefficients are:" << endl;
	//								if(steadyIndex+2<word.size() && steadyIndex-2>=0)
	//								{
	//									for(int i=steadyIndex-2;i<=steadyIndex+2;i++)
	//									{
	//										compressFile << i+1 << "th window coefficients are:" << endl;
	//						
	//										pair<vector<double>,vector<double>> windowCoefficients = CalculateandFindCoefficients(word[i]);

	//										vector<double> R = windowCoefficients.first;
	//										compressFile << "R's are:" << endl;
	//										for(int j=0;j<13;j++)
	//										{
	//											compressFile << R[j] << " ";
	//										}
	//										compressFile << endl;
	//										//R.clear();

	//										vector<double> alpha = windowCoefficients.second;
	//										if(alpha.size()==0)
	//										{
	//											compressFile << "Alpha's can't be calculated as R[0] is 0" << endl;
	//										}
	//										else
	//										{
	//											compressFile << "Alpha's are:" << endl;
	//											for(int j=1;j<alpha.size();j++)
	//											{
	//												compressFile << " " << alpha[j];
	//											}
	//											compressFile << endl;
	//										}

	//										vector<long double> ceps;
	//										cepstralcoefficients(R,alpha,ceps);

	//										compressFile << "Ceps's are:" << endl;
	//										for(int j=1;j<ceps.size();j++)
	//										{
	//											compressFile << " " << ceps[j];
	//										}
	//										compressFile << endl;
	//										compressFile << endl;
	//									}
	//								}
	//								else
	//								{
	//									compressFile << "Steady Windows out of range" << endl;
	//								}
	//							}
	//							else
	//							{
	//								outfile << "word was having less than 10 windows" << endl;
	//							}
	//						outfile << "word ended last word incompletely recorded" << endl;
	//						word.clear();
	//						wordStarted = 0;
	//					}
	//					break;
	//				}

	//				else
	//				{
	//					if(wordStarted==0)
	//					{
	//						if(meanSqWindow > thresholdSound)            // word started
	//						{
	//							wordStarted = 1;
	//							wordCount++;
	//							word.push_back(windowSamples);
	//							outfile<< endl;
	//							outfile<< "word started" << endl;
	//						}
	//					}
	//					else
	//					{
	//						if(meanSqWindow < thresholdSound)
	//						{
	//							if(word.size()>10)
	//							{
	//								int steadyIndex = findSteadyIndex(word);
	//								compressFile << wordCount << "th word windows are " << word.size() << " and Coefficients are:" << endl;
	//								if(steadyIndex+2<word.size() && steadyIndex-2>=0)
	//								{
	//									int windowNo = 0;
	//									for(int i=steadyIndex-2;i<=steadyIndex+2;i++)
	//									{
	//										compressFile << i+1 << "th window coefficients are:" << endl;
	//						
	//										pair<vector<double>,vector<double>> windowCoefficients = CalculateandFindCoefficients(word[i]);

	//										vector<double> R = windowCoefficients.first;
	//										compressFile << "R's are:" << endl;
	//										for(int j=0;j<13;j++)
	//										{
	//											compressFile << R[j] << " ";
	//										}
	//										compressFile << endl;
	//										//R.clear();

	//										vector<double> alpha = windowCoefficients.second;
	//										if(alpha.size()==0)
	//										{
	//											compressFile << "Alpha's can't be calculated as R[0] is 0" << endl;
	//										}
	//										else
	//										{
	//											compressFile << "Alpha's are:" << endl;
	//											for(int j=1;j<alpha.size();j++)
	//											{
	//												compressFile << " " << alpha[j];
	//											}
	//											compressFile << endl;
	//										}

	//										vector<long double> ceps;
	//										cepstralcoefficients(R,alpha,ceps);

	//										compressFile << "Ceps's are:" << endl;
	//										for(int j=1;j<ceps.size();j++)
	//										{
	//											compressFile << " " << ceps[j];
	//										}
	//										compressFile << endl;

	//										vector<long double> raisedCeps(Q+1,0);

	//										vector<long double> weights = RaisedSineWeights();

	//										raisedCeps[0] = ceps[0]*weights[0]; 

	//										compressFile << "Raised Ceps's are:" << endl;
	//										for(int j=1;j<ceps.size();j++)
	//										{
	//											raisedCeps[j] = ceps[j]*weights[j];
	//											compressFile << " " << raisedCeps[j];
	//										}
	//										compressFile << endl;

	//										cepsTableForEachVowel[windowNo].push_back(raisedCeps);

	//										windowNo++;

	//										compressFile << endl;
	//									}

	//								}
	//								else
	//								{
	//									compressFile << "Steady Windows out of range" << endl;
	//								}
	//							}
	//							else
	//							{
	//								outfile << "word was having less than 10 windows" << endl;
	//							}
	//							outfile<< "word ended" << endl;
	//							outfile<< endl;
	//							word.clear();
	//							wordStarted = 0;
	//						}
	//						else
	//						{
	//							word.push_back(windowSamples);
	//						}
	//					}
	//					windowNumber++;
	//					outfile << windowNumber << "th window:- " << meanSqWindow << endl;
	//					windowSamples.clear();
	//				}
	//		
	//			}
	//			infile.close();	
	//		}
	//		outfile.close();
	//		compressFile.close();
	//	}
	//	vector<long double> averageCepsforParticularVowel[COUNT_STEADY_WINDOW];
	//	cout << cepsTableForEachVowel[0].size() << endl;

	//	for(int k=0;k<COUNT_STEADY_WINDOW;k++)
	//	{
	//		for(int j=1;j<=Q;j++)
	//		{
	//			long double temp = 0;
	//			for(int i=0;i<10;i++)
	//			{
	//				temp += cepsTableForEachVowel[k][i][j]/10;
	//			}
	//			averageCepsforParticularVowel[k].push_back(temp);
	//		}
	//	}

	//	for(int k=0;k<COUNT_STEADY_WINDOW;k++)
	//	{
	//		for(int i=0;i<10;i++)
	//		{
	//			for(int j=1;j<=Q;j++)
	//			{
	//				cout << cepsTableForEachVowel[k][i][j] << " ";
	//			}
	//			cout << endl;
	//		}
	//		cout << endl;
	//	}
	//	//string averageCeps = AVERAGE_CEPS;

	//	string suffix = "";
	//	suffix = suffix + vowels[x] + ".txt";

	//	string averageCepsFileName = AVERAGE_CEPS;

	//	averageCepsFileName = averageCepsFileName + suffix;
	//	//cout << inFileName << endl;
	//	ofstream averageCepsEachVowel;							
	//	averageCepsEachVowel.open(averageCepsFileName);
	//	for(int i=0;i<COUNT_STEADY_WINDOW;i++)
	//	{
	//		for(int j=0;j<averageCepsforParticularVowel[i].size();j++)
	//		{
	//			averageCepsEachVowel << averageCepsforParticularVowel[i][j]<< " ";
	//		}
	//		averageCepsEachVowel << endl;
	//	}
	//	averageCepsEachVowel.close();
	//}

	ofstream result;
	result.open(RESULT_FILE);
	
	for(int x=0;x<vowelCount;x++)
	{
		for(int y=1;y<=10;y++)
		{
			vector<long double> cepsForEachTestVowel[COUNT_STEADY_WINDOW];

			double dcShift = CalculateDCShift(vowels[x],y,PRE_PATH_TEST);

			pair<double,double> p = CalculateMaxAmpandSilenceSTE(dcShift,vowels[x],y,PRE_PATH_TEST);
	
			double steSilence = p.first;
			double thresholdSound = steSilence*THRESH_FACTOR;

			double maxAmp = p.second;
			double normalFactor = MAX_AMP/maxAmp;

			//cout << steSilence << " " << thresholdSound << " " << maxAmp << " " << normalFactor << endl;
			string str = to_string(long long(y));
			string inFileName = PRE_PATH_TEST;

			string suffix = "";
			suffix = suffix + vowels[x] + "_" + str + ".txt";

			inFileName = inFileName + suffix;
			//cout << inFileName << endl;
			ifstream infile;							
			infile.open(inFileName);
			/*ifstream infile;
			infile.open(SAMPLE_FILE);*/

			string outFileName = PRE_PATH_TEST_OUTPUT;
			outFileName = outFileName + suffix;
			ofstream outfile;
			outfile.open(outFileName);

			/*ofstream outfile;
			outfile.open("outputMultipleWords.txt");*/


			string compressFileName = PRE_PATH_TEST_COFF;
			compressFileName = compressFileName + suffix;
			ofstream compressFile;
			compressFile.open(compressFileName);

			/*ofstream compressFile;
			compressFile.open(COMPRESSED_FILE);*/

			if(infile.is_open())
			{
				int countIgnore =0;
				while(!infile.eof() && countIgnore<IGNORE_SAMPLES)
				{
					int sampleValue;
					infile >> sampleValue;
					countIgnore++;
				}

				vector<vector<double>> word;							//word variable store the sample of word
				int wordStarted=0;
				int wordCount = 0;
				int windowNumber = 0;
				while(!infile.eof())
				{
					double meanSqWindow=0;
					vector<double> windowSamples;

					while(windowSamples.size()<WINDOW_SIZE && !infile.eof())
					{
						double sampleValue;
						infile >> sampleValue;

						sampleValue = ShiftandNormalize(sampleValue,dcShift,normalFactor);

						windowSamples.push_back(sampleValue);
						meanSqWindow += (sampleValue*sampleValue)/double(WINDOW_SIZE);
					}

					if(windowSamples.size() < WINDOW_SIZE)
					{
						if(wordStarted==1)
						{
							if(word.size()>10)
								{
									int steadyIndex = findSteadyIndex(word);
									compressFile << wordCount << "th word windows are " << word.size() << " and Coefficients are:" << endl;
									if(steadyIndex+2<word.size() && steadyIndex-2>=0)
									{
										for(int i=steadyIndex-2;i<=steadyIndex+2;i++)
										{
											compressFile << i+1 << "th window coefficients are:" << endl;
							
											pair<vector<double>,vector<double>> windowCoefficients = CalculateandFindCoefficients(word[i]);

											vector<double> R = windowCoefficients.first;
											compressFile << "R's are:" << endl;
											for(int j=0;j<13;j++)
											{
												compressFile << R[j] << " ";
											}
											compressFile << endl;
											//R.clear();

											vector<double> alpha = windowCoefficients.second;
											if(alpha.size()==0)
											{
												compressFile << "Alpha's can't be calculated as R[0] is 0" << endl;
											}
											else
											{
												compressFile << "Alpha's are:" << endl;
												for(int j=1;j<alpha.size();j++)
												{
													compressFile << " " << alpha[j];
												}
												compressFile << endl;
											}

											vector<long double> ceps;
											cepstralcoefficients(R,alpha,ceps);

											compressFile << "Ceps's are:" << endl;
											for(int j=1;j<ceps.size();j++)
											{
												compressFile << " " << ceps[j];
											}
											compressFile << endl;
											compressFile << endl;
										}
									}
									else
									{
										compressFile << "Steady Windows out of range" << endl;
									}
								}
								else
								{
									outfile << "word was having less than 10 windows" << endl;
								}
							outfile << "word ended last word incompletely recorded" << endl;
							word.clear();
							wordStarted = 0;
						}
						break;
					}

					else
					{
						if(wordStarted==0)
						{
							if(meanSqWindow > thresholdSound)            // word started
							{
								wordStarted = 1;
								wordCount++;
								word.push_back(windowSamples);
								outfile<< endl;
								outfile<< "word started" << endl;
							}
						}
						else
						{
							if(meanSqWindow < thresholdSound)
							{
								if(word.size()>10)
								{
									int steadyIndex = findSteadyIndex(word);
									compressFile << wordCount << "th word windows are " << word.size() << " and Coefficients are:" << endl;
									if(steadyIndex+2<word.size() && steadyIndex-2>=0)
									{
										int windowNo = 0;
										for(int i=steadyIndex-2;i<=steadyIndex+2;i++)
										{
											compressFile << i+1 << "th window coefficients are:" << endl;
							
											pair<vector<double>,vector<double>> windowCoefficients = CalculateandFindCoefficients(word[i]);

											vector<double> R = windowCoefficients.first;
											compressFile << "R's are:" << endl;
											for(int j=0;j<13;j++)
											{
												compressFile << R[j] << " ";
											}
											compressFile << endl;
											//R.clear();

											vector<double> alpha = windowCoefficients.second;
											if(alpha.size()==0)
											{
												compressFile << "Alpha's can't be calculated as R[0] is 0" << endl;
											}
											else
											{
												compressFile << "Alpha's are:" << endl;
												for(int j=1;j<alpha.size();j++)
												{
													compressFile << " " << alpha[j];
												}
												compressFile << endl;
											}

											vector<long double> ceps;
											cepstralcoefficients(R,alpha,ceps);

											compressFile << "Ceps's are:" << endl;
											for(int j=1;j<ceps.size();j++)
											{
												compressFile << " " << ceps[j];
											}
											compressFile << endl;

											vector<long double> raisedCeps(Q+1,0);

											vector<long double> weights = RaisedSineWeights();

											raisedCeps[0] = ceps[0]*weights[0]; 

											compressFile << "Raised Ceps's are:" << endl;
											for(int j=1;j<ceps.size();j++)
											{
												raisedCeps[j] = ceps[j]*weights[j];
												compressFile << " " << raisedCeps[j];
											}
											compressFile << endl;

											cepsForEachTestVowel[windowNo] = raisedCeps;

											windowNo++;

											compressFile << endl;
										}

									}
									else
									{
										compressFile << "Steady Windows out of range" << endl;
									}
								}
								else
								{
									outfile << "word was having less than 10 windows" << endl;
								}
								outfile<< "word ended" << endl;
								outfile<< endl;
								word.clear();
								wordStarted = 0;
							}
							else
							{
								word.push_back(windowSamples);
							}
						}
						windowNumber++;
						outfile << windowNumber << "th window:- " << meanSqWindow << endl;
						windowSamples.clear();
					}
			
				}
				infile.close();	
			}
			outfile.close();
			compressFile.close();

			/*for(int i=0;i<COUNT_STEADY_WINDOW;i++)
			{
				for(int j=1;j<=Q;j++)
				{
					cout << cepsForEachTestVowel[i][j] << " ";
				}
				cout << endl;
			}
			cout << endl;*/
			int flag =0;
			long double minDistance;
			char vow;
			for(int l=0;l<vowelCount;l++)
			{
				string suffix = "";
				suffix = suffix + vowels[l] + ".txt";

				string averageCepsFileName = AVERAGE_CEPS;

				averageCepsFileName = averageCepsFileName + suffix;
				//cout << averageCepsFileName << endl;
				ifstream averageCepsEachVowel;							
				averageCepsEachVowel.open(averageCepsFileName);

				long double totalDist =0;
				cout <<  "comparing " << vowels[l] << " with " << inFileName << endl;
				if(averageCepsEachVowel.is_open())
				{
					vector<long double> averageCepsforParticularVowel[COUNT_STEADY_WINDOW];
					for(int i=0;i<COUNT_STEADY_WINDOW;i++)
					{
						for(int j=0;j<Q;j++)
						{
							long double temp;
							averageCepsEachVowel >> temp;
							averageCepsforParticularVowel[i].push_back(temp);
						}
					}
					
					for(int i=0;i<COUNT_STEADY_WINDOW;i++)
					{
						long double currWindowDist = 0;
						for(int j=0;j<Q;j++)
						{
							currWindowDist += Tokhura_Weights[j]*((cepsForEachTestVowel[i][j+1] -averageCepsforParticularVowel[i][j])*(cepsForEachTestVowel[i][j+1] - averageCepsforParticularVowel[i][j]));
						}
						//cout << i+1 << "th window tokhura distance " << currWindowDist << endl; 
						totalDist += currWindowDist;
					}
					/*for(int i=0;i<COUNT_STEADY_WINDOW;i++)
					{
						for(int j=0;j<Q;j++)
						{
							cout << averageCepsforParticularVowel[i][j] << " ";
						}
						cout << endl;
					}*/
				}
				cout << "total distance from " <<  vowels[l] << " by summing all window is " << totalDist << endl;
				if(flag==0)
				{
					minDistance = totalDist;
					vow = vowels[l];
					flag = 1;
				}
				else
				{
					if(minDistance > totalDist)
					{
						cout << " entered min" << endl;
						minDistance = totalDist;
						vow = vowels[l];
					}
				}
				cout << minDistance << endl;
				//cout << vowelCount << endl;
				//cout << i << endl;
			}
			cout << inFileName << " vowel: " << vow << " minDistance: " << minDistance << endl;
			cout << endl;
			result << inFileName << " vowel: " << vow << " minDistance: " << minDistance << endl;
			//cepsForEachTestVowel[windowNo]

		}
		
		//vector<long double> averageCepsforParticularVowel[COUNT_STEADY_WINDOW];
		//cout << cepsTableForEachVowel[0].size() << endl;

		//for(int k=0;k<COUNT_STEADY_WINDOW;k++)
		//{
		//	for(int j=1;j<=Q;j++)
		//	{
		//		long double temp = 0;
		//		for(int i=0;i<10;i++)
		//		{
		//			temp += cepsTableForEachVowel[k][i][j]/10;
		//		}
		//		averageCepsforParticularVowel[k].push_back(temp);
		//	}
		//}

		//for(int k=0;k<COUNT_STEADY_WINDOW;k++)
		//{
		//	for(int i=0;i<10;i++)
		//	{
		//		for(int j=1;j<=Q;j++)
		//		{
		//			cout << cepsTableForEachVowel[k][i][j] << " ";
		//		}
		//		cout << endl;
		//	}
		//	cout << endl;
		//}
		////string averageCeps = AVERAGE_CEPS;

		//string suffix = "";
		//suffix = suffix + vowels[x] + ".txt";

		//string averageCepsFileName = AVERAGE_CEPS;

		//averageCepsFileName = averageCepsFileName + suffix;
		////cout << inFileName << endl;
		//ofstream averageCepsEachVowel;							
		//averageCepsEachVowel.open(averageCepsFileName);
		//for(int i=0;i<COUNT_STEADY_WINDOW;i++)
		//{
		//	for(int j=0;j<averageCepsforParticularVowel[i].size();j++)
		//	{
		//		averageCepsEachVowel << averageCepsforParticularVowel[i][j]<< " ";
		//	}
		//	averageCepsEachVowel << endl;
		//}
		//averageCepsEachVowel.close();
	}

	result.close();













	system("pause");
	return 0;
}


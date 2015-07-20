/*
 * loadDL_LZ.cpp
 *
 *  Created on: 25-03-2015
 *      Author: hector
 */

#include "LZDocList64.h"
#include <ConfigFile.h>

bool PRINT = false;			// true: print all details for console
bool TEST_IND = true;		// true: apply exhaustive test
bool RUN_EXP = false;
uint REPET = 5;				// number of repetitions for each experiment
uint MAX_M = 10;			// maximum length pattern value to run test

// Structure with all globals parameters program
typedef struct {
	string configFile;		// properties file

	uchar *seq;				// original sequence (1 byte for symbol)
	ulong n;				// Length of generalize Text = T1$T2$...TD$
	char cutDoc;			// symbol to separate documents
	bool lowerCPatt;		// 1: transform the patterns to lowercase

	LZDocList64 *index;

	ulong* EndDocs;			// this store the final phrase number for each document. This is no include in the final structure

	bool pattFromFile;		// 0: random patterns, 1: from file (as in todoCL)
	char dirStore[300];		// directory to save/load the data structure (files *.tk)
	char dirResult[300];	// directory to save tables

	ulong* patterns;
	uchar **patt;
	suint* lenPat;

	char fileTodoCL1W[300];
	char fileTodoCL2W[300];
} ParamProgram;

void runExperiments(ParamProgram *par);
void runExperimentsTodocl(ParamProgram *par);
void testDL(ParamProgram *par);
bool searchPattInFMI(ParamProgram *par, uchar* pat, uint m, ulong *occt1, uint *V1, uint *d1, ulong *occt2, uint *V2, uint *d2, ulong *occt3, uint *V3, uint *d3);
bool searchPattRevTrie(ParamProgram *par, uchar* pat, uint m, long int *x);

int main(int argc, char *argv[]) {
	ParamProgram *par = new ParamProgram();
	char fileName[300];

	if(argc != 2){
		cout << "ERRORR !! " << endl;
		cout << "loadDL_LZ.cpp's usage requires the Properties File as parameter !! " << endl;
		cout << "Example for the file 'config.txt' in the same directory: ./loadDL_LZ.cpp config.txt" << endl;
		exit(1);
	}
	par->configFile = string(argv[1]);
	cout << "Congif file: " << par->configFile << endl;
	ConfigFile cf(par->configFile);

	PRINT = cf.Value("GLOBALS","TRACE");
	TEST_IND = cf.Value("GLOBALS","TEST");
	REPET = cf.Value("GLOBALS","N_REP");
	RUN_EXP = cf.Value("GLOBALS","RUN_EXP");
	MAX_M = cf.Value("GLOBALS","MAX_M");

	strcpy(par->dirStore, ((string)(cf.Value("DL","dirStore"))).c_str());
	strcpy(par->dirResult, ((string)(cf.Value("DL","dirResult"))).c_str());
	par->pattFromFile = cf.Value("DL","pattFromFile");	// 0:random
	par->lowerCPatt = cf.Value("DL","lowercase");
	if(par->pattFromFile){
		strcpy(par->fileTodoCL1W, ((string)(cf.Value("DL","fileTodoCL1W"))).c_str());
		strcpy(par->fileTodoCL2W, ((string)(cf.Value("DL","fileTodoCL2W"))).c_str());
	}

	cout << "loadAppTopkLZ parameters..." << endl;
	cout << "dirStore: " << par->dirStore << endl;
	cout << "dirResult: " << par->dirResult << endl;
	cout << "patterns type: " << par->pattFromFile << endl;
	cout << "lowercase patterns: " << par->lowerCPatt << endl;
	if(par->pattFromFile){
		cout << "fileTodoCL1W : " << par->fileTodoCL1W << endl;
		cout << "fileTodoCL2W : " << par->fileTodoCL2W << endl;
	}

	par->index = new LZDocList64(par->dirStore);
	par->n = par->index->n;
	par->cutDoc = par->index->cutDoc;

	cout << "____________________________________________________" << endl;
	cout << "***  SizeUpRange = " << par->index->sizeUpRange << " bytes = " << (float)par->index->sizeUpRange*8.0/(float)par->n << " bpc" << endl;
	cout << "***  Index size " << par->index->sizeDS << " bytes = " << (float)par->index->sizeDS*8.0/(float)par->n << " bpc" << endl;
	cout << "====================================================" << endl;

	if (TEST_IND){
		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "sequence.test");
		ifstream is(fileName, ios::binary);
		par->seq = new uchar[par->n];
		is.read((char*)par->seq, par->n*sizeof(uchar));
		is.close();

		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "sep_rrr.test");
		load_from_file(par->index->sep_rrr, fileName);

		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "sep_rank.test");
		load_from_file(par->index->sep_rank, fileName);
		util::init_support(par->index->sep_rank, &par->index->sep_rrr);

		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "EndDocs.test");
		ifstream is2(fileName, ios::binary);
		par->EndDocs = new ulong[par->index->nDocs];
		is2.read((char*)par->EndDocs, par->index->nDocs*sizeof(ulong));
		is2.close();
		par->index->EndDocs = par->EndDocs;

		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "fmi.test");
		load_from_file(par->index->fmi, fileName);
		cout << "Load FMI complete, FMI requires " << size_in_mega_bytes(par->index->fmi) << " MiB." << endl;

		cout << "Test Index..." << endl;
		testDL(par);
		cout << "Test Index OK !!" << endl;
	}

	if (RUN_EXP){
		if(par->pattFromFile){
			runExperimentsTodocl(par);
		}else{
			cout << "Experiments..." << endl;
			runExperiments(par);
			cout << "Experiments OK !!" << endl;
		}
	}

	cout << "$$$$$$$$$$$$$$$$$$$$$" << endl;
	return 0;
}

void runExperimentsOneTodocl(ParamProgram *par){
	uint m;
	uchar *patron;
	uint *occ = new uint[par->index->nDocs];
	double t1, t2, t3, avgTime1, avgTime2, avgTime3;
	uint i, k, nDoccTot, nDocs1, nDocs2, nDocs3;
	float avgnOcc1, avgnOcc2, avgnOcc3, avgM;
	char aFile[300];

	avgTime1 = avgTime2 = avgTime3 = avgM =0.0;
	avgnOcc1 = avgnOcc2 = avgnOcc3 = 0;
	for (i=0; i<par->index->nDocs; i++)
		par->index->V[i] = 0;

	for (k=0; k<REPET; k++){
		nDocs1 = nDocs2 = nDocs3 =0;
		patron = par->patt[k];
		m = par->lenPat[k];
		avgM += m;
		t1 = getTime_ms();
		par->index->documentListingTypeOne(patron, m, occ, &nDocs1);
		t1 = getTime_ms() - t1;
		avgTime1 += t1/(double)REPET;
		nDoccTot = nDocs1;
		avgnOcc1 += nDocs1;

		par->index->documentListingTypeTwoThree(patron, m, occ, &nDoccTot, &nDocs2, &nDocs3, &t2, &t3);
		avgTime2 += t2/(double)REPET;
		avgTime3 += t3/(double)REPET;
		avgnOcc2 += nDocs2;
		avgnOcc3 += nDocs3;

		// It cleans the V[] vector or unmarks the documents reported !!
		for (i=0; i<nDoccTot; i++)
			par->index->V[occ[i]] = 0;
	}
	avgnOcc1 /= (float)REPET;
	avgnOcc2 /= (float)REPET;
	avgnOcc3 /= (float)REPET;
	avgM /= (float)REPET;
	cout << " TODOCL... " << endl;
	cout << "Average CPU time for execution, Search type 1 : " << avgTime1*1000.0 << " Microseconds" << endl;
	cout << "Average CPU time for execution, Search type 2 : " << avgTime2*1000.0 << " Microseconds" << endl;
	cout << "Average CPU time for execution, Search type 3 : " << avgTime3*1000.0 << " Microseconds" << endl;
	cout << "Average CPU total time for execution : " << (avgTime1+avgTime2+avgTime3)*1000.0 << " Microseconds" << endl;
	cout << "Average nDocs type 1 found : " << avgnOcc1 << endl;
	cout << "Average nDocs type 2 found : " << avgnOcc2 << endl;
	cout << "Average nDocs type 3 found : " << avgnOcc3 << endl;
	cout << "Size : " << par->index->sizeDS*8.0/(float)par->n << endl;
	cout << "Average length pattern : " << avgM << endl;
	cout << "____________________________________________________" << endl;

	strcpy(aFile, par->dirResult);
	strcat(aFile, "TodoclLZDL.txt");
	FILE *fp = fopen(aFile, "a+" );
	//[avgM] [size without Range] [total size] [n_occ1] [t_occ1] [n_occ2] [t_occ2] [n_occ3] [t_occ3] [n_occTotal] [t_occTotal]
	fprintf(fp, "%f %f %f %f %G %f %G %f %G %f %G\n", avgM, par->index->sizeUpRange*8.0/(float)par->n, par->index->sizeDS*8.0/(float)par->n, avgnOcc1, avgTime1*1000.0, avgnOcc2, avgTime2*1000.0, avgnOcc3, avgTime3*1000.0, avgnOcc1+avgnOcc2+avgnOcc3, (avgTime1+avgTime2+avgTime3)*1000.0);
	fclose(fp);
}

void runExperimentsOne(ParamProgram *par, uint m){
	uint *occ = new uint[par->index->nDocs];
	double t1, t2, t3, avgTime1, avgTime2, avgTime3;
	uint i, k, nDoccTot, nDocs1, nDocs2, nDocs3;
	float avgnOcc1, avgnOcc2, avgnOcc3;
	char aFile[300];

	cout << "____________________________________________________" << endl;
	cout << "Start DL for patterns of length m = " << m << endl;
	avgTime1 = avgTime2 = avgTime3 = 0.0;

	avgnOcc1 = avgnOcc2 = avgnOcc3 = 0;
	for (i=0; i<par->index->nDocs; i++)
		par->index->V[i] = 0;

	for (k=0; k<REPET; k++){
		nDocs1 = nDocs2 = nDocs3 =0;
		t1 = getTime_ms();
		par->index->documentListingTypeOne(par->seq+par->patterns[k], m, occ, &nDocs1);
		t1 = getTime_ms() - t1;
		avgTime1 += t1/(double)REPET;
		nDoccTot = nDocs1;
		avgnOcc1 += nDocs1;

		par->index->documentListingTypeTwoThree(par->seq+par->patterns[k], m, occ, &nDoccTot, &nDocs2, &nDocs3, &t2, &t3);
		avgTime2 += t2/(double)REPET;
		avgTime3 += t3/(double)REPET;
		avgnOcc2 += nDocs2;
		avgnOcc3 += nDocs3;

		// It cleans the V[] vector or unmarks the documents reported !!
		for (i=0; i<nDoccTot; i++)
			par->index->V[occ[i]] = 0;
	}
	avgnOcc1 /= (float)REPET;
	avgnOcc2 /= (float)REPET;
	avgnOcc3 /= (float)REPET;
	cout << "Average CPU time for execution, Search type 1 : " << avgTime1*1000.0 << " Microseconds" << endl;
	cout << "Average CPU time for execution, Search type 2 : " << avgTime2*1000.0 << " Microseconds" << endl;
	cout << "Average CPU time for execution, Search type 3 : " << avgTime3*1000.0 << " Microseconds" << endl;
	cout << "Average CPU total time for execution : " << (avgTime1+avgTime2+avgTime3)*1000.0 << " Microseconds" << endl;
	cout << "Average nDocs type 1 found : " << avgnOcc1 << endl;
	cout << "Average nDocs type 2 found : " << avgnOcc2 << endl;
	cout << "Average nDocs type 3 found : " << avgnOcc3 << endl;
	cout << "Size : " << par->index->sizeDS*8.0/(float)par->n << endl;
	cout << "____________________________________________________" << endl;

	strcpy(aFile, par->dirResult);
	strcat(aFile, "LZDL.txt");
	FILE *fp = fopen(aFile, "a+" );
	// [m] [size without Range] [total size] [n_occ1] [t_occ1] [n_occ2] [t_occ2] [n_occ3] [t_occ3] [n_occTotal] [t_occTotal]
	fprintf(fp, "%d %f %f %f %G %f %G %f %G %f %G\n", m, par->index->sizeUpRange*8.0/(float)par->n, par->index->sizeDS*8.0/(float)par->n, avgnOcc1, avgTime1*1000.0, avgnOcc2, avgTime2*1000.0, avgnOcc3, avgTime3*1000.0, avgnOcc1+avgnOcc2+avgnOcc3, (avgTime1+avgTime2+avgTime3)*1000.0);
	fclose(fp);
}

void runExperimentsTwo(ParamProgram *par, uint m){
	double t, avgTime;
	uint i, k, lev, nDocs;
	float avgOcc;
	uint *occ = new uint[par->index->nDocs];
	char aFile[300], str[30];
	ulong sizeRMQ_Range = 0;

	strcpy(aFile, par->dirResult);
	sprintf(str, "LZDL_m%d_RangeLevels", m);
	strcat(aFile, str);
	FILE *fp = fopen(aFile, "a+" );

	cout << "____________________________________________________" << endl;
	cout << "Start DL Range-Levels for patterns of length m = " << m << endl;

	for (i=0; i<par->index->nDocs; i++)
		par->index->V[i] = 0;
	for (lev=0; lev<par->index->h; lev++)
		sizeRMQ_Range += par->index->Range[lev].rmqB->getSize();

	cout << " ********* sizeDS = " << par->index->sizeDS << endl;
	cout << " ********* sizeDS - sizeRMQ_Range = " << par->index->sizeDS - sizeRMQ_Range << endl;

	for (lev=0; lev<=par->index->h; lev++){
		if(lev%2 == 0){
			avgOcc = 0;
			avgTime = 0.0;
			for (k=0; k<REPET; k++){
				nDocs = 0;
				if (lev == par->index->h)
					par->index->levRMQ = par->index->h-1;
				else
					par->index->levRMQ = lev;
				t = getTime_ms();
				par->index->documentListing(par->seq+par->patterns[k], m, occ, &nDocs);
				t = getTime_ms() - t;
				avgTime += t/(double)REPET;
				avgOcc += nDocs;

				// It cleans the V[] vector (it unmarks the documents reported)
				for (i=0; i<nDocs; i++)
					par->index->V[occ[i]] = 0;
			}
			avgOcc /= (float)REPET;
			cout << "Summary with " << lev << " levels of Range with RMQs" << endl;
			cout << "Average CPU time for execution: " << avgTime*1000.0 << " Microseconds" << endl;
			cout << "Average nDocs found : " << avgOcc << endl;
			cout << "Size : " << (par->index->sizeDS-sizeRMQ_Range)*8.0/(float)par->n << endl;
			cout << "____________________________________________________" << endl;

			// [lev] [size] [time] [nDocs]
			fprintf(fp, "%d %f %G %f\n", lev, (par->index->sizeDS-sizeRMQ_Range)*8.0/(float)par->n, avgTime*1000.0, avgOcc);
		}
		if (lev < par->index->h)
			sizeRMQ_Range -= par->index->Range[lev].rmqB->getSize();
	}
	fclose(fp);
}

void createPatterns(ParamProgram *par, uint m){
	ulong i, j, k;
	bool eq;

	cout << "Creating " << REPET << " patterns of length m=" << m << endl;
	for(k=0; k<REPET; k++){
		eq = true;
		while(eq){
			i = (rand() % (par->n-(m+1)))+1;
			for(j=i; j<i+(m+1); j++){
				if (par->seq[j] <= par->cutDoc){
					i = (rand() % (par->n-(m+1)))+1;
					j = i-1;
				}
				if(j==0) break;
				else{
					if (j>i && par->seq[j-1] != par->seq[j])
						eq = false;
				}
			}
		}
		par->patterns[k] = i;
		//cout << "["<<k<<"] " << i << endl;
	}
}

void loadPatternsAll(ParamProgram *par, char *fileQueries){
	uint i, j, k, len, totLines;
	string line;
	ifstream input (fileQueries);
	assert(input.good()); 				// check open
	if (input.is_open()){
		totLines = len = 0;
		while (input.good()){
			totLines++;
			getline(input,line);
			k = line.size();
			if (k > 99 || k < 3) continue;
			len++;
		}
		input.close();
	} else{
		cout << "Unable to open file " << fileQueries << endl;
		exit(0);
	}

	cout << "Total lines: " << totLines << ", patterns selected: " << len << endl;
	REPET = len;
	par->patt = new uchar*[len];
	par->lenPat = new suint[len];
	ifstream input2(fileQueries);// open file
	if (input2.is_open()){
		i = 0;
		while (input2.good() && i < len){
			getline(input2,line);
			k = line.size();
			if (k > 99 || k < 3) continue;
			par->patt[i] = new uchar[k+1];
			par->lenPat[i] = k;
			for(j=0; j<k; j++){
				if (par->lowerCPatt)
					par->patt[i][j] = (uchar)tolower(line.at(j));
				else
					par->patt[i][j] = (uchar)(line.at(j));
			}
			par->patt[i][k] = '\0';
			i++;
			//cout << par->patt[i] << endl;
		}
		input2.close();
	}
	cout << i << " patterns loaded !!" << endl;
}

void runExperimentsTodocl(ParamProgram *par){
	cout << " Todocl for patterns with 1 word..." << endl;
	loadPatternsAll(par, par->fileTodoCL1W);
	runExperimentsOneTodocl(par);
	delete [] par->patt;
	delete [] par->lenPat;

	cout << " Todocl for patterns with 2 words..." << endl;
	loadPatternsAll(par, par->fileTodoCL2W);
	runExperimentsOneTodocl(par);
	delete [] par->patt;
	delete [] par->lenPat;
}

void runExperiments(ParamProgram *par){
	uint m;
	cout << "====================================================" << endl;

	par->patterns = new ulong[REPET];
	m = 6;
	createPatterns(par, m);
	cout << "Patterns created for m = " << m << endl;
	runExperimentsOne(par, m);
	runExperimentsTwo(par, m);

	m = 10;
	createPatterns(par, m);
	cout << "Patterns created for m = " << m << endl;
	runExperimentsOne(par, m);
	runExperimentsTwo(par, m);

	delete [] (par->patterns);
}

void testDL(ParamProgram *par){
	uint nDocs1, nDocs2, nDocs3, d1, d2, d3, m;
	ulong nOcc1, nOcc2, nOcc3, occt1, occt2, occt3;
	ulong k, l, t, pv, x;
	float av1, av2, av3;
	bool isInRev;
	uchar *pat = new uchar[MAX_M+1];

	uint *occ1 = new uint[par->index->nDocs];
	uint *occ2 = new uint[par->index->nDocs];
	uint *occ3 = new uint[par->index->nDocs];

	uint *V1 = new uint[par->index->nDocs];
	uint *V2 = new uint[par->index->nDocs];
	uint *V3 = new uint[par->index->nDocs];

	uint *auxV2 = new uint[par->index->nDocs];
	bool *auxV3 = new bool[par->index->nDocs];

	for (m=1; m<=MAX_M; m++){
		av1 = av2 = av3 = 0.0;

		for (t=0; t<REPET; t++){
			k = (rand() % (par->n-m-1));
			for(l=k; l<k+m; l++){
				if (par->seq[l] == par->cutDoc){
					k = (rand() % (par->n-m-1));
					l = k-1;
				}
			}
			for (l=0; l<m; l++)
				pat[l] = par->seq[k+l];
			pat[m] = '\0';
			// ERROR type2, [e="Twel], m=7, k=21197788, test=1, nDocs2 = 19 != d2 20
			/*pat[0] = 'e';
			pat[1] = '=';
			pat[2] = '"';
			pat[3] = 'T';
			pat[4] = 'w';
			pat[5] = 'e';
			pat[6] = 'l';
			PRINT = true;*/

			nOcc1 = nOcc2 = nOcc3 = 0;
			nDocs1 = nDocs2 = nDocs3 = 0;
			occt1 = occt2 = occt3 = 0;
			d1 = d2 = d3 = 0;
			for (l=0; l<par->index->nDocs; l++){
				par->index->V[l] = 0;
				V1[l] = V2[l] = V3[l] = auxV2[l] = auxV3[l] = 0;
				occ1[l] = occ2[l] = occ3[l] = 0;
			}

			// search in CSA...
			searchPattInFMI(par, pat, m, &occt1, V1, &d1, &occt2, V2, &d2, &occt3, V3, &d3);
			av1 += occt1;
			av2 += occt2;
			av3 += occt3;
			if (PRINT){
				cout << "pat=" << pat;
				cout << ", m=" << m<< ", k:" << k << ", test=" << t  << endl;
				cout << "occ type_1 = " << occt1 << ", nDocs1 = " << d1 << endl;
				cout << "occ type_2 = " << occt2 << ", nDocs2 = " << d2 << endl;
				cout << "occ type_3 = " << occt3 << ", nDocs3 = " << d3 << endl;
			}

			// TEST OCCURRENCES TYPE 1...
			if(true){
				x = 1;
				isInRev = par->index->searchPattern_Rev(pat, m, &x, &pv);
				if (isInRev){
					par->index->documentListingTypeOne_test(pat, m, occ1, &nOcc1, &nDocs1);
					/*cout << "OUTPUT ";
					for (l=0; l<nOcc1; l++)
						cout << occ1[l];
					cout << endl;*/
					if (nOcc1 != occt1){
						cout << "ERROR type1b, [" << pat << "], m=" << m << ", k=" << k << ", test=" << t << ", has " << nOcc1 << " occ type_1 != " << occt1 << endl;
						exit(1);
					}
					if (nDocs1 != d1){
						cout << "ERROR type1b, [" << pat << "], m=" << m << ", k=" << k << ", test=" << t << ", nDocs1 = " << nDocs1 << " != d1 " << d1 << endl;
						exit(1);
					}
				}
			}

			if (m>1){
				// TEST OCCURRENCES TYPE 2 and 3...
				par->index->mark = 0;
				par->index->documentListingTypeTwoThree_test(pat, m, occ2, &nOcc2, &nDocs2, occ3, &nOcc3, &nDocs3, auxV2, auxV3);

				if (nOcc2 != occt2){
					cout << "ERROR type2, [" << pat << "], m=" << m << ", k=" << k << ", l=" << l << ", test=" << t << ", has " << nOcc2 << " occ type_2 != " << occt2 << endl;
					exit(1);
				}
				if (nDocs2 != d2){
					cout << "ERROR type2, [" << pat << "], m=" << m << ", k=" << k << ", test=" << t << ", nDocs2 = " << nDocs2 << " != d2 " << d2 << endl;
					exit(1);
				}

				if (m>2){
					// TEST OCCURRENCES TYPE 3...
					if (nOcc3 != occt3){
						cout << "ERROR type3, [" << pat << "], m=" << m << ", k=" << k << ", l=" << l << ", test=" << t << ", has " << nOcc3 << " occ type_3 != " << occt3 << endl;
						exit(1);
					}
					if (nDocs3 != d3){
						cout << "ERROR type3, [" << pat << "], m=" << m << ", k=" << k << ", test=" << t << ", nDocs3 = " << nDocs3 << " != d3 " << d3 << endl;
						exit(1);
					}
				}
			}

			//cout << "OK" << endl;exit(0);

		}
		av1 /= (float)REPET;
		av2 /= (float)REPET;
		av3 /= (float)REPET;
		//cout << "  AVG occ1 = " << av1 << ", occ2 = " << av2 << ", occ3 = " << av3 << endl;
	}
}

bool searchPattInFMI(ParamProgram *par, uchar* pat, uint m, ulong *occt1, uint *V1, uint *d1, ulong *occt2, uint *V2, uint *d2, ulong *occt3, uint *V3, uint *d3){
	ulong o1, o2, o3, doc;
	uint64_t lb=0, rb=par->index->fmi.size()-1;
	uint64_t l, r;
	string query = string((char *)pat);

	//size_t occs = sdsl::count(par->index->fmi, query.begin(), query.begin()+m);
	backward_search(par->index->fmi, lb, rb, query.begin(), query.begin()+m, l, r);

	//cout << "Occ Interval = [" << l << ", " << r << "]" << endl;
	size_t occs = r-l+1;
	if (!occs)
		return false;

	auto locations = locate(par->index->fmi, query.begin(), query.begin()+m);
	//cout << "Total occurrences found with FMI : " << occs << endl;

	o1 = o2 = o3 = 0;
	//cout << "locations..." << endl;
	for(ulong i=0; i<occs; i++){
		//cout << locations[i] << endl;
		//cout << " : " << par->sep_rrr[locations[i]] << par->sep_rrr[locations[i]+1] << par->sep_rrr[locations[i]+2] << endl;
		//cout << " : " << par->index->seq[locations[i]] << par->index->seq[locations[i]+1] << par->index->seq[locations[i]+2] << endl;
		r = par->index->sep_rank.rank(locations[i]+m) - par->index->sep_rank.rank(locations[i]+1);
		ulong rank = par->index->sep_rank.rank(locations[i]+1);
		doc = par->index->searchDocument(rank);

		if (r==0){
			o1++;
			if (V1[doc]==0){
				V1[doc]=1;
				(*d1)++;
			}
			//cout << locations[i] << ",(t1)doc:" << doc << endl;
		}else{
			if (r==1){
				o2++;
				if (V2[doc]==0){
					V2[doc]=1;
					(*d2)++;
				}
				//cout << locations[i] << ",(t2)doc:" << doc << endl;
			}else{
				o3++;
				if (V3[doc]==0){
					V3[doc]=1;
					(*d3)++;
				}
				//cout << locations[i] << ",(t3)doc:" << doc << endl;
			}
		}
	}
	//cout << endl;
	*occt1 = o1;
	*occt2 = o2;
	*occt3 = o3;

	return true;
}

bool searchPattRevTrie(ParamProgram *par, uchar* pat, uint m, long int *x){
	LZ78Tries64::RevNode* p = par->index->tries->revTrie;
	LZ78Tries64::RevNode *node = p->fChildRev;
	m--;
	while(node){
		if (pat[m] == node->symbol){
			if (m == 0){
				*x =node->idNode;
				return true;
			}else{
				m--;
				if(node->fict == true && node->uPath == true){
					// check unary path...
					uint i;
					for (i=0; i<node->lenUPath && m>0; i++, m--){
						if (pat[m] != node->uSymbols[i])
							return false;
					}
					if (m == 0){
						if (i<node->lenUPath){
							if (pat[0] == node->uSymbols[i]){
								*x =node->idNode;
								return true;
							}else
								return false;
						}
					}else{
						if(i < node->lenUPath)
							return false;
					}
				}
				node = node->fChildRev;
			}
		}else{
			if (pat[m] > node->symbol)
				node = node->nextSibRev;
			else
				return false;
		}

	}
	return false;
}

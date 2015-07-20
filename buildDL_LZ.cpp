//============================================================================
// Name        : buildDL_LZ.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "LZDocList64.h"
#include <ConfigFile.h>

bool TRACE = false;			// true: print all details for console
bool TEST = false;			// true: apply exhaustive test
uint N_REP = 5;

// Structure with all globals parameters program
typedef struct {
	uchar *seq;				// original sequence (1 byte for symbol)
	ulong n;				// Length of generalize Text = T1$T2$...TD$
	char cutDoc;			// symbol to separate documents
	LZDocList64 *index;

	ulong* EndDocs;			// this store the final phrase number for each document. This is no include in the final structure

	string configFile;		// properties file
	char inputFile[200];	// list of files
	bool lowercase;			// 1: transform to lowercase
	bool filesInList;		// 1: list of files, 0: Unique file
	char boundSymbol;		// original symbol delimiter of documents when we read all documents in 1 file.
	char dirStore[200];		// directory to save/load the data structure (files *.tk)
	uint levRMQ_range;		// The first levRMQ levels of Range will store RMQ structure


	// The following data structure are only for test !!
	ulong* patterns;
} ParamProgram;

void testSearchPatt_DL_CSA_LZ(ParamProgram *par);
bool searchPatternInFMI(ParamProgram *par, uchar* pat, uint m, ulong *occt1, uint *V1, uint *d1, ulong *occt2, uint *V2, uint *d2, ulong *occt3, uint *V3, uint *d3);
bool searchPatternRevTrie(ParamProgram *par, uchar* pat, uint m, long int *x);

int main(int argc, char *argv[]) {
	ParamProgram *par = new ParamProgram();
	char fileName[300];

	if(argc != 2){
		cout << "ERRORR !! " << endl;
		cout << "buildDL_LZ's usage requires the Properties File as parameter !! " << endl;
		cout << "Example for the file 'config.txt' in the same directory: ./buildDL_LZ config.txt" << endl;
		exit(1);
	}
	par->configFile = string(argv[1]);
	cout << "Congif file: " << par->configFile << endl;
	ConfigFile cf(par->configFile);

	TRACE = cf.Value("GLOBALS","TRACE");
	TEST = cf.Value("GLOBALS","TEST");
	N_REP = cf.Value("GLOBALS","N_REP");

	strcpy(par->inputFile, ((string)(cf.Value("DL","inputFile"))).c_str());
	par->filesInList = cf.Value("DL","filesInList");
	par->boundSymbol = cf.Value("DL","boundSymbol");
	par->cutDoc = cf.Value("DL","cutDoc");
	par->lowercase = cf.Value("DL","lowercase");
	par->levRMQ_range = cf.Value("DL","levRMQ");
	strcpy(par->dirStore, ((string)(cf.Value("DL","dirStore"))).c_str());

	cout << "buildDL_LZ config parameters..." << endl;
	cout << "Input File            : " << par->inputFile << endl;
	cout << "File in list          : " << par->filesInList << endl;
	cout << "Boundary Symbol(code) : " << (int)par->boundSymbol << endl;
	cout << "Cut Doc. Symbol(code) : " << (int)par->cutDoc << endl;
	cout << "Lowercase             : " << par->lowercase << endl;
	cout << "Range levels with RMQ : " << par->levRMQ_range << endl;
	cout << "Store Folder          : " << par->dirStore << endl;

	LZDocList64::TRACE = TRACE;
	LZDocList64::TEST = TEST;
	par->index = new LZDocList64(par->levRMQ_range, par->inputFile, par->filesInList, par->cutDoc, par->lowercase, par->dirStore, par->boundSymbol, false);
	par->n = par->index->n;
	par->EndDocs = par->index->EndDocs;
	cout << "____________________________________________________" << endl;
	cout << "***  Index size " << par->index->sizeDS << " bytes = " << (float)par->index->sizeDS*8.0/(float)par->n << " bpc" << endl;
	cout << "***  SizeUpRange = " << par->index->sizeUpRange << " bytes = " << (float)par->index->sizeUpRange*8.0/(float)par->n << " bpc" << endl;
	cout << "====================================================" << endl;

	par->index->saveDS(true);

	{
		// create the FMI to execute test...
		cout << "____________________________________________________" << endl;
		cout << " Make the FMI ..." << endl;
		strcpy(fileName, "");
		strcpy(fileName, par->inputFile);
		cout << " Reading... " << fileName << endl;
		strcat(fileName, "_copy.txt");
		construct(par->index->fmi, fileName, 1); // generate index
		cout << " **  FMI size " << size_in_bytes(par->index->fmi) << " bytes = " << (float)size_in_bytes(par->index->fmi)/(float)par->n << "|T|" << endl;
		cout << " **  FMI length " << par->index->fmi.size() << endl;
		if (par->index->fmi.size() != par->n){
			cout << "ERROR. FMI length != n = " << par->n << endl;
			exit(1);
		}
		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "fmi.test");
		store_to_file(par->index->fmi, fileName);
	}

	if (TEST){
		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "EndDocs.test");
		ofstream os2 (fileName, ios::binary);
		os2.write((const char*)par->EndDocs, par->index->nDocs*sizeof(ulong));					// save LbRev[]
		if(TRACE) cout << " .- EndDocs[] " << par->index->nDocs*sizeof(ulong) << " Bytes" << endl;
		os2.close();

		par->index->sep_rrr = rrr_vector<127>(par->index->sepPhrase_b);
		par->index->sep_rank = rrr_vector<127>::rank_1_type(&par->index->sep_rrr);
		std::cout << " sep_rrr uses " << size_in_bytes(par->index->sep_rrr) << " Bytes" << endl;

		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "sep_rrr.test");
		store_to_file(par->index->sep_rrr, fileName);

		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "sep_rank.test");
		store_to_file(par->index->sep_rank, fileName);

		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "EndDocs.test");
		ofstream os3 (fileName, ios::binary);
		os3.write((const char*)par->EndDocs, par->index->nDocs*sizeof(ulong));					// save LbRev[]
		if(TRACE) cout << " .- EndDocs[] " << par->index->nDocs*sizeof(ulong) << " Bytes" << endl;
		os3.close();

		// load Sequence...
		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "sequence.test");
		ifstream is(fileName, ios::binary);
		par->seq = new uchar[par->n];
		is.read((char*)par->seq, par->n*sizeof(uchar));
		is.close();

		cout << "Running test searchPattern.." << endl;
		testSearchPatt_DL_CSA_LZ(par);
		cout << "Test searchPattern OK !!" << endl;
	}

	par->index->~LZDocList64();

	cout << "$$$$$$$$$$$$$$$$$$$$$" << endl;
	return 0;
}

void testSearchPatt_DL_CSA_LZ(ParamProgram *par){
	uint *occ1 = new uint[par->index->nDocs];
	uint *occ2 = new uint[par->index->nDocs];
	uint *occ3 = new uint[par->index->nDocs];
	uint nDocs1, nDocs2, nDocs3;
	ulong nOcc1, nOcc2, nOcc3;

	ulong k, l, t, pv, x;
	uint m, M=10;
	float av1, av2, av3;
	bool found, isInRev;
	uchar *pat = new uchar[M+1];
	long int nod1;

	uint *V1 = new uint[par->index->nDocs];
	uint *V2 = new uint[par->index->nDocs];
	uint *V3 = new uint[par->index->nDocs];

	uint *auxV2 = new uint[par->index->nDocs];
	bool *auxV3 = new bool[par->index->nDocs];

	uint d1, d2, d3;
	ulong occt1, occt2, occt3;

	if (par->index->tries->nPhra < N_REP){
		N_REP = par->index->tries->nPhra;
		M = 6;
	}

	for (m=1; m<=M; m++){
		av1 = av2 = av3 = 0.0;

		for (t=0; t<N_REP; t++){
			k = (rand() % (par->n-m-1));
			for(l=k; l<k+m; l++){
				if (par->seq[l] <= 31){
					k = (rand() % (par->n-m-1));
					l = k-1;
				}
			}

			for (l=0; l<m; l++)
				pat[l] = par->seq[k+l];
			pat[m] = '\0';

			nOcc1 = nOcc2 = nOcc3 = 0;
			nDocs1 = nDocs2 = nDocs3 = 0;
			occt1 = occt2 = occt3 = 0;
			d1 = d2 = d3 = 0;
			for (l=0; l<par->index->nDocs; l++){
				par->index->V[l] = 0;
				V1[l] = V2[l] = V3[l] = auxV2[l] = 0;
				auxV3[l] = 0;
				occ1[l] = occ2[l] = occ3[l] = 0;
			}

			// search in CSA...
			searchPatternInFMI(par, pat, m, &occt1, V1, &d1, &occt2, V2, &d2, &occt3, V3, &d3);
			av1 += occt1;
			av2 += occt2;
			av3 += occt3;
			if (TRACE){
				cout << "pat=" << pat;
				cout << ", m=" << m<< ", k:" << k << ", test=" << t  << endl;
				cout << "occ type_1 = " << occt1 << ", nDocs1 = " << d1 << endl;
				cout << "occ type_2 = " << occt2 << ", nDocs2 = " << d2 << endl;
				cout << "occ type_3 = " << occt3 << ", nDocs3 = " << d3 << endl;
			}

			// TEST OCCURRENCES TYPE 1...
			if(true){
				nod1 = x = 1;
				found = searchPatternRevTrie(par, pat, m, &nod1);
				isInRev = par->index->searchPattern_Rev(pat, m, &x, &pv);
				if (found != isInRev){
					if (found)
						cout << "ERROR type1a, [" << pat << "], m=" << m << ", k=" << k << ", test=" << t << ", is in the node " << nod1 << " and NOT FOUND in RMM representation!!" << endl;
					else
						cout << "ERROR type1a, [" << pat << "], m=" << m << ", k=" << k << ", test=" << t << ", NOT FOUND in RevTrie, but IS FOUND in RMM representation!!" << endl;
					exit(1);
				}else{
					if (isInRev){
						par->index->documentListingTypeOne_test(pat, m, occ1, &nOcc1, &nDocs1);

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
			}

			if (m>1){
				// TEST OCCURRENCES TYPE 2 and 3...
				par->index->mark = 1;
				par->index->documentListingTypeTwoThree_test(pat, m, occ2, &nOcc2, &nDocs2, occ3, &nOcc3, &nDocs3, auxV2, auxV3);

				if (nOcc2 != occt2){
					cout << "ERROR type2, [" << pat << "], m=" << m << ", k=" << k << ", l=" << l << ", test=" << t << ", has " << nOcc2 << " occ type_2 != " << occt2 << endl;
					exit(1);
				}
				if (nDocs2 != d2){
					cout << "ERROR documents type2, [" << pat << "], m=" << m << ", k=" << k << ", test=" << t << ", nDocs2 = " << nDocs2 << " != d2 " << d2 << endl;
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
		}
		av1 /= (float)N_REP;
		av2 /= (float)N_REP;
		av3 /= (float)N_REP;
		//cout << "  AVG occ1 = " << av1 << ", occ2 = " << av2 << ", occ3 = " << av3 << endl;
	}
}

void createPatterns(ParamProgram *par, uint m, ulong repeat){
	ulong i, j, k;
	par->patterns = new ulong[repeat];
	bool eq;

	//cout << "Creating patterns of length m=" << m << endl;
	for(k=0; k<repeat; k++){
		eq = true;
		while(eq){
			i = (rand() % (par->n-(m+1)))+1;
			for(j=i; j<i+(m+1); j++){
				if (par->seq[j] == par->cutDoc){
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
	//cout << "Patterns created !!" << endl;
}

bool searchPatternInFMI(ParamProgram *par, uchar* pat, uint m, ulong *occt1, uint *V1, uint *d1, ulong *occt2, uint *V2, uint *d2, ulong *occt3, uint *V3, uint *d3){
	ulong r, o1, o2, o3, doc;
	string query = string((char *)pat);
	size_t occs = sdsl::count(par->index->fmi, query.begin(), query.begin()+m);
	auto locations = locate(par->index->fmi, query.begin(), query.begin()+m);
	//cout << "Total occurrences found with FMI : " << occs << endl;

	o1 = o2 = o3 = 0;
	//cout << "locations..." << endl;
	for(ulong i=0; i<occs; i++){
		//cout << locations[i] << endl;
		//cout << " : " << par->index->sep_rrr[locations[i]] << par->index->sep_rrr[locations[i]+1] << par->index->sep_rrr[locations[i]+2] << endl;
		//cout << " : " << par->index->seq[locations[i]] << par->index->seq[locations[i]+1] << par->index->seq[locations[i]+2] << endl;
		r = par->index->sep_rank.rank(locations[i]+m) - par->index->sep_rank.rank(locations[i]);
		ulong rank = par->index->sep_rank.rank(locations[i]+1);
		doc = par->index->searchDocument(rank);

		if(par->index->sep_rrr[locations[i]]){
			if (r==1){
				o1++;
				if (V1[doc]==0){
					V1[doc]=1;
					(*d1)++;
				}
				//cout << locations[i] << ",(t1)doc:" << doc << endl;
			}else{
				if (r==2){
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
		}else{
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
	}
	//cout << endl;
	*occt1 = o1;
	*occt2 = o2;
	*occt3 = o3;

	return true;
}

bool searchPatternRevTrie(ParamProgram *par, uchar* pat, uint m, long int *x){
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

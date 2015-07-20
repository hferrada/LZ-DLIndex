/*
 * LZDocList64.h
 *
 *  Created on: 04-09-2014
 *      Author: hector
 */

#ifndef LZDOCLIST64_H_
#define LZDOCLIST64_H_

#include "LZ78Tries64.h"
#include <RangeMMTree64.h>
#include <RMQRMM64.h>

using namespace sdsl;
using namespace std;
using namespace drf64;

const uint W32 = 32;

class LZDocList64 {

	// structure for wt of Range, each level has a biststring B_rrr (of legth n'-1) and a rmq for its respective pointer of documents
	typedef struct {
		rrr_vector<127> B_rrr;
		rrr_vector<127>::rank_1_type B_rank1;
		rrr_vector<127>::rank_0_type B_rank0;
		RMQRMM64 *rmqB;
	} LevelRange;

	typedef struct {
		bit_vector B_b;
		long int *C;
		long int *Aux;
	} BaseRange;

private:
	uchar lgD;				// Binary logarithm (ceiling) of nDocs
	char boundSymbol;		// original symbol delimiter of documents when we read all documents in 1 file.
	ulong* charTf;			// store the tf for each symbol. This is no include in the final structure

	ulong nod;				// number of phrase nodes nod. Here nod = n'+1, n' of LZIndex
	uint lgNod;				// ceiling binary log of nod
	ulong nF;				// number of fictitious nodes in RevTrie, nF only include one fictitious node for each unary path of fictitious nodes
	ulong nEF;				// this is the difference between the total fictitious nodes in the original RevTrie an the RevTrie represented where we join the fictitious nodes in unary path
	uint nEFHead;			// number of fictitious nodes are head of a unary path
	ulong nodRev;			// nod + nF
	ulong *PLz;				// LZTrie's DFUDS sequence representation in 2*nod bits. The size is include in the RMM representation
	ulong nLz;				// length of PLz DFUDS sequences, these have exactly nLz = 2*nod nodes
	ulong nRev;				// length of PRev DFUDS sequences, these have exactly nRev = 2(nod+nF) nodes
	ulong *PRev;			// RevTrie's DFUDS sequence representation in 2(nod+nF) bits. The size is include in the RMM representation

	// Labels...
	uchar *LbRev;			// RevTrie's labels, for the n' phrase nodes, in DFUDS representation
	bit_vector fRev_b;		// this marks the fictitious nodes in LbRev
	rrr_vector<127> fRev_rrr;
	rrr_vector<127>::rank_1_type fRev_rank;
	uchar *LbRevF;			// RevTrie's labels for fictitious nodes with unary paths
	bit_vector fictU_b;		// this marks the fictitious nodes with unary path in FRev
	rrr_vector<127> fictU_rrr;
	rrr_vector<127>::rank_1_type fictU_rank;
	bit_vector listF_b;		// Delimiter for each unary sequence in LbRevF
	rrr_vector<127> listF_rrr;
	rrr_vector<127>::select_1_type lisF_sel;

	RangeMMTree64 *treeRev;	// range min-max tree representations for PRev

	ulong *DocLZ;			// Document array to phrases with length nod-1
	ulong *Node;			// This array give the preorder value in LZTrie, that is: Range[j] = i <-> the node with
							// preorder j in RevTrie corresponds to the node with preorder i in LZTrie

	ulong nDocRev;			// lenght of LDocRev_b
	bit_vector LDocRev_b;
	rrr_vector<127> LDocRev_rrr;
	rrr_vector<127>::rank_1_type LDocRev_rank1;
	rrr_vector<127>::select_1_type LDocRev_sel1;

	RMQRMM64 *rmqC;			// RMQ structure for Document array for all occ type 1, with length n

	void createStructure(uint *DocRev, ulong *ArrRange);
	void genSeqExtraFictRevTrie(LZ78Tries64::RevNode *node, ulong *posExtF, ulong *pfRev, ulong *posFU);
	void genDFUDSSeqRevTrie(LZ78Tries64::RevNode *node, ulong *pos, ulong *posRev);
	void genDFUDSSeqLZTrieNotLbl(LZ78Tries64::LZNode *node, ulong *pos);
	void setIntervalDocRev(LZ78Tries64::LZNode* nod, uint *DocRev, ulong *posDR);
	void genDocRevArrays(LZ78Tries64::RevNode* nod, ulong *ArrRange, ulong *preRev, uint *DocRev, ulong *posDR);
	void genDocArray(LZ78Tries64::LZNode* nod, ulong *DocLZ, ulong *ini);
	void createSupportRMQDocRev(uint *DocRev);
	void genNodoWTRange(BaseRange *netRange, uint level, ulong *ArrRange, ulong start, ulong length, ulong lowSymbol);
	ulong searchLZCoord(uint lev, ulong lcur, ulong rcur, ulong obj);
	void createRange(ulong *ArrRange);

	// read the list of input file "inputFile", count the documents of these files
	// and store it in a the sequence seq[]. in Seq[] the files are separated by "cutDoc" symbol.
	// If bitLowercase==true --> the symbols are converted to lowercase
	void readListFiles(char *inputFile, bool bitLowercase);

	// read the input file "inputFile", count the documents in this file which are separated by "boundS" symbol
	// and store it in a the sequence seq[]. in Seq[] the files are separated by "cutDoc" symbol
	// If bitLowercase==true --> the symbols are converted to lowercase
	void readUniqueFile(char *inputFile, bool bitLowercase);

	// makes priorities for hash tables
	void makePriorities(ulong* cPosTb);

public:
	static bool TRACE;		// true: print all details for console
	static bool TEST;		// true: print all details for console

	ulong sizeDS;			// size in bytes for all data structure (this include sizeRMQRange)
	ulong sizeUpRange;		// size in bytes for the data structure without include Range strunture
	uchar *seq;				// original sequence (1 byte for symbol)
	ulong n;				// Length of generalize Text = T1$T2$...TD$
	uint nDocs;				// Number of distinct documents D
	uint h;					// h = lg(n'-1) = lg(nod-2). Binary logarithm
	uint levRMQ;			// The first levRMQ levels of Range will store RMQ structure. If levRMQ > h --> we set levRMQ = h

	char cutDoc;			// symbol to separate documents
	char dirStore[200];		// directory to save/load the data structure (files *.tk)
	LZ78Tries64 *tries;		// structure to store the original LZTrie and RevTrie.

	LevelRange *Range;	    // Wavelet tree represented as an array of (h-1) Bitmaps of length (n'-1), with support rank/select. Range[1..h-1][i..n'-1]
							// it represents the pairs (j,i) <--> If node j in RevTrie has id=k, then node i in LZTrie has id=k+1
	ulong* EndDocs;			// this store the final phrase number for each document. This is used only in time construction and it will be no include in the final structure

	// this bit_vector is only for test search of patterns and not belong to the index!!
	bit_vector sepPhrase_b;	// This marks the beginning of each LZ-phrase
	rrr_vector<127> sep_rrr;
	rrr_vector<127>::rank_1_type sep_rank;
	rrr_vector<127>::select_1_type sep_sel;
	sdsl::csa_wt<wt_huff<rrr_vector<127> >, 4, 4> fmi;

	suint *V;					// In this we marking the reported documents.
	suint mark;					// V and mark are used to check reported documents

	LZDocList64(char dirSaveLoad[400]);
	LZDocList64(uint diffLev, char *inputFile, bool filesInList, char cutDocCode, bool bitLowercase, char dirSaveLoad[200], char boundS, bool showSize);
	virtual ~LZDocList64();

	uint searchDocument(ulong idNode);

	// [lcur, rcur] is the current bitstring in the level 'lev' (B[lcur, rcur]). [lobj, robj] is the target interval to search in LZTrie, and [lrmq, rrmq] is the interval to query by RMQ in RevTrie
	void searchIntervalInRange(uint lev, ulong lcur, ulong rcur, ulong lobj, ulong robj, ulong lrmq, ulong rrmq, uint* occ, uint* nDocs);

	// this stores in mat[i] the preorder of the node that matches with a phrase p_{i..m} (find all phrase i for this m)
	void searchPhase_Rev(uchar* pat, uint m, ulong *mat);

	// return true if the pattern is in the RevTrie and in '*pv' the respective preorder value in the tree,
	// in '*x' stores its DFUDS bit position. Return false if the pattern does not exist
	bool searchPattern_Rev(uchar* pat, uint m, ulong *x, ulong *pv);

	// check if obj (LZ) matches with some point in [lRev... rRev]
	bool searchRevRange(uint lev, ulong lcur, ulong rcur, ulong obj, ulong lRev, ulong rRev);

	void documentListingTypeOne(uchar* pat, uint m, uint *occ, uint *nDocs1);

	void reportDocOcc1(ulong l, ulong r, uint* occ, uint* nDocs);

	void reportDocOcc2(uint lev, ulong lcur, ulong rcur, ulong l, ulong r, uint* occ, uint* nDocs);

	// apply sequential search in the range [l,r]...
	void reportDocOcc2_levRange(uint lev, ulong lcur, ulong rcur, ulong l, ulong r, ulong lobj, ulong robj, uint* occ, uint* nDocs);

	void documentListingTypeTwoThree(uchar* pat, uint m, uint* occ, uint *nDocs, uint* nDocs2, uint* nDocs3, double* t2, double* t3);

	// DL for pattern p of length m
	void documentListing(uchar* pat, uint m, uint* occ, uint* nDocs);

	// save the Data Structure in folder 'dirStore'
	void saveDS(bool showSize);

	// load the Data Structure from the file 'fileName'
	void loadDS(bool showSize);


	// THE NEXT METHODS ARE ONLY FOR TEST:
	void documentListingTypeOne_test(uchar* pat, uint m, uint *occ, ulong *nOcc, uint *nDocs1);
	void documentListingTypeTwoThree_test(uchar* pat, uint m, uint* occ2, ulong* nOcc2, uint* nDocs2, uint* occ3, ulong* nOcc3, uint* nDocs3, uint* auxV2, bool* auxV3);
	void reportDocOcc2_levRangeOcc(uint lev, ulong lcur, ulong rcur, ulong l, ulong r, ulong lobj, ulong robj, uint* occ, uint* nDocs, uint* nDocs2);
	void reportDocOcc2_levRange_test(uint lev, ulong lcur, ulong rcur, ulong l, ulong r, ulong lobj, ulong robj, uint* occ2, ulong* nOcc2, uint* nDocs2, uint* auxV2);
	void reportDocOcc2Occ(uint lev, ulong lcur, ulong rcur, ulong l, ulong r, uint* occ, uint* nDocs, uint *nDocs2);
	void reportDocOcc2_test(uint lev, ulong lcur, ulong rcur, ulong l, ulong r, uint* occ2, ulong* nOcc2, uint* nDocs2, uint* V2_aux);
	void searchIntervalInRangeOcc(uint lev, ulong lcur, ulong rcur, ulong lobj, ulong robj, ulong lrmq, ulong rrmq, uint* occ, uint* nDocs, uint* nDocs2);
	void searchIntervalInRange_test(uint lev, ulong lcur, ulong rcur, ulong lobj, ulong robj, ulong lrmq, ulong rrmq, uint* occ, ulong* nOcc2, uint* nDocs, uint* auxV2);

};

#endif /* LZDOCLIST64_H_ */

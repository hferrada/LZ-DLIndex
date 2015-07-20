/*
 * LZDocList64.cpp
 *
 *  Created on: 04-09-2014
 *      Author: hector
 */

#include "LZDocList64.h"
bool LZDocList64::TRACE = false;
bool LZDocList64::TEST = false;
uchar auxPath[10000];	// this is used to store the labels of unary path at construction

LZDocList64::LZDocList64(char dirSaveLoad[400]) {
	bool showSize=true;
	strcpy(dirStore, dirSaveLoad);
	loadDS(showSize);
}

LZDocList64::LZDocList64(uint levRange, char *inputFile, bool filesInList, char cutDocCode, bool bitLowercase, char dirSaveLoad[300], char boundS, bool showSize) {
	ulong i;

	this->sizeUpRange = 0;

	// size for variables
	this->sizeDS = 7*sizeof(ulong) + 5*sizeof(uint) + sizeof(char) + sizeof(uchar);	// size for variables
	cout << " ** Size of Variables" << sizeDS << " = " << sizeDS << " bpc" << endl;
	this->levRMQ = levRange;
	this->cutDoc = cutDocCode;
	this->boundSymbol = boundS;
	strcpy(dirStore, dirSaveLoad);
	this->nDocs = this->n = this->nEFHead = 0;

	charTf = new ulong[LZ78Tries64::SIGMA];
	for(i=0; i<LZ78Tries64::SIGMA; i++)
		charTf[i] = 0;
	if (filesInList){
		readListFiles(inputFile, bitLowercase);
	}else
		readUniqueFile(inputFile, bitLowercase);

	ulong* cPosTb = new ulong[LZ78Tries64::SIGMA]; // hasTable position for each character
	makePriorities(cPosTb);

	this->EndDocs = new ulong[nDocs];
	this->lgD = ceilingLog64(nDocs, 2);
	if(!lgD) lgD=1;

	this->V = new suint[nDocs];
	sizeDS += nDocs*sizeof(suint);
	cout << " ** size of V[1..nDocs] " << nDocs*sizeof(suint) << " Bytes = " << nDocs*sizeof(suint)*8.0/(float)this->n << " bpc" << endl;

	//cout << "  Create Original LZTrie and RevTrie..." << endl;
	tries = new LZ78Tries64(nDocs, cutDoc, cPosTb);
	this->sepPhrase_b = bit_vector(n, 0);

	tries->createTriesFromText(seq, n, &Node, EndDocs, &sepPhrase_b);
	cout << "  LZTrie with " << tries->nPhra << " nodes (phrases)" << endl;
	sizeDS += tries->sizeNode;
	cout << " ** Size of Node array " << tries->sizeNode << " = " << tries->sizeNode*8.0/(float)n << " bpc" << endl;

	{	// save seq[] for test...
		char fileName[300];
		strcpy(fileName, "");
		strcpy(fileName, dirSaveLoad);
		strcat(fileName, "sequence.test");
		ofstream os (fileName, ios::binary);
		os.write((const char*)seq, n*sizeof(uchar));
		if(TRACE) cout << " .- Seq[] " << n*sizeof(uchar) << " Bytes" << endl;
		os.close();
		delete []seq;
	}

	tries->countFictNodRevTrie(tries->revTrie);	// this set tries->nFictU
	cout << "  RevTrie with " << tries->nPhra << " nodes, " << tries->nFictU << " nFictU nodes and " << tries->nExtraFict << " nExtraFict nodes" << endl;

	this->nod = tries->nPhra;
	this->lgNod = ceilingLog64(nod, 2);

	ulong preRev, posDR;
	preRev = posDR = 0;
	ulong *ArrRange = new ulong[tries->nPhra];	// It stores the points of range. Range[j] = i <--> If node j is in RevTrie has id=k, then node i in

	// LZTrie has id=k+1 (this is needed only in construction time)
	uint *DocRev = new uint[n];		// Document array for all occ ty.pe 1 (this is needed only in construction time), with length n
	this->LDocRev_b = bit_vector(n, 0);

	genDocRevArrays(tries->revTrie, ArrRange, &preRev, DocRev, &posDR);
	this->nDocRev = posDR;
	delete [] tries->IdDocPreLZ;

	LDocRev_rrr = rrr_vector<127>(LDocRev_b);
	sizeDS += size_in_bytes(LDocRev_rrr);
	cout << " ** size of LDocRev_rrr: " << size_in_bytes(LDocRev_rrr) << " = " << size_in_bytes(LDocRev_rrr)*8.0/(float)n << " bpc" << endl;
	LDocRev_rank1 = rrr_vector<127>::rank_1_type(&LDocRev_rrr);
	LDocRev_sel1 = rrr_vector<127>::select_1_type(&LDocRev_rrr);

	if (TRACE){
		cout << " LDocRev_rrr[0.." << LDocRev_rrr.size() << "]:" << endl;
		cout << LDocRev_b[0];
		for(i=1; i<LDocRev_rrr.size(); i++)
			cout << LDocRev_b[i];
		cout << endl;

		cout << endl << " LZTrie with n'= " << tries->nPhra << " nodes..." << endl;
		tries->listLZTrie();
		cout << "====================================================" << endl;
		cout << endl << " RevTrie with n'= " << tries->nPhra << " phrase nodes, nFictU = " << tries->nFictU << " fictitious nodes, and " << tries->nExtraFict << " extra fictitious nodes (nExtraFict)... " << endl;
		tries->listRevTrie();

		cout << " ArrRange[0.." << tries->nPhra-1 << "]:" << endl;
		cout << ArrRange[0] << " ";
		for(i=1; i<tries->nPhra; i++){
			if(i%10 == 0)
				cout << "- ";
			cout << ArrRange[i] << " ";
		}
		cout << endl;

		cout << " Node[0.." << tries->nPhra-1 << "]:" << endl;
		cout << getNum64(Node, 0, tries->lgPhr) << " ";
		uint j=tries->lgPhr;
		for(i=1; i<tries->nPhra; i++, j+=tries->lgPhr){
			if(i%10 == 0)
				cout << "- ";
			cout << getNum64(Node, j, tries->lgPhr) << " ";
		}
		cout << endl;

		ulong len = sepPhrase_b.bit_size();
		cout << " separator_b[0.." << len-1 << "]:" << endl;
		cout << sepPhrase_b[0];
		for(i=1; i<len; i++)
			cout << sepPhrase_b[i];
		cout << endl;
	}

	this->nF = tries->nFictU;
	this->nEF = tries->nExtraFict;
	this->nodRev = nod+nF;
	this->nRev = 2*nodRev;
	this->nLz = 2*nod;

	cout << "  DFUFS RevTrie seq with " << nRev << " bits" << endl;
	cout << "  DFUFS LZTrie seq with " << nLz << " bits" << endl;

	ulong size = nLz/W64;
	if (nLz%W64)
		size++;
	this->PLz = new ulong[size];	// The size is include in the RMM representation of LZTrie

	size = nRev/W64;
	if (nRev%W64)
		size++;
	this->PRev = new ulong[size]; 	// The size is include in the RMM representation of RevTrie

	this->LbRev = new uchar[nodRev];
	//cout << " length of LbRev  " << nodRev << endl;
	size = nodRev*sizeof(uchar);
	sizeDS += size;
	cout << " ** size of LbRev[]: " << size << " = " << size*8.0/(float)n << " bpc" << endl;
	//cout << " LbRev lenght " << nodRev << endl;

	this->LbRevF = new uchar[nEF];
	//cout << "## len of LbRevF  " << nEF << endl;
	size = nEF*sizeof(uchar);
	sizeDS += size;
	cout << " ** size of LbRevF[]: " << size << " = " << size*8.0/(float)n << " bpc" << endl;
	//cout << " LbRevF lenght " << nEF << endl;

	size = nod*lgD;
	if(size%W64)
		size = size/W64 + 1;
	else
		size /= W64;
	this->DocLZ = new ulong[size];
	size *= sizeof(ulong);
	sizeDS += size;
	cout << " ** size of DocLZ[1..n']: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;
	cout << "====================================================" << endl << endl;

	// set initial open parenthesis in DFUDS sequences
	setBit64(PLz, 0);
	setBit64(PRev, 0);

	createStructure(DocRev, ArrRange);
}

void LZDocList64::createStructure(uint *DocRev, ulong *ArrRange){
	ulong i, j, pos, posRev, posF, posExtF, posFU;
	cout << "  Create LZDL Structure from RevTrie..." << endl;

	//[1] Create the DFUDS sequences of parentheses PRev, the bitstring FRev to mark the fictitious nodes, and the label vectors LbRev and LbRevF...
	fRev_b = bit_vector(nRev/2, 0);
	//cout << "lenght of fRev_b (posF) = " << fRev_b.bit_size() << endl;
	fictU_b = bit_vector(nF, 0);
	//cout << "lenght of fictU_b (posFU) = " << fictU_b.bit_size() << endl;
	listF_b = bit_vector(nEF, 0);
	//cout << "lenght of listF_b (posExtF) = " << listF_b.bit_size() << endl;

	posF = posExtF = posFU = 0;
	genSeqExtraFictRevTrie(tries->revTrie, &posExtF, &posF, &posFU);
	//cout << "posExtF " << posExtF << ", posF " << posF << ", posFU " << posFU << endl;
	//cout << "genSeqExtraFictRevTrie OK" << endl;

	pos = 1; posRev = 0;
	genDFUDSSeqRevTrie(tries->revTrie, &pos, &posRev);
	//cout << "pos " << pos << ", posRev (LbRev) " << posRev << endl;
	pos = 1;
	genDFUDSSeqLZTrieNotLbl(tries->lzTrie, &pos);
	//cout << "pos (Plz) " << pos << " =? nLz = " << nLz <<endl;

	if(TEST){ // test for balanced sequences
		long int sum = 0;
		ulong r=0;

		for (; r<nLz; r++){
			if(readBit64(PLz, r))
				sum++;
			else sum--;
		}
		if(sum != 0){
			cout << " PLz DFUDS is not a balanced sequence of parentheses !! " << endl;
			exit(1);
		}
		//else cout << " Test for PLz. DFUDS is a well balanced sequence of parentheses !! " << endl;

		sum = r = 0;
		for (; r<nRev; r++){
			if(readBit64(PRev, r))
				sum++;
			else sum--;
		}
		if(sum != 0){
			cout << " PRev DFUDS is not a balanced sequence of parentheses !! " << endl;
			exit(1);
		}
		//else cout << " Test for PRev. DFUDS is a well balanced sequence of parentheses !! " << endl;
	}

	fRev_rrr = rrr_vector<127>(fRev_b);
	sizeDS += size_in_bytes(fRev_rrr);
	cout << " ** size of fRev_rrr: " << size_in_bytes(fRev_rrr) << " = " << size_in_bytes(fRev_rrr)*8.0/(float)n << " bpc" << endl;
	fRev_rank = rrr_vector<127>::rank_1_type(&fRev_rrr);

	fictU_rrr = rrr_vector<127>(fictU_b);
	sizeDS += size_in_bytes(fictU_rrr);
	cout << " ** size of fictU_rrr: " << size_in_bytes(fictU_rrr) << " = " << size_in_bytes(fictU_rrr)*8.0/(float)n << " bpc" << endl;
	fictU_rank = rrr_vector<127>::rank_1_type(&fictU_rrr);

	listF_rrr = rrr_vector<127>(listF_b);
	sizeDS += size_in_bytes(listF_rrr);
	cout << " ** size of listF_rrr: " << size_in_bytes(listF_rrr) << " = " << size_in_bytes(listF_rrr)*8.0/(float)n << " bpc" << endl;
	lisF_sel = rrr_vector<127>::select_1_type(&listF_rrr);

	if(TRACE){
		cout << "RevTrie with " << nod << " nodes, " << tries->nFictU << " nFictU nodes and " << tries->nExtraFict << " nExtraFict nodes" << endl;
		tries->listRevTrie();

		cout << " DFUDFS PLz[0.." << nLz-1 << "]:" << endl;
		cout << readBit64(PLz, 0);
		for(i=1; i<nLz; i++){
			if(i%10 == 0)
				cout << "-";
			cout << readBit64(PLz, i);
		}
		cout << endl;

		cout << " DFUDFS PRev[0.." << nRev-1 << "]:" << endl;
		cout << readBit64(PRev, 0);
		for(i=1; i<nRev; i++){
			if(i%10 == 0)
				cout << "-";
			cout << readBit64(PRev, i);
		}
		cout << endl;
		cout << " fRev[0.." << fRev_b.bit_size()-1 << "]:" << endl;
		cout << fRev_b[0];
		for(i=1; i<fRev_b.bit_size(); i++){
			if(i%10 == 0)
				cout << "-";
			cout << fRev_b[i];
		}
		cout << endl;
		cout << " fictU[0.." << fictU_b.bit_size()-1 << "]:" << endl;
		cout << fictU_b[0];
		for(i=1; i<fictU_b.bit_size(); i++){
			if(i%10 == 0)
				cout << "-";
			cout << fictU_b[i];
		}
		cout << endl;
		cout << " LbRev[0.." << nodRev-1 << "]:" << endl;
		for(i=0; i<nodRev; i++)
			cout << LbRev[i];

		cout << endl;
		if (nEF){
			cout << " LbRevF[0.." << nEF-1 << "]:" << endl;
			for(i=0; i<nEF; i++)
				cout << LbRevF[i];
			cout << endl;
			cout << " listF[0.." << listF_b.bit_size()-1 << "]:" << endl;
			for(i=0; i<listF_b.bit_size(); i++)
				cout << listF_b[i];
			cout << endl;
		}else
			cout << " LbRevF does not have any extra labels (unary path of fictitious node with length > 1)" << endl;

		cout << " EndDocs[0.." << nDocs-1 << "]:" << endl;
		for(i=0; i<nDocs; i++)
			cout << EndDocs[i] << " ";
		cout << endl;

		cout << " Node array[0.." << nod-1 << "]:" << endl;
		cout << getNum64(Node, 0, lgNod) << " ";
		for(i=1, j=lgNod; i<nod; i++, j+=lgNod){
			if(i%10 == 0)
				cout << "- ";
			cout << getNum64(Node, j, lgNod) << " ";
		}
		cout << endl;
	}
	cout << "  Create the representation with Range min-max of PLz and PRev DFUDS..." << endl;

	treeRev = new RangeMMTree64(PRev, nRev, LbRev, false);
	sizeDS += treeRev->sizeRMM;
	cout << " ** size of treeRev: " <<  treeRev->sizeRMM << " = " <<  treeRev->sizeRMM*8.0/(float)n << " bpc" << endl;

	//cout << "  Create DocLZ..." << endl;
	pos = 0;
	// making DocLZ in LZ78Tries
	genDocArray(tries->lzTrie, DocLZ, &pos);
	if(TRACE){
		cout << "  Document array[0.." << nod-1 << "]:" << endl;
		for(i=j=0; i<nod; i++, j+=lgD){
			if(i%10 == 0) cout << "-";
			cout << getNum64(DocLZ, j, lgD);
		}
		cout << endl;
	}

	cout << "====================================================" << endl;
	cout << "Creating RMQ structure on C array..." << endl;
	createSupportRMQDocRev(DocRev);
	delete [] DocRev;

	// it is the size of the structure without Range
	sizeUpRange = sizeDS;

	cout << "====================================================" << endl;
	cout << "Creating Range structure..." << endl;
	createRange(ArrRange);
	if(TRACE){
		cout << " ## Range[0.." << h-1 << "]x[0.." << nod-1 << "] ..." << endl;
		for(i=0; i<h; i++){
			cout << Range[i].B_rrr[0];
			for(j=1; j<nod; j++){
				if(j%10 == 0) cout << "-";
				cout << Range[i].B_rrr[j];
			}
			cout << endl;
		}
		cout << endl;
	}
}

// set extra fictitious labels...
void LZDocList64::genSeqExtraFictRevTrie(LZ78Tries64::RevNode *node, ulong *posExtF, ulong *pfRev, ulong *posFU){
	LZ78Tries64::RevNode *child;
	uint i;

	if (node->fict){
		fRev_b[*pfRev] = 1;
		(*posFU)++;
		if(node->nChildren == 1){
			child = node->fChildRev;
			if (child->fict){
				nEFHead++;
				fictU_b[*posFU - 1] = 1;
				listF_b[*posExtF] = 1;
				LbRevF[*posExtF] = child->symbol;
				auxPath[0] = child->symbol;
				uint unary = 1;
				(*posExtF)++;

				// search the last node in this unary path...
				while (child->nChildren == 1 && child->fChildRev->fict){
					LbRevF[*posExtF] = child->fChildRev->symbol;
					auxPath[unary] = child->fChildRev->symbol;
					unary++;
					(*posExtF)++;
					child = child->fChildRev;
				}
				node->nChildren = child->nChildren;
				node->fChildRev = child->fChildRev;
				node->uPath = true;
				node->lenUPath = unary;
				node->uSymbols = new uchar[unary];
				if (unary > 10000){
					cout << "ERROR unary = " << unary << " > 10000" << endl;
					exit(1);
				}

				for (i=0; i<unary; i++)
					node->uSymbols[i] = auxPath[i];
			}
		}
	}
	(*pfRev)++;
	if (*pfRev > fRev_b.bit_size()){
		cout << "ERRRRR: *pfRev = " << *pfRev << " > " << "fRev_b.bit_size() = " << fRev_b.bit_size() << endl;
		exit(1);
	}

	// recursive call for all children
	if (node->nChildren){
		child = node->fChildRev;
		for(i=0; i<node->nChildren; i++){
			genSeqExtraFictRevTrie(child, posExtF, pfRev, posFU);
			child = child->nextSibRev;
		}
	}
}

//	create 	PRev, LbRev and LbRevF structures
void LZDocList64::genDFUDSSeqRevTrie(LZ78Tries64::RevNode *node, ulong *pos, ulong *posRev){
	LZ78Tries64::RevNode *child;
	uint i;

	// set open parentheses
	for(i=0, child = node->fChildRev; i<node->nChildren; i++, (*pos)++, child=child->nextSibRev, (*posRev)++){
		LbRev[*posRev] = child->symbol;
		setBit64(PRev, *pos);
		//cout << "(";
	}
	if (child){
		cout << "ERRR.. nod " << node->idNode << " has " << node->nChildren << " children, but there is also: " << child->idNode << endl;
		exit(1);
	}

	// set close parenthesis
	cleanBit64(PRev, *pos);
	(*pos)++;
	//cout << ")";

	// recursive call for all children
	child = node->fChildRev;
	for(i=0; i<node->nChildren; i++){
		genDFUDSSeqRevTrie(child, pos, posRev);
		child = child->nextSibRev;
	}
	if (child){
		cout << "ERRR.. nod " << node->idNode << " has " << node->nChildren << " children, but there is also: " << child->idNode << endl;
		exit(1);
	}
}

//	create 	PLZ	:	LZTrie's DFUDS representation without Labels
void LZDocList64::genDFUDSSeqLZTrieNotLbl(LZ78Tries64::LZNode *node, ulong *pos){
	LZ78Tries64::LZNode *child;
	uint i;

	// set open parentheses
	for(i=0, child = node->fChildLZ; i<node->nChildren; i++, (*pos)++, child=child->nextSibLZ)
		setBit64(PLz, *pos);

	if (child){
		cout << "ERRR.. nod " << node->idNode << " has " << node->nChildren << " children, but there is also: " << child->idNode << endl;
		exit(1);
	}

	// set close parenthesis
	cleanBit64(PLz, *pos);
	(*pos)++;

	// recursive call for all children
	child = node->fChildLZ;
	for(i=0; i<node->nChildren; i++){
		genDFUDSSeqLZTrieNotLbl(child, pos);
		child = child->nextSibLZ;
	}
	if (child){
		cout << "ERRR.. nod " << node->idNode << " has " << node->nChildren << " children, but there is also: " << child->idNode << endl;
		exit(1);
	}
}

void LZDocList64::setIntervalDocRev(LZ78Tries64::LZNode* nod, uint *DocRev, ulong *posDR){
	LZ78Tries64::LZNode *p = nod->fChildLZ;

	if (nod){
		DocRev[*posDR] = searchDocument(nod->idNode);
		(*posDR)++;
		if(*posDR < this->n)
			LDocRev_b[*posDR] = false;
		for(uint i=0; i<nod->nChildren; i++){
			setIntervalDocRev(p, DocRev, posDR);
			p = p->nextSibLZ;
		}
	}
}

void LZDocList64::genDocRevArrays(LZ78Tries64::RevNode* nod, ulong *ArrRange, ulong *preRev, uint *DocRev, ulong *posDR){
	LZ78Tries64::RevNode* p = nod->fChildRev;
	if (nod->idNode == 0){
		setNum64(Node, 0, lgNod, 0);
		ArrRange[0] = tries->IdDocPreLZ[1];
		(*preRev)++;
	}else{
		if(nod->fict == false){
			setNum64(Node, (*preRev)*lgNod, lgNod, nod->nodeLZ->preorder);

			// store a new point for range structure
			if (nod->idNode == (int)(tries->nPhra-1))
				ArrRange[*preRev] = tries->IdDocPreLZ[0];
			else
				ArrRange[*preRev] = tries->IdDocPreLZ[nod->idNode+1];

			LDocRev_b[*posDR] = true;
			setIntervalDocRev(nod->nodeLZ, DocRev, posDR);
			(*preRev)++;
		}
	}

	for(uint i=0; i<nod->nChildren; i++){
		genDocRevArrays(p, ArrRange, preRev, DocRev, posDR);
		p = p->nextSibRev;
	}
}

// determine the doc Id for the phrase idNode
uint LZDocList64::searchDocument(ulong idNode){
	uint m, ini=0, end=nDocs-1;

	// binary search in the interval endDocs[0..nPhra-1]
	while(ini < end){
		m = ini+(end-ini)/2;
		if (idNode > EndDocs[m])
			ini = m+1;
		else{
			if(m){
				if (idNode > EndDocs[m-1])
					return m;
				else
					end = m-1;
			}else
				return 0;
		}
	}

	return ini;
}

// create the document array for LZTrie's phrases
void LZDocList64::genDocArray(LZ78Tries64::LZNode* nod, ulong *DocLZ, ulong *ini){
	LZ78Tries64::LZNode* p = nod->fChildLZ;
	ulong doc = searchDocument(nod->idNode);

	setNum64(DocLZ, *ini, lgD, doc);
	(*ini) += lgD;
	for(ulong i=0; i<nod->nChildren; i++){
		genDocArray(p, DocLZ, ini);
		p = p->nextSibLZ;
	}
}

void LZDocList64::createSupportRMQDocRev(uint *DocRev){
	long int *C = new long int[n];
	long int *Aux = new long int[nDocs];
	ulong i, j;

	for(i=0; i<nDocs; i++) Aux[i] = -1;
	for(i=0; i<n; i++){
		C[i] = Aux[DocRev[i]];
		Aux[DocRev[i]] = i;
	}
	delete [] Aux;
	if (TRACE){
		cout << "C = ";
		for(i=0; i<n; i++) cout << C[i] << " ";
		cout << endl;
	}

	rmqC = new RMQRMM64(C,n);
	cout << " ** size of RMQ_optimal on C[1.." << n << "]: " << rmqC->getSize() << " bytes = " << rmqC->getSize()*8.0/(float)n << " bpc" << endl;
	sizeDS += rmqC->getSize();

	if (TEST && false){
		cout << "Exhaustive test for RMQ on C..." << endl;
		uint rmq_min, min, k;

		for (i=0; i<n; i+=100000){
			for (j=i+1; j<n; j+=100000){
				rmq_min = rmqC->queryRMQ(i,j);

				min = i;
				for (k=i+1; k<=j; k++){
					if (C[k] < C[min])
						min = k;
				}

				if (rmq_min < i || rmq_min > j){
					cout << "ERROR... rmq_min = " << rmq_min << " out of range [" << i << ", " << j << " ]" << endl;
					exit(1);
				}else{
					if (C[rmq_min] != C[min]){
						cout << "ERROR... (" << i << ", " << j << " ) = " << rmq_min << " != " << min << endl;
						exit(1);
					}
					//else{ if (listingTrace) cout << "rmq(" << beg << "," << end << "):" << rmq_min << endl;}
				}
			}
		}

		cout << "RMQ on C test OK!!" << endl;
	}

	delete [] C;
}

void LZDocList64::genNodoWTRange(BaseRange *netRange, uint level, ulong *ArrRange, ulong start, ulong length, ulong lowSymbol){
	uint doc;
	ulong i, lenR = length>>1;
	ulong medSymbol = lowSymbol + lenR;
	if (length%2)
		medSymbol++;

	for(i=start; i<(start+length); i++){
		if(ArrRange[i] < medSymbol)
			netRange[level].B_b[i] = 0;
		else
			netRange[level].B_b[i] = 1;
	}

	if (level<levRMQ){
		for(i=start; i<(start+length); i++){
			doc = getNum64(DocLZ, (ArrRange[i])*lgD, lgD);
			netRange[level].C[i] = netRange[level].Aux[doc];
			netRange[level].Aux[doc] = i;
		}
	}

	if (false){
		cout << "C ponter..."<< endl;
		for(i=start; i<(start+length); i++)
			cout << netRange[level].C[i] << " ";
		cout << endl;
	}

	if (length>2){
		ulong l, r, *auxR = new ulong[lenR];
		for(i=l=start, r=0; i<(start+length); i++){
			if (netRange[level].B_b[i]){
				auxR[r] = ArrRange[i];
				r++;
			}else{
				ArrRange[l] = ArrRange[i];
				l++;
			}
		}
		for(i=0; i<lenR; i++)
			ArrRange[start+i+length-lenR] = auxR[i];

		delete [] auxR;

		genNodoWTRange(netRange, level+1, ArrRange, start, length-lenR, lowSymbol);
		if(lenR>1)
			genNodoWTRange(netRange, level+1, ArrRange, medSymbol, lenR, medSymbol);
	}
}

// search the coordinate y that match with obj. => (y, obj) is in Range
ulong LZDocList64::searchLZCoord(uint lev, ulong lcur, ulong rcur, ulong obj){
	if (lev == h-1){	// this is a leaf !!
		if (Range[lev].B_rrr[obj])
			return rcur;
		else
			return lcur;
	}

	ulong med = lcur + ((rcur-lcur)>>1);
	if (Range[lev].B_rrr[obj] == 0){
		// we go down to the left...
		obj = Range[lev].B_rank0.rank(obj+1);
		if (lcur)
			obj -= Range[lev].B_rank0.rank(lcur);
		return searchLZCoord(lev+1, lcur, med, lcur+obj-1);
	}else{
		//we go down to the right...
		obj = Range[lev].B_rank1.rank(obj+1);
		if (lcur)
			obj -= Range[lev].B_rank1.rank(lcur);
		return searchLZCoord(lev+1, med+1, rcur, med+obj);
	}

	return 0;
}

// create wavelet tree for range represented as an array of BitSequences...
void LZDocList64::createRange(ulong *ArrRange){
	ulong i;
	ulong *Rcopy=NULL;
	this->h = ceilingLog64(nod,2);

	if (h < levRMQ) levRMQ = h;
	//cout << "n' " << nod << ", h " << h << ", levRMQ " << levRMQ << endl;

	if (TEST){
		Rcopy = new ulong[nod];
		for(i=0; i<nod; i++)
			Rcopy[i] = ArrRange[i];
	}

	BaseRange *netRange = new BaseRange[h];
	for (i=0; i<h; i++){
		netRange[i].B_b = bit_vector(nod, 0);
		if (i<levRMQ){
			netRange[i].C = new long int[nod];
			netRange[i].Aux = new long int[nDocs];
			for(uint j=0; j<nDocs; j++){
				netRange[i].C[j] = j;
				netRange[i].Aux[j] = -1;
			}
		}
	}
	Range = new LevelRange[h];
	sizeDS += h*sizeof(LevelRange);
	cout << " ** size of Range[]'s pointers " << h*sizeof(LevelRange) << " = " << ( h*sizeof(LevelRange))*8.0/(float)this->n << " bpc" << endl;
	genNodoWTRange(netRange, 0, ArrRange, 0, nod, 0);

	if (TEST == false)
		delete []ArrRange;

	if(false){
		cout << " ## C for Range..." << endl;
		for(i=h-1; i<h; i++){
			cout << netRange[i].C[0] << " ";
			for(uint j=1; j<nod; j++){
				if (j%10 == 0) cout << "||";
				cout << netRange[i].C[j] << " ";

			}
			cout << endl;
		}
		cout << endl;
	}

	ulong sizeRMQsRange, sizeBitMapsRange;
	sizeRMQsRange = sizeBitMapsRange = 0;
	for (i=0; i<h; i++){
		Range[i].B_rrr = rrr_vector<127>(netRange[i].B_b);
		Range[i].B_rank1 = rrr_vector<127>::rank_1_type(&(Range[i].B_rrr));
		Range[i].B_rank0 = rrr_vector<127>::rank_0_type(&(Range[i].B_rrr));
		sizeBitMapsRange += size_in_bytes(Range[i].B_rrr);

		netRange[i].B_b.~int_vector();
		if (i<levRMQ){
			Range[i].rmqB = new RMQRMM64(netRange[i].C, nod);
			if (TEST && false){
				cout << "Test RMQ on netRange[" <<i<< "].C..." << endl;
				uint min, rmq_min, k, t, j;

				for (t=0; t<nod/2; t+=10000){
					for (j=nod/2+t+1; j<nod; j+=10000){
						rmq_min = Range[i].rmqB->queryRMQ(t,j);
						min = t;
						for (k=t+1; k<=j; k++){
							if (netRange[i].C[k] < netRange[i].C[min])
								min = k;
						}

						if (rmq_min < t || rmq_min > j){
							cout << "ERROR... rmq_min = " << rmq_min << " out of range [" << t << ", " << j << " ]" << endl;
							exit(1);
						}else{
							if (netRange[i].C[rmq_min] != netRange[i].C[min]){
								cout << "ERROR... (" << t << ", " << j << " ) = " << rmq_min << " != " << min << endl;
								exit(1);
							}
							//else{ if (listingTrace) cout << "rmq(" << beg << "," << end << "):" << rmq_min << endl;}
						}
					}
				}
				cout << " RMQ OK!!" << endl;
			}
			delete [] netRange[i].C;
			sizeRMQsRange += Range[i].rmqB->getSize();
		}
	}

	cout << " ** size of bitVectors (RRR) of Range: " << sizeBitMapsRange << " = " << sizeBitMapsRange*8.0/(float)this->n << " bpc" << endl;
	cout << " ** size of RMQs of Range: " << sizeRMQsRange << " = " << sizeRMQsRange*8.0/(float)this->n << " bpc" << endl;

	sizeDS += sizeBitMapsRange;
	sizeDS += sizeRMQsRange;

	if (TEST){
		//uint doc;
		ulong x;
		cout << " Test searchLZCoord..." << endl;
		//cout << " ArrRange Docs: " << endl;
		for (i=0; i<nod; i++){
			x = searchLZCoord(0, 0, nod-1, i);
			if (Rcopy[i] != x){
				cout << "ERROR!! points[" << i << "] = " << ArrRange[i] << " != searchLZCoord = " << x << endl;
				exit(0);
			}
			//doc = getNum64(DocLZ, (ArrRange[i])*lgD, lgD);
			//if (i%10 == 0) cout << " - ";
			//cout << doc << " ";
		}
		cout << " Test Ok!" << endl;
		delete []Rcopy;
		delete []ArrRange;
	}
}

void LZDocList64::readListFiles(char *inputFile, bool lowercase){
	ulong i, len, lenText;
	char fileName[300];

	std::ifstream in(inputFile);
	string line;
	std::getline(in,line);
	while(in){
	    strcpy(fileName,line.c_str());
	    //cout << "File: " << fileName << endl;
	    std::ifstream input(fileName);
		assert(input.good());
		input.seekg(0, ios_base::end);
		len = (size_t)input.tellg();
		if(len > 1){
			n += len;
			nDocs++;
		}
		input.close();
		std::getline(in,line);
	}
	in.close();

	cout << "Length of generalize text(n): " << n << ", in " << nDocs << " Documents" << endl;
	// allocate to memory for text...
	seq = new uchar[n];
	char *aux;
	std::ifstream in2(inputFile);
	lenText = 0;

	for(ulong texts=0; texts < nDocs;){
		std::getline(in2,line);
		strcpy(fileName,line.c_str());
		std::ifstream input(fileName); 			// open file
		//cout << "... File: " << fileName << endl;
		assert(input.good()); 				// check open
		input.seekg(0, ios_base::end);			// move to the end
		len = (size_t)input.tellg();			// add the final symbol (smallest)
		if(len > 1){
			aux = new char[len];
			//if (TRACE)
			//cout << "To read " << fileName << " pos: " << lenText << "..." << lenText+len-1 << endl;

			input.seekg(0, ios_base::beg);		// move back to the beginning
			if (input.good()){
				input.read(aux, len);

				len--;
				aux[len] = cutDoc;
				//cout << aux << endl;
				for (i=0; i<len; i++, lenText++){
					if((uchar)(aux[i]) <= cutDoc)
						seq[lenText] = ' ';
					else{
						if (lowercase)
							seq[lenText] = (uchar)(tolower(aux[i]));
						else
							seq[lenText] = (uchar)(aux[i]);
					}
					(charTf[seq[lenText]])++;
				}
				seq[lenText] = cutDoc;
				(charTf[(ulong)cutDoc])++;
				lenText++;
				assert(!input.fail());
				input.close();
			}else{
				cout << "Can't to open the file: <" << fileName << ">";
				exit(1);
			}
			delete [] aux;
			texts++;
		}
	}
	seq[n-1] = '\0';
	in2.close();
	char fileCpy[300];
	strcpy(fileCpy,inputFile);
	strcat(fileCpy, "_copy.txt");
	ofstream myfile;
	myfile.open (fileCpy);
	myfile << seq;
	myfile.close();
	seq[n-1] = cutDoc;

	if(TEST){
		uint DD = 0;
		for(i=0; i<n; i++){
			if (seq[i] == cutDoc)
				DD++;
		}
		if(nDocs != DD){
			cout << "ERROR with nDocs in Sequence !! " << endl;
			cout << "nDocs = " << nDocs << " != " << DD << endl;
			exit(1);
		}
	}
	if(TRACE){		// listing original sequence
		cout << endl << "T[0.." << n-1 << "]:" << endl;
		for(i=0; i<n; i++){
			if (seq[i] == cutDoc)
				cout << "$";
			else
				cout << seq[i];
		}
		cout << endl;

		cout << "charTf[0..255]: ";
		for(i=0; i<LZ78Tries64::SIGMA; i++)
			if (charTf[i])
				cout << i << "(" << charTf[i] << ") ";
		cout << endl;
	}
}

void LZDocList64::readUniqueFile(char *inputFile, bool lowercase){
	ulong i, j, len;
	n = nDocs = 0;
	ifstream input(inputFile);			// open file
	assert(input.good()); 				// check open
	input.seekg(0, ios_base::end);		// move to the end
	n = (size_t)input.tellg();
	seq = new uchar[n];
	char *aux = new char[n];

	input.seekg(0, ios_base::beg);		// move back to the beginning
	if (input.good()){
		input.read(aux, n-1);
		for (i=0; i<n-1; i++){
			if((uchar)(aux[i]) == boundSymbol){
				nDocs++;
				while((uchar)(aux[i+1]) == boundSymbol){
					seq[i] = '\n';
					(charTf['\n'])++;
					i++;
				}
				seq[i] = cutDoc;
			}else{
				if((uchar)(aux[i]) <= cutDoc)
					seq[i] = '\n';
				else{
					if (lowercase)
						seq[i] = (uchar)(tolower(aux[i]));
					else
						seq[i] = (uchar)aux[i];
				}
			}
			(charTf[seq[i]])++;
		}
		if(seq[n-2] == cutDoc){
			seq[n-2] = '\n';
			(charTf['\n'])++;
			nDocs--;
		}
		seq[n-1] = cutDoc;
		nDocs++;
		//cout << seq << endl;
		assert(!input.fail());
		input.close();
	}else{
		cout << "Can't to open the file: <" << inputFile << ">";
		exit(1);
	}
	input.close();
	delete [] aux;
	cout << "Length of generalize text(n): " << n << ", in " << nDocs << " Documents" << endl;

	for (i=j=len=0; i<n; i++){
		if(seq[i] == cutDoc){
			len=0;
			j++;
		}else
			len++;
	}
	if (j != nDocs){
		cout << "Error cutDocs Symbols = " << j << " != nDocs = " << nDocs << endl;
		exit(1);
	}

	seq[n-1] = '\0';
	char fileCpy[300];
	strcpy(fileCpy,inputFile);
	strcat(fileCpy, "_copy.txt");
	ofstream myfile;
	myfile.open (fileCpy);
	myfile << seq;
	myfile.close();
	seq[n-1] = cutDoc;

	if(TEST){
		uint DD = 0;
		for(i=0; i<n; i++){
			if (seq[i] == cutDoc)
				DD++;
		}
		if(nDocs != DD){
			cout << "ERROR with nDocs in Sequence !! " << endl;
			cout << "nDocs = " << nDocs << " != " << DD << endl;
			exit(1);
		}
	}

	if(TRACE){		// listing original sequence
		cout << "T[0.." << n-1 << "]:" << endl;
		for(i=0; i<n; i++){
			if (seq[i] == cutDoc)
				cout << "$";
			else
				cout << seq[i];
		}
		cout << endl;
		cout << "charTf[0..255]: ";
		for(i=0; i<LZ78Tries64::SIGMA; i++)
			if (charTf[i])
				cout << i << "(" << charTf[i] << ") ";
		cout << endl;
	}
}

// make priorities for hash table by probability
void LZDocList64::makePriorities(ulong* cPosTb){
	uint i, j, sigma = 1;
	uchar *prior = new uchar[LZ78Tries64::SIGMA];
	prior[0] = 0; // --> charTf[0];

	bool tr = TRACE;
	TRACE = false;

	for (i=1; i<LZ78Tries64::SIGMA; i++){
		if (charTf[i] > 0)
			sigma++;
		for (j=i; j>0 && charTf[prior[j-1]] < charTf[i]; j--){
			prior[j] = prior[j-1];
		}
		prior[j] = i; // --> charTf[i];
	}

	cout << " ## sigma = " << sigma << endl;

	if(TRACE){
		cout << "prior[0..255]: ";
		for(i=0; i<LZ78Tries64::SIGMA; i++){
			cout << i << "(" << (uint)prior[i] << ") ";
		}
		cout << endl;
	}

	uint rows = LZ78Tries64::SIGMA / LZ78Tries64::LENHASH;
	if (LZ78Tries64::SIGMA % LZ78Tries64::LENHASH)
		rows++;
	uint code = 0;
	for(uint fil=0; fil<rows; fil++){
		for(i=0; i<LZ78Tries64::LENHASH && code<LZ78Tries64::SIGMA; i++){
			if(fil%2 == 0)
				cPosTb[prior[code]] = i;
			else
				cPosTb[prior[code]] = LZ78Tries64::LENHASH-1-i;
			code++;
		}
	}

	if(TRACE){
		cout << "charPosTb[0..255]: ";
		for(i=0; i<LZ78Tries64::SIGMA; i++){
			cout << i << "(" << cPosTb[i] << ") ";
		}
		cout << endl;
	}
	delete [] prior;
	TRACE = tr;
}

// [lcur, rcur] is the current bitstring in the level 'lev' (B[lcur, rcur]). [lobj, robj] is the target interval to search in LZTrie, and [lrmq, rrmq] is the interval to query by RMQ in RevTrie
void LZDocList64::searchIntervalInRange(uint lev, ulong lcur, ulong rcur, ulong lobj, ulong robj, ulong lrmq, ulong rrmq, uint* occ, uint* nDocs){
	if (lev == h-1){	// this is a leaf !!
		if (lrmq < rrmq){
			if (lobj <= lcur){
				uint doc = getNum64(DocLZ, lcur*lgD, lgD);
				if (V[doc] == 0){
					occ[*nDocs] = doc;
					(*nDocs)++;
					//cout << " hoja -> doc " << doc << endl;
				}
				V[doc] = mark;
			}
			if (rcur <= robj){
				uint doc = getNum64(DocLZ, rcur*lgD, lgD);
				if (V[doc] == 0){
					occ[*nDocs] = doc;
					(*nDocs)++;
					//cout << " hoja -> doc " << doc << endl;
				}
				V[doc] = mark;
			}
		}else{ // then lcur == rcur
			if (Range[lev].B_rrr[lrmq] == 0 && lobj <= lcur){
				uint doc = getNum64(DocLZ, lcur*lgD, lgD);
				if (V[doc] == 0){
					occ[*nDocs] = doc;
					(*nDocs)++;
					//cout << " hoja -> doc " << doc << endl;
				}
				V[doc] = mark;
			}else{
				if (Range[lev].B_rrr[lrmq] && rcur <= robj){
					uint doc = getNum64(DocLZ, rcur*lgD, lgD);
					if (V[doc] == 0){
						occ[*nDocs] = doc;
						(*nDocs)++;
						//cout << " hoja -> doc " << doc << endl;
					}
					V[doc] = mark;
				}
			}
		}
	}else{
		if(lobj <= lcur && rcur <= robj){			// we consider the complete node --> apply RMQ...  the node is totally inside of objective
			if (lev<levRMQ){	// is there a RMQ structure in this level ?
				mark++;
				reportDocOcc2(lev, lcur, rcur, lrmq, rrmq, occ, nDocs);
			}else
				reportDocOcc2_levRange(lev, lcur, rcur, lrmq, rrmq, lobj, robj, occ, nDocs);
		}else{
			ulong med = lcur + ((rcur-lcur)>>1);

			if(robj <= med){
				// we go down only to the left...
				ulong befL0=0, befR0=0, aux;
				aux = Range[lev].B_rank0.rank(lcur);
				befL0 = Range[lev].B_rank0.rank(lrmq) - aux;
				befR0 = Range[lev].B_rank0.rank(rrmq+1) - aux;
				if (lev<h && befR0>befL0)
					searchIntervalInRange(lev+1, lcur, med, lobj, robj, lcur+befL0, lcur+befR0-1, occ, nDocs);
			}else{
				if(lobj > med){
					// we go down only to the right...
					ulong befL1=0, befR1=0, aux;
					aux = Range[lev].B_rank1.rank(lcur);
					befL1 = Range[lev].B_rank1.rank(lrmq) - aux;
					befR1 = Range[lev].B_rank1.rank(rrmq+1) - aux;
					if (befR1>befL1)
						searchIntervalInRange(lev+1, med+1, rcur, lobj, robj, med+befL1+1, med+befR1, occ, nDocs);
				}else{
					// first, we go down to the left...
					ulong befL0=0, befR0=0, aux;
					aux = Range[lev].B_rank0.rank(lcur);
					befL0 = Range[lev].B_rank0.rank(lrmq) - aux;
					befR0 = Range[lev].B_rank0.rank(rrmq+1) - aux;
					if (lev<h && befR0>befL0)
						searchIntervalInRange(lev+1, lcur, med, lobj, robj, lcur+befL0, lcur+befR0-1, occ, nDocs);

					// second, we go down to the right...
					ulong befL1=0, befR1=0;
					aux = Range[lev].B_rank1.rank(lcur);
					befL1 = Range[lev].B_rank1.rank(lrmq) - aux;
					befR1 = Range[lev].B_rank1.rank(rrmq+1) - aux;
					if (lev<h && befR1>befL1)
						searchIntervalInRange(lev+1, med+1, rcur, lobj, robj, med+befL1+1, med+befR1, occ, nDocs);
				}
			}
		}
	}
}

// this stores in mat[i] the preorder of the node that matches with a phrase p_{i..m} (find all phrase i for this m)
void LZDocList64::searchPhase_Rev(uchar* pat, uint m, ulong *mat){
	ulong d, w1, w2, sig0, pre=0, pos=1;

	while (m > 0){
		if (fRev_rrr[pre]){						// is it a fictitious node ?
			w1 = fRev_rank.rank(pre+1)-1;
			if (fictU_rrr[w1]){					// is it a fictitious node with unary path of other fictitious nodes ?
				w2 = fictU_rank.rank(w1+1);
				uint r,l=lisF_sel.select(w2);
				if (w2 < nEFHead)
					r = lisF_sel.select(w2+1);
				else
					r = nEF;
				if ((r-l) > m)
					return;
				for(;l<r && LbRevF[l]==pat[m];m--,l++);
				if (l<r || m==0)
					return;
				else{
					pre++;
					if (fRev_rrr[pre]){
						if (fictU_rrr[w1+1]){
							sig0 = treeRev->selectNext0(pos);
							if (sig0 > nRev) sig0 = nRev;
							else if (pos>=sig0) return;
							if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
								if(w2 > 1){
									d=1;
									treeRev->fwd_search(sig0-w2+1, &d, &pos);
									pos++;
								}else
									pos = sig0+1;
							}else return;
							pre = pos-treeRev->rank_1(pos-1);
							if (fRev_rrr[pre]==0)
								mat[m] = pre;
						}else m++;
					}else m++;
				}
			}else{									// one fictitious node
				sig0 = treeRev->selectNext0(pos);
				if (sig0 > nRev) sig0 = nRev;
				else if (pos>=sig0) return;
				if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
					if(w2 > 1){
						d=1;
						treeRev->fwd_search(sig0-w2+1, &d, &pos);
						pos++;
					}else
						pos = sig0+1;
				}else return;
				pre = pos-treeRev->rank_1(pos-1);
				if (fRev_rrr[pre]==0)
					mat[m] = pre;
			}
		}else{										// phrase node
			sig0 = treeRev->selectNext0(pos);
			if (sig0 > nRev) sig0 = nRev;
			else if (pos>=sig0) return;
			if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
				if(w2 > 1){
					d=1;
					treeRev->fwd_search(sig0-w2+1, &d, &pos);
					pos++;
				}else
					pos = sig0+1;
			}else return;
			pre = pos-treeRev->rank_1(pos-1);
			if (fRev_rrr[pre]==0)
				mat[m] = pre;
		}
		m--;
	}
}

// return true if the pattern is in the RevTrie and in '*pv' the respective preorder value in the tree,
// in '*x' stores its DFUDS bit position. Return false if the pattern does not exist
bool LZDocList64::searchPattern_Rev(uchar* pat, uint m, ulong *x, ulong *pv){
	ulong d, w1, w2, sig0, pre=0, pos=*x;

	while (m > 0){
		m--;
		if (fRev_rrr[pre]){						// is it a fictitious node ?
			w1 = fRev_rank.rank(pre);

			if (fictU_rrr[w1]){					// is it a fictitious node with unary path of other fictitious nodes ?
				w2 = fictU_rank.rank(w1+1);
				ulong r,l=lisF_sel.select(w2);
				if (w2 < nEFHead)
					r = lisF_sel.select(w2+1);
				else
					r = nEF;
				if(m){
					while (l<r && LbRevF[l]==pat[m]){
						m--; l++;
						if (m==0) break;
					}
					if (l<r){
						if (m==0){
							if (LbRevF[l]!=pat[0])return false;
						}else return false;
					}else{
						pre++;
						if (fRev_rrr[pre]){
							if (fictU_rrr[w1+1]){
								sig0 = treeRev->selectNext0(pos);
								if (sig0 > nRev) sig0 = nRev;
								else if (pos>=sig0) return false;

								if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
									if(w2 > 1){
										d=1;
										treeRev->fwd_search(sig0-w2+1, &d, &pos);
										pos++;
									}else
										pos = sig0+1;
								}else return false;

								pre = pos-treeRev->rank_1(pos-1);
							}else m++;
						}else m++;
					}
				}else{
					if (LbRevF[l]==pat[0]){
						*x = pos;
						*pv = pre;
						return true;
					}
					else return false;
				}
			}else{									// normal case
				sig0 = treeRev->selectNext0(pos);
				if (sig0 > nRev) sig0 = nRev;
				else if (pos>=sig0) return false;

				if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
					if(w2 > 1){
						d=1;
						treeRev->fwd_search(sig0-w2+1, &d, &pos);
						pos++;
					}else
						pos = sig0+1;
				}else return false;
				pre = pos-treeRev->rank_1(pos-1);
			}
		}else{										// normal case
			sig0 = treeRev->selectNext0(pos);
			if (sig0 > nRev) sig0 = nRev;
			else if (pos>=sig0) return false;

			if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
				if(w2 > 1){
					d=1;
					treeRev->fwd_search(sig0-w2+1, &d, &pos);
					pos++;
				}else
					pos = sig0+1;
			}else return false;
			pre = pos-treeRev->rank_1(pos-1);
		}
	}
	*x = pos;
	*pv = pre;
	return true;
}

// check if obj (LZ) matches with some point in [lRev... rRev]
bool LZDocList64::searchRevRange(uint lev, ulong lcur, ulong rcur, ulong obj, ulong lRev, ulong rRev){
	if(lRev <= lcur && rRev >= rcur)
		return true;

	if (lev == h-1){	// this is a leaf !!
		if (Range[lev].B_rrr[lRev] == 0 && obj == lcur)
			return true;
		else{
			if (Range[lev].B_rrr[lRev] && obj == rcur)
				return true;
		}
	}else{
		ulong med = lcur + ((rcur-lcur)>>1);
		if(obj <= med){	// we go down to the left...
			ulong befL0=0, befR0=0, aux;
			aux = Range[lev].B_rank0.rank(lcur);
			befL0 = Range[lev].B_rank0.rank(lRev) - aux;
			befR0 = Range[lev].B_rank0.rank(rRev+1) - aux;
			if (lev<h && befR0>befL0)
				return (searchRevRange(lev+1, lcur, med, obj, lcur+befL0, lcur+befR0-1));
		}else{
			// we go down to the right...
			ulong befL1=0, befR1=0, aux;
			aux = Range[lev].B_rank1.rank(lcur);
			befL1 = Range[lev].B_rank1.rank(lRev) - aux;
			befR1 = Range[lev].B_rank1.rank(rRev+1) - aux;
			if (befR1>befL1)
				return (searchRevRange(lev+1, med+1, rcur, obj, med+befL1+1, med+befR1));
		}
	}
	return false;
}

void LZDocList64::reportDocOcc1(ulong l, ulong r, uint *occ, uint *nDocs){
	if (l > r)
		return;
	ulong min, preRev, preLZ;
	uint doc;

	if (l==r) min = l;
	else min = rmqC->queryRMQ(l,r);
	preRev = LDocRev_rank1.rank(min+1);
	preLZ = getNum64(Node, preRev*lgNod, lgNod) + min - LDocRev_sel1.select(preRev);
	doc = getNum64(DocLZ, preLZ*lgD, lgD);

	if (V[doc] == 0){
		V[doc] = 1;
		occ[*nDocs] = doc;
		(*nDocs)++;
		//cout << " Occ 1 -> pos " << min << ", doc " << doc << endl;
		if (min > 0)
			reportDocOcc1(l, min-1, occ, nDocs);
		if (min < n-1)
			reportDocOcc1(min+1, r, occ, nDocs);
	}
}

// In occ we will be put the nOcc1 occurrences type 1 of document listing for pat[1..m]
void LZDocList64::documentListingTypeOne(uchar* pat, uint m, uint *occ, uint *nDocs1){
	ulong pv, len, fictBef, x=1, l, r;

	if (searchPattern_Rev(pat, m, &x, &pv)){
		len = treeRev->subTreeSize(x);
		fictBef = fRev_rank.rank(pv);
		len -= fRev_rank.rank(pv+len) - fictBef;
		pv -= fictBef;
		l = LDocRev_sel1.select(pv);
		if (pv+len < nod)
			r = LDocRev_sel1.select(pv+len)-1;
		else r = nDocRev-1;
		reportDocOcc1(l, r, occ, nDocs1);
	}
}

// apply RMQ on the bitstring segment B[l,r] in a level 'lev'. In the level 'lev' there is a RMQ structure !!
void LZDocList64::reportDocOcc2(uint lev, ulong lcur, ulong rcur, ulong l, ulong r, uint* occ, uint* nDocs){
	if (l > r)
		return;

	uint doc;
	ulong med, i, min, t=lev, lq=lcur, rq=rcur;
	if (l==r) min = i = l;
	else min = i = Range[lev].rmqB->queryRMQ(l,r);

	while (t < h-1){
		med = lq + (rq-lq)/2;
		if(Range[t].B_rrr[i]){			// go down to the right...
			i = med+Range[t].B_rank1.rank(i+1);
			i -= Range[t].B_rank1.rank(lq);
			lq = med+1;
		}else{								// go down to the left...
			i = lq+Range[t].B_rank0.rank(i+1)-1;
			i -= Range[t].B_rank0.rank(lq);
			rq = med;
		}
		t++;
	}
	// now, i is the preorder in lz of the occurrence type 2 to report...
	if (Range[t].B_rrr[i])
		i = rq;
	else
		i = lq;
	doc = getNum64(DocLZ, i*lgD, lgD);
	if (V[doc] != mark){
		if (V[doc] == 0){
			occ[*nDocs] = doc;
			(*nDocs)++;
			//cout << " C Occ 2 -> doc " << doc << endl;
		}//else cout << " * C Occ 2 -> doc " << doc << endl;
		V[doc] = mark;
		if (min > 0)
			reportDocOcc2(lev, lcur, rcur, l, min-1, occ, nDocs);
		if (min < nod-2)
			reportDocOcc2(lev, lcur, rcur, min+1, r, occ, nDocs);
	}
}

// apply sequential search in the range [l,r]...
void LZDocList64::reportDocOcc2_levRange(uint lev, ulong lcur, ulong rcur, ulong l, ulong r, ulong lobj, ulong robj, uint* occ, uint* nDocs){
	uint auxLev, doc;
	ulong i, obj, auxL, auxR, med;

	for(i=l; i<=r; i++){
		for(obj=i, auxLev=lev, auxL=lcur, auxR=rcur; auxLev<h-1; auxLev++){
			med = auxL + ((auxR-auxL)>>1);
			if (Range[auxLev].B_rrr[obj] == 0){
				// we go down to the left...
				obj = auxL+Range[auxLev].B_rank0.rank(obj+1);
				obj -= Range[auxLev].B_rank0.rank(auxL)+1;
				auxR = med;
			}else{
				//we go down to the right...
				obj = med+Range[auxLev].B_rank1.rank(obj+1);
				obj -= Range[auxLev].B_rank1.rank(auxL);
				auxL = med+1;
			}
		}

		if (auxLev == h-1){
			// this is a leaf --> report...
			if (Range[auxLev].B_rrr[obj])
				obj = auxR;
			else
				obj = auxL;
			doc = getNum64(DocLZ, obj*lgD, lgD);
			if (V[doc] == 0){
				occ[*nDocs] = doc;
				(*nDocs)++;
				//cout << " R Occ 2 -> doc " << doc << endl;
			}//else cout << " * R Occ 2 -> doc " << doc << endl;
			V[doc] = 1;
		}
	}
}

void LZDocList64::documentListingTypeTwoThree(uchar* pat, uint m, uint* occ, uint *nDocs, uint* nDocs2, uint* nDocs3, double* t2, double* t3){
// =========================================================================================================================
// occurrences type 2...
// =========================================================================================================================
	*t2 = getTime_ms();
	// for each partition, pref[k] stores the prefix P_{0..k}. [0]=preorder-Rev, [1]=subtree size of [0]
	// suff[k] stores the suffix P_{k+1..m-1}.: [0]:preorder-LZ, [1]:subtree size of [0].

	uint i, j, k, l, doc;
	ulong r, pv, len, fictBef, x=1;
	ulong pref[m-1][2], suff[m-1][2], M[m-1][m-1][3];
	ulong *mat = new ulong[m-1];
	for (i=0; i<(m-1); i++){
		mat[i] = pref[i][0] = pref[i][1] = suff[i][0] = suff[i][1] = 0;
		for (j=0; j<(m-1); j++)
			M[i][j][0] = M[i][j][1] = 0;
	}

	// [2] occ type 2...
	mark = 2;
	for (l=1; l<m; l++){
		x=1;
		if (searchPattern_Rev(pat, l, &x, &pv)){				// search prefix P_{0..l-1}...
			len = treeRev->subTreeSize(x);
			fictBef = fRev_rank.rank(pv);
			pref[l][0] = pv - fictBef;
			pref[l][1] = len - fRev_rank.rank(pv+len) + fictBef;

			x = 1;
			if (searchPattern_Rev(pat+l, m-l, &x, &pv)){		// search suffix P_{l..m-1}...
				if (fRev_rrr[pv]){
					suff[l-1][0] = -1;
				}else{
					pv -= fRev_rank.rank(pv+1);
					x = getNum64(Node, pv*lgNod, lgNod);					// preorder in LZTrie
					if (pv < nod-1)
						len = LDocRev_sel1(pv+1) - LDocRev_sel1(pv);
					else
						len = n - LDocRev_sel1(pv);
					suff[l-1][0] = x;
					suff[l-1][1] = len;
					searchIntervalInRangeOcc(0, 0, nod-1, x, x+len-1, pref[l][0], pref[l][0]+pref[l][1]-1, occ, nDocs, nDocs2);
				}
			}
		}else{
			x=1;
			if (searchPattern_Rev(pat+l, m-l, &x, &pv)){		// search suffix P_{l..m-1}...
				if (fRev_rrr[pv]){
					suff[l-1][0] = -1;
				}else{
					pv -= fRev_rank.rank(pv+1);
					x = getNum64(Node, pv*lgNod, lgNod);					// preorder in LZTrie
					if (pv < nod-1)
						len = LDocRev_sel1(pv+1) - LDocRev_sel1(pv);
					else
						len = n - LDocRev_sel1(pv);
					suff[l-1][0] = x;
					suff[l-1][1] = len;
				}
			}
		}
	}
	*t2 = getTime_ms() - *t2;

	//cout << "Found " << nOcc2 << " occ type 2" << endl;
	//cout << "type1, Prefixes and Suffixes..." << endl;
	//for (i=0; i<(m-1); i++)
	//	cout << i << ":(" << pref[i][0] << "," << pref[i][1] << "," << suff[i][0] << "," << suff[i][1]<< ") " << endl;

// =========================================================================================================================
// [3] occ type 3...
// =========================================================================================================================
	// M[i][j] represents to the phrase P_{i..j}, where:
	// M[i][j][0] stores the preorderLZ for the node that matches with p_{i..j}
	// M[i][j][1] stores the M's index r of phrase(K+1) which is stored in M[j+1][r]=P_{j+1..r} that join to the right with phrase(K)=P_{i..j}
	// M[i][j][2] stores the preorderLZ of phrase(K+1) that join to the right with phrase(K) = P_{i..j}

	if (m>2){
		*t3 = getTime_ms();
		for (j=m-2; j>0; j--){
			searchPhase_Rev(pat, j, mat);
			for (i=1; i<=j; i++){
				if (mat[i]){
					pv = mat[i] - fRev_rank.rank(mat[i]);//FRev->rank1(mat[i]-1);
					M[i][j][0] = getNum64(Node, pv*lgNod, lgNod);
					// to search preorderLZ of phrase(K+1) that join with phrase(K) = P_{i..j}
					M[i][j][2] = searchLZCoord(0, 0, nod-1, pv);
					for (l=r=j+1; r<m-1; r++){ 	// search concatenations...
						if(M[l][r][0] == M[i][j][2]){
							M[i][j][1]=r;
							break;
						}
					}
				}
			}
		}

		/*cout << "M (preLZ, next M's index (phrase concatenated), preLZ of next phrase)...";
		for (i=1; i<m-1; i++){
			cout << endl << "[" << i << "] ";
			for (j=i; j<m-1; j++)
				cout << j << ":(" << M[i][j][0] << "," << M[i][j][1] << "," << M[i][j][2]<< ") ";
		}
		cout << endl;*/

		// *** HERE THE MODEL FOR P_{0..m-1} IS SPLIT P IN P_{0..i-1} + P_{i..j} + P_{j+1..m-1},
		// where P_{0..i-1} is in pref[i], P_{i..j} is in M[i][j] and P_{j+1..m-1} is in suff[j]
		// and we also consider to split P_{i..j} in more phrases stored in M[][]
		for (i=1; i<m-1; i++){
			if(pref[i][0]){
				for (j=i; j<m-1; j++){
					if(M[i][j][0]){
						doc = getNum64(DocLZ, M[i][j][0]*lgD, lgD);
						if (V[doc] == 0){
							// check pref[i]+M[i][j]...
							if (searchRevRange(0, 0, nod-1, M[i][j][0], pref[i][0], pref[i][0]+pref[i][1]-1)){
								if (suff[j][0] > 0){
									// can be 'M phase' and 'suffix' concatenated ?
									if (M[i][j][2] >= suff[j][0] && M[i][j][2] <= (suff[j][0]+suff[j][1]-1)){
										// REPORT OCCURRENCE
										V[doc] = 1;
										occ[*nDocs] = doc;
										(*nDocs)++;
										(*nDocs3)++;
									}
								}

								// to search concatenations of full phrase...
								l=i;r=j;
								while(M[l][r][1]){
									k=l;
									l = r+1;
									r = M[k][r][1];
									if (suff[r][0] > 0){
										// can be 'M phase' and 'suffix' concatenated ?
										if (M[l][r][2] >= suff[r][0] && M[l][r][2] <= (suff[r][0]+suff[r][1]-1)){
											// REPORT OCCURRENCE
											V[doc] = 1;
											occ[*nDocs] = doc;
											(*nDocs)++;
											(*nDocs3)++;
										}
									}
								}
							}
						}
					}
				}
			}
		}
		*t3 = getTime_ms() - *t3;
	}else
		*t3 = 0;

	//cout << "Found " << nOcc3 << " occ type 3" << endl;
}

void LZDocList64::documentListing(uchar* pat, uint m, uint* occ, uint* nDocs){
	ulong pv, len, fictBef, x=1, l, r;
	uint i, j, k, doc;

	//mark = 1;	// this is the mark for occurrences type one that we will mark in V[]
	if (searchPattern_Rev(pat, m, &x, &pv)){
		len = treeRev->subTreeSize(x);
		fictBef = fRev_rank.rank(pv);
		len -= fRev_rank.rank(pv+len) - fictBef;
		pv -= fictBef;
		l = LDocRev_sel1.select(pv);
		if (pv+len < nod)
			r = LDocRev_sel1.select(pv+len)-1;
		else r = nDocRev-1;
		reportDocOcc1(l, r, occ, nDocs);
	}

	// =========================================================================================================================
	// occurrences type 2...
	// =========================================================================================================================
	// for each partition, pref[k] stores the prefix P_{0..k}. [0]=preorder-Rev, [1]=subtree size of [0]
	// suff[k] stores the suffix P_{k+1..m-1}.: [0]:preorder-LZ, [1]:subtree size of [0].

	ulong pref[m-1][2], suff[m-1][2], M[m-1][m-1][3];
	ulong *mat = new ulong[m-1];
	for (i=0; i<(m-1); i++){
		mat[i] = pref[i][0] = pref[i][1] = suff[i][0] = suff[i][1] = 0;
		for (j=0; j<(m-1); j++)
			M[i][j][0] = M[i][j][1] = 0;
	}

	// [2] occ type 2...
	for (l=1; l<m; l++){
		x=1;
		if (searchPattern_Rev(pat, l, &x, &pv)){				// search prefix P_{0..l-1}...
			len = treeRev->subTreeSize(x);
			fictBef = fRev_rank.rank(pv);
			pref[l][0] = pv - fictBef;
			pref[l][1] = len - fRev_rank.rank(pv+len) + fictBef;

			x = 1;
			if (searchPattern_Rev(pat+l, m-l, &x, &pv)){		// search suffix P_{l..m-1}...
				if (fRev_rrr[pv]){
					suff[l-1][0] = -1;
				}else{
					pv -= fRev_rank.rank(pv+1);
					x = getNum64(Node, pv*lgNod, lgNod);					// preorder in LZTrie
					if (pv < nod-1)
						len = LDocRev_sel1(pv+1) - LDocRev_sel1(pv);
					else
						len = n - LDocRev_sel1(pv);
					suff[l-1][0] = x;
					suff[l-1][1] = len;
					mark = 2;
					searchIntervalInRange(0, 0, nod-1, x, x+len-1, pref[l][0], pref[l][0]+pref[l][1]-1, occ, nDocs);
				}
			}
		}else{
			x = 1;
			if (searchPattern_Rev(pat+l, m-l, &x, &pv)){		// search suffix P_{l..m-1}...
				if (fRev_rrr[pv]){
					suff[l-1][0] = -1;
				}else{
					pv -= fRev_rank.rank(pv+1);
					x = getNum64(Node, pv*lgNod, lgNod);					// preorder in LZTrie
					if (pv < nod-1)
						len = LDocRev_sel1(pv+1) - LDocRev_sel1(pv);
					else
						len = n - LDocRev_sel1(pv);
					suff[l-1][0] = x;
					suff[l-1][1] = len;
				}
			}
		}
	}
	/*cout << "Found " << *nOcc2 << " occ type 2" << endl;
	cout << "type1, Prefixes and Suffixes..." << endl;
	for (i=0; i<(m-1); i++)
		cout << i << ":(" << pref[i][0] << "," << pref[i][1] << "," << suff[i][0] << "," << suff[i][1]<< ") " << endl;
	 */

// =========================================================================================================================
// [3] occ type 3...
// =========================================================================================================================
	// M[i][j] represents to the phrase P_{i..j}, where:
	// M[i][j][0] stores the preorderLZ for the node that matches with p_{i..j}
	// M[i][j][1] stores the M's index r of phrase(K+1) which is stored in M[j+1][r]=P_{j+1..r} that join to the right with phrase(K)=P_{i..j}
	// M[i][j][2] stores the preorderLZ of phrase(K+1) that join to the right with phrase(K) = P_{i..j}
	if (m>2){
		for (j=m-2; j>0; j--){
			searchPhase_Rev(pat, j, mat);
			for (i=1; i<=j; i++){
				if (mat[i]){
					pv = mat[i] - fRev_rank.rank(mat[i]);//FRev->rank1(mat[i]-1);
					M[i][j][0] = getNum64(Node, pv*lgNod, lgNod);
					//cout << "M[i][j][0] = " << M[i][j][0] << endl;
					// to search preorderLZ of phrase(K+1) that join with phrase(K) = P_{i..j}
					M[i][j][2] = searchLZCoord(0, 0, nod-1, pv);
					//cout << "M[i][j][2] = " << M[i][j][2] << endl;
					for (l=r=j+1; r<m-1; r++){ 	// search concatenations...
						if(M[l][r][0] == M[i][j][2]){
							M[i][j][1]=r;
							break;
						}
					}
				}
			}
		}

		/*cout << "M (preLZ, next M's index (phrase concatenated), preLZ of next phrase)...";
		for (i=1; i<m-1; i++){
			cout << endl << "[" << i << "] ";
			for (j=i; j<m-1; j++)
				cout << j << ":(" << M[i][j][0] << "," << M[i][j][1] << "," << M[i][j][2]<< ") ";
		}
		cout << endl;*/


		// *** HERE THE MODEL FOR P_{0..m-1} IS SPLIT P IN P_{0..i-1} + P_{i..j} + P_{j+1..m-1},
		// where P_{0..i-1} is in pref[i], P_{i..j} is in M[i][j] and P_{j+1..m-1} is in suff[j]
		// and we also consider to split P_{i..j} in more phrases stored in M[][]
		for (i=1; i<m-1; i++){
			if(pref[i][0]){
				for (j=i; j<m-1; j++){
					if(M[i][j][0]){
						doc = getNum64(DocLZ, M[i][j][0]*lgD, lgD);

						// check pref[i]+M[i][j]...
						if (searchRevRange(0, 0, nod-1, M[i][j][0], pref[i][0], pref[i][0]+pref[i][1]-1)){
							if (suff[j][0] > 0){
								// can be 'M phase' and 'suffix' concatenated ?
								if (M[i][j][2] >= suff[j][0] && M[i][j][2] <= (suff[j][0]+suff[j][1]-1)){
									if (V[doc] == 0){
										// REPORT OCCURRENCE
										V[doc] = 1;
										occ[*nDocs] = doc;
										(*nDocs)++;
									}
								}
							}

							// to search concatenations of full phrase...
							l=i;r=j;
							while(M[l][r][1]){
								k=l;
								l = r+1;
								r = M[k][r][1];
								if (suff[r][0] > 0){
									// can be 'M phase' and 'suffix' concatenated ?
									if (M[l][r][2] >= suff[r][0] && M[l][r][2] <= (suff[r][0]+suff[r][1]-1)){
										if (V[doc] == 0){
											// REPORT OCCURRENCE
											V[doc] = 1;
											occ[*nDocs] = doc;
											(*nDocs)++;
										}
									}
								}
							}

						}
					}
				}
			}
		}
	}
	//cout << "Found " << nOcc << " occ" << endl;
}

// save the Data Structure in file 'fileName'
void LZDocList64::saveDS(bool showSize){
	char str[100];
	char *fileName = new char[300];
	cout << "Save data structure in folder " << dirStore << endl;

	strcpy(fileName, dirStore);
	strcat(fileName, "dataStructures.dl");
	ofstream os (fileName, ios::binary);
	cout << "   Saving. Data structure size: " << sizeDS << endl;

	os.write((const char*)&n, sizeof(ulong));
	os.write((const char*)&h, sizeof(uint));
	os.write((const char*)&levRMQ, sizeof(uint));
	os.write((const char*)&nDocs, sizeof(uint));
	os.write((const char*)&lgD, sizeof(uchar));
	os.write((const char*)&cutDoc, sizeof(char));
	os.write((const char*)&nod, sizeof(ulong));
	os.write((const char*)&lgNod, sizeof(uint));
	os.write((const char*)&nF, sizeof(ulong));
	os.write((const char*)&nEF, sizeof(ulong));
	os.write((const char*)&nEFHead, sizeof(uint));
	os.write((const char*)&nodRev, sizeof(ulong));
	os.write((const char*)&nLz, sizeof(ulong));
	os.write((const char*)&nRev, sizeof(ulong));

	ulong sizeDSav = 7*sizeof(ulong) + 5*sizeof(uint) + sizeof(char) + sizeof(uchar);	// size for variables
	cout << " .- Variables " << sizeDSav << " Bytes = " << sizeDSav << " bpc" << endl;

	sizeDSav += nDocs*sizeof(suint);
	cout << " .- V[1..nDocs] " << nDocs*sizeof(suint) << " Bytes = " << nDocs*sizeof(suint)*8.0/(float)this->n << " bpc" << endl;

	os.write((const char*)LbRev, nodRev*sizeof(uchar));					// save LbRevF[]
	sizeDSav += nodRev*sizeof(uchar);
	if(showSize) cout << " .- LbRev[] " << nodRev*sizeof(uchar) << " Bytes" << endl;

	os.write((const char*)LbRevF, nEF*sizeof(uchar));					// save LbRevF[]
	sizeDSav += nEF*sizeof(uchar);
	if(showSize) cout << " .- LbRevF[] " << nEF*sizeof(uchar) << " Bytes" << endl;

	ulong size = nod*lgD/W64;
	if ((nod*lgD)%W64)
		size++;
	os.write((const char*)DocLZ, size*sizeof(ulong));				// save DocLZ[]
	sizeDSav += size*sizeof(ulong);
	if(showSize) cout << " .- DocLZ[] " << size*sizeof(ulong) << " Bytes" << endl;

	size = nod*lgNod/W64;
	if ((nod*lgNod)%W64)
		size++;
	os.write((const char*)Node, size*sizeof(ulong));				// save Node[]
	sizeDSav += size*sizeof(ulong);
	if(showSize) cout << " .- Node[] " << size*sizeof(ulong) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fRev_rrr.dl");
	store_to_file(fRev_rrr, fileName);
	sizeDSav += size_in_bytes(fRev_rrr);
	if(showSize) cout << " .- fRev_rrr " << size_in_bytes(fRev_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fRev_rank.dl");
	store_to_file(fRev_rank, fileName);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fictU_rrr.dl");
	store_to_file(fictU_rrr, fileName);
	sizeDSav += size_in_bytes(fictU_rrr);
	if(showSize) cout << " .- fictU_rrr " << size_in_bytes(fictU_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fictU_rank.dl");
	store_to_file(fictU_rank, fileName);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "listF_rrr.dl");
	store_to_file(listF_rrr, fileName);
	sizeDSav += size_in_bytes(listF_rrr);
	if(showSize) cout << " .- listF_rrr " << size_in_bytes(listF_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "lisF_sel.dl");
	store_to_file(lisF_sel, fileName);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "LDocRev_rrr.dl");
	store_to_file(LDocRev_rrr, fileName);
	sizeDSav += size_in_bytes(LDocRev_rrr);
	if(showSize) cout << " .- LDocRev_rrr " << size_in_bytes(LDocRev_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "LDocRev_rank1.dl");
	store_to_file(LDocRev_rank1, fileName);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "LDocRev_sel1.dl");
	store_to_file(LDocRev_sel1, fileName);

	// Range
	sizeDSav += h*sizeof(LevelRange);
	cout << "  .- Range[]'s pointers " << h*sizeof(LevelRange) << " = " << (h*sizeof(LevelRange))*8.0/(float)this->n << " bpc" << endl;
	ulong sizeRMQsRange, sizeBitMapsRange;
	sizeRMQsRange = sizeBitMapsRange = 0;

	cout << "levRMQ = " << levRMQ << endl;
	cout << "h = " << h << endl;

	for(uint i=0; i<this->h; i++){
		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcpy(str, "");
		sprintf(str, "Range_%d.B_rrr.dl", i);
		strcat(fileName, str);
		store_to_file(Range[i].B_rrr, fileName);
		sizeBitMapsRange += size_in_bytes(Range[i].B_rrr);

		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcpy(str, "");
		sprintf(str, "Range_%d.B_rank0.dl", i);
		strcat(fileName, str);
		store_to_file(Range[i].B_rank0, fileName);

		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcpy(str, "");
		sprintf(str, "Range_%d.B_rank1.dl", i);
		strcat(fileName, str);
		store_to_file(Range[i].B_rank1, fileName);

		if (i<levRMQ){
			strcpy(fileName, "");
			strcpy(fileName, dirStore);
			strcpy(str, "");
			sprintf(str, "Range_%d.rmq.dl", i);
			strcat(fileName, str);
			Range[i].rmqB->saveDS(fileName);
			sizeRMQsRange += Range[i].rmqB->getSize();
		}
	}

	if(showSize){
		cout << " .- Range[]'s BitVector(RRR) " << sizeBitMapsRange << " Bytes" << endl;
		cout << " .- Range[]'s RMQs " << sizeRMQsRange << " Bytes" << endl;
	}
	sizeDSav += sizeBitMapsRange;
	sizeDSav += sizeRMQsRange;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "RevTrie.dl");
	cout << "***   Saving RevTrie... " << endl;
	sizeDSav += treeRev->saveDS(fileName, showSize);;
	if(showSize) cout << " .- treeRev " << treeRev->sizeRMM << " Bytes =? " << treeRev->saveDS(fileName, showSize) << endl;

	// RMQ
	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "rmq.dl");
	rmqC->saveDS(fileName);
	sizeDSav += rmqC->getSize();
	if(showSize) cout << " .- rmqC " << rmqC->getSize() << " Bytes" << endl;

	cout << "   Total bytes saved from data structure DL_LZ = " << sizeDSav << endl;
}

// load the Data Structure from the file 'fileName'
void LZDocList64::loadDS(bool showSize){
	char str[100];
	cout << " Load data structure from " << dirStore << endl;
	char *fileName = new char[300];

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "dataStructures.dl");
	ifstream is(fileName, ios::binary);

	is.read((char*)&n, sizeof(ulong));
	is.read((char*)&h, sizeof(uint));
	is.read((char*)&levRMQ, sizeof(uint));
	is.read((char*)&nDocs, sizeof(uint));
	is.read((char*)&lgD, sizeof(uchar));
	is.read((char*)&cutDoc, sizeof(char));
	is.read((char*)&nod, sizeof(ulong));
	is.read((char*)&lgNod, sizeof(uint));
	is.read((char*)&nF, sizeof(ulong));
	is.read((char*)&nEF, sizeof(ulong));
	is.read((char*)&nEFHead, sizeof(uint));
	is.read((char*)&nodRev, sizeof(ulong));
	is.read((char*)&nLz, sizeof(ulong));
	is.read((char*)&nRev, sizeof(ulong));

	sizeDS = 7*sizeof(ulong) + 5*sizeof(uint) + sizeof(char) + sizeof(uchar);	// size for variables
	if(showSize) cout << " ** Size of Variables " << sizeDS << " Bytes = " << sizeDS << " bpc" << endl;

	this->V = new suint[nDocs];
	sizeDS += nDocs*sizeof(suint);
	if(showSize) cout << " ** Size of V[1..nDocs] " << nDocs*sizeof(suint) << " Bytes = " << nDocs*sizeof(suint)*8.0/(float)this->n << " bpc" << endl;

	LbRev = new uchar[nodRev];
	is.read((char*)LbRev, nodRev*sizeof(uchar));
	sizeDS += nodRev*sizeof(uchar);
	if(showSize) cout << " ** size of LbRev[] " << nodRev*sizeof(uchar) << " Bytes" << endl;

	LbRevF = new uchar[nEF];
	is.read((char*)LbRevF, nEF*sizeof(uchar));
	sizeDS += nEF*sizeof(uchar);
	if(showSize) cout << " ** size of LbRevF[] " << nEF*sizeof(uchar) << " Bytes" << endl;

	ulong size = nod*lgD/W64;
	if ((nod*lgD)%W64)
		size++;
	DocLZ = new ulong[size];
	is.read((char*)DocLZ, size*sizeof(ulong));
	sizeDS += size*sizeof(ulong);
	if(showSize) cout << " ** size of DocLZ[] " << size*sizeof(ulong) << " Bytes" << endl;

	size = nod*lgNod/W64;
	if ((nod*lgNod)%W64)
		size++;
	Node = new ulong[size];
	is.read((char*)Node, size*sizeof(ulong));
	sizeDS += size*sizeof(ulong);
	if(showSize) cout << " ** size of Node[] " << size*sizeof(ulong) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fRev_rrr.dl");
	load_from_file(fRev_rrr, fileName);
	sizeDS += size_in_bytes(fRev_rrr);
	if(showSize) cout << " ** size of fRev_rrr " << size_in_bytes(fRev_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fRev_rank.dl");
	load_from_file(fRev_rank, fileName);
	util::init_support(fRev_rank, &fRev_rrr);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fictU_rrr.dl");
	load_from_file(fictU_rrr, fileName);
	sizeDS += size_in_bytes(fictU_rrr);
	if(showSize) cout << " ** size of fictU_rrr " << size_in_bytes(fictU_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fictU_rank.dl");
	load_from_file(fictU_rank, fileName);
	util::init_support(fictU_rank, &fictU_rrr);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "listF_rrr.dl");
	load_from_file(listF_rrr, fileName);
	sizeDS += size_in_bytes(listF_rrr);
	if(showSize) cout << " ** size of listF_rrr " << size_in_bytes(listF_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "lisF_sel.dl");
	load_from_file(lisF_sel, fileName);
	util::init_support(lisF_sel, &listF_rrr);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "LDocRev_rrr.dl");
	load_from_file(LDocRev_rrr, fileName);
	sizeDS += size_in_bytes(LDocRev_rrr);
	if(showSize) cout << " ** size of LDocRev_rrr " << size_in_bytes(LDocRev_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "LDocRev_rank1.dl");
	load_from_file(LDocRev_rank1, fileName);
	util::init_support(LDocRev_rank1, &LDocRev_rrr);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "LDocRev_sel1.dl");
	load_from_file(LDocRev_sel1, fileName);
	util::init_support(LDocRev_sel1, &LDocRev_rrr);

	sizeUpRange = sizeDS;

	// Range
	Range = new LevelRange[h];
	sizeDS += h*sizeof(LevelRange);
	if(showSize) cout << " ** Size of Range[]'s pointers " << h*sizeof(LevelRange) << " = " << (h*sizeof(LevelRange))*8.0/(float)this->n << " bpc" << endl;

	ulong sizeRMQsRange, sizeBitMapsRange;
	sizeRMQsRange = sizeBitMapsRange = 0;

	for(uint i=0; i<h; i++){
		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcpy(str, "");
		sprintf(str, "Range_%d.B_rrr.dl", i);
		strcat(fileName, str);
		load_from_file(Range[i].B_rrr, fileName);
		sizeBitMapsRange += size_in_bytes(Range[i].B_rrr);

		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcpy(str, "");
		sprintf(str, "Range_%d.B_rank0.dl", i);
		strcat(fileName, str);
		load_from_file(Range[i].B_rank0, fileName);
		util::init_support(Range[i].B_rank0, &(Range[i].B_rrr));

		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcpy(str, "");
		sprintf(str, "Range_%d.B_rank1.dl", i);
		strcat(fileName, str);
		load_from_file(Range[i].B_rank1, fileName);
		util::init_support(Range[i].B_rank1, &(Range[i].B_rrr));

		if (i<levRMQ){
			strcpy(fileName, "");
			strcpy(fileName, dirStore);
			strcpy(str, "");
			sprintf(str, "Range_%d.rmq.dl", i);
			strcat(fileName, str);
			Range[i].rmqB = new RMQRMM64(fileName);
			cout << " ..... Range[i].rmqB->getSize  " << Range[i].rmqB->getSize() << " Bytes" << endl;
			sizeRMQsRange += Range[i].rmqB->getSize();
		}
	}
	is.close();
	if(showSize){
		cout << " .- Range[]'s BitVector(RRR) " << sizeBitMapsRange << " Bytes" << endl;
		cout << " .- Range[]'s RMQs " << sizeRMQsRange << " Bytes" << endl;
	}
	sizeDS += sizeBitMapsRange;
	sizeDS += sizeRMQsRange;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "RevTrie.dl");
	treeRev = new RangeMMTree64(fileName, false);
	treeRev->labels = LbRev;
	sizeDS += treeRev->sizeRMM;
	sizeUpRange += treeRev->sizeRMM;
	if(showSize) cout << " ** size of treeRev " << treeRev->sizeRMM << " Bytes" << endl;

	// RMQ
	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "rmq.dl");
	rmqC = new RMQRMM64(fileName);
	sizeDS += rmqC->getSize();
	sizeUpRange += rmqC->getSize();
	if(showSize) cout << " ** size of rmqC " << rmqC->getSize() << " Bytes" << endl;

	cout << sizeDS << " Bytes Loaded !!" << endl;
	cout << "sizeUpRange = " << sizeUpRange << endl;
}

LZDocList64::~LZDocList64() {
	sizeDS = n = nDocs = lgD = cutDoc = nod = lgNod = 0;
	nF = nEF = nEFHead = nodRev = nLz = nRev = 0;

	delete []EndDocs;
	delete []charTf;
	delete []PLz;
	delete []PRev;
	delete []LbRev;
	delete []LbRevF;
	delete []DocLZ;
	delete []Node;

	for(uint i=0; i<h; i++){
		if (i<levRMQ)
			Range[i].rmqB->~RMQRMM64();
	}
	delete []Range;


	tries->~LZ78Tries64();

	cout << "destroyeing treeRev..." << endl;
	treeRev->~RangeMMTree64();
	cout << "destroyeing rmqC..." << endl;
	rmqC->~RMQRMM64();
}

// ======================================================================================================================
// THE NEXT METHODS ARE ONLY FOR TEST:
// ======================================================================================================================

// In occ we will be put the nOcc1 occurrences type 1 of document listing for pat[1..m]
void LZDocList64::documentListingTypeOne_test(uchar* pat, uint m, uint *occ, ulong *nOcc, uint *nDocs1){
	ulong pv, len, fictBef, x=1, l, r;

	if (searchPattern_Rev(pat, m, &x, &pv)){
		len = treeRev->subTreeSize(x);
		fictBef = fRev_rank.rank(pv);
		len -= fRev_rank.rank(pv+len) - fictBef;
		pv -= fictBef;
		l = LDocRev_sel1.select(pv);
		if (pv+len < nod)
			r = LDocRev_sel1.select(pv+len)-1;
		else r = nDocRev-1;

		//cout << "Interval occ in LZDocList64 = [" << l << ", " << r << "]" << endl;

		reportDocOcc1(l, r, occ, nDocs1);

		*nOcc = r-l+1;
	}else
		*nOcc = 0;
}

void LZDocList64::documentListingTypeTwoThree_test(uchar* pat, uint m, uint* occ2, ulong* nOcc2, uint* nDocs2, uint* occ3, ulong* nOcc3, uint* nDocs3, uint* auxV2, bool* auxV3){
// =========================================================================================================================
// occurrences type 2...
// =========================================================================================================================
	// for each partition, pref[k] stores the prefix P_{0..k}. [0]=preorder-Rev, [1]=subtree size of [0]
	// suff[k] stores the suffix P_{k+1..m-1}.: [0]:preorder-LZ, [1]:subtree size of [0].

	uint i, j, k, l, doc;
	ulong r, pv, len, fictBef, x=1;
	ulong pref[m-1][2], suff[m-1][2], M[m-1][m-1][3];
	ulong *mat = new ulong[m-1];
	for (i=0; i<(m-1); i++){
		mat[i] = pref[i][0] = pref[i][1] = suff[i][0] = suff[i][1] = 0;
		for (j=0; j<(m-1); j++)
			M[i][j][0] = M[i][j][1] = 0;
	}

	// [2] occ type 2...
	mark = 2;
	for (l=1; l<m; l++){
		x=1;
		if (searchPattern_Rev(pat, l, &x, &pv)){				// search prefix P_{0..l-1}...
			len = treeRev->subTreeSize(x);
			fictBef = fRev_rank.rank(pv);
			pref[l][0] = pv - fictBef;
			pref[l][1] = len - fRev_rank.rank(pv+len) + fictBef;

			x = 1;
			if (searchPattern_Rev(pat+l, m-l, &x, &pv)){		// search suffix P_{l..m-1}...
				if (fRev_rrr[pv]){
					suff[l-1][0] = -1;
				}else{
					pv -= fRev_rank.rank(pv+1);
					x = getNum64(Node, pv*lgNod, lgNod);					// preorder in LZTrie
					if (pv < nod-1)
						len = LDocRev_sel1(pv+1) - LDocRev_sel1(pv);
					else
						len = n - LDocRev_sel1(pv);
					suff[l-1][0] = x;
					suff[l-1][1] = len;

					//cout << "search in ["<<x<<","<<x+len-1<<"]x["<<pref[l][0]<<","<<pref[l][0]+pref[l][1]-1<<"]"<<endl;
					searchIntervalInRange_test(0, 0, nod-1, x, x+len-1, pref[l][0], pref[l][0]+pref[l][1]-1, occ2, nOcc2, nDocs2, auxV2);
				}
			}
		}else{
			x = 1;
			if (searchPattern_Rev(pat+l, m-l, &x, &pv)){		// search suffix P_{l..m-1}...
				if (fRev_rrr[pv]){
					suff[l-1][0] = -1;
				}else{
					pv -= fRev_rank.rank(pv+1);//FRev->rank1(pv);
					x = getNum64(Node, pv*lgNod, lgNod);					// preorder in LZTrie
					if (pv < nod-1)
						len = LDocRev_sel1(pv+1) - LDocRev_sel1(pv);
					else
						len = n - LDocRev_sel1(pv);
					suff[l-1][0] = x;
					suff[l-1][1] = len;
				}
			}
		}
	}
	//cout << "Found " << *nOcc2 << " occ type 2" << endl;
	/*cout << "type1, Prefixes and Suffixes..." << endl;
	for (i=0; i<(m-1); i++)
		cout << i << ":(" << pref[i][0] << "," << pref[i][1] << "," << suff[i][0] << "," << suff[i][1]<< ") " << endl;
	 */


// =========================================================================================================================
// [3] occ type 3...
// =========================================================================================================================
	// M[i][j] represents to the phrase P_{i..j}, where:
	// M[i][j][0] stores the preorderLZ for the node that matches with p_{i..j}
	// M[i][j][1] stores the M's index r of phrase(K+1) which is stored in M[j+1][r]=P_{j+1..r} that join to the right with phrase(K)=P_{i..j}
	// M[i][j][2] stores the preorderLZ of phrase(K+1) that join to the right with phrase(K) = P_{i..j}
	if (m>2){
		for (j=m-2; j>0; j--){
			searchPhase_Rev(pat, j, mat);
			for (i=1; i<=j; i++){
				if (mat[i]){
					pv = mat[i] - fRev_rank.rank(mat[i]);//FRev->rank1(mat[i]-1);
					M[i][j][0] = getNum64(Node, pv*lgNod, lgNod);
					//cout << "M[i][j][0] = " << M[i][j][0] << endl;
					// to search preorderLZ of phrase(K+1) that join with phrase(K) = P_{i..j}
					M[i][j][2] = searchLZCoord(0, 0, nod-1, pv);
					//cout << "M[i][j][2] = " << M[i][j][2] << endl;
					for (l=r=j+1; r<m-1; r++){ 	// search concatenations...
						if(M[l][r][0] == M[i][j][2]){
							M[i][j][1]=r;
							break;
						}
					}
				}
			}
		}

		/*cout << "M (preLZ, next M's index (phrase concatenated), preLZ of next phrase)...";
		for (i=1; i<m-1; i++){
			cout << endl << "[" << i << "] ";
			for (j=i; j<m-1; j++)
				cout << j << ":(" << M[i][j][0] << "," << M[i][j][1] << "," << M[i][j][2]<< ") ";
		}
		cout << endl;*/


		// *** HERE THE MODEL FOR P_{0..m-1} IS SPLIT P IN P_{0..i-1} + P_{i..j} + P_{j+1..m-1},
		// where P_{0..i-1} is in pref[i], P_{i..j} is in M[i][j] and P_{j+1..m-1} is in suff[j]
		// and we also consider to split P_{i..j} in more phrases stored in M[][]
		for (i=1; i<m-1; i++){
			if(pref[i][0]){
				for (j=i; j<m-1; j++){
					if(M[i][j][0]){
						doc = getNum64(DocLZ, M[i][j][0]*lgD, lgD);

						// check pref[i]+M[i][j]...
						if (searchRevRange(0, 0, nod-1, M[i][j][0], pref[i][0], pref[i][0]+pref[i][1]-1)){
							if (suff[j][0] > 0){
								// can be 'M phase' and 'suffix' concatenated ?
								if (M[i][j][2] >= suff[j][0] && M[i][j][2] <= (suff[j][0]+suff[j][1]-1)){
									(*nOcc3)++;
									if (auxV3[doc] == 0){
										// REPORT OCCURRENCE
										auxV3[doc] = 1;
										occ3[*nDocs3] = doc;
										(*nDocs3)++;
										//cout << "Doc reported type3: " << doc << endl;
									}
								}
							}

							// to search concatenations of full phrase...
							l=i;r=j;
							while(M[l][r][1]){
								k=l;
								l = r+1;
								r = M[k][r][1];
								if (suff[r][0] > 0){
									// can be 'M phase' and 'suffix' concatenated ?
									if (M[l][r][2] >= suff[r][0] && M[l][r][2] <= (suff[r][0]+suff[r][1]-1)){
										(*nOcc3)++;
										if (auxV3[doc] == 0){
											// REPORT OCCURRENCE
											auxV3[doc] = 1;
											occ3[*nDocs3] = doc;
											(*nDocs3)++;
											//cout << "Doc reported type3: " << doc << endl;
										}
									}
								}
							}

						}
					}
				}
			}
		}
	}
	//cout << "Found " << nOcc3 << " occ type 3" << endl;
}

// apply sequential search in the range [l,r]...
void LZDocList64::reportDocOcc2_levRangeOcc(uint lev, ulong lcur, ulong rcur, ulong l, ulong r, ulong lobj, ulong robj, uint* occ, uint* nDocs, uint* nDocs2){
	uint auxLev, doc;
	ulong i, obj, auxL, auxR, med;

	for(i=l; i<=r; i++){
		for(obj=i, auxLev=lev, auxL=lcur, auxR=rcur; auxLev<h-1; auxLev++){
			med = auxL + ((auxR-auxL)>>1);
			if (Range[auxLev].B_rrr[obj] == 0){
				// we go down to the left...
				obj = auxL+Range[auxLev].B_rank0.rank(obj+1);
				obj -= Range[auxLev].B_rank0.rank(auxL)+1;
				auxR = med;
			}else{
				//we go down to the right...
				obj = med+Range[auxLev].B_rank1.rank(obj+1);
				obj -= Range[auxLev].B_rank1.rank(auxL);
				auxL = med+1;
			}
		}

		if (auxLev == h-1){
			// this is a leaf --> report...
			if (Range[auxLev].B_rrr[obj])
				obj = auxR;
			else
				obj = auxL;
			doc = getNum64(DocLZ, obj*lgD, lgD);
			if (V[doc] == 0){
				occ[*nDocs] = doc;
				(*nDocs)++;
				(*nDocs2)++;
				//cout << " R Occ 2 -> doc " << doc << endl;
			}//else cout << " * R Occ 2 -> doc " << doc << endl;
			V[doc] = 1;
		}
	}
}

// apply sequential search in the range [l,r]...
void LZDocList64::reportDocOcc2_levRange_test(uint lev, ulong lcur, ulong rcur, ulong l, ulong r, ulong lobj, ulong robj, uint* occ2, ulong* nOcc2, uint* nDocs2, uint* auxV2){
	uint auxLev, doc;
	ulong i, obj, auxL, auxR, med;

	for(i=l; i<=r; i++){
		for(obj=i, auxLev=lev, auxL=lcur, auxR=rcur; auxLev<h-1; auxLev++){
			med = auxL + ((auxR-auxL)>>1);
			if (Range[auxLev].B_rrr[obj] == 0){
				// we go down to the left...
				obj = auxL+Range[auxLev].B_rank0.rank(obj+1);
				obj -= Range[auxLev].B_rank0.rank(auxL)+1;
				auxR = med;
			}else{
				//we go down to the right...
				obj = med+Range[auxLev].B_rank1.rank(obj+1);
				obj -= Range[auxLev].B_rank1.rank(auxL);
				auxL = med+1;
			}
		}

		if (auxLev == h-1){
			// this is a leaf --> report...
			if (Range[auxLev].B_rrr[obj])
				obj = auxR;
			else
				obj = auxL;

			doc = getNum64(DocLZ, obj*lgD, lgD);
			//cout << " chequear levR Occ 2 -> doc " << doc << endl;
			if (auxV2[doc] == 0){
				occ2[*nDocs2] = doc;
				(*nDocs2)++;
				auxV2[doc] = 1;
				//cout << " R Occ 2 -> doc " << doc << endl;
			}//else cout << " DESCARTADO... R Occ 2 -> doc " << doc << endl;
		}
	}
}

// apply RMQ on the bitstring segment B[l,r] in a level 'lev'. In the level 'lev' there is a RMQ structure !!
void LZDocList64::reportDocOcc2Occ(uint lev, ulong lcur, ulong rcur, ulong l, ulong r, uint* occ, uint* nDocs, uint *nDocs2){
	if (l > r)
		return;

	uint doc;
	ulong med, i, min, t=lev, lq=lcur, rq=rcur;
	if (l==r) min = i = l;
	else min = i = Range[lev].rmqB->queryRMQ(l,r);

	while (t < h-1){
		med = lq + ((rq-lq)>>1);
		if(Range[t].B_rrr[i]){			// go down to the right...
			i = med+Range[t].B_rank1.rank(i+1);
			i -= Range[t].B_rank1.rank(lq);
			lq = med+1;
		}else{								// go down to the left...
			i = lq+Range[t].B_rank0.rank(i+1)-1;
			i -= Range[t].B_rank0.rank(lq);
			rq = med;
		}
		t++;
	}

	// now, i is the preorder in LZ of the occurrence type 2 to report...
	if (Range[t].B_rrr[i])
		i = rq;
	else
		i = lq;
	doc = getNum64(DocLZ, i*lgD, lgD);
	if (V[doc] != mark){
		if (V[doc] == 0){
			occ[*nDocs] = doc;
			(*nDocs)++;
			(*nDocs2)++;
			//cout << " C Occ 2 -> doc " << doc << endl;
		}
		V[doc] = mark;
		if (min > 0)
			reportDocOcc2Occ(lev, lcur, rcur, l, min-1, occ, nDocs, nDocs2);
		if (min < nod-2)
			reportDocOcc2Occ(lev, lcur, rcur, min+1, r, occ, nDocs, nDocs2);
	}
}

// apply RMQ on the bitstring segment B[l,r] in a level 'lev'. In the level 'lev' there is a RMQ structure !!
void LZDocList64::reportDocOcc2_test(uint lev, ulong lcur, ulong rcur, ulong l, ulong r, uint* occ2, ulong* nOcc2, uint* nDocs2, uint* V2_aux){
	if (l > r)
		return;

	uint doc;
	ulong med, i, min, t=lev, lq=lcur, rq=rcur;
	if (l==r) min = i = l;
	else min = i = Range[lev].rmqB->queryRMQ(l,r);

	while (t < h-1){
		med = lq + (rq-lq)/2;
		if(Range[t].B_rrr[i]){			// go down to the right...
			i = med+Range[t].B_rank1.rank(i+1);
			i -= Range[t].B_rank1.rank(lq);
			lq = med+1;
		}else{								// go down to the left...
		i = lq+Range[t].B_rank0.rank(i+1)-1;
		i -= Range[t].B_rank0.rank(lq);
			rq = med;
		}
		t++;
	}
	// now, i is the preorder in lz of the occurrence type 2 to report...
	if (Range[t].B_rrr[i])
		i = rq;
	else
		i = lq;
	doc = getNum64(DocLZ, i*lgD, lgD);
	//cout << " to report Occ 2 ? -> doc " << doc << endl;

	if (V2_aux[doc] != mark){
		if (V2_aux[doc] == 0){
			occ2[*nDocs2] = doc;
			(*nDocs2)++;
			//cout << " Occ 2 for doc " << doc << endl;
		}//else cout << " DESCARTADo C Occ 2 -> doc " << doc << endl;
		V2_aux[doc] = mark;
		if (min > 0)
			reportDocOcc2_test(lev, lcur, rcur, l, min-1, occ2, nOcc2, nDocs2, V2_aux);
		if (min < nod-1)
			reportDocOcc2_test(lev, lcur, rcur, min+1, r, occ2, nOcc2, nDocs2, V2_aux);
	}//else cout << " DESCARTADOOO V2_aux[doc] != mark, doc " << doc << endl;
}

// [lcur, rcur] is the current bitstring in the level 'lev' (B[lcur, rcur]). [lobj, robj] is the target interval to search in LZTrie, and [lrmq, rrmq] is the interval to query by RMQ in RevTrie
void LZDocList64::searchIntervalInRangeOcc(uint lev, ulong lcur, ulong rcur, ulong lobj, ulong robj, ulong lrmq, ulong rrmq, uint* occ, uint* nDocs, uint *nDocs2){
	if (lev == h-1){	// this is a leaf !!
		if (lrmq < rrmq){
			if (lobj <= lcur){
				uint doc = getNum64(DocLZ, lcur*lgD, lgD);
				if (V[doc] == 0){
					occ[*nDocs] = doc;
					(*nDocs)++;
					(*nDocs2)++;
					//cout << " hoja -> doc " << doc << endl;
				}
				V[doc] = mark;
			}
			if (rcur <= robj){
				uint doc = getNum64(DocLZ, rcur*lgD, lgD);
				if (V[doc] == 0){
					occ[*nDocs] = doc;
					(*nDocs)++;
					(*nDocs2)++;
					//cout << " hoja -> doc " << doc << endl;
				}
				V[doc] = mark;
			}
		}else{ // then lcur == rcur
			if (Range[lev].B_rrr[lrmq] == 0 && lobj <= lcur){
				uint doc = getNum64(DocLZ, lcur*lgD, lgD);
				if (V[doc] == 0){
					occ[*nDocs] = doc;
					(*nDocs)++;
					(*nDocs2)++;
					//cout << " hoja -> doc " << doc << endl;
				}
				V[doc] = mark;
			}else{
				if (Range[lev].B_rrr[lrmq] && rcur <= robj){
					uint doc = getNum64(DocLZ, rcur*lgD, lgD);
					if (V[doc] == 0){
						occ[*nDocs] = doc;
						(*nDocs)++;
						(*nDocs2)++;
						//cout << " hoja -> doc " << doc << endl;
					}
					V[doc] = mark;
				}
			}
		}
	}else{
		if(lobj <= lcur && rcur <= robj){			// we consider the complete node --> apply RMQ...  the node is totally inside of objective
			if (lev<levRMQ){	// is there in this level a RMQ structure ?
				mark++;
				reportDocOcc2Occ(lev, lcur, rcur, lrmq, rrmq, occ, nDocs, nDocs2);
			}else
				reportDocOcc2_levRangeOcc(lev, lcur, rcur, lrmq, rrmq, lobj, robj, occ, nDocs, nDocs2);
		}else{
			ulong med = lcur + ((rcur-lcur)>>1);

			if(robj <= med){
				// we go down only to the left...
				ulong befL0=0, befR0=0, aux;
				aux = Range[lev].B_rank0.rank(lcur);
				befL0 = Range[lev].B_rank0.rank(lrmq) - aux;
				befR0 = Range[lev].B_rank0.rank(rrmq+1) - aux;
				if (lev<h && befR0>befL0)
					searchIntervalInRangeOcc(lev+1, lcur, med, lobj, robj, lcur+befL0, lcur+befR0-1, occ, nDocs, nDocs2);
			}else{
				if(lobj > med){
					// we go down only to the right...
					ulong befL1=0, befR1=0, aux;
					aux = Range[lev].B_rank1.rank(lcur);
					befL1 = Range[lev].B_rank1.rank(lrmq) - aux;
					befR1 = Range[lev].B_rank1.rank(rrmq+1) - aux;
					if (befR1>befL1)
						searchIntervalInRangeOcc(lev+1, med+1, rcur, lobj, robj, med+befL1+1, med+befR1, occ, nDocs, nDocs2);
				}else{
					// first, we go down to the left...
					ulong befL0=0, befR0=0, aux;
					aux = Range[lev].B_rank0.rank(lcur);
					befL0 = Range[lev].B_rank0.rank(lrmq) - aux;
					befR0 = Range[lev].B_rank0.rank(rrmq+1) - aux;
					if (lev<h && befR0>befL0)
						searchIntervalInRangeOcc(lev+1, lcur, med, lobj, robj, lcur+befL0, lcur+befR0-1, occ, nDocs, nDocs2);

					// second, we go down to the right...
					ulong befL1=0, befR1=0;
					aux = Range[lev].B_rank1.rank(lcur);
					befL1 = Range[lev].B_rank1.rank(lrmq) - aux;
					befR1 = Range[lev].B_rank1.rank(rrmq+1) - aux;
					if (lev<h && befR1>befL1)
						searchIntervalInRangeOcc(lev+1, med+1, rcur, lobj, robj, med+befL1+1, med+befR1, occ, nDocs, nDocs2);
				}
			}
		}
	}
}

// [lcur, rcur] is the current bitstring in the level 'lev' (B[lcur, rcur]). [lobj, robj] is the target interval to search in LZTrie, and [lrmq, rrmq] is the interval to query by RMQ in RevTrie
void LZDocList64::searchIntervalInRange_test(uint lev, ulong lcur, ulong rcur, ulong lobj, ulong robj, ulong lrmq, ulong rrmq, uint* occ2, ulong* nOcc2, uint* nDocs2, uint* auxV2){
	if (lev == h-1){	// this is a leaf !!
		if (lrmq < rrmq){
			if (lobj <= lcur){
				uint doc = getNum64(DocLZ, lcur*lgD, lgD);
				(*nOcc2)++;
				if (auxV2[doc] == 0){
					occ2[*nDocs2] = doc;
					(*nDocs2)++;
					//cout << " hoja -> doc " << doc << ", *nOcc2 = " << *nOcc2 << endl;
					auxV2[doc] = mark;
				}//else cout << "DESCARTADA hoja -> doc " << doc << ", *nOcc2 = " << *nOcc2 << endl;
			}
			if (rcur <= robj){
				uint doc = getNum64(DocLZ, rcur*lgD, lgD);
				(*nOcc2)++;
				if (auxV2[doc] == 0){
					occ2[*nDocs2] = doc;
					(*nDocs2)++;
					//cout << " hoja -> doc " << doc << ", *nOcc2 = " << *nOcc2 << endl;
					auxV2[doc] = mark;
				}//else cout << "DESCARTADA hoja -> doc " << doc << ", *nOcc2 = " << *nOcc2 << endl;
			}
		}else{ // then lcur == rcur
			if (Range[lev].B_rrr[lrmq] == 0 && lobj <= lcur){
				uint doc = getNum64(DocLZ, lcur*lgD, lgD);
				(*nOcc2)++;
				if (auxV2[doc] == 0){
					occ2[*nDocs2] = doc;
					(*nDocs2)++;
					//cout << " hoja -> doc " << doc << ", *nOcc2 = " << *nOcc2 << endl;
					auxV2[doc] = mark;
				}//else cout << "DESCARTADA hoja -> doc " << doc << ", *nOcc2 = " << *nOcc2 << endl;
			}else{
				if (Range[lev].B_rrr[lrmq] && rcur <= robj){
					uint doc = getNum64(DocLZ, rcur*lgD, lgD);
					(*nOcc2)++;
					if (auxV2[doc] == 0){
						occ2[*nDocs2] = doc;
						(*nDocs2)++;
						//cout << " hoja -> doc " << doc << ", *nOcc2 = " << *nOcc2 << endl;
						auxV2[doc] = mark;
					}//else cout << "DESCARTADA hoja -> doc " << doc << ", *nOcc2 = " << *nOcc2 << endl;
				}
			}
		}
	}else{
		if(lobj <= lcur && rcur <= robj){			// we consider the complete node --> apply RMQ...  the node is totally inside of objective
			(*nOcc2)+= rrmq-lrmq+1;

			if (lev<levRMQ){	// is there in this level a RMQ structure ?
				mark++;
				//cout << " Apply reportDocOcc2_test in lev " << lev << ", +" << rrmq-lrmq+1 << ", *nOcc2 = " << *nOcc2 <<  ", MARK = " << mark << endl;

				reportDocOcc2_test(lev, lcur, rcur, lrmq, rrmq, occ2, nOcc2, nDocs2, auxV2);
			}else{
				//cout << " reportDocOcc2_levRange_test in lev " << lev << ", +" << rrmq-lrmq+1 << ", *nOcc2 = " << *nOcc2 << endl;
				reportDocOcc2_levRange_test(lev, lcur, rcur, lrmq, rrmq, lobj, robj, occ2, nOcc2, nDocs2, auxV2);
			}
		}else{
			ulong med = lcur + ((rcur-lcur)>>1);
			if(robj <= med){
				// we go down only to the left...
				ulong befL0=0, befR0=0, aux;
				aux = Range[lev].B_rank0.rank(lcur);//(lcur-1);
				befL0 = Range[lev].B_rank0.rank(lrmq) - aux;//(lrmq-1) - aux;
				befR0 = Range[lev].B_rank0.rank(rrmq+1) - aux;//(rrmq) - aux;
				if (lev<h && befR0>befL0)
					searchIntervalInRange_test(lev+1, lcur, med, lobj, robj, lcur+befL0, lcur+befR0-1, occ2, nOcc2, nDocs2, auxV2);
			}else{
				if(lobj > med){
					// we go down only to the right...
					ulong befL1=0, befR1=0, aux;
					aux = Range[lev].B_rank1.rank(lcur);//(lcur-1);
					befL1 = Range[lev].B_rank1.rank(lrmq) - aux;//(lrmq-1) - aux;
					befR1 = Range[lev].B_rank1.rank(rrmq+1) - aux;//(rrmq) - aux;
					if (befR1>befL1)
						searchIntervalInRange_test(lev+1, med+1, rcur, lobj, robj, med+befL1+1, med+befR1, occ2, nOcc2, nDocs2, auxV2);
				}else{
					// first, we go down to the left...
					ulong befL0=0, befR0=0, aux;
					aux = Range[lev].B_rank0.rank(lcur);//(lcur-1);
					befL0 = Range[lev].B_rank0.rank(lrmq) - aux;//(lrmq-1) - aux;
					befR0 = Range[lev].B_rank0.rank(rrmq+1) - aux;//(rrmq) - aux;
					if (lev<h && befR0>befL0)
						searchIntervalInRange_test(lev+1, lcur, med, lobj, robj, lcur+befL0, lcur+befR0-1, occ2, nOcc2, nDocs2, auxV2);

					// second, we go down to the right...
					ulong befL1=0, befR1=0;
					aux = Range[lev].B_rank1.rank(lcur);//(lcur-1);
					befL1 = Range[lev].B_rank1.rank(lrmq) - aux;//(lrmq-1) - aux;
					befR1 = Range[lev].B_rank1.rank(rrmq+1) - aux;//(rrmq) - aux;
					if (lev<h && befR1>befL1)
						searchIntervalInRange_test(lev+1, med+1, rcur, lobj, robj, med+befL1+1, med+befR1, occ2, nOcc2, nDocs2, auxV2);
				}
			}
		}
	}
}


/*
	cout << " n = " << n << endl;
	cout << " h = " << h << endl;
	cout << " levRMQ = " << levRMQ << endl;
	cout << " nDocs = " << nDocs << endl;
	cout << " lgD = " << lgD << endl;
	cout << " cutDoc = " << cutDoc << endl;
	cout << " nod = " << nod << endl;
	cout << " lgNod = " << lgNod << endl;
	cout << " nF = " << nF << endl;
	cout << " nEF = " << nEF << endl;
	cout << " nEFHead = " << nEFHead << endl;
	cout << " nodRev = " << nodRev << endl;
	cout << " nLz = " << nLz << endl;
	cout << " nRev = " << nRev << endl;
 */

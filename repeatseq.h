//*----------------------------------------------------------------------------------------*
//*RepeatSeq is available through the Virginia Tech non-commerical license.                *
//*For more details on the license and use, see license.txt included in this distribution. *
//*----------------------------------------------------------------------------------------*

#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//from bamtools:
#include <api/BamReader.h>
#include <api/BamWriter.h>

//from fastahack:
#include "fastahack/Fasta.h"
#include "fastahack/disorder.h"
#include "fastahack/split.h"

using namespace std;
using namespace BamTools;       // use BamTools classes/methods

//STRUCTURES:
struct Sequences {
	bool insertions;
	string preSeq;
	string alignedSeq;
	string postSeq;

	Sequences();
	Sequences(string, string, string, bool);
};

struct STRING_GT {
	string print;
	Sequences reads;
	int GT;
	bool paired;
	bool reverse;
	int MapQ;
	int minFlank;
	double avgBQ;

	STRING_GT(string, Sequences, int, bool, int, int, bool, double);
	STRING_GT();
	bool operator<(const STRING_GT &other) const; //for sorting
};

struct SETTINGS_FILTERS {
	int MAX_READ_SIZE;
	int LR_CHARS_TO_PRINT;
	int mode;
	bool emitAll;
	bool multi;
	bool properlyPaired;
	bool makeRepeatseqFile;
	bool makeCallsFile;
	int readLengthMin;
	int readLengthMax;
	int consLeftFlank;
	int consRightFlank;
	int MapQuality;
	string paramString;
	
	SETTINGS_FILTERS(){
		LR_CHARS_TO_PRINT = 8;
		MAX_READ_SIZE = 200;
		mode = 2;
		emitAll = false;
		multi = false;
		properlyPaired = false;
		makeRepeatseqFile = false;
		makeCallsFile = false;
		readLengthMin = 0;
		readLengthMax = 0;
		consLeftFlank = 3;
		consRightFlank = 3;
		MapQuality = 0;
		paramString = "";
	}
};

//to be used for printing GT: information to header line:
struct GT {
	int readlength;
	int occurrences;
	int reverse;
	int avgMinFlank;
	double avgBQ;

	GT(int rl, int oc, int rev, int minF, double avg);
    static bool sortByReadLength(const GT & a, const GT & b) { return (a.readlength > b.readlength); }
};

//counter struct is used in array for table files:
struct counter {
	int numGT;              //number of repeats that have a GT
	int numRepeats;         //number of total repeats
	double numRepeats2;     //number of repeats with 2 or more reads present
	double tallyC;          //a running sum of C:'s to be later divided by numRepeats
	
	counter();
};

//structure to be passed to VCF-writing function:
struct VCF_INFO {
	string chr;
	int start;
	string unit;
	int length;
	int purity;
	int depth;
	bool emitAll;
};

//class for parsing region argument:
class Region {
public:
	string startSeq;
	int startPos;
	int stopPos;
	
	Region(string& region);
	int length(void);
};

//function declarations:
float fact(int);
double retSumFactOverIndFact(int, int, int);
string getVCF(vector<string>, string, string, int, char, VCF_INFO, map<pair<int,int>,double> &);
double PhredToFloat(char);
string setToCD (string);
bool fileCheck(string);
void buildFastaIndex(string);
void printHeader(ofstream&);
void parseSettings(char**, int, SETTINGS_FILTERS&, string&, string&, string&);
void printArguments();
vector<int> printGenoPerc(vector<GT>, int, int, double&, int, map<pair<int, int>, double> &);
bool fileCheck(string);
void buildFastaIndex(string);
void print_output(string, FastaReference*, ofstream&, ofstream&, ofstream&,  const SETTINGS_FILTERS&, BamReader&);

inline bool vectorGTsort(GT a, GT b) { return (a.occurrences > b.occurrences); }


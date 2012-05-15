// Custom structure & class function definitions (Version7a)
//

#include "repeatseq.h"

STRING_GT::STRING_GT(string a, Sequences b, int c, bool d, int e, int f, bool g, double h){
	GT = c;
	reads = Sequences(b);
	print = a;
	paired = d;
	MapQ = e;
	minFlank = f;
	reverse = g;
	avgBQ = h;
}

STRING_GT::STRING_GT(){
	print = "";
	GT = 0;
	paired = 0;
	MapQ = 0;
	minFlank = 0;
	reverse = 0;
}

//overloaded for sorting
bool STRING_GT::operator<(const STRING_GT &other) const {
	if (minFlank >= other.minFlank) return false;
	else return true;
}

Sequences::Sequences(string a, string b, string c, bool d){
	preSeq = a;
	alignedSeq = b;
	postSeq = c;
	insertions = d;
}

Sequences::Sequences(){
	preSeq = "";
	alignedSeq = "";
	postSeq = "";
	insertions = 0;
}

GT::GT(int rl, int oc, int rev, int minF, double avgbq){
	avgMinFlank = minF;
	readlength = rl;
	occurrences = oc;
	reverse = rev;
	avgBQ = avgbq;
};

counter::counter(){
	numGT = 0;
	numRepeats = 0;
	tallyC = 0;
	numRepeats2 = 0;
}

Region::Region(string& region) {
	startPos = -1;
	stopPos = -1;
	size_t foundFirstColon = region.find(":");
	// we only have a single string, use the whole sequence as the target
	if (foundFirstColon == string::npos) {
		startSeq = region;
	}
	else {
		startSeq = region.substr(0, foundFirstColon);
		size_t foundRangeDots = region.find("-", foundFirstColon);
		if (foundRangeDots == string::npos) {
			startPos = atoi(region.substr(foundFirstColon + 1).c_str());
			stopPos = startPos; // just print one base if we don't give an end
		} else {
			startPos = atoi(region.substr(foundFirstColon + 1, foundRangeDots - foundRangeDots - 1).c_str());
			stopPos = atoi(region.substr(foundRangeDots + 1).c_str()); // to the start of this chromosome
		}
	}
}

int Region::length(void) {
	if (stopPos > 0) {
		return stopPos - startPos + 1;
	} else {
		return 1;
	}
}

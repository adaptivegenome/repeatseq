//*----------------------------------------------------------------------------------------*
//*RepeatSeq is available through the Virginia Tech non-commerical license.                *
//*For more details on the license and use, see license.txt included in this distribution. *
//*----------------------------------------------------------------------------------------*

/*
 repeatseq.cpp - main source code file for RepeatSeq
 
 See "repeatseq.h" for function & custom data structure declarations
 
 This .cpp contains functions: 
  (1) main() - Parse command line options, iterate line by line through the TRF file (calling 
               print_output() on each region).
               
  (2) print_output() - This function is called for each repeat in the repeat file, and 
		       handles the calling of other functions to determine genotype and print 
                       data to files.

  (3) parseCigar() - Uses CIGAR sequence to align read with reference sequence.
 
  (4) printGenoPerc() - perform statistical analysis to determine most likely genotype and its 
                        likelihood.
 
  (5) getVCF() - print variant record to VCF file.
*/

#include "repeatseq.h"
#include <algorithm>
#include <pthread.h>
#include <unistd.h>

double log_factorial[100000] = {};
string VERSION = "0.5.3";

typedef struct worker_data {
    worker_data(const SETTINGS_FILTERS & settings, const vector<string> & regions)
    : settings(settings)
    , regions(regions)
    {}
    FastaReference * fr;
    ofstream vcfFile, oFile, callsFile;
    string vcfFilename, oFilename, callsFilename;
    const SETTINGS_FILTERS & settings;
    const vector<string> & regions;
    size_t region_start, region_stop;
    pthread_t thread;
    BamReader reader;
} worker_data_t;

void * worker_thread(void * pdata) {
    worker_data_t & worker_data = *((worker_data_t *) pdata);
    
    worker_data.vcfFile.open(worker_data.vcfFilename.c_str());
    if (worker_data.settings.makeRepeatseqFile)
        worker_data.oFile.open(worker_data.oFilename.c_str());
    if (worker_data.settings.makeCallsFile)
        worker_data.callsFile.open(worker_data.callsFilename.c_str());
    
    for(size_t i = worker_data.region_start; i != worker_data.region_stop; i++)
        print_output(worker_data.regions[i], worker_data.fr, worker_data.vcfFile, worker_data.oFile, worker_data.callsFile, worker_data.settings, worker_data.reader);
    
    worker_data.vcfFile.close();
    if (worker_data.settings.makeRepeatseqFile)
        worker_data.oFile.close();
    if (worker_data.settings.makeCallsFile)
        worker_data.callsFile.close();

    return NULL;
}

int main(int argc, char* argv[]){
    
    ofstream oFile, callsFile, vcfFile;
	try{
		SETTINGS_FILTERS settings;	
		srand( time(NULL) );
		string bam_file = "", fasta_file = "", position_file = "", region;
		
		//load log_factorial vector
		for (int i=1,val=0 ; i < 100000; ++i){ 
			val += log(i);
			log_factorial[i] = val;
		}

		//parse arguments, store in settings:
		parseSettings(argv, argc, settings, bam_file, fasta_file, position_file);
		if (bam_file == "") { throw "NO BAM FILE"; }
		if (fasta_file == "") { throw "NO FASTA FILE"; }
		if (position_file == "") { throw "NO POSITION FILE"; }
		
		//create index filepaths & output filepaths (ensuring output is to current directory):
		string fasta_index_file = fasta_file + ".fai";
		string bam_index_file = bam_file + ".bai";
		string output_filename = setToCD(bam_file + settings.paramString + ".repeatseq");
		string calls_filename = setToCD(bam_file + settings.paramString + ".calls");
		string vcf_filename = setToCD(bam_file + settings.paramString + ".vcf");
		
		//open FastaReference object (creating fasta index file if needed):
		if (!fileCheck(fasta_index_file)) {
			cout <<  "Fasta index file not found, creating...";
			buildFastaIndex(fasta_file);
		}

		//open input & output filestreams:
		if (settings.makeRepeatseqFile){ oFile.open(output_filename.c_str()); }
	 	if (settings.makeCallsFile){ callsFile.open(calls_filename.c_str()); }
		vcfFile.open(vcf_filename.c_str());
		ifstream range_file(position_file.c_str());
		if (!range_file.is_open()) { throw "Unable to open input range file."; }
		
		//print VCF header information:
		printHeader(vcfFile);
		
        //read in the region file
        vector<string> regions;
		while(getline(range_file,region))
            regions.push_back(region);
        
        long num_threads = sysconf(_SC_NPROCESSORS_ONLN);
        vector<worker_data_t *> thread_worker_data;
        
        //set up threads to actually print the output
        for(int thread = 0; thread != num_threads; thread++) {
            thread_worker_data.push_back(new worker_data_t(settings, regions));
            worker_data_t & data = *(thread_worker_data.back());
            if (!data.reader.Open(bam_file)){ throw "Could not open BAM file.."; }
            if (!data.reader.OpenIndex(bam_index_file)){ throw "Could not open BAM index file.."; }

            data.fr = new FastaReference();
            data.fr->open(fasta_file);

            char filename_buffer[32];
            sprintf(filename_buffer, "repeatseq_vcf_temp_%03d.vcf", thread);
            data.vcfFilename = filename_buffer;
            sprintf(filename_buffer, "repeatseq_o_temp_%03d.o", thread);
            data.oFilename = filename_buffer;
            sprintf(filename_buffer, "repeatseq_calls_temp_%03d.calls", thread);
            data.callsFilename = filename_buffer;
            
            data.region_start = thread * (regions.size() / num_threads);
            if(thread == num_threads - 1)
                data.region_stop = regions.size();
            else
                data.region_stop = (thread+1) * (regions.size() / num_threads);
        }
        
        //start worker threads
        for(int thread = 0; thread != num_threads; thread++) {
            if(0 != pthread_create(&thread_worker_data[thread]->thread, NULL, worker_thread, thread_worker_data[thread]))
                perror("Error starting worker thread");
        }
        
        //wait for all workers to finish
        for(int thread = 0; thread != num_threads; thread++) {
            if(0 != pthread_join(thread_worker_data[thread]->thread, NULL))
                perror("Error closing worker thread");
        }
        
        //consolidate results from the worker threads
        for(int thread = 0; thread != num_threads; thread++) {
            worker_data_t & data = *thread_worker_data[thread];

            ifstream in_vcf(data.vcfFilename.c_str(), ios::binary);
            vcfFile << in_vcf.rdbuf();
            remove(data.vcfFilename.c_str());

            if (settings.makeRepeatseqFile) {
                ifstream in_o(data.oFilename.c_str(), ios::binary);
                oFile << in_o.rdbuf();
                remove(data.oFilename.c_str());
            }

            if (settings.makeCallsFile) {
                ifstream in_calls(data.callsFilename.c_str(), ios::binary);
                callsFile << in_calls.rdbuf();
                remove(data.callsFilename.c_str());
            }
        }
	}
	catch(const char* exOutput) {
		cout << endl << exOutput << endl;
		printArguments();
		return 0;
	}	
}

inline string parseCigar(stringstream &cigarSeq, string &alignedSeq, string &QS, vector<string> & insertions, int alignStart, int refStart, int LR_CHARS_TO_PRINT, double &avgBQ){
	int reserveSize = alignedSeq.length() + 500;

	//reserve sufficient space (so iterators remain valid)
	alignedSeq.reserve(reserveSize);
	string tempInsertions = "";
	tempInsertions.reserve(reserveSize);
	
	//iterators & other variables
	string::iterator it=alignedSeq.begin();
	string::iterator START;
	char cigChar;
	int cigLength;	
	bool STARTset = 0;
	int posLeft = refStart - alignStart;
	int posLeftINS = refStart - alignStart - LR_CHARS_TO_PRINT;
	bool firstRun = true;
	
	//determine average base quality:
	avgBQ = 0;
	for (int i=0; i<QS.length(); ++i){ avgBQ += PhredToFloat(QS[i]); }
	avgBQ /= QS.length();

	cigarSeq.clear();
	while (!cigarSeq.eof()) {
		//parse:
		cigLength = -1;
		cigarSeq >> cigLength;
		cigarSeq >> cigChar;
		
		//Perform operations on aligned seq:
		switch(cigChar) {
			case 'M':                   //MATCH to the reference
				for (int i = cigLength; i>0; i--) {
					if (posLeft>0) {
						posLeft--;
						posLeftINS--;
					}
					else if (!STARTset) {
						START = it;
						STARTset = 1;
					}
					it++;
				}
				
				break;
				
			case 'I':                  //INSERTION to the reference
				tempInsertions = "";
				*(it-1) += 32;	//convert previous letter to lower case (to mark the following insertion)
				
				for (int i = cigLength; i>0; i--) {
					tempInsertions += *it + 1;
					*it = 'd';         //convert to d's to remove later
					++it;
				}
				if (posLeftINS <= 0) { insertions.push_back(tempInsertions); }
				
				break;
				
			case 'D':                       //DELETION from the reference
				for (int i = cigLength; i>0; i--) {
					alignedSeq.insert(it, 1, '-');	//inserts - into the aligned seq
					posLeft--; 
					posLeftINS--;
					if (posLeft < 0 && !STARTset) {
						START = it;
						STARTset = 1;
					}
					++it;
				}
				break;
				
			case 'N':       //SKIPPED region from the reference
				return "";	//fail the read (return null string)
				
			case 'S':                       //SOFT CLIP on the read (clipped sequence present in <seq>)
				if (firstRun && !STARTset) posLeft+=cigLength;
				else if (firstRun) START -= cigLength;
				
				for (int i = cigLength; i>0; i--) {
					if (posLeft>0) {
						posLeft--;
						posLeftINS--;
					}
					else if (!STARTset) {
						START = it;
						STARTset = 1;
					}
					
					*it = 'S';				//mark as soft-clipped
					++it;
				}
				break;
				
			case 'H':   //HARD CLIP on the read (clipped sequence NOT present in <seq>)
				break;
				
			case 'P':   //PADDING (silent deletion from the padded reference sequence)
				if (posLeft>0) {
					posLeft--;
					posLeftINS--;
				}
				else if (!STARTset) {
					START = it;
					STARTset = 1;
				}
				
				it += cigLength;
				break;
		}       //end case
		firstRun = 0;
	}
	
	int offset = alignStart - refStart;
	if (!STARTset) {
		START = alignedSeq.begin();
		while ((offset--) > 0) { alignedSeq.insert(alignedSeq.begin(),1,'x'); }
	}
	int numD = 0;
	for (string::iterator ii = START; ii > START - LR_CHARS_TO_PRINT && ii >= alignedSeq.begin(); --ii) {
		if (*ii == 'd') numD++;
	}
	
	string temp = "";
	temp.reserve(500);
	
	string::iterator ii = START;
	for (int i = 0; i < numD + LR_CHARS_TO_PRINT; ++i) {
		if (ii > alignedSeq.begin()) {
			--ii;
			temp.insert(0,1,*ii);
		}
		else {
			temp.insert(0,1,'x');
		}
	}
	
	ii = START;
	for (int i = 0; i < alignStart - refStart; ++i){ temp += 'x'; }
	while (ii < alignedSeq.end()) { temp += *(ii++); }
	
	return temp; //return modified string
}

inline void print_output(string region,FastaReference* fr, ofstream &vcf,  ofstream &oFile, ofstream &callsFile, const SETTINGS_FILTERS &settings, BamReader & reader){
	
	vector<string> insertions;
	string sequence;                // holds reference sequence
	string secondColumn;            // text string to the right of tab
	int unitLength;
	double purity;
	
	// parse region argument:
	secondColumn = region.substr(region.find('\t',0)+1,-1);
	if (secondColumn == "") cout << "missing information after the tab in region file for " << region << ".\ncontinuing..." << endl;
	region = region.substr(0,region.find('\t',0));          //erases all of region string after tab
	
	// parse secondColumn:
	if (int(secondColumn.find('_',0)) == -1) {
		cout << "improper second column found for " << region << ".\ncontinuing with next region..." << endl;
		return;
	}
	unitLength = atoi(secondColumn.substr(0,secondColumn.find('_',0)).c_str());
	string UnitSeq = secondColumn.substr(secondColumn.rfind('_')+1);	
	
	int pos = 0;
	for (int i = 0; i < 3; ++i) pos = secondColumn.find('_',pos + 1);
	++pos; //increment past fourth '_'
	purity = atof(secondColumn.substr(pos,secondColumn.find('_',pos)).c_str());
	
	Region target(region);
	if (target.startPos > target.stopPos) throw "Invalid input file...";
	
	//ensure target doesn't overrun end of chromosome
	if (target.startPos+target.length() > fr->sequenceLength(target.startSeq)+1) throw "Target range is outside of chromosome.\n exiting..";
	
	//if asked to print entire sequence:
	if (target.startPos == -1) sequence = fr->getSequence(target.startSeq);
	
	//when start position is exactly at the beginning of the chromosome:
	else if (target.startPos == 1)
		sequence = " "
		+ fr->getSubSequence(target.startSeq, target.startPos - 1, target.length())
		+ " "
		+ fr->getSubSequence(target.startSeq, target.startPos - 1 + target.length(), settings.LR_CHARS_TO_PRINT);
	
	//when start position is within 20 of beginning of chromosome:
	else if (target.startPos < 1 + settings.LR_CHARS_TO_PRINT)
		sequence = fr->getSubSequence(target.startSeq, 0, target.startPos - 1)
		+ " "
		+ fr->getSubSequence(target.startSeq, target.startPos - 1, target.length())
		+ " "
		+ fr->getSubSequence(target.startSeq, target.startPos - 1 + target.length(), settings.LR_CHARS_TO_PRINT);
	
	
	//when end position is exactly at the end of chromosome:
	else if (target.startPos+target.length() ==  fr->sequenceLength(target.startSeq)+1)
		sequence = fr->getSubSequence(target.startSeq, target.startPos - 1 - settings.LR_CHARS_TO_PRINT, settings.LR_CHARS_TO_PRINT)
		+ " "
		+ fr->getSubSequence(target.startSeq, target.startPos - 1, target.length())
		+ " ";
	
	//when end position is within 20 of end of chromosome:
	else if (target.startPos+target.length()+settings.LR_CHARS_TO_PRINT > fr->sequenceLength(target.startSeq)+1)
		sequence = fr->getSubSequence(target.startSeq, target.startPos - 1 - settings.LR_CHARS_TO_PRINT, settings.LR_CHARS_TO_PRINT)
		+ " "
		+ fr->getSubSequence(target.startSeq, target.startPos - 1, target.length())
		+ " "
		+ fr->getSubSequence(target.startSeq, target.startPos - 1 + target.length(), fr->sequenceLength(target.startSeq)-target.startPos-target.length()+1);
	
	//all other cases:
	else sequence = fr->getSubSequence(target.startSeq, target.startPos - 1 - settings.LR_CHARS_TO_PRINT, settings.LR_CHARS_TO_PRINT)
		+ " "
		+ fr->getSubSequence(target.startSeq, target.startPos - 1, target.length())
		+ " "
		+ fr->getSubSequence(target.startSeq, target.startPos - 1 + target.length(), settings.LR_CHARS_TO_PRINT);
	
	int firstSpace = sequence.find(' ',0);
	int secondSpace = sequence.find(' ',firstSpace+1);
	
	string leftReference, centerReference, rightReference;
	if (firstSpace != 0) leftReference = sequence.substr(0,firstSpace);
	else leftReference = "";
	centerReference = sequence.substr(firstSpace+1,secondSpace-firstSpace-1);
	if (secondSpace != -1) rightReference = sequence.substr(secondSpace+1,-1);
	else rightReference = "";
	
	// ensure reference is all caps (for matching purposes):
	std::transform(leftReference.begin(), leftReference.end(), leftReference.begin(), ::toupper);	
	std::transform(centerReference.begin(), centerReference.end(), centerReference.begin(), ::toupper);	
	std::transform(rightReference.begin(), rightReference.end(), rightReference.begin(), ::toupper);	
	
	// define our region of interest:
	// debug-cout << "region: " << target.startSeq << ":" << target.startPos-1 << "-" << target.stopPos-1 << endl;
	BamRegion bamRegion(reader.GetReferenceID(target.startSeq), target.startPos - 1,reader.GetReferenceID(target.startSeq), target.stopPos - 1);
	reader.SetRegion(bamRegion);
	
	// prep for getting alignment info
	BamAlignment al;
	stringstream ssPrint;                   //where data to print will be stored
	string PreAlignedPost = "";             //contains all 3 strings to be printed
	double concordance = 0;
	int totalOccurrences = 0;
	int majGT = 0;
	int occurMajGT = 0;
	int depth = 0;
	int numReads = 0;
	int numStars = 0;
	
	vector<GT> vectorGT;
	vectorGT.reserve(100);
	vector<STRING_GT> toPrint;
	toPrint.reserve(100);
	
	string vcfPrint;
	
	//cout << "trying " << target.startSeq << ":" << target.startPos - 1 << "-" << target.stopPos - 1 << endl;
	// iterate through alignments in this region,
	while (reader.GetNextAlignment(al)) {
		//cout << " found\n";
		insertions.clear();
		ssPrint.str("");
		stringstream cigarSeq;
		int gtBonus = 0;
		
		if (al.CigarData.begin()==al.CigarData.end()) {
			numStars++;
			continue;
			//if CIGAR is not there, it's * case..
			//so increment numStars and get next alignment
		}
		
		//load cigarSeq
		for ( vector<BamTools::CigarOp>::const_iterator it=al.CigarData.begin(); it < al.CigarData.end(); it++ ) {
			cigarSeq << it->Length;
			cigarSeq << it->Type;
		}
		
		//run parseCigar:
		double avgBQ;
		PreAlignedPost = parseCigar(cigarSeq, al.QueryBases, al.Qualities, insertions, al.Position + 1, target.startPos, settings.LR_CHARS_TO_PRINT, avgBQ);
		if (PreAlignedPost == ""){ 
			//If an 'N' or other problem was found
			cout << "N found-- Possible Error!\n";
			continue; 
		} 
		
		//adjust for d's
		for (int a = PreAlignedPost.find('d',0); a!=-1; a=PreAlignedPost.find('d',0)) {
			if ( (a + 1) > settings.LR_CHARS_TO_PRINT && (a + 1) < settings.LR_CHARS_TO_PRINT + target.length()) gtBonus+=1;
			PreAlignedPost.erase(a,1);
		}
		
		//set strings to print based off of value input
		string PreSeq, AlignedSeq, PostSeq;
		
		//if there's not enough characters to make it through PreSeq, skip read
		if (PreAlignedPost.length() < settings.LR_CHARS_TO_PRINT+1) continue;
		
		//Split PreAlignedPost into 3 substrings
		PreSeq = PreAlignedPost.substr(0,settings.LR_CHARS_TO_PRINT);
		AlignedSeq = PreAlignedPost.substr(settings.LR_CHARS_TO_PRINT, target.length());
		if (AlignedSeq.length() < target.length()) AlignedSeq.resize(target.length(),'x');
		else PostSeq = PreAlignedPost.substr(settings.LR_CHARS_TO_PRINT + target.length(), settings.LR_CHARS_TO_PRINT);
		PostSeq.resize(settings.LR_CHARS_TO_PRINT,'x');
		
		if (AlignedSeq[target.length()/2] != 'x') ++depth;      //increment depth (if middle character is NOT an x)
		int numMatchesL = 0, numMatchesR = 0;
		int minflank = 0;

		// if first and last characters of sequence range are present in read, print it's information:
		if (AlignedSeq[0] != ' ' && AlignedSeq[0] != 'x' && AlignedSeq[0] != 'X' && AlignedSeq[0] != 'S') {
			if (AlignedSeq[AlignedSeq.length()-1]!= 'x' && AlignedSeq[AlignedSeq.length()-1]!= ' ' && AlignedSeq[AlignedSeq.length()-1]!='X' && AlignedSeq[AlignedSeq.length()-1] != 'S') {
				string toprintPre = string(PreSeq);
				string toprintAligned = string(AlignedSeq);
				string toprintPost = string(PostSeq);
				
				bool hasinsertions = (! insertions.empty());
				if (hasinsertions){
					//PROCESS SEQUENCE:
					//put insertions back in pre-sequence (as lower case) here
					for (int i = 0; i < toprintPre.length();){
						if (toprintPre[i] > 96 && toprintPre[i] != 'x'){	//is lowercase
							toprintPre[i++] -= 32;							//convert to uppercase
							if (i == toprintPre.length()) toprintAligned = insertions.front() + toprintAligned;
							else toprintPre.insert(i,insertions.front());
							insertions.erase(insertions.begin());
						}
						else ++i;
					}
					//put insertions back in Aligned-sequence (as lower case) here
					for (int i = 0; i < toprintAligned.length();){
						if (toprintAligned[i] > 96 && toprintAligned[i] != 'x'){	//is lowercase
							toprintAligned[i++] -= 32;							//convert to uppercase
							if (i == toprintAligned.length()) toprintPost = insertions.front() + toprintPost;
							else toprintAligned.insert(i,insertions.front());
							insertions.erase(insertions.begin());
						}
						else ++i;
					}
					//put insertions back in Post-sequence (as lower case) here
					for (int i = 0; i < toprintPost.length();){
						if (toprintPost[i] > 96 && toprintPost[i] != 'x'){	//is lowercase
							toprintPost[i++] -= 32;							//convert to uppercase
							if (i == toprintPost.length()) toprintPost += insertions.front();
							else toprintPost.insert(i,insertions.front());
							insertions.erase(insertions.begin());
						}
						else ++i;
					}
				}
				
				ssPrint << " " << (al.Position + 1) << " ";   //start position
				
				//Determine & print read size information:
				int readSize = 0;
				for (vector<BamTools::CigarOp>::const_iterator it=al.CigarData.begin(); it < al.CigarData.end(); it++){
					if (it->Type == 'M' || it->Type == 'I' || it->Type == 'S' || it->Type == '=' || it->Type == 'X'){
						readSize += it->Length;         //increment readsize by the length
					}
				}
				ssPrint << readSize << " ";      //read size
				
				//FILTER based on min/max read length restrictions:
				if (settings.readLengthMin && readSize < settings.readLengthMin){ continue; }
				if (settings.readLengthMax && readSize > settings.readLengthMax){ continue; }
			
				//Determine consecutive matching flanking bases (LEFT):
				string::iterator i = PreSeq.end()-1;
				string::iterator i2 = leftReference.end()-1;
				bool consStreak = 1;
				numMatchesL = 0;
				for (int ctr = 0; ctr < PreSeq.length(); ++ctr ) {      //-1 compensates for matching null character @ end of all strings
					if ((*i != *i2) && (*i != *i2 + 32)) {
						consStreak = 0;
						if (ctr < 3){
							if (*i == 'x' || *i == 'S' || (*i2 != '-' && *i == '-') || (*i2 == '-' && *i != '-' )){ 
								continue; //fail the read
							}
						}
					}
					else if (consStreak){ ++numMatchesL;}
					--i; --i2;
				}
				
				//Determine consecutive matching flanking bases (RIGHT):
				i = PostSeq.begin();
				i2 = rightReference.begin();
				consStreak = 1; 
				numMatchesR = 0;
				for (int ctr = 0; ctr < PostSeq.length(); ctr++) { 
					if ((*i != *i2) && (*i != *i2 + 32)){
						consStreak = 0;
						if (ctr < 3){ 
							if (*i == 'x' || *i == 'S' || (*i2 != '-' && *i == '-') || (*i2 == '-' && *i != '-' )){
								continue; //fail the read
							}
						}
					}
					else{
						if (consStreak) ++numMatchesR;
					}
					++i; ++i2;
				}
				
				// Set minflank & print matching # of consecutive bases to the left/right of repeat
				if (numMatchesR < minflank) minflank = numMatchesR;
				else { minflank = numMatchesL; }
				ssPrint << numMatchesL << " " << numMatchesR << " ";  
				
				//FILTER based on consecutive flank bases
				if (numMatchesL < settings.consLeftFlank) continue;
				if (numMatchesR < settings.consRightFlank) continue;
				
				//Print avgBQ:
				ssPrint << "B:" << float(int(10000*avgBQ))/10000 << " ";

				//FILTER based on MapQ, then print MapQ
				if (al.MapQuality < settings.MapQuality) continue;  //MapQuality Filter
				ssPrint << "M:" << al.MapQuality << " ";
				
				//PRINT FLAG STRING:
				ssPrint << "F:";
				if (al.IsPaired()) ssPrint << 'p';
				if (al.IsProperPair()) ssPrint << 'P';
				if (!al.IsMapped()) ssPrint << 'u';
				if (!al.IsMateMapped()) ssPrint << 'U';
				if (al.IsReverseStrand()) ssPrint << 'r';
				if (al.IsMateReverseStrand()) ssPrint << 'R';
				if (al.IsFirstMate()) ssPrint << '1';
				if (al.IsSecondMate()) ssPrint << '2';
				if (!al.IsPrimaryAlignment()) ssPrint << 's';
				if (al.IsFailedQC()) ssPrint << 'f';
				if (al.IsDuplicate()) ssPrint << 'd';
				
				//print CIGAR string:
				ssPrint << " C:";
				for (vector<BamTools::CigarOp>::const_iterator it=al.CigarData.begin(); it < al.CigarData.end(); it++) {
					ssPrint << it->Length;
					ssPrint << it->Type;
				}
				
				//-MULTI filter (check for XT:A:R tag):
				string stringXT;
				al.GetTag("XT",stringXT);
				if (settings.multi && stringXT.find('R',0) != -1) continue;  //if stringXT contains R, ignore read
				
				//-PP filter (check if read is properly paired):
				if (settings.properlyPaired && !al.IsProperPair()){ continue; }
				
				ssPrint << " ID:" << al.Name << endl;
				
				toPrint.push_back( STRING_GT(ssPrint.str(), Sequences(toprintPre, toprintAligned, toprintPost, hasinsertions), AlignedSeq.length() + gtBonus, al.IsProperPair(), al.MapQuality, minflank, al.IsReverseStrand(), avgBQ) );
			}
		}        //end if statements
		
	} //end while loop
	
	numReads = toPrint.size();
	
	//push reference sequences into vectors for expansion & printing:
	toPrint.insert( toPrint.begin(), STRING_GT("\n", Sequences(leftReference, centerReference, rightReference, 0), 0, 0, 0, 0, 0, 0.0) );
	
	// If any of the reads have insertions, expand the reads without inserted bases so all reads are fully printed:
	bool skip = 1;
	for (vector<STRING_GT>::iterator it=toPrint.begin(); it < toPrint.end(); it++){
		if (it->reads.insertions) skip = 0;
	}
	if (true){
		//Handle PRE-SEQ:
		for (int index = 0, limit = settings.LR_CHARS_TO_PRINT + 1; index < limit; ++index){
			for (vector<STRING_GT>::iterator jt=toPrint.begin(); jt < toPrint.end(); jt++){
				if (index >= jt->reads.preSeq.length()) continue;
				if (jt->reads.preSeq[index] == 'B' || jt->reads.preSeq[index] == 'U' || jt->reads.preSeq[index] == 'D' || jt->reads.preSeq[index] == 'H' || jt->reads.preSeq[index] == 'O'){
					limit++;
					
					vector<STRING_GT>::iterator pt;
					if (jt + 1 == toPrint.end()) pt = toPrint.begin();
					else pt = jt+1;
					
					while(pt != jt){
						if (pt->reads.preSeq[index] == 'B' || pt->reads.preSeq[index] == 'U' || pt->reads.preSeq[index] == 'D' || pt->reads.preSeq[index] == 'H' || pt->reads.preSeq[index] == 'O') pt->reads.preSeq[index] -=1;
						else pt->reads.preSeq.insert(index,"-");
						
						if (pt + 1 == toPrint.end()) pt = toPrint.begin();
						else pt++;
					}
					
					jt->reads.preSeq[index] -= 1;
					
				}
			}
		}
		//Handle ALIGNED-SEQ:
		for (int index = 0, limit = target.length() + 1; index < limit; ++index){
			for (vector<STRING_GT>::iterator jt=toPrint.begin(); jt < toPrint.end(); jt++){
				if (index >= jt->reads.alignedSeq.length()) continue;
				if (jt->reads.alignedSeq[index] == 'B' || jt->reads.alignedSeq[index] == 'U' || jt->reads.alignedSeq[index] == 'D' || jt->reads.alignedSeq[index] == 'H' || jt->reads.alignedSeq[index] == 'O'){
					limit++;
					
					vector<STRING_GT>::iterator pt;
					if (jt + 1 == toPrint.end()) pt = toPrint.begin();
					else pt = jt+1;
					
					while(pt != jt){
						if (pt->reads.alignedSeq[index] == 'B' || pt->reads.alignedSeq[index] == 'U' || pt->reads.alignedSeq[index] == 'D' || pt->reads.alignedSeq[index] == 'H' || pt->reads.alignedSeq[index] == 'O') pt->reads.alignedSeq[index] -=1;
						else pt->reads.alignedSeq.insert(index,"-");
						
						if (pt + 1 == toPrint.end()) pt = toPrint.begin();
						else pt++;
					}
					
					jt->reads.alignedSeq[index] -= 1;
					
				}
			}
		}
		//Handle POST-SEQ:
		for (int index = 0, limit = settings.LR_CHARS_TO_PRINT + 1; index < limit; ++index){
			for (vector<STRING_GT>::iterator jt=toPrint.begin(); jt < toPrint.end(); jt++){
				if (index >= jt->reads.postSeq.length()) continue;
				if (jt->reads.postSeq[index] == 'B' || jt->reads.postSeq[index] == 'U' || jt->reads.postSeq[index] == 'D' || jt->reads.postSeq[index] == 'H' || jt->reads.postSeq[index] == 'O'){
					limit++;
					
					vector<STRING_GT>::iterator pt;
					if (jt + 1 == toPrint.end()) pt = toPrint.begin();
					else pt = jt+1;
					
					while(pt != jt){
						if (pt->reads.postSeq[index] == 'B' || pt->reads.postSeq[index] == 'U' || pt->reads.postSeq[index] == 'D' || pt->reads.postSeq[index] == 'H' || pt->reads.postSeq[index] == 'O') pt->reads.postSeq[index] -=1;
						else pt->reads.postSeq.insert(index,"-");
						
						if (pt + 1 == toPrint.end()) pt = toPrint.begin();
						else pt++;
					}
					
					jt->reads.postSeq[index] -= 1;
					
				}
			}
		}
		
		// fix for insertions/deletions immediately following repeat:
		int index = 0;
		while(toPrint.begin()->reads.postSeq[index] == '-'){ ++index; }
		for (vector<STRING_GT>::iterator jt=toPrint.begin(); jt < toPrint.end(); jt++){
			jt->reads.alignedSeq += jt->reads.postSeq.substr(0, index);
			jt->reads.postSeq.erase(0,index);
			
			if (jt->GT){ //if it's not the reference..
				string repeat = jt->reads.alignedSeq;
				repeat.erase(std::remove (repeat.begin(), repeat.end(), '-'), repeat.end());
				
				jt->GT = repeat.length();
			}
		}
	}

	
	// Build VectorGT from toPrint:
	for (vector<STRING_GT>::iterator tP=toPrint.begin(); tP < toPrint.end(); ++tP) {
		if (tP->GT == 0) continue; //ignore reference
		if (!vectorGT.empty()) {
			for (vector<GT>::iterator it = vectorGT.begin(); 1; ) {
				if (it->readlength == tP->GT) {
					it->occurrences += 1;
					it->avgBQ += tP->avgBQ;
					it->avgMinFlank += tP->minFlank;
					if (tP->reverse) it->reverse += 1;
					break;
				}
				else {
					++it;
					if (it == vectorGT.end()) {
						GT a(tP->GT, 1, tP->reverse, tP->minFlank, tP->avgBQ);
						vectorGT.insert(vectorGT.end(), a);
						break;
					}
				}
			}
		}
		else {
			GT a(tP->GT, 1, tP->reverse, tP->minFlank, tP->avgBQ);
			vectorGT.insert(vectorGT.end(), a);
		}
	}
	//average out BQs & flanks
	for (vector<GT>::iterator it = vectorGT.begin(); it < vectorGT.end(); ++it) {
		it->avgBQ /= it->occurrences;
		it->avgMinFlank /= it->occurrences;
	}

	int tallyMapQ = 0;
	int occ = 0;
	double avgMapQ;
	for (vector<STRING_GT>::iterator it=toPrint.begin(); it < toPrint.end(); ++it) {
		tallyMapQ = tallyMapQ + it->MapQ;
		++occ;
	}
	if (occ) avgMapQ = double(tallyMapQ)/double(occ);
	else avgMapQ = -1;
	
	//sort vectorGT by occurrences..
	sort(vectorGT.begin(), vectorGT.end(), vectorGTsort);
	
	//output header line
	oFile << "~" << region << " ";
	oFile << secondColumn;
	oFile << " REF:" << target.length();
	oFile << " A:";
	if (!vectorGT.size()) {
		oFile << "NA ";
		concordance = -1;
		majGT = -1;
	}
	else {
		if (vectorGT.size() == 1) {
			if (numReads == 1) {
				concordance = -1;
				majGT = vectorGT.begin()->readlength;
				oFile << "NA ";
			}
			else {
				concordance = 1;
				majGT = vectorGT.begin()->readlength;
				oFile << vectorGT.begin()->readlength << " ";
				//oFile << "(" << vectorGT.begin()->avgBQ << ")"; //temp
			}
		}
		else {
			for (vector<GT>::iterator it=vectorGT.begin(); it < vectorGT.end(); it++) {
				oFile << it->readlength << "[" << it->occurrences << "] " ;
				//oFile << "(" << it->avgBQ << ") ";
				if (it->occurrences >= occurMajGT) {
					occurMajGT = it->occurrences;
					if (it->readlength > majGT) majGT = it->readlength;
				}
				totalOccurrences += it->occurrences;
			}
			concordance = double(double(occurMajGT)-1.00) / double(double(totalOccurrences)-1.00);
		}
	}
	
	//concordance = # of reads that support the majority GT / total number of reads
	if (concordance < 0) oFile << "C:NA";
	else oFile << " C:" << concordance;
	
	oFile << " D:" << depth << " R:" << numReads << " S:" << numStars;
	if (avgMapQ >= 0) oFile << " M:" << float(int(100*avgMapQ))/100;
	else oFile << " M:NA";
	
	oFile << " GT:";
	callsFile << region << "\t" << secondColumn << "\t";
	vector<int> vGT;
	double conf = 0;
	if (vectorGT.size()  == 0) {
		oFile << "NA L:NA" << endl;
		callsFile << "NA\tNA\n";
	}
	else { 
		vGT = printGenoPerc(vectorGT, target.length(), unitLength, conf, settings.mode); 
		if (numReads <= 1){ conf = 0; }
		//write genotypes to calls & repeats file
		if (vGT.size() == 0) { throw "vGT.size() == 0.. ERROR!\n"; }
		else if (vGT.size() == 1 && conf > 3.02) { oFile << vGT[0] << " L:" << conf << "\n"; callsFile << vGT[0] << '\t' << conf << '\n'; }
		else if (vGT.size() == 2 && conf > 3.02) { oFile << vGT[0] << "h" << vGT[1] << " L:" << conf << "\n"; callsFile << vGT[0] << "h" << vGT[1] << '\t' << conf << '\n'; }
		else{ oFile << "NA L:" << conf << endl; callsFile << "NA\tNA\n"; }
	}
	
	// Set info for printing VCF file
	VCF_INFO INFO;
	INFO.chr = target.startSeq;
	INFO.start = target.startPos + 1;
	INFO.unit = UnitSeq; 
	INFO.length = target.length();
	INFO.purity = purity;
	INFO.depth = numReads;
	INFO.confidence = conf;
	
	// GO THROUGH VECTOR AND PRINT ALL REMAINING
	if (toPrint.size()>1){ //if there are reads present..
		string REF = toPrint[0].reads.alignedSeq;
		bool homo = false;
		if (vGT.size() == 1) homo = true;
		
		for (vector<STRING_GT>::iterator it=toPrint.begin(); it < toPrint.end(); it++) {
			// print .repeats file:
			oFile << it->reads.preSeq << " " << it->reads.alignedSeq << " " << it->reads.postSeq << it->print;
			
			// finished printing to .repeats file.
			if (vGT.size() != 0 && conf > 3.02){
				if (vGT.size() > 1 || vGT[0] != target.length() /*there's been a mutation*/){
					// print .vcf file:
					vector<int>::iterator tempgt = std::find(vGT.begin(), vGT.end(), it->GT);
					if (tempgt != vGT.end() && it->GT != target.length()){
						//debug vcf << "VCF record for " << REF << " --> " << it->reads.alignedSeq << "..\n";
						
						// the read represents one of our genotypes..
						string vcfRecord = getVCF(it->reads.alignedSeq, REF, target.startSeq, target.startPos, *(leftReference.end()-1), homo, INFO);
						vcf << vcfRecord;
						
						//remove the genotype from the genotype list..
						vGT.erase( tempgt ); 
					}	
					// finished printing to .vcf file. 
				}
			}
			
			// continue iterating through each read..
		}
	}
	
	return;
}

inline int nCr (int n, int r){
    return fact(n)/fact(r)/fact(n-r);
}

class tagAndRead{
public:
    string m_name;
    float m_pX;
    tagAndRead(string a, float b){
        m_name = a;
        m_pX = b;
    }
};

inline bool compareTAR(tagAndRead a, tagAndRead b){
    return (a.m_pX > b.m_pX);
}

inline double retBetaMult(int* vector){
	double value = 1, sum = 0;
        for (int i = 0; i < 3; ++i) {
		value += log_factorial[vector[i]-1];
		sum += vector[i];
	}
        value -= log_factorial[int(sum) - 1];
	return value;
	
}

inline vector<int> printGenoPerc(vector<GT> vectorGT, int ref_length, int unit_size, double &confidence, int mode){
	if (ref_length > 70) ref_length = 70;
	if (unit_size > 5) unit_size = 5;
	else if (unit_size < 1) unit_size = 1;
	for (vector<GT>::iterator it = vectorGT.begin(); it < vectorGT.end(); ++it){
		it->avgBQ = -30*log10(it->avgBQ);
		if (it->avgBQ < 0){ it->avgBQ = 0; }
		else if (it->avgBQ > 4){ it->avgBQ = 4; }
	}

	vector<tagAndRead> pXarray;
	vector<int> gts;
	stringstream toReturn;
	extern int PHI_TABLE[5][5][5][2]; 
	
	vectorGT.push_back(GT(0,0,0,0,0.0)); //allows locus to be considered homozygous
    	double pXtotal = 0;
    	string name;
    	
	// Calculate LOCAL_PHI 
    	int mostCommon = 0, secondCommon = 0; double totalSum = 0;
	for (vector<GT>::iterator it = vectorGT.begin(); it < vectorGT.end(); ++it){
	    	int tempOccur = it->occurrences;
		if (tempOccur > mostCommon) {secondCommon = mostCommon; mostCommon = tempOccur;}
        	else if (tempOccur > secondCommon) {secondCommon = tempOccur;}
        	totalSum += tempOccur;
    	}
	double LOCAL_PHI, ERROR;
	if (mode == 1){ LOCAL_PHI = float(totalSum-mostCommon)/totalSum; }
    	else if (mode == 2){ LOCAL_PHI = float(totalSum-(mostCommon+secondCommon))/totalSum; }
	else{ LOCAL_PHI = 0; }
	
	if (vectorGT.size() == 2){ // reference has already been added so 2
		confidence = 100;
		gts.push_back(vectorGT.begin()->readlength);
		return gts;
	}
	
    for (vector<GT>::iterator it = vectorGT.begin(); it < vectorGT.end(); ++it){
        for (vector<GT>::iterator jt = it+1; jt < vectorGT.end(); ++jt){
	    int alleles = 1, errorOccurrences = 0;
	    double CHANCE_error = 1;
            for (vector<GT>::iterator errt = vectorGT.begin(); errt < vectorGT.end(); ++errt){
	    	if (errt != jt && errt != it) { errorOccurrences += errt->occurrences; }
	    }
	    
	    stringstream tempss;
            if (jt->occurrences != 0) {
                tempss << it->readlength <<  "h"  << jt->readlength;
                alleles = 2;
                name = tempss.str();
            }
            else {
                tempss << it->readlength;
                name = tempss.str();
	    }

		//if haploid mode is enabled, ensure only one allele is present:
		if (mode == 1 && alleles == 2) continue;	
			
		//determine likelihood:
		int* ERROR_TABLE_1 = PHI_TABLE[unit_size-1][ref_length/15][int(it->avgBQ)];
		int* ERROR_TABLE_2 = PHI_TABLE[unit_size-1][ref_length/15][int(jt->avgBQ)];
		
		ERROR_TABLE_1[0] = 1;
		ERROR_TABLE_1[1] = 15;
		ERROR_TABLE_2[0] = 1;
		ERROR_TABLE_2[1] = 15;
		int ERROR_1[2]; 
		int ERROR_2[2];

		if (it->occurrences == 0){ ERROR_1[0] = 0; ERROR_1[1] = 0;}
		else { ERROR_1[0] = ERROR_TABLE_1[0]; ERROR_1[1] = ERROR_TABLE_1[1];}
		
		if (jt->occurrences == 0){ ERROR_2[0] = 0; ERROR_2[1] = 0;}
		else { ERROR_2[0] = ERROR_TABLE_2[0]; ERROR_2[1] = ERROR_TABLE_2[1];}


		int v_numerator[3];
		v_numerator[0] = 1 + ERROR_1[1] - ERROR_1[0] + it->occurrences;
		v_numerator[1] = 1 + ERROR_2[1] - ERROR_2[0] + jt->occurrences;
		v_numerator[2] = 1 + ERROR_1[0] + ERROR_2[0] + errorOccurrences; 
		//cout << v_numerator[0] << " ";
		//cout << v_numerator[1] << " ";
		//cout << v_numerator[2] << "\n";

		int v_denom[3];
		v_denom[0] = 1 + ERROR_1[1] - ERROR_1[0];
		v_denom[1] = 1 + ERROR_2[1] - ERROR_2[0];
		v_denom[2] = 1 + ERROR_1[0] + ERROR_2[0]; 
		//cout << v_denom[0] << " ";
		//cout << v_denom[1] << " ";
		//cout << v_denom[2] << "\n";
		
	    	double NUMERATOR = retBetaMult(v_numerator);
		double DENOM = retBetaMult(v_denom);
		//cout << log(retSumFactOverIndFact(it->occurrences,jt->occurrences,errorOccurrences)) << " + " << NUMERATOR << " - " << DENOM << endl;
		//cout << exp(log(retSumFactOverIndFact(it->occurrences,jt->occurrences,errorOccurrences))+NUMERATOR-DENOM) << endl;
		//cout << name << " --> " << "retSum(" << it->occurrences << ", " << jt->occurrences << ", " << errorOccurrences << ")" << endl << endl; // --> " << retSumFactOverIndFact(it->occurrences,jt->occurrences,errorOccurrences) << " * " << CHANCE_allele1 << " * " << CHANCE_allele2 << " * " << CHANCE_error << " == " << retSumFactOverIndFact(it->occurrences,jt->occurrences,errorOccurrences)*CHANCE_allele1*CHANCE_allele2*CHANCE_error << endl;
		
		//add genotype & likelihood to pXarray:
		tagAndRead temp = tagAndRead(name, exp(log(retSumFactOverIndFact(it->occurrences,jt->occurrences,errorOccurrences))+NUMERATOR-DENOM));
		pXarray.push_back(temp);
            	pXtotal += temp.m_pX;
        }
    }
	
    for (vector<tagAndRead>::iterator it = pXarray.begin(); it < pXarray.end(); ++it){
		it->m_pX /= pXtotal;
    		//cout << it->m_name << " --> " << it->m_pX << endl;
	}
    
	// sort, based on likelihood
	sort(pXarray.begin(), pXarray.end(), compareTAR);
	
	// set gts, based on sorted pXarray
	int hpos = pXarray.begin()->m_name.find('h');
	if (hpos == -1){
		//homozygous..
		gts.push_back( atoi(pXarray.begin()->m_name.c_str()) );
	}
	else {
		//heterozygous...
		gts.push_back( atoi(pXarray.begin()->m_name.substr(0, hpos).c_str()) );
		gts.push_back( atoi(pXarray.begin()->m_name.substr(hpos+1, -1).c_str()) );
	}
	
	// set confidence value
	confidence = -10*log10(1-pXarray.begin()->m_pX);
	//cout << confidence << endl << endl;
	if (confidence > 100) confidence = 100; //impose our upper bound to confidence..
	
	//check for NaN --> set to 0
	if (confidence != confidence) {	confidence = 0;	}

	return gts;
}



float fact ( int n ){
    float fact = 1;
    while (n > 1) fact *= n--;
    return fact;
}

// retSumFactOverIndFact(a,b,c) returns fact(a+b+c))/(fact(a)*fact(b)*fact(c))
// (avoids overflow better by avoiding directly computing fact(a+b+c)..)
double retSumFactOverIndFact(int a, int b, int c){
	double val = 1;
	
	// set max/min values
	int max = a, min1 = b, min2 = c;
	if (b > max && b > c){
		max = b;
		min1 = a;
		min2 = c;
	}
	else if (c > max) {
		max = c;
		min1 = a;
		min2 = b;
	}
	
	// determine value to return..
	for (int i = 1; i <= min1; ++i){
		++max;
		val *= max;
		val /= i;
	}
	for (int i = 1; i <= min2; ++i){
		++max;
		val *= max;
		val /= i;
	}
	
	return val;
}


string getVCF(string alignment, string reference, string chr, int start, char precBase, bool homozygous, VCF_INFO info){
	stringstream vcf;
	int begin, end = -1, bothInsOffset = 0;
	
	//assumes alignment & reference are the same length..
	for(int index = alignment.length()-1; index >= 0; --index){
		if ((alignment[index] != reference[index])){
			if (alignment[index] == '-' || reference[index] == '-'){
				end = index;
				while (((alignment[index] != reference[index]) || alignment[index] == '-') && index >= 0){ index -= 1; }
				begin = index;
				while (index >= 0){
					if (reference[index] == '-' && alignment[index] != '-'){ begin -= 1; }
					else if (reference[index] == '-' && alignment[index] == '-'){ bothInsOffset += 1; }
					index -= 1;
				}
				break;
			}
		}
	}
	
	if (end == -1) { return ""; } // no difference was found..
	
	start += begin;
	if (begin == -1) {
		reference = precBase + reference.substr(0, end+1);
		alignment = precBase + alignment.substr(0, end+1);
	}
	else{
		//call getVCF to recursively get VCF for any indels earlier in the sequences...
		vcf << getVCF(alignment.substr(0,begin), reference.substr(0,begin), chr, start-begin, precBase, homozygous, info);
		
		reference = reference.substr(begin, end-begin+1);
		alignment = alignment.substr(begin, end-begin+1);
	}
	
	//remove -'s
	reference.erase( std::remove(reference.begin(), reference.end(), '-'), reference.end() );
	alignment.erase( std::remove(alignment.begin(), alignment.end(), '-'), alignment.end() );
	
	vcf << chr << '\t';
	vcf << start-bothInsOffset << '\t';
	vcf << "." << '\t'; //ID
	vcf << reference << '\t';
	vcf << alignment << '\t';
	vcf << info.confidence << '\t'; //qual -> -10log_10 prob(call in ALT is wrong)
	if (info.confidence > 0.8) vcf << "PASS\t"; //filter
	else vcf << ".\t";
	vcf << "RU=" << info.unit << ";DP=" << info.depth << ";RL=" << info.length << "\t"; //info 
	vcf << "GT\t"; //format
	if (homozygous){ vcf << "1/1\n"; }
	else{ vcf << "1/0\n"; }
	
	return vcf.str();
}

//function to convert phred score to probability score
double PhredToFloat(char chr){
	// p_right-base = 1 - 10^(-Q/10)
	double temp = chr - 33; 
	return (1 - pow(10,temp/-10));
}

//function to ensure filepath is in the current directory
string setToCD (string filepath){
	if (filepath.rfind('/') != -1){ filepath = filepath.substr( filepath.rfind('/') + 1, -1); }
	return filepath;
}

//function to check if a file exists:
bool fileCheck(string strFilename) {
	//attempt to open the file with filename strFilename:
	ifstream ifile;
	ifile.open(strFilename.data(),ios::in);
	//return statement will cast the object to TRUE if the file exists, otherwise false:
	return ifile;
}

//build the fasta index file
void buildFastaIndex(string fastaFileName){
	FastaIndex* fai = new FastaIndex();
	fai->indexReference(fastaFileName.c_str());
	fai->writeIndexFile(fastaFileName + fai->indexFileExtension());
}	

void printHeader(ofstream &vcf){
	vcf << "##fileformat=VCFv4.0" << endl;
	vcf << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
	vcf << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" << endl;
	vcf << "##INFO=<ID=RU,Number=1,Type=String,Description=\"Repeated Unit\">" << endl;
	vcf << "##INFO=<ID=RL,Number=1,Type=Integer,Description=\"Reference Length of Microsatellite\">" << endl;
	vcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" << endl;
}



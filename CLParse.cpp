// Command-line parsing & error reporting module
// 

#include "repeatseq.h"

int PHI_TABLE[5][5][5][2] = { { { { 177 , 7 } , { 234 , 15 } , { 991 , 76 } , { 425 , 39 } , { 234 , 28 } } ,
{ { 324 , 23 } , { 587 , 59 } , { 250 , 28 } , { 123 , 17 } , { 764 , 110 } } ,
{ { 268 , 21 } , { 375 , 40 } , { 139 , 20 } , { 66 , 6 } , { 21 , 2 } } ,
{ { 142 , 7 } , { 25 , 1 } , { 6 , 0 } , { 7 , 0 } , { 2 , 0 } } ,
{ { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } }
} ,
{ { { 948 , 35 } , { 120 , 5 } , { 458 , 29 } , { 174 , 11 } , { 75 , 7 } } ,
{ { 226 , 10 } , { 405 , 22 } , { 159 , 9 } , { 622 , 49 } , { 296 , 18 } } ,
{ { 512 , 33 } , { 111 , 8 } , { 410 , 32 } , { 170 , 17 } , { 61 , 6 } } ,
{ { 636 , 37 } , { 141 , 7 } , { 48 , 3 } , { 12 , 1 } , { 1 , 0 } } ,
{ { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } }
} ,
{ { { 364 , 11 } , { 471 , 40 } , { 176 , 28 } , { 73 , 7 } , { 33 , 8 } } ,
{ { 351 , 12 } , { 637 , 28 } , { 295 , 14 } , { 116 , 16 } , { 62 , 3 } } ,
{ { 460 , 18 } , { 80 , 0 } , { 39 , 5 } , { 13 , 0 } , { 4 , 0 } } ,
{ { 70 , 1 } , { 5 , 1 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } } ,
{ { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } }
} ,
{ { { 337 , 13 } , { 287 , 34 } , { 102 , 18 } , { 21 , 9 } , { 13 , 4 } } ,
{ { 266 , 14 } , { 277 , 37 } , { 90 , 11 } , { 28 , 5 } , { 17 , 3 } } ,
{ { 886 , 68 } , { 106 , 12 } , { 40 , 1 } , { 17 , 1 } , { 7 , 1 } } ,
{ { 170 , 7 } , { 26 , 1 } , { 13 , 0 } , { 2 , 0 } , { 0 , 0 } } ,
{ { 52 , 2 } , { 8 , 0 } , { 1 , 1 } , { 1 , 0 } , { 0 , 0 } }
} ,
{ { { 304 , 6 } , { 302 , 24 } , { 105 , 19 } , { 38 , 7 } , { 7 , 7 } } ,
{ { 192 , 10 } , { 172 , 23 } , { 50 , 11 } , { 23 , 3 } , { 9 , 1 } } ,
{ { 300 , 11 } , { 54 , 13 } , { 15 , 4 } , { 1 , 5 } , { 1 , 0 } } ,
{ { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } } ,
{ { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } , { 0 , 0 } }
}
};


void parseSettings(char *argv[], int argc, SETTINGS_FILTERS &settings, string &bam_file, string &fasta_file, string &position_file){
	// repeatseq [options] <in.bam> <in.fasta> <in.regions>
	if (argc < 4) { throw "Not enough arguments given. Exiting.."; }
	bam_file = argv[argc - 3];
	fasta_file = argv[argc - 2];
	position_file = argv[argc - 1];
	
	for (int i = 1; i < argc-3; ++i) {
		string sw = argv[i];
		
		//SETTINGS:
		if (sw == "-r") {
			//read length select setting
			++i;
			string range = argv[i];
			size_t firstColon = range.find(":");
			if (firstColon == string::npos){
				//no colon was found:
				settings.readLengthMin = atoi(argv[i]);
				settings.readLengthMax = atoi(argv[i]);
				settings.paramString += ".readlength";
				settings.paramString += argv[i];
			}
			else { 
				// a colon was found:
				string first = range.substr(0, firstColon);
				string second = range.substr(firstColon+1,string::npos);
				settings.readLengthMin = atoi(first.c_str());
				settings.readLengthMax = atoi(second.c_str());
				settings.paramString += ".readlength";
				settings.paramString += first;
				settings.paramString += "-";
				settings.paramString += second;
			}
		}
		else if (sw == "-t") {
			//append tag to paramString
			++i;
			settings.paramString += ".";
			settings.paramString += argv[i];
		}
		else if (sw == "-o") {
			//set number of L/R chars to print
			++i;
			settings.LR_CHARS_TO_PRINT = atoi(argv[i]);
		}
		else if (sw == "-error") {
			//set PHI_TABLE
			float constError = 0.05;
			if (string(argv[i+1]).find('-') == -1){
				++i;
				constError = atof(argv[i]);
			}
			for (int w=0;w<5;++w){
				for (int x=0; x<5; ++x){
					for (int y=0; y<5; ++y){
						PHI_TABLE[w][x][y][0] = int(100*constError);
						PHI_TABLE[w][x][y][1] = int(100-100*constError);
					}
				}
			}
		}
		else if (sw == "-haploid") {
			//set haploid/diploid mode
			++i;
			settings.mode = 1;
		}			

		//FILTERS:
		else if (sw == "-pp") {
			//Properly Paired Filter
			settings.properlyPaired = true;
			settings.paramString += ".ProperlyPaired";
		}
		else if (sw == "-L") {
			//Characters that much consecutively match to the LEFT for a read to be used
			++i;
			settings.consLeftFlank = atoi(argv[i]);
			settings.paramString += ".L";
			settings.paramString += argv[i];
		}
		else if (sw == "-R") {
			//Characters that much consecutively match to the RIGHT for a read to be used
			++i;
			settings.consRightFlank = atoi(argv[i]);
			settings.paramString += ".R";
			settings.paramString += argv[i];
		}
		else if (sw == "-M") {
			//MINIMUM MapQuality Score
			++i;
			settings.MapQuality = atoi(argv[i]);
			settings.paramString += ".M";
			settings.paramString += argv[i];
		}
		else if (sw == "-multi") {
			//MULTI Filter (exclude read if XT:A:R tag is present)
			settings.paramString += ".multi";
			settings.multi = true;
		}
		else if (sw == "-repeatseq") { 
			settings.makeRepeatseqFile = true;
		}
		else if (sw == "-calls") {
			settings.makeCallsFile = true;
		}
		else throw "IMPROPER COMMAND LINE ARGUMENT. Exiting..";
	}
}



void printArguments(){
	extern string VERSION;

	cout << endl << "-----------------------------------------------------------\n\n";
	cout << "RepeatSeq v" << VERSION << "\n\n";
	cout << "Usage:\t repeatseq [options] <in.bam> <in.fasta> <in.regions>\n\n";
	cout << "Options:";
	cout << " -r\t\tuse only a specific read length or range of read lengths (e.g. LENGTH or MIN:MAX)";
	cout << "\n\t -L\t\trequired number of reference matching bases BEFORE the repeat [3]";
	cout << "\n\t -R\t\trequired number of reference matching bases AFTER the repeat [3]";
	cout << "\n\t -M\t\tminimum mapping quality for a read to be used for allele determination";
	cout << "\n\t -multi\t\texclude reads flagged with the XT:A:R tag";
	cout << "\n\t -pp\t\texclude reads that are not properly paired (for PE reads only)";
	cout << "\n";
	cout << "\n\t -error\t\tmanually override the RepeatSeq error model and set a constant error rate [0.05]";
	cout << "\n\t -haploid\tassume a haploid rather than diploid genome";
	cout << "\n";	
	cout << "\n\t -repeatseq\twrite .repeatseq file containing additional information about reads";
	cout << "\n\t -calls\t\twrite .calls file";
	cout << "\n\t -t\t\tinclude user-defined tag in the output filename";
	cout << "\n\t -o\t\tnumber of flanking bases to output from each read";
	cout << "\n";
	cout << endl << "-----------------------------------------------------------" << endl;
}







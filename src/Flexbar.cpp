/*==================================================

   Flexbar - flexible barcode and adapter removal

   Version 3.5.0

   BSD 3-Clause License

   uses SeqAn library release 2.4.0
   and TBB library 4.0 or later


             Developer:  Johannes Roehr

   Former contributors:  Matthias Dodt
                         Benjamin Menkuec
                         Sebastian Roskosch


   https://github.com/seqan/flexbar

===================================================*/

#include "Flexbar.h"


int main(int argc, const char* argv[]){

	using namespace std;
	using namespace seqan;

	const string version = "3.5.0";
	const string date    = "May 2019";

	ArgumentParser parser("flexbar");

	defineOptions(parser, version, date);
	parseCmdLine(parser, version, argc, argv);

	Options o;

	initOptions(o, parser);
	loadOptions(o, parser);

        //TODO extract Reads
        //additionally make sure they are unique
        //TODO std::move??
        {
            std::vector<BamAlignmentRecord > recordstable = extractReads(o);
            removeCDNA(o,recordstable);
        }
        std::cout << "TargetName: " << o.targetName << "\n";



        //TODO removeCDNA
        //tbb::concurrent_vector<flexbar::TBar> bars;

        //P1 alignment

        //use adapter para
	startComputation(o);


        splitReads(o);

        //TODO use own parameters //overwrite adapter parameters in options
        //o.out = &o.fstrmOut;
        startComputation(o);

        //o.out = &o.fstrmOut;
        //TODO //overwrite adapter parameters in options
        startComputation(o);

	return 0;
}

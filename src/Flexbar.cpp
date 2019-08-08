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

	startComputation(o);




    //o.readsFile     for extractReads
//     o.readsFile = "";


    getOptionValue(o.adapterFile, parser, "adapters");
    *out << "Adapter file:          " << o.adapterFile << endl;
    o.adapRm = NORMAL;
    o.useAdapterFile = true;


    o.bundleSize = 20;
	o.a_match = o.barcode_match;
    o.a_mismatch = o.barcode_mismatch;
    o.a_gapCost = o.barcode_gapCost;
	o.a_errorRate = o.barcode_errorRate;
    o.a_end = LTAIL;
    rcMode    = RCOFF;

    //comp
    o.a_end = RTAIL;
    o.rcMode = RCONLY;

	return 0;
}

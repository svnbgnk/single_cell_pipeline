#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <iostream>
#include <seqan/align.h>
#include <seqan/bam_io.h>
#include<iostream>
#include<fstream>
#include <seqan/find.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <map>
#include <iostream>
#include <fstream>

using namespace seqan;
using namespace std;



std::vector<BamAlignmentRecord > extractReads(Options &o)
{
//     o.regionsFile and
//     o.readsFile
/*
    ArgumentParser parser("Extract Reads");

    addOption(parser, ArgParseOption("b", "bam", "Path to the bam file", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "bam");

    addOption(parser, ArgParseOption("g", "gtf", "Path to the gtf file", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "gtf");


    addOption(parser, ArgParseOption("su", "suffix", "Suffix concatenated to ouput files (beside .bam)",
                                     ArgParseOption::STRING));

    addOption(parser, ArgParseOption("o", "output", "Path to the output prefix", ArgParseArgument::INPUT_FILE, "IN"));

    addOption(parser, ArgParseOption("ov", "overlap", "Extend each region by that amount", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("s", "step", "Number of exons in single bam", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("t", "threads", "Number of threads", ArgParseArgument::INTEGER, "INT"));




    CharString bamPath, gtfPath, outputPathPrefix = "reads", suffix = "";
    int batchSize1 = 100000;

    int barcode_umi_length = 30;
    int threshold = 5;
    int first = 9999999;
    int threads = 1;
    int step = 1;
    int overlapI = 0;

    getOptionValue(bamPath, parser, "bam");
    getOptionValue(gtfPath, parser, "gtf");
    getOptionValue(overlapI, parser, "overlap");
    getOptionValue(outputPathPrefix, parser, "output");
    getOptionValue(suffix, parser, "suffix");
    getOptionValue(step, parser, "step");
    getOptionValue(threads, parser, "threads");
//     getOptionValue(barcodeLength, parser, "barcodeL");
//     getOptionValue(first, parser, "first");
//     getOptionValue(threshold, parser, "threshold");
//     bool verbose = isSet(parser, "verbose");
  */
    int const overlap = 0;

    //Read gtf file
    ifstream inputFile(o.regionsFile);
    string line;
    vector<std::tuple<CharString, uint32_t, uint32_t> > table;
    if(inputFile.is_open()){
//         bool head = true;
        while (getline(inputFile, line, '\n'))
        {
            if(line.length() < 10)
                continue;
//             if(head){
//                 head = false;
//                 continue;
//             }
            std::tuple<CharString, uint32_t, uint32_t> row;
            std::istringstream sline(line);
            string element;
            std::vector<string> tmprow;
            int k = 0;
            while (getline(sline, element, '\t')){
//                 std::cout << "k: " << k << "\t" << element << "\n";
                if(k == 0 || k == 3 || k == 4)
                    tmprow.push_back(element);
                ++k;
            }
            CharString refid = tmprow[0];
//             std::cout << tmprow[0] << "\t" << tmprow[1] << "\t" << tmprow[2] << "\n";
            int s = std::stoi(tmprow[1]);
            int e = std::stoi(tmprow[2]);
            table.push_back(std::make_tuple(refid, s, e));
        }
//         std::cout << "Finished reading\n";
        /*
        for (unsigned i = 0; i < table.size(); i++){
                std::cout << std::get<0>(table[i]) << "\t" << std::get<1>(table[i]) << "\t" << std::get<2>(table[i]) << "\n";
        }
        std::cout << "\n";*/
        inputFile.close();
    }
    else
    {
        std::cout << "Could not open file: " << o.regionsFile << "\n";
    }

    // Open input file, BamFileIn can read SAM and BAM files.
    BamFileIn bamFileIn;
    if (!open(bamFileIn, o.readsFile))
    {
        std::cerr << "ERROR: Could not open " << o.readsFile << std::endl;
        return 1;
    }

    BamHeader header;
    std::vector<std::vector<BamAlignmentRecord > > recordtable;
    try
    {
        // Copy header.
        readHeader(header, bamFileIn);
        // Copy records.
        BamAlignmentRecord record;
        recordtable.resize(table.size());

        string lastContig = "___";
        int st = 0;
        int end = 0;
        bool stNew = true;
        bool endNew = false;
        bool skip = false;
        uint32_t k = 0;
        while (!atEnd(bamFileIn))
        {
            readRecord(record, bamFileIn);
            //use map to jump to correct chromosom //use start pos and length
            if(hasFlagUnmapped(record))
                continue;
            string recordContig = toCString(getContigName(record, bamFileIn));
            uint32_t recordBegin = record.beginPos;
            uint32_t recordEnd = recordBegin + length(record.seq);
            //determine search range
//             bool same = true;
            if(lastContig.compare(recordContig) != 0){
                std::cout << "Calc new Range for: \n";
                if(!skip)
                {
                    for(int i = st; i < end; ++i)
                        table.erase(table.begin() + st);
                }
//                 same = false;
                skip = false;
                stNew = true;
                string rowContig;
                for(int i = 0; i < table.size(); ++i){
                    rowContig = toCString(std::get<0>(table[i]));
                    if(recordContig.compare(rowContig) == 0 && stNew){
                        stNew = false;
                        st = i;
                    }
                    if(!stNew && recordContig.compare(rowContig) != 0)
                    {
                        end = i;
                        break;
                    }
                }

                // in case last element matches
                if(recordContig.compare(rowContig) == 0)
                    end = table.size();

                lastContig = recordContig;
                std::cout << "Chrom: " << lastContig << "\n";
                std::cout << "Start: " << st << "\tEnd: " << end << "\n";
                std::cout << "Record: " << k << "\n";

                // in case no element matches
                if(stNew){
                    std::cout << "Skip not in gtf file: " << lastContig << "\n";
                    skip = true;
                }

            }
            ++k;

            if(skip){
                st = table.size();
                end = table.size();
                continue;
            }
            //TODO assume sorted
            #pragma omp parallel for num_threads(threads) schedule(static)
            for(int i = st; i < end; ++i){
                //check row
                string rowContig = toCString(std::get<0>(table[i]));
                uint32_t rowBegin = std::get<1>(table[i]);
                uint32_t rowEnd = std::get<2>(table[i]);
                //std::cout << "recordBegin: " << recordBegin << "\t" << recordEnd << "\trow: " << rowBegin << "\t" << rowEnd << "\n";
                if((recordBegin + overlap >= rowBegin && recordBegin <= rowEnd + overlap) || (recordEnd + overlap >= rowBegin && recordEnd <= rowEnd + overlap)
                    || (recordBegin <= rowBegin + overlap && recordEnd + overlap >= rowEnd))
                {
                    //last record did not match sorted coordinates therefore next does not also
//                     if(same && st < i)
//                         st = i;

                    recordtable[i].push_back(record);
                }
            }
        }
/*
        string prefix = toCString(outputPathPrefix);
        #pragma omp parallel for schedule(dynamic) num_threads(threads)
        for(int b = 0; b < recordtable.size(); b += step)
        {
            int cumLength = 0;
            for(int j = 0; j < step && (b + j) < recordtable.size(); ++j){
                cumLength += recordtable[b + j].size();
            }
            if(cumLength == 0)
                continue;

            ofstream mybamstream;
            string bamName = prefix + to_string(b) + toCString(suffix) + ".bam";
            mybamstream.open(bamName);
            // Open output file, BamFileOut accepts also an ostream and a format tag.
            BamFileOut bamFileOut(context(bamFileIn), mybamstream, Bam());
            writeHeader(bamFileOut, header);
            for(int j = 0; j < step && (b + j) < recordtable.size(); ++j){
                for(int i = 0; i < recordtable[b + j].size(); ++i){
                    writeRecord(bamFileOut, recordtable[b + j][i]);
                }
            }
            close(bamFileOut);
        }*/
    }


    catch (Exception const & e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    std::vector<BamAlignmentRecord > uniqueExtractedReads;
    //post processing making reads unique
    map<TString, short> idMap;
    uint32_t dups = 0;
    uint32_t empty_dups = 0;

    //TODO use move int
    for(int i = 0; i < recordtable.size(); ++i){
        for(int j = 0; j < recordtable[i].size(); ++j){
            auto & record = recordtable[i][j];
            if(length(record.seq) == 0){
                ++empty_dups;

            }
            else
            {
                if(!idMap.count(record.qName) == 1){
                    idMap[record.qName] = 1;
                    uniqueExtractedReads.push_back(std::move(recordtable[i][j]));
                }
                else
                {
                    ++dups;
                }
            }
        }
    }

    std::cout "Finished extraction. Removed " << (int)dups << " duplicates." << "Duplicates with no sequence: " << (int)empty_dups << "\n";

    return uniqueExtractedReads;
}
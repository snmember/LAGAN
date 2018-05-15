#include<iostream>
#include<seqan/file.h>
#include<seqan/sequence.h>
#include<seqan/arg_parse.h>

using namespace std;
using namespace seqan;


int main(int argc, char* argv[])
{

    if (argc < 2) {
      cerr << "Bitte zwei Dateien angeben!" << endl;
      return 1;
    }
  
    CharString seqFileName1 = getAbsolutePath(argv[0]);
    CharString seqFileName2 = getAbsolutePath(argv[1]);
    CharString id1;
    Dna5String seq1;
    CharString id2;
    Dna5String seq2;

    SeqFileIn seqFileIn;
    if (!open(seqFileIn1, toCString(seqFileName1)))
    {
        std::cerr << "ERROR: Could not open the first file.\n";
        return 1;
    }

    if (!open(seqFileIn2, toCString(seqFileName2)))
    {
        std::cerr << "ERROR: Could not open the second file.\n";
        return 1;
    }

    try
    {
        readRecord(id, seq, seqFileIn1);
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR in 1: " << e.what() << std::endl;
        return 1;
    }

    try
    {
        readRecord(id, seq, seqFileIn2);
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR in 2: " << e.what() << std::endl;
        return 1;
    }

    std::cout << id << '\t' << seq << '\n';

    return 0;
}

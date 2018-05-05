#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include <algorithm>
#include <string>
using namespace seqan;


typedef std::tuple<unsigned, size_t> QGram; //Qgram Tupel
typedef Iterator<String<QGram>, Rooted>::Type QIter; //Rooted Iterator
typedef Seed<Simple> TSeed; //Seed Klasse
typedef SeedSet<TSeed> TSeedSet; //Seed Set
typedef Iterator<TSeedSet, Rooted>::Type SIter;

String<QGram> qGram1;
String<QGram> qGram2; //Strings, welche die HashWerte des QGram Index beinhalten sollen.

int main(int argc, char const *argv[]) {

   seqan::ArgumentParser parser("SimpleLagan");
   addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "Path1"));
   addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "Path2"));
   addOption(parser, seqan::ArgParseOption("q", "seed", "Set the length of the seed",seqan::ArgParseArgument::INTEGER, "INT"));
   setDefaultValue(parser, "seed", "10");
   seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

   if (res != seqan::ArgumentParser::PARSE_OK)
       return res == seqan::ArgumentParser::PARSE_ERROR;

   unsigned q = 0;
   getOptionValue(q, parser, "seed");
   std::cout << q << std::endl;
   CharString path1;
   CharString path2;

   getArgumentValue(path1, parser, 0);
   getArgumentValue(path2, parser, 1);

   CharString id1;
   CharString id2;
   Dna5String seq1;
   Dna5String seq2;


   SeqFileIn seqFileIn(toCString(path1));
   readRecord(id1, seq1, seqFileIn);

   SeqFileIn seqFileIn2(toCString(path2));
   readRecord(id2, seq2, seqFileIn2);
   //einlesen der Fasta Dateien


   Shape<Dna, UngappedShape<10> > myShape; //erstmals festgelegte Seedlänge, soll noch als optionale Eingabe modifiziert werden.


   hashInit(myShape, begin(seq1));
   for (size_t i = 0; i < length(seq1) - length(myShape) + 1; ++i)
       appendValue(qGram1, std::make_tuple(hashNext(myShape, begin(seq1) + i), i)) ;

   hashInit(myShape, begin(seq2));
   for (size_t i = 0; i < length(seq2) - length(myShape) + 1; ++i)
       appendValue(qGram2, std::make_tuple(hashNext(myShape, begin(seq2) + i), i)) ; //QGram Index und Position in Tupel speichern.

   std::sort(begin(qGram1), end(qGram1), [] (QGram const & elem1, QGram const & elem2) // sortiert QGram Liste anhand des Hashwerts
   {
       return std::get<0>(elem1) < std::get<0>(elem2);
   });

   std::sort(begin(qGram2), end(qGram2), [] (QGram const & elem1, QGram const & elem2) // sortiert QGram Liste anhand des Hashwerts
   {
       return std::get<0>(elem1) < std::get<0>(elem2);
   });


   TSeedSet seedSet;
   QIter it1 = begin(qGram1);
   QIter it2 = begin(qGram2); //setze Iteratoren jeweils an den Anfang meiner Strings

   std::cout << std::get<0>(*it1) << std::endl << std::get<0>(*it2); //Das war ein Test, ob es tatsächlich Hash Werte ausgibt.

   while (!atEnd(it1) || !atEnd(it2)) //solange keiner der Iteratoren hinter dem letzten Element landet soll Loop ausgeführt werden
   {
       if(std::get<0>(*it1) > std::get<0>(*it2)) { goNext(it2); }

       if(std::get<0>(*it1) < std::get<0>(*it2)) { goNext(it1); } //Iteratoren entsprechend ihrem Größenverhältnis bewegt.

       else
       {
           while(std::get<0>(*it1) == std::get<0>(*it2))
           {
               addSeed(seedSet, TSeed(std::get<1>(*it1),std::get<1>(*it2), 9),Single());
               goNext(it2);
           }
           --it2;
           while(std::get<0>(*it1) == std::get<0>(*it2))
           {
               addSeed(seedSet,TSeed(std::get<1>(*it1),std::get<1>(*it2), 9),Single());
               goNext(it1);
           }// den Teil muss ich noch etwas bearbeiten, da mehrere Fälle abgedeckt werden müssen.
           goNext(it2);
       }
   }

/*
   Score<int, Simple> scoringScheme(1, -1, -1);
   for (SIter it3 = begin(seedSet, Standard()); it3 != end(seedSet, Standard()); ++it3)
       extendSeed(*it3, seq1, seq2, EXTEND_BOTH, scoringScheme, 5, GappedXDrop()); Der Teil funktioniert noch nicht
*/




  return 0;
}

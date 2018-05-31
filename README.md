# LAGAN
Software Praktikum: Seqan

LAGAN is an efficient algorithm to generate a gobal sequence alignment by generating seeds and running a banded sequence alignment. LAGAN was published 2013 by Michael Brudo, Chuong B.\,Do, Gregory M.\,Cooper, et al. in 2003.
This LAGAN program was implemented using the SeqAn library in C++.

To use the lagan.cpp you have to follow the SeqAn Userguide on http://seqan.readthedocs.io/en/master/Infrastructure/Use/index.html#infra-use

Then you have to download the lagan.cpp and include it into your CMakeLists.txt

After compiling you can use LAGAN by calling ./lagan -q --insert seed length (optional)-- pathtofile1 pathtofile2

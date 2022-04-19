#include <sdsl/csa_bitcompressed.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/utility.hpp>

#include <sdsl/lcp.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <stdio.h>

// buildsa: Input
// The input consists of one optional parameter, as well as 2 required arguments, given in this order:

// --preftab <k> - if the option --preftab is passed to the buildsa executable (with the parameter k), then a prefix table will be built atop the suffix array, capable of jumping to the suffix array interval corresponding to any prefix of length k.
// reference - the path to a FASTA format file containing the reference of which you will build the suffix array. Note
// output - the program will write a single binary output file to a file with this name, that contains a serialized version of the input string and the suffix array.

// buildsa --preftab <k> <reference> <output>

// buildsa <reference> <output>
using namespace sdsl;
using namespace std;
using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

int main(int argc, char **argv)
{
    int refIdx = 1;
    int outputIdx = 2;
    // cout << argc << endl;
    if (strcmp(argv[1],"--preftab") == 0)
    {
        refIdx += 2;
        outputIdx += 2;
    }

    string ref_fname = argv[refIdx];
    string out_fname = argv[refIdx + 1];
    std::ifstream file(ref_fname);
    bool first = true;
    string ref;
    if (file.is_open())
    {
        std::string line;
        while (std::getline(file, line))
        {
            // using printf() in all tests for consistency
            if (!first)
            {
                ref = ref + line;
            }
            first = false;
        }
        file.close();
    }
    else
    {
        cout << "error" << endl;
    }

    string tname = "tmp";
    std::ofstream tfile(tname);
    tfile << ref.c_str();
    tfile.close();
    // printf("%s", ref.c_str());

    cache_config cc(false); // do not delete temp files after csa construction
    csa_wt<> csa;
    
    auto start = std::chrono::steady_clock::now();
    construct(csa, tname, 1);

    auto end = std::chrono::steady_clock::now();
    cout << "csa const time for " << ref_fname << " " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << endl;
    unordered_map<string, pair<int, int>> preft = unordered_map<string, pair<int, int>>();
    if (strcmp(argv[1],"--preftab") == 0) {
        int k = stoi(argv[2]);
        string prev = ref.substr(csa[1], k);
        cout << "building preftab with k=" <<k<< endl;
        start = std::chrono::steady_clock::now();
        preft[prev].first = 1;
        string curr;
        for (int i = 2; i < csa.size(); i++)
        {
            curr = ref.substr(csa[i], k);
            if (prev.compare(curr) != 0)
            {
                preft[prev].second = i - 1;
                preft[curr] = {i, i};
            }
            prev = curr;
            // cout << i << " " << ref.substr(csa[i]) << endl;
        }
        end = std::chrono::steady_clock::now();
        cout << "preftabe const time for " << ref_fname << " " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << endl;
    }
    std::ofstream ofile(out_fname);
    cereal::BinaryOutputArchive oarchive(ofile); // Create an output archive

    oarchive(preft);
    oarchive(ref);
    csa.serialize(ofile);
}
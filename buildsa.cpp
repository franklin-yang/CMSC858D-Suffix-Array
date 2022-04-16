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
#include "SA.hpp"

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

void build_sa(string ref_fname, string out_fname)
{
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
        cout << 34 << endl;
    }

    string tname = "tmp";
    std::ofstream tfile(tname);
    tfile << ref.c_str();
    tfile.close();
    printf("%s", ref.c_str());
    SA sufa;

    cache_config cc(false); // do not delete temp files after csa construction
    csa_wt<> csa;
    construct(csa, tname, 1);

    cc.delete_files = true; // delete temp files after lcp construction
    lcp_wt<> lcp;
    construct(lcp, tname, 1);
    // cout << ref.length();
    // sufa.ref = ref;
    int k = 2;
    unordered_map<string, pair<int, int>> preft = unordered_map<string, pair<int, int>>();
    string prev = ref.substr(csa[1], k);
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

    // for (auto const &pair : preft)
    // {
    //     std::cout << "{" << pair.first << ": " << pair.second.first << " " << pair.second.second << "}\n";
    // }
    // if (csa.size() < 1000)
    // {
    //     cout << csa << endl;
    //     cout << "-------" << endl;
    //     cout << lcp << endl;
    // }
    // }

    std::ofstream ofile(out_fname);
    cereal::BinaryOutputArchive oarchive(ofile); // Create an output archive

    // oarchive(sufa);
    // sa.serialize(ofile);
    // ofile.close();
    oarchive(preft);
    csa.serialize(ofile);
    lcp.serialize(ofile);
}

int main(int argc, char **argv)
{
    int refIdx = 1;
    int outputIdx = 2;

    if (argc == 4)
    {
        refIdx += 2;
        outputIdx += 2;
    }

    // std::string line;
    // std::ifstream myfile (argv[refIdx]);
    // bool first = true;
    // if (myfile.is_open())
    // {
    //     while ( getline (myfile,line) )
    //     {
    //     // std::cout << line << '\n';
    //     first = false;
    //     if (!first) {

    //     }
    //     }
    //     myfile.close();
    // }
    build_sa(argv[refIdx], argv[refIdx + 1]);

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
        cout << 34 << endl;
    }

    string tname = "tmp";
    std::ofstream tfile(tname);
    tfile << ref.c_str();
    tfile.close();
    printf("%s", ref.c_str());
    SA sufa;

    cache_config cc(false); // do not delete temp files after csa construction
    csa_wt<> csa;
    construct(csa, tname, 1);

    cc.delete_files = true; // delete temp files after lcp construction
    lcp_wt<> lcp;
    construct(lcp, tname, 1);
    // cout << ref.length();
    // sufa.ref = ref;
    int k = 2;
    unordered_map<string, pair<int, int>> preft = unordered_map<string, pair<int, int>>();
    string prev = ref.substr(csa[1], k);
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

    // for (auto const &pair : preft)
    // {
    //     std::cout << "{" << pair.first << ": " << pair.second.first << " " << pair.second.second << "}\n";
    // }
    if (csa.size() < 1000)
    {
        // cout << csa << endl;
        // cout << "-------" << endl;
        // cout << lcp << endl;
    }
    // }

    std::ofstream ofile(out_fname);
    cereal::BinaryOutputArchive oarchive(ofile); // Create an output archive

    oarchive(preft);
    oarchive(ref);
    csa.serialize(ofile);
}
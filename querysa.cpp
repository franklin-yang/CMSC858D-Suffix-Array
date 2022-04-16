#include <sdsl/csa_bitcompressed.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <sdsl/suffix_arrays.hpp>

#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/utility.hpp>
// querysa: Input
// index - the path to the binary file containing your serialized suffix array (as written by buildsa above).
// queries - the path to an input file in FASTA format containing a set of records. You will need to care about both the name and sequence of these fasta records, as you will report the output using the name that appears for a record. Note, query sequences can span more than one line (headers will always be on one line).
// query mode - this argument should be one of two strings; either naive or simpaccel. If the string is naive you should perform your queries using the naive binary search algorithm. If the string is simpaccel you should perform your queries using the “simple accelerant” algorithm we covered in class. Note: If you are reading a serialized input file with no prefix lookup table, then these algorithms are to be run on each query on the full suffix array, and if you are reading in a prefix lookup table as well, then this is the algorithm that should be used on the relevant interval for each query.
// output - the name to use for the resulting output.

using namespace sdsl;
using namespace std;
using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

int lcp(string s1, string s2)
{
    int i = 0;
    for (i = 0; i < min(s1.length(),s2.length()); i++) {
        if (s1.at(i) != s2.at(i)) {
            break;
        }
    }
    return i;
}
int main(int argc, char **argv)
{
    string idx_fname = argv[1];
    std::ifstream ifile(idx_fname);
    cereal::BinaryInputArchive oarchive(ifile); // Create an output archive
    unordered_map<string, pair<int, int>> preft = unordered_map<string, pair<int, int>>();
    string ref;
    oarchive(preft);
    oarchive(ref);
    
    csa_wt<> csa;
    csa.load(ifile);

    // for (int i = 0; i < csa.size(); i++)
    // {
    //     cout << i << " " << ref.substr(csa[i]) << endl;
    // }

    string q_fname = argv[2];
    cout << q_fname << endl;
    std::ifstream file(q_fname);
    
    unordered_map<string,string> queries = unordered_map<string,string>();
    unordered_map<string,pair<int,int>> matches = unordered_map<string,pair<int,int>>();
    string prev_name;
    if (file.is_open())
    {
    // cout << q_fname << endl;
        std::string line;
        while (std::getline(file, line))
        {
            // using printf() in all tests for consistency
            if (line.at(0) == '>' ) {
                // refname
                prev_name = line.substr(1);
                // cout << qname << endl;
                // prev_name = 
                queries[line.substr(1)] = "";
            } else {
                queries[prev_name].append(line);
            }
        }
        file.close();
    }
    else
    {
        cout << 34 << endl;
    }


    for (auto const &pair : queries)
    {
        int least_idx;
        string p = pair.second;
        string name = pair.first;
        matches[name] = {-1,-2};

        int left_idx = 0;
        int right_idx = csa.size();
        int mid = (left_idx + right_idx) / 2;
        string sac;
        bool found = false;
        int cmp;
        while (right_idx >= left_idx) {
            mid = (left_idx + right_idx) / 2;
            sac = ref.substr(csa[mid],p.length());
            cmp = p.compare(sac);
            if (cmp < 0) {
                right_idx = mid - 1;
            } else if (cmp > 0) {
                left_idx = mid + 1;
            } else if (cmp == 0) {
                if (mid == 0 || p.compare(ref.substr(csa[mid - 1],p.length())) > 0) {
                    least_idx = mid;
                    found = true;
                    break;
                }
                right_idx = mid - 1;
            }
            cout << mid << endl;
            // sac = ref.substr(mid);
        }

        int most_idx;
         left_idx = least_idx;
         right_idx = csa.size() - 1;
         mid = (left_idx + right_idx) / 2;
        while (right_idx >= left_idx) {
            mid = (left_idx + right_idx) / 2;
            sac = ref.substr(csa[mid],p.length());
            cmp = p.compare(sac);
            if (cmp < 0) {
                right_idx = mid - 1;
            } else if (cmp > 0) {
                left_idx = mid + 1;
            } else if (cmp == 0) {
                if (mid == csa.size()-1 || p.compare(ref.substr(csa[mid + 1],p.length())) < 0) {
                    most_idx = mid;
                    break;
                }
                left_idx = mid + 1;
            }
        }
        if (found) {
            matches[name] = {least_idx,most_idx};
        }
        std::cout << "{" << pair.first << ": " << pair.second << "}\n" << least_idx << " " << most_idx;
        // return 0;
        // break;
    }

    string out_fname = argv[4];
    cout << out_fname;
    ofstream ofile(out_fname);
    int nm;
    for (auto const &pair : matches)
    {
        ofile << pair.first << "\t";
        nm = pair.second.second - pair.second.first + 1;
        ofile << nm ;
        for (int j = pair.second.first; j <= pair.second.second; j++) {
            ofile << "\t" << csa[j];
        }
        ofile << endl;
        std::cout << "{" << pair.first << ": " << pair.second.first << " " << pair.second.second << "}\n";
    }


    // std::string name = "13 dna:chromosome chromosome:GRCh38:13:1:114364328:1 REF";
}
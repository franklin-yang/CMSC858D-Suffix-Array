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

int lcp(string &s1, int offset, string &s2)
{
    int i = 0;
    for (i = 0; i < min(s1.length()-offset, s2.length()); i++)
    {
        if (s1.at(i+offset) != s2.at(i))
        {
            break;
        }
    }
    return i;
}


int lcp2(string &s1, int offset, string &s2, int offset2)
{
    int i = 0;
    for (i = 0; i < min(s1.length()-offset, s2.length()-offset2); i++)
    {
        if (s1.at(i+offset) != s2.at(i+offset2))
        {
            break;
        }
    }
    return i;
}

int subcmp(string &s1, int offset, string &s2, int offset2)
{
    int i=0;
    for ( i; i < min(s1.length()-offset, s2.length()-offset2); i++)
    {
        if (s1.at(i+offset) != s2.at(i+offset2))
        {
            // return s2.at(i+offset2) -  s1.at(i+offset);
            if ( s2.at(i+offset2) -  s1.at(i+offset) > 0) {
                return i+1;
            } else {
                return -i-1;
            }
        }
    }
    // cout << "i: " << i << "offset" << offset2 << "s2l" << s2.length() << endl;
    if ( (i+offset2) == s2.length()){
        return 0;
    }
   if (s1.length()-offset < s2.length()-offset2)
   {
       return 1;
   } else if (s1.length()-offset > s2.length()-offset2) {
       return -1;
   } else {
       return 0;
   }
   return 0;
}

pair<unordered_map<string, pair<int, int>>, double> naive(unordered_map<string, string> queries, csa_wt<> csa, string ref)
{
    unordered_map<string, pair<int, int>> matches;
    auto start = std::chrono::steady_clock::now();
    for (auto const &pair : queries)
    {
        int least_idx;
        string p = pair.second;
        string name = pair.first;
        matches[name] = {-1, -2};

        int left_idx = 0;
        int right_idx = csa.size();
        int mid = (left_idx + right_idx) / 2;
        string sac;
        bool found = false;
        int cmp;
        while (right_idx >= left_idx)
        {
            mid = (left_idx + right_idx) / 2;
            // sac = ref.substr(csa[mid], p.length());
            // cmp = p.compare(sac);
            cmp = subcmp(ref, csa[mid], p,0);
            if (cmp < 0)
            {
                right_idx = mid - 1;
            }
            else if (cmp > 0)
            {
                left_idx = mid + 1;
            }
            else if (cmp == 0)
            {
                if (mid == 0 || p.compare(ref.substr(csa[mid - 1], p.length())) > 0)
                {
                    least_idx = mid;
                    found = true;
                    break;
                }
                right_idx = mid - 1;
            }
            // cout << mid << endl;
            // sac = ref.substr(mid);
        }

        int most_idx;
        left_idx = least_idx;
        right_idx = csa.size() - 1;
        mid = (left_idx + right_idx) / 2;
        while (right_idx >= left_idx)
        {
            mid = (left_idx + right_idx) / 2;
            // sac = ref.substr(csa[mid], p.length());
            // cmp = p.compare(sac);
            
            cmp = subcmp(ref, csa[mid], p,0);
            if (cmp < 0)
            {
                right_idx = mid - 1;
            }
            else if (cmp > 0)
            {
                left_idx = mid + 1;
            }
            else if (cmp == 0)
            {
                if (mid == csa.size() - 1 || p.compare(ref.substr(csa[mid + 1], p.length())) < 0)
                {
                    most_idx = mid;
                    break;
                }
                left_idx = mid + 1;
            }
        }
        if (found)
        {
            matches[name] = {least_idx, most_idx};
        }
    }
    auto end = std::chrono::steady_clock::now();
    return {matches, std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()};
}

pair<unordered_map<string, pair<int, int>>, double> naivepreft(unordered_map<string, string> queries, csa_wt<> csa, string ref, unordered_map<string, pair<int, int>> preft, int k)
{
    string p;
    string name;
    string pref;
    unordered_map<string, pair<int, int>> matches;
    auto start = std::chrono::steady_clock::now();
    for (auto const &pair : queries)
    {
        int least_idx;
        p = pair.second;
        name = pair.first;
        pref = p.substr(0, k);
        matches[name] = {-1, -2};

        int left_idx = preft[pref].first;
        int right_idx = preft[pref].second;
        int mid = (left_idx + right_idx) / 2;
        string sac;
        bool found = false;
        int cmp;
        while (right_idx >= left_idx)
        {
            mid = (left_idx + right_idx) / 2;
            sac = ref.substr(csa[mid], p.length());
            cmp = p.compare(sac);
            if (cmp < 0)
            {
                right_idx = mid - 1;
            }
            else if (cmp > 0)
            {
                left_idx = mid + 1;
            }
            else if (cmp == 0)
            {
                if (mid == 0 || p.compare(ref.substr(csa[mid - 1], p.length())) > 0)
                {
                    least_idx = mid;
                    found = true;
                    break;
                }
                right_idx = mid - 1;
            }
            // cout << mid << endl;
            // sac = ref.substr(mid);
        }

        int most_idx;
        left_idx = least_idx;
        right_idx = csa.size() - 1;
        mid = (left_idx + right_idx) / 2;
        while (right_idx >= left_idx)
        {
            mid = (left_idx + right_idx) / 2;
            sac = ref.substr(csa[mid], p.length());
            cmp = p.compare(sac);
            if (cmp < 0)
            {
                right_idx = mid - 1;
            }
            else if (cmp > 0)
            {
                left_idx = mid + 1;
            }
            else if (cmp == 0)
            {
                if (mid == csa.size() - 1 || p.compare(ref.substr(csa[mid + 1], p.length())) < 0)
                {
                    most_idx = mid;
                    break;
                }
                left_idx = mid + 1;
            }
        }
        if (found)
        {
            matches[name] = {least_idx, most_idx};
        }
    }

    auto end = std::chrono::steady_clock::now();
    return {matches, std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()};
}

pair<unordered_map<string, pair<int, int>>, double> simpaccel(unordered_map<string, string> queries, csa_wt<> csa, string ref)
{
    unordered_map<string, pair<int, int>> matches;
    int least_idx,most_idx,left_idx,right_idx,mid,llcp,mlcp,rlcp,cmp,tlcp;
    bool found;
    string sac,p,name;
    auto start = std::chrono::steady_clock::now();
    // cout << "ref: " << ref.length() << endl;
    // int tlcp;
    for (auto const &pair : queries)
    {
        p = pair.second;
        // cout << "searching for " << p << endl;
        name = pair.first;
        matches[name] = {-1, -2};

        left_idx = 0;
        right_idx = csa.size();
        // mid = (left_idx + right_idx) / 2;
        llcp = lcp(ref,csa[left_idx],p);
        rlcp = lcp(ref,csa[right_idx],p);
        found = false;
        while (right_idx >= left_idx)
        {      
            
            // llcp = lcp(ref,csa[left_idx],p);
            // rlcp = lcp(ref,csa[right_idx],p);
            mlcp = min(llcp,rlcp);
            // cout << mlcp << endl;
            mid = (left_idx + right_idx) / 2;
            // sac = ref.substr(csa[mid]+mlcp, p.length()-mlcp);
            // cmp = (p.substr(mlcp)).compare(sac);
            // int neww = subcmp(ref,csa[mid]+mlcp,p,mlcp);
            cmp = subcmp(ref,csa[mid]+mlcp,p,mlcp);
            // if (cmp != neww)
            //     cout << "og" << cmp << "new" << neww << "diff" << cmp-neww << "s1" << ref.substr(csa[mid]+mlcp,p.length()+2) << "s2" << p.substr(mlcp)<< endl;
            if (cmp < 0)
            {
//                 if (mid == left_idx+1) {
//                     least_idx=mid;
//                     if (p.compare(ref.substr(csa[mid], p.length())) == 0)
// {                        found = true;
// }                    break;
//                 }
                right_idx = mid -1;
                // right_idx = mid - 1;
                // tlcp = min(llcp,abs(cmp)-1);
                // rlcp = lcp2(ref,csa[right_idx],p,0);
                // rlcp = lcp(ref.substr(csa[right_idx],p.length()),p);
                // rlcp = lcp(ref,csa[right_idx],p);
                rlcp = abs(cmp) - 1;
                // rlcp = abs(cmp)-1;
            }
            else if (cmp > 0)
            {
                // if (mid == right_idx-1) {
                //     least_idx=right_idx;
                //     if (p.compare(ref.substr(csa[right_idx], p.length())) == 0){
                //         found = true;}
                //     break;
                // }
                
                left_idx = mid + 1;
                // left_idx = mid+1 ;
                // llcp = abs(cmp)-1;
                // tlcp = min(rlcp,abs(cmp)-1);
                // llcp = lcp(ref,csa[left_idx],p);
                llcp = abs(cmp) - 1;
                
                // llcp = lcp2(ref,csa[left_idx]+0,p,0);
                // llcp = lcp(ref.substr(csa[left_idx],p.length()),p);
            }
            else if (cmp == 0)
            {
                // cout << "here" << endl;
                if (mid == 0 || p.compare(ref.substr(csa[mid - 1], p.length())) > 0)
                {
                    least_idx = mid;
                    found = true;
                    break;
                }
                
                // if (mid == right_idx-1) {
                //     if (p.compare(ref.substr(csa[right_idx], p.length())) == 0)
                //         found = true;
                //     least_idx=right_idx;
                //     break;
                // }
                right_idx = mid - 1;
                rlcp = p.length();
                // right_idx = mid-1 ;
                // rlcp = lcp(ref,csa[right_idx],p);
                
                // tlcp = min(llcp,abs(cmp)-1);
                // rlcp = lcp2(ref,csa[right_idx]+0,p,0);
                // rlcp = lcp(ref.substr(csa[right_idx],p.length()),p);
            }
            // cout << "li" << left_idx << "m" << mid << "ri" << right_idx << endl;
            // sac = ref.substr(mid);
        }

        if (!found) {
            continue;
        }

        int most_idx;
        left_idx = least_idx;
        right_idx = csa.size() - 1;
        llcp = lcp(ref,csa[left_idx],p);
        rlcp = lcp(ref,csa[right_idx],p);
        // mlcp = min(llcp,rlcp);
        mid = (left_idx + right_idx) / 2;
        while (right_idx >= left_idx)
        {   
            
            mlcp = min(llcp,rlcp);
            mid = (left_idx + right_idx) / 2;
            // sac = ref.substr(csa[mid], p.length());
            // cmp = p.compare(sac);
            cmp = subcmp(ref,csa[mid]+mlcp,p,mlcp);

            if (cmp < 0)
            {
                right_idx = mid - 1;
                rlcp = lcp(ref,csa[right_idx],p);
                
                // tlcp = min(llcp,abs(cmp)-1);
                // rlcp = lcp2(ref,csa[right_idx]+tlcp,p,tlcp);
            }
            else if (cmp > 0)
            {
                left_idx = mid + 1;
                llcp = lcp(ref,csa[left_idx],p);
            }
            else if (cmp == 0)
            {
                if (mid == csa.size() - 1 || p.compare(ref.substr(csa[mid + 1], p.length())) < 0)
                {
                    most_idx = mid;
                    break;
                }
                left_idx = mid + 1;
                llcp = lcp(ref,csa[left_idx],p);
            }
        }
        if (found)
        {
            matches[name] = {least_idx, most_idx};
        }
        // return{matches,0};
        
        if (most_idx < least_idx) {
            cout << "most_idx: " << most_idx << "least_idx: " << least_idx << endl;
        }
    }
    auto end = std::chrono::steady_clock::now();
    return {matches, std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()};
}

pair<unordered_map<string, pair<int, int>>, double> simpaccelpreft(unordered_map<string, string> queries, csa_wt<> csa, string ref, unordered_map<string, pair<int, int>> preft, int k)
{
    unordered_map<string, pair<int, int>> matches;
    int least_idx,most_idx,left_idx,right_idx,mid,llcp,mlcp,rlcp,cmp;
    bool found;
    string sac,p,name,pref;
    auto start = std::chrono::steady_clock::now();
    // cout << "ref: " << ref.length() << endl;
    for (auto const &pair : queries)
    {
        p = pair.second;
        // cout << "searching for " << p << endl;
        name = pair.first;
        matches[name] = {-1, -2};
        
        pref = p.substr(0, k);

        left_idx = preft[pref].first;
        right_idx = preft[pref].second;

        // left_idx = 0;
        // right_idx = csa.size();
        mid = (left_idx + right_idx) / 2;
         
        // cout << "li" << left_idx << "sidx" << csa[left_idx] << "str" << ref.substr(csa[left_idx]) << endl;
        // cout << "ri:" << right_idx << "sidx" << csa[right_idx] << "str"<< ref.substr(csa[right_idx]) << endl;
        llcp = lcp(ref,csa[left_idx],p);
        rlcp = lcp(ref,csa[right_idx],p);
        //  rlcp = lcp(ref.substr(csa[right_idx],p.length()),p);
         mlcp = min(llcp,rlcp);
        // string sac;
         found = false;
         cmp;
        while (right_idx >= left_idx)
        {      
            
            llcp = lcp(ref,csa[left_idx],p);
            rlcp = lcp(ref,csa[right_idx],p);
            mlcp = min(llcp,rlcp);
            // cout << mlcp << endl;
            mid = (left_idx + right_idx) / 2;
            sac = ref.substr(csa[mid]+mlcp, p.length()-mlcp);
            cmp = (p.substr(mlcp)).compare(sac);
            if (cmp < 0)
            {
                right_idx = mid - 1;
                // rlcp = lcp(ref.substr(csa[right_idx],p.length()),p);
            }
            else if (cmp > 0)
            {
                left_idx = mid + 1;
                // llcp = lcp(ref.substr(csa[left_idx],p.length()),p);
            }
            else if (cmp == 0)
            {
                if (mid == 0 || p.compare(ref.substr(csa[mid - 1], p.length())) > 0)
                {
                    least_idx = mid;
                    found = true;
                    break;
                }
                right_idx = mid - 1;
                // rlcp = lcp(ref.substr(csa[right_idx],p.length()),p);
            }
            // cout << mid << endl;
            // sac = ref.substr(mid);
        }

        int most_idx;
        left_idx = least_idx;
        right_idx = csa.size() - 1;
        mid = (left_idx + right_idx) / 2;
        while (right_idx >= left_idx)
        {
            mid = (left_idx + right_idx) / 2;
            sac = ref.substr(csa[mid], p.length());
            cmp = p.compare(sac);
            if (cmp < 0)
            {
                right_idx = mid - 1;
            }
            else if (cmp > 0)
            {
                left_idx = mid + 1;
            }
            else if (cmp == 0)
            {
                if (mid == csa.size() - 1 || p.compare(ref.substr(csa[mid + 1], p.length())) < 0)
                {
                    most_idx = mid;
                    break;
                }
                left_idx = mid + 1;
            }
        }
        if (found)
        {
            matches[name] = {least_idx, most_idx};
        }
        // return{matches,0};
    }
    auto end = std::chrono::steady_clock::now();
    return {matches, std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()};
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

    int k =0;
    if (preft.size() > 0)
    {
        k = preft.begin()->first.length();
    }
    // auto it = preft.begin();
    // string samp_key = it->first;
    // int k = samp_key.length();
    // for (int i = 0; i < csa.size(); i++)
    // {
    //     cout << i << " " << ref.substr(csa[i]) << endl;
    // }

    string q_fname = argv[2];
    cout << q_fname << endl;
    std::ifstream file(q_fname);
    string mode = argv[3];

    unordered_map<string, string> queries = unordered_map<string, string>();
    string prev_name;
    if (file.is_open())
    {
        // cout << q_fname << endl;
        std::string line;
        while (std::getline(file, line))
        {
            // using printf() in all tests for consistency
            if (line.at(0) == '>')
            {
                // refname
                prev_name = line.substr(1);
                // cout << qname << endl;
                // prev_name =
                queries[line.substr(1)] = "";
            }
            else
            {
                queries[prev_name].append(line);
            }
        }
        file.close();
    }
    else
    {
        cout << 34 << endl;
    }
    // unordered_map<string,pair<int,int>> matches = naive(queries,csa,ref);
pair<unordered_map<string, pair<int, int>>, double> res;
    if (mode.compare("naive") == 0)
    {
        if (preft.size() > 0)
        {
            cout << "naive preft" << endl;
            res = naivepreft(queries, csa, ref, preft, k);
        }
        else
        {
            
            cout << "naive no preft" << endl;
            res = naive(queries, csa, ref);
        }
    }
    else
    {
        if (preft.size() > 0)
        {
            cout << "simpaccel preft" << endl;
            // res = naivepreft(queries, csa, ref, preft, k);
        }
        else
        {
            
            cout << "simpaccel no preft" << endl;
            res = simpaccel(queries, csa, ref);
        }
    }

    string out_fname = argv[4];
    cout << out_fname << endl;
    cout << res.second << endl;
    ofstream ofile(out_fname);
    int nm;
    for (auto const &pair : res.first)
    {
        ofile << pair.first << "\t";
        nm = pair.second.second - pair.second.first + 1;
        ofile << nm;
        for (int j = pair.second.first; j <= pair.second.second; j++)
        {
            ofile << "\t" << csa[j];
        }
        ofile << endl;
        // std::cout << "{" << pair.first << ": " << pair.second.first << " " << pair.second.second << "}\n";
    }

    // std::string name = "13 dna:chromosome chromosome:GRCh38:13:1:114364328:1 REF";
}
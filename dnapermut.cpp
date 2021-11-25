#include <cstdio> //printf
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <inttypes.h>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace std;
using namespace seqan3;

const int MOD = 100003;

enum Braid { Single, Double, TripleLinear, TripleBranched, Longer };

void read_fasta(const char *filename , vector<string> * sequences)
{
    ifstream infile(filename);
    string line;
    string s1, s2;
    while (getline(infile, line)) {
        s1 = line.substr(1);
        getline(infile, line);
        s2 = line;
        sequences->push_back(s2);
    }
    infile.close();
}

void find_substrings(vector<string> * sequences, set<unsigned> * dependencies, set<unsigned> * relationships, set<unsigned> * expanded)
{
    fm_index index{*sequences};
    unsigned current_string = 0;
    

    // Iterate every string
    vector<string>::iterator it;
    for(it = sequences->begin(); it != sequences->end(); ++it)
    {
        // Match it against the index (exact matching)
        auto search_results = search(*it, index);

        // Add every match to the list of relationships 
        for (auto && hit : search_results)
        {
            // If its not a self-match or a multi-match with the same string
            if(hit.reference_id() != current_string && 
                std::find(relationships[current_string].begin(), relationships[current_string].end(), hit.reference_id()) == relationships[current_string].end() )
            {
                // dependencies holds the directed edges
                dependencies[current_string].emplace(hit.reference_id());

                // relationships holds dependencies and their reflexive relations
                relationships[current_string].emplace(hit.reference_id());
                relationships[hit.reference_id()].emplace(current_string);

                // expanded is a copy of relationships that will later include all the elements that belong to independent branches
                expanded[current_string].emplace(hit.reference_id());
                expanded[hit.reference_id()].emplace(current_string);

            }
            
        
        }
        ++current_string;
    }
}

void rec_expansion_sets(set<unsigned>::iterator it, set<unsigned> * current_row, 
                        set<unsigned> * relationships, set<unsigned> * visited, 
                        set<unsigned> * expanded, unsigned original)
{
    if(it == current_row->end()) return;

    if(*it != original) expanded[original].emplace(*it);
    
    if(!visited->contains(*it)){
        visited->emplace(*it);
        rec_expansion_sets(relationships[*it].begin(), &relationships[*it], relationships, visited, expanded, original);
    }
    rec_expansion_sets(++it, current_row, relationships, visited, expanded, original);
}

Braid detect_braid_type(unsigned element, set<unsigned> * dependencies, set<unsigned> * expanded)
{

    unsigned count = expanded[element].size();

    switch(count)
    {
        case 0: return Single;
        break;
        case 1: return Double;
        break;
        case 2: { 
                    // if 2 dependencies => tripleBranched, if 3 dependencies (transitive property) => TripleLinear
                    
                    unsigned intracount = dependencies[element].size();
                    set<unsigned>::iterator expit;
                    for(expit = expanded[element].begin(); expit != expanded[element].end(); ++expit)
                        intracount += dependencies[*expit].size();
                    
                    if(intracount == 2) return TripleBranched;
                    else if(intracount == 3) return TripleLinear;
                    else { cerr << "Found inconsistency in triple braid\n"; exit(-1); }
                }
        break;
        default: return Longer;
        break;
    }
}

uint64_t binomial(uint64_t n, uint64_t k)
{
    if(k > n)
        { cerr << "Bad input in binomial coefficient calculation\n"; exit(-1); }
    if(k == 0)
        return 1;
    if(k > n/2)
        return binomial(n, n-k);
    return n * binomial(n-1, k-1) / k;
}

uint64_t calculate_permutations(unsigned n_seqs, set<unsigned> * dependencies, set<unsigned> * expanded)
{
    uint64_t t_perm = 1;
    uint64_t remaining = (uint64_t) n_seqs;
    set<unsigned> * visited = new set<unsigned>;

    for(unsigned i=0; i<n_seqs; i++)
    {
        if(!visited->contains(i))
        {
            Braid b = detect_braid_type(i, dependencies, expanded);

            if(b == Single)
                { t_perm *= binomial(remaining, 1); remaining -= 1; }

            else if(b == Double)
                { t_perm *= binomial(remaining, 2); remaining -= 2; }

            else if(b == TripleLinear)
                { t_perm *= binomial(remaining, 3); remaining -= 3; }

            else if(b == TripleBranched)
                { t_perm *= 2 * binomial(remaining, 3); remaining -= 3; }

            else
            {
                // Generate permutations and filter them
            }

            // Add nodes from braid to visited so that we do not count them again
            visited->emplace(i);
            set<unsigned>::iterator it;
            for(it = expanded[i].begin(); it != expanded[i].end(); ++it)
                visited->emplace(*it);
        }

    }

    return t_perm;
}

void heap_rec_permutation(unsigned * a, unsigned size, unsigned n)
{
    if(size == 1)
    {
        //for(int i=0; i<n; i++) cout << a[i] << ",";
        //cout << endl;
        return;
    }

    for(unsigned i = 0; i < size; i++)
    {
        heap_rec_permutation(a, size - 1, n);

        if(size % 2 == 1)
            swap(a[0], a[size - 1]);
        else
            swap(a[i], a[size - 1]);
    }
}

void heap_iterate_permutation(unsigned * a, unsigned n)
{

    unsigned c[n];

    for(unsigned i=0; i<n; i++)
    {
        c[i] = 0;
        //cout << a[i] << ",";
    }
    //cout << endl;
    

    unsigned i = 0;
    while(i < n)
    {
        if(c[i] < i)
        {
            if(i % 2 == 0)
                { unsigned x = a[0]; a[0] = a[i]; a[i] = x; }
            else
                { unsigned x = a[c[i]]; a[c[i]] = a[i]; a[i] = x; }
            
            //for(unsigned j=0; j<n; j++) cout << a[j] << ",";
            //cout << endl;
            ++c[i];
            i = 0;
        }
        else
        {
            c[i] = 0;
            ++i;
        }
    }
}

int main(int argc, char **argv) {


    vector<string> sequences;
    if(argc != 2) { cerr << " Error: Use " << argv[0] << " <fasta>\n"; exit(-1); }
    read_fasta(argv[1], &sequences);


    // Calculate dependencies (substrings) between strings using an FM-index
    set<unsigned> * dependencies  = new set<unsigned>[sequences.size()];
    set<unsigned> * relationships = new set<unsigned>[sequences.size()];
    set<unsigned> * expanded      = new set<unsigned>[sequences.size()];
    find_substrings(&sequences, dependencies, relationships, expanded);

    for(unsigned i=0; i<sequences.size(); i++){
        set<unsigned>::iterator it;
        cout << "### Seq " << i << " is related to:\n";
        for(it = relationships[i].begin(); it != relationships[i].end(); ++it){
            cout << *it << " ";
        }
        cout << "\n";
    }


    // Expand relationships to find lose braids
    for(unsigned i=0; i<sequences.size(); i++)
    {
        set<unsigned> * visited = new set<unsigned>;
        set<unsigned>::iterator it = relationships[i].begin();
        visited->emplace(i);
        rec_expansion_sets(it, &relationships[i], relationships, visited, expanded, i);
    }

    cout << "Expansions: \n";
    for(unsigned i=0; i<sequences.size(); i++){
        set<unsigned>::iterator it;
        cout << "### Seq " << i << " is expanded to:\n";
        for(it = expanded[i].begin(); it != expanded[i].end(); ++it){
            cout << *it << " ";
        }
        cout << "\n";
    }
    
    cout << "\n\n\nBRAIDS\n";
    
    for(unsigned i=0; i<sequences.size(); i++){
        cout << "Seq " << i << " is type: " << detect_braid_type(i, dependencies, expanded) << endl;
    }
    
    uint64_t t_perm = calculate_permutations(sequences.size(), dependencies, expanded);

    cout << "Total permutations: " << t_perm << endl;

    delete [] relationships;
    delete [] dependencies;
    delete [] expanded;

    /*
    vector<sequence> sequences;
    int ans;
    read_fasta(argv[1], sequences);
    // Good luck and have fun!
    // ans = your code
    printf("Number of different ways to sort the fasta file (mod 100003) is: %d\n", ans);
    */
    return 0;
}
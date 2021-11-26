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

const uint64_t MOD = 100003;

// Enum identifying each of the possible groupings within braids
enum Braid { Single, Double, TripleLinear, TripleBranched, Longer };

/* 
* Function read_fasta
* Reads an input fasta file where each sequence has one numeric identifier and one string without endlines
* @param filename          String name containing the path of the fasta file
* @param sequences         A pointer to an std::vector in which only the dna strings will be stored
* @returns                 Nothing
*/
void read_fasta(const char * filename , vector<string> * sequences)
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

/* 
* Function find_substrings
* Uses an FM-index to compute matches and partial matches (substrings) between query and reference. In this
* particular case, the query and the reference are the same data.
* @param sequences         Pointer to std::vector containing the sequences
* @param dependencies      Pointer to an array of std::set of size n_seqs, such that each element contains the id of the strings for which it is a substring
* @param relationships     Same as above but including reflexive dependencies (if seq 0 comes before seq 2, then 2 is also related to 0)
* @param expanded          Copy of relationships for posterior use
* @returns                 Nothing
*/
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

/* 
* Function rec_expansion_sets
* A recursive function that expands dependencies, i.e. given a list of reflexive relationships, rec_expansion_sets will
* connect together all elements (transitive property) that share a dependency, such that if seq 0 comes before 2 and 3, 
* and seq 2 comes before 4, then seq 0, 2, 3 and 4 are related.
* @param it                An iterator to the currently in-expansion std::set
* @param current_row       Point to a row (with reflexive dependencies) where the current iterator is traversing
* @param relationships     Pointer to the table containing all reflexive dependencies
* @param visited           A pointer to an std::set containing which sequences have already been visited (prevents cycles)
* @param expanded          Pointer to the table containing all expanded relationships (transitive)
* @param original          The current sequence for which the expansion has been launched
* @returns                 Nothing
*/
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

/* 
* Function detect_braid_type
* Given an expanded relationships table, this function returns the type of braid (or grouping) for a given sequence. 
* For instance, if the sequence is related to 2 other sequences, detect_braid_type will return the type of braid which can
* be either 3 sequences with linear relationships (TripleLinear) or 2 sequences related to 1, or viceversa (TripleBranched)
* @param element           The current sequence for which the grouping or braid is being detected
* @param dependencies      Pointer to the table containing all dependencies
* @param expanded          Pointer to the table containing all expanded relationships (reflexive and transitive)
* @returns                 The Braid type up to 4 different types
*/
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
                    // All cases are straight-forward except for the triple one
                    // Which can be linear (1->2->3) or branch-like (1->2, 1->3 or 1->3, 2->3)
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

/* 
* Function mod_of_mult
* A recursive function to calculate the modulo of a large multiplication ((a*b) % m). It proceeds by splitting each
* multiplication into smaller parts that can be added together while performing the modulo at the same time.
* @param a                 The first multiplication operand
* @param b                 The second multiplication operand
* @param m                 The modulo operand
* @returns                 The modulo m of the multiplication of a times b
*/
uint64_t mod_of_mult(uint64_t a, uint64_t b, uint64_t m) 
{ 
    uint64_t res = 0;
    a = a % m;
    while(b > 0) 
    { 
        if (b % 2 == 1) 
            res = (res + a) % m; 
        a = (a * 2) % m; 
        b /= 2; 
    } 
    return res % m; 
}

/* 
* Function binomial
* A recursive function to calculate the binomial coefficient efficiently. It takes advantage of the recursive formula:
* (n, k) = (n / k) * (n-1, k-1) which enables to calculate combinations much quicker while avoiding overflow.
* @param n                 The number of items to choose from
* @param k                 The number of items chosen
* @returns                 The number of combinations
*/
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

/* 
* Function filter_permutation
* A function that filters permutations (when generated with Heap's algorithm or from a different source) 
* which only requires O(2p) time to accept or discard a permutation (where p is the length of the permutation)
* @param a                 Pointer to an array containing the permutation
* @param n                 The length of the permutation
* @param dependencies      Pointer to the table containing all dependencies
* @param n_seqs            The number of sequences (used for translation into lookup table)
* @returns                 True or False depending if the permutation evaluated is accepted or discarded
*/
bool inline filter_permutation(unsigned * a, unsigned n, set<unsigned> * dependencies, unsigned n_seqs)
{
    // Create key-value dictionary indexed by the string where value is the position in the sequence
    // This way we only have to traverse the permutation twice: once now and the other to check dependencies
    unsigned keypos[n_seqs]; // Note: right now translation into key-pos is done using O(n_seqs) memory. 
    for(unsigned i=0; i<n; i++) keypos[a[i]] = i;

    for(unsigned i=0; i<n; i++)
    {
        unsigned current = a[i];
        set<unsigned>::iterator deps;
        for(deps = dependencies[current].begin(); deps != dependencies[current].end(); ++deps)
            if(keypos[current] > keypos[*deps])
                return false;
    }
    return true;
}

/* 
* Function heap_iterate_permutation_filtering
* A function that generates all possible permutations of a given array and that filters them accordingly to a set
* of dependencies. This function is to be used only when the braid can not be computed analyitically.
* @param a                 Pointer to an array containing the permutation
* @param n                 The length of the permutation
* @param dependencies      Pointer to the table containing all dependencies
* @param n_seqs            The number of sequences (used for translation into lookup table)
* @returns                 The number of accepted permutations
*/
uint64_t heap_iterate_permutation_filtering(unsigned * a, unsigned n, set<unsigned> * dependencies, unsigned n_seqs)
{

    unsigned c[n];
    uint64_t t_perm = 0;

    for(unsigned i=0; i<n; i++)
        c[i] = 0;

    if(filter_permutation(a, n, dependencies, n_seqs))
        ++t_perm;

    unsigned i = 0;
    while(i < n)
    {
        if(c[i] < i)
        {
            if(i % 2 == 0)
                { unsigned x = a[0]; a[0] = a[i]; a[i] = x; }
            else
                { unsigned x = a[c[i]]; a[c[i]] = a[i]; a[i] = x; }
            
            if(filter_permutation(a, n, dependencies, n_seqs))
                ++t_perm;
            
            ++c[i];
            i = 0;
        }
        else
        {
            c[i] = 0;
            ++i;
        }
    }
    return t_perm;
}

/* 
* Function calculate_permutations
* This function computes the number of permutations by making use of the principle of multiplication. This means
* that all braids will be calculated independently and then multiplied, thus reducing complexity from O(n!) to 
* O(k1! + k2! + k3! + ... kj!). Additionally, braids of size < 4 will be computed analytically in closed form without
* requiring generating any permutation at all. 
* @param n_seqs            The number of sequences
* @param dependencies      Pointer to the table containing all dependencies
* @param expanded          Pointer to the table containing all expanded relationships (reflexive and transitive)
* @returns                 The total number of permutations modulo 100003
*/
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
            {
                t_perm = mod_of_mult(t_perm, binomial(remaining, 1), MOD);
                remaining -= 1; 
            }

            else if(b == Double)
            {
                t_perm = mod_of_mult(t_perm, binomial(remaining, 2), MOD); 
                remaining -= 2; 
            }

            else if(b == TripleLinear)
            {
                t_perm = mod_of_mult(t_perm, binomial(remaining, 3), MOD); 
                remaining -= 3; 
            }

            else if(b == TripleBranched)
            {
                t_perm = mod_of_mult(t_perm, 2 * binomial(remaining, 3), MOD); 
                remaining -= 3; 
            }

            else
            {
                // Generate permutations and filter them
                unsigned n = (unsigned) expanded[i].size() + 1;
                unsigned a[n];
                a[0] = i;
                set<unsigned>::iterator it;
                unsigned j = 1;
                for(it = expanded[i].begin(); it != expanded[i].end(); ++it)
                {
                    a[j] = *it;
                    ++j;
                }
                t_perm = mod_of_mult(t_perm, heap_iterate_permutation_filtering(a, n, dependencies, n_seqs), MOD);
                remaining -= (uint64_t) n;
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

/* 
* Function main
* This program takes as input a fasta file and generates all permutations following the substring constraint, i.e.
* a sequence that is subsequence of another one, must appear before in all permutations. To do so, this program 
* first identifies which sequences are related to others and isolates them in "braids" or "groupings". Once so, 
* each braid or group can be calculated independently (thus reducing complexity by a factor of the size of each group).
* In fact, some of these groups are computed in closed form, meaning that permutations do not need to be calculated.
* @param argc              System arg - number of arguments
* @param argv              System arg - actual arguments
* @returns                 0 if program exited correctly, -1 otherwise
*/
int main(int argc, char **argv) {


    vector<string> sequences;
    if(argc != 2) { cerr << " Error: Use " << argv[0] << " <fasta>\n"; exit(-1); }
    read_fasta(argv[1], &sequences);


    // Calculate dependencies (substrings) between strings using an FM-index
    set<unsigned> * dependencies  = new set<unsigned>[sequences.size()];
    set<unsigned> * relationships = new set<unsigned>[sequences.size()];
    set<unsigned> * expanded      = new set<unsigned>[sequences.size()];
    find_substrings(&sequences, dependencies, relationships, expanded);

    // Writes dependencies in graph format (V,V)
    unsigned n_substrings = 0;
    for(unsigned i=0; i<sequences.size(); i++){
        set<unsigned>::iterator it;
        for(it = dependencies[i].begin(); it != dependencies[i].end(); ++it){
            //cout << i << " " << *it << endl; // Enable this print fro graph-like output
            ++n_substrings;
        }
    }

    // Expand relationships to find lose braids
    for(unsigned i=0; i<sequences.size(); i++)
    {
        set<unsigned> * visited = new set<unsigned>;
        set<unsigned>::iterator it = relationships[i].begin();
        visited->emplace(i);
        rec_expansion_sets(it, &relationships[i], relationships, visited, expanded, i);
    }

    // Finds longest chain and prints expansion table
    unsigned longest_chain = 0;
    for(unsigned i=0; i<sequences.size(); i++){
        if((expanded[i].size() + 1) > longest_chain)
            longest_chain = expanded[i].size() + 1;
    }
    
    // Prints the type of braids detected (commented for clarity)
    /*
    cout << "Braid type: \n";
    for(unsigned i=0; i<sequences.size(); i++){
        cout << "Sequence " << i << " is of type: " << detect_braid_type(i, dependencies, expanded) << endl;
    }
    */
    
    uint64_t t_perm = calculate_permutations(sequences.size(), dependencies, expanded);

    cout << "Number of substrings:              " << n_substrings << endl;
    cout << "Longest substring chain:           " << longest_chain << endl;
    cout << "Total permutations (mod 100003):   " << t_perm << endl;

    delete [] relationships;
    delete [] dependencies;
    delete [] expanded;
    
    return 0;
}
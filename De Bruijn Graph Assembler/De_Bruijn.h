//De Bruijn header file
//Jacob Gerace

using namespace std;
#include <stdio.h>
#include <iostream>
#include <unordered_map>
#include <list>
#include <vector>
#include <set>
#ifndef De_Bruijn_H
#define De_Bruijn_H

//Structure for nodes in the graph
struct DBNode {
        string kmer;
        set<string> out;
        set<string> in;
        bool branch;
        bool marker;
        int count;
};

//comments for each funciton found in the .cpp implementation
class De_Bruijn
{
        public:
                De_Bruijn(int new_ksize);
                ~De_Bruijn();
                void add_sequence(string sequence);
                bool good_read(string sequence);
                vector<string> assemble_contigs();
                int get_num_nodes();

        private:
                vector<string> make_kmers(string sequence);
                void add_node(string kmer, string out, string in,
                              bool branch, bool marker);

                DBNode* make_DBNode(string key, string out, string in,
                                    bool branch, bool marker, int count);
                void mark_branching();
                vector<DBNode*> find_start_nodes();
                vector<DBNode*> find_post_branch_nodes();
                string get_linear_seq(string key);
                string make_string(vector<string> kmerlist);
                void print();

                int ksize;
                int num_nodes;
                //mapped by the key: a string, and a pointer to the node
                unordered_map<string, DBNode*> db_graph;
};
#endif

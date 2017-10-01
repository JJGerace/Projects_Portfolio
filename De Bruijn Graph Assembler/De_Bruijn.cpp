//implementation of de bruijn graph
//Jacob Gerace
#include "De_Bruijn.h"

//construct graph without nodes
De_Bruijn::De_Bruijn(int new_ksize)
{
        num_nodes = 0;
        ksize = new_ksize;
}

//only have to delete the DBNodes, the rest will free automatically
//Due to the C++ standard library
De_Bruijn::~De_Bruijn()
{
        for (auto it = db_graph.begin(); it != db_graph.end(); it++)
                delete it->second;
}

//Add a single sequence to the graph
void De_Bruijn::add_sequence(string sequence)
{
        //Breaks into kmers
        vector<string> kmerarray = make_kmers(sequence);
        unsigned length = kmerarray.size();

        //add each kmer to the graph
        if (length == 1) {
                add_node(kmerarray[0], "", "", false, false);
        }
        else if (length > 1) {
                add_node(kmerarray[0], kmerarray[1], "", false, false);
                for (unsigned i = 1; i < kmerarray.size() - 1; i++) {
                        add_node(kmerarray[i], kmerarray[i + 1],
                                 kmerarray[i - 1], false, false);
                }
                add_node(kmerarray[length - 1], "", kmerarray[length - 2],
                         false, false);
        }
}

//Checks to see if each sequence is a good read
//The sequence must have previously been added to the graph
//Good reads are those that have no nodes with a multiplicity
//of 1.
bool De_Bruijn::good_read(string sequence)
{
        if (sequence.length() < (unsigned)ksize)
                return false;


        //make the kmers in order to look up the nodes in the graph
        vector<string> kmerarray = make_kmers(sequence);
        DBNode* lookup;
        unordered_map<string, DBNode*>::iterator it;

        for (unsigned i = 0; i < kmerarray.size(); i++) {
                it = db_graph.find(kmerarray[i]);
                if (it != db_graph.end()) { //exists in graph
                        lookup = it->second;
                        //if multiplicity of 1, then bad read
                        if(lookup->count == 1) {
                                return false;
                        }
                }
        }
        return true;
}

//Heart and soul of the project:
//Returns an array with all of the linear reigons of the graph (contigs)
vector<string> De_Bruijn::assemble_contigs() {
        //Identify the branching nodes
        mark_branching();
        vector<DBNode*> start_nodes = find_start_nodes();
        vector<string> assembled;
        string temp = "";

        //for every startnode, walk along the graph using
        //get_linear_seq and add each string ot the array
        for(unsigned i = 0; i < start_nodes.size(); i++) {
                if(start_nodes[i]->marker == false) {
                        temp = get_linear_seq(start_nodes[i]->kmer);
                        if (temp != "")
                                assembled.push_back(temp);
                }
        }
        temp = "";

        //For each node that comes from a branch node, also walk
        //along the graph using get_linear_seq and add that string
        //to the array
        vector<DBNode*> post_branch_nodes = find_post_branch_nodes();
        for(unsigned i = 0; i < post_branch_nodes.size(); i++) {
                if(post_branch_nodes[i]->marker == false) {
                        temp = get_linear_seq(post_branch_nodes[i]->kmer);
                        if (temp != "")
                                assembled.push_back(temp);
                }
        }
        return assembled;
}

//accessor
int De_Bruijn::get_num_nodes()
{
        return num_nodes;
}

//add node with the relevant info to the graph
void De_Bruijn::add_node(string kmer, string out, string in,
                         bool branch, bool marker)
{
        unordered_map<string, DBNode*>::iterator it = db_graph.find(kmer);
        if(it == db_graph.end()) {//idiom for key not found
                DBNode* new_node = make_DBNode(kmer, out, in,
                                               branch, marker, 1);
                db_graph.insert({kmer, new_node});
        }
        else {
                DBNode *lookup = db_graph.at(kmer);
                if (out != "")
                        lookup->out.insert(out);
                if (in != "")
                        lookup->in.insert(in);

                lookup->count++;
        }
        num_nodes++;

}

//return "" if invalid string formed
//marks all it visits
string De_Bruijn::get_linear_seq(string key)
{
        string kmer = key;
        unordered_map<string, DBNode*>::iterator it;
        set<string>::iterator set_it;
        vector<string> kmerlist;
        DBNode* lookup;

        //infinite loop, but every case in the loop has a exit condition
        while(true){
                //lookup the current key
                it = db_graph.find(kmer);
                if (it == db_graph.end()) {
                        //something went wrong or an incorrect key was given
                        return "";
                }
                else {
                        lookup = it->second;
                        kmer = lookup->kmer;
                        //stop when you reach branched
                        if (lookup->branch) {
                                lookup->marker = true;
                                kmerlist.push_back(kmer);
                                return make_string(kmerlist);
                        }
                        //discard when you reach marker
                        else if (lookup->marker) {
                                return "";
                        }
                        //walk along graph otherwise
                        else {
                                lookup->marker = true;
                                kmerlist.push_back(kmer);
                                set_it = lookup->out.begin();

                                //reached end of list
                                if (set_it == lookup->out.end())
                                        return make_string(kmerlist);

                                kmer = *(lookup->out.begin());
                        }
                }
        }
}

//Given a list of adjacent kmers, compress them into a string
string De_Bruijn::make_string(vector<string> kmerlist)
{
        (void)kmerlist;
        if (kmerlist.size() == 0)
                return "";

        string contig = kmerlist[0];
        for (vector<string>::iterator it = kmerlist.begin();
                                      it != kmerlist.end(); ++it) {
                if (it != kmerlist.begin()) {
                        contig += (*it)[ksize - 1];
                }

        }

        return contig;
}

//Make a DBNode from given info. Completely abstract from the graph
//implementation
DBNode* De_Bruijn::make_DBNode(string key, string out, string in,
                               bool branch, bool marker, int count)
{
        DBNode* new_node = new DBNode;
        new_node->kmer = key;
        if (out != "")
                new_node->out.insert(out);

        if (in != "")
                new_node->in.insert(in);

        new_node->branch = branch;
        new_node->marker = marker;
        new_node->count = count;
        return new_node;
}

//Split a sequence up into kmers
vector<string> De_Bruijn::make_kmers(string sequence)
{
        vector<string> new_kmers;
        unsigned start;
        for (start = 0; start < ((sequence.length() - ksize) + 1); start++) {
                new_kmers.push_back(sequence.substr(start, ksize));
        }
        return new_kmers;
}

//Iterate through all nodes on the graph, if the indegree or outdegree
//is greater than one then mark them as branching
void De_Bruijn::mark_branching()
{
        unordered_map<string, DBNode*>::iterator graph_it;
        DBNode* lookup;

        for (graph_it = db_graph.begin();
             graph_it != db_graph.end(); ++graph_it) {
                lookup = graph_it->second;
                if ((lookup->out.size() > 1) || (lookup->in.size() > 1))
                        lookup->branch = true;
                else
                        lookup->branch = false;
        }
}

//Iterate through all nodes on the graph, return an array of all the nodes
//That have indegree = 0
vector<DBNode* > De_Bruijn::find_start_nodes()
{
        vector<DBNode*> nodes;
        unordered_map<string, DBNode*>::iterator graph_it;
        DBNode* lookup;

        for (graph_it = db_graph.begin();
             graph_it != db_graph.end(); ++graph_it) {
                lookup = graph_it->second;
                if (lookup->in.size() == 0) {
                        nodes.push_back(lookup);
                }
        }
        return nodes;
}

//Iterate through all nodes on the graph, return an array of all the nodes
//That come directly after a branching node, and only a branching node
vector<DBNode*> De_Bruijn::find_post_branch_nodes()
{
        vector<DBNode*> nodes;
        unordered_map<string, DBNode*>::iterator graph_it;
        DBNode* lookup;
        DBNode* previous;
        string previous_key;

        for (graph_it = db_graph.begin();
             graph_it != db_graph.end(); ++graph_it) {
                //get each node with indegree 1 whose incoming edge comes from
                //a branching node
                lookup = graph_it->second;
                if (lookup->in.size() == 1) {
                        previous_key = *lookup->in.begin();
                        previous =  db_graph.find(previous_key)->second;
                        if (previous->branch == true) {
                                nodes.push_back(lookup);
                        }
                }
        }
        return nodes;


}

//print function, for debug purposes
void De_Bruijn::print()
{
        unordered_map<string, DBNode*>::iterator graph_it;

        for (graph_it = db_graph.begin();
             graph_it != db_graph.end(); ++graph_it) {

                DBNode* lookup = graph_it->second;
                cout << "key: " << lookup->kmer << endl;
                cout << "       branch: " << lookup->branch << endl;
                cout << "       marker: " << lookup->marker << endl;
                cout << "       Outgoing:";

                for (set<string>::iterator it = lookup->out.begin();
                                            it != lookup->out.end(); ++it) {
                        cout << endl << "               " << *it;
                }
                cout << endl;
                cout << "       Incoming:";


                for (set<string>::iterator it = lookup->in.begin();
                                            it != lookup->in.end(); ++it) {
                        cout << endl << "               " << *it;
                }

                cout << endl;
                cout << "       count: " << lookup->count << endl;
        }
}


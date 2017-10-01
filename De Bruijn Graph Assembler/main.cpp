//main.cpp
//Jacob Gerace
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include "De_Bruijn.h"
#include <algorithm>
using namespace std;
const int ksize = 31;
const int print_threshold = 100;

//for using std:sort at the end, note that sorting is optional
//with this code
void print(vector<string> list, vector<string> reads);
bool comparison_function (int i, int j) { return (i > j); }
int calulate_n50(vector<int> sorted_lengths);

//Main method, reads in data from standard in, uses the two
//de_bruijn graphs to get the contig scaffolds
int main(int argc, char* argv[])
{

        /*read in from stdin or a file */
        vector<string> list;
        ifstream in_stream;
        string line;
        if (argc == 1) {
                while (!cin.eof()) {
                        cin >> line;
                        list.push_back(line);
                }
        }
        else {
                in_stream.open(argv[1]);
                while (!in_stream.eof()) {
                        in_stream >> line;
                        list.push_back(line);
                }
        }
        in_stream.close();
        De_Bruijn* first_graph = new De_Bruijn(ksize);
        /* Put all read sequences into first DeBrujin graph*/
        for (unsigned i = 0; i < (list.size() - 1); i++) {
                if (list[i].length() < ksize) {
                        cerr << "ignoring input line: "
                             << list[i]
                             << " it is below the ksize of "
                             << ksize
                             << endl;
                        list[i] = "";
                }
                else
                        first_graph->add_sequence(list[i]);
        }
        if (first_graph->get_num_nodes() == 0)
                return 1;


        //erase all sequences that are bad reads
        for (unsigned i = 0; i < (list.size() - 1); i++) {
                if (list[i] != ""
                    && !first_graph->good_read(list[i])) {
                        list[i] = "";
                }
        }
        delete first_graph;

        //reconstruct graph with good reads
        De_Bruijn* second_graph = new De_Bruijn(ksize);
        for (unsigned i = 0; i < (list.size() - 1); i++) {
                if (list[i] != "") {
                        second_graph->add_sequence(list[i]);
                }
        }
        //get and report output
        vector<string> assembled_seqs = second_graph->assemble_contigs();
        delete second_graph;
        print(assembled_seqs, list);
}

//prints output_contigs from contigs
//prints contig_lengths from contigs
//prints good_reads     from reads
void print(vector<string> contigs, vector<string> reads)
{
        //print assembled graph to output_contigs
        //print lengths to contig_lengths
        //print good reads to good_reads
        ofstream file;
        file.open("output_contigs");
        for(unsigned i = 0; i < contigs.size(); i++) {
                if (contigs[i].length() >= print_threshold)
                        file << contigs[i] << endl;
        }
        file.close();
        file.open("contig_lengths");
        /*IF YOU WANT TO SORT, UNCOMMENT THIS SECTION
        vector<int> sorted_lengths;

        for(unsigned i = 0; i < contigs.size(); i++) {
                if (contigs[i].length() >= print_threshold)
                        sorted_lengths.push_back((int)contigs[i].length());
        }
        std::sort(sorted_lengths.begin(), sorted_lengths.end(),
                  comparison_function);


        for(unsigned i = 0; i < sorted_lengths.size(); i++) {
                file << sorted_lengths[i] << endl;
        }

        UNCOMMENT ENDING HERE*/

         /*Can automatically calculate N50 if you want:
        file << "N50 size: " << calulate_n50(sorted_lengths) << endl;
         */


        /* COMMENT OUT STARTING HERE TO SORT */
        for(unsigned i = 0; i < contigs.size(); i++) {
                if (contigs[i].length() >= print_threshold)
                        file << contigs[i].length() << endl;
        }
        /*COMMENT OUT ENDING HERE TO SORT */
        file.close();

        file.open("good_reads");
        for(unsigned i = 0; i < (reads.size() - 1); i++) {
                if (reads[i] != "")
                        file << reads[i] << endl;
        }
        file.close();

}

//calculates the n50 value, part of the optional output
int calulate_n50(vector<int> sorted_lengths)
{
        int sum = 0;
        for(unsigned i = 0; i < sorted_lengths.size(); i++) {
                sum += sorted_lengths[i];
        }

        int sum_largest = 0;
        double fraction = 0;
        for(unsigned i = 0; i < sorted_lengths.size(); i++) {
                sum_largest += sorted_lengths[i];
                fraction = (double)sum_largest / (double)sum;
                if (fraction >= .5)
                        return sorted_lengths[i];
        }
        //something went wrong
        return -1;
}



#include <iostream>
#include <set>
#include <string>
#include <fstream>
#include <cstdlib>
using namespace std;

// written by Pierre Lindenbaum (https://www.biostars.org/p/10353/)

int main(int argc,char **argv)
    {
    string line;
    set<string> names;
    /* reads your names */
    ifstream in(argv[1],ios::in);
    if(!in.is_open()) return EXIT_FAILURE;
    while(getline(in,line,'\n'))
        {
        names.insert(line);
        }
    in.close();
    /* reads the fastq */
    string name;
    string seq;
    string name2;
    string qual;
    set<string>::iterator r;
    /* assuming 4 lines / read */
    while(!names.empty())
        {
        if(!getline(cin,name,'\n')) break;
        if(!getline(cin,seq,'\n')) break;
        if(!getline(cin,name2,'\n')) break;
        if(!getline(cin,qual,'\n')) break;
        r=names.find(name);
        if(r==names.end()) continue;
        names.erase(r);//names are unique we don't need this one anymore
        cout << name << "\n" << seq << "\n" << name2 << "\n" << qual << endl;
        }
    return 0;
    }

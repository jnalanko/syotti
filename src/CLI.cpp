#include <iostream>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <iterator>
#include <cassert>
#include <algorithm>
#include <set>
#include <map>
#include <sstream>
#include <fstream>
#include "commands.hh"

using namespace std;

static vector<string> commands = {"index", "design", "examine", "fill-gaps", "subset"};

void print_help(int argc, char** argv){
    (void) argc; // Unused parameter
    cerr << "Available commands: " << endl;
    for(string S : commands) cerr << "   " << argv[0] << " " << S << endl;
    cerr << "Running a command without arguments prints the usage instructions for the command." << endl;
}

int main(int argc, char** argv){

    if(argc == 1){
        print_help(argc, argv);
        return 1;
    }

    string command = argv[1];
    if(command == "--help" || command == "-h"){
        print_help(argc, argv);
        return 1;
    }

    // Drop the first element of argv
    for(int64_t i = 1; i < argc; i++) argv[i-1] = argv[i];
    argc--;

    try{
        if(command == "index") return index_main(argc, argv);
        else if(command == "design") return design_main(argc, argv);
        else if(command == "examine") return examine_main(argc, argv);
        else if(command == "fill-gaps") return fill_gaps_main(argc, argv);
        else if(command == "subset") return subset_main(argc, argv);
        else{
            throw std::runtime_error("Invalid command: " + command);
            return 1;
        }
    } catch (const std::runtime_error &e){
        std::cerr << "Runtime error: " << e.what() << '\n';
        return 1;
    } catch(const std::exception& e){
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }

}
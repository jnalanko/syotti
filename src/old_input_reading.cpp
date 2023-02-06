#include "SeqIO.hh"

// Vector of (read, header) pairs
std::vector<std::pair<std::string, std::string> > parse_FASTA(std::string filename){
    throwing_ifstream input(filename);

    std::vector<std::pair<std::string,std::string> > reads;

    std::string line;
    while(input.getline(line)){
        while(line.size() > 0 && isspace(line.back()))
            line.pop_back(); // Trim trailing whitespace just in case

        if(line.size() == 0 && !(input.stream.eof())) continue; // Ignore empty lines, just in case

        if(line[0] == '>')
            reads.push_back({"", line}); // Start new read
        else 
            reads.back().first += line; // Append base pairs to read
    }

    //Erase invalid reads
    //reads.erase(std::remove_if(reads.begin(), reads.end(), is_invalid_read), reads.end());
    return reads;
}


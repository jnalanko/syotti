#pragma once

/*
  Buffered reading for FASTA and FASTQ files.
  Authors: Jarno Alanko & Simon Puglisi
*/

#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <memory>
#include "throwing_streams.hh"
#include "buffered_streams.hh"

using namespace std;

namespace SeqIO{

enum Format {FASTA, FASTQ};

struct FileFormat{
    Format format;
    bool gzipped;
    string extension; // Includes the possible .gz extension
};

FileFormat figure_out_file_format(string filename);

class NullStream : public std::ostream {
public:
  NullStream() : std::ostream(nullptr) {}
};

template <class T>
const NullStream &operator<<(NullStream &&os, const T &value) { 
  return os;
}

void reverse_complement_c_string(char* S, int64_t len);

template<typename ifstream_t = Buffered_ifstream<std::ifstream>> // The underlying file stream.
class Reader {

// The class is used like this:
// Sequence_Reader_Buffered sr;
// while(true) { 
//   int64_t len = sr.get_next_read_to_buffer();
//   if(len == 0) break;
//   do something with sr.read_buf
//}
//
// or (slow):
// while(true) { 
//   read = sr.get_next_read()
//   if(read.size() == 0) break;
//}

private:

Reader(const Reader& temp_obj) = delete; // No copying
Reader& operator=(const Reader& temp_obj) = delete;  // No copying

std::unique_ptr<ifstream_t> stream;
int64_t mode;
int64_t read_buf_cap;
int64_t header_buf_cap;

bool reverse_complements = false; // Whether reverse complements are enabled
bool return_rc_next = false; // If reverse complements are enabled, this flag is used internally to manage the process
string filename;

vector<char> rc_buf; // Internal buffer for reverse complements

public:

    // These buffers are intended to be read from outside the class
    char* read_buf; // Stores a sequence read
    char* header_buf; // Stores the header of a read (without the '>' or '@')

    void read_first_char_and_sanity_check(){
        
        char c = 0; stream->get(&c);
        if(mode == FASTA && c != '>')
            throw runtime_error("ERROR: FASTA file does not start with '>'");
        if(mode == FASTQ && c != '@')
            throw runtime_error("ERROR: FASTQ file does not start with '@'");

        // This leaves the input stream pointer after the first character, but
        // get_next_read_to_buffer is written such that it's ok.
    }

    // mode should be FASTA_MODE or FASTQ_MODE
    // Note: FASTQ mode does not support multi-line FASTQ
    Reader(string filename, int64_t mode) : mode(mode), filename(filename) {
        stream = std::make_unique<ifstream_t>(filename, ios::binary);
        if(mode != FASTA && mode != FASTQ)
            throw std::invalid_argument("Unkown sequence format");
        
        read_buf_cap = 256;
        read_buf = (char*)malloc(read_buf_cap);

        header_buf_cap = 256;
        header_buf = (char*)malloc(header_buf_cap);

        read_first_char_and_sanity_check();
    }

    Reader(string filename) : filename(filename) {
        stream = std::make_unique<ifstream_t>(filename, ios::binary);
        SeqIO::FileFormat fileformat = figure_out_file_format(filename);
        if(fileformat.format == FASTA) mode = FASTA;
        else if(fileformat.format == FASTQ) mode = FASTQ;
        else throw(runtime_error("Unknown file format: " + filename));

        read_buf_cap = 256;
        read_buf = (char*)malloc(read_buf_cap);

        header_buf_cap = 256;
        header_buf = (char*)malloc(header_buf_cap);

        read_first_char_and_sanity_check();
    }

    void enable_reverse_complements(){
        reverse_complements = true;
        return_rc_next = false;
    }


    ~Reader(){
        free(read_buf);
        free(header_buf);
    }

    void rewind_to_start(){
        // Create a new stream
        stream = std::make_unique<ifstream_t>(filename, ios::binary);

        read_first_char_and_sanity_check();
        return_rc_next = false;
    }

    int64_t get_mode() const {return mode;}

    void fix_alphabet(char* buf, int64_t len){
        for(int64_t i = 0; i < len; i++){
            char c = buf[i];
            if(c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N'){
                cerr << "Warning: non-ACGTN-character found: '" << c << "'. Replacing with 'N'" << endl; 
                buf[i] = 'N';
            }
        }
    }

    // Returns length of read, or zero if no more reads.
    // The read is null-terminated.
    // The read is stored in the member pointer `read_buffer`
    // The header is stored in the member pointer `header buffer`
    // When called, the read that is currently in the buffer is overwritten
    int64_t get_next_read_to_buffer() {
        if(reverse_complements){
            if(return_rc_next){
                // Copy from rc_buf to read_buf
                strcpy(read_buf, rc_buf.data());
                rc_buf.clear();
                return_rc_next = false;
                return strlen(read_buf);
            } else {
                if(!stream->eof()) return_rc_next = true;
            }
        }
        
        if(stream->eof()) return 0;

        int64_t header_length = 0;
        if(mode == FASTA){
            char c = 0;

            while(c != '\n'){
                // Make space if needed
                if(header_length >= header_buf_cap) {
                    header_buf_cap *= 2;
                    header_buf = (char*)realloc(header_buf, header_buf_cap);
                }

                // Read next character to buffer
                stream->get(&c); 
                header_buf[header_length++] = c;
            }
            header_buf[header_length-1] = '\0'; // Overwrite the newline with a null terminator

            int64_t buf_pos = 0;
            while(true){
                if(!stream->get(&c)) break; // Last read end
                else {
                    if(c == '\n') continue;
                    else if(c == '>') break;
                    else {
                        if(buf_pos + 1 >= read_buf_cap) { // +1: space for null terminator
                            read_buf_cap *= 2;
                            read_buf = (char*)realloc(read_buf, read_buf_cap);
                        }
                        read_buf[buf_pos++] = toupper(c);
                    }
                }
            }
            if(buf_pos == 0) throw std::runtime_error("Error: empty sequence in FASTA file.");
            read_buf[buf_pos] = '\0';
            if(reverse_complements){
                // Store the reverse complement for later
                for(int64_t i = 0; i < buf_pos+1; i++) // +1: also copy the null
                    rc_buf.push_back(read_buf[i]);
                reverse_complement_c_string(rc_buf.data(), buf_pos);
            }
            fix_alphabet(read_buf, buf_pos);
            return buf_pos;
        } else if(mode == FASTQ){
            char c = 0;
            while(c != '\n'){
                // Make space if needed
                if(header_length >= header_buf_cap) {
                    header_buf_cap *= 2;
                    header_buf = (char*)realloc(header_buf, header_buf_cap);
                }

                // Read next character to buffer
                stream->get(&c); 
                header_buf[header_length++] = c;
            }
            header_buf[header_length-1] = '\0'; // Overwrite the newline with a null terminator

            int64_t buf_pos = 0;
            while(true){
                stream->get(&c);
                if(c == '\n') break; // End of read
                if(buf_pos + 1 >= read_buf_cap) { // +1: space for null terminator
                    read_buf_cap *= 2;
                    read_buf = (char*)realloc(read_buf, read_buf_cap);
                }
                read_buf[buf_pos++] = toupper(c);
            }
            read_buf[buf_pos] = '\0';

            c = 0;
            while(c != '\n') stream->get(&c); // Skip '+'-line

            c = 0;
            while(c != '\n') stream->get(&c); // Skip quality line

            stream->get(&c); // Consume the '@' of the next read. If no more reads left, sets the eof flag.
            if(buf_pos == 0) throw std::runtime_error("Error: empty sequence in FASTQ file.");

            if(reverse_complements){
                // Store the reverse complement for later
                for(int64_t i = 0; i < buf_pos+1; i++) // +1: also copy the null
                    rc_buf.push_back(read_buf[i]);
                reverse_complement_c_string(rc_buf.data(), buf_pos);
                fix_alphabet(rc_buf.data(), buf_pos);
            }

            fix_alphabet(read_buf, buf_pos);
            return buf_pos;
        } else{
            throw std::runtime_error("Should not come to this else-branch");
        }
    }

    // Slow
    string get_next_read(){
        int64_t len = get_next_read_to_buffer();
        string read = (len > 0 ? string(read_buf) : "");
        return read;
    }

};

// Produces reads from multiple files like it was a single file
template<typename reader_t = SeqIO::Reader<>>
class Multi_File_Reader{

    public:

    char* read_buf; // Does not own this memory
    char* header_buf; // Does not own this memory

    vector<string> filenames;
    int64_t current_file_idx;
    std::unique_ptr<reader_t> reader;
    bool reverse_complements = false;

    Multi_File_Reader(const vector<string>& filenames) : filenames(filenames), current_file_idx(0){
        if(filenames.size() > 0){
            reader = make_unique<reader_t>(filenames[0]);
        }
    }

    int64_t get_next_read_to_buffer(){
        if(current_file_idx == filenames.size()) return 0; // All files processed

        int64_t len = reader->get_next_read_to_buffer();
        while(len == 0){ // End of file
            current_file_idx++;
            if(current_file_idx == filenames.size()) return 0; // All files processed
            reader = make_unique<reader_t>(filenames[current_file_idx]);
            if(reverse_complements) reader->enable_reverse_complements();
            len = reader->get_next_read_to_buffer();
        }

        this->read_buf = reader->read_buf; // Update pointer in case there was a realloc
        this->header_buf = reader->header_buf; // Update pointer in case there was a realloc
        return len;
    }

    void enable_reverse_complements(){
        reverse_complements = true;
        if(filenames.size() > 0) reader->enable_reverse_complements();
    }

    void rewind_to_start(){
        current_file_idx = 0;
        if(filenames.size() > 0){
            reader = make_unique<reader_t>(filenames[0]);
            if(reverse_complements) reader->enable_reverse_complements();
        }
    }
};


template<typename ofstream_t = Buffered_ofstream<std::ofstream>> // The underlying file stream.
class Writer{

    string fasta_header = ">\n";
    string fastq_header = "@\n";
    string newline = "\n";
    string plus = "+";

    public:

    ofstream_t out;
    int64_t mode;

    // Tries to figure out the format based on the file extension.
    Writer(string filename) : out(filename) {
        SeqIO::FileFormat fileformat = figure_out_file_format(filename);
        if(fileformat.format == FASTA) mode = FASTA;
        else if(fileformat.format == FASTQ) mode = FASTQ;
        else throw(runtime_error("Unknown file format: " + filename));
    }

    void write_sequence(const char* seq, int64_t len){
        if(mode == FASTA){
            // FASTA format
            out.write(fasta_header.c_str(), 2);
            out.write(seq, len);
            out.write(newline.c_str(), 1);
        } else{
            // FASTQ
            out.write(fastq_header.c_str(), 2);
            out.write(seq, len);
            out.write(newline.c_str(), 1);
            out.write(plus.c_str(), 1);
            out.write(newline.c_str(), 1);
            out.write(seq, len); // Use the read again for the quality values
            out.write(newline.c_str(), 1);
        }
    }

    // Flush the stream. The stream is also automatically flushed when the object is destroyed.
    void flush(){
        out.flush();
    }
};

int64_t count_sequences(const string& filename);

} // Namespace SeqIO

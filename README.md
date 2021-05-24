# Compiling

```
git submodule init
git submodule update

cd sdsl-lite
sh install.sh
cd ..

make toolkit
```

# Usage

For example, to compute a set of baits of length 120, allowing up to 40 mismatches, use the following:

```
./bin/greedy -L 120 -d 40 -s input.fna -o ./my_result_path
```

This writes three output files: 

| File name     | Explanation |
| ------------- | ------------- |
| ./my_result_path-baits.fna            | The baits in fasta format.  |
| ./my_result_path-cover-fractions.txt  | A text file containing one line for each bait. The i-th line is the fraction of characters covered by baits 0..i.  |
| ./my_result_path-cover-marks.txt      | A text file containing one line per input sequence. Each line contains an ascii string of zeroes and ones, with length equal to the length of the corresponding string. Character i of the j-th input sequence is covered by the bait set if the bit at position i of the j-th line is a one (indexing starts from zero). This is not and if-and-only-if condition! A position might be covered even if we have marked it with a zero if the FM-index seed-and-extend approximate search missed it. |

The full set of input parameter is listed below.

```
Computes a greedy bait cover.
Usage:
  ./bin/greedy [OPTION...]

  -L, --bait-len arg          Length of the baits. (default: 120)
  -d, --hamming-distance arg  Number of allowed mismatches in the baits.
                              (default: 40)
  -s, --sequences arg         Path to a fasta file of the input sequences.
                              (default: "")
      --fm-index-out arg      The algorithm is based on FM-index, which we
                              build at the start. Building the index can take a
                              lot of time and memory. Use this option to save
                              the FM-index to disk so that you can later run
                              the algorithm with different parameters re-using
                              the same FM-index. (optional). (default: "")
  -f, --fm-index arg          Path to a previously saved FM-index on disk
                              (--fm-index-out). This option loads the FM index
                              from disk instead of building it again. (default:
                              "")
  -o, --out arg               Filename prefix for the output files. (default:
                              "")
  -r, --randomize             Randomize the processing order the sequences in
                              the greedy algorithm.
  -t, --n-threads arg         Maximum number of parallel threads. The program
                              is not very well optimized for parallel
                              processing, so don't expect much of a speedup here.
                              (default: 1)
  -c, --cutoff arg            Stop the greedy algorithm after this fraction
                              of positions is covered. For example: 0.99.
                              (default: 1)
      --no-rev-comp           By default, a string also matches to its
                              reverse complement. Use this option to turn off
                              reverse complement matching.
  -g, --seed-len arg          The length of the seeds in the FM-index
                              seed-and-extend approximate string search subroutine. A
                              lower value will find more matches, but will be
                              slower. (default: 20)
  -h, --help                  Print instructions.

```

To compute various statistics on any bait set, run for example the following:

```
./bin/examine -d 40 -s input.fna -b baits.fna -o ./my_examine_output_path
```

This writes the following output files:

| File name     | Explanation |
| ------------- | ------------- |
| ./my_examine_output_path-cover_marks.txt            | Same as the cover marks file output by the greedy algorithm. |
| ./my_examine_output_path-cover-fractions.txt"       | Same as the cover fractions file output by the greedy algorithm. |
| ./my_examine_output_path-crossings.txt              | A text file with a number of lines equal to the lengths of the baits. The i-th file has a pair (x,y), were y is the fraction of baits that cross a branch in the de-bruijn graph of order x of the input sequences. |
| ./my_examine_output_path-gaps.txt | A text file consisting of one line of space-separated integers representing the lengths of the runs of zeroes in the cover marks file. The cover has no gaps the line is empty. |

The full set of input parameters is listed below

```
Computes various statistics on a given bait set.
Usage:
  ./bin/examine [OPTION...]

  -d, --hamming-distance arg  Number of allowed mismatches in the baits.
                              (default: 40)
  -s, --sequences arg         Path to the fasta file of the input sequences.
                              (default: "")
  -b, --baits arg             Path to the fasta file of the baits (default:
                              "")
  -o, --out-prefix arg        Filename prefix for the output files. (default:
                              "")
      --fm-index-out arg      The algorithm is based on FM-index, which we
                              build at the start. Building the index can take a
                              lot of time and memory. Use this option to save
                              the FM-index to disk so that you can later run
                              the algorithm with different parameters re-using
                              the same FM-index. (optional). (default: "")
  -f, --fm-index arg          Path to a previously saved FM-index on disk
                              (--fm-index-out). This option loads the FM index
                              from disk instead of building it again. (default:
                              "")
  -t, --n-threads arg         Maximum number of parallel threads. The program
                              is not very well optimized for parallel
                              processing, so don't expect much of a speedup here.
                              (default: 1)
  -g, --seed-len arg          The length of the seeds in the FM-index
                              seed-and-extend approximate string search subroutine. A
                              lower value will find more matches, but will be
                              slower. (default: 20)
  -v, --verbose               Print more debug output
  -h, --help                  Print instructions.

```

# Compiling the tests (optional)

```
cd googletest
mkdir build
cd build
cmake ..
make
cd ../..

make tests
```

The test executables will be written to ./bin
#!/bin/bash
set -e

# Activate Holy Build Box environment.
source /hbb_exe/activate

set -x

# Clone and enter source
git clone https://github.com/jnalanko/syotti
cd syotti

# Compile

git submodule init
git submodule update

cd sdsl-lite
sh install.sh
cd ..

make syotti

# Copy result to host
cp bin/syotti /io/

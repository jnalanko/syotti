# Compiling

```
git submodule init
git submodule update

cd sdsl-lite
sh install.sh
cd ..

cd googletest
mkdir build
cd build
cmake ..
make

cd ../..
make toolkit
make tests # Optional
```


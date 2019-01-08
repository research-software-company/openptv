# Build liboptv
set -ev
cd liboptv
mkdir _build && cd _build
# ./configure --prefix=$PREFIX --with-charset=utf
cmake ../
make
make verify
make
cp src/liboptv.dylib ../../py_bind/
 
# Build liboptv Python wrappers
cd ../../py_bind
python setup.py build_ext -I/usr/local/include -L.
python setup.py install --prefix=$PREFIX
export PATH=$PATH:/usr/local/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/usr/local/lib
echo $LD_LIBRARY_PATH
echo $PATH
python -c "import sys; print(sys.path)"
cd test
nosetests -v
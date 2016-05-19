# Speed comparison between naive implementation and OpenBLAS
This is a simple speedtest with matrix matrix multiplication and matrix
vector multiplication with OpenBLAS and a naive implementation. You can change
the value `rounds` to any desired number to alter the time used for the calculation.

### Usage
If you installed [OpenBLAS](http://www.openblas.net/) with the standard directory, you
will need to use:
* `export LD_LIBRARY_PATH=/opt/OpenBLAS/lib/` (or to your directory with OpenBLAS)
* `g++ -o benchmark main.cpp -I /opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas` for compiling

### Notes
The time is used from the CPU clock and may differ from realtime but since most of the
time is running on the CPU, I prefer this way. If you have a different opinion, feel
free to use `time()`.
Please note the different representation for the matrices and the vector for OpenBLAS.
libtnp
======

arbitrary-order automatic differentiation library

libtnp offers an implementation for the AD method described in [1]

Build & Installation
--------------------

Requirements:

        * C99 compiler
        * C++11 compiler
        * Boost http://www.boost.org/
        * cmake
        * make
        
If you are running on linux or a decent unix, you should be able to install all required dependencies using your package manager, e.g. 

```sh
apt-get install build-essential cmake libboost-all-dev
```
Note, that both gcc and clang should work well.

In any other case (i.e. MS Windows) you will have to download and install all required tools yourself.

Since cmake favours out-of-tree builds, the next step is to create a dedicated build directory and run the build process in it. 

```sh
mkdir build
cd build
cmake .. # add your specific cmake flags here
make && make install
```

In case of non-standard installations, specific installation targets etc. you have to set the appropriate cmake flags.

Footnotes
---------

[1] Operational Semantics for a Modular Equation Language, C.Höger, 2013, Proceedings of the 4th Analytic Virtual Integration of Cyber-Physical Systems Workshop, December 3, Vancouver, Canada


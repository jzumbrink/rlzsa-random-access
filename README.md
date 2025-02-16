rlzsa-random-access
===============
A modification of the r-index and the r-index-rlzsa.

R-Index Author: Nicola Prezza (nicola.prezza@gmail.com)
Joint work with Travis Gagie and Gonzalo Navarro

rlzsa has been added according to [1] and [2]


cite as:

Gagie T, Navarro G, Prezza N. Optimal-time text indexing in BWT-runs bounded space. In Proceedings of the Twenty-Ninth Annual ACM-SIAM Symposium on Discrete Algorithms, SODA 2018, New Orleans, NA, USA, January 7-10 2017.

Travis Gagie, Gonzalo Navarro, and Nicola Prezza. Fully Functional Suffix Trees and Optimal Text Searching in BWT-Runs Bounded Space. J. ACM 67, 1, Article 2 (April 2020)

### Brief description

RLZ-compressed sufffix array. Only offers support for random access.
### Download

To clone the repository, run:

> git clone https://github.com/jzumbrink/rlzsa-random-access

### Compile

We use cmake to generate the Makefile. Create a build folder in the main r-index folder:

> mkdir build

run cmake:

> cd build; cmake ..

and compile:

> make

### Run

After compiling, run 

>  rlzsa-build input

This command will create the rlz compressed suffix array of the text file "input" and will store it as "input.rlzsa". Use option -o to specify a different basename for the index file. 

Run

> rlzsa-ra index.rlzsa indices

to perform random access with indices specified in indices.

### Funding

Nicola Prezza has been supported by the project Italian MIUR-SIR CMACBioSeq ("Combinatorial methods for analysis and compression of biological sequences") grant n.~RBSI146R5L, PI: Giovanna Rosone. Link: http://pages.di.unipi.it/rosone/CMACBioSeq.html

### References
[1] Simon J. Puglisi and Bella Zhukova. Smaller RLZ-Compressed Suffix Arrays. In 31st Data Compression Conference (DCC), 2021 ([pdf](https://ieeexplore.ieee.org/document/9418726))

[2] Bella Zhukova, New space-time trade-offs for pattern matching with compressed indexes, PhD Thesis 2024 ([pdf](https://helda.helsinki.fi/items/a672a30c-0611-408c-abeb-36aa069df2a1))
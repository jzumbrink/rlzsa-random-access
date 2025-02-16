// Copyright (c) 2017, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>

#include "internal/rlzsa.hpp"
#include "utils.hpp"

using namespace rlz;
using namespace std;
using namespace sdsl;

string out_basename=string();
string input_file=string();
string integer_file=string();
int sa_rate = 512;
bool sais=true;
ulint T = 0;//Build fast index with SA rate = T

void help(){
	cout << "rlzsa-build: builds the rlzsa-index. Extension .rlzsa is automatically added to output index file" << endl << endl;
	cout << "Usage: rlzsa-build [options] <input_file_name>" << endl;
	cout << "   <input_file_name>    input text file." << endl;
	cout << "   -o <basename>        use 'basename' as prefix for all index files. Default: basename is the specified input_file_name"<<endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if(s.compare("-o")==0){

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -o option." << endl;
			help();
		}

		out_basename = string(argv[ptr]);
		ptr++;

	}else if(s.compare("-divsufsort")==0){

		sais = false;

	} else {
		cout << "Error: unrecognized '" << s << "' option." << endl;
		help();
	}

}

int main(int argc, char** argv){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

    auto t1 = high_resolution_clock::now();

	//parse options

    out_basename=string();
    input_file=string();
	int ptr = 1;

	if(argc<2) help();

	while(ptr<argc-1)
		parse_args(argv, argc, ptr);

	input_file = string(argv[ptr++]);

	if(out_basename.compare("")==0)
		out_basename = string(input_file);

	string idx_file = out_basename;
	idx_file.append(".rlzsa");


	cout << "Building rlzsa-index of input file " << input_file << endl;
	cout << "Index will be saved to " << idx_file << endl;

	string input;

	{

		std::ifstream fs(input_file);
		std::stringstream buffer;
		buffer << fs.rdbuf();

		input = buffer.str();

	}

	string path = string(out_basename).append(".rlzsa");
	std::ofstream out(path);

	auto idx = rlzsa<>(input);

	idx.serialize(out);


	auto t2 = high_resolution_clock::now();
	ulint total = duration_cast<duration<double, std::ratio<1>>>(t2 - t1).count();
	cout << "Build time : " << get_time(total) << endl;


	out.close();

	idx.log_data_structure_sizes();

	ulint time_ms = duration_cast<chrono::milliseconds>(t2 - t1).count();
	cout << "RESULT"
		<< " algo=rlzsa_build"
		<< " time_ms=" << time_ms
		<< " n=" << idx.n
		<< " idx_size=" << idx.size_in_bytes()
		<< endl;
}
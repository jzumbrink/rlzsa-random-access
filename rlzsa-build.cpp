// Copyright (c) 2017, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>

#include "internal/r_index.hpp"
#include "utils.hpp"
#include "internal/r_index.hpp"

using namespace ri_rlzsa;
using namespace std;

string out_basename=string();
string input_file=string();
string integer_file=string();
int sa_rate = 512;
bool sais=true;
ulint T = 0;//Build fast index with SA rate = T
bool fast = false;//build fast index
bool hyb = false; //use hybrid bitvectors instead of sd_vectors?
ulint reference_size = 0;

void help(){
	cout << "rlzsa-build: builds the rlzsa-index. Extension .rlzsa is automatically added to output index file" << endl << endl;
	cout << "Usage: rlzsa-build [options] <input_file_name> <reference size> <integer file>" << endl;
	cout << "   <input_file_name>    input text file." << endl;
	cout << "   <reference size>     the size of the reference string used for rlz (typically around 5*r)." << endl;
	cout << "   <integer file>       path to file which contains an integer value at each line." << endl;
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

	if(argc<4) help();

	while(ptr<argc-3)
		parse_args(argv, argc, ptr);

	input_file = string(argv[ptr++]);
	reference_size = std::stoll(argv[ptr++]);
	integer_file = string(argv[ptr++]);

	// load indices
	cout << "Loading indices" << endl;
	std::ifstream arbitrary_indices_in(integer_file);
	string last_index = string();
	vector<int64_t> arbitrary_indices;
	while(arbitrary_indices_in >> last_index) {
		arbitrary_indices.push_back(std::stoll(last_index));
	}

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

	//save flag storing whether index is fast or small
	out.write((char*)&fast,sizeof(fast));


	auto idx = r_index<>();

	int64_t n = input.size();
	auto sa = std::make_unique<int64_t[]>(n);
	libsais64((uint8_t const*) input.data(), sa.get(), n, 0, nullptr);

	vector<int64_t> SA_d;
	SA_d.reserve(n);
	if (n > 0) SA_d.push_back(sa[0]);
	for (int i = 1; i < n; ++i) {
		SA_d.push_back(sa[i] + n + 1 - sa[i - 1]);
	}
	std::cout << "SA_d size=" << SA_d.size() << std::endl;
	std::cout << "using rlz reference_size" << reference_size << std::endl;

	//std::function<uint64_t(uint64_t)> SA_function = [std::move(sa)](uint64_t i) { return sa[i]; };
	auto SA_function = [&](uint64_t i){return sa[i];};

	idx.build_rlzsa<int64_t>(SA_d, SA_function, reference_size);

	idx.serialize(out);


	auto t2 = high_resolution_clock::now();
	ulint total = duration_cast<duration<double, std::ratio<1>>>(t2 - t1).count();
	cout << "Build time : " << get_time(total) << endl;


	out.close();

	auto ra_start = high_resolution_clock::now();
	// measure random access time
	int64_t a = 0;
	for (int64_t i = 0; i < arbitrary_indices.size(); i++) {
		a += idx.access_SA(arbitrary_indices[i]);
	}
	auto ra_end = high_resolution_clock::now();
	cout << a << endl;

	cout << "Managed " << arbitrary_indices.size() << " random accesses in " << duration_cast<chrono::milliseconds>(ra_end - ra_start).count() << " ms" << endl;

}
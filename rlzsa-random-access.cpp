// Copyright (c) 2017, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>

#include "internal/rlzsa.hpp"
#include "utils.hpp"

using namespace rlz;
using namespace std;
using namespace sdsl;

int sa_rate = 512;
ulint T = 0;//Build fast index with SA rate = T

void help(){
	cout << "rlzsa-ra: does random access on the rlzsa" << endl << endl;
	cout << "Usage: rlzsa-build <rlzsa file> <integer file>" << endl;
    cout << "   <rlzsa file>         path to file which contains the RLZ-compressed suffix array." << endl;
	cout << "   <integer file>       path to file which contains an integer value at each line." << endl;
	exit(0);
}

int main(int argc, char** argv){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

    auto t1 = high_resolution_clock::now();
    if(argc<3) help();

    string rlzsa_file=string(argv[1]);
    string integer_file=string(argv[2]);

	// load indices
	cout << "Loading indices" << endl;
	std::ifstream arbitrary_indices_in(integer_file);
	string last_index = string();
	vector<int64_t> arbitrary_indices;
	while(arbitrary_indices_in >> last_index) {
		arbitrary_indices.push_back(std::stoll(last_index));
	}

	auto idx = rlzsa<>();
    idx.load_from_file(rlzsa_file);


	auto t2 = high_resolution_clock::now();
	ulint total = duration_cast<duration<double, std::ratio<1>>>(t2 - t1).count();
	cout << "Load time : " << get_time(total) << endl;

	auto ra_start = high_resolution_clock::now();
	// measure random access time
	int64_t a = 0;
	for (int64_t i = 0; i < arbitrary_indices.size(); i++) {
		a += idx.access_SA(arbitrary_indices[i]);
	}
	auto ra_end = high_resolution_clock::now();

    ulint time_ms = duration_cast<chrono::milliseconds>(ra_end - ra_start).count();
    cout << "Managed " << arbitrary_indices.size() << " random accesses in " << time_ms << " ms" << endl;
    idx.log_data_structure_sizes();
    cout << "RESULT"
        << " algo=rlzsa_ra"
        << " tp_ms=" << (double) arbitrary_indices.size() / (double) time_ms
        << " time_ms=" << time_ms
        << " idx_size=" << idx.size_in_bytes()
        << " a=" << a
        << " bits_per_symbol=" << (double) idx.size_in_bytes() / idx.n * 8
        << endl;
}
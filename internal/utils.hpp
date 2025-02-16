// Copyright (c) 2017, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#pragma once

#include <sstream>

#include <malloc_count.h>

using namespace std;

namespace rlz{

using ulint = uint64_t;

string get_time(uint64_t time){

	stringstream ss;

	if(time>=3600){

		uint64_t h = time/3600;
		uint64_t m = (time%3600)/60;
		uint64_t s = (time%3600)%60;

		ss  << time << " seconds. ("<< h << "h " << m << "m " << s << "s" << ")";

	}else if (time>=60){

		uint64_t m = time/60;
		uint64_t s = time%60;

		ss << time << " seconds. ("<< m << "m " << s << "s" << ")";

	}else{

		ss << time << " seconds.";

	}

	return ss.str();

}

inline uint8_t char_to_uchar(char c) {
    return *reinterpret_cast<uint8_t*>(&c);
}

uint8_t bitsize(uint64_t x){

	if(x==0) return 1;
	return 64 - __builtin_clzll(x);

}

//parse pizza&chilli patterns header:
void header_error(){
	cout << "Error: malformed header in patterns file" << endl;
	cout << "Take a look here for more info on the file format: http://pizzachili.dcc.uchile.cl/experiments.html" << endl;
	exit(0);
}

ulint get_number_of_patterns(string header){

	ulint start_pos = header.find("number=");
	if (start_pos == std::string::npos or start_pos+7>=header.size())
		header_error();

	start_pos += 7;

	ulint end_pos = header.substr(start_pos).find(" ");
	if (end_pos == std::string::npos)
		header_error();

	ulint n = std::atoi(header.substr(start_pos).substr(0,end_pos).c_str());

	return n;

}

ulint get_patterns_length(string header){

	ulint start_pos = header.find("length=");
	if (start_pos == std::string::npos or start_pos+7>=header.size())
		header_error();

	start_pos += 7;

	ulint end_pos = header.substr(start_pos).find(" ");
	if (end_pos == std::string::npos)
		header_error();

	ulint n = std::atoi(header.substr(start_pos).substr(0,end_pos).c_str());

	return n;

}

template <typename T>
bool contains(std::vector<T> vec, T val) {
    return std::find(vec.begin(),vec.end(),val) != vec.end();
}

std::chrono::steady_clock::time_point now() {
    return std::chrono::steady_clock::now();
}

std::string format_time(uint64_t ns) {
    std::string time_str;

    if (ns > 10000000000) {
        time_str = std::to_string(ns/1000000000) + " s";
    } else if (ns > 10000000) {
        time_str = std::to_string(ns/1000000) + " ms";
    } else if (ns > 10000) {
        time_str = std::to_string(ns/1000) + " us";
    } else {
        time_str = std::to_string(ns) + " ns";
    }

    return time_str;
}

uint64_t time_diff_ns(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
}

std::chrono::steady_clock::time_point log_runtime(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2) {
    std::cout << ", in ~ " << format_time(time_diff_ns(t1,t2)) << std::endl;
    return std::chrono::steady_clock::now();
}

std::chrono::steady_clock::time_point log_runtime(std::chrono::steady_clock::time_point t) {
    return log_runtime(t,std::chrono::steady_clock::now());
}

inline char uchar_to_char(uint8_t c) {
    return *reinterpret_cast<char*>(&c);
}

uint64_t time_diff_min(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2) {
    return std::chrono::duration_cast<std::chrono::minutes>(t2-t1).count();
}

std::string format_size(uint64_t B) {
    std::string size_str;

    if (B > 10000000000) {
        size_str = std::to_string(B/1000000000) + " GB";
    } else if (B > 10000000) {
        size_str = std::to_string(B/1000000) + " MB";
    } else if (B > 10000) {
        size_str = std::to_string(B/1000) + " KB";
    } else {
        size_str = std::to_string(B) + " B";
    }
    
    return size_str;
}

}
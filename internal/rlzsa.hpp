// Copyright (c) 2017, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 * rlzsa.hpp
 *
 *  Created on: Apr 13, 2017
 *      Author: nico
 *
 * Small version of the r-index: O(r) words of space, O(log(n/r)) locate time per occurrence
 *
 */

#pragma once

#include "definitions.hpp"
#include "rle_string.hpp"
#include "utils.hpp"
#include <sdsl/csa_wt.hpp>
#include <sdsl/csa_bitcompressed.hpp>
#include <sdsl/suffix_array_algorithm.hpp>
#include <ips4o.hpp>
#include <ankerl/unordered_dense.h>
#include <libsais.h>
#include <libsais64.h>

using namespace sdsl;

namespace rlz{

template<class rle_string_t = rle_string_sd>
class rlzsa{

public:

	using triple = std::tuple<range_t, ulint, ulint>;

	rlzsa(){}

	static constexpr ulint a = 64;

	/*
	 * Build index
	 */
	rlzsa(string &input, bool log = true) {
		auto time = now();

		n = input.size() + 1;
		if (log) cout << "n = " << n << endl;

		if (log) cout << "building BWT and computing SA samples (SE-SAIS)";

		auto bwt_and_sa = sufsort(input);

		string& bwt_s = get<0>(bwt_and_sa);
		int_vector_buffer<>& SA = *get<1>(bwt_and_sa);
		cache_config cc_sa = get<2>(bwt_and_sa);

		auto bwt = rle_string_sd(bwt_s);
		r = bwt.number_of_runs();

		int log_r = bitsize(uint64_t(r));
		int log_n = bitsize(uint64_t(n));

		if (log) {
			cout << "Number of BWT equal-letter runs (used for determination of reference size): r = " << r << endl;
		}

		if (n + 1 <= std::numeric_limits<int32_t>::max()) {
			build_rlzsa_internal<int32_t>([&](ulint i){return SA[i];}, log);
		} else {
			build_rlzsa_internal<int64_t>([&](ulint i){return SA[i];}, log);
		}

	    sdsl::remove(cache_file_name(conf::KEY_SA, cc_sa));
	}

	public:
	template <typename sad_t>
	void build_rlzsa_internal(std::function<ulint(ulint)> SA, bool log = false) {
		auto time = now();
		if (log) cout << "building SA^d" << flush;

		vector<sad_t> SA_d;
		SA_d.reserve(n + 1);
		SA_d.resize(n);
		SA_d[0] = SA(0);

		for (ulint i=1; i<n; i++) {
			SA_d[i] = (SA(i)+n)-SA(i-1);
		}
		
		if (log) time = log_runtime(time);
		ulint size_R_target = std::min(std::max<ulint>(1,n/3),5*r);
		build_rlzsa<sad_t>(SA_d, SA, size_R_target, log);
	}

	public:
	template <typename sad_t = int32_t>
	void build_rlzsa(vector<sad_t>& SA_d, std::function<ulint(ulint)> SA, ulint reference_size, bool log = true) {
		/*std::cout << "SA:";
		for (int i = 0; i < 100; i++) {
			std::cout << SA(i) << ", ";
		}
		std::cout << std::endl;
		std::cout << "SA^d:";
		for (int i = 0; i < 100; i++) {
			std::cout << SA_d[i] << ", ";
		}
		std::cout << std::endl;*/

		n = SA_d.size();
		auto time = now();
		if (log) cout << "sorting SA^d" << flush;

		vector<ulint> alphabet; // alphabet of SA^d
		alphabet.resize(n);

		for (ulint i=0; i<n; i++) {
			alphabet[i] = SA_d[i];
		}

		ips4o::sort(alphabet.begin(),alphabet.end()); // sort the values in SA^d
	
		if (log) time = log_runtime(time);
		if (log) cout << "building SA^d alphabet map" << flush;

		auto it = std::unique(alphabet.begin(),alphabet.end()); // remove duplicates in SA^d
		alphabet.resize(std::distance(alphabet.begin(),it));
		ulint sigma_sad = alphabet.size();
		ankerl::unordered_dense::map<ulint,ulint> alphabet_map; // map to the effective alphabet of SA^d
		int_vector<> alphabet_map_inv; // map from the effective alphabet of SA^d
		alphabet_map_inv.width(ceil(log2((2*n))/8.0)*8);
		alphabet_map_inv.resize(sigma_sad+1);
		alphabet_map_inv[0] = 0;

		{
			ulint j = 1;

			// build alphabet_map and alphabet_map_inv
			for (ulint i=0; i<sigma_sad; i++) {
				alphabet_map.insert({alphabet[i],j});
				alphabet_map_inv[j] = alphabet[i];
				j++;
			}
		}

		if (log) time = log_runtime(time);
		if (log) cout << "mapping SA^d to its effective alphabet" << flush;

		// map SA^d to its effective alphabet
		for (ulint i=0; i<n; i++) {
			SA_d[i] = alphabet_map[SA_d[i]];
		}

		if (log) time = log_runtime(time);
		if (log) cout << "building suffix array of SA^d" << flush;

		sad_t fs = 6 * (sigma_sad + 1);
		vector<sad_t> SA_sad;
		SA_sad.resize(n + fs);

		if constexpr (std::is_same_v<sad_t, int32_t>) {
			libsais_int(SA_d.data(), SA_sad.data(), n, sigma_sad + 1, fs);
		} else {
			libsais64_long(SA_d.data(), SA_sad.data(), n, sigma_sad + 1, fs);
		}

		SA_sad.resize(n);
		if (log) time = log_runtime(time);
		if (log) cout << "building inverse suffix array of SA^d" << flush;

		int_vector<> ISA_sad;
		ISA_sad.width(ceil(log2(n)/8.0)*8);
		ISA_sad.resize(n);

		for (ulint i=0; i<n; i++) {
			ISA_sad[SA_sad[i]] = i;
		}

		if (log) time = log_runtime(time);
		if (log) cout << "counting frequencies of k-mers in SA^d" << flush;

		ulint size_R_target = reference_size;
		
		std::cout << "size_R_target=" << size_R_target << std::endl;

		ulint segment_size = std::min<ulint>(4096,size_R_target);
		ulint num_segments = n/segment_size;
		vector<ulint> selected_segments;
		ulint k = std::min<ulint>(8,size_R_target);
		ulint num_groups;
		
		int_vector<> group;
		int_vector<> freq;
		
		group.width(ceil(log2(n)/8.0)*8);
		freq.width(ceil(log2(n)/8.0)*8);

		group.resize(n);
		freq.resize(n);

		{
			ulint i = 1;
			ulint g = 0;

			while (i < n) {
				ulint j = 1;
				bool equal = true;

				while (equal && i+j < n) {
					for (ulint l=0; l<k; l++) {
						if (
							SA_sad[i+j]+l == n ||
							SA_sad[i+j-1]+l == n ||
							SA_d[SA_sad[i+j]+l] != SA_d[SA_sad[i+j-1]+l]
						) {
							equal = false;
							break;
						}
					}
					
					if (equal) {
						j++;
					}
				}

				for (ulint p=0; p<j; p++) {
					group[i+p] = g;
					freq[i+p] = j;
				}

				g++;
				i += j;
			}

			num_groups = g;
		}

		if (log) time = log_runtime(time);
		if (log) cout << "counting frequencies of k-mers per segment of SA^d" << flush;

		{
			ulint num_segments_to_select = size_R_target/segment_size;
			std::vector<double> segment_freqs;
			segment_freqs.resize(num_segments+1,0);

			for (ulint s=0; s<num_segments; s++) {
				sdsl::bit_vector is_group_processed;
				is_group_processed.resize(num_groups+1);

				for (ulint i=s*segment_size; i<(s+1)*segment_size-k; i++) {
					ulint p = ISA_sad[i];

					if (is_group_processed[group[p]] == 0) {
						segment_freqs[s] += std::sqrt((double)freq[p]);
						is_group_processed[group[p]] = 1;
					}
				}
			}

			if (log) time = log_runtime(time);
			if (log) cout << "building R" << flush;

			selected_segments.resize(num_segments_to_select);
			sdsl::bit_vector is_segment_selected;
			is_segment_selected.resize(num_segments+1);
			sdsl::bit_vector is_group_processed;
			is_group_processed.resize(num_groups+1);

			for (ulint s=0; s<num_segments_to_select; s++) {
				ulint best_segment = 0;
				
				for (ulint i=1; i<num_segments; i++) {
					if (is_segment_selected[i] == 0 && segment_freqs[i] > segment_freqs[best_segment]) {
						best_segment = i;
					}
				}
				
				selected_segments[s] = best_segment;
				is_segment_selected[best_segment] = 1;
				segment_freqs[best_segment] = 0;

				for (ulint i=best_segment*segment_size; i<(best_segment+1)*segment_size-k; i++) {
					ulint p = ISA_sad[i];

					if (p < n && is_group_processed[group[p]] == 0) {
						ulint g = group[p];
						double freq_sqrt = std::sqrt((double)freq[p]);

						while (p > 0 && group[p-1] == g) {
							p--;
						}

						sdsl::bit_vector is_segment_processed;
						is_segment_processed.resize(num_segments+1);
						is_segment_processed[best_segment] = 1;

						do {
							ulint seg = SA_sad[p]/segment_size;

							if (
								SA_sad[p] < (seg+1)*segment_size-k &&
								is_segment_processed[seg] == 0 &&
								is_segment_selected[seg] == 0
							) {
								segment_freqs[seg] -= freq_sqrt;
								is_segment_processed[seg] = 1;
							}
							
							p++;
						} while (p < n && group[p] == g);

						is_group_processed[g] = 1;
					}
				}
			}

			ips4o::sort(selected_segments.begin(),selected_segments.end());
		}

		{
			R.width(ceil(log2((2*n))/8.0)*8);
			R.resize(segment_size*selected_segments.size());
			ulint j = 0;

			for (ulint i=0; i<selected_segments.size(); i++) {
				ulint segment_start = selected_segments[i]*segment_size;

				for (ulint l=0; l<segment_size; l++) {
					R[j] = alphabet_map_inv[SA_d[segment_start+l]];
					j++;
				}
			}
		}

		selected_segments.clear();
		selected_segments.shrink_to_fit();

		selected_segments.clear();
		selected_segments.shrink_to_fit();

		group.resize(0);

		freq.resize(0);

		SA_sad.clear();
		SA_sad.shrink_to_fit();

		ISA_sad.resize(0);

		if (log) time = log_runtime(time);
		if (log) cout << "building FM-Index for R" << flush;

		int_vector<> R_rev;
		R_rev.width(ceil(log2((2*n))/8.0)*8);
		R_rev.resize(R.size());

		for (ulint i=0; i<R.size(); i++) {
			R_rev[i] = R[R.size()-i-1];
		}

		sdsl::csa_bitcompressed<sdsl::int_alphabet<>> FM_rrev;
		sdsl::construct_im(FM_rrev,R_rev,0);
		R_rev.resize(0);

		if (log) time = log_runtime(time);
		if (log) cout << "building rlzsa" << flush;

		S.width(ceil(log2((n))/8.0)*8);
		PS.width(ceil(log2((n))/8.0)*8);

		S.resize(1);
		PS.resize(1);
		PL.resize(1);
		S[0] = SA(0);
		PS[0] = 0;
		PL[0] = 0;
		ulint z_l = 1;

		{
			ulint i = 1;
			ulint z = 1;
			int_vector<> SA_d_val;
			SA_d_val.resize(1);
			
			// rlz-encode SA^d
			while (i < n) {
				ulint b = 0;
				ulint e = FM_rrev.size()-1;
				ulint b_last,e_last,pos_sa;
				ulint l = 0;

				// find the longest prefix SA^d[i,i+d) of SA^d[i,n) that occurs in R
				while (true) {
					b_last = b;
					e_last = e;
					SA_d_val[0] = alphabet_map_inv[SA_d[i+l]];

					if (sdsl::backward_search(FM_rrev, b, e, SA_d_val.begin(), SA_d_val.end(), b, e) == 0) {
						break;
					}

					l++;

					if (l == UINT16_MAX || i+l == n) {
						b_last = b;
						e_last = e;
						break;
					}
				}

				z++;
				PL.resize(z);
				PL[z-1] = l;
				S.resize(z);

				if (z % a == 1) {
					PS.resize(PS.size()+1);
					PS[PS.size()-1] = i;
				}

				if (l == 0) {
					// literal phrase
					S[z-1] = SA(i);
					i++;
					z_l++;
				} else {
					// copy-phrase
					S[z-1] = (R.size()-FM_rrev[b_last])-l;
					i += l;

					// add a literal phrase after each copy-phrase
					if (i < n) {
						z++;
						PL.resize(z);
						PL[z-1] = 0;
						S.resize(z);
						S[z-1] = SA(i);

						if (z % a == 1) {
							PS.resize(PS.size()+1);
							PS[PS.size()-1] = i;
						}
						
						i++;
						z_l++;
					}
				}
			}

			PL.resize(z+1);
			PL[z] = 1;
		}

		if (log) time = log_runtime(time);

		if (log) {
			std::cout << std::endl;
			std::cout << "z: " << PL.size() << std::endl;
			std::cout << "z_l/z: " << z_l/(double)PL.size() << std::endl;
			std::cout << "z/r: " << PL.size()/(double)std::max<ulint>(r,1) << std::endl;
			std::cout << "peak memory usage: " << format_size(malloc_count_peak()) << std::endl;
		}

		std::cout << "Build RLZSA finished" << std::endl;
		/*
		std::cout << "SA: ";
		for (ulint i=0; i<n; i++) {
			std::cout << SA(i) << ", ";
		}

		std::cout << std::endl << "SA^d: ";
		std::cout << SA(0) << ", ";
		for (ulint i=1; i<n; i++) {
			std::cout << sad_t{alphabet_map_inv[SA_d[i]]}-sad_t{n} << ", ";
		}

		std::cout << std::endl << "R: ";
		for (ulint i=0; i<R.size(); i++) {
			std::cout << sad_t{R[i]}-sad_t{n} << ", ";
		}
		
		std::cout << std::endl << "PL: ";
		for (ulint i=0; i<PL.size(); i++) {
			std::cout << PL[i] << ", ";
		}

		std::cout << std::endl << "S: ";
		for (ulint i=0; i<S.size(); i++) {
			std::cout << S[i] << ", ";
		}

		std::cout << std::endl << "PS: ";
		for (ulint i=0; i<PS.size(); i++) {
			std::cout << PS[i] << ", ";
		}
		std::cout << std::endl;
		*/

	}

	ulint access_SA(ulint i) {
		ulint x_ps; // index in PS of the last phrase starting before or at position i
		{
			ulint b = 0;
			ulint e = PS.size()-1;
			ulint m;

			while (b != e) {
				m = b+(e-b)/2+1;

				if (PS[m] <= i) {
					b = m;
				} else {
					e = m-1;
				}
			}

			x_ps = b;
		}

		ulint s_cp = PS[x_ps]; // starting position of the current phrase
		ulint x_p = a * x_ps;
		
		// find the phrase containing i
		while (s_cp + std::max<ulint>(1,PL[x_p]) <= i) {
			s_cp += std::max<ulint>(1,PL[x_p]);
			x_p++;
		}

		if (PL[x_p] == 0) {
			// the x_p-th phrase is a literal phrase
			return S[x_p];
		} else {
			// the x_p-th phrase is a copy phrase, hence the x_p-1-th phrase is a
			// literal phrase, because there is a literal phrase before each copy-phrase;
			// so store its SA-value in s and decode the x_p-th phrase up to position i

			ulint j = s_cp;
			ulint s = S[x_p-1];
			ulint p_r = S[x_p]; // current position in R

			while (j <= i) {
				s += R[p_r];
				s -= n;
				p_r++;
				j++;
			}

			return s;
		}
	}

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	ulint serialize(std::ostream& out){
		ulint w_bytes = 0;

		w_bytes += R.serialize(out);
		w_bytes += PS.serialize(out);
		w_bytes += PL.serialize(out);
		w_bytes += S.serialize(out);

		out.write((char*)&n, sizeof(ulint));
		w_bytes += sizeof(ulint);

		std::cout << "Wrote " << w_bytes << " to disk." << std::endl;

		return w_bytes;
	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {
		R.load(in);
		PS.load(in);
		PL.load(in);
		S.load(in);

		in.read((char*)&n, sizeof(ulint));
	}

	/*
	 * save the structure to the path specified.
	 * \param path_prefix prefix of the index files. suffix ".ri" will be automatically added
	 */
	void save_to_file(string path_prefix){
		string path = string(path_prefix).append(".rlzsa");

		std::ofstream out(path);
		serialize(out);
		out.close();
	}

	/*
	 * load the structure from the path specified.
	 * \param path: full file name
	 */
	void load_from_file(string path){
		std::ifstream in(path);
		load(in);
		in.close();
	}
	
    uint64_t size_in_bytes() {
        uint64_t size = 0;

		size += sdsl::size_in_bytes(R);
		size += sdsl::size_in_bytes(PS);
		size += sdsl::size_in_bytes(PL);
		size += sdsl::size_in_bytes(S);            

        return size;
    }

    void log_data_structure_sizes() {
        std::cout << "index size: " << format_size(size_in_bytes()) << std::endl;

        std::cout << "R: " << format_size(sdsl::size_in_bytes(R)) << std::endl;
        std::cout << "PS: " << format_size(sdsl::size_in_bytes(PS)) << std::endl;
        std::cout << "PL: " << format_size(sdsl::size_in_bytes(PL)) << std::endl;
        std::cout << "S: " << format_size(sdsl::size_in_bytes(S)) << std::endl;
    }

	/*
	 * returns a triple containing BWT of input string
	 * (uses 0x1 character as terminator), text positions corresponding
	 * to first letters in BWT runs (plus their ranks from 0 to R-1), and text positions corresponding
	 * to last letters in BWT runs (in BWT order)
	 */
	tuple<string, int_vector_buffer<>*, cache_config> sufsort(string &s){

		string bwt_s;

	    cache_config cc;

	    int_vector<8> text(s.size());
	    assert(text.size()==s.size());

	    for(ulint i=0;i<s.size();++i)
	    	text[i] = (uchar)s[i];

	    assert(text.size()==s.size());

	    append_zero_symbol(text);

	    store_to_cache(text, conf::KEY_TEXT, cc);

	    construct_config::byte_algo_sa = SE_SAIS;
	    construct_sa<8>(cc);

	    //now build BWT from SA
	    int_vector_buffer<>* sa = new int_vector_buffer<>(cache_file_name(conf::KEY_SA, cc));

		{

	        for (ulint i=0; i<sa->size(); i++){
	            auto x = (*sa)[i];

	            assert(x<=text.size());

	            if ( x > 0 )
	            	bwt_s.push_back((uchar)text[x-1]);
	            else
	            	bwt_s.push_back(TERMINATOR);
	        }

	    }

	    sdsl::remove(cache_file_name(conf::KEY_TEXT, cc));
	    //sdsl::remove(cache_file_name(conf::KEY_SA, cc));

	    return tuple<string, int_vector_buffer<>*, cache_config>(bwt_s, sa, cc);

	}

	static const uchar TERMINATOR = 1;

	ulint n = 0;
	ulint r = 0;//number of BWT runs

	sdsl::int_vector<> R; // rlzsa reference
	sdsl::int_vector<> PS; // rlzsa phrase starting positions (every a-th phrase is sampled)
	sdsl::int_vector<16> PL; // rlzsa phrase lengths
	sdsl::int_vector<> S; // starting positions of the rlzsa phrases in R
};

}
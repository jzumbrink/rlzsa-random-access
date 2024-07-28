#include <gtest/gtest.h>
#include "internal/r_index.hpp"

using namespace ri_rlzsa;

TEST(test_move_r,fuzzy_test) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::lognormal_distribution<double> avg_input_rep_length_distrib(4.0,2.0);
    std::uniform_real_distribution<double> prob_distrib(0.0,1.0);
    std::uniform_int_distribution<uint8_t> alphabet_size_distrib(1,26);
    std::uniform_int_distribution<uint8_t> uchar_distrib('a','z');
    std::uniform_int_distribution<uint32_t> input_size_distrib(1,200000);

    uint32_t input_size;
    uint8_t alphabet_size;
    std::string input;
    std::string input_reverted;
    uint32_t max_pattern_length;
    uint32_t num_queries;

    auto start_time = now();

    // generate random input strings and test all operations until one hour has passed
    while (time_diff_min(start_time,now()) < 60) {
        // choose a random input length
        input_size = input_size_distrib(gen);
        input.resize(input_size);

        // choose a random alphabet size
        alphabet_size = alphabet_size_distrib(gen);

        // choose a random alphabet
        std::vector<uint8_t> alphabet;
        uint8_t uchar;
        for (uint32_t i=0; i<alphabet_size; i++) {
            do {uchar = uchar_distrib(gen);} while (contains(alphabet,uchar));
            alphabet.push_back(uchar);
        }

        // choose a random input based on the alphabet
        std::uniform_int_distribution<uint8_t> char_idx_distrib(0,alphabet_size-1);
        double avg_input_rep_length = 1.0+avg_input_rep_length_distrib(gen);
        uint8_t cur_uchar = alphabet[char_idx_distrib(gen)];
        for (uint32_t i=0; i<input_size; i++) {
            if (prob_distrib(gen) < 1/avg_input_rep_length) cur_uchar = alphabet[char_idx_distrib(gen)];
            input[i] = uchar_to_char(cur_uchar);
        }

        // build the index
        r_index<> index(input,true);

        // generate patterns from the input and test count- and locate queries
        std::uniform_int_distribution<uint32_t> pattern_pos_distrib(0,input_size-1);
        max_pattern_length = std::min<uint32_t>(10000,std::max<uint32_t>(100,input_size/1000));
        std::uniform_int_distribution<uint32_t> pattern_length_distrib(1,max_pattern_length);
        num_queries = std::min<uint32_t>(10000,std::max<uint32_t>(1000,input_size/100));
        uint32_t pattern_pos;
        uint32_t pattern_length;
        std::string pattern;
        std::vector<uint64_t> correct_occurrences;
        std::vector<uint64_t> occurrences;
        bool match;
        for(uint32_t cur_query=0; cur_query<num_queries; cur_query++) {
            pattern_pos = pattern_pos_distrib(gen);
            pattern_length = std::min<uint32_t>(input_size-pattern_pos,pattern_length_distrib(gen));
            pattern.resize(pattern_length);
            for (uint32_t i=0; i<pattern_length; i++) pattern[i] = input[pattern_pos+i];
            for (uint32_t i=0; i<=input_size-pattern_length; i++) {
                match = true;
                for (uint32_t j=0; j<pattern_length; j++) {
                    if (input[i+j] != pattern[j]) {
                        match = false;
                        break;
                    }
                }
                if (match) correct_occurrences.emplace_back(i);
            }
            auto range = index.count(pattern);
            EXPECT_EQ(range.second-range.first+1,correct_occurrences.size());
            occurrences = index.locate_all(pattern);
            std::sort(occurrences.begin(),occurrences.end());
            EXPECT_EQ(occurrences,correct_occurrences);
            correct_occurrences.clear();
        }
    }
}
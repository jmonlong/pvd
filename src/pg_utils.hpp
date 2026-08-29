#ifndef PVD_PG_UTILS_HPP_INCLUDED
#define PVD_PG_UTILS_HPP_INCLUDED

#include <iostream>
#include <sstream>
#include <vector>
#include "types_and_structs.hpp"
#include <map>
#include <zlib.h>

// splits a string by a provided separator into a vector
std::vector<std::string> split_str(std::string input_str, char sep='\t');

// parses the INFO fields with values of a VCF into a map
std::map<std::string, std::string> parse_info(std::string info_str);

// parses an AT like >1>2>3 into [>1, >2, >3]
std::vector<std::string> parse_at(std::string at_str);

// joins a string vector range into a string
std::string join_str(std::vector<std::string> str_vec, int start_pos, int end_pos);

// joins a string vector with a separator
std::string join_str_sep(std::vector<std::string> str_vec, std::string sep);

// get the reference position one base after a node (e.g. '>1234')
size_t get_var_pos(std::string node, std::map<size_t, node_t>& node_info);

// get the sequence of a node, e.g. '>1234'
std::string get_sequence(std::string node, std::map<size_t, node_t>& node_info);

// get the padding base for a node (i.e. the last base of its sequence)
std::string get_padding(std::string node, std::map<size_t, node_t>& node_info);

// reads a gzipped GFA file to extract node information: node ID, sequence,
// reference position (if on the reference)
void read_nodes(const std::string& filename, std::map<size_t, node_t>& node_info, bool read_full_gfa=false);

#endif

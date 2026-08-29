#include "pg_utils.hpp"

// splits a string by a provided separator into a vector
std::vector<std::string> split_str(std::string input_str, char sep) {
    std::vector<std::string> fields;
    std::istringstream iss(input_str);
    std::string field;
    while (std::getline(iss, field, sep)) {
        fields.push_back(field);
    }
    return (fields);
}

// parses the INFO fields with values of a VCF into a map
std::map<std::string, std::string> parse_info(std::string info_str) {
    std::map<std::string, std::string> infos;
    std::vector<std::string> infos_v = split_str(info_str, ';');
    for (std::string info: infos_v) {
        std::vector<std::string> info_tag = split_str(info, '=');
        if (info_tag.size() == 2) {
            infos[info_tag[0]] = info_tag[1];
        }
    }
    return (infos);
}

// parses an AT like >1>2>3 into [>1, >2, >3]
std::vector<std::string> parse_at(std::string at_str) {
    std::vector<std::string> al_trav;
    // traverse the string and keep track of latest </> position
    int last_sep_pos = 0;
    for (int str_ii=1; str_ii<at_str.size(); ++str_ii) {
        if (at_str[str_ii] == '<' || at_str[str_ii] == '>'){
            // add current node and update sep position
            al_trav.push_back(at_str.substr(last_sep_pos, str_ii - last_sep_pos));
            last_sep_pos = str_ii;
        }
    }
    // include last ongoing node
    al_trav.push_back(at_str.substr(last_sep_pos, at_str.size() - last_sep_pos));
    return al_trav;
}

// joins a string vector range into a string
std::string join_str(std::vector<std::string> str_vec, int start_pos, int end_pos) {
    std::string joined_str = "";
    for (int idx=start_pos; idx<end_pos; ++idx) {
        joined_str += str_vec[idx];
    }
    return (joined_str);
}

// joins a string vector with a separator
std::string join_str_sep(std::vector<std::string> str_vec, std::string sep) {
    std::string joined_str = "";
    for (int idx=0; idx<str_vec.size(); ++idx) {
        if (idx > 0){
            joined_str += sep;
        }
        joined_str += str_vec[idx];
    }
    return (joined_str);
}

// get the reference position one base after a node (e.g. '>1234')
size_t get_var_pos(std::string node, std::map<size_t, node_t>& node_info) {
    size_t node_id = std::stoll(node.substr(1, node.size()-1));
    size_t node_size = node_info[node_id].seq.size();
    size_t refpos = node_info[node_id].refpos;
    return (refpos + node_size + 1);
}

// get the sequence of a node, e.g. '>1234'
std::string get_sequence(std::string node, std::map<size_t, node_t>& node_info) {
    size_t node_id = std::stoll(node.substr(1, node.size()-1));
    std::string seq = node_info[node_id].seq;
    if (node[0] == '>') {
        return (seq);
    } else {
        // reverse complement
        std::string rc_seq = seq;
        rc_seq.reserve(seq.size());
        for (int ii=0; ii<seq.size(); ++ii) {
            if (seq[ii] == 'A') {
                rc_seq[seq.size()-ii-1] = 'T';
            } else if (seq[ii] == 'T') {
                rc_seq[seq.size()-ii-1] = 'A';
            } else if (seq[ii] == 'C') {
                rc_seq[seq.size()-ii-1] = 'C';
            } else if (seq[ii] == 'G') {
                rc_seq[seq.size()-ii-1] = 'G';
            }
        }
        return (rc_seq);
    }     
}

// get the padding base for a node (i.e. the last base of its sequence)
std::string get_padding(std::string node, std::map<size_t, node_t>& node_info) {
    std::string seq = get_sequence(node, node_info);
    return (seq.substr(seq.size() - 1, 1));
}

// reads a gzipped GFA file to extract node information: node ID, sequence,
// reference position (if on the reference)
void read_nodes(const std::string& filename, std::map<size_t, node_t>& node_info,
                bool read_full_gfa) {
    gzFile file = gzopen(filename.c_str(), "rb"); // Open the gzip file in read-binary mode
    if (!file) {
        std::cerr << "Error opening file: " << filename << '\n';
        return;
    }

    // read a buffer buffer allowing for very long lines
    char buffer_char[4096];
    std::string buffer_str;
    bool keep_going = true;
    int n_non_s = -1;
    while (keep_going) {
        int bytes_read = gzread(file, buffer_char, 4096);
        // if nothing was read, we've reached the end (or there's an error)
        if (bytes_read <= 0) {
            break;
        }
        // append the new characters to the string buffer
        buffer_str.append(buffer_char, bytes_read);
        // loop through all the lines in that buffer
        size_t buffer_pos;
        while ((buffer_pos = buffer_str.find('\n')) != std::string::npos) {
            // extract the next line in the buffer
            std::string line = buffer_str.substr(0, buffer_pos);
            // split by tab
            std::vector<std::string> values = split_str(line, '\t');
            // is it an S line?
            if (values[0] == "S") {
                // create a node object
                node_t node;
                node.id = std::stoll(values[1]);
                node.seq = values[2];
                // init reference position
                node.refpos = -1;
                // look for the SO tag
                for (size_t idx = 3; idx < values.size(); idx++) {
                    std::vector<std::string> tag_info = split_str(values[idx], ':');
                    if (tag_info[0] == "SO") {
                        // found the SO tag
                        node.refpos = std::stoll(tag_info[2]);
                        break;
                    }
                }
                node_info[node.id] = node;
                n_non_s = 0;
            } else {
                // should we stop?
                if (!read_full_gfa && (n_non_s > 1000 || values[0] == "W" || values[0] == "P")) {
                    // yes if we haven't seen a S in a while or we reached the path part
                    keep_going = false;
                    std::cerr << "Assuming all S lines have been read. Stopping.\n";
                }
                if (n_non_s >= 0) {
                    n_non_s++;
                }
            }

            // remove that line from the buffer
            buffer_str.erase(0, buffer_pos + 1);
        }
    }

    gzclose(file);
}


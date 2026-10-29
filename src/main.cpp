#include <iostream>
#include <vector>
#include <map>
#include <zlib.h>
#include <ctime>
#include <chrono>
#include <getopt.h>

#include "types_and_structs.hpp"
#include "pg_utils.hpp"
#include "variant_site.hpp"

void print_help() {
    std::cerr << "Usage: pvd -v <variants.vcf.gz> -g <pangenome.gfa.gz> [options]\n\n"
              << "  -g, --gfa FILE            Path to the gzipped GFA file (only S lines really used)\n"
              << "  -v, --vcf FILE            Path to the gzipped VCF file\n"
              << "  -f, --full-gfa            Read the full GFA when extracting node information. E.g. if S lines are not all together.\n"
              << "  -h, --help                      Print this help message\n";
}


int main(int argc, char* argv[]) {

    std::string input_vcf_fn;
    std::string input_gfa_fn;
    bool read_full_gfa = false;
    
    // Parse arguments
    int c;

    static struct option long_options[] = {
        {"vcf", required_argument, 0, 'v'},
        {"gfa", required_argument, 0, 'g'},
        {"full-gfa", no_argument, 0, 'f'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    while ((c = getopt_long(argc, argv, "v:g:fh", long_options, nullptr)) != -1) {
        switch (c) {
        case 'v': input_vcf_fn = optarg; break;
        case 'g': input_gfa_fn = optarg; break;
        case 'f': read_full_gfa = true; break;
        case '?':
            print_help();
            return EXIT_FAILURE;
        case 'h': 
            print_help(); 
            return EXIT_SUCCESS; 
        default:
            print_help();
            return EXIT_FAILURE;
        }
    }

    // no positional arguments
    if (optind < argc) {
        std::cerr << "Invalid argument.\n";
        print_help();
        return EXIT_FAILURE;
    }
    
    if (argc == 1) {
        print_help();
        return EXIT_FAILURE;
    }
    
    time_t now = time(0);
    std::cerr << ctime(&now) << " - Reading node information" << std::endl;
    
    // read node information
    std::map<size_t, node_t> node_info;
    read_nodes(input_gfa_fn, node_info);
    now = time(0);
    std::cerr << ctime(&now) << " - Information for " << node_info.size() << " nodes saved." << std::endl;

    // start reading the VCF file
    gzFile file = gzopen(input_vcf_fn.c_str(), "rb"); // Open the gzip file in read-binary mode    
    if (!file) {
        std::cerr << "Error opening file: " << input_vcf_fn << '\n';
        return 1;
    }

    // keep track of the current variant and sample names
    std::map<int, std::string> sample_names;
    Snarl* cur_snarl;
    std::string cur_id;

    // have we already inserted our VCF headers?
    bool is_header_done = false;
    
    // read buffer objects
    size_t buffer_size = 4096;
    // char buffer_char[buffer_size];
    char* buffer_char = new char[buffer_size];
    std::string buffer_str;
    bool expand_buffer = false;
    auto start_timer = std::chrono::high_resolution_clock::now();
    auto end_timer = std::chrono::high_resolution_clock::now();
    size_t n_was_simple = 0;
    while (true) {
        // read up to a whole buffer of character
        if (expand_buffer) {
            buffer_size = buffer_size * 2;
            delete[] buffer_char;
            buffer_char = new char[buffer_size];
        }
        if (!expand_buffer && buffer_size > 4096) {
            buffer_size = 4096;
            delete[] buffer_char;
            buffer_char = new char[buffer_size];
        }
        int bytes_read = gzread(file, buffer_char, buffer_size);
        expand_buffer = true;
        // if nothing was read, we've reached the end (or there's an error)
        if (bytes_read <= 0) {
            break;
        }
        // append the new characters to the string buffer
        buffer_str.append(buffer_char, bytes_read);
        // loop through all the lines in that buffer
        size_t buffer_pos;
        while ((buffer_pos = buffer_str.find('\n')) != std::string::npos) {
            expand_buffer = false;
            // extract the next line in the buffer
            std::string current_line = buffer_str.substr(0, buffer_pos);
            // process the VCF line
            // write headers
            if (current_line[0] == '#') {
                std::vector<std::string> line_v_eq = split_str(current_line, '=');
                if (line_v_eq[0] == "##INFO" || line_v_eq[0] == "##FORMAT") {
                    // skip INFOs and FORMAT, we'll make our own
                    // remove that line from the buffer
                    buffer_str.erase(0, buffer_pos + 1);
                    continue;
                }
                std::vector<std::string> line_v = split_str(current_line, '\t');
                if ((line_v_eq[0] == "##contig"  || line_v[0] == "##CHROM") && !is_header_done) {
                    // add our headers just before the contig definition or the last header line
                    std::cout << "##INFO=<ID=PS,Number=1,Type=String,Description=\"ID of variant corresponding to parent snarl\">" << std::endl <<
                        "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">" << std::endl;
                    // add out FORMAT header
                    std::cout << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
                    is_header_done = true;
                }
                if (line_v[0] == "#CHROM"){
                    // save the sample names
                    for (int idx = 9; idx < line_v.size(); ++idx){
                        sample_names[idx] = line_v[idx];
                    }
                }
                std::cout << current_line << std::endl;
                // remove that line from the buffer
                buffer_str.erase(0, buffer_pos + 1);
                continue;
            }

            // if not a header
            std::vector<std::string> line_v = split_str(current_line, '\t');
            std::map<std::string, std::string> infos = parse_info(line_v[7]);
            
            // skip if not LV=0
            if (infos.find("LV") != infos.end() && infos["LV"] != "0") {
                // remove that line from the buffer
                buffer_str.erase(0, buffer_pos + 1);
                continue;
            }
            
            // if we're seeing a different variant ID,
            // it means we're getting a new snarl
            // std::cerr << line_v[2] << std::endl;
            if (cur_id != line_v[2]) {
                // TODO don't do anything if AT too short (all <=3)
                // eventually write the current variant
                if (!cur_id.empty()) {
                    // check if SNP here?
                    if (cur_snarl->is_simple()) {
                        // easy, no need to look for breakpoints
                        n_was_simple++;
                        cur_snarl->write_out_simple();
                    } else {
                        // compare the traversals to look for breakpoints
                        cur_snarl->find_breakpoints();
                        cur_snarl->write_out(node_info);
                    }
                    end_timer = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration<double>(end_timer - start_timer).count();
                    if (duration > 10) {
                        std::cerr << cur_id << " took " << std::to_string(duration) + " s!" << std::endl;
                    }
                    delete(cur_snarl);
                }
                // start a new variant
                start_timer = std::chrono::high_resolution_clock::now();
                cur_snarl = new Snarl(line_v, infos["AT"]);
            } else {
                // same variant ID, just update with new ATs
                cur_snarl->add_ats(infos["AT"]);
            }

            // add genotype for the samples
            std::vector<std::string> fmt = split_str(line_v[8], ':');
            bool found_gt = false;
            for (int gt_idx=0; gt_idx<fmt.size(); ++gt_idx) {
                if (fmt[gt_idx] == "GT") {
                    for (int samp_idx=9; samp_idx < line_v.size(); ++samp_idx) {
                        std::string gt = split_str(line_v[samp_idx], ':')[gt_idx];
                        cur_snarl->add_sample(sample_names[samp_idx], gt);
                    }
                    found_gt = true;
                    break;
                }
            }
            if (!found_gt) {
                std::cerr << "Warning: no GT in FORMAT field" << std::endl;
            }
            
            // update the current ID
            cur_id = line_v[2];
            
            // remove that line from the buffer
            buffer_str.erase(0, buffer_pos + 1);
        }
    }

    // close the input gzip file
    gzclose(file); 
    
    // finish the last record
    if (!cur_id.empty()) {
        cur_snarl->find_breakpoints();
        cur_snarl->write_out(node_info);
        delete(cur_snarl);
    }

    now = time(0);
    std::cerr << ctime(&now) << " - Done." << std::endl;
    std::cerr << n_was_simple << " were simple snarls." << std::endl;

    return 0;
}

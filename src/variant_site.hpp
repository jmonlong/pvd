#ifndef PVD_VARIANT_SITE_HPP_INCLUDED
#define PVD_VARIANT_SITE_HPP_INCLUDED

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "pg_utils.hpp"
#include "types_and_structs.hpp"

// Class to save information about a top-level snarl (one or two LV=0 lines in the input VCF)
class Snarl {
public:
    Snarl() {};

    // fill basic information (chromosome, quality, filter) and a first set
    // of allele traversals
    Snarl(std::vector<std::string> vcf_line, std::string info_ats);

    // add more ATs to the Snarl object, e.g. when a second VCF record for the
    // same top-level snarl is read
    void add_ats(std::string new_ats_str);

    // add genotype information for samples, read from the GT in the VCF record
    void add_sample(std::string samp_name, std::string gt);

    // compare ATs that are seen in at least one individual and identify boundary
    // reference nodes
    void find_breakpoints();

    // check if this is a simple situation with nothing to do (simple SNP/indel)
    bool is_simple();
    
    void write_out(std::map<size_t, node_t>& node_info);
    void write_out_simple();
    
protected:
    // save contig name and other info that we want to write out
    std::string seqn;
    std::string qual;
    std::string filter;
    std::string id;
    std::vector<std::string> first_vcf_rec;
    // save the reference AT
    std::string ref_at;
    // to save memory on huge snarls, assign an index to ATs and then work
    // with that index instead of the very long string
    std::map<int, std::string> id_to_at;
    std::map<std::string, int> at_to_id;
    // the original ATs' index and their vector representation (node, orientation)
    // e.g. 1 (AT index for ">10>11>12" for example) -> ['>10', '>11', '>12']
    std::map<int, std::vector<std::string>> ats;
    // current ATs to use when adding the next samples
    // helps match an allele index from the GT field to its allele index in this object
    std::vector<int> current_ats;
    // sample order, used when writing decomposed VCF record
    std::vector<std::string> samp_order;
    std::map<std::string, bool> phased_samps;
    // how many times have we seen this allele index across all samples
    // e.g. 1 (allele index, see id_to_at/at_to_id) -> 4
    std::map<int, int> tot_at_counts;
    // allele count for each sample
    // e.g. "samp3" -> [3, 1]
    // -1 means missing (".")
    std::map<std::string, std::vector<int>> samp_at_gts;
    // breakpoints and both their reference and alternate alleles
    // a breakpoint is defined by two nodes forming a subpath in
    // the reference AT (e.g. ">10>12"). For each breakpoint ID, we save the
    // different alleles (subpaths like '>10>13>12' or 'ref') and their
    // corresponding full-length allele index
    std::map<std::string, std::map<std::string, std::vector<int>>> bkpts;
    // save the reference ATs for the nested variants (breakpoint ID)
    std::map<std::string, std::vector<std::string>> ref_sub_ats;
    // save if that variant was split into multiple VCF records
    // if that's the case, we shouldn't assume '.' means an actual missing copy anymore
    // because it could be just due to splitting the snarl into multiple records
    bool multi_record_variant;
};

#endif

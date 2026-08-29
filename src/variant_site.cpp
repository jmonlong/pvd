#include "variant_site.hpp"

// fill basic information (chromosome, quality, filter) and a first set
// of allele traversals
Snarl::Snarl(std::vector<std::string> vcf_line, std::string info_ats) {
    seqn = vcf_line[0];
    qual = vcf_line[5];
    filter = vcf_line[6];
    id = vcf_line[2];
    // parse the ATs
    std::vector<std::string> orig_ats = split_str(info_ats, ',');
    ref_at = orig_ats[0];
    current_ats.clear();
    // saves which allele index corresponds to which AT in the current VCF record
    int at_idx = 0;
    for (std::string at: orig_ats) {
        id_to_at[at_idx] = at;
        at_to_id[at] = at_idx;
        ats[at_idx] = parse_at(at);
        current_ats.push_back(at_idx);
        at_idx++;
    }
    // for now only one record was used on this variants
    multi_record_variant = false;
}

// add more ATs to the Snarl object, e.g. when a second VCF record for the
// same top-level snarl is read
void Snarl::add_ats(std::string new_ats_str) {
    // if we're doing this, it means that the info was split in at least two records
    multi_record_variant = true;
    // split the ATs and update current ATs to use when adding the next samples
    std::vector<std::string> new_ats = split_str(new_ats_str, ',');
    // add potentially new AT from a new VCF record
    if (new_ats[0] != ref_at) {
        std::cerr << "the reference AT is different than the previous variant!" << std::endl;
    }
    // add new AT or assigns their index
    int at_idx = id_to_at.size();
    current_ats.clear();
    for (std::string new_at: new_ats) {
        if (at_to_id.find(new_at) == at_to_id.end()){
            // new AT
            id_to_at[at_idx] = new_at;
            at_to_id[new_at] = at_idx;
            ats[at_idx] = parse_at(new_at);
            current_ats.push_back(at_idx);
            at_idx++;
        } else {
            // we've seen this one in the previous VCF record
            current_ats.push_back(at_to_id[new_at]);
        }
    }
}

// add genotype information for samples, read from the GT in the VCF record
void Snarl::add_sample(std::string samp_name, std::string gt) {
    // potentially we'll want to skip this VCF record
    bool skip_if_missing = false; 
    // first update the sample list/order (if a new samples/first VCF record)
    if (samp_at_gts.find(samp_name) == samp_at_gts.end()) {
        samp_order.push_back(samp_name);
    } else {
        // we have this sample, so we're combining a second VCF record
        // check if the old GT was ./. and if so replace it
        bool all_missing = true;
        for (int at_idx: samp_at_gts[samp_name]) {
            if (at_idx != -1){
                all_missing = false;
                break;
            }
        }
        if (all_missing) {
            // remove the information from the previous VCF record
            phased_samps[samp_name] = false;
            samp_at_gts[samp_name].clear();
        } else {
            // skip this new VCF record but check that it is indeed ./.
            skip_if_missing = true;
        }
    }
    // for each allele of the genotype, update the AT genotypes
    char gt_sep = '/';
    if (gt.find('|') != std::string::npos) {
        phased_samps[samp_name] = true;
        gt_sep = '|';
    }
    std::vector<std::string> gt_als = split_str(gt, gt_sep);
    if (skip_if_missing) {
        bool all_missing = true;
        for (std::string al: gt_als) {
            if (al != ".") {
                all_missing = false;
            }
        }
        if (all_missing) {
            // expected, all good, just skip it
            return;
        } else {
            // potentially a problem...
            std::cerr << "Two VCFs record have GT information for " << samp_name << " " << id << std::endl;
        }
    }
    for (std::string al: gt_als) {
        if (al == ".") {
            // if missing, no need to check the AT and other info
            // just add a "fake" index to remember the (potential) phasing
            samp_at_gts[samp_name].push_back(-1);
            continue;
        }
        // convert to the AT index of the allele
        int at_idx = current_ats[std::stoi(al)];
        // update the AT genotype for that sample
        samp_at_gts[samp_name].push_back(at_idx);
        // update the total allele counts
        if (tot_at_counts.find(at_idx) == tot_at_counts.end()) {
            tot_at_counts[at_idx] = 0;
        }
        tot_at_counts[at_idx]++;
    }
}

// compare ATs that are seen in at least one individual and identify boundary
// reference nodes
void Snarl::find_breakpoints() {
    // first prepare reference nodes as a dict matching node to
    // its preceding node in the reference path
    std::map<std::string, std::string> ref_nodes;
    std::string prev_node = "";
    for (std::string node: ats[at_to_id[ref_at]]) {
        ref_nodes[node] = prev_node;
        prev_node = node;
    }
    // split alternate ATs at every discordant reference node
    for (auto [at_idx, at_count]: tot_at_counts) {
        // note: no need to require at_count>0 because the AT won't be in tot_at_counts
        if (id_to_at[at_idx] == ref_at) {
            // no need to compare the reference path to itself
            continue;
        }
        // work with the parsed AT ([>1, >2, <3, ...]
        std::vector<std::string> atp = ats[at_idx];
        // compare alternate AT with reference AT
        prev_node = "";
        int bk_start = 0;
        int bk_end = 0;
        for (int node_idx=0; node_idx<atp.size(); ++node_idx) {
            std::string node = atp[node_idx];
            if (ref_nodes.find(node) != ref_nodes.end()) {
                // current node is in the ref AT
                // was the preceding node different in the ref AT?
                if (prev_node != ref_nodes[node]) {
                    // yes then we've just found the second breakpoint
                    bk_end = node_idx;
                }
                // check if we were within a variant
                // we don't count (for now) cases where there is a variant but it's a cycle/duplication
                if (bk_start != bk_end && atp[bk_start] != atp[bk_end]) {
                    // both breakpoints are different, make an ID
                    std::string bkid = atp[bk_start] + atp[bk_end];
                    // save/update potential alleles between these breakpoints
                    if (bkpts.find(bkid) == bkpts.end()) {
                        bkpts[bkid]["ref"]; // not sure it's the right way to init
                    }
                    // concatenate string
                    // TODO make this more efficient? like reserve memory first?
                    std::string alt_at = join_str(atp, bk_start, bk_end + 1);
                    if (bkpts[bkid].find(alt_at) == bkpts[bkid].end()) {
                        bkpts[bkid][alt_at]; // not sure it's the right way to init
                    }
                    bkpts[bkid][alt_at].push_back(at_idx);
                }
                // we're (maybe now) following the reference path
                bk_start = node_idx;
                bk_end = node_idx;
                prev_node = node;
            } else {
                // if we're on non-reference nodes, move bk_end
                bk_end = node_idx;
                prev_node = node;
            }
        }
    }
    // to fill up the AT corresponding to the "reference" allele of each
    // identified breakpoint ('ref' in bkpts) we need to check each one to make sure they go
    // from on boundary node to the other using only reference edges
    // first index reference edges in all full ats
    std::map<std::string, std::vector<int>> ref_edge_to_at_idx;
    for (auto [at_idx, al_c]: tot_at_counts) {
        if (id_to_at[at_idx] == ref_at) {
            // no need to compare the reference path to itself
            continue;
        }
        std::vector<std::string> atp = ats[at_idx];
        for (int atp_ii=1; atp_ii<atp.size(); ++atp_ii) {
            // save if node is a ref node and its previous node was the previous
            // node on the reference AT
            if(ref_nodes.find(atp[atp_ii]) != ref_nodes.end() && ref_nodes[atp[atp_ii]] == atp[atp_ii-1]) {
                ref_edge_to_at_idx[atp[atp_ii-1] + atp[atp_ii]].push_back(at_idx);
            }
        }
    }
    // check the reference alleles for each breakpoints
    for (auto [bkid, bkpt_o]: bkpts){
        // parse the breakpoint id into two nodes (e.g. ['>10', '>12'])
        std::vector<std::string> bk_nd = parse_at(bkid);
        // first look for the reference sub-AT
        // start from second node and try to reach first node
        std::vector<std::string> ref_sub_at;
        ref_sub_at.push_back(bk_nd[1]);
        if (ref_nodes.find(ref_sub_at[0]) == ref_nodes.end()) {
            std::cerr << ref_sub_at[0] << " not in ref_nodes" << std::endl;
        }
        while (ref_nodes[ref_sub_at[0]] != "" && ref_sub_at[0] != bk_nd[0]) {
            ref_sub_at.insert(ref_sub_at.begin(), ref_nodes[ref_sub_at[0]]);
        }
        // try the other way if we didn't find a reference path
        if (ref_sub_at[0] != bk_nd[0]) {
            ref_sub_at.clear();
            ref_sub_at.push_back(bk_nd[0]);
            if (ref_nodes.find(ref_sub_at[0]) == ref_nodes.end()) {
                std::cerr << ref_sub_at[0] << " not in ref_nodes" << std::endl;
            }
            while (ref_nodes[ref_sub_at[0]] != "" && ref_sub_at[0] != bk_nd[1]) {
                ref_sub_at.insert(ref_sub_at.begin(), ref_nodes[ref_sub_at[0]]);
            }
            if (ref_sub_at[0] != bk_nd[1]) {
                std::cerr <<  "Error: no ref path for breakpoint " << bkid << std::endl;
            }
        }
        // add ref sub-AT in bkpts (ugly...)
        // doesn't work anymore... but did we actually use it?
        // ?? don't understand the comment above...
        ref_sub_ats[bkid] = ref_sub_at;
        bkpts[bkid]["ref"].push_back(at_to_id[ref_at]);
        // list ref-to-ref edges involving the variant boundaries
        std::string ref_edge_1 = ref_sub_at[0] + ref_sub_at[1];
        std::string ref_edge_2 = ref_sub_at[ref_sub_at.size()-2] + ref_sub_at[ref_sub_at.size()-1];
        // find alleles going through that ref sub-path
        std::vector<int> als_2 = ref_edge_to_at_idx[ref_edge_2];
        for (int at_idx: ref_edge_to_at_idx[ref_edge_1]) {
            if (std::find(als_2.begin(), als_2.end(), at_idx) != als_2.end()) {
                // this allele is present in both
                bkpts[bkid]["ref"].push_back(at_idx);
            }
        }
    }
}

void Snarl::write_out(std::map<size_t, node_t>& node_info) {
    // TODO don't change if not decomposition was needed, for example if only two original AT of length 3 (simple SNP case)
    // save the records to write in a list first that will be sorted by
    // position to make sure the output VCF is still sorted
    std::vector<std::pair<size_t, std::string>> vcf_recs;
    for (auto [bkid, bkinfo]: bkpts) {
        // write VCF record for this breakpoint
        // position from first ref node
        std::vector<std::string> ref_sub_at_p = ref_sub_ats[bkid];
        std::string ref_node_s = ref_sub_at_p[0];
        int pos = get_var_pos(ref_node_s, node_info);
        // and allele sequences
        std::string ref_seq = "";
        for (int node_idx=1; node_idx < ref_sub_at_p.size()-1; ++node_idx) {
            ref_seq += get_sequence(ref_sub_at_p[node_idx], node_info);
        }
        bool add_padding = ref_seq == "";
        std::vector<std::string> alt_seqs;
        // also prepare the new INFO fields
        std::string info_at = join_str(ref_sub_at_p, 0, ref_sub_at_p.size());
        // and save the allele ID to use for the sample GTs for each AT
        std::map<int, int> alt_at_o;
        int alt_ii = 1;
        // loop through each "allele" and their corresponding AT
        for (auto [alt_at, alt_als]: bkinfo) {
            if (alt_at == "ref") {
                // reference ATs
                for (int orig_at_idx: alt_als) {
                    alt_at_o[orig_at_idx] = 0;
                }
                // we skip the rest which is about ALT ATs and sequences
                continue;
            }
            std::vector<std::string> alt_at_p = parse_at(alt_at);
            // find ALT sequence
            std::string alt_seq = "";
            for (int node_idx=1; node_idx < alt_at_p.size()-1; ++node_idx) {
                alt_seq += get_sequence(alt_at_p[node_idx], node_info);
            }
            if (alt_seq == "") {
                add_padding = true;
            }
            alt_seqs.push_back(alt_seq);
            // remember which AT correspond to that allele index (for GT)
            for (int orig_at_idx: alt_als) {
                alt_at_o[orig_at_idx] = alt_ii;
            }
            // update the AT INFO field
            info_at += "," + alt_at;
            // increment ALT index
            alt_ii++;
        }
        // do we need a padding base?
        if (add_padding) {
            // get the last nucleotide of the first ref node
            std::string pad_seq = get_padding(ref_node_s, node_info);
            // add it to all sequences
            ref_seq = pad_seq + ref_seq;
            for (int alt_idx=0; alt_idx<alt_seqs.size(); ++alt_idx) {
                alt_seqs[alt_idx] = pad_seq + alt_seqs[alt_idx];
            }
            // decrement position
            pos += -1;
        }
        // combine ALT allele sequences
        std::string alt_seqs_j = join_str_sep(alt_seqs, ",");
        // genotypes for all samples
        std::string gts = "";
        for (std::string sampname: samp_order) {
            std::string gt = "";
            int ngt = 0;
            for (int samp_at_idx: samp_at_gts[sampname]) {
                // add appropriate GT seperator
                if (ngt > 0) {
                    if (phased_samps.find(sampname) != phased_samps.end()) {
                        gt += "|";
                    } else {
                        gt += "/";
                    }
                }
                // add appropriate allele index
                if (alt_at_o.find(samp_at_idx) != alt_at_o.end()) {
                    gt += std::to_string(alt_at_o[samp_at_idx]);
                } else {
                    gt += ".";
                }
                ngt++;
            }
            // // TODO ideally record sample ploidy and add just enough missing below, when needed
            // while (ngt < 2 && multi_record_variant) {
            //     if (ngt > 0) {
            //         gt += "/";
            //     }
            //     gt += ".";
            //     ngt++;
            // }
            if (ngt > 2) {
                std::cerr << "More than 2 alleles in sample " << sampname << " " << bkid << ": " << gt << std::endl;
            }
            if (gts != ""){
                gts += "\t" + gt;
            } else {
                gts += gt;
            }
        }
        // prepare record
        std::string info_val = "AT=" + info_at + ";PS=" + id;
        std::string out_vcf_rec = seqn + "\t" + std::to_string(pos) + "\t" + bkid + "\t" + ref_seq + "\t" + alt_seqs_j + "\t" + qual + "\t" + filter + "\t" + info_val + "\tGT\t" + gts;
        vcf_recs.push_back(std::make_pair(pos, out_vcf_rec));
    }

    std::sort(vcf_recs.begin(), vcf_recs.end(), [&] (std::pair<int, std::string> pa, std::pair<int, std::string> pb) {
        return (pa.first < pb.first);
    });
        
    for (auto [pos, vcf_out]: vcf_recs) {
        std::cout << vcf_out << std::endl;
    }
}

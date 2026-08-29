#ifndef types_and_structs
#define types_and_structs

// stores ID, sequence, and reference position of a node
struct node_t {
    size_t id;
    std::string seq;
    size_t refpos; // -1 if not on the reference path
};

#endif

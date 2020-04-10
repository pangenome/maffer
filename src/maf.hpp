#include <iostream>
#include "xg.hpp"

namespace maffer {

using namespace handlegraph;

struct segment_t {
    handle_t start;
    handle_t end;
};

struct path_trav_t {
    bool is_rev = false;
    uint64_t start = 0;
    uint64_t end = 0;
};

void write_maf(std::ostream& out, const xg::XG& graph);

}

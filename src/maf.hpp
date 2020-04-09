#include <iostream>
#include "xg.hpp"

namespace maffer {

struct segment_t {
    handle_t start;
    handle_t end;
};

void write_maf(std::ostream& out, const xg::XG& graph);

}

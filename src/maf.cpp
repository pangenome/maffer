#include "maf.hpp"

namespace maffer {

void write_maf(std::ostream& out, const xg::XG& graph) {
    // the algorithm is really simple
    // we find the segments
    // these are the maximal ranges in the vectorized pangenome space where we don't see any paths with non-monotonic positional changes
    std::vector<segment_t> segments;
    segment_t* curr = nullptr;
    // our path trajectories
    std::unordered_map<path_handle_t, std::pair<bool, uint64_t>> path_traj_pos;
    graph.for_each_handle(
        [&](const handle_t& h) {
            // starting case
            if (curr == nullptr) {
                segments.emplace_back();
                curr = &segments.back();
                curr.start = h;
            }
            // determine if we should break
            uint64_t handle_length = graph.get_length(h);
            bool should_break = false;
            graph.for_each_step_on_handle(
                h,
                [&](const step_handle_t& step) {
                    path_handle_t path = graph.get_path_handle_of_step(step);
                    uint64_t pos = graph.get_position_of_step(step);
                    bool is_rev = (graph.get_handle_of_setp(step) != h);
                    auto f = path_traj_pos.find(path);
                    if (f != path_traj_pos.end()) {
                        auto& traj_pos = f->second;
                        if (traj_pos.first == is_rev && traj_pos.second == pos) {
                            traj_pos.second += handle_length;
                        } else {
                            should_beak = true;
                        }
                    } else {
                        path_traj_pos[path] = std::make_pair(is_rev, pos + handle_length);
                    }
                });

            // XXX something is wrong with the update, or we might end up with an incomplete range
            if (should_break) {
                segments.emplace_back();
                curr = &segments.back();
                curr.start = h;
            }
        });
    
    // we write the segments

    // node_vector_offset is our friend
    
    // for each, extract the MAF record
    // then compress it use a heuristic to eliminate SNPs and maybe MNPs
    
}

}

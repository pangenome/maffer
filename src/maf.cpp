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
    //handle_t last_handle;
    graph.for_each_handle(
        [&](const handle_t& h) {
            std::cerr << "id = " << graph.get_id(h) << std::endl;
            // starting case
            if (curr == nullptr) {
                segments.emplace_back();
                curr = &segments.back();
                curr->start = h;
                curr->end = h;
            }
            // determine if we should break
            uint64_t handle_length = graph.get_length(h);
            std::cerr << "handle length " << handle_length << std::endl;
            for (auto& p : path_traj_pos) {
                std::cerr << "traj_pos " << graph.get_path_name(p.first) << " " << p.second.first << " " << p.second.second << std::endl;
            }
            bool should_break = false;
            graph.for_each_step_on_handle(
                h,
                [&](const step_handle_t& step) {
                    path_handle_t path = graph.get_path_handle_of_step(step);
                    uint64_t pos = graph.get_position_of_step(step);
                    bool is_rev = (graph.get_handle_of_step(step) != h);
                    std::cerr << "path " << graph.get_path_name(path) << " " << is_rev << " " << pos << std::endl;
                    if (is_rev) {
                        pos += handle_length;
                    }
                    auto f = path_traj_pos.find(path);
                    if (f != path_traj_pos.end()) {
                        auto& traj_pos = f->second;
                        if (traj_pos.first == is_rev && traj_pos.second == pos) {
                            if (traj_pos.first) {
                                traj_pos.second -= handle_length;
                            } else {
                                traj_pos.second += handle_length;
                            }
                        } else {
                            std::cerr << "breaking at " << graph.get_path_name(path) << " got traj_pos "
                                      << traj_pos.first << " " << traj_pos.second << " but wanted " << is_rev << " " << pos << std::endl;
                            path_traj_pos.erase(f);
                            path_traj_pos[path] = std::make_pair(is_rev, (is_rev ? pos - handle_length : pos + handle_length));
                            should_break = true;
                        }
                    } else {
                        path_traj_pos[path] = std::make_pair(is_rev, (is_rev ? pos - handle_length : pos + handle_length));
                    }
                });
            if (should_break) {
                //path_traj_pos.clear();
                //path_traj_pos[path] = std::make_pair(is_rev, (is_rev ? pos - handle_length : pos + handle_length));
                //
                segments.emplace_back();
                curr = &segments.back();
                curr->start = h;
            }
            curr->end = h;
            //last_handle = h;
            return true;
        });
    //curr->end = last_handle;
    
    // we write the segments
    uint64_t i = 0;
    for (auto& segment : segments) {
        std::cerr << "segment " << i++ << " "
                  << graph.get_id(segment.start) << "@"
                  << graph.node_vector_offset(graph.get_id(segment.start))
                  << " - "
                  << graph.get_id(segment.end) << "@"
                  << graph.node_vector_offset(graph.get_id(segment.end)) + graph.get_length(segment.end) << std::endl;
    }

    // node_vector_offset is our friend
    
    // for each, extract the MAF record
    // then compress it use a heuristic to eliminate SNPs and maybe MNPs
    
}

}

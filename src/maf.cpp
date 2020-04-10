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
                    if (is_rev) { pos += handle_length; }
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
    uint64_t j = 0;
    // node_vector_offset is our friend
    for (auto& segment : segments) {
        // hard assumption: our id space is contiguous
        nid_t start_id = graph.get_id(segment.start);
        nid_t end_id = graph.get_id(segment.end);
        std::cerr << "segment " << j++ << " "
                  << start_id << "@"
                  << graph.node_vector_offset(start_id)
                  << " - "
                  << end_id << "@"
                  << graph.node_vector_offset(end_id) + graph.get_length(segment.end) << std::endl;
        // collect the path set
        std::unordered_map<path_handle_t, std::unordered_map<uint64_t, path_trav_t>> path_limits;
        for (nid_t i = start_id; i <= end_id; ++i) {
            handle_t h = graph.get_handle(i);
            uint64_t handle_length = graph.get_length(h);
            graph.for_each_step_on_handle(
                h,
                [&](const step_handle_t& step) {
                    path_handle_t path = graph.get_path_handle_of_step(step);
                    uint64_t pos = graph.get_position_of_step(step);
                    bool is_rev = (graph.get_handle_of_step(step) != h);
                    if (is_rev) { pos += handle_length; }
                    auto f = path_limits.find(path);
                    if (f != path_limits.end()) {
                        auto& path_travs = f->second;
                        auto q = path_travs.find(pos);
                        if (q != path_travs.end()) {
                            // update
                            auto path_trav = q->second;
                            path_trav.end += (is_rev ? -handle_length : handle_length);
                            path_travs.erase(q);
                            path_travs[path_trav.end] = path_trav;
                        } else {
                            uint64_t end_pos = (is_rev ? pos - handle_length : pos + handle_length);
                            path_travs[end_pos] = { is_rev, pos, end_pos };
                        }
                    } else {
                        uint64_t end_pos = (is_rev ? pos - handle_length : pos + handle_length);
                        path_limits[path][end_pos] = { is_rev, pos, end_pos };
                    }
                });
        }
        for (auto& p : path_limits) {
            auto& path = p.first;
            auto& travs = p.second;
            for (auto& t : travs) {
                auto& trav = t.second;
                std::cerr << graph.get_path_name(path)
                          << " " << (trav.is_rev ? "-" : "+") << " "
                          << trav.start << " " << trav.end << std::endl;
            }
        }
        // find the limits of each path
        // and its orientation in the range
    }
    
    // for each, extract the MAF record
    // then compress it use a heuristic to eliminate SNPs and maybe MNPs
    
}

}

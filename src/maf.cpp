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
            std::unordered_set<path_handle_t> seen_paths;
            graph.for_each_step_on_handle(
                h,
                [&](const step_handle_t& step) {
                    path_handle_t path = graph.get_path_handle_of_step(step);
                    bool seen_path = false;
                    if (seen_paths.count(path)) {
                        std::cerr << "breaking becaus" << std::endl;
                        seen_path = true;
                        should_break = true;
                    } else {
                        seen_paths.insert(path);
                    }
                    if (!seen_path) {
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
        uint64_t pangenome_start = graph.node_vector_offset(start_id);
        uint64_t pangenome_end = graph.node_vector_offset(end_id) + graph.get_length(segment.end);
        std::cerr << "segment " << j++ << " "
                  << start_id << "@" << pangenome_start
                  << " - "
                  << end_id << "@" << pangenome_end << std::endl;
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
                std::string gapped;
                auto& trav = t.second;
                std::cerr << graph.get_path_name(path)
                          << " " << (trav.is_rev ? "-" : "+") << " "
                          << trav.start << " " << trav.end << std::endl;
                // print the gapped sequence against the pangenome
                // find the pangenome position of trav.start
                if (!trav.is_rev) {
                    uint64_t path_pos = trav.start;
                    uint64_t last_seq_pos = pangenome_start;
                    while (path_pos < trav.end) {
                        handle_t handle = graph.get_handle_of_step(graph.get_step_at_position(path, path_pos));
                        uint64_t seq_pos = graph.node_vector_offset(graph.get_id(handle));
                        if (seq_pos >= last_seq_pos) {
                            gapped.append(std::string(seq_pos - last_seq_pos, '-'));
                        } else {
                            std::cerr << "looping middle " << seq_pos << " " << last_seq_pos << std::endl;
                        }
                        gapped.append(graph.get_sequence(handle));
                        uint64_t handle_length = graph.get_length(handle);
                        path_pos += handle_length;
                        last_seq_pos = seq_pos + handle_length;
                    }
                    if (pangenome_end >= last_seq_pos) {
                        gapped.append(std::string(pangenome_end - last_seq_pos, '-'));
                    } else {
                        std::cerr << "looping end " << pangenome_end << " " << last_seq_pos << std::endl;
                    }
                    std::cerr << gapped << std::endl;
                } else {
                    // blah
                }
                
                // pad with - from our segment start to there
                // walk the path, padding as we need
                // pad with - to the end of our segment

            }
        }
        // find the limits of each path
        // and its orientation in the range
    }
    
    // for each, extract the MAF record
    // then compress it use a heuristic to eliminate SNPs and maybe MNPs
    
}

}

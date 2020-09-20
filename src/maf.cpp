#include "maf.hpp"

namespace maffer {

//#define debug_maf

void write_maf(std::ostream& out, const xg::XG& graph, const char* filename) {
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
#ifdef debug_maf
            std::cerr << "id = " << graph.get_id(h) << std::endl;
#endif

            // starting case
            if (curr == nullptr) {
                segments.emplace_back();
                curr = &segments.back();
                curr->start = h;
                curr->end = h;
            }
            // determine if we should break
            uint64_t handle_length = graph.get_length(h);

#ifdef debug_maf
            std::cerr << "handle length " << handle_length << std::endl;
            for (auto& p : path_traj_pos) {
                std::cerr << "traj_pos " << graph.get_path_name(p.first) << " " << p.second.first << " " << p.second.second << std::endl;
            }
#endif

            bool should_break = false;
            std::unordered_set<path_handle_t> seen_paths;
            graph.for_each_step_on_handle(
                h,
                [&](const step_handle_t& step) {
                    path_handle_t path = graph.get_path_handle_of_step(step);
                    bool seen_path = false;

                    if (seen_paths.count(path)) {
#ifdef debug_maf
                        std::cerr << "breaking at " << graph.get_path_name(path)
                        << ": encountered the same node several times in the path" << std::endl;
#endif

                        seen_path = true;
                        should_break = true;
                    } else {
                        seen_paths.insert(path);
                    }
                    if (!seen_path) {
                        uint64_t pos = graph.get_position_of_step(step);
                        bool is_rev = (graph.get_handle_of_step(step) != h);

#ifdef debug_maf
                        std::cerr << "path " << graph.get_path_name(path) << " " << is_rev << " " << pos << std::endl;
#endif

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
#ifdef debug_maf
                                std::cerr << "breaking at " << graph.get_path_name(path) << ": got traj_pos "
                                << traj_pos.first << " " << traj_pos.second << " but wanted " << is_rev << " "
                                << pos << std::endl;
#endif

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

                segments.emplace_back();
                curr = &segments.back();
                curr->start = h;
            }
            curr->end = h;
            //last_handle = h;

#ifdef debug_maf
            std::cerr << std::endl;
#endif

            return true;
        });
    //curr->end = last_handle;

    // maf header
    out << "##maf version=1" << std::endl;
    out << "# maffer" << std::endl;
    out << "# input=" << filename << std::endl;
    out << std::endl;
    
    // we write the segments
    uint64_t j = 0;
    // node_vector_offset is our friend
    for (auto& segment : segments) {
        // hard assumption: our id space is contiguous
        nid_t start_id = graph.get_id(segment.start);
        nid_t end_id = graph.get_id(segment.end);
        uint64_t pangenome_start = graph.node_vector_offset(start_id);
        uint64_t pangenome_end = graph.node_vector_offset(end_id) + graph.get_length(segment.end);

#ifdef debug_maf
        std::cerr << "segment " << j++ << " "
                  << start_id << "@" << pangenome_start
                  << " - "
                  << end_id << "@" << pangenome_end << std::endl;
#endif

        // collect the path set
        bool contains_loops = false;
        std::unordered_map<path_handle_t, std::unordered_map<uint64_t, path_trav_t>> path_limits;
        for (nid_t i = start_id; i <= end_id; ++i) {
            handle_t h = graph.get_handle(i);
            uint64_t handle_length = graph.get_length(h);
            std::vector<path_pos_t> steps;
            std::unordered_set<path_handle_t> seen_paths;
            graph.for_each_step_on_handle(
                h,
                [&](const step_handle_t& step) {
                    path_handle_t path = graph.get_path_handle_of_step(step);
                    if (seen_paths.count(path)) {
                        contains_loops = true;
                    } else {
                        seen_paths.insert(path);
                    }
                    uint64_t pos = graph.get_position_of_step(step);
                    bool is_rev = (graph.get_handle_of_step(step) != h);
                    if (is_rev) { pos += handle_length; }
                    steps.push_back({path, pos, is_rev});
                });
            // in order to correctly extend, we need to process our steps in the correct path-relative order
            std::sort(steps.begin(), steps.end(),
                      [&graph](const path_pos_t& a,
                               const path_pos_t& b) {
                          auto& a_i = as_integer(a.path);
                          auto& b_i = as_integer(b.path);
                          return a_i < b_i || a_i == b_i && (a.is_rev ? a.pos > b.pos : a.pos < b.pos);
                      });
            for (auto& s : steps) {
                auto& path = s.path;
                auto& pos = s.pos;
                auto& is_rev = s.is_rev;
                auto f = path_limits.find(path);
                if (f != path_limits.end()) {
                    auto& path_travs = f->second;
                    auto q = path_travs.find(pos);
                    if (q != path_travs.end()) {
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
            }
        }
        
        // build our MAF records
        std::vector<maf_record_t> records;
        for (auto& p : path_limits) {
            auto& path = p.first;
            auto& travs = p.second;
            for (auto& t : travs) {
                auto& trav = t.second;
                records.emplace_back();
                auto& record = records.back();
                record.src = graph.get_path_name(path);

                uint64_t record_src_size = graph.get_path_length(path);

                // If the strand field is "-" then this is the start relative to the reverse-complemented source sequence
                record.start = std::to_string(trav.is_rev ? record_src_size - trav.start : trav.start);

                record.size = std::to_string(trav.is_rev ? trav.start - trav.end : trav.end - trav.start);
                record.is_rev = (trav.is_rev ? "-" : "+");
                record.src_size = std::to_string(record_src_size);
                // print the gapped sequence against the pangenome
                // find the pangenome position of trav.start
                if (!trav.is_rev) {
                    uint64_t path_pos = trav.start;
                    uint64_t last_seq_pos = pangenome_start;
                    while (path_pos < trav.end) {
                        handle_t handle = graph.get_handle_of_step(graph.get_step_at_position(path, path_pos));
                        uint64_t seq_pos = graph.node_vector_offset(graph.get_id(handle));
                        if (seq_pos >= last_seq_pos) {
                            record.text.append(std::string(seq_pos - last_seq_pos, '-'));
                        } else {

#ifdef debug_maf
                            std::cerr << "looping middle " << seq_pos << " " << last_seq_pos << std::endl;
#endif
                        }
                        record.text.append(graph.get_sequence(handle));
                        uint64_t handle_length = graph.get_length(handle);
                        path_pos += handle_length;
                        last_seq_pos = seq_pos + handle_length;
                    }
                    if (pangenome_end >= last_seq_pos) {
                        record.text.append(std::string(pangenome_end - last_seq_pos, '-'));
                    } else {
#ifdef debug_maf
                        std::cerr << "looping end " << pangenome_end << " " << last_seq_pos << std::endl;
#endif
                    }
                    //std::cout << gapped << std::endl;
                } else {
                    uint64_t path_pos = trav.start;
                    uint64_t last_seq_pos = pangenome_start;

#ifdef debug_maf
                    std::cerr << "rev trav " << trav.start << " " << trav.end << std::endl;
#endif
                    while (path_pos > trav.end) {
                        // evil -1 hack
                        handle_t handle = graph.flip(graph.get_handle_of_step(graph.get_step_at_position(path, path_pos-1)));
                        uint64_t seq_pos = graph.node_vector_offset(graph.get_id(handle));
#ifdef debug_maf
                        std::cerr << "last_seq_pos " << last_seq_pos << std::endl;
                        std::cerr << "seq_pos " << seq_pos << std::endl;
#endif
                        if (seq_pos >= last_seq_pos) {
                            record.text.append(std::string(seq_pos - last_seq_pos, '-'));
                        } else {
#ifdef debug_maf
                            std::cerr << "looping middle " << seq_pos << " " << last_seq_pos << std::endl;
#endif
                        }
                        record.text.append(graph.get_sequence(handle));
                        uint64_t handle_length = graph.get_length(handle);
                        path_pos -= handle_length;
                        last_seq_pos = seq_pos + handle_length;
                    }
                    if (pangenome_end >= last_seq_pos) {
                        record.text.append(std::string(pangenome_end - last_seq_pos, '-'));
                    } else {
#ifdef debug_maf
                        std::cerr << "looping end " << pangenome_end << " " << last_seq_pos << std::endl;
#endif
                    }
                    //std::cout << gapped << std::endl;
                }
            }
        }

        // determine block length
        size_t block_length = 0;
        for (auto& record : records) {
            block_length = std::max(record.text.size(), block_length);
        }
        if (!contains_loops) {
            for (auto& record : records) {
                //assert(record.text.size() == block_length);
            }
        } else {
            // pad looping blocks
            for (auto& record : records) {
                record.text.append(std::string(block_length - record.text.size(), '-'));
            }
        }

        uint64_t last_block_length = block_length;
        // compress SNPs in the MAF records
        // run through each, finding places where all sequences have a indel of a given size alternating
        uint64_t compress_iter_max = 10;
        for (uint64_t iter = 0; iter < compress_iter_max; ++iter) {
            std::vector<bool> gap_drops(block_length);
            for (uint64_t i = 0; i < block_length; ++i) {
                gap_drops[i] = false;
            }
            if (block_length >= 2) {
                for (uint64_t i = 0; i < block_length - 1; ++i) {
                    // simple filter; do all records have gaps in the current or next column?
                    uint64_t snp_gap_drop_count = 0;
                    for (auto& record : records) {
                        auto& text = record.text;
                        snp_gap_drop_count += (text.at(i) == '-' ^ text.at(i+1) == '-');
                    }
                    if (snp_gap_drop_count == records.size()) {
                        gap_drops[i] = true;
                        ++i;
                        gap_drops[i] = true;
                        ++i;
                    }
                }
            }
            // remove gaps
            for (auto& record : records) {
                auto& text = record.text;
                std::string ungapped_text;
                for (uint64_t i = 0; i < text.size(); ++i) {
                    if (gap_drops[i] && text[i] == '-') {
                    } else {
                        ungapped_text.push_back(text[i]);
                    }
                }
                text = ungapped_text;
            }
            block_length = 0;
            for (auto& record : records) {
                block_length = std::max(record.text.size(), block_length);
            }
            for (auto& record : records) {
                if (record.text.size() != block_length) {
                    std::cerr << "[maffer] failure to compress block" << std::endl;
                    assert(false);
                    //iter = compress_iter_max;
                }
            }
            if (block_length == last_block_length) {
#ifdef debug_maf
                std::cerr << "done" << std::endl;
#endif
                break;
            } else {
                last_block_length = block_length;
            }
        }
        
        // pad and write the MAF records
        // determine output widths for everything
        size_t max_src_length = 0;
        size_t max_start_length = 0;
        size_t max_size_length = 0;
        size_t max_is_rev_length = 0;
        size_t max_src_size_length = 0;
        size_t max_text_length = 0;
        for (auto& record : records) {
            max_src_length = std::max(max_src_length, record.src.size());
            max_start_length = std::max(max_start_length, record.start.size());
            max_size_length = std::max(max_size_length, record.size.size());
            max_is_rev_length = std::max(max_is_rev_length, record.is_rev.size());
            max_src_size_length = std::max(max_src_size_length, record.src_size.size());
            max_text_length = std::max(max_text_length, record.text.size());
        }
        // write and pad them
        out << "a " << "loops=" << (contains_loops?"true":"false") << std::endl;
        for (auto& record : records) {
            out << "s "
                << record.src << std::string(max_src_length - record.src.size(), ' ')
                << std::setw(max_start_length+1) << record.start //string(max_start_length-record.start.size(), ' ') << record.start
                << std::setw(max_size_length+1) << record.size
                << std::setw(max_is_rev_length+1) << record.is_rev
                << std::setw(max_src_size_length+1) << record.src_size
                << " " << record.text
                << std::endl;
        }
        out << std::endl;
        // find the limits of each path
        // and its orientation in the range
    }
    
    // for each, extract the MAF record
    // then compress it use a heuristic to eliminate SNPs and maybe MNPs
    
}

}

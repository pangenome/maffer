/** \file version_main.cpp
 *
 * Defines the "vg version" subcommand, which evaluates graphs and alignments.
 */


#include <omp.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include "args.hxx"
#include "sdsl/bit_vectors.hpp"
#include "xg.hpp"
#include "maf.hpp"

using namespace std;
using namespace xg;
using namespace maffer;

int main(int argc, char** argv) {
    args::ArgumentParser parser("maffer: extract MSAs from variation graphs");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<std::string> gfa_in(parser, "FILE", "index the graph in this GFA file", {'g', "gfa-in"});
    args::ValueFlag<std::string> xg_out(parser, "FILE", "write the resulting xg index to this file", {'o', "out"});
    args::ValueFlag<std::string> xg_in(parser, "FILE", "read the xg index from this file", {'i', "in"});
    args::ValueFlag<std::string> base(parser, "BASE", "use this basename for temporary files during build", {'b', "base"});
    args::Flag gfa_out(parser, "FILE", "write the graph in GFA to stdout", {'G', "gfa-out"});
    //args::ValueFlag<std::string> structure_html(parser, "FILE", "store the sdsl structure description", {'S', "structure"});
    args::ValueFlag<uint64_t> num_threads(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
    args::Flag validate(parser, "validate", "validate construction", {'V', "validate"});
    args::Flag maf_out(parser, "write-maf", "write MAF output to stdout", {'m', "write-maf"});
    //args::Flag debug(parser, "debug", "enable debugging", {'d', "debug"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    size_t n_threads = args::get(num_threads);
    if (n_threads) {
        omp_set_num_threads(args::get(num_threads));
    } else {
        omp_set_num_threads(1);
    }

    XG graph;
    if (!args::get(xg_in).empty()) {
        std::ifstream in(args::get(xg_in));
        graph.deserialize(in);
    } else if (!args::get(gfa_in).empty()) {
        graph.from_gfa(args::get(gfa_in), args::get(validate),
                       args::get(base).empty() ? args::get(gfa_in) : args::get(base));
    }

    if (!args::get(xg_out).empty()) {
        std::ofstream out(args::get(xg_out));
        unique_ptr<sdsl::structure_tree_node> structure;
        structure = unique_ptr<sdsl::structure_tree_node>(new sdsl::structure_tree_node("name", "type"));
        graph.serialize_and_measure(out, structure.get(), "xg");
    }

    if (args::get(gfa_out)) {
        graph.to_gfa(std::cout);
    } else if (args::get(maf_out) || !args::get(xg_in).empty()) {
        // default: MAF output if we get -i input, or -m flag
        write_maf(std::cout, graph);
    }

    return 0;
}

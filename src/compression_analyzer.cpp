#include <iostream>
#include <string>
#include <getopt.h>
#include <vector>
#include <time.h>
#include <cmath>
#include <atomic>

#include "bioparser/fastq_parser.hpp"
#include "biosoup/sequence.hpp"

#define VERSION "v0.1.9"

static int help_flag = 0;         /* Flag set by �--help�.    */
static int version_flag = 0;      /* Flag set by �--version�. */

std::atomic<std::uint32_t> biosoup::Sequence::num_objects{0};

void printFragmentsInfo(const std::vector<std::unique_ptr<biosoup::Sequence>>& fragments) {
    uint64_t length_sum = 0;
    std::vector<size_t> lengths(fragments.size());
    for (int i = 0; i < int(fragments.size()); i++) {
        lengths[i] = fragments[i]->data.size();
        length_sum += lengths[i];
    }
    sort(lengths.begin(), lengths.end(), std::greater<size_t>());

    uint64_t N50 = -1, tmp_sum = 0;
    for (int i = 0; i < int(fragments.size()); i++) {
        tmp_sum += lengths[i];
        if (tmp_sum * 2 >= length_sum) {
            N50 = lengths[i];
            break;
        }
    }
    std::cerr << "FASTQ fragments:\n";
    std::cerr << "Number of fragments: " << fragments.size() << '\n';
    std::cerr << "Average length: " << length_sum * 1.0 / fragments.size() << '\n';
    std::cerr << "N50 length: " << N50 << '\n';
    std::cerr << "Minimal length: " << lengths.back() << '\n';
    std::cerr << "Maximal length: " << lengths.front() << "\n\n";
}

/* Modificiran primjer https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html */
/* Pojasnjenje primjera https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Options.html */

const std::string HELP_MESSAGE = "compression_analyzer usage: \n\n"
                                 "flags: \n"
                                 "-h or --help     prints help message \n"
                                 "-v or --version  prints version      \n"
                                 "\ncompression_analyzer takes one FASTQ filename as a command line argument.\n";

int main (int argc, char **argv) {
    srand (time(NULL)); /* initialize random seed: */
    int c;              /* result variable for getopt_long function */

    while (1) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                {"help",    no_argument, &help_flag,    1},
                {"version", no_argument, &version_flag, 1},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long(argc, argv, "hv",
                        long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c) {
        case 0:
            break;

        case 'h':
            help_flag = 1;
            break;

        case 'v':
            version_flag = 1;
            break;
        
        case '?':
            /* getopt_long already printed an error message. */
            break;

        default:
            abort ();
        }
    }

    // For fast I/O
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if (help_flag) {
        std::cout << HELP_MESSAGE;
    } 
    if (version_flag) {
        std::cout << VERSION << std::endl;
    }
    
    if (optind < argc) {
        auto x = biosoup::Sequence();
        std::cout << "hej" << std::endl;
        auto fragment_parser = bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastqParser>(argv[optind]);
        std::cout << "hejxxxx:" << std::endl;

        //parse in chunks
        std::vector<std::unique_ptr<biosoup::Sequence>> fragments;
        std::uint32_t chunk_size = 500 * 1024 * 1024;  // 500 MB
        for (auto t = fragment_parser->Parse(chunk_size); !t.empty(); t = fragment_parser->Parse(chunk_size)) {
            fragments.insert(
                fragments.end(),
                std::make_move_iterator(t.begin()),
                std::make_move_iterator(t.end()));
        }
        printFragmentsInfo(fragments);
        std::cout << "heyyyy:" << std::endl;
    }

    return 0;
}
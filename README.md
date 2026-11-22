> Originally developed in 2021. Cleaned up and refactored in 2026 with AI assistance (Claude), mainly for readability, documentation, and minor bug fixes. Core design and algorithms for all problems are from the original 2021 work.

# Rosalind Bioinformatics Stronghold Solutions

My solutions to problems from [Rosalind](https://rosalind.info/), mostly the **Bioinformatics Stronghold**
(string algorithms, sequence alignment, genome assembly, phylogeny, probability), plus the hidden Markov model
problems `ba10a`–`ba10k` from the Bioinformatics Textbook Track.

- **Python** — one module per topic in [`python/`](python/), one or more functions per problem.
- **Java** — four alignment problems (`laff`, `smgb`, `oap`, `ksim`) in [`java/`](java/), written for speed on the
  full-size datasets (see [Why some problems are in Java](#why-some-problems-are-in-java)).

## Repository layout

```
.
├── python/
│   ├── strings.py            # string, counting and motif problems   (dna … mprt)
│   ├── alignment.py          # sequence alignment                    (lcsq … osym)
│   ├── assembly.py           # genome assembly                       (grph, long, corr, dbru, pcov, gasm, grep, asmq)
│   ├── rna_structure.py      # RNA secondary structure               (pmch, mmch, cat, motz, rnas)
│   ├── phylogeny.py          # trees, distance matrices, sets        (pdst, nwck, nkew, tree, inod, seto)
│   ├── hmm.py                # hidden Markov models                  (ba10a … ba10k)
│   └── utils.py              # FASTA reader, output formatters, BLOSUM62/PAM250 loaders, HMM matrix writer
├── java/
│   ├── Laff.java             # local alignment with affine gap penalty (BLOSUM62)
│   ├── Smgb.java             # semiglobal alignment
│   ├── Oap.java              # overlap alignment
│   └── Ksim.java             # all substrings of s within edit distance k of t
└── data/
    ├── inputs/               # rosalind_<id>.txt  -- problem datasets
    ├── scoring_matrices/     # BLOSUM62.txt, PAM250.txt
    └── outputs/              # results written by the programs (e.g. laff_rslt.txt)
```

**All file paths in the code are relative to the repository root, so run everything from there.**

## Requirements

| For | You need |
| --- | --- |
| `python/` | Python 3, `numpy`, `pandas`, `scipy`, `biopython` |
| `lcsm` (`dna_lcs`) | [`suffix_trees`](https://pypi.org/project/suffix-trees/) |
| `mprt` (`glycosylation_motif`) | `requests`, `beautifulsoup4`, and internet access (protein sequences are downloaded from UniProt) |
| `java/` | A JDK. Use **Java 11 or newer**: `Smgb.java` calls `String.repeat`, which does not exist in Java 8 |

`suffix_trees`, `requests` and `beautifulsoup4` are only imported inside the two functions that use them, so the rest of
`strings.py` works without them.

## Usage

### Python

Each module is a set of solution functions. At the bottom is an `if __name__ == "__main__":` block with the demo calls
for that module (on a sample string, or on the matching file in `data/inputs/`). Importing a module does not run them.

```bash
cd Rosalind_Bioinfo_Stronghold_Solutions

python python/alignment.py          # run one module's demo calls

# use a function from your own code or a REPL (PYTHONPATH makes the modules importable)
PYTHONPATH=python python -c "from alignment import edit_distance; edit_distance('PLEASANTLY', 'MEANLY')"   # prints 5
```

A module's demos run top to bottom and stop at the first one that fails, for example a demo whose dataset is not in the
repository. See [Known issues](#known-issues).

`utils.py` reads the scoring matrices when it is imported, so `data/scoring_matrices/` must be reachable from the
current directory (another reason to run from the repository root).

### Java

```bash
cd Rosalind_Bioinfo_Stronghold_Solutions
mkdir -p build && javac -d build java/*.java

java -cp build Laff                 # reads data/inputs/rosalind_laff.txt,  writes data/outputs/laff_rslt.txt
java -cp build Smgb                 # reads data/inputs/rosalind_smgb.txt,  writes data/outputs/smgb_rslt.txt
java -cp build Oap                  # reads data/inputs/rosalind_oap.txt,   writes data/outputs/oap_rslt.txt
java -cp build Ksim                 # reads data/inputs/rosalind_ksim.txt,  writes data/outputs/ksim_rslt.txt
```

Each program first runs a tiny built-in test case, then the real dataset, and prints the alignment time.
`Ksim` on the full dataset can take several minutes. Compiled `.class` files are not tracked (see `.gitignore`).

Output formats: `Laff` — score, then the two aligned local substrings; `Smgb` / `Oap` — score, then the two aligned
strings; `Ksim` — one `start length` pair per line.

### Why some problems are in Java

`laff`, `smgb` and `oap` are also solved in Python (`alignment.py`), which is fine for the Rosalind sample inputs but too
slow for the real datasets, where each sequence has thousands of residues (the docstrings in `alignment.py` note this
for proteins above about 1,000 aa and DNA of 1,000–10,000 bp). The Java programs are equivalent implementations written
for speed. For example, on the 8,388 × 8,356 aa `laff` dataset the Java program takes about 3 s (measured on an idle
node with Java 8; the program is single-threaded and the time varies with machine load), while the Python function is
estimated to need about 26–27 minutes (extrapolated from timings on 200 × 200, 400 × 400 and 800 × 800 slices, all at
roughly 22–23 µs per cell). Java is therefore several hundred times faster here. `ksim` exists only in Java.

## Data

- `data/inputs/` holds the datasets I downloaded from Rosalind, named `rosalind_<id>.txt`. Rosalind generates a
  different dataset for every download, so these will not match what you get for the same problem.
  `rosalind_ksim_test.txt` is a small extra test input for `Ksim` (no program reads it by default).
- The Python modules also refer to datasets that are **not** in this repository:
  `asmq`, `corr`, `grep`, `long` (as `rosalind_long-2.txt`), `mprt`, and `ba10a`–`ba10k`.
  Download them from Rosalind and put them in `data/inputs/` under the same names.
- Python demo calls write scratch output files (`test.txt`, `result.txt`, `*_rslt.txt`, `ba10g_input1.txt`) into
  `data/outputs/`.

## Problem index

Titles are the official Rosalind titles. Function names are in the module named in each section heading.

### Strings, counting and combinatorics — `python/strings.py`

| Problem | Title | Function(s) |
| --- | --- | --- |
| [dna](https://rosalind.info/problems/dna/) | Counting DNA Nucleotides | `count_base` |
| [rna](https://rosalind.info/problems/rna/) | Transcribing DNA into RNA | `transcribe` |
| [revc](https://rosalind.info/problems/revc/) | Complementing a Strand of DNA | `reverse_complement` |
| [fib](https://rosalind.info/problems/fib/) | Rabbits and Recurrence Relations | `fibonacci` |
| [fibd](https://rosalind.info/problems/fibd/) | Mortal Fibonacci Rabbits | `mortal_fibonacci` |
| [gc](https://rosalind.info/problems/gc/) | Computing GC Content | `gc_content`, `get_highest_gc` |
| [hamm](https://rosalind.info/problems/hamm/) | Counting Point Mutations | `hamming_dist` |
| [perm](https://rosalind.info/problems/perm/) | Enumerating Gene Orders | `rearrange` |
| [sign](https://rosalind.info/problems/sign/) | Enumerating Oriented Gene Orderings | `signed_permutation` |
| [lexf](https://rosalind.info/problems/lexf/) | Enumerating k-mers Lexicographically | `lexicograph` |
| [kmer](https://rosalind.info/problems/kmer/) | k-Mer Composition | `kmer_composition` |
| [iprb](https://rosalind.info/problems/iprb/) | Mendel's First Law | `mendel_prob` |
| [lia](https://rosalind.info/problems/lia/) | Independent Alleles | `prob_lia` |
| [rstr](https://rosalind.info/problems/rstr/) | Matching Random Motifs | `prob_rstr` |

### Translation and motif finding — `python/strings.py`

| Problem | Title | Function(s) |
| --- | --- | --- |
| [prot](https://rosalind.info/problems/prot/) | Translating RNA into Protein | `translate` |
| [mrna](https://rosalind.info/problems/mrna/) | Inferring mRNA from Protein | `infer_mrna` |
| [orf](https://rosalind.info/problems/orf/) | Open Reading Frames | `find_orfs` |
| [splc](https://rosalind.info/problems/splc/) | RNA Splicing | `exon_protein` |
| [subs](https://rosalind.info/problems/subs/) | Finding a Motif in DNA | `find_motif` |
| [cons](https://rosalind.info/problems/cons/) | Consensus and Profile | `consensus_string` |
| [kmp](https://rosalind.info/problems/kmp/) | Speeding Up Motif Finding | `failure_array`, plus a Knuth–Morris–Pratt search (`find_lps`, `kmp`) |
| [sseq](https://rosalind.info/problems/sseq/) | Finding a Spliced Motif | `find_spliced_motif` |
| [revp](https://rosalind.info/problems/revp/) | Locating Restriction Sites | `find_palindrome`; brute-force version `find_palindrome_brute_force` |
| [lcsm](https://rosalind.info/problems/lcsm/) | Finding a Shared Motif | `dna_lcs` (generalized suffix tree) |
| [mprt](https://rosalind.info/problems/mprt/) | Finding a Protein Motif | `glycosylation_motif` (N-glycosylation motif, UniProt download) |

### RNA secondary structure — `python/rna_structure.py`

| Problem | Title | Function(s) |
| --- | --- | --- |
| [pmch](https://rosalind.info/problems/pmch/) | Perfect Matchings and RNA Secondary Structures | `perfect_matching` |
| [mmch](https://rosalind.info/problems/mmch/) | Maximum Matchings and RNA Secondary Structures | `maximum_matching` |
| [cat](https://rosalind.info/problems/cat/) | Catalan Numbers and RNA Secondary Structures | `noncrossing` (memoized recursion) |
| [motz](https://rosalind.info/problems/motz/) | Motzkin Numbers and RNA Secondary Structures | `noncrossing_all` |
| [rnas](https://rosalind.info/problems/rnas/) | Wobble Bonding and RNA Secondary Structures | `noncrossing_with_wobble` |

### Sequence alignment — `python/alignment.py`

| Problem | Title | Function(s) |
| --- | --- | --- |
| [lcsq](https://rosalind.info/problems/lcsq/) | Finding a Shared Spliced Motif | `lcs_with_length_table` (keeps a table of lengths next to a table of strings) |
| [scsp](https://rosalind.info/problems/scsp/) | Interleaving Two Motifs | `shortest_common_supersequence`, built on `lcs_string_table` (a single table of strings) |
| [edit](https://rosalind.info/problems/edit/) | Edit Distance | `edit_distance` |
| [edta](https://rosalind.info/problems/edta/) | Edit Distance Alignment | `edit_reconstruct` |
| [ctea](https://rosalind.info/problems/ctea/) | Counting Optimal Alignments | `count_optimal_alignments_memo` (memoized recursion) and `count_optimal_alignments_dp` (dynamic programming) |
| [glob](https://rosalind.info/problems/glob/) | Global Alignment with Scoring Matrix | `global_alignment_score_linear` (BLOSUM62) |
| [gaff](https://rosalind.info/problems/gaff/) | Global Alignment with Scoring Matrix and Affine Gap Penalty | `global_alignment_affine` |
| [gcon](https://rosalind.info/problems/gcon/) | Global Alignment with Constant Gap Penalty | `global_alignment_score_constant` |
| [loca](https://rosalind.info/problems/loca/) | Local Alignment with Scoring Matrix | `local_alignment_linear` (Smith–Waterman, PAM250) |
| [laff](https://rosalind.info/problems/laff/) | Local Alignment with Affine Gap Penalty | `local_alignment_affine`; Java: `java/Laff.java` |
| [smgb](https://rosalind.info/problems/smgb/) | Semiglobal Alignment | `semiglobal_alignment_linear_dna`; Java: `java/Smgb.java` |
| [oap](https://rosalind.info/problems/oap/) | Overlap Alignment | `overlap_alignment_linear_dna`; Java: `java/Oap.java` |
| [sims](https://rosalind.info/problems/sims/) | Finding a Motif with Modifications | `find_modified_motif_linear_dna` |
| [ksim](https://rosalind.info/problems/ksim/) | Finding All Similar Motifs | Java only: `java/Ksim.java` |
| [mult](https://rosalind.info/problems/mult/) | Multiple Alignment | `generalized_needleman_wunsch` (n-dimensional DP) |
| [mgap](https://rosalind.info/problems/mgap/) | Maximizing the Gap Symbols of an Optimal Alignment | `max_gap` |
| [osym](https://rosalind.info/problems/osym/) | Isolating Symbols in Alignments | `align_symbols_cost` |

### Genome assembly — `python/assembly.py`

| Problem | Title | Function(s) |
| --- | --- | --- |
| [grph](https://rosalind.info/problems/grph/) | Overlap Graphs | `overlap_graph` |
| [long](https://rosalind.info/problems/long/) | Genome Assembly as Shortest Superstring | `scs` (overlap graph + topological sort); exact subset-DP version `shortest_common_superstring` |
| [corr](https://rosalind.info/problems/corr/) | Error Correction in Reads | `correct_reads` |
| [dbru](https://rosalind.info/problems/dbru/) | Constructing a De Bruijn Graph | `de_bruijn` |
| [pcov](https://rosalind.info/problems/pcov/) | Genome Assembly with Perfect Coverage | `construct_perfect_cycle` |
| [gasm](https://rosalind.info/problems/gasm/) | Genome Assembly Using Reads | `construct_from_reads` |
| [grep](https://rosalind.info/problems/grep/) | Genome Assembly with Perfect Coverage and Repeats | `output_seqs` (`kmers_repeat`, `construct_all_seqs`) |
| [asmq](https://rosalind.info/problems/asmq/) | Assessing Assembly Quality with N50 and N75 | `nstat`, `output_nstats` |

### Phylogeny and sets — `python/phylogeny.py`

| Problem | Title | Function(s) |
| --- | --- | --- |
| [pdst](https://rosalind.info/problems/pdst/) | Creating a Distance Matrix | `p_dist_matrix` |
| [nwck](https://rosalind.info/problems/nwck/) | Distances in Trees | `output_dist` (Biopython `Phylo`) |
| [nkew](https://rosalind.info/problems/nkew/) | Newick Format with Edge Weights | `output_dist` (same function as `nwck`) |
| [tree](https://rosalind.info/problems/tree/) | Completing a Tree | `edges_to_add` (connected components − 1) |
| [inod](https://rosalind.info/problems/inod/) | Counting Phylogenetic Ancestors | `internal_nodes` |
| [seto](https://rosalind.info/problems/seto/) | Introduction to Set Operations | `set_operation`, `output_set_op` |

### Hidden Markov models (Bioinformatics Textbook Track) — `python/hmm.py`

| Problem | Title | Function(s) |
| --- | --- | --- |
| [ba10a](https://rosalind.info/problems/ba10a/) | Compute the Probability of a Hidden Path | `markov_prob` |
| [ba10b](https://rosalind.info/problems/ba10b/) | Compute the Probability of an Outcome Given a Hidden Path | `markov_prob_2` |
| [ba10c](https://rosalind.info/problems/ba10c/) | Implement the Viterbi Algorithm | `viterbi` |
| [ba10d](https://rosalind.info/problems/ba10d/) | Compute the Probability of a String Emitted by an HMM | `hmm_likelihood` (forward algorithm) |
| [ba10e](https://rosalind.info/problems/ba10e/) | Construct a Profile HMM | `profile_hmm_matrix` |
| [ba10f](https://rosalind.info/problems/ba10f/) | Construct a Profile HMM with Pseudocounts | `profile_hmm_matrix` (with pseudocounts) |
| [ba10g](https://rosalind.info/problems/ba10g/) | Perform a Multiple Sequence Alignment with a Profile HMM | `msa_hmm`, `ViterbiGraph` |
| [ba10h](https://rosalind.info/problems/ba10h/) | Estimate the Parameters of an HMM | `hmm_param_est` |
| [ba10i](https://rosalind.info/problems/ba10i/) | Implement Viterbi Learning | `viterbi_learning` |
| [ba10j](https://rosalind.info/problems/ba10j/) | Solve the Soft Decoding Problem | `soft_decoding` (`forward_algo`, `backward_algo`) |
| [ba10k](https://rosalind.info/problems/ba10k/) | Implement Baum-Welch Learning | `baum_welch` |

## Known issues

- The datasets listed under [Data](#data) are not in the repository, so their demo calls (`mprt`, `long`, `corr`,
  `grep`, `asmq`, `ba10a`–`ba10k`) raise `FileNotFoundError` until you add them to `data/inputs/`.

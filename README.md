# Miniature Telegram
A math and science utility for assisting with different calculations written in and for OCaml. 

## bioinformatics.ml
This portion of the project is dedicated to solving bioinformatics problems. The current functions include: 
* `type nucleotide = A | C | G | T | B` A new type definition representing a single nucleotide. A `B` represents a blank space when doing an alignment. 
* `generate_sequence str` Generate a sequence of nucleotides from a string `str`. Unknown characters will be converted to a `B`. 
* `string_of_nucleotide n` Convert a nucleotide `n` to a corresponding string. 
* `print_nucleotide n` Print a nucleotide `n` to the screen. 
* `print_sequence s` Print the sequence `s` to the screen.
* `print_sequence_list l` Print a list of sequences `l` to the screen. They will be lined up virtically so corresponding nucleotides are aligned. 
* `score_alignment (s, t)` Score the alignment of two sequences of nucleotides `s` and `t`. We will use a standard substitution cost matrix and a linear gap penalty to calculate the score.
* `show_alignment (s, t)` Print sequences `s` and `t` to the screen while highlighting matched nucleotides as well as transitions and transversions. 
*  `align_dfs (s, t)` Align the two sequences `s` and `t` such that the score is as good as possible. This is done with a depth first search of all possible alignments and can only be reasonably completed with sequences of length <10.
*  `generate_random_sequence n` Generate a random sequence of nucleotides of length n.

### Usage examples. 
This example is being used in the OCaml toplevel.
```ocaml 
# let s = generate_random_sequence 10;; 
val s : nucleotide list = [G; C; A; A; A; T; C; G; C; A]

# let t = generate_random_sequence 10;; 
val t : nucleotide list = [C; A; G; C; A; T; G; C; C; C]

# let (s', t') = align_dfs (s, t);;
val s' : nucleotide list = [G; C; A; A; B; A; T; C; G; C; A; B]
val t' : nucleotide list = [B; C; A; G; C; A; T; B; G; C; C; C]

# score_alignment (s, t);;
- : int = -10

# score_alignment (s', t');;
- : int = -5

# show_alignment (s,t);;
G C A A A T C G C A 
. . | . | | . . | . 
C A G C A T G C C C
- : unit = ()

# show_alignment (s',t');;
G C A A - A T C G C A - 
  | | |   | |   | | .   
- C A G C A T - G C C C
- : unit = ()
```

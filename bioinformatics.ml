open Printf
open ANSITerminal

(* The alphabet that we will be using for DNA sequences (B is a blank space.) *)
type nucleotide = A | C | G | T | B

(* A function that will generate a sequence from a string of characters. *)
let generate_sequence str = 
    let str = String.lowercase_ascii str in 
    let l = List.init (String.length str) (String.get str) in 
    let convert c = match c with 
        | 'a' -> A 
        | 'c' -> C
        | 'g' -> G 
        | 't' -> T 
        | _   -> B 
    in 
    List.map convert l 

(* Convert a nucleotide to a char. *)
let string_of_nucleotide n = match n with 
    | A -> "A"
    | C -> "C"
    | G -> "G"
    | T -> "T"
    | B -> "B"

(* Print a nucleotide with a specific color. *)
let print_nucleotide n = match n with 
    | A -> print_string [red; on_default] "A"
    | C -> print_string [blue; on_default] "C"
    | G -> print_string [green; on_default] "G"
    | T -> print_string [yellow; on_default] "T"
    | B -> print_string [black; on_default] "-"

(* A function to print a nucleotide sequince to the screen. *)
let print_sequence s = 
    let x, y = ANSITerminal.size() in 
    let rec ps n s' = 
        if n = 0 then print_string [default; on_default] "...\n"  
        else
            match s' with 
                | []   -> print_string [default; on_default] "\n"
                | c::l -> print_nucleotide c; print_string [default; on_default] " "; ps (n-1) l
    in
    ps (x-5) s
;;

(* A function to print a list of sequences in paralell. *) 
let rec print_sequence_list l = match l with 
    | []   -> ()
    | c::l -> print_sequence c; print_sequence_list l
;;

(* A function that will calculate and return the complementary DNA sequence. *)
let complement s = 
    let c n = match n with 
        | A -> T
        | T -> A
        | C -> G
        | G -> C
        | _ -> B
    in List.map c s
;;


(* A function that will score how well matched two sequences are. *)
let score_alignment (s, t) = 
    let break_cost = -2 in 
    let transition_cost = -1 in 
    let transversion_cost = -2 in 
    let alignment_score = 1 in 
    let rec sa n s' t' = match (s', t') with 
        | ([], _) | (_, [])    -> n
        | ((h1::t1), (h2::t2)) -> match (h1, h2) with 
                                    | (B, _) | (_, B) -> sa (n + break_cost) t1 t2 
                                    | (A, C) | (C, A) 
                                    | (A, T) | (T, A) 
                                    | (G, C) | (C, G) 
                                    | (T, G) | (G, T) -> sa (n + transversion_cost) t1 t2 
                                    | (A, G) | (G, A) 
                                    | (T, C) | (C, T) -> sa (n + transition_cost) t1 t2 
                                    | (_, _)          -> sa (n + alignment_score) t1 t2
    in 
    sa 0 s t
;;

(* A function to print the comparison between two sequences. *)
let show_alignment (s, t) = 
    print_sequence s;
    let x, y = ANSITerminal.size() in 
    let print_diff = 
        let min a b = if a<b then a else b in 
        let l = (min (List.length s) (List.length t)) in
        let l = min l x in 
        for i = 0 to (l-1) do 
            if (List.nth s i) = (List.nth t i) then (print_string [green; on_default] "| ")
            else (
                match ((List.nth s i), (List.nth t i)) with 
                | (A, G) | (G, A)
                | (T, C) | (C, T) -> print_string [blue; on_default] "| "
                | (_, B) | (B, _) -> print_string [default; on_default] "  "
                | (_, _)          -> print_string [blue; on_default] ". "
                )
        done
    in 
    print_string [default; on_default] "\n";
    print_diff;
    print_sequence t; 
;;

(* A function that will find the optimal alignment of two DNA sequences. *)
(* This method uses an exhaustive depth first search and is NOT efficient AT ALL. *)
(* DO NOT use with large inputs.(Like >10.) *)
let align_dfs (s, t) = 
    let max a b c = 
        if a > b then (if a>c then a else c)
        else (if b>c then b else c)
    in
    let rec mat (s', t') = match (s', t') with 
        | ([], [])             -> ([], [])
        | ([], h::t)           -> let (a, b) = (mat ([], t)) in ((B::a), (h::b))
        | (h::t, [])           -> let (a, b) = (mat ([], t)) in ((h::a), (B::b))
        | ((h1::t1), (h2::t2)) -> (
            let o1 = score_alignment ([h1], [h2]) + score_alignment (mat (t1, t2))
            and o2 = score_alignment ([h1], [B]) + score_alignment (mat (t1, t'))
            and o3 = score_alignment ([B], [h2]) + score_alignment (mat (s', t2)) in 
            let m = max o1 o2 o3 in 
            if m = o1 then let (a, b) = (mat (t1, t2)) in (h1::a, h2::b)
            else if m = o2 then let (a, b) = (mat (t1, t')) in (h1::a, B::b)
            else let (a, b) = (mat (s', t2)) in (B::a, h2::b)
        )
    in
    mat (s, t)
;;

(* A function that will generate a random DNA sequence of length n. *)
let generate_random_sequence n = 
    let rec r n' s = 
        if n' = 0 then s
        else (
            match Random.int (4) with 
            | 0 -> r (n'-1) (A::s)
            | 1 -> r (n'-1) (C::s)
            | 2 -> r (n'-1) (G::s)
            | _ -> r (n'-1) (T::s)
        )
    in 
    r n []
;;
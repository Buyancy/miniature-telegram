open Complex;; 
open Float;;
open Gnuplot;; 
module Gp = Gnuplot;;

"A mathematics library for OCaml.";;
"Made by Sean Moore <sean.moore3@mail.mcgill.ca>";;

(* Mathematical constants. *)
let e = 2.7182818284590452353602874;;
let pi = 3.14159265358979323846;;

(* A function to compute a definate integral. *)
let def_int (f : float -> float) (a, b) = 
  let segs = 10000.0 in 
  let dx = (b -. a) /. segs in
  let rec eval x n acc dx' = 
    if n <= 0.0 then acc 
    else eval (x +. dx') (n -. 1.0) (acc +. ((f x) *. dx')) dx'
  in 
  eval a segs 0.0 dx
;;

(* A function to take the derivitive of a function  f at x. *)
let derivitive (f : float -> float) x = 
  let dx = 0.0000000001 in 
  ((f (x +. dx)) -. (f x)) /. dx
;;

(* Graph a function with real (float) inputs and outputs. *)
let graph (f : float -> float) (a, b) (title: string) (x_label: string)  (y_label: string) = 
  let segs = 100.0 in 
  let dx = (b -. a) /. segs in 
  let points = 
    let rec p' n x = 
      if n<=0. then (x, f x) :: []
      else (x, f x) :: (p' (n -. 1.0) (x +. dx)) 
    in 
    p' segs a
  in
  let gp = Gp.create () in 
  Gp.plot gp ~title:title ~labels:(Gp.Labels.create ~x:x_label ~y:y_label ()) ~range:(Gp.Range.X (a, b)) (Gp.Series.lines_xy (points));
  let _ = read_line () in 
  Gp.close gp; 
;;

(* This is the Nieve implementation. We can do better. *)
let is_prime n = 
    let rec p' i n = 
        if i*i > n then true 
        else if (n mod i) = 0 then false
        else p' (i+1) n 
    in
    p' 2 n
;;

(* Simply the oposite of is prime. *)
let is_composit n = not (is_prime n);;

let square x = pow x 2.0;;

(* The factorial function for integers. *)
let fact n = 
  let rec f' n' acc = 
    if n' <= 0 then acc 
    else f' (n'-1) (n'*acc)
    in 
  f' n 1
;;

(* Stirling's approximation for factorial/gamma function. *)
let stirling x = (sqrt (2. *. pi *. x)) *. (pow (x /. e) x);; 

(* The Gamma function appriximated using stirling's approximation. *)
let gamma x = stirling x;;

(* The choose function for integers. *)
let choose n k = (fact n) / ((fact k)*(fact (n-k)));; 

(* The probability mass function for a binomial distribution. *)
let binomial_pmf n p k = 
  let q' = 1. -. p in 
  let n' = float_of_int n in 
  let k' = float_of_int k in 
  (float_of_int (choose n k)) *. (pow p k') *. (pow q' (n' -. k'))
;;
(* The probability density function for a exponential distribution. *)
let exponential_pdf x l  = l *. exp (-1.0 *. l *. x);;

(* The pdf for a gamma distribution. *)
(* a = shape paramater    b = scale parameter *)
let gamma_pdf x a b = ((Float.pow b a) /. (gamma a)) *. (pow x (a -. 1.0)) *. (exp ((-1.0) *. b *. x));;

(* The probability density function / cumulative distribution function for a normal distribution. *)
let normal_pdf x m s = (1.0 /. (s *. (sqrt (2.0 *. pi)))) *. (exp (-0.5 *. square ((x -. m) /. s)));;

let normal_cdf x m s = 
  (* The gauss error function calculated using a taylor series. *)
  let erf z = 
    let term n = 
      (Complex.div 
        (Complex.mul (Complex.pow { re = -1.0; im = 0.0 } { re = n; im = 0.0 }) (Complex.pow z { re= (2. *. n) +. 1.0; im = 0.0 }))
        {re=((float_of_int (fact (int_of_float n))) *. (2. *. n +. 1.)); im=0.0})
    in 
    let rec sum n max acc = 
      if n >= max then acc
      else sum (n +. 1.) max (acc +. (term n).re)
    in 
    (2.0 /. (sqrt pi)) *. (sum 0.0 10. 0.0)
  in 
  0.5 *. (1.0 +. erf({ re=(x -. m) /. (s *. sqrt 2.0); im=0.0 }))
;; 

(* The probability mass function for a poisson distribution. *)
let poisson_pmf k l =  (pow l (float_of_int k)) *. (exp (-1.0 *. l)) /. float_of_int(fact k);;
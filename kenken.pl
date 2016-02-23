% kenken/3
% solves an NxN KenKen puzzle with arithmetic constraints C on the NxN grid T
% Solves the 6x6 given testcase in 0.014 sec (CPU time)
% Solves a 4x4 puzzle in 0.000 sec (CPU time) using statistic/0. Evidently VERY fast.
% Memory usage:
% Memory          limit       in use     free
%
% trail  stack    16383 Kb    9 Kb       16374 Kb
% cstr   stack    16383 Kb    29 Kb      16354 Kb
% global stack    32767 Kb    7 Kb       32760 Kb
% local  stack    16383 Kb    6 Kb       16377 Kb
% atom   table    32768 atoms 1781 atoms 30987 atoms
kenken(N,C,T):-
    length(T,N), maplist(check_col_len(N), T), % make sure T is an NxN grid
    maplist(set_row_domain(N), T),
    maplist(fd_all_different, T),            % ensures each row of T contains unique values
    transpose(T, Transpose),  maplist(fd_all_different, Transpose),% ensures each column of T contains unique values
    maplist(arithmetic_constraints(T), C),   % match T with the constraints in C
    maplist(fd_labeling, T).

%%%%%%
% BEGIN HELPER FUNCTIONS
%%%%%%

% reverse syntax to fit maplist syntax
check_col_len(Len, Something) :- length(Something, Len).

% sets fd_domain for each row in T with maplist
set_row_domain(Max, Row) :- fd_domain(Row, 1, Max).

% grabs the Ith row of T and checks its Jth element against Val
matrix_element(I-J, T, Val) :- nth(I, T, Row),
			       nth(J, Row, Val). 

% transposes a matrix. pretty straightforward
% very heavily influenced from the SWI-prolog clpfd module on GitHub
transpose([], []).
transpose([F|Fs], Ts) :-
    transpose(F, [F|Fs], Ts).

transpose([], _, []).
transpose([_|Rs], Ms, [Ts|Tss]) :-
    lists_firsts_rests(Ms, Ts, Ms1),
    transpose(Rs, Ms1, Tss).

lists_firsts_rests([], [], []).
lists_firsts_rests([[F|Os]|Rest], [F|Fs], [Os|Oss]) :-
    lists_firsts_rests(Rest, Fs, Oss).
    

%%%%%%
% END HELPER FUNCTIONS
%%%%%%


% Create arithmetic constraint wrapper
arithmetic_constraints(T, C) :- wrapper_constraint(T,C).
% Match all 4 cases for elements of C
% Each constraint needs knowledge of the entire grid so we can check
% the elements of L (or J, K for / or -) against the grid itself.
wrapper_constraint(T, +(S, L)) :- addition(T, S, L, 0).
wrapper_constraint(T, *(P, L)) :- multiplication(T, P, L, 1).
wrapper_constraint(T, -(D, J, K)) :- subtraction(T, D, J, K).
wrapper_constraint(T, /(Q, J, K)) :- division(T, Q, J, K). 

% addition/4
% addition(Grid, sum, list of elems, accumulator).
% Takes an NxN grid, the sum of the integers in the list, a list of
% squares to fill, and an accumulator. Fills in the squares of the grid
% recursively according to the goal sum and the accumulator.
addition(_, S, [], S).
addition(T, S, [Hd|Tl], Acc) :-
    matrix_element(Hd, T, Val),
    Rec_acc #= Acc + Val,
    addition(T, S, Tl, Rec_acc).

% multiplication/4
% multiplication(grid, product, list of elems, accumulator).
% Similar fashion of addition but with multiplication for the accumulator
multiplication(_, P, [], P).
multiplication(T, P, [Hd|Tl], Acc) :-
    matrix_element(Hd, T, Val),
    Rec_acc #= Acc * Val,
    multiplication(T, P, Tl, Rec_acc).

% subtraction/4
% subtraction(grid, difference, J, K) where J,K  === i-j
% Subtraction and division in KenKen only deals with two squares,
% J and K. We only need to check against the grid for those two squares
% and make sure we can make the difference with J - K or K - J
% Doesn't need an explicit accumulator but uses one in the calculation
subtraction(_, D, _, _, D). %Accumulator base fact
% J - K rule
subtraction(T, D, J, K) :-
    matrix_element(J, T, Val_1),
    matrix_element(K, T, Val_2),
    Rec_acc #= Val_1 - Val_2,
    subtraction(T, D, J, K, Rec_acc).
% K - J rule
subtraction(T, D, J, K) :-
    matrix_element(J, T, Val_1),
    matrix_element(K, T, Val_2),
    Rec_acc #= Val_2 - Val_1,
    subtraction(T, D, J, K, Rec_acc).
    
% division/4
% Same logic and reasoning as subtraction/4 but uses division instead of subtraction
% Also need to make sure the division has remainder 0.
division(_, Q, _, _, Q). %Accumulator base fact
% J - K rule
division(T, Q, J, K) :-
    matrix_element(J, T, Val_1),
    matrix_element(K, T, Val_2),
    Rec_acc #= Val_1 / Val_2,
    division(T, Q, J, K, Rec_acc).
% K - J rule
division(T, Q, J, K) :-
    matrix_element(J, T, Val_1),
    matrix_element(K, T, Val_2),
    Rec_acc #= Val_2 / Val_1,
    division(T, Q, J, K, Rec_acc).
    
    





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN PLAIN_KENKEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plain_kenken/3
% exactly like kenken/3 but does not use the GNU Prolog finite domain solver
% since we don't have the finite domain to work with, we have to create our
% own range of permutations for the rows and columns of T to exist in. This
% creates an O(n!) runtime that makes it much much worse than kenken/3.
% Takes 1.045 seconds of CPU time to solve a 4x4.
plain_kenken(N,C,T) :-
    length(T, N), maplist(check_col_len(N), T),
    set_domain_plain(N, L), maplist(permutation(L), T),
    transpose(T, Transpose), maplist(plain_all_different, Transpose),
    maplist(arithmetic_constraints_plain(T), C).

%%%%%%%%%%%%%%%%%%%
% BEGIN PLAIN HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%
% creates a list L containing the range of 1 to N
set_domain_plain(N,L) :-
    findall(Num, between(1,N,Num), L).

% O(n^2) naive approach to checking if a list has unique elements.
plain_all_different([]).
plain_all_different([Hd|Tl]) :-
    \+(member(Hd, Tl)),
    plain_all_different(Tl).

%%%%%%%%%%%%%%%%%
% END PLAIN HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%


% Winter 16 Piazza post 174

arithmetic_constraints_plain(T, C) :- wrapper_constraint_plain(T, C).

wrapper_constraint_plain(T, +(S, L)) :- addition_plain(T, S, L, 0).
wrapper_constraint_plain(T, *(P, L)) :- multiplication_plain(T, P, L, 1).
wrapper_constraint_plain(T, -(D, J, K)) :- subtraction_plain(T, D, J, K).
wrapper_constraint_plain(T, /(Q, J, K)) :- division_plain(T, Q, J, K). 

addition_plain(_, S, [], S).
addition_plain(T, S, [Hd|Tl], Acc) :-
    matrix_element(Hd, T, Val),
    Rec_acc is Acc + Val,
    addition_plain(T, S, Tl, Rec_acc).

multiplication_plain(_, P, [], P).
multiplication_plain(T, P, [Hd|Tl], Acc) :-
    matrix_element(Hd, T, Val),
    Rec_acc is Acc * Val,
    multiplication_plain(T, P, Tl, Rec_acc).

subtraction_plain(_, D, _, _, D). %Accumulator base fact
% J - K rule
subtraction_plain(T, D, J, K) :-
    matrix_element(J, T, Val_1),
    matrix_element(K, T, Val_2),
    Rec_acc is Val_1 - Val_2,
    subtraction_plain(T, D, J, K, Rec_acc).
% K - J rule
subtraction_plain(T, D, J, K) :-
    matrix_element(J, T, Val_1),
    matrix_element(K, T, Val_2),
    Rec_acc is Val_2 - Val_1,
    subtraction_plain(T, D, J, K, Rec_acc).
    
division_plain(_, Q, _, _, Q). %Accumulator base fact
% J - K rule
division_plain(T, Q, J, K) :-
    matrix_element(J, T, Val_1),
    matrix_element(K, T, Val_2),
    Val_1 rem Val_2 =:= 0,
    Rec_acc is Val_1 // Val_2,
    division_plain(T, Q, J, K, Rec_acc).
% K - J rule
division_plain(T, Q, J, K) :-
    matrix_element(J, T, Val_1),
    matrix_element(K, T, Val_2),
    Val_2 rem Val_1 =:= 0,
    Rec_acc is Val_2 // Val_1,
    division_plain(T, Q, J, K, Rec_acc).
    



%%%%%%%%%%%%%%%%%%%%
% END PROLOG CODE
%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% no-op KenKen proposed API
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/*
* let the predicate be
* no_op_kenken/4
* (N, C, OperationList, T)
* Identical to kenken/3 but it will attempt to unify the operation list along with the solution T
* In this case, elements of C do not contain the operator but still contain the caged cells and 
* expected outcome.
* 
* I'm going to use a naive high level approach here.
* Instead of trying to match each arithmetic constraint, it will go through each NxN permutation of
* T that at least satifies the unique row and column constraints. Then for each permutation, it can try
* each possible arithmetic predicate and attempt to match one. Since each constraint is disjoint from
* the rest of the matrix elements, we can safely run this.
* 
* Here is an example call:
* no_op_kenken(3,
*     [(5, [1-1, 1-2, 2-2]),
*      (1, [3-3, 2-1]),
*      (8, [1-2, 1-3, 2-3]),
*      (1, [1-1, 2-2, 3-3])],
*     O,
*     T).
*  ____________________
*  
* O = [+,-,+,*]
* T = [[1,3,2],
*      [2,1,3],
*      [3,2,1]]
* 
* no
*/

# Amino-Acid-Project-Public
The public version of my ongoing efforts to investigate knowledge representation and information encoding in amino acids.

Generalized_framework_V2.py contains the code for testing individual within-codon Interpretation Frameworks.

Generalized_framework_V3.py contains the code for testing all (assuming a few restrictions) within-codon Interpretation Frameworks at once.

[Song List](https://docs.google.com/document/d/1H2WqFh46SxRQdmDHeAnxCJIUzvgVl9-fe-KuK2qX6wM/edit?usp=sharing) has the music I listened to on repeat while working on this project.

### Introduction ###
[This](https://docs.google.com/document/d/1fJm4VFJXzaQFAEqZuBsYYX4VEwT2wIqQ7nArpJahXto/edit?usp=sharing) paper contains a more in-depth look at the thought process. Unfortunately, it is nowhere near publishable and no doubt contains numerous errors as well as out-of-date information. Most notably, I've dropped using Carbon as a fundamental unit of organization (though the general logic still applies). 

Finally, my background is in Economics -- not Biology, Computer Science, or Chemistry -- so please let me know if you see problems with my approach or have clarifying comments or questions. You can email me [here](Michael.ar.campbell@gmail.com) or on [LinkedIn](https://www.linkedin.com/in/michael-campbell-73159b5a/).

### Problem Statement (non-Technical) ###

We're looking to see if there is a specific *type* of instruction set that nature uses when building an amino acid. If you're reading a set of instructions and come upon a passage that says, "A: put down three blocks in order 1, 2, 3", that's simple enough. However, what if the instruction set read, "A, if the next instruction is G: put down three blocks in order..." and right below that it read, "A, if the next instruction is C: put down *four* blocks in order..."?

If we solve the math problem stated above (using real world values for amino acids), we can determine if the instruction set cares what the letter (or letters) after the first is. This lets us know the "Interpretation Framework" (how information is stored) that nature uses to build amino acids.

### Problem Statement (Technical) ###

Let X<sub>i</sub> be a three dimensional vector with positions x<sub>1</sub>, x<sub>2</sub>, x<sub>3</sub>. Each x<sub>i</sub> may take one of four possible values: [A, C, U, G]. Let A<sub>t</sub> be an undirected heterogeneous hypergraph of finite size at time t. Find the set of graph updates rules R that takes A<sub>0</sub> -> B<sub>i</sub> for every possible X<sub>i</sub> in precisely three time steps, where B<sub>i</sub> is an undirected heterogeneous hypergraph of finite size and known topology. Assume A<sub>0</sub> is empty, such that A<sub>0</sub> = [ ].

X<sub>i</sub> - codon associated with an amino acid

[A, C, U, G] - the possible base pairs

A<sub>0</sub> - the presumed initial conditions (I, perhaps incorrectly, assume all acids start construction from the same empty graph)

B<sub>i</sub> - the final graph structure of the amino acid

#### Key and Potentially Problematic Assumptions ####
* I use physiological pH and the associated amino acid arrangements therein
* Currently, I restrict Interpretation Frameworks to using a single codon (cannot reference existence of neighbors)
* No backwards time steps. Effectively, only atoms present in the final graph state are considered for the analysis. This means that if time steps 1 and 2 use atoms that are not present at the end of 3, we would *not* pick up on that. Consider, for example, Aspartic Acid uses Xenon (for some reason), even though it doesn't appear in the final graph structure.



### Runtime ###

Code is runable on a personal computer (though V3 may take up to 3 hours -- I promise to work on further optimization, eventually...*) 
Most IFs seem to take 10k - 15k runs and ~10 minutes (though I haven't counted) with a i7-6700K running at the base 4GHz.

"*" - Truthfully, I'll only get to this when there are less important things to work on and knowing how life works -- that might be never. Given the shape of the search space, it's far more important to test assumptions (no environmental transformations, no backwards time steps, etc.) than optimize something where runtime isn't really an issue. That said, I feel like the way I handle DFs is wasteful and rewrite that *would* be a good way to improve. Fingers crossed that I find the answer and never get to work on that :P 

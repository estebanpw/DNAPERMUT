# DNAPERMUT

DNAPERMUT is a method to compute all partially-ordered permutations of `n` sequences with a reduced complexity factor. Instead of computing the permutations in `O(n!)` it will take `O(k1! + k2! + ... + kj!)` where each `ki` is an individual group of sequences, e.g. if three sequences are related to only each other but not the rest (via substring operator) then we can compute them separately and then multiply together the total number of permutations.

## Install

First, clone this repository by issuing `git clone https://github.com/estebanpw/DNAPERMUT`. Then, in the same folder, proceed to install `seqAn3` as follows.

DNAPERMUT requires [seqAn3](https://github.com/seqan/seqan3). Make sure to install it in the same folder as the cloned repository of `DNAPERMUT`. You can find installation instructions [here](https://docs.seqan.de/seqan/3-master-user/setup.html). 

Once `seqAn3`is downloaded, your directory structure should look like this:

```
├── CMakeLists.txt
├── dnapermut.cpp
├── README.md
└── seqan3
    ├── build_system
    ├── CHANGELOG.md
    ...
    └── test
```

Now simply run the following command to compile DNAPERMUT:

`cmake -DCMAKE_BUILD_TYPE=Release .`

If all dependencies were met, proceed to compile `DNAPERMUT` by firing:

`make`

And thats it!

## Use

The method works for any number of sequences, however, keep in mind that highly tangled sequences (i.e. if there are big groups of substrings) will require O(k!) computation.

`./dnapermut <fasta>`

This will generate an output such as:

```
Number of substrings:              32
Longest substring chain:           5
Total permutations (mod 100003):   88285
```

Where the `Number of substrings` corresponds to how many input strings (sequences) from the FASTA dataset are substring from another sequence in the very same FASTA. The `Longest substring chain` is the biggest group of sequences that are related by substring dependencies.

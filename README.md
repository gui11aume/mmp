## MEM Mapper Prototype: mapping short reads, faithfully ##
---
## Contents: ##
    1. Prototype? What do you mean?
    3. Compilation and installation.
    4. Running mmp

---
## I. Prototype? What do you mean? ##

`mmp` is a proper mapper, we use it in the lab and you are more than
welcome to use it as well. When we say ‘prototype’, we mean that
our purpose was not really to secure a stable user base, but rather to
show that faithful mapping is achievable.

Below are the main things to consider before using `mmp`:

    1. It works only with single-end reads.
    2. The whole read is aligned, there is no trimming.
    3. There is no way to change the sensitivity.

For now, `mmp` is single core, but we'll fix that soon.

Otherwise, `mmp` is pretty good. The most important feature is that it
is approximately faithful, so you should be able to trust the mapping
quality. You will notice that some reads have a mapping quality over 150.
It is not a bug, we mean it: the location of those reads is extremely
likely to be correct.

Finally, `mmp` is quite fast but it uses a little more memory than
many mainstream mappers.

You can find more information about `mmp` in the
[bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2020.02.10.942599v1)
that describes the mapping strategy.


II. Compilation and installation
---------------------------------

To install `mmp`, open a terminal and clone this git repository:

    git clone https://github.com/gui11aume/mmp

The files should be downloaded in a directory named `mmp`. Go to this new
directory and use `make` to compile:

    cd mmp
    make

This will create a binary file `mmp`.

III. Running `mmp`
--------------------

`mmp` runs on Linux.

### Indexing:

You first have to index a genome file before you can map reads in it. If
the genome of interest is saved in a file called `genome.fasta`, you
index it with the following command:

    ./mmp --index genome.fasta
  
This creates a bunch of files that all start with `genome.fasta.`. This
is the index.
      
### Mapping:

Once the genome is indexed, you can start mapping reads in it. For this
you must know the approximate error rate of the sequencer. By default,
`mmp` assumes that it is 1%.

If your reads are in a file called `reads.fastq`, you map them with the
following command:

    ./mmp genome.fasta reads.fastq

This produces a sam file printed on `stdout`. You can change the error
rate with the option `-e`. For instance, if you know that it is 2%,
then the command should be:

    ./mmp -e 0.02 genome.fasta reads.fastq


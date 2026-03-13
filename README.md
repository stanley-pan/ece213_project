## ECE213 Project: CUDA FastSP

After you clone, please unzip and untar the data files.

Please look at the run-commands.sh to choose what files to run
```
./build/FastSP -r data/test.true -e data/test.estimated
./build/FastSP -r data/10000/true_aligned.fa -e data/10000/twilight.aln
./build/FastSP -r data/20000/true_aligned.fa -e data/20000/twilight.aln
./build/FastSP -r data/RNASim_100K/true_aligned.fa -e data/RNASim_100K/twilight.aln
```

To run with gpu, add `-cuda` to the end of the command e.g.
`./build/FastSP -r data/test.true -e data/test.estimated -cuda`

then run with docker image -- example:
`/opt/launch-sh/bin/launch.sh -v a30 -c 8 -g 1 -m 8 -i yatisht/ece213-wi26:latest -f ./ece213_project/run-commands.sh`

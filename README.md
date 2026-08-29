# PVD - Pangenome VCF Decomposer

# Compile

```sh
mkdir build
cd build
cmake ..
make
cd ..
```

Then use `./bin/pvd`

> Add the `bin` directory to your *PATH* with: `export PATH=$PATH:/path/to/pvd/directory/bin`.

# Inputs

## Gzipped GFA

For now, the node information (sequence, position on the reference path) are extracted from the *S* line of a gzipped GFA.

To speed up the reading of that GFA, it might be useful to first make a subset of the full GFA with *S* lines only:

```sh
zcat pg.gfa.gz | awk '{if($1=="S"){print}}' | gzip > pg.Sonly.gfa.gz
```

If you have a GBZ, use something like:

```sh
vg convert -t 1 -Hf -Q refPathPrefix --no-translation pg.gbz | awk '{if($1=="S"){print}}' | gzip > pg.Sonly.gfa.gz

## or if you know/assume the S lines are all in one block around the beginning (they should, here)
## this is faster because it won't read all the other W/P/L lines
## also faster to parallelize the gzip compression with pigz
vg convert -t 1 -Hf -Q refPathPrefix --no-translation pg.gbz | awk '{if($1=="S"){E=1; print} else if(E){exit 0}}' | pigz -p 5 > pg.Sonly.gfa.gz
```

> *SO:* tags are necessary to inform the position of the nodes on the reference path. Hence the `-Q refPathPrefix` in the `vg convert` command above.

> `--no-translation` is to make sure that the nodes in the output GFA exactly match the GBZ, and that chopped nodes are not translated back (assuming that is what we want).

## Pangenome VCF

Could be a VCF made by `vg deconstruct`.
Eventually `vcfbub` could also help remove huge top-level snarls.

Input VCF could also be a `bcftool merge -m id` of single-sample VCF genotyped with `vg call`.
A same snarl can be present in up to two lines (because of the padding base) but as long as they are next to each other, they will be combined properly during decomposition.

# Usage

```sh
pvd input.vcf.gz pg.gfa.gz > input.dec.vcf
```


# Example for the HPRC v2.0 pangenome

```sh
## download pangenome and prepare GFA with S lines
wget https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/release2/minigraph-cactus/hprc-v2.0-mc-chm13.gbz
vg convert -t 4 -Hf -Q CHM13 --no-translation hprc-v2.0-mc-chm13.gbz | awk '{if($1=="S"){E=1; print} else if(E){exit 0}}' | gzip > hprc-v2.0-mc-chm13.Sonly.gfa.gz

## download deconstructed VCF
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.0/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.vcf.gz

bcftools view -o hprc-v2.0-mc-chm13.region.vcf.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.0/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.vcf.gz chr14:50000000-53000000

pvd hprc-v2.0-mc-chm13.region.vcf.gz hprc-v2.0-mc-chm13.Sonly.gfa.gz | bgzip > hprc-v2.0-mc-chm13.region.pvd.vcf.gz
```


# Potential future plan

- Option to input a VG/PG/XG/GBZ graph instead of GFA
- Add a biallelic output
- Recursively apply the decomposition on non-reference nested traversals.

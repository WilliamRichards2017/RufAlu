# Rufalu
**Author** - Will Richards

Alu mobile element detection plugin for  kmer based denovo variant detection software [rufus](https://github.com/jandrewrfarrell/RUFUS). 


#  Run with rufus (recommended)

 This plugin uses intermediate files generated from rufus as input, and writes out called Alu mobile elements to the final vcf file produced by rufus.  Rufalu is installed as an external dependency inside rufus, and the best way to use this tool to call alus is to run rufus directly, by follow the install and running instructions [here](https://github.com/jandrewrfarrell/RUFUS)
# Running as standalone app

This plugin is primarily ment to be run inside of rufus, although it can be run as a standalone app if a user has the appropriate files.

### Dependencies
[Cmake](https://cmake.org/download/) - c++ build management

### Install
```
cd bin
cmake ../
make
```

### Run 

To run Rufalu, you must run the aluDetect script with the following positional paramers:
```
./src/aluDetect $1={{path-to-bam-file}} $2={{path-to-contig file}} $3={{path-to-alu-list}} $4={{path-to-reference-file}} $5={{path-to-fastahack}} $6={{path-to-jellyfish}}
```


**Note:** For standalone use, fastahack and jellyfish are installed as external dependencies.  These will be generated during the  build, and can be found in bin/externals once cmke and make have been run

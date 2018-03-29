#include "fastqParse.h"
#include <string.h>
#include <list>
#include <stdio.h>
#include "libs/polyATail.h"
#include "libs/discordantAndChimericReads.h"
#include "libs/knownAlus.h"

int main(void){

  const char * bamPath = "../test_data/big-test.bam";
  
  const char * contigFilePath = "../../RUFUS/scripts/Family1.child.bam.generator.V2.overlap.hashcount.fastq";
  const char * aluFilePath = "../test_data/primate_non-LTR_Retrotransposon.fasta";
  const char * aluIndexPath =  "../test_data/primate_non-LTR_Retrotransposon.fasta.fai";
  const char * refPath = "/uufs/chpc.utah.edu/common/home/marth-ucgdstor/resources/references/human/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa";
  const char * refIndexPath = "/uufs/chpc.utah.edu/common/home/marth-ucgdstor/resources/references/human/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.fai";

  DACReads *dacReads = new DACReads(bamPath);
  KnownAlus *knownAlus = new KnownAlus(contigFilePath, aluFilePath, aluIndexPath, refPath, refIndexPath);

  return 0;
}


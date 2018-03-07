#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "fastqParse.h"
#include "discordantAndChimericReads.h"



void DACReads::getAllDiscordantReads(std::string inputFile, std::string outputFile){
  std::cout << "inside get all discordant reads" << std::endl;


  BamTools::BamReader reader;
  if (!reader.Open(inputFile)){
    std::cout << "Could not open input Bam file" << std::endl;
    return;
  }
  
  
  const BamTools::SamHeader header = reader.GetHeader();
  const BamTools::RefVector references = reader.GetReferenceData();

  BamTools::BamWriter writer;
  if(!writer.Open(outputFile, header, references)){
    std::cout << "Could not open output Bam file" << std::endl;
  }

  BamTools::BamAlignment al;
  while(reader.GetNextAlignment(al)){
    //Discordant reads
    if(al.HasTag("S")){
      writer.SaveAlignment(al);
      std::cout << "found discordant read" << std::endl;
      std::cout << al.QueryBases << std::endl;
      std::cout << "Alignment flag: " << al.AlignmentFlag << std::endl;

    }
  }

}

void DACReads::getAllChimericReads(std::string inputFile, std::string outputFile){
  std::cout << "inside get all discordant reads" << std::endl;


  BamTools::BamReader reader;
  if (!reader.Open(inputFile)){
    std::cout << "Could not open input Bam file" << std::endl;
    return;
  }


  const BamTools::SamHeader header = reader.GetHeader();
  const BamTools::RefVector references = reader.GetReferenceData();

  BamTools::BamWriter writer;
  if(!writer.Open(outputFile, header, references)){
    std::cout << "Could not open output Bam file" << std::endl;
  }

  BamTools::BamAlignment al;
  while(reader.GetNextAlignment(al)){
    //Chimeric reads
    if(al.HasTag("SA")){
      writer.SaveAlignment(al);
      std::cout << "found discordant read" << std::endl;
      std::cout << al.QueryBases << std::endl;
      std::cout << "Alignment flag: " << al.AlignmentFlag << std::endl;

    }
  }

}

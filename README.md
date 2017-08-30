Usage:

1. Dependencies

       LAST, Biopython, Samtools, UCSC Kent Utilities, ROAST  
   
2. Install CNSpipeline

       git clone https://github.com/liangpingping/CNSpipeline.git

       mv CNSpipeline /to/your/working/path/
   
3. Run programs

     Step1: Prepare reference genome for subsequent alignment.

             python del-reference.py reference_genome.fa

     Step2: Pairwise alignment and polish the alignment.

             python del-query.py query_genome.fa number_of_thread

     Step3: Multiple pairwise.

             roast X=0 E=reference-species species-guid-tree maf-source destination 

     Step4: Calculating the score of each base.

             python score-for-bases.py multiple_pairwise num_species
			 
	 Step5: Give a cutoff value to get CNSs according to the coverage of CDS in reference genome.
	 
	         python get-cns.py score-file CDS-file cutoff
			 

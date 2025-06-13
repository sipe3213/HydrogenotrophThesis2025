### MolMikro protocol

All credit to molecular microbiology course and Ian Marshall

3. Run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to check the quality of your reads.
    
    ```
    mkdir fastqc_output
    fastqc -o fastqc_output -t 2 *.fastq.gz
    ```


    
6. Modify the command below to trim your sequences appropriately. The backslashes `\` at the end of each line denote a single command that continues along multiple lines. This will take approximately 20 minutes.

For Samanthas files:

```
java -jar /usr/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
-threads 10 \
Genome-H9_S124_L001_R1_001.fastq.gz Genome-H9_S124_L001_R2_001.fastq.gz \
S9R1.paired.fastq.gz \
S9R1.unpaired.fastq.gz \
S9R2.paired.fastq.gz \
S9R2.unpaired.fastq.gz \
CROP:240 \
HEADCROP:20 \
SLIDINGWINDOW:4:20 \
MINLEN:100 \
ILLUMINACLIP:/usr/bbmap/resources/adapters.fa:2:40:15

```
EG files:

```
screen -L
java -jar /usr/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
-threads 8 \
EGP025_24_012_R1_001.fastq.gz EGP025_24_012_R2_001.fastq.gz \
EG12R1.paired.fastq.gz \
EG12R1.unpaired.fastq.gz \
EG12R2.paired.fastq.gz \
EG12R2.unpaired.fastq.gz \
ILLUMINACLIP:/usr/bbmap/resources/adapters.fa:3:30:10 MINLEN:36

```


  * `java -jar /usr/Trimmomatic-0.39/trimmomatic-0.39.jar PE` - This is the location of the trimmomatic program and a parameter `PE` indicating we're working with paired end reads.
  * `-threads 2` - This specifies the number of processors that trimmomatic can use. We have 12 processors on our server, so using two CPUs will allow all groups to work simultaneously.
  * `Course-X_SXXX_L001_R1_001.fastq.gz` 
  `Course-X_SXXX_L001_R2_001.fastq.gz` - This is your input data, the filenames of the files you have been given from MiSeq sequencing. Rename these to match your own data.
  * `Course-X_SXXX_L001_R1_001.paired.fastq.gz` 
  `Course-X_SXXX_L001_R1_001.unpaired.fastq.gz` 
  `Course-X_SXXX_L001_R2_001.paired.fastq.gz` 
  `Course-X_SXXX_L001_R2_001.unpaired.fastq.gz` - These are your output filenames. The `paired` output is where both paired reads have been retained following trimming and filtering, the `unpaired` is where one of the reads has been removed. Rename these to match your own filenames.
  * `CROP:240` - This crops (cuts) the sequence to a certain maximum length. You determine this number in step 5.
  * `HEADCROP:20` - This crops the sequence at its "HEAD", or at the start. You determine this number in step 5.
  * `SLIDINGWINDOW:4:20` - This parameter trims based on quality. A "sliding window" finds the average sequence quality across a given window size (in this case 4) and cuts off the sequence when the mean quality drops below a certain value (in this case 20). You shouldn't need to change this value.
  * `MINLEN:100` - This specifies a minimum length for sequences following trimming (in this case 100 bases). You shouldn't need to change this value.
  * `ILLUMINACLIP:/usr/bbmap/resources/adapters.fa:2:40:15` - This trims off certain specific sequences (Illumina adapters) that may be left over following other trimming. You shouldn't need to change this value.


7. Re-run FastQC (step 3) and look at your trimmed sequences, making sure to make a new output folder with a different name than the first time around. If you are unsatisfied with your trimming (i.e. the read quality is still too low, or you suspect that trimming was too strict and too many reads were removed) then try tweaking the parameters in step 6 and repeating the process.

    ```
    mkdir fastqc_output_post_trimming
    fastqc -o fastqc_output_post_trimming -t 2 *.paired.fastq.gz
    ```
    ```
    mkdir fastqc_EGoutput_post_trimming
    fastqc -o fastqc_EGoutput_post_trimming -t 2 *.paired.fastq.gz
    ```

    8. For technical reasons, the genome of the PhiX phage is often used as a spike-in during Illumina sequencing. This is usually removed due to different index sequences, but errors in reading the index sequences can occasionally allow some PhiX sequences to leak into your genome libraries. To make sure this isn't the case we will scan PhiX sequences and remove them using the following command:

    ```
    bbduk.sh \
    in=S10R1.paired.fastq.gz \
    in2=S10R2.paired.fastq.gz \
    out=S10R1.paired.nophix.fastq.gz \
    out2=S10R2.paired.nophix.fastq.gz \
    ref=/usr/bbmap/resources/phix174_ill.ref.fa.gz \
    k=31 \
    threads=10 \
    hdist=1
    ```
    ```
    bbduk.sh \
    in=EG11R1.paired.fastq.gz \
    in2=EG11R2.paired.fastq.gz \
    out=EG11R1.paired.nophix.fastq.gz \
    out2=EG11R2.paired.nophix.fastq.gz \
    ref=/usr/bbmap/resources/phix174_ill.ref.fa.gz \
    k=31 \
    threads=10 \
    hdist=1
    ```

    * `in=Course-X_SXXX_L001_R1_001.paired.fastq.gz` 
    `in2=Course-X_SXXX_L001_R2_001.paired.fastq.gz` - This indicates the name of your input fastq files containing forward reads (`in`) and reverse reads (`in2`). Trimmed files are used here.
    * `out=Course-X_SXXX_L001_R1_001.paired.nophix.fastq.gz` `out2=Course-X_SXXX_L001_R2_001.paired.nophix.fastq.gz` - Indicates the names of your output fastq files.
    * `ref=/usr/bbmap/resources/phix174_ill.ref.fa.gz` - A fasta file containing the reference PhiX genome that your reads will be compared against.
    * `k=31` - Specifying kmer length as 31. The longer this is the slower it will run but the more accurate it will be, 31 is recommended for this operation by the author of BBTools.
    * `threads = 2` - This specifies the number of processors that bbduk can use. We have 12 processors on our server, so using two CPUs will allow all groups to work simultaneously.
    * `hdist = 1` - hdist stands for "hamming distance" - the number of allowable errors per kmer

    Look at the output from `bbduk.sh` - you probably won't have phiX contamination - in this case, proceed to step 9 using your non-decontaminated trimmed reads as inputs. If you do have contamination, you will need to use the `nophix` decontaminated output that `bbduk.sh` has produced.

9. Assemble your genome using [SPAdes](http://cab.spbu.ru/software/spades/). This will take about two hours. If you won't have time to keep your laptop connected for a the full time period then be sure to prepend the command below with `screen -L`. 

    cd ~/genome/samantha_genomes/S_trimmed_fastq
    ```
    screen -L spades.py \
    --pe1-1 S10R1.paired.fastq.gz \
    --pe1-2 S10R2.paired.fastq.gz \
    -o spades_output_S10 \
    --careful \
    -k 21,33,55,77,99,127 \
    --cov-cutoff auto \
    --threads 8
    ```
    cd ~/genome/samantha_genomes/EGtrimmedfastq/Nophix
    ```
    screen -L spades.py \
    --pe1-1 EG11R1.paired.nophix.fastq.gz \
    --pe1-2 EG11R2.paired.nophix.fastq.gz \
    -o spades_output_EG11 \
    --careful \
    -k 21,33,55,77,99,127 \
    --cov-cutoff auto \
    --threads 10
    ```
  * `spades.py --pe1-1 Course-X_SXXX_L001_R1_001.paired.fastq.gz --pe1-2 Course-X_SXXX_L001_R2_001.paired.fastq.gz` - This specifies the input file name for the paired end forward read file (`--pe1-1`) and reverse read file (`--pe1-2`). Change this to match your filenames.
  * `-o spades_output_directory` - This specifies the name of the directory where the assembly will be stored. You should change this if this isn't the first time you run `spades`.
  * `--careful` - Careful mode tries to reduce number of mismatches and short indels.
  * `-k 21,33,55,77,99,127` - k-mer sizes - you shouldn't have to change this value.
  * `--cov-cutoff auto` - automatically determines a minimum necessary coverage.
  * `--threads 2` - This specifies the number of processors that SPAdes can use. You shouldn't have to change this value.


  3. To obtain basic statistics about the assembly, run Quast.
    
    ```
    cd genome/samantha_genomes/Scaffolds
    quast.py EG12_high_GC_decontaminated_genome_assembly.fasta
    ```
    ```
    cd genome/samantha_genomes/Scaffolds
    quast.py EG11_scaffolds.fasta
    ```

4. Take a look at a text report of the Quast results using `cat`, a small program that prints the contents of any file quickly to the command line. You will now see some basic statistics about your genome assembly, including the total number of contigs, total length, and N50. Keep in mind that although Quast is listing "contigs", these are in fact "scaffolds" produced by SPAdes. The `quast_results/latest` directory also contains the report in several other formats, including pdf, that you can download and view on your own computer.

    ```
    cd ~/genome/samantha_genomes/Scaffolds

    cat EG12_quast_report.txt
    ```
    cat S3_quast_report.txt
    ```
Reflect on how your genome statistics (especially total length and GC content) compare with closely related genomes you found in NCBI or IMG. Download the `report.txt` file to your own computer for later use in the Genome Report.

5. Calculate the mean coverage of your genome. To do this, divide the total number of bases used to make your assembly by the `Total length (>= 0 bp)` from your Quast output. Find the total number of bases in your fastq input we will use a small program called `fastq-stats` that's part of the [**ea-utils**](http://expressionanalysis.github.io/ea-utils/) package. You will need the coverage for your Genome Report, so be sure to make a note of these numbers.

    ```
    cd ~/genome/samantha_genomes/S_trimmed_fastq
    fastq-stats S1R1.paired.fastq.gz
    fastq-stats S1R2.paired.fastq.gz
    ```
    The final line of the output called `total bases` is what we're interested in here. Add the number of bases together for your forward and reverse paired reads for the total number of bases you used in your assembly. If you also used unpaired reads you will need to run `fastq-stats` on those too. Now take that number and divide it by the total assembly length to get the coverage.

6.  Return to your home directory (`~`) and then run **CheckM2** to estimate the completeness and degree of contamination of your genome. Note that CheckM2 is designed to operate on all files ending with a given file extension (in this case `.fasta`) within a given folder, so you don't need to specify the input file directly, just the directory containing the input file(s). This will take about three minutes to run.

cd genome/samantha_genomes
    ```
    checkm2 predict \
    --input cekm \
    --extension fasta \
    --output-directory new_checkm2_output \
    --threads 10
    ```

    * `checkm2 predict` - The command that calls the program CheckM2 and tells it to `predict` the completeness and contamination of a genome.
    * `--input genome`- This specifies the input directory, in this case called `genome`. This directory can contain one or more fasta files including all the scaffolds for a genome.
    * `--extension fasta`- The input directory may contain other files that are not genome fasta files. How can CheckM2 tell the difference? By the file extension - in this case, any file ending in `.fasta` will be considered a genome.
    * `--output-directory checkm2_output`- This is the directory where checkm2 will store its output. Note that checkm2 will not work if this directory already exists.
    * `--threads 2` - This parameter specifies the number of threads that the program will operate across - the more threads, the faster it should run.

7. Look at your **CheckM** results using `cat`. Note the "Completeness" and "Contamination" percentages. Make a note of these numbers, or save `quality_report.tsv` on your own computer, for use in your Genome Report. Note that you might need to open this tsv (tab-separated value) file in Excel to get the headers to line up with each column correctly.

    ```
    cat new_checkm2_output/quality_report.tsv
    ```

8. Determine the coverage for each scaffold in your genome assembly using **bbmap**. You will need to tell bbmap where your reads files are (your `.fastq.gz` files) and where your genome assembly `.fasta` file is. This will take about 30 minutes.

    ```
    cd genome/samantha_genomes

    bbmap.sh \
    ref=S9_scaffolds.fasta \
    in=S9R1.paired.fastq.gz \
    in2=S9R2.paired.fastq.gz \
    covstats=S9_bbmap_coverage_statistics.tsv \
    threads=10 \
    -Xmx6g nodisk
    ```
    cd genome/samantha_genomes

    bbmap.sh \
    ref=EG12_scaffolds.fasta \
    in=EG12R1.paired.nophix.fastq.gz \
    in2=EG12R2.paired.nophix.fastq.gz \
    covstats=EG12_bbmap_coverage_statistics.tsv \
    threads=10 \
    -Xmx6g nodisk
    ```
    * `bbmap.sh` - This calls the program **bbmap**.
    * `ref=meaningfulname_scaffolds.fasta` - This is the name of your genome assembly file (your scaffolds).
    * `in=Course-X_SXXX... in2=Course-X_SXXX...` - `in` and `in2` specify the filenames for the forward and reverse read `fastq` files used to make your genome assembly
    * `covstats=bbmap_coverage_statistics.tsv` - This specifies the name of your output file. You may want to change this to help keep track of your files.
    * `threads=2` - This specifies the number of threads you will use (2, in this case).
    * `-Xmx12g nodisk` - This specifies the maximum amount of memory that bbmap is permitted to use (12 gigabytes, the available memory for each of five groups when dividing 64GB of RAM between them) and `nodisk` specifies that the reference database should not be written to the disk.
    *Note:* bbmap is weird in that there are no spaces in between flags and most parameters (i.e. it's `threads=2`, not `threads= 2`). So don't put spaces where they shouldn't be.


    18. Re-run **Quast** and **CheckM** on the decontaminated genome assembly, changing the input filename for `quast` accordingly. Has the decontamination procedure changed the result? Have you reduced the amount of contamination detected by CheckM? Have you removed too much of your genome? You can repeat this procedure and try a different lassoing pattern in step 15 if you are unhappy with the result.


    19. Find taxonomy of potential decontaminated fasta files.


First make copy of scaffold fasta files or rename them to fasta after.
Fasta to fna

cd ~/genome/samantha_genomes/Genometotax
"
for file in *.fasta; do mv "$file" "${file%.fasta}.fna"; done
"
"
gtdbtk classify_wf --genome_dir samantha_genomes/Genometotax --out_dir classify_wf_out --cpus 10 --skip_ani_screen
"
or
"
gtdbtk classify_wf --genome_dir ~/genome/samantha_genomes/Genometotax --out_dir ~/genome/samantha_genomes/Genometotax/deviants_classify_wf_out --cpus 10 --mash_db ~/genome/samantha_genomes/Genometotax/Mashlocation
"
Then check the summary
"
cat classify_wf_out/gtdbtk.bac120.summary.tsv
"


20. For pairwise genome alignment, you can use fastANI to compaire two genomes:

cd ~/genome/samantha_genomes/Scaffolds

"
fastANI -q EG3_scaffolds.fasta -r EG3_scaffolds.fasta -o eg3vseg3.txt
"



7. Once you have uploaded annotated genome sequences for approximately 5 relative genomes, unzip these on the server and then concatenate protein sequences for all genomes into a single file using the `cat` command:

    ```
    cd relative_genomes
    gunzip *genomic.gbff.gz
    cat *genomic.gbff > relative_genomes.gb
    ```

8. You are now ready to annotate your genome. Use the `cd` command to change to the directory where your assembled contigs fasta file is and then use `prokka` to annotate your genome. This will take about ten minutes to run.

    ```
    cd ../genome

    prokka \
    --outdir prokka_annotation \
    --prefix output_filename_prefix \
    --kingdom Bacteria \
    --proteins ~/relative_genomes/relative_genomes.gb \
    --rawproduct \
    --mincontiglen 500 \
    --cpus 2 \
    decontaminated_genome_assembly.fasta
    ```
    cd ~/genome/samantha_genomes/Scaffolds

    prokka --outdir S1_annotation --prefix S1_Desulfolegal --genus Desulfovibrio --species legallii --kingdom Bacteria --cpus 10 --locustag DESULF --mincontiglen 500 S1_scaffolds.fasta




   * `--outdir prokka_annotation` - Prokka will make a new folder to dump its output into - in this case `prokka_annotation`, but you can change this if you like. Note that if you run `prokka` for a second time you will have to change this - `prokka` will (wisely) not overwrite a previous annotation.
   * `--prefix output_filename_prefix` - Here it's `output_filename_prefix`, but this is the string of letters you want to be on the beginning of every output file (i.e. output_filename_prefix.gbk, output_filename_prefix.fna etc.), so it should be the name of your organism - a logical name for your organism would be the genus name then, `MM2023` for Molecular Microbiology 2023, followed by your group number, i.e. `Bacillus_MM2023_1` for group 1.
   * `--kingdom Bacteria` - This specifies the kingdom of life your genome belongs to - this should be changed if you're annotating a viral, archaeal, or mitochondrial genome.
   * `--proteins ~/relative_genomes/relative_proteins.faa` - This specifies your protein database to annotate from - Prokka will first run BLAST against this database for each protein before moving on to its larger databases. This is the file you prepared in steps 2--7.
   * `--rawproduct` - By default, prokka "cleans up" your protein annotations to make it easier to submit them to NCBI/ENA. This cleaning process sometimes removes useful information from the annotation, so here we are turning this feature off to just return the "raw product". When we submit these genomes to NCBI we will reannotate them using NCBI's annotation pipeline, so there is no need to worry about compliance with database standards at this point.
   * `--mincontiglen 500` - This sets them minimum contig (scaffold, in our case) length to annotate (500 base pairs here).
   * `--cpus 2` - This restricts the number of CPUs (threads) you will use to two.
   * `decontaminated_genome_assembly.fasta` - This is the file with your assembled scaffolds - change this to whatever you called your decontaminated scaffolds file (if you carried out any decontamination).
   You will see there are other command line parameters for Prokka in the manual or by typing `prokka -h`.

11. If you `cd` into the output directory (specified by you – `prokka_annotation` from above) you will see a number of files that make up your annotated genome:
   * `err` – errors and warnings from the annotation procedure
   * `faa` – fasta file with amino acid sequences for all proteins in genome
   * `ffn` – fasta file with nucleotide sequences for all genes in genome
   * `fna` – fasta file with all scaffolds
   * `fsa` – same as fna file, but with extra information required for NCBI submission 
   * `gbk` – full annotated genome, genbank format
   * `gff` – full annotated genome, gff format
   * `log` – prokka’s output from when it was run
   * `sqn` – file for Genbank submission
   * `tbl` – feature table file (also related to Genbank submission)
   * `txt` – summary of the genome - here you will find the total number of protein-coding sequences (CDS) and various kinds of RNA genes (rRNA, tRNA, etc.), you will need this numbers for your report

   Download the `prokka_annotation` folder to your own computer to carry out subsequent steps.

12. If you have a Sanger-sequenced 16S rRNA gene for your genome you should confirm that your genome 16S rRNA sequence is near-identical to the Sanger sequence. Open the `ffn` file on your own computer in **Sublime Text** or another text editor, and search for "16S ribosomal RNA" (ensure that it is the ribosomal RNA gene you have found, and not a sequence encoding a ribosomal protein). If you have multiple 16S rRNA genes in your `ffn` file you should check them all to make sure they all come from the same organism. Copy the description line and gene sequence to a new file, and copy your Sanger description line and sequence to the same file to make a `fasta`-formatted file with the two 16S rRNA gene sequences. Open a web browser and go to the EBI muscle alignment web tool [https://www.ebi.ac.uk/Tools/msa/muscle/](https://www.ebi.ac.uk/Tools/msa/muscle/). Copy and paste the two sequences into the input box and click *Submit*. Are the sequences nearly identical? Some differences, especially towards the 3' end of the Sanger, are an acceptable consequence of degrading sequence quality towards the end of the read and not necessarily cause for concern.



14. It is sometimes useful to try an alternative annotation tool to determine protein function. We will use *Interproscan*, one of many only annotation tools. This tool takes a long time to run, so run it with `screen` so it will keep running after you log off the server. Make sure you're in your prokka output directory, then use the following command:

    ```
    screen -L interproscan.sh -i output_filename_prefix.faa \
    -f TSV \
    -cpu 2
    ```

    * `-i output_filename_prefix.faa` - This is the fasta file containing the protein sequences inferred from your genome. Change `output_filename_prefix` to whatever you used in step 10.
    * `-f TSV` - This specifies that your output will be in tab-separated value (tsv) format, a way of making a table in a text file.
    * `-cpu 2` - This specifies the number of CPU threads to run in parallel.

    The output file will appear with the name `output_filename_prefix.tsv` - download this file to your own computer and open it in a text editor or Excel to see the matches.

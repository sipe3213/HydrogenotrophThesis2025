### Lucas thesis 2025
### Bioinformatics Protocol
# All codes in the order I would run it. 


# Getting reads fixed up

* 1. Collect nanoplot fastq files to one
    
        cat *.fastq.gz > x_reads.fastq.gz

* 2. Nanoplot

        screen -L \
        NanoPlot \
        --loglength \
        --threads 10 \
        --outdir pretrimming_nanoplot_output \
        --fastq 46_reads.fastq.gz

* 3. Filtering

        screen -L \
        filtlong \
        --target_bases 3000 \
        barcode_8_reads.fastq.gz \
        | gzip > b8_filtered.fastq.gz
        
        or
        
        filtlong \
        --min_length 4000 \
        --keep_percent 90 \
        --target_bases 250000000 \
        55_reads.fastq.gz \
        | gzip > 55_filtered.fastq.gz

* 4. Trimming

        screen -L \
        porechop \
        -i PC8_reads.fastq.gz \
        -o pc8_trimmed.fastq.gz \
        --threads 10

* Run nanoplot again on newly trimmed and filtered fastq files

# Assembly

* 1. Flye

        screen -L \
        flye \
        --nano-hq filtertrimmed/pc2_filtered_trimmed.fastq.gz \
        --genome-size 5m \
        --threads 10 \
        --out-dir pc2_flyeplish2ndtry_assembly
        

* 2. autocycler 

        screen -L autocycler subsample \
        --reads inputfile_trimmed_filtered.fastq.gz \
        --out_dir subsampled_reads \
        --genome_size Xm
        ```

        screen -L
        mkdir assemblies
        for assembler in flye raven; do
            for i in 01 02 03 04; do
                $assembler.sh \
                    subsampled_reads/sample_$i.fastq \
                    assemblies/"$assembler"_"$i" 10 Xm
            done
        done
        exit

    * Next we will compress these assemblies to simplify the assembly of assemblies computationally.

        screen -L autocycler compress -i assemblies -a autocycler_out

    * Next we will cluster the assemblies together to find the contigs that match from each assembly and put them in the same cluster.

        screen -L autocycler cluster -a autocycler_out


    * Trim and resolve clusters:

        for c in autocycler_out/clustering/qc_pass/cluster_*; do
            autocycler trim -c "$c"
            autocycler resolve -c "$c"
        done

    * Combine clusters into final assembly:

        autocycler combine -a autocycler_out \
        -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa
    


# Polish

* 1. Polypolish
    * Polish nanopore assembly with illumina short-read files.

        cd main_fastq

    * First run this

        bwa index b13_flye_assembly/assembly.fasta
        bwa mem -t 10 -a b13_flye_assembly/assembly.fasta illuminashortreadforpolish/S2R1.paired.fastq.gz > b13_alignments_1.sam
        bwa mem -t 10 -a b13_flye_assembly/assembly.fasta illuminashortreadforpolish/S2R2.paired.fastq.gz > b13_alignments_2.sam
        polypolish filter --in1 b13_alignments_1.sam --in2 b13_alignments_2.sam --out1 b13_filtered_1.sam --out2 b13_filtered_2.sam

    * Then run this

        polypolish polish b13_flye_assembly/assembly.fasta b13_filtered_1.sam b13_filtered_2.sam > b13_polished.fasta
        rm *.amb *.ann *.bwt *.pac *.sa *.sam

    * If mixed sample (here two different genomes)

        bwa index b11_meta_flye_assembly/assembly.fasta
        bwa mem -t 10 -a b11_meta_flye_assembly/assembly.fasta illuminashortreadforpolish/EG10R1.paired.nophix.fastq.gz > b11_alignments_a_1.sam
        bwa mem -t 10 -a b11_meta_flye_assembly/assembly.fasta illuminashortreadforpolish/EG10R2.paired.nophix.fastq.gz > b11_alignments_a_2.sam
        bwa mem -t 10 -a b11_meta_flye_assembly/assembly.fasta illuminashortreadforpolish/EG7R1.paired.nophix.fastq.gz > b11_alignments_b_1.sam
        bwa mem -t 10 -a b11_meta_flye_assembly/assembly.fasta illuminashortreadforpolish/EG7R2.paired.nophix.fastq.gz > b11_alignments_b_2.sam
        polypolish filter --in1 b11_alignments_a_1.sam --in2 b11_alignments_a_2.sam --out1 b11_filtered_a_1.sam --out2 b11_filtered_a_2.sam
        polypolish filter --in1 b11_alignments_b_1.sam --in2 b11_alignments_b_2.sam --out1 b11_filtered_b_1.sam --out2 b11_filtered_b_2.sam

    * then 

        polypolish polish b11_meta_flye_assembly/assembly.fasta b11_filtered_*.sam > b11_polished.fasta


* Flye polish happened already as it assembled the genome, but for more iterations i used `--iterations 10`  for none-short read polish.
        
        screen -L \
        flye \
        --nano-hq filtertrimmed/pc2_filtered_trimmed.fastq.gz \
        --genome-size 5m \
        --threads 10 \
        --out-dir pc2_flyeplish2ndtry_assembly
        --iterations 10

# Contig configering 

* 1. Split contigs into individual files:
        
        awk '/^>/ {OUT="<output_folder_name>/" substr($0,2) ".fasta"}; {print >> OUT; close(OUT)}' <input_assembly_name>.fasta

    * Rename contigs appropriately (must not contain spacing " ")
        
        cd <output_folder_name>
        for file in *.fasta ; do mv "$file" "<id_prefix>_${file}" ; done
        cd .. 
    
    *** <id_prefix> = the first part of all contig names, such as strain id, which is not to be altered.

* 2. Check contig completeness
        
        checkm taxonomy_wf \
        --threads 10 \
        --file <output_results>.txt \
        --extension fasta \
        domain Bacteria \ 
        <input_folder> \
        <output_folder>

        cat <output_results>.txt
    
# Classification (gtdbtk)

* 1. First make copy of scaffold fasta files which you will convert fasta to fna:
        cd ~/<selected_contigs_folder>
        for file in *.fasta; do mv "$file" "${file%.fasta}.fna"; done
        cd ..

* 2. Then run gtdbtk on them:
        gtdbtk classify_wf --genome_dir <selected_contigs_folder> --out_dir <output_folder> --cpus 10 --mash_db </location/preferred_mash_location>

    * Then check the summary
        cat <output_folder>/gtdbtk.bac120.summary.tsv



# Annotation

* 1. Closely-related strains for annotation reference (downloaded from NCBI genomic database) processed sequencially:
        
    * from .gbff into .faa:

        prokka-genbank_to_fasta_db --idtag=TAG --format=gbff <relative_genome_name_ncbi>.gbff > <relative_genome_name_ncbi>.faa
    
    * To get rid of the redundancy we'll use a program called **cd-hit** to make clusters of similar proteins and then pick a representative sequence.

        cdhit -i <relative_genome_name_ncbi>.faa -o <relative_genome_name_ncbi_cdhit_output> -T 6 -M 8000 -g 1 -s 0.8 -c 0.9

    * Convert new cdhit processed file into database files    
        
        makeblastdb -dbtype prot -in <relative_genome_name_ncbi_cdhit_output>

    * Import into prokkadb
        
        cp <relative_genome_name_ncbi_cdhit_output>.p* /usr/prokka/db/genus/

        prokka --listdb


* 2. Prokka Annotation NANOPORE
    
        prokka \
        --outdir </output_folder> \
        --prefix <prefix> \
        --locustag <TAG> \
        --genus <genus> --species <species> --strain <strain> \
        --usegenus \
        <input_polished_assembly>.fasta \
        --mincontiglen 500 --centre A \
        --cpus 10

* 3. Illumina 

        cat ~/allnanopore/relatives/Paenibacilluschibensis_ncbi.gbff > ~/allnanopore/relatives/Paenibacilluschibensis_ncbi.gbk

        prokka \
        --prefix S3 \
        --mincontiglen 500 \
        --kingdom Bacteria \
        --cpus 6 \
        fna/S3_scaffolds.fna

# HMMscan

* 1. First, align each HydDB group fasta sequences
        
        Usage: `muscle` -in <inputfile> -out <outputfile>

        ```    
        muscle \
        -in ~/allnanopore/hmm/hyddb/NiFe4i.fasta \
        -out ~/allnanopore/hmm/msa/NiFe4i.clw \
        -clw

        ```

* 2. Then build hhm-profile
        
        Usage: `hmmbuild` [-options] <hmmfile_out> <msafile>
        
        ```
        hmmbuild  \
        --amino  \
        --cpu 10 -n NiFe1d  \
        -o ~/allnanopore/hmm/hmmbuildout/summary/NiFe1d.txt  \
        ~/allnanopore/hmm/hmmbuildout/NiFe1d.hmm  \
        ~/allnanopore/hmm/msa/NiFe1d.clw 

        ```
    * combine into a hmm database
        
        ```
        cat *.hmm > HydDBgroups.hmm

        ```



* 3. then hmmpress
        
        Usage: hmmpress [-options] <hmmfile>
        
        hmmpress  \
        ~/allnanopore/hmm/hmmbuildout/HydDBgroups.hmm

* 4. Then hmmscan
        Usage: hmmscan [-options] <hmmdb> <seqfile>

        hmmscan \
        -E 0.00001 \
        -o ~/allnanopore/hmm/hmmresults/76HydDBgroups.txt \
        --tblout ~/allnanopore/hmm/hmmresults/76HydDBgroups_sum.txt \
        --cpu 6 \
        ~/allnanopore/hmm/hmmbuildout/HydDBgroups.hmm \
        ~/allnanopore/annotation/S9.faa




* 5. For converting hmm results to a combined excel format.

        cd ~/allnanopore/hmm/txt_to_csv

        python multi_txt_to_excel_converter.py hmmscan_output.xlsx 




# kofam scan  
Version 1.3 https://github.com/takaram/kofam_scan
Inspiration from Taylor Reiters guide https://taylorreiter.github.io/2019-05-11-kofamscan/ 

cd ~/working_directory

* 1. Downloaded the databases and executables:
        wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz		# download the ko list 
        wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz 		# download the hmm profiles
        wget ftp://ftp.genome.jp/pub/tools/kofamscan/kofamscan.tar.gz	# download kofamscan tool
        wget ftp://ftp.genome.jp/pub/tools/kofamscan/README.md		# download README

    * if they dont work, just download manually from website.

* 2. Unzipping and untarring relevant files:
        gunzip ko_list.gz
        tar xf profiles.tar.gz
        tar xf kofam_scan-1.3.0.tar.gz


* 3. Edit template of config.yml to your preference and place it with exec_annotation file. Then run it:

        exec_annotation -f mapper -o </location/genome_name> </location/genome_name>.faa

    * Combine all output txt files to one:
        cat *.txt > collectedkofamscangenomeswall.txt

* 4. Then run Keggdecoder with a python prompt (I used anaconda) to either make .svg (picuture) or .list (table format) file  
        
        conda activate keggdecoder

        KEGG-decoder -i <C:\location\exec_annotation_output>.txt -o <C:\location\exec_annotation_output>.svg -v static

        KEGG-decoder -i <C:\location\exec_annotation_output>.txt -o <C:\location\exec_annotation_output>.list -v static


# If spacing in fasta headers (The reference genomes from NCBI, for example):

    * Remove space in header of fasta:
        
        sed 's, ,_,g' -i <genome_name>.faa

    * Rename header of fasta contigs:
        
        rename.sh ignorejunk=t <in=genome_name>.faa <out=genome_name>.faa prefix=<PREFIX>


In order to blast our sequences and get the closest relatives to the test sequences we created a blast database in the following steps

1. From the NCBI we obtained all the available sequences for microbiota of Anopheles sinensis to form the refence sequences.
2. We next combined these sequences into a single FASTA file named `my_reference.fasta` . 
3. Downloaded the blast executables `wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz` and unzipped it on the command-line using the code `tar -zxvf ncbi-blast-2.13.0+-x64-linux.tar.gz`
4. On the command-line, we then used the reference file to create a nucloetide data, by calling `"makeblast"`:

    `makeblastdb -in <my_reference.fasta> -dbtype nucl -parse_seqids -out <database_name> -title "Database title"`
5. Now to blast the sequences from the three colonies, against the sequences already in the local database we run, the code;
   `blastn -query <testfile.fasta> -db <database_name> -out <blasted-seq.fasta>`

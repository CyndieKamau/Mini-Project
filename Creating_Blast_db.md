# **Creating a local BLAST database**

1. Download all the sequences to be used in the database in NCBI. Example; 'Microbiota of _A. sinensis_.'

2. Combine all sequences into a single FASTA file named eg. `my_reference.fasta` .


3. Download the blast executables 
 
 `wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz` 


4. Unzip it on the command-line using the code 


 `tar -zxvf ncbi-blast-2.13.0+-x64-linux.tar.gz`.


5. A directory called `ncbi-blast-2.13.0+` will be formed.

 
 
6. To export the PATH, use this command: 

`export PATH=$PATH:$HOME/ncbi-blast-2.13.0+/bin`



7. On the command-line,use the reference files to create a nucleotide/protein database, by calling `"makeblast"`:


`makeblastdb -in <my_reference.fasta> -dbtype nucl -parse_seqids -out <database_name> -title "Database title"`
    

8. If the run is successful output should be as follows:

```
Building a new DB, current time: 08/25/2022 20:04:17
    
New DB name:   /home/icipe/Desktop/mini-project/blasting/database/ref-db

New DB title:  Sequence database

Sequence type: Nucleotide

Keep MBits: T

Maximum file size: 1000000000B

Adding sequences from FASTA; added 1524388 sequences in 49.1813 seconds.
```

9. Blast your query sequences against the local database using this code;


`blastn -query <testfile.fasta> -db <database_name> -out <blasted-seq.fasta>`

# Kmer2SNP: 

**Kmer2SNP** is a fast and accurate variants calling tool for next-generation sequencing data. It's a reference-free and alignment free method.
 
# Downloading Kmer2SNP

To download **Kmer2SNP**, you have to install the following software
<pre><code>
(1) DSK: 
https://github.com/GATB/dsk
(2) findGSE: 
https://github.com/schneebergerlab/findGSE/blob/master/INSTALL
(3) Python package networkx
https://networkx.github.io/documentation/networkx-1.11/install.html
 </code></pre>
  
Then clone the **Kmer2SNP** repository to your machine.
<pre><code> git clone https://github.com/yanboANU/Kmer2SNP.git </code></pre>

After you download Kmer2SNP, you need to change path in the runkmercalling.sh.

# Input Format

The input of **Kmer2SNP**  can be fasta, fastq, either gzipped or not. 

# Ouput Format

See Kmer2SNP/example, k_31_pair.snp stores isolated SNP kmer pair, k_31_pair.non and k_31_pair.non.sep store non-isolated SNP kmer pair. 
k_31_pair.non.sep puts one SNP in the middle position. 
k_31_pair.non stores the result after merge two non-isolated SNP kmer pair. 

# Example Usage (first version)

<pre><code>  sh runkmercalling.sh 31 60 single_chr22_60x_data.fq  </code></pre>

31 is the size of k-mer,  60 is reads coverage (or homozygous coverage ), single_chr22_60x_data.fq is reads fastq file.

For more than one fastq file
<pre><code>  sh runkmercalling.sh 31 60 reads1.fq,reads2.fq,reads3.fq  </code></pre>

Input fasta file
<pre><code>  sh runkmercalling.sh 31 60 single_chr22_60x_data.fasta  </code></pre>

# Example Usage (version > 1)

<pre><code> python3 kmer2snp.py single --k 31 --c 60 --fastaq example.fasta   </code></pre>

To see more input parameters,
<pre><code> python3 kmer2snp.py single --help </code></pre>

Ideally, users will have multiple samples for comparisons. Therefore we have implemented Kmer2SNP such that multiple samples can be studied collectively in population mode.

Typical command for population mode is as follows:

<pre><code> python3 kmer2snp.py population --k 31 --cfile sample.coverage.txt --faqfile
[file with per sample fastq/fasta file] </code></pre>

--cfile: a file contains homozygous coverage information for each sample. Each row of this file must contain information about one sample; the first column must list the sample identifier (no spaces allowed) and the second column lists the numeric homozygous coverage value. 
--faqfile: a text file containing list of NGS sequence reads file of all samples. Each row of the file must contain information about one sample; the first column must list the sample identifier (no spaces allowed) and the second column must contain path to at least one sequence read file in fastq/fasta format. If the sample has multiple files, file paths can be separated by a comma (”,”). Sample identifiers must match with sample identifiers supplied to the --cfile.


To see more input parameters,
<pre><code> python3 kmer2snp.py population --help </code></pre>





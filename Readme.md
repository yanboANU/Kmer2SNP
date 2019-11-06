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

  
Then clone the **Kmer2SNP** repository to your machine.
<pre><code> git clone https://github.com/yanboANU/Kmer2SNP.git </code></pre>

After you download Kmer2SNP, you need to change path in the runkmercalling.sh.

# Input Format

The input of **Kmer2SNP**  can be fasta, fastq, either gzipped or not. 

# Ouput Format

See Kmer2SNP/example, k_31_pair.snp stores isolated SNP kmer pair, k_31_pair.non and k_31_pair.non.sep store non-isolated SNP kmer pair. 
k_31_pair.non.sep puts one SNP in the middle position. 
k_31_pair.non stores the result after merge two non-isolated SNP kmer pair. 

# Example Usage

<pre><code>  sh runkmercalling.sh 31 60 single_chr22_60x_data.fq  </code></pre>

31 is the size of k-mer,  60 is reads coverage (or homozygous coverage ) .



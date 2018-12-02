# WillItCode

Compute the coding potential of RNA sequences.

## Overview

WillItCode is tool to classify RNA sequences into protein-coding (mRNA) or non-coding (lncRNA). It computes a series of features that describe the sequence composition and similarity to protein motifs and then classify them using a gradient boosted trees model.

## Preparing the work environment

1. Install the [Conda](https://conda.io/) package and environment manager. This can be achieved by the installation of [Miniconda](https://conda.io/miniconda.html) or [Anaconda](https://www.anaconda.com/download/#linux).
2. Create the WillItCode environment: `conda env create -f willitcode_environment.yml`
3. Activate the WillItCode environment: `source activate willitcode`
4. Download the latest version of the [Pfam database](http://pfam.xfam.org/) into the `Pfam` directory and unpack it.
5. Prepare the HMM database: `hmmpress Pfam-A.hmm`

## Scripts
- `willitcode_make_hexamer_table.py`:
  - Outputs a file containing the hexamer frequencies in the ORFs of coding and non-coding transcripts.
  - Take two FASTA files as input, one of the training mRNA sequences and the other containing the training lncRNA.
- `willitcode_train_classification_model.py`:
  - Outputs the trained classification model file.
  - Take the two training FASTA files and the hexamer frequency file as input.
- `willitcode_classify_fasta.py`:
  - Outputs the classification of the sequences in the target FASTA and (optionally) the values of each feature.
  - Take the hexamer frequency file, the classification model file and the target FASTA as input.
- `willitcode_make_feature_matrix.py`:
  - Outputs a file containing the values of the features for each sequence in the target FASTA file.
  - Take the hexamer frequency file and the target FASTA as input.

## Example

```
$ git clone https://github.com/apcamargo/WillItCode.git
$ cd WillItCode/
$ conda env create -f willitcode_environment.yml
$ source activate willitcode
$ cd Pfam/
$ wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz && gunzip Pfam-A.hmm.gz && hmmpress Pfam-A.hmm
$ cd ../
$ python willitcode_make_hexamer_table.py ../mRNA_training.fa ../lncRNA_training.fa hexamer_frequencies.pickle
$ python willitcode_train_classification_model.py --hmmer_cpu 8 ../mRNA_training.fa ../lncRNA_training.fa hexamer_frequencies.pickle model.pickle
$ python willitcode_classify_fasta.py --output_file classification_result.tsv --output_features --hmmer_cpu 8 ../test_transcripts.fa model.pickle hexamer_frequencies.pickle
```

The ouput of `willitcode_classify_fasta.py` contains the computed features and the classification of each transcript:

```
sequence_id	orf_integrity	log_sequence_length	log_orf_length	orf_ratio	fickett_score	gc_content	gc_bias	hexamer_bias	hexamer_bias_distance	protein_pi	snr	log_hmmer_score	coding_probability	prediction
HOTAIR-201	1.0	7.792348924113037	5.7745515455444085	0.1325898389095415	0.9457	60.74766355140187	5.607476635514018	0.07669675570903184	0.14526636980902713	11.83282470703125	0.9328568042965191	0.0	0.05229618772864342	Non-coding
GAPDH-201	1.0	7.53689712956617	6.9167150203536085	0.5376	1.2926	55.257936507936506	7.142857142857146	0.2934898815462439	0.5006784307968564	8.56610107421875	3.3520498724216043	5.454038241544812	0.9477038383483887	Coding
```

## Features

## Features

- `orf_integrity`: Whether a stop codon was found.
- `log_sequence_length`: Log-transformed transcript length.
- `log_orf_length`: Log-transformed longest ORF length.
- `orf_ratio`: Ratio between the longest ORF length and the transcript length.
- `fickett_score`: Score computed using the ORF nucleotide composition as described by JW Fickett (1982).
- `gc_content`: GC content of the longest ORF.
- `gc_bias`: Difference between the largest GC content among the reading frames and the ORF GC content.
- `hexamer_bias`: Score computed by measusing how similar the transcript's hexamer composition is to the average hexamer composition of mRNAs and lncRNAs.
- `hexamer_bias_distance`: Average difference between the value of the largest `hexamer_bias` among the reading frames and the values of the other two frames.
- `protein_pi`: Isoelectric point of the protein translated from the longest ORF.
- `snr`: A quantification of signal-to-noise ratio in the ORF as described by Pian, Cong, et al. (2016). Uses the Discrete Fourier Transform to detect a period-3 behaviour and then computes how strong is this signal.
- `log_hmmer_score`: Log-transformed score of the first hit found by `hmmsearch` in the Pfam database. Uses the protein sequence translated from the longest ORF.

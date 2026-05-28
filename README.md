# Genomic Fasta to Music

Genomic Fasta to Music is a Python-based creative bioinformatics prototype that converts DNA sequences into sound. The script reads a genomic FASTA file, cleans the nucleotide sequence, translates codons into amino acids, calculates GC content, and maps amino acids to musical notes.

The project explores how biological sequence information can be transformed into audio using basic computational biology and sound synthesis concepts.

## Overview

This prototype takes a DNA sequence and generates a short audio composition based on:

* Codon-to-amino acid translation
* GC content
* Amino acid chemical properties
* Major and minor musical scales
* ADSR envelope modulation
* Sine wave audio synthesis

Depending on the GC content of the sequence, the program selects either a major or minor scale. Amino acids are then mapped to specific note frequencies, while amino acid properties such as hydrophobic, polar, or charged influence the sound envelope.

## Features

* Reads genomic FASTA files
* Removes FASTA headers and non-nucleotide characters
* Supports lead strand, reverse strand, and reverse complement strand
* Translates DNA codons into amino acid sequences
* Calculates GC content
* Maps amino acids to musical frequencies
* Uses ADSR modulation to shape each note
* Plays generated audio directly using Python

## Requirements

This project uses:

```bash
numpy
sounddevice
```

Install the required packages with:

```bash
pip install numpy sounddevice
```

## How It Works

1. The script loads a genomic FASTA file.
2. The DNA sequence is cleaned to keep only A, T, C, and G bases.
3. The user selects which strand to use:

   * Lead strand
   * Reverse complement strand
   * Reverse strand
4. The selected DNA sequence is translated into amino acids.
5. GC content is calculated.
6. Amino acids are converted into musical notes.
7. Amin acid properties modify ADSR sound parameters.
8. The final audio sequence is played through the computer speakers.

## Biological Logic

The script uses the standard genetic code to translate codons into amino acids. GC content is used as a simple biological feature to influence musical scale selection.

Amino acids are grouped into broad chemical categories:

* Hydrophobic / non-polar
* Polar
* Charged

These categories affect sound parameters such as sustain, decay, release, and attack.

## Audio Logic

Each amino acid is mapped to a note frequency. The note is generated as a sine wave and then shaped using an ADSR envelope:

* Attack: how quickly the sound starts
* Decay: how quickly the sound drops after the initial peak
* Sustain: the stable volume level of the note
* Release: how quickly the note fades out

## Usage

Update the file path in the script:

```python
file_1 = open("path/to/your/genomic_file.fna", "r")
```

Then run:

```bash
python protein_music.py
```

Choose a strand option when prompted:

```bash
Press 0 for Lead strand | Press 1 for Reverse Complement strand | Press 2 for Reverse strand:
```

## Project Status

This project is an early-stage prototype. It is intended as a creative exploration of bioinformatics, DNA translation, and audio synthesis.

Future improvements may include:

* Command-line arguments for input files
* Exporting generated audio as WAV files
* Support for multiple reading frames
* More biologically meaningful ORF detection
* Better sound design and waveform options
* Visualization of GC content and amino acid composition

## Disclaimer

This project is for educational and creative purposes. It is not intended for biological diagnosis, formal genomic annotation, or clinical interpretation.

## Author

Enrique Garcia-Mireles
M.S. Bioinformatics student
Biotechnology Engineer

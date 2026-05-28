import numpy as np
import sounddevice as sd
import re

with open("example_genome.fna", "r") as file_1:
text = file_1.read()

# Remove FASTA header
text = re.sub(r'^>.*\n?', '', text)

# Keep only nucleotide letters
sequence = re.sub(r'[^ATCGatcg]', '', text).upper()

aminoacids = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y",
    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C",
    "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "TAA": "*", "TAG": "*", "TGA": "*"
}

print("Original sequence (first 200 bases):\n", sequence[0:200])
print("Length:", len(sequence), "bases")

seq_rev = sequence[::-1]

rev_complement = ""

for base in seq_rev:
    if base == "A":
        rev_complement += "T"
    elif base == "T":
        rev_complement += "A"
    elif base == "G":
        rev_complement += "C"
    elif base == "C":
        rev_complement += "G"


class protein2music:
    def __init__(self, using_seq, aminoacids):
        self.using_seq = using_seq
        self.aminoacids = aminoacids
        self.protein_seq = ""

    def seq2prot(self):
        protein = ""

        for i in range(0, len(self.using_seq) - 2, 3):
            codon = self.using_seq[i:i+3]

            if codon in self.aminoacids:
                protein += self.aminoacids[codon]
            else:
                protein += "X"

        self.protein_seq = protein
        return protein

    def gc_content(self):
        gcount = 0
        ccount = 0

        for base in self.using_seq:
            if base == "G":
                gcount += 1
            elif base == "C":
                ccount += 1

        gc_cont = ((gcount + ccount) / len(self.using_seq)) * 100
        return gc_cont

    def modulate_adsr(self, note, attack, sustain, release, decay, fs):
        attack_len = int(attack * fs)
        decay_len = int(decay * fs)
        release_len = int(release * fs)

        sustain_len = len(note) - attack_len - decay_len - release_len

        if sustain_len < 0:
            sustain_len = 0

        attack_env = np.linspace(0, 1, attack_len)
        decay_env = np.linspace(1, sustain, decay_len)
        sustain_env = np.ones(sustain_len) * sustain
        release_env = np.linspace(sustain, 0, release_len)

        envelope = np.concatenate([
            attack_env,
            decay_env,
            sustain_env,
            release_env
        ])

        if len(envelope) < len(note):
            extra = np.ones(len(note) - len(envelope)) * sustain
            envelope = np.concatenate([envelope, extra])

        envelope = envelope[:len(note)]

        return note * envelope

    def prot2notes(self, protein, gc_cont, attack, sustain, release, decay):
        fs = 44100
        duration = 0.12
        t = np.linspace(0, duration, int(fs * duration), endpoint=False)

        if gc_cont > 50:
            notes = {
                "G": 261.63,   # C4
                "A": 293.66,   # D4
                "V": 311.13,   # Eb4
                "L": 349.23,   # F4
                "I": 392.00,   # G4
                "M": 415.30,   # Ab4
                "F": 466.16,   # Bb4
                "W": 523.25,   # C5
                "P": 587.33,   # D5

                "S": 622.25,   # Eb5
                "T": 698.46,   # F5
                "N": 783.99,   # G5
                "Q": 830.61,   # Ab5
                "C": 932.33,   # Bb5
                "Y": 1046.50,  # C6

                "D": 1174.66,  # D6
                "E": 1244.51,  # Eb6
                "K": 1396.91,  # F6
                "R": 1567.98,  # G6
                "H": 1661.22   # Ab6
            }
        else:
            notes = {
                "G": 261.63,   # C4
                "A": 293.66,   # D4
                "V": 329.63,   # E4
                "L": 349.23,   # F4
                "I": 392.00,   # G4
                "M": 440.00,   # A4
                "F": 493.88,   # B4
                "W": 523.25,   # C5
                "P": 587.33,   # D5

                "S": 659.25,   # E5
                "T": 698.46,   # F5
                "N": 783.99,   # G5
                "Q": 880.00,   # A5
                "C": 987.77,   # B5
                "Y": 1046.50,  # C6

                "D": 1174.66,  # D6
                "E": 1318.51,  # E6
                "K": 1396.91,  # F6
                "R": 1567.98,  # G6
                "H": 1760.00   # A6
            }

        if protein in notes:
            frequency = notes[protein]
            note = np.sin(2 * np.pi * frequency * t)
            note = self.modulate_adsr(note, attack, sustain, release, decay, fs)
        else:
            note = np.zeros_like(t)

        return note


using_seq = ""

decision_user = input("Press 0 for Lead strand | Press 1 for Reverse Complement strand | Press 2 for Reverse strand: ")

if decision_user == "1":
    using_seq = rev_complement
elif decision_user == "0":
    using_seq = sequence
elif decision_user == "2":
    using_seq = seq_rev
else:
    print("Invalid digit. Using lead strand.")
    using_seq = sequence

obj = protein2music(using_seq, aminoacids)

gc_cont = obj.gc_content()
print("GC content:", gc_cont)

protein_seq = obj.seq2prot()
print("First 50 amino acids:", protein_seq[:50])

song = []

hydrophobic = {"A", "V", "I", "L", "M", "F", "W", "Y"}
polar = {"S", "T", "N", "Q", "C"}
charged = {"D", "E", "K", "R", "H"}

attack = 0.01
sustain = 0.3
release = 0.02
decay = 0.02

for aa in protein_seq[:500]:

    if aa in hydrophobic:
        sustain += 0.02
        decay += 0.002
    else:
        sustain -= 0.015
        decay -= 0.001

    if aa in polar:
        release += 0.002
    else:
        release -= 0.001

    if aa in charged:
        attack -= 0.003
    else:
        attack += 0.001

    sustain = np.clip(sustain, 0.1, 0.9)
    attack = np.clip(attack, 0.002, 0.03)
    decay = np.clip(decay, 0.002, 0.04)
    release = np.clip(release, 0.005, 0.04)

    note = obj.prot2notes(aa, gc_cont, attack, sustain, release, decay)
    song.append(note)

song = np.concatenate(song)

sd.play(0.3 * song, 44100)
sd.wait()

file_1.close()

from nucleotides import *
import random
import sys

def func(seq):
    with open(seq) as file1:
        t1 = ""
        for i in file1.readlines():
            if i.startswith(">"):
                pass
            else:
                t1 += i.rstrip("\n")
    return t1


class bio_materials:
    """Biological material analysis"""
    """Sequence should be string type"""
    def __init__(self, seq, data_type):
        self.seq = seq
        self.data_type = data_type

        if self.data_type == "DNA":
            self.valid = self.__DNAvalidate()
            assert self.valid, f"data that you entered is not correct {data_type} type. Try to change data_type or input the correct letters"
        elif self.data_type == "RNA":
            self.valid = self.___RNAvalidate()
            assert self.valid, f"data that you entered is not correct {data_type} type. Try to change data_type or input the correct letters"
    def __DNAvalidate(self):
        return set(NUCLEOTIDE_BASE["DNA"]).issuperset(self.seq)
    def ___RNAvalidate(self):
        return set(NUCLEOTIDE_BASE["RNA"]).issuperset(self.seq)

    def nuc_numbers(self, freq = False):
        """Get the frequency of each nucleotide types"""
        frequency = dict()
        for i in self.seq:
            if i in frequency:
                frequency[i] += 1
            else:
                frequency[i] = 1
        if freq:
            for i in frequency.keys():
                frequency[i] = round(frequency[i]/len(self.seq), 2)
            return frequency
        return frequency

    def transcription(self):
        """Transcribe DNA nucleotides to RNA nucleotides"""
        if self.data_type != "DNA":
            return "just DNA can be transcribed"
        self.seq = self.seq.replace("T", "U")
        self.data_type = "RNA"
        return self

    def reverse_transcription(self):
        """Transcribe RNA to DNA sequence"""
        if self.data_type != "RNA":
            return "Just RNA can be transcribe reversely"
        self.data_type = "DNA"
        return self

    def info(self):
        """Brief information about provided sequence"""
        return f"Sequence type: {self.data_type}\nSequence length: {len(self.seq)}"

    def generate_rnd_seq(self, data_type= "DNA",length=100):
        """Delete provided sequence and generate a new random sequence"""
        seq = ''
        for j in random.choices(NUCLEOTIDE_BASE[data_type], k=length):
            seq += j
        self.__init__(seq, data_type)

        return self

    def reverse_complementary(self):
        """Take the DNA sequence and find the reverse and complementary sequence"""
        if self.data_type != "DNA" and self.data_type != "RNA":
            return f"Just DNA and RNA can have a reverse complementary sequence"
        seq = self.seq[::-1]
        new_seq = ''
        x = eval("Complementary_"+ self.data_type)
        for i in seq:
            new_seq += x[i]
        self.seq = new_seq

        # for i in seq:
        #     if i == "G":
        #         new_seq += "C"
        #     if i == "C":
        #         new_seq += "G"
        #     if i == "U" or i == "T":
        #         new_seq += "A"
        #     if i == "A" and self.data_type == "DNA":
        #         new_seq += "T"
        #     if i == "A" and self.data_type == "RNA":
        #         new_seq += "U"
        return self

    def colorize(self):
        """Colorize the nucleotides"""
        col = {
            'A': '\033[33m',
            'U': '\033[34m',
            'G': '\033[35m',
            'T': '\033[36m',
            'C': '\033[37m',
            'reset': '\033[0m'
        }
        StrCol = ""

        for i in self.seq:
            if i in col:
                StrCol += col[i] + i
            else:
                StrCol += col['reset'] + i

        return StrCol + '\033[0m'

    def GC_content(self):
        """Find the GC content of nucleotide sequence"""
        t = 0
        for i in self.seq:
            if i == 'G' or i == 'C':
                t += 1
        return 10 * (t / len(self.seq))

    def GC_content_sub(self, k=5):
        """"Calculate GC content of a sequence by dividing k length"""
        seq = []
        if len(self.seq) % k == 0:

            for i in range(0, len(self.seq), k):
                rat = self.seq[i:i + k]
                seq.append((rat.count("G")+rat.count("C"))/len(rat))
        elif k > len(self.seq):
            return sys.stderr.write("K is bigger than X \n")
            sys.stderr.flush()

        else:
            for i in range(0, len(self.seq) - k + 1, k):
                rat = self.seq[i:i + k]
                seq.append((rat.count("G")+rat.count("C"))/len(rat))
        return seq

    def translation(self, start_point=0):
        """Translate DNA or RNA to Amino Acid sequence"""
        seq = ''
        x = eval(self.data_type + "_Codons")
        if len(self.seq) % 3 == 0:
            for i in range(start_point, len(self.seq), 3):
                seq += x[self.seq[i:i + 3]]
        else:
            for i in range(start_point, len(self.seq)-2, 3):
                seq += x[self.seq[i:i + 3]]
        return seq

    def codon_frequency(self, letter = "L"):
        """Get the codons that translated as Amino Acids which is indicated as letter"""
        """Input should be RNA or DNA"""
        if self.data_type != "DNA" and self.data_type != "RNA":
            return "Input should be RNA or DNA"
        letters = {}
        letter = letter.upper()
        x = eval(self.data_type + "_Codons")
        if len(self.seq) % 3 == 0:
            for i in range(0, len(self.seq), 3):
                if (letter == x[self.seq[i:i + 3]]) and (self.seq[i:i + 3] in letters):
                    letters[self.seq[i:i + 3]] += 1
                elif (letter == x[self.seq[i:i + 3]]) and (self.seq[i:i + 3] not in letters):
                    letters[self.seq[i:i + 3]] = 1
                else:
                    pass
        else:
            for i in range(0, len(self.seq) - 2, 3):
                if (letter == x[self.seq[i:i + 3]]) and (self.seq[i:i + 3] in letters):
                    letters[self.seq[i:i + 3]] += 1
                elif (letter == x[self.seq[i:i + 3]]) and (self.seq[i:i + 3] not in letters):
                    letters[self.seq[i:i + 3]] = 1
                else:
                    pass
        total = sum(letters.values())
        for i, j in letters.items():
            letters[i] = j / total
        return letter, letters

    def open_reading_frame(self):
        """Get the 6 possible Amino Acids sequence from DNA or RNA"""
        str1 = ''
        str1 += (self.translation(0) + "\n\n")
        str1 += (self.translation(1) + "\n\n")
        str1 += (self.translation(2) + "\n\n")

        str1 += (self.reverse_complementary().translation(0) + "\n\n")
        str1 += (self.reverse_complementary().translation(1) + "\n\n")
        str1 += (self.reverse_complementary().translation(2) + "\n\n")

        return str1

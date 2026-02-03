import re

DNA_ALPHABET = {"A", "C", "G", "T"}
IUPAC = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"}, "W": {"A", "T"},
    "K": {"G", "T"}, "M": {"A", "C"},
    "B": {"C", "G", "T"}, "D": {"A", "G", "T"}, "H": {"A", "C", "T"},
    "V": {"A", "C", "G"}, "N": {"A", "C", "G", "T"}, ".": {"A", "C", "G", "T"},
}

def parse_bracket_class(s, i):
    assert s[i] == "["
    j = i + 1
    allowed = set()
    if j >= len(s):
        raise ValueError("Unclosed '[' in motif.")
    while j < len(s) and s[j] != "]":
        ch = s[j].upper()
        if ch in IUPAC:
            allowed |= IUPAC[ch]
        elif ch in DNA_ALPHABET:
            allowed.add(ch)
        else:
            raise ValueError(f"Unsupported character '{s[j]}' in bracket class.")
        j += 1
    if j >= len(s) or s[j] != "]":
        raise ValueError("Unclosed '[' in motif.")
    return allowed, j + 1


def parse_motif_to_symbol_sets(motif):
    s = motif.strip().upper()
    i = 0
    pieces = []
    if not s:
        raise ValueError("Motif is empty.")
    while i < len(s):
        ch = s[i]
        if ch == "[":
            allowed, i = parse_bracket_class(s, i)
            pieces.append(allowed)
        elif ch in IUPAC:
            pieces.append(IUPAC[ch].copy())
            i += 1
        else:
            raise ValueError(f"Unsupported motif character '{s[i]}' at position {i}.")
    return pieces

def read_fasta(text):
    entries = []
    header = None
    seq_lines = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                entries.append((header, "".join(seq_lines).upper()))
            header = line[1:].strip()
            seq_lines = []
        else:
            seq_lines.append(re.sub(r"\s+", "", line))
    if header is not None:
        entries.append((header, "".join(seq_lines).upper()))
    return entries

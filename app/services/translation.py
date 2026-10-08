"""DNA to protein translation with the standard genetic code (NCBI table 1).

Rules (from the specification):

* frame 1/2/3 starts at nucleotide 1/2/3 of the supplied strand: no start-codon
  search, no reverse complement;
* only complete codons are translated; 1-2 trailing nucleotides are ignored;
* any codon containing N becomes ``X``, even if the amino acid could be inferred
  (e.g. ``GCN`` is always Ala, but is still reported as ``X``);
* translation stops *before* the first unambiguous stop codon (TAA, TAG, TGA).
  ``*`` is never part of the protein.
"""

from __future__ import annotations

from dataclasses import dataclass

VALID_FRAMES = (1, 2, 3)
UNKNOWN_AMINO_ACID = "X"
STOP = "*"

# The standard code in the conventional TCAG order: codon (b1, b2, b3) maps to
# _AMINO_ACIDS[16 * i1 + 4 * i2 + i3], where i is the index of each base in "TCAG".
_BASES = "TCAG"
_AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

CODON_TABLE: dict[str, str] = {
    b1 + b2 + b3: _AMINO_ACIDS[16 * i1 + 4 * i2 + i3]
    for i1, b1 in enumerate(_BASES)
    for i2, b2 in enumerate(_BASES)
    for i3, b3 in enumerate(_BASES)
}
STOP_CODONS = frozenset(codon for codon, aa in CODON_TABLE.items() if aa == STOP)


class InvalidFrameError(ValueError):
    pass


@dataclass(frozen=True, slots=True)
class Translation:
    frame: int
    protein: str
    has_stop_codon: bool

    @property
    def length(self) -> int:
        return len(self.protein)


def translate_codon(codon: str) -> str:
    """Translate one codon; ``X`` if it contains N, ``*`` for a stop codon."""
    if len(codon) != 3:
        raise ValueError(f"codon must have 3 bases, got {codon!r}")
    if "N" in codon:
        return UNKNOWN_AMINO_ACID
    return CODON_TABLE[codon]


def codons(sequence: str, frame: int = 1) -> list[str]:
    """Split ``sequence`` into the complete codons of the given frame."""
    _check_frame(frame)
    start = frame - 1
    return [sequence[i : i + 3] for i in range(start, len(sequence) - 2, 3)]


def translate(sequence: str, frame: int = 1) -> Translation:
    """Translate an upper-case DNA sequence until the first stop codon."""
    protein: list[str] = []
    for codon in codons(sequence, frame):
        amino_acid = translate_codon(codon)
        if amino_acid == STOP:
            return Translation(frame=frame, protein="".join(protein), has_stop_codon=True)
        protein.append(amino_acid)
    return Translation(frame=frame, protein="".join(protein), has_stop_codon=False)


def translate_full(sequence: str, frame: int = 1) -> str:
    """Translate every complete codon, keeping ``*`` for stops (used for annotation)."""
    return "".join(translate_codon(codon) for codon in codons(sequence, frame))


def _check_frame(frame: int) -> None:
    if frame not in VALID_FRAMES:
        raise InvalidFrameError(f"frame must be 1, 2 or 3, got {frame!r}")

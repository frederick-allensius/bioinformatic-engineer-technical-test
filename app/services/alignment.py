"""Global pairwise alignment (Needleman-Wunsch), implemented from scratch.

Scoring
-------
The default scheme is the one used in the task's example response, so the example
(patient_02, score 294.0) reproduces exactly:

=====================  =====
column                 score
=====================  =====
identical A/C/G/T      +1
different A/C/G/T      -1
any column with an N    0   (unknown base: neither rewarded nor penalised)
gap character          -2   (linear: a gap of length k costs 2k)
=====================  =====

Tie-breaking
------------
With a linear gap penalty, many alignments share the optimal score, including some
where one 3 bp deletion is split into a 2 bp and a 1 bp gap. To make the result
deterministic *and* biologically sensible, ties are broken in this order:

1. **Fewest gap openings.** Among all alignments with the optimal score, keep the one
   with the fewest separate gap runs, so consecutive indel bases stay grouped. This is
   done with a three-state (Gotoh) DP where every score is scaled by ``K = n + m + 1``
   and each gap opening costs an extra ``-1``. Because an alignment can never have
   ``K`` or more gap openings, this secondary term can only decide between alignments
   whose primary scores are equal; it never changes the optimal score itself.
2. **Fixed traceback preference:** match/mismatch, then deletion, then insertion.
   Traceback starts at the 3' end, so preferring the diagonal pushes gaps as far
   5' (left) as possible. In repeats this gives left-aligned indels, the same
   convention VCF normalisation uses.

The returned score is always recomputed from the alignment columns, so it is exactly
the primary score of the alignment that is returned.
"""

from __future__ import annotations

from dataclasses import dataclass

GAP = "-"
UNAMBIGUOUS = frozenset("ACGT")

# Sentinel for impossible DP states. Plain ints keep the inner loop fast and exact.
_NEG = -(10**15)
# Traceback states.
_MATCH, _DEL, _INS = 0, 1, 2


@dataclass(frozen=True, slots=True)
class ScoringScheme:
    match: int = 1
    mismatch: int = -1
    gap: int = -2
    n: int = 0

    def pair(self, a: str, b: str) -> int:
        """Score of aligning base ``a`` against base ``b`` (no gaps)."""
        if a == "N" or b == "N":
            return self.n
        return self.match if a == b else self.mismatch


DEFAULT_SCORING = ScoringScheme()


@dataclass(frozen=True, slots=True)
class Alignment:
    reference_aligned: str
    query_aligned: str
    score: int

    @property
    def length(self) -> int:
        return len(self.reference_aligned)

    @property
    def identity(self) -> float:
        """Identical A/C/G/T columns divided by *all* columns (gaps and N included)."""
        return identity(self.reference_aligned, self.query_aligned)


def identity(reference_aligned: str, query_aligned: str) -> float:
    if not reference_aligned:
        return 0.0
    identical = sum(
        1
        for r, q in zip(reference_aligned, query_aligned, strict=True)
        if r == q and r in UNAMBIGUOUS
    )
    return identical / len(reference_aligned)


def score_alignment(
    reference_aligned: str, query_aligned: str, scoring: ScoringScheme = DEFAULT_SCORING
) -> int:
    """Score an existing alignment column by column."""
    total = 0
    for r, q in zip(reference_aligned, query_aligned, strict=True):
        if r == GAP and q == GAP:
            raise ValueError("an alignment column cannot be gap/gap")
        total += scoring.gap if GAP in (r, q) else scoring.pair(r, q)
    return total


def global_align(reference: str, query: str, scoring: ScoringScheme = DEFAULT_SCORING) -> Alignment:
    """Optimal global alignment of ``query`` against ``reference``.

    Runs in O(n*m) time and memory (n, m <= 1000 here, about 1M cells per matrix
    at most). Both sequences are expected to be upper-case.
    """
    n, m = len(reference), len(query)
    k = n + m + 1  # strictly greater than the maximum possible number of gap openings
    gap = scoring.gap * k
    # Substitution scores precomputed per reference base, so the inner loop does no
    # function calls.
    substitution = {base: [scoring.pair(base, q) * k for q in query] for base in set(reference)}

    # M: ends in ref[i-1]/query[j-1]; D: ends in ref[i-1]/gap; I: ends in gap/query[j-1].
    M = [[_NEG] * (m + 1) for _ in range(n + 1)]
    D = [[_NEG] * (m + 1) for _ in range(n + 1)]
    I = [[_NEG] * (m + 1) for _ in range(n + 1)]  # noqa: E741 (standard DP notation)
    M[0][0] = 0
    for i in range(1, n + 1):
        D[i][0] = gap * i - 1
    for j in range(1, m + 1):
        I[0][j] = gap * j - 1

    for i in range(1, n + 1):
        sub = substitution[reference[i - 1]]
        m_row, d_row, i_row = M[i], D[i], I[i]
        m_up, d_up, i_up = M[i - 1], D[i - 1], I[i - 1]
        for j in range(1, m + 1):
            # (mis)match: diagonal move from any state
            best = m_up[j - 1]
            if d_up[j - 1] > best:
                best = d_up[j - 1]
            if i_up[j - 1] > best:
                best = i_up[j - 1]
            m_row[j] = best + sub[j - 1]
            # deletion: vertical move; opening a gap (from M or I) costs the extra 1
            best = m_up[j] - 1
            if d_up[j] > best:
                best = d_up[j]
            if i_up[j] - 1 > best:
                best = i_up[j] - 1
            d_row[j] = best + gap
            # insertion: horizontal move
            best = m_row[j - 1] - 1
            if d_row[j - 1] - 1 > best:
                best = d_row[j - 1] - 1
            if i_row[j - 1] > best:
                best = i_row[j - 1]
            i_row[j] = best + gap

    ref_out, query_out = _traceback(reference, query, M, D, I, substitution, gap)
    return Alignment(ref_out, query_out, score_alignment(ref_out, query_out, scoring))


def _traceback(reference, query, M, D, I, substitution, gap):  # noqa: E741
    i, j = len(reference), len(query)
    matrices = (M, D, I)
    state = _pick(M[i][j], D[i][j], I[i][j], target=max(M[i][j], D[i][j], I[i][j]))
    ref_out: list[str] = []
    query_out: list[str] = []

    while i > 0 or j > 0:
        value = matrices[state][i][j]
        if state == _MATCH:
            prev = value - substitution[reference[i - 1]][j - 1]
            ref_out.append(reference[i - 1])
            query_out.append(query[j - 1])
            i, j = i - 1, j - 1
            state = _pick(M[i][j], D[i][j], I[i][j], target=prev)
        elif state == _DEL:
            prev = value - gap
            ref_out.append(reference[i - 1])
            query_out.append(GAP)
            i -= 1
            state = _pick(M[i][j] - 1, D[i][j], I[i][j] - 1, target=prev)
        else:
            prev = value - gap
            ref_out.append(GAP)
            query_out.append(query[j - 1])
            j -= 1
            state = _pick(M[i][j] - 1, D[i][j] - 1, I[i][j], target=prev)

    return "".join(reversed(ref_out)), "".join(reversed(query_out))


def _pick(from_match: int, from_del: int, from_ins: int, *, target: int) -> int:
    """Return the predecessor state that produced ``target``, preferring M > D > I."""
    if from_match == target:
        return _MATCH
    if from_del == target:
        return _DEL
    if from_ins == target:
        return _INS
    raise AssertionError("traceback inconsistency")  # pragma: no cover

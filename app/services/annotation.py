"""Variant effect annotation in reference frame 1 (bonus ``/api/analyze``).

Effects
-------
``synonymous`` / ``missense`` / ``nonsense`` / ``stop_lost``
    SNVs, judged on the whole codon. Several SNVs in the **same codon** are evaluated
    together (the codon is read straight from the alignment), so two changes that
    would each be missense but jointly synonymous are reported as synonymous.
``frameshift`` / ``inframe_insertion`` / ``inframe_deletion``
    Indels, decided by length (``len % 3``). Each indel is applied to the reference
    *on its own* and the two translations are compared to describe the protein change.
``undetermined``
    The effect cannot be decided reliably: the codon contains an N, the codon is also
    hit by an indel, or the variant lies outside the last complete reference codon.

Each variant is annotated against the reference frame (like VEP does). Haplotype-level
interactions are only flagged in ``notes``: a variant downstream of a frameshift or of
a premature stop, or a later indel that restores the frame. The ``translation`` field
of the analyze response shows the protein actually encoded by the sample.

``protein_change`` uses simplified one-letter HGVS-like notation (``p.L82*``,
``p.N36del``, ``p.P6_L7insI``, ``p.K77Xfs*12``). No ``dup``/``ext`` detection beyond
``stop_lost``.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, replace
from typing import Literal, get_args

from app.services.alignment import GAP, Alignment
from app.services.translation import STOP, UNKNOWN_AMINO_ACID, translate_codon, translate_full
from app.services.variants import Variant

Effect = Literal[
    "synonymous",
    "missense",
    "nonsense",
    "stop_lost",
    "frameshift",
    "inframe_insertion",
    "inframe_deletion",
    "undetermined",
]
EFFECTS: tuple[Effect, ...] = get_args(Effect)


@dataclass(frozen=True, slots=True)
class AnnotatedVariant:
    variant: Variant
    effect: Effect
    codon: int | None = None  # 1-based codon number in reference frame 1
    ref_codon: str | None = None
    alt_codon: str | None = None
    ref_aa: str | None = None
    alt_aa: str | None = None
    protein_change: str | None = None
    notes: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class _AlignedView:
    """The alignment projected onto reference coordinates."""

    query_at: list[str | None]  # query base aligned to reference position p (index p, 1-based)
    inserted_after: frozenset[int]  # reference positions followed by an insertion

    @classmethod
    def build(cls, alignment: Alignment, reference_length: int) -> _AlignedView:
        query_at: list[str | None] = [None] * (reference_length + 1)
        inserted_after: set[int] = set()
        ref_pos = 0
        for r, q in zip(alignment.reference_aligned, alignment.query_aligned, strict=True):
            if r == GAP:
                inserted_after.add(ref_pos)
            else:
                ref_pos += 1
                query_at[ref_pos] = q
        return cls(query_at, frozenset(inserted_after))


def annotate_variants(
    reference: str, alignment: Alignment, variants: Sequence[Variant]
) -> list[AnnotatedVariant]:
    view = _AlignedView.build(alignment, len(reference))
    snv_positions = {v.position for v in variants if v.type == "SNV"}
    ref_protein = translate_full(reference)

    annotated = [
        _annotate_snv(v, reference, view, snv_positions)
        if v.type == "SNV"
        else _annotate_indel(v, reference, ref_protein)
        for v in variants
    ]
    return _add_context_notes(annotated)


# --------------------------------------------------------------------------- SNVs


def _annotate_snv(
    variant: Variant, reference: str, view: _AlignedView, snv_positions: set[int]
) -> AnnotatedVariant:
    codon_no = (variant.position - 1) // 3 + 1
    start = 3 * codon_no - 2  # 1-based first base of the codon
    if start + 2 > len(reference):
        return AnnotatedVariant(
            variant, "undetermined", notes=("Outside the last complete reference codon.",)
        )

    span = range(start, start + 3)
    ref_codon = reference[start - 1 : start + 2]
    if any(view.query_at[p] == GAP for p in span) or view.inserted_after & {start, start + 1}:
        return AnnotatedVariant(
            variant,
            "undetermined",
            codon_no,
            ref_codon=ref_codon,
            notes=("Codon is also affected by an insertion or deletion.",),
        )

    alt_codon = "".join(view.query_at[p] or "" for p in span)
    ref_aa, alt_aa = translate_codon(ref_codon), translate_codon(alt_codon)
    notes: list[str] = []
    partners = [p for p in span if p != variant.position and p in snv_positions]
    if partners:
        joined = ", ".join(map(str, partners))
        notes.append(f"Evaluated together with SNV(s) at {joined} in the same codon.")

    if UNKNOWN_AMINO_ACID in (ref_aa, alt_aa):
        notes.append("Codon contains N, so the amino acid cannot be determined.")
        return AnnotatedVariant(
            variant,
            "undetermined",
            codon_no,
            ref_codon=ref_codon,
            alt_codon=alt_codon,
            notes=tuple(notes),
        )

    effect = _classify_substitution(ref_aa, alt_aa)
    if codon_no == 1 and ref_aa == "M" and alt_aa != "M":
        notes.append("Changes the start codon (start_lost).")
    return AnnotatedVariant(
        variant,
        effect,
        codon_no,
        ref_codon,
        alt_codon,
        ref_aa,
        alt_aa,
        _substitution_hgvs(effect, ref_aa, alt_aa, codon_no),
        tuple(notes),
    )


def _classify_substitution(ref_aa: str, alt_aa: str) -> Effect:
    if ref_aa == alt_aa:
        return "synonymous"
    if alt_aa == STOP:
        return "nonsense"
    if ref_aa == STOP:
        return "stop_lost"
    return "missense"


def _substitution_hgvs(effect: Effect, ref_aa: str, alt_aa: str, codon_no: int) -> str:
    if effect == "synonymous":
        return f"p.{ref_aa}{codon_no}="
    if effect == "stop_lost":
        return f"p.*{codon_no}{alt_aa}ext*?"
    return f"p.{ref_aa}{codon_no}{alt_aa}"


# --------------------------------------------------------------------------- indels


def _annotate_indel(variant: Variant, reference: str, ref_protein: str) -> AnnotatedVariant:
    if variant.length % 3:
        effect: Effect = "frameshift"
    else:
        effect = "inframe_insertion" if variant.type == "insertion" else "inframe_deletion"

    # Insertion after base p starts affecting codon p // 3 + 1; deletion at p affects
    # the codon containing p.
    pos = variant.position
    codon_no = pos // 3 + 1 if variant.type == "insertion" else (pos - 1) // 3 + 1
    if codon_no > len(reference) // 3:
        return AnnotatedVariant(
            variant, effect, notes=("Outside the last complete reference codon.",)
        )

    alt_protein = translate_full(_apply(reference, variant))
    ref_aa, alt_aa, change, notes = _describe_protein_change(
        ref_protein, alt_protein, frameshift=effect == "frameshift"
    )
    return AnnotatedVariant(
        variant,
        effect,
        codon_no,
        ref_aa=ref_aa,
        alt_aa=alt_aa,
        protein_change=change,
        notes=notes,
    )


def _apply(reference: str, variant: Variant) -> str:
    """Reference with this single indel applied."""
    if variant.type == "insertion":
        return reference[: variant.position] + variant.alt + reference[variant.position :]
    start = variant.position - 1
    return reference[:start] + reference[start + len(variant.ref) :]


def _describe_protein_change(
    ref_protein: str, alt_protein: str, *, frameshift: bool
) -> tuple[str | None, str | None, str | None, tuple[str, ...]]:
    """Compare full reference/alt translations; return (ref_aa, alt_aa, hgvs, notes)."""
    ref_stop = ref_protein.find(STOP)
    coding_end = len(ref_protein) if ref_stop == -1 else ref_stop + 1
    prefix = _common_prefix(ref_protein, alt_protein)
    if prefix >= coding_end:
        return None, None, "p.=", ("The encoded protein is unchanged.",)

    if frameshift:
        ref_aa = ref_protein[prefix]
        alt_aa = alt_protein[prefix] if prefix < len(alt_protein) else None
        if alt_aa == STOP:
            return ref_aa, alt_aa, f"p.{ref_aa}{prefix + 1}*", ("Creates an immediate stop codon.",)
        new_stop = alt_protein.find(STOP, prefix)
        tail = f"*{new_stop - prefix + 1}" if new_stop != -1 else "*?"
        return ref_aa, alt_aa, f"p.{ref_aa}{prefix + 1}{alt_aa or '?'}fs{tail}", ()

    # In-frame: trim the common prefix (maximal, giving the HGVS 3' rule) then suffix.
    ref_rest, alt_rest = ref_protein[prefix:], alt_protein[prefix:]
    suffix = _common_prefix(ref_rest[::-1], alt_rest[::-1])
    suffix = min(suffix, len(ref_rest), len(alt_rest))
    ref_mid = ref_rest[: len(ref_rest) - suffix]
    alt_mid = alt_rest[: len(alt_rest) - suffix]
    if not ref_mid and not alt_mid:
        return None, None, "p.=", ("The encoded protein is unchanged.",)
    first = prefix + 1  # 1-based index of the first changed residue
    notes = ("Introduces a premature stop codon.",) if STOP in alt_mid else ()

    if not alt_mid:  # pure deletion
        last = prefix + len(ref_mid)
        span = f"{ref_mid[0]}{first}" + (f"_{ref_mid[-1]}{last}" if len(ref_mid) > 1 else "")
        return ref_mid, "-", f"p.{span}del", notes
    if not ref_mid:  # pure insertion between residues `prefix` and `prefix + 1`
        if prefix == 0 or prefix >= len(ref_protein):
            return "-", alt_mid, None, notes
        left, right = ref_protein[prefix - 1], ref_protein[prefix]
        return "-", alt_mid, f"p.{left}{prefix}_{right}{prefix + 1}ins{alt_mid}", notes
    last = prefix + len(ref_mid)
    span = f"{ref_mid[0]}{first}" + (f"_{ref_mid[-1]}{last}" if len(ref_mid) > 1 else "")
    return ref_mid, alt_mid, f"p.{span}delins{alt_mid}", notes


def _common_prefix(a: str, b: str) -> int:
    n = 0
    for x, y in zip(a, b, strict=False):
        if x != y:
            break
        n += 1
    return n


# --------------------------------------------------------------------------- context


def _add_context_notes(annotated: list[AnnotatedVariant]) -> list[AnnotatedVariant]:
    """Flag variants whose reference-frame effect is overshadowed by an upstream event."""
    out: list[AnnotatedVariant] = []
    frame_offset = 0
    shifted_at: int | None = None  # position of the frameshift currently in effect
    stop_codon: int | None = None  # first premature stop codon created by an SNV

    for item in annotated:
        notes = list(item.notes)
        v = item.variant
        if stop_codon is not None and (item.codon is None or item.codon > stop_codon):
            notes.append(
                f"Downstream of the premature stop at codon {stop_codon}; "
                "not translated in this sample."
            )
        elif shifted_at is not None:
            notes.append(
                f"Downstream of the frameshift at position {shifted_at}; the "
                "reference-frame effect may not reflect the translated protein."
            )

        if v.type != "SNV":
            was_shifted = frame_offset % 3 != 0
            frame_offset += v.length if v.type == "insertion" else -v.length
            if was_shifted and frame_offset % 3 == 0:
                notes.append(f"Restores the reading frame shifted at position {shifted_at}.")
                shifted_at = None
            elif not was_shifted and frame_offset % 3 != 0:
                shifted_at = v.position
        if item.effect == "nonsense" and stop_codon is None:
            stop_codon = item.codon

        out.append(replace(item, notes=tuple(notes)))
    return out

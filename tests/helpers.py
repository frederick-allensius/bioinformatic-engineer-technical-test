"""Shared test oracles."""


def apply_variants(reference: str, variants: list[dict]) -> str:
    """Rebuild the query from the reference and the reported variants (test oracle)."""
    out, cursor = [], 0  # cursor = number of reference bases already emitted
    for v in sorted(variants, key=lambda v: (v["position"], v["type"] != "deletion")):
        pos = v["position"]
        if v["type"] == "insertion":
            out.append(reference[cursor:pos])
            out.append(v["alt"])
            cursor = pos
        else:
            out.append(reference[cursor : pos - 1])
            out.append("" if v["type"] == "deletion" else v["alt"])
            cursor = pos - 1 + len(v["ref"])
    out.append(reference[cursor:])
    return "".join(out)

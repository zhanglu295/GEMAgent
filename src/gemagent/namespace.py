"""Model namespace adapters for the reconstruction families GEMAgent supports.

Automated reconstruction tools disagree on three things that every review skill
needs to know: how a species id is spelled, how a compartment is encoded, and
how an exchange reaction is named.  Historically each of those differences was
handled ad hoc (a CarveMe ``C_e`` special case here, an AGORA ``M_`` prefix fix
there), which is exactly the source of silent no-ops such as "every ADD was
refused because the id did not exist in the model".

This module centralises the differences in one small adapter:

================  ==========================  ==============  =====================
family            species id                  compartment     exchange id
================  ==========================  ==============  =====================
``seed``          ``M_cpd00001_c0``           ``c0``/``e0``   ``R_EX_cpd00001_e0``
(gapseq,          ``M_cpd00001_c0``           ``c0``/``e0``   ``R_EX_cpd00001_e0``
 ModelSEED)
``bigg``          ``M_10fthf_c``              ``C_c``/``C_e`` ``R_EX_2ameph_e``
(CarveMe)
``agora``         ``M_10fthf__91__c__93__``   ``c``/``e``     ``R_EX_10fthf__40__e__41__``
(AGORA2)
================  ==========================  ==============  =====================

Usage::

    from gemagent.namespace import detect, for_base

    ns = detect(sbml_text)          # or for_base("carveme")
    ns.split_species("M_10fthf_c")  # -> ("10fthf", "c")
    ns.cytosolic_partner("M_10fthf_c")  # -> ""

    ns = for_base("gapseq")
    ns.cytosolic_partner("M_cpd00001_e0")  # -> "M_cpd00001_c0"

The adapter never guesses silently: an unknown layout raises ``UnknownNamespace``
with the sampled evidence, so a new tool fails loudly instead of producing a
model full of refused edits.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

DECODE = (
    ("__40__", "("), ("__41__", ")"), ("__91__", "["), ("__93__", "]"),
    ("__45__", "-"), ("__95__", "_"), ("__46__", "."), ("__58__", ":"),
)

SPECIES_PREFIX = "M_"
REACTION_PREFIX = "R_"

#: logical compartment -> suffix used by each family
#: logical compartment -> accepted suffixes (first entry is written back)
SUFFIXES = {
    "seed": {"c": ("_c0",), "e": ("_e0",), "p": ("_p0",)},
    "bigg": {"c": ("_c",), "e": ("_e",), "p": ("_p",)},
    # AGORA spells species with square brackets and the exchange with parens
    "agora": {"c": ("__91__c__93__",),
              "e": ("__91__e__93__", "__40__e__41__"),
              "p": ("__91__p__93__",)},
}

#: logical compartment -> compartment attribute used by the model
COMPARTMENTS = {
    "seed": {"c": "c0", "e": "e0", "p": "p0"},
    "bigg": {"c": "C_c", "e": "C_e", "p": "C_p"},
    "agora": {"c": "c", "e": "e", "p": "p"},
}


def decode(value: str) -> str:
    """Turn the SBML escape sequences back into readable characters."""
    text = str(value)
    for encoded, plain in DECODE:
        text = text.replace(encoded, plain)
    return text


def norm(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


class UnknownNamespace(ValueError):
    """Raised when a model does not match any known reconstruction family."""


@dataclass(frozen=True)
class Namespace:
    """Adapter for one reconstruction family.

    ``name`` is what the user asked for (``gapseq``, ``modelseed``, ...);
    ``family`` is the id convention actually in use.  gapseq and ModelSEED share
    the SEED namespace, so they share one adapter but keep distinct labels in
    reports and decision tables.
    """

    name: str
    family: str

    # -- species ---------------------------------------------------------
    @property
    def suffixes(self) -> dict[str, tuple[str, ...]]:
        return SUFFIXES[self.family]

    def strip_species_prefix(self, species_id: str) -> str:
        text = str(species_id)
        return text[len(SPECIES_PREFIX):] if text.startswith(SPECIES_PREFIX) else text

    def split_species(self, species_id: str) -> tuple[str, str]:
        """``M_cpd00001_e0`` -> ``("cpd00001", "e")``.

        Both spellings are accepted: the raw SBML form (``__91__c__93__``,
        ``M_`` prefix) and the decoded form cobra uses in memory (``[c]``, no
        prefix).  Returns ``(decoded_id, "")`` when no known compartment suffix
        matches, which lets callers decide whether to skip or raise.
        """
        text = decode(self.strip_species_prefix(species_id))
        for comp, suffixes in self.suffixes.items():
            for suffix in sorted(suffixes, key=len, reverse=True):
                for candidate in {suffix, decode(suffix)}:
                    if text.endswith(candidate):
                        return text[: -len(candidate)], comp
        return text, ""

    def is_external(self, species_id: str) -> bool:
        return self.split_species(species_id)[1] == "e"

    def compartment_attr(self, comp: str = "c") -> str:
        """Compartment attribute the family writes into its SBML."""
        return COMPARTMENTS[self.family][comp]

    def species_id(self, compound: str, comp: str = "c") -> str:
        return SPECIES_PREFIX + str(compound) + self.suffixes[comp][0]

    def species_variants(self, compound: str, comp: str = "c") -> list[str]:
        """Every spelling a model may use for one compound/compartment."""
        encoded = self.species_id(compound, comp)
        plain = self.strip_species_prefix(encoded)
        return [encoded, plain, decode(encoded), decode(plain)]

    def cytosolic_partner(self, species_id: str, comp: str = "c") -> str:
        """Return the same compound in the ``comp`` compartment, or ``""``."""
        compound, current = self.split_species(species_id)
        if not compound or current == comp:
            return ""
        return self.species_id(compound, comp)

    # -- exchanges -------------------------------------------------------
    def is_exchange(self, reaction_id: str) -> bool:
        rid = str(reaction_id)
        if rid.startswith(REACTION_PREFIX):
            rid = rid[len(REACTION_PREFIX):]
        return rid.startswith("EX_")

    def exchange_compound(self, reaction_id: str) -> str:
        """``R_EX_cpd00001_e0`` -> ``cpd00001`` (family suffix removed)."""
        rid = str(reaction_id)
        if rid.startswith(REACTION_PREFIX):
            rid = rid[len(REACTION_PREFIX):]
        rid = re.sub(r"^EX_", "", rid)
        rid = decode(rid)
        candidates = [decode(suffix)
                      for suffixes in self.suffixes.values() for suffix in suffixes]
        for suffix in sorted(set(candidates), key=len, reverse=True):
            if rid.endswith(suffix):
                rid = rid[: -len(suffix)]
                break
        return rid

    def new_exchange_id(self, compound: str) -> str:
        """Exchange id written by the expand skill (kept stable per family)."""
        token = re.sub(r"[^A-Za-z0-9_]", "_", str(compound)) or "claim"
        ext = "e0" if self.family == "seed" else "e"
        return f"EX_{token}_{ext}"


#: canonical labels; every accepted ``--base`` value maps onto one adapter
BASE_ALIASES = {
    "agora": "agora", "agora2": "agora", "agora-2": "agora",
    "gapseq": "seed", "gap": "seed", "seed": "seed",
    "modelseed": "seed", "ms": "seed", "model-seed": "seed",
    "carveme": "bigg", "carve": "bigg", "bigg": "bigg",
}


def for_base(base: str) -> Namespace:
    """Build the adapter for a user-supplied base name."""
    key = str(base).strip().lower()
    family = BASE_ALIASES.get(key)
    if family is None:
        raise UnknownNamespace(
            f"unknown base {base!r}; choose one of: "
            + ", ".join(sorted(BASE_ALIASES)))
    return Namespace(name=key, family=family)


def detect(source: "str | Path") -> Namespace:
    """Infer the adapter from an SBML file, an SBML string, or a cobra model."""
    if hasattr(source, "metabolites"):  # cobra.Model
        species = [m.id for m in source.metabolites]
        text = ""
    else:
        text = str(source)
        if "<" not in text:
            candidate = Path(text)
            try:
                if candidate.is_file():
                    text = candidate.read_text(encoding="utf-8", errors="replace")
            except OSError:
                pass
        species = re.findall(r'<species\b[^>]*\bid="([^"]+)"', text)
    joined = " ".join(species)

    if "__91__" in joined or "__40__" in text:
        return Namespace(name="agora2", family="agora")
    if re.search(r'compartment="C_', text):
        return Namespace(name="carveme", family="bigg")
    if re.search(r"(?:_c0|_e0)(?:\"|$)", joined) or re.search(
            r'compartment="[ce]0"', text):
        return Namespace(name="gapseq", family="seed")
    raise UnknownNamespace(
        "cannot infer the reconstruction family; pass --gapseq/--carveme/"
        "--modelseed/--agora2 explicitly. "
        f"sampled species ids: {species[:3]}")


# -- helpers shared by the review skills ---------------------------------

def model_metabolite_id(model, species_id: str, ns: Namespace | None = None) -> str:
    """Map a raw SBML species id onto the in-memory cobra id.

    cobra strips the ``M_`` prefix when it reads an SBML it wrote itself, so
    ``M_a_e0`` in the file is ``a_e0`` on the model.  Both spellings resolve.
    """
    sid = str(species_id)
    if sid in model.metabolites:
        return sid
    stripped = ns.strip_species_prefix(sid) if ns else sid
    for candidate in (stripped, sid):
        if candidate in model.metabolites:
            return candidate
    target = norm(decode(stripped))
    for metabolite in model.metabolites:
        if norm(decode(metabolite.id)) == target:
            return metabolite.id
    return ""


def cytosolic_partner(model, species_id: str, ns: Namespace) -> str:
    """Resolve the cytosolic partner of an external metabolite in ``model``."""
    external = model_metabolite_id(model, species_id, ns)
    if not external:
        return ""
    compound, _ = ns.split_species(external)
    if not compound:
        return ""
    # the model may keep the M_ prefix, use decoded brackets, or both
    for candidate in ns.species_variants(compound, "c"):
        if candidate in model.metabolites:
            return candidate
    target = norm(decode(ns.strip_species_prefix(ns.species_id(compound, "c"))))
    for metabolite in model.metabolites:
        if norm(decode(ns.strip_species_prefix(metabolite.id))) == target:
            return metabolite.id
    return ""

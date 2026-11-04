"""Command-line interface for gemagent-core."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from . import __version__, evaluation, workflow
from .skills import SKILLS
from . import config


def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument('--db-dir', default='', help='GEMAgent-curated data directory')
    p.add_argument('--gapseq-bin', default='', help='official gapseq executable')
    p.add_argument('--rscript', default='', help='Rscript executable')
    p.add_argument('--eval-root', default='', help='optional evaluation root for legacy scoring and exchange-map discovery')
    p.add_argument('--force', action='store_true', help='rerun even when outputs exist')


def cmd_eval_us(args: argparse.Namespace) -> int:
    calls = pd.read_csv(args.calls, sep="\t", dtype=str)
    if args.derived_from_decisions:
        calls = evaluation.decision_to_predicted(calls, args.decision_column, args.applied_column)
    truth = pd.read_csv(args.truth, sep="\t", dtype=str)
    metrics, merged = evaluation.evaluate_frames(truth, calls,
                                                 predicted_column=args.predicted_column)
    print(json.dumps(metrics, indent=2))
    if args.out:
        out = Path(args.out)
        out.mkdir(parents=True, exist_ok=True)
        (out / "metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
        merged.to_csv(out / "evaluated_calls.tsv", sep="\t", index=False)
    return 0
def cmd_namespace(args: argparse.Namespace) -> int:
    """Inspect (or assert) the id namespace of a reconstructed model."""
    from . import namespace as nslib

    try:
        ns = nslib.for_base(args.base) if args.base else nslib.detect(args.model)
    except nslib.UnknownNamespace as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    import contextlib
    import io

    import cobra

    # cobra announces every boundary metabolite it turns into an exchange
    # reaction; keep that out of the report unless asked for.
    if args.verbose:
        model = cobra.io.read_sbml_model(args.model)
    else:
        with contextlib.redirect_stdout(io.StringIO()):
            model = cobra.io.read_sbml_model(args.model)
    exchanges = [r for r in model.reactions if ns.is_exchange(r.id)]
    unresolved: list[tuple[str, str]] = []
    for reaction in exchanges:
        metabolites = list(reaction.metabolites)
        if len(metabolites) != 1:
            unresolved.append((reaction.id, "not a single-metabolite exchange"))
            continue
        if not nslib.cytosolic_partner(model, metabolites[0].id, ns):
            unresolved.append((reaction.id, "no cytosolic partner"))
    unknown_comp = [m.id for m in model.metabolites if not ns.split_species(m.id)[1]]

    report = {
        "model": str(args.model),
        "base": ns.name,
        "family": ns.family,
        "reactions": len(model.reactions),
        "metabolites": len(model.metabolites),
        "exchanges": len(exchanges),
        "unresolved_exchanges": len(unresolved),
        "metabolites_without_known_compartment": len(unknown_comp),
    }
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(f"base: {ns.name} (family {ns.family})")
        print(f"  reactions={report['reactions']} metabolites={report['metabolites']} "
              f"exchanges={report['exchanges']}")
        print(f"  unresolved exchanges: {len(unresolved)}")
        print(f"  metabolites without a known compartment suffix: {len(unknown_comp)}")
        for reaction in exchanges[: args.samples]:
            metabolites = list(reaction.metabolites)
            if len(metabolites) != 1:
                continue
            external = nslib.model_metabolite_id(model, metabolites[0].id, ns)
            partner = nslib.cytosolic_partner(model, metabolites[0].id, ns)
            print(f"  {reaction.id} -> {ns.exchange_compound(reaction.id)} "
                  f"external={external} cytosolic={partner or '-'}")
        for rid, reason in unresolved[:5]:
            print(f"  ! {rid}: {reason}")
    if args.strict and unresolved:
        return 1
    return 0


def cmd_build(args: argparse.Namespace) -> int:
    result = workflow.build(
        args.genome, args.out, recipe=args.recipe,
        org_name=args.org_name, taxonomy=args.taxonomy,
        medium=args.medium, eval_medium=args.eval_medium,
        claims_dir=args.claims_dir, us_microbe_id=args.us_microbe_id,
        kofam_tsv=args.kofam_tsv, db_dir=args.db_dir,
        gapseq_bin=args.gapseq_bin, rscript=args.rscript,
        eval_root=args.eval_root, threads=args.threads, force=args.force,
        from_draft=args.from_draft, biomass_review=not args.no_biomass_review,
        gapfill_review=not args.no_gapfill_review,
        gpr_check=not args.no_gpr_check, us_review=not args.no_us_review,
        us_prune=not args.no_us_prune, dry_run=args.dry_run, gpr_min_or=args.gpr_min_or,
        gtdb_summary=args.gtdb_summary, mimedb_tsv=args.mimedb_tsv,
        exchange_map=args.exchange_map,
    )
    print(json.dumps(result, indent=2))
    return 0


def cmd_skill_list(_: argparse.Namespace) -> int:
    for name in sorted(SKILLS):
        spec = SKILLS[name]
        print(f"{name}\t{spec.description}\trequires={','.join(spec.requires)}")
    return 0


def cmd_skill_run(args: argparse.Namespace) -> int:
    result = workflow.run_skill(
        args.skill, args.workspace, org_name=args.org_name,
        taxonomy=args.taxonomy, genome_id=args.genome_id, medium=args.medium,
        kofam_tsv=args.kofam_tsv, claims_dir=args.claims_dir,
        microbe_id=args.microbe_id, db_dir=args.db_dir,
        eval_root=args.eval_root, force=args.force, gpr_min_or=args.gpr_min_or,
        gtdb_summary=args.gtdb_summary, mimedb_tsv=args.mimedb_tsv,
        features_only=args.features_only,
        exchange_map=args.exchange_map, base=args.base,
    )
    print(json.dumps(result, indent=2))
    return 0


def cmd_recipe_list(_: argparse.Namespace) -> int:
    root = Path(__file__).resolve().parent / 'recipes'
    for p in sorted(root.glob('*.json')):
        data = json.loads(p.read_text(encoding='utf-8'))
        print(f"{data.get('name', p.stem)}\t{data.get('description', '')}")
    return 0


def cmd_recipe_run(args: argparse.Namespace) -> int:
    return cmd_build(argparse.Namespace(
        genome=args.genome, out=args.out, recipe=args.recipe,
        org_name=args.org_name, taxonomy=args.taxonomy, medium=args.medium,
        eval_medium=args.eval_medium, claims_dir=args.claims_dir,
        us_microbe_id=args.us_microbe_id, kofam_tsv=args.kofam_tsv,
        db_dir=args.db_dir, gapseq_bin=args.gapseq_bin, rscript=args.rscript,
        eval_root=args.eval_root, threads=args.threads, force=args.force,
        from_draft=args.from_draft, no_biomass_review=args.no_biomass_review,
        no_gapfill_review=args.no_gapfill_review, no_gpr_check=args.no_gpr_check,
        no_us_review=args.no_us_review, no_us_prune=args.no_us_prune,
        dry_run=args.dry_run, gpr_min_or=args.gpr_min_or,
        gtdb_summary=args.gtdb_summary, mimedb_tsv=args.mimedb_tsv,
        exchange_map=args.exchange_map,
    ))


def parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog='gemagent', description='Official gapseq base + modular review skills')
    ap.add_argument('--version', action='version', version=f'gemagent {__version__}')
    sub = ap.add_subparsers(dest='command', required=True)

    b = sub.add_parser('build', help='build one GEM')
    b.add_argument('genome')
    b.add_argument('--out', required=True)
    b.add_argument('--recipe', default='full', choices=['full', 'base-only'])
    b.add_argument('--org-name', default='')
    b.add_argument('--taxonomy', default='')
    b.add_argument('--medium', default='', help='explicit gapfill medium CSV')
    b.add_argument('--eval-medium', default='', help='benchmark evaluation medium CSV')
    b.add_argument('--claims-dir', default='')
    b.add_argument('--us-microbe-id', default='')
    b.add_argument('--kofam-tsv', default='')
    b.add_argument('--threads', type=int, default=1)
    b.add_argument('--from-draft', action='store_true')
    b.add_argument('--gpr-min-or', type=int, default=2)
    b.add_argument('--no-biomass-review', action='store_true')
    b.add_argument('--no-gapfill-review', action='store_true')
    b.add_argument('--no-gpr-check', action='store_true')
    b.add_argument('--no-us-review', action='store_true')
    b.add_argument('--no-us-prune', action='store_true')
    b.add_argument('--exchange-map', default='', help='exchange_map.tsv or its data directory')
    b.add_argument('--gtdb-summary', default='')
    b.add_argument('--mimedb-tsv', default='')
    b.add_argument('--dry-run', action='store_true')
    _add_common(b)
    b.set_defaults(func=cmd_build)

    s = sub.add_parser('skill', help='list or run a review skill')
    ss = s.add_subparsers(dest='skill_command', required=True)
    sl = ss.add_parser('list')
    sl.set_defaults(func=cmd_skill_list)
    sr = ss.add_parser('run')
    sr.add_argument('skill', choices=sorted(SKILLS))
    sr.add_argument('--workspace', required=True)
    sr.add_argument('--genome-id', default='')
    sr.add_argument('--org-name', default='')
    sr.add_argument('--taxonomy', default='')
    sr.add_argument('--medium', default='')
    sr.add_argument('--kofam-tsv', default='')
    sr.add_argument('--claims-dir', default='')
    sr.add_argument('--microbe-id', default='')
    sr.add_argument('--gpr-min-or', type=int, default=2)
    sr.add_argument('--gtdb-summary', default='')
    sr.add_argument('--mimedb-tsv', default='')
    sr.add_argument('--features-only', action='store_true')
    sr.add_argument('--exchange-map', default='')
    skill_fam = sr.add_mutually_exclusive_group()
    skill_fam.add_argument('--base', default='',
                           help='model family: gapseq/carveme/modelseed/agora2')
    skill_fam.add_argument('--gapseq', dest='base', action='store_const', const='gapseq')
    skill_fam.add_argument('--carveme', dest='base', action='store_const', const='carveme')
    skill_fam.add_argument('--modelseed', dest='base', action='store_const', const='modelseed')
    skill_fam.add_argument('--agora2', dest='base', action='store_const', const='agora2')
    _add_common(sr)
    sr.set_defaults(func=cmd_skill_run)

    r = sub.add_parser('recipe', help='list or run a named recipe')
    rs = r.add_subparsers(dest='recipe_command', required=True)
    rl = rs.add_parser('list')
    rl.set_defaults(func=cmd_recipe_list)
    rr = rs.add_parser('run')
    rr.add_argument('recipe', choices=['full', 'base-only'])
    rr.add_argument('genome')
    rr.add_argument('--out', required=True)
    rr.add_argument('--org-name', default='')
    rr.add_argument('--taxonomy', default='')
    rr.add_argument('--medium', default='')
    rr.add_argument('--eval-medium', default='')
    rr.add_argument('--claims-dir', default='')
    rr.add_argument('--us-microbe-id', default='')
    rr.add_argument('--kofam-tsv', default='')
    rr.add_argument('--threads', type=int, default=1)
    rr.add_argument('--from-draft', action='store_true')
    rr.add_argument('--gpr-min-or', type=int, default=2)
    rr.add_argument('--no-biomass-review', action='store_true')
    rr.add_argument('--no-gapfill-review', action='store_true')
    rr.add_argument('--no-gpr-check', action='store_true')
    rr.add_argument('--no-us-review', action='store_true')
    rr.add_argument('--no-us-prune', action='store_true')
    rr.add_argument('--exchange-map', default='')
    rr.add_argument('--gtdb-summary', default='')
    rr.add_argument('--mimedb-tsv', default='')
    rr.add_argument('--dry-run', action='store_true')
    _add_common(rr)
    rr.set_defaults(func=cmd_recipe_run)
    n = sub.add_parser('namespace', help='inspect or assert a model id namespace')
    nn = n.add_subparsers(dest='namespace_command', required=True)
    ni = nn.add_parser('inspect', help='detect or assert the reconstruction family')
    ni.add_argument('model')
    fam = ni.add_mutually_exclusive_group()
    fam.add_argument('--base', default='', help='explicit base: gapseq/carveme/modelseed/agora2')
    fam.add_argument('--gapseq', dest='base', action='store_const', const='gapseq')
    fam.add_argument('--carveme', dest='base', action='store_const', const='carveme')
    fam.add_argument('--modelseed', dest='base', action='store_const', const='modelseed')
    fam.add_argument('--agora2', dest='base', action='store_const', const='agora2')
    ni.add_argument('--samples', type=int, default=5)
    ni.add_argument('--json', action='store_true')
    ni.add_argument('--verbose', action='store_true',
                    help='show cobra''s SBML loading messages')
    ni.add_argument('--strict', action='store_true',
                    help='exit 1 when an exchange cannot be resolved')
    ni.set_defaults(func=cmd_namespace)

    e = sub.add_parser('eval-us', help='evaluate final uptake/secretion calls against truth')
    e.add_argument('--truth', required=True, help='TSV with observed and exchange_id/direction')
    e.add_argument('--calls', required=True, help='TSV containing calls or prune decisions')
    e.add_argument('--predicted-column', default='predicted')
    e.add_argument('--derived-from-decisions', action='store_true')
    e.add_argument('--decision-column', default='decision')
    e.add_argument('--applied-column', default='applied')
    e.add_argument('--out', default='')
    e.set_defaults(func=cmd_eval_us)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    return args.func(args)


if __name__ == '__main__':
    raise SystemExit(main())

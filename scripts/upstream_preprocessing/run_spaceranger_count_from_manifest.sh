#!/usr/bin/env bash
# Run Space Ranger count commands from config/spaceranger_count_manifest_from_archived_scripts.csv.
# This is a portable wrapper around the archived server-specific scripts.
# Edit SPACERANGER_BIN, WORKDIR, and paths in the manifest before running on another system.

set -euo pipefail

MANIFEST=${1:-config/spaceranger_count_manifest_from_archived_scripts.csv}
WORKDIR=${WORKDIR:-results/spaceranger}
SPACERANGER_BIN=${SPACERANGER_BIN:-spaceranger}
THREADS_DEFAULT=${THREADS_DEFAULT:-16}
MEM_DEFAULT=${MEM_DEFAULT:-128}
RUN_IDS=${RUN_IDS:-}
DRY_RUN=${DRY_RUN:-true}

mkdir -p "${WORKDIR}"

python3 - "$MANIFEST" "$WORKDIR" "$SPACERANGER_BIN" "$THREADS_DEFAULT" "$MEM_DEFAULT" "$RUN_IDS" "$DRY_RUN" <<'PY'
import csv, os, shlex, subprocess, sys
manifest, workdir, spaceranger, threads_default, mem_default, run_ids, dry_run = sys.argv[1:]
selected = set([x.strip() for x in run_ids.split(',') if x.strip()])
with open(manifest, newline='') as f:
    for row in csv.DictReader(f):
        run_id = row.get('id','').strip()
        if not run_id:
            continue
        if selected and run_id not in selected:
            continue
        args = [spaceranger, 'count', f'--id={run_id}']
        for csv_key, cli_key in [
            ('transcriptome','transcriptome'), ('fastqs','fastqs'), ('sample','sample'),
            ('image','image'), ('slide','slide'), ('area','area'), ('loupe_alignment','loupe-alignment'),
            ('probe_set','probe-set'), ('filter_probes','filter-probes'), ('create_bam','create-bam')
        ]:
            val = row.get(csv_key, '').strip()
            if val:
                args.append(f'--{cli_key}={val}')
        args.append(f"--localcores={row.get('localcores','').strip() or threads_default}")
        args.append(f"--localmem={row.get('localmem','').strip() or mem_default}")
        if row.get('nosecondary','').strip().lower() == 'true':
            args.append('--nosecondary')
        cmd = ' '.join(shlex.quote(a) for a in args)
        print('\n#', run_id, 'from', row.get('source_script',''))
        print('cd', shlex.quote(workdir), '&&', cmd)
        if dry_run.lower() not in {'1','true','yes','y'}:
            subprocess.check_call(args, cwd=workdir)
PY

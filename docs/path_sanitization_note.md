# Path sanitization note

Private absolute server paths in archived provenance scripts and manifests have been replaced with generic placeholders before public release. These placeholders preserve the directory logic without exposing local usernames, server mount points, or institutional storage topology.

Common placeholders used in this repository include:

- `/path/to/local_home`
- `/path/to/local_storage/raid2/project`
- `/path/to/local_storage/raid3/project`
- `/path/to/GEO_or_local_project`
- `/path/to/collaborator_home`

To re-run archived upstream scripts, edit these placeholders to point to the local download directory containing the GEO-linked raw and intermediate files. Public-facing scripts in `scripts/` are written to use relative repository paths or configurable environment variables when possible.

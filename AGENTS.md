# ACCLIM source-of-truth and editing protocols

## Authoritative scripts

The following Drive files are the authoritative versions of the ACCLIM web-app scripts:

- `G:\My Drive\Model\ACCLIM\ACCLIM_GEE_web_app_script_Public.txt`
- `G:\My Drive\Model\ACCLIM\ACCLIM_GEE_web_app_script_Private.txt`

The local repository, workspace output files, and GitHub are synchronized copies, not substitutes for the current Drive versions. Do not change this rule without explicit user authorization.

## Before editing

1. Read the current authoritative Drive scripts and compare each with its corresponding local copy. Record content hashes of the Drive files to detect subsequent edits.
2. Preserve all user changes, including experimental parameter overrides such as `C4frac = ee.Number(0);`. ACCLIM is research code; do not restore model defaults or remove experiments unless requested.
3. Reconcile local differences against Drive before implementing the requested change. Preserve intentional differences between the public and private scripts.
4. If conflicting changes have unclear intent, ask the user before resolving them. If Drive cannot be read, do not assume a local copy is authoritative; report the limitation before proceeding with dependent edits.

## During editing

- Apply only the requested changes to the reconciled versions.
- Do not overwrite one script wholesale with the other to mirror an update. Apply the requested change to each while preserving its existing content and settings.
- Review the diff against the Drive baseline and run checks appropriate to the change.

## Before saving to Drive

1. Re-read and hash each Drive file immediately before replacing it. Compare with the baseline captured before editing.
2. If a file changed in the meantime, reconcile the new edits and repeat the review and checks before saving. Ask about ambiguous conflicts.
3. Save a uniquely dated backup of each current Drive file before replacing it. Do not overwrite existing backups. If backup creation fails, do not replace the original.
4. Write the reconciled result and verify the saved contents match the intended result.

Never overwrite Drive directly from a potentially stale local copy.

## Synchronization and reporting

- Update the repository and workspace mirrors from the reconciled Drive versions. Synchronize GitHub within the user's authorized scope.
- Keep public/private differences intact; do not publish private-only content to the public repository.
- State which files were updated, what was tested, and any remaining deployment step.
- Updating Drive or GitHub does not republish an Earth Engine app. Do not claim the live app is updated without verifying deployment.

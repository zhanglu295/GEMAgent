# Security policy

## Secrets

Do not commit API keys, access tokens, passwords, SSH keys, solver licenses, or
private data. GEMAgent reads credentials from environment variables or a local
`.env` file, which is excluded by `.gitignore`.

If a credential is committed, revoke it immediately and remove it from the
repository history before publishing the repository.

## Reporting a vulnerability

Do not open a public issue for a security-sensitive problem. Contact the
corresponding author at `ericluzhang@hkbu.edu.hk` with:

- a description of the issue;
- the affected version or commit;
- minimal reproduction steps;
- the potential impact.

Please avoid including active credentials or private research data in the
initial report.

# Contributing

Thanks for helping. Bug reports, receiver logs that do not decode, and pull requests are all welcome.

## Reporting a problem

Open an [issue](../../issues/new/choose). Say what you ran, what you expected and what happened, and
include:

- the `gnss --version` or `geomaticist_gnss_py.version()` output and your OS;
- the receiver make and model and the file type (RINEX version, T04, MDB, ...);
- the full error text or the JSON warnings.

If you can share the data file, attach it or link it. Only share data you are allowed to share. A short
trimmed piece (`gnss rinex trim`) is often enough.

Questions about how to use the engine belong in [Discussions](../../discussions). Security problems: see
[SECURITY.md](SECURITY.md), not a public issue.

## Pull requests

1. Open an issue first for anything larger than a small fix, so we can agree on the approach.
2. Keep each pull request to one change, with a test.
3. Run `cargo xtask ci` before you push: formatting, clippy (warnings are errors), tests, dependency
   rules and the licence check. CI runs the same.
4. Add a line to `CHANGELOG.md` under "Unreleased".

Rules of the code base: no `unsafe`, no `unwrap` outside tests, no network access from the engine (it
takes files and parsed products handed to it), and no new dependency without a reason; licences must pass
`cargo deny` (permissive only).

## Contributor License Agreement

Before your first pull request is merged you sign the [Contributor License Agreement](CLA.md). The CLA
Assistant bot asks you on the pull request; you sign by posting the sentence it gives you as a comment, once,
and that covers all your later contributions.

In short: you keep the copyright in your contribution. You give the project owner a permanent licence to
use, change and distribute it, including under other terms (for example a commercial licence) as well as
under this project's MIT OR Apache-2.0 licence. You confirm the work is your own, or that you have the
right to submit it, and you grant a patent licence for it. Your contribution is always published here
under MIT OR Apache-2.0, so everyone can use it on the same terms as the rest of the project.

## Licence

By contributing you agree that your contribution is licensed under MIT OR Apache-2.0, as the project is,
in addition to the terms of the CLA.

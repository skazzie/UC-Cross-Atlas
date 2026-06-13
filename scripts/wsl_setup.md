# WSL setup — human gate for the laptop MAGMA arm

**Why this exists.** The MAGMA binary at `data/reference/magma` is an
ELF 64-bit Linux x86-64 executable (verified: BuildID
`a49e0123dbe15bb65734ba3bc99d2d8f88cd1330`, GNU/Linux 3.2.0, MAGMA
v1.10). It can't run natively on Windows. WSL2 is installed
(`wsl --version` reports 2.7.3.0) but no Linux distribution is
present, so this is the one-time setup step that unblocks Track A
of the local handoff.

**This step is a hard human gate.** The handoff doc explicitly forbids
unattended `wsl --install` because Windows requires a reboot and
prompts for an Ubuntu username/password during first-launch. So,
please run this yourself.

## Step 1 — install Ubuntu under WSL

In **Windows PowerShell as Administrator**:

```powershell
wsl --install -d Ubuntu
```

This downloads Ubuntu 24.04 LTS (~500 MB), registers it with WSL,
and may prompt for a reboot. **Reboot if prompted.**

After reboot, Windows will auto-open an Ubuntu terminal asking for
a Linux username + password. Pick any (these are Linux-only; not
linked to Windows). Once you see the `$` prompt, Ubuntu is ready.

If `wsl --install -d Ubuntu` returns "Ubuntu is already registered,"
just run `wsl -d Ubuntu` to drop into it.

## Step 2 — verify MAGMA runs

From the Ubuntu prompt (you'll be at `~`, i.e. `/home/<username>`):

```bash
cd /mnt/c/Users/muska/UC-Cross-Atlas
chmod +x data/reference/magma
data/reference/magma --version
```

Expected output:

```
MAGMA, version 1.10 (custom build, 2022-01-10)
...
```

If `--version` prints the MAGMA banner, **Track A is unblocked**.

### Troubleshooting

- **`Permission denied`**: the `chmod +x` step is required because
  Windows-side file permissions don't flag the bit as executable.
  Re-run the `chmod` line.
- **`No such file or directory`** despite the path existing: the
  Linux loader may be missing libc6 on a minimal Ubuntu install. The
  MAGMA binary is statically linked (verified above) so this
  shouldn't fire, but if it does: `sudo apt update && sudo apt
  install libc6` and retry.
- **Wrong binary or corrupted**: replacement download URL is in
  `scripts/download_refs.sh` (the `vu.data.surf.nl` mirror that
  `_fetch reference/magma.zip` already pulled cleanly once).
- **Performance**: WSL2 file IO on `/mnt/c/...` is ~3-5× slower than
  native. For MAGMA's ~minute-scale runs this is fine; if you ever
  see it crawling, copy the LD ref + munged sumstats to
  `~/uc-cross-atlas-data/` on the WSL filesystem and re-run from
  there. Not needed for v1.

## Step 3 — confirm back in Windows

After step 2 prints the MAGMA banner, you can:

- Just leave the Ubuntu window open (it stays alive for further
  WSL commands).
- Or close it; `wsl -d Ubuntu` re-enters when needed.

Either way, Track A1-A3 in the handoff plan can now proceed.
Claude Code's next session can run the three MAGMA jobs (de Lange,
Yengo, SCZ) inside WSL via:

```powershell
wsl -d Ubuntu -- bash -c 'cd /mnt/c/Users/muska/UC-Cross-Atlas && export PATH=$PWD/data/reference:$PATH && bash code/01_magma/run_magma.sh ...'
```

(Claude will fill the exact `...` args when it actually runs.)

## What this WSL install does NOT enable

- LDSC (Track B) needs Python 2.7 in its own conda env under WSL.
  Separate setup; documented in `code/01b_ldsc/README.md` once that
  module is scaffolded.
- Anything touching the big atlases (TAURUS / HCA Gut / Pan-GI) —
  those stay HB-pinned regardless of WSL availability.

#!/usr/bin/env bash
# Native assembly comparison tools; no Docker and no raw-read assembly stack.
set -euo pipefail
venv_dir=${1:-${VIRTUAL_ENV:-"$HOME/isoscout-pi-venv"}}
if [[ $# -gt 1 || ! -f "$venv_dir/pyvenv.cfg" ]]; then
    echo 'Usage: bash install_native_alignment.sh EXISTING_VENV_DIRECTORY' >&2
    exit 2
fi
venv_dir=$(cd -- "$venv_dir" && pwd)
case "$(uname -s):$(uname -m)" in
    Linux:aarch64|Linux:arm64)
        platform=linux-aarch64; datasets_arch=arm64
        mamba_sha=7154d66112a623d8d930f8aab467347971d5c1403a84b1c25d19298f4fb9da03
        datasets_sha=2d817c9621bf3e5ea1e1f2c2d751ec4134ef4dd6c514bf37f8ab7b4476695b64 ;;
    Linux:x86_64)
        platform=linux-64; datasets_arch=amd64
        mamba_sha=ffc3cb8d52d4d6b354bdbb979c407719c485392b74e462cbd50811aa88e58f85
        datasets_sha=b569d60ca45d389776582da8529dda7ae5e013db2a16c7d2b19abedcab74a4fa ;;
    *) echo 'Supported: 64-bit Linux ARM64 or x86_64.' >&2; exit 2 ;;
esac
for tool in curl sha256sum unzip; do command -v "$tool" >/dev/null || { echo "Install $tool first." >&2; exit 1; }; done
root="$venv_dir/native-alignment"
mkdir -p "$root"
for tool in minimap2 paftools.js k8 datasets; do
    if [[ -e "$venv_dir/bin/$tool" || -L "$venv_dir/bin/$tool" ]]; then
        case "$(readlink "$venv_dir/bin/$tool" || true)" in
            "$root/"*) ;;
            *) echo "Refusing to replace an unmanaged tool: $venv_dir/bin/$tool" >&2; exit 1 ;;
        esac
    fi
done
work=$(mktemp -d "$root/download.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
curl --fail --location --retry 3 -o "$work/micromamba" \
    "https://github.com/mamba-org/micromamba-releases/releases/download/2.3.2-0/micromamba-$platform"
printf '%s  %s\n' "$mamba_sha" "$work/micromamba" | sha256sum -c -
install -m 755 "$work/micromamba" "$root/micromamba"
# Isolate native libraries from the Python virtual environment. Micromamba is
# used only during installation; no background service runs during benchmarks.
action=create
[[ ! -d "$root/env/conda-meta" ]] || action=install
"$root/micromamba" --no-rc --root-prefix "$root/mamba" "$action" --yes \
    --prefix "$root/env" --override-channels --strict-channel-priority \
    --channel conda-forge --channel bioconda 'minimap2=2.31' 'k8=1.2'
curl --fail --location --retry 3 -o "$work/datasets.zip" \
    "https://github.com/ncbi/datasets/releases/download/v18.35.0/linux-$datasets_arch.cli.package.zip"
printf '%s  %s\n' "$datasets_sha" "$work/datasets.zip" | sha256sum -c -
unzip -q "$work/datasets.zip" datasets -d "$work"
install -m 755 "$work/datasets" "$root/datasets"
for tool in minimap2 paftools.js k8; do
    test -x "$root/env/bin/$tool"
    ln -sfn "$root/env/bin/$tool" "$venv_dir/bin/$tool"
done
ln -sfn "$root/datasets" "$venv_dir/bin/datasets"
export PATH="$venv_dir/bin:$PATH"
minimap2 --version
k8 -e 'print("k8 works")'
datasets version
"$root/micromamba" --no-rc --root-prefix "$root/mamba" list --prefix "$root/env" --explicit > "$root/packages-explicit.txt"
printf 'micromamba 2.3.2-0\nminimap2 2.31\nk8 1.2\ndatasets 18.35.0\nplatform %s\n' "$platform" > "$root/versions.txt"
echo "Native assembly tools ready. Activate: source \"$venv_dir/bin/activate\""

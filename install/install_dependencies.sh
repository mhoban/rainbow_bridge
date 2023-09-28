#!/usr/bin/env bash

pkg_installed() {
  dpkg-query -W -f='${Status}' $1 2>/dev/null >/dev/null
}

confirm() {
  local msg=$1
  while true; do
    read -p "$msg (Y/N/c) " yn
    case $yn in
      [Yy]* ) return 0;;
      [Nn]* ) return 1;;
      [Cc]* ) exit;;
      * ) echo "Please answer y, n, or cancel.";;
    esac
  done
}
 

get_singularity() {
  latest_url="https://github.com/sylabs/singularity/releases/latest"
  url=$(curl -Ls -o /dev/null -w %{url_effective} $latest_url)
  os=$(lsb_release -c | awk '{print $2}') 
  version=$(basename $url)
  # https://github.com/sylabs/singularity/releases/download/v4.0.0/singularity-ce_4.0.0-jammy_amd64.deb
  pkg=$(printf "https://github.com/sylabs/singularity/releases/download/%s/singularity-ce_%s-%s_amd64.deb" $version ${version##v} $os)
  echo $pkg
}

get_nextflow() {
  latest_url="https://github.com/nextflow-io/nextflow/releases/latest"
  url=$(curl -Ls -o /dev/null -w %{url_effective} $latest_url)
  version=$(basename $url)
  # https://github.com/sylabs/singularity/releases/download/v4.0.0/singularity-ce_4.0.0-jammy_amd64.deb
  pkg=$(printf "https://github.com/nextflow-io/nextflow/releases/download/%s/nextflow-%s-all" $version ${version##v})
  echo $pkg
}

if [[ "$(whoami)" != "root" ]]; then
  echo "This will work better/at all if you run this script with sudo"
fi

tmpdir=$(mktemp -d)

if ! [ -x singularity ]; then
  # get singularity deb package
  sing=$(get_singularity)
  if code=$(curl -w %{http_code} -fsL "$sing" -o "$tmpdir/singularity.deb"); then
    echo "Installing singularity..."
    apt install "$tmpdir/singularity.deb"
  else
    echo "Singularity package failed to download: $code!"
    exit 1
  fi
fi

if ! pkg_installed "default-jre"; then
  if confirm "Java runtime isn't installed, install it?"; then
    apt install default-jre
  fi
fi

# get nextflow executable
if ! [ -x nextflow ]; then
  nf=$(get_nextflow)
  if code=$(curl -w %{http_code} -fsL "$nf" -o "$tmpdir/nextflow"); then
    echo "Installing nextflow"
    chmod 755 "$tmpdir/nextflow"
    mv "$tmpdir/nextflow" /usr/local/bin 
  else
    echo "Nextflow failed to download: $code!"
  fi
fi

rm -rf "$tmpdir"

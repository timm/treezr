#!/bin/bash
# ell - Essential shell config

# If executed (not sourced), start bash with this as init file
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    exec bash --init-file "${BASH_SOURCE[0]}" -i
fi

# Core setup
Here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
Xdir="$(basename "$(dirname "$Here")")"
export BASH_SILENCE_DEPRECATION_WARNING=1
export PATH="$Here:$PATH"

# Colors & prompt
bold=$(tput bold) col0=$(tput sgr0) col1=$(tput setaf 6) col2=$(tput setaf 3)
branch() { git branch 2>/dev/null | awk '/^\*/ {print $2}'; }
PROMPT_COMMAND='PS1="🌴 ${bold}${col1}${Xdir}${col0} ${col1}$(basename "$(dirname "$PWD")")/$(basename "$PWD")${col0} ${col2}$(branch)${col0} ▶ "'

# Essential aliases
alias ..='cd ..' c='clear' Q='exit'
alias l='ls -lh' la='ls -la' t='tree -L 1' ls="\ls --color"
alias gs='git status -sb' ga='git add .' gc='git commit -m' gp='git push' gl='git log --oneline --graph --decorate'
alias vi="nvim -u '$Here/nvim_init.lua'" h='history'
alias reload="source '$Here/ell' && echo ✅"

# Mega useful extras
alias grep='grep --color=auto' less='less -R'
alias ports='lsof -i -P -n | grep LISTEN'
alias myip='curl -s ipinfo.io/ip'
mkcd() { mkdir -p "$1" && cd "$1"; }

exit() { builtin exit "$@"; }

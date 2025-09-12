#!/usr/bin/env bash
src() { cat <<'EOF'
BEGIN {
  a.i.j = 42
  print a.i.j
  a.i=223
  array(a.i)
  print .123 }

#--------------------------------------------------------------------
EOF
}

prep() { gawk '
  BEGIN { print "# add 2 blank lines to fix line numbers (in errors)\n"  } 
        { print gensub(/\.([a-zA-Z_][a-zA-Z0-9_]*)/, "[\"\\1\"]", "g")}'; }

gawk -f <(src | prep) "$@"

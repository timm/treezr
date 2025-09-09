#!/usr/bin/env bash
src() { cat <<'EOF'
BEGIN {
  a.i.j = 42
  print a.i.j
  print a.123 }

EOF
}

prep() { gawk '
  BEGIN { print "# correct any error messages by adding 2 blank lines\n"  } 
        { print gensub(/\.([a-zA-Z_][a-zA-Z0-9_]*)/, "[\"\\1\"]", "g")}'; }

gawk -f <(src | prep) "$@"

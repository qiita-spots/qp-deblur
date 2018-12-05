import re
import string
import sys
from json import load


with open('support_files/sepp/tmpl_tiny-revnamemap.json') as f:
    revnamemap = load(f)


def relabel_newick(newick_string):
    pattern = re.compile("(UQrYOlnDN[^(,:)<>]+)")
    invalidChars = set(string.punctuation).union(set(string.whitespace))

    def replace_func(m):
        repl = m.group(1)
        if m.group(1) in revnamemap:
            repl = revnamemap[m.group(1)]
            if any(char in invalidChars for char in repl):
                repl = "'%s'" % repl
        else:
            repl = m.group(1)

        return repl
    t = pattern.sub(replace_func, newick_string)
    return t


for l in sys.stdin.readlines():
        sys.stdout.write(relabel_newick(l))

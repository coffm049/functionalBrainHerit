"""Validate that the qmd's R chunks are syntactically well formed.

Extracts each ```{r ...} chunk and parses it. Knitr inline expressions are
skipped (they are not standalone R). Prints the offending chunk if any fails.
"""
import io, re, sys

qmd = 'summary/06_results_catalog.qmd'
t = io.open(qmd, encoding='utf-8').read()

FENCE = '```'
pat = re.compile(re.escape(FENCE) + r'\{r[^\n]*\n(.*?)\n' + re.escape(FENCE),
                 re.DOTALL)
chunks = pat.findall(t)
print('found %d R chunks' % len(chunks))

# Minimal R parser via a subprocess-free route is unreliable; use Rscript if
# available, else fall back to a brace/paren balance check.
try:
    import subprocess
    have_r = subprocess.run(['Rscript', '--version'], capture_output=True)
    have_r = have_r.returncode == 0
except Exception:
    have_r = False
print('Rscript available:', have_r)

bad = 0
for i, c in enumerate(chunks):
    open('_chunk.R', 'w', encoding='utf-8').write(c)
    if have_r:
        import subprocess
        r = subprocess.run(['Rscript', '-e', 'invisible(parse("_chunk.R"))'],
                           capture_output=True, text=True)
        if r.returncode != 0:
            bad += 1
            print('--- chunk %d FAILED to parse ---' % i)
            print(c[:400])
            print(r.stderr[:400])
    else:
        # balance check
        for op, cl in [('{', '}'), ('(', ')'), ('[', ']')]:
            if c.count(op) != c.count(cl):
                bad += 1
                print('--- chunk %d unbalanced %s%s (%d vs %d) ---'
                      % (i, op, cl, c.count(op), c.count(cl)))
                break

import os
if os.path.exists('_chunk.R'):
    os.remove('_chunk.R')

print('BAD chunks:', bad)
sys.exit(1 if bad else 0)
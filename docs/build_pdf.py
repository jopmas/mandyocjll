"""Builds bugs_mandyoc_cpp.pdf from the .md, using only Python and xelatex: python3 docs/build_pdf.py"""
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

src = Path(__file__).parent / "bugs_mandyoc_cpp.md"
_, head, body = src.read_text().split("---\n", 2)
meta = dict(re.findall(r'(?m)^(\w+): "?(.*?)"?$', head))
SPECIAL = {"\\": r"\textbackslash{}", "{": r"\{", "}": r"\}", "_": r"\_", "%": r"\%", "&": r"\&",
           "#": r"\#", "$": r"\$", "^": r"\^{}", "~": r"\~{}"}


def esc(s):
    return "".join(SPECIAL.get(c, c) for c in s)


def code(s):  # allow line breaks in long code
    return r"\texttt{" + re.sub(r"(\\_|/|:|\(|,)", r"\1\\allowbreak{}", esc(s)) + "}"


def inline(s):
    out = []
    for p in re.split(r"(`[^`]*`|\$[^$]+\$|\*\*.+?\*\*|\*[^*]+\*)", s):
        if p.startswith("`") and len(p) > 1:
            out.append(code(p[1:-1]))
        elif p.startswith("$") and len(p) > 2:
            out.append(p)
        elif p.startswith("**") and len(p) > 4:
            out.append(r"\textbf{" + inline(p[2:-2]) + "}")
        elif p.startswith("*") and len(p) > 2:
            out.append(r"\textit{" + inline(p[1:-1]) + "}")
        else:
            out.append(esc(p).replace('"', "''"))
    return "".join(out)


L = body.split("\n")
o = []
i = 0
while i < len(L):
    l = L[i]
    if l.startswith("```"):
        j = i + 1
        while not L[j].startswith("```"):
            j += 1
        o += [r"\begin{Verbatim}[fontsize=\footnotesize,frame=leftline,framesep=3mm,samepage=true]"] + L[i + 1:j] + [r"\end{Verbatim}"]
        i = j
    elif l.startswith("# "):
        o.append(r"\section*{" + inline(l[2:]) + "}")
    elif l.startswith("|"):
        rows = []
        while i < len(L) and L[i].startswith("|"):
            rows.append([c.strip() for c in L[i].strip().strip("|").split("|")])
            i += 1
        i -= 1
        sep = rows.pop(1)
        k = len(sep)
        w = [max(3, len(c)) for c in sep]
        if sum(w) < 60:  # numeric table: width from the content
            w = [max(6, min(26, max(len(r[n]) for r in rows))) for n in range(k)]
        w = [x / sum(w) * (0.98 - 0.025 * k) for x in w]
        spec = "".join(r">{\raggedright\arraybackslash}p{%.3f\linewidth}" % x for x in w)
        o.append(r"\needspace{%d\baselineskip}" % min(len(rows) + 3, 12))
        o.append(r"{\small\begin{longtable}{" + spec + r"}\toprule")
        o.append(" & ".join(r"\textbf{" + inline(c) + "}" for c in rows[0]) + r"\\\midrule\endhead\bottomrule\endlastfoot")
        o += [" & ".join(inline(c) for c in r) + r"\\" for r in rows[1:]] + [r"\end{longtable}}"]
    elif re.match(r"(\d+\.|-) ", l):
        env = "itemize" if l[0] == "-" else "enumerate"
        o.append(r"\begin{%s}" % env)
        while i < len(L) and re.match(r"(\d+\.|-) ", L[i]):
            o.append(r"\item " + inline(re.sub(r"^(\d+\.|-) ", "", L[i])))
            i += 1
        o.append(r"\end{%s}" % env)
        i -= 1
    elif l.startswith("\\newpage"):
        o.append(l)
    else:
        o.append(inline(l))
    i += 1

tex = r"""\documentclass[10pt]{article}
\usepackage[a4paper,margin=2.2cm]{geometry}
\usepackage{fontspec,longtable,booktabs,array,fancyvrb,titlesec,needspace}
\setmainfont{Helvetica Neue}\setmonofont[Scale=0.82]{Menlo}
\setlength{\parindent}{0pt}\setlength{\parskip}{5pt}\raggedbottom\emergencystretch=3em
\titleformat{\section}{\needspace{6\baselineskip}\large\bfseries}{}{0pt}{}
\begin{document}
\begin{center}{\LARGE @TITLE@}\\[4pt]{\large @SUB@}\\[4pt]@DATE@\end{center}
@BODY@
\end{document}
"""
for key, val in (("@TITLE@", esc(meta["title"])), ("@SUB@", esc(meta["subtitle"])), ("@DATE@", esc(meta["date"])), ("@BODY@", "\n".join(o))):
    tex = tex.replace(key, val)
with tempfile.TemporaryDirectory() as d:
    Path(d, "a.tex").write_text(tex)
    for _ in range(2):
        r = subprocess.run(["xelatex", "-interaction=nonstopmode", "a.tex"], cwd=d, capture_output=True, text=True)
    if not Path(d, "a.pdf").exists():
        raise SystemExit(r.stdout[-2000:])
    for line in r.stdout.splitlines():
        if line.startswith("!") or "Missing character" in line:
            print(line)
    shutil.copy(Path(d, "a.pdf"), src.with_suffix(".pdf"))
print("ok")

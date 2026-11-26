'''Relatório Zarr x ASCII (Mandyoc): tabelas CSV + PDF.'''
import csv, json, os, sys
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

S = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
L = json.load(open(f"{S}/local_bench.json"))

ASCII_C, ZARR_C = "#eb6834", "#2a78d6"          # categorical slots 2 / 1 (dataviz reference palette)
INK, INK2, GRID, SURF = "#1f1f1e", "#5f5e58", "#e4e3de", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "text.color": INK, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.edgecolor": GRID, "figure.facecolor": "white"})

VARS = ["density", "viscosity", "pressure", "strain", "strain_rate", "temperature", "heat", "velocity"]

# ---------------- craton_colision_lwr (hypatia, 96 ranks, 2501 x 351), saída 3300 ----------------
lwr_t = {"strain": 0.8286, "temperature": 2.3716, "heat": 0.7549, "viscosity": 0.9103, "velocity": 1.8096,
         "strain_rate": 0.8615, "pressure": 0.9181, "density": 0.9186}          # médias de 68 saídas (FD.out)
lwr_t_other = {"particles": 0.11, "lithology": 0.18}                          # sequência de mtimes, média de 67
lwr_ascii = {"density": 5167444, "heat": 2981252, "pressure": 10252239, "strain": 7645844, "strain_rate": 10438102,
             "temperature": 6684006, "velocity": 21809258, "viscosity": 9817299, "surface": 2600026,
             "particles": 39269220, "lithology": 364397012}
lwr_zarr = {"density": 93379, "pressure": 791666, "temperature": 1372376, "viscosity": 4415772, "strain": 5372759,
            "strain_rate": 6636355, "heat": 172256, "velocity": 13257862, "particles": 6912788,
            "lithology": 2825689, "surface": 391634}
ORDER = VARS + ["surface", "particles", "lithology"]
NAMES = {"density": "density", "viscosity": "viscosity", "pressure": "pressure", "strain": "strain",
         "strain_rate": "strain_rate", "temperature": "temperature", "heat": "heat", "velocity": "velocity (vx, vz)",
         "surface": "surface", "particles": "partículas (step_*)", "lithology": "litologia (litho_*)"}
OUT_SO_FAR, MYR_SO_FAR, HOURS_SO_FAR, TMAX = 68, 26.2, 28.5, 375.0
N_FULL = round(OUT_SO_FAR * TMAX / MYR_SO_FAR)

mb = lambda b: b / 1e6
rows1 = []
for v in ORDER:
    a, z = lwr_ascii[v], lwr_zarr[v]
    t = lwr_t.get(v, lwr_t_other.get(v, ""))
    rows1.append([NAMES[v], t, round(mb(a), 2), round(mb(z), 2), round(a / z, 1)])
ta, tz = sum(lwr_ascii.values()), sum(lwr_zarr.values())
t_ascii_out = sum(lwr_t.values()) + sum(lwr_t_other.values())
rows1.append(["TOTAL por saída", round(t_ascii_out, 2), round(mb(ta), 1), round(mb(tz), 1), round(ta / tz, 1)])

# ---------------- continental_rift (local, 2 ranks, 801 x 151), 2 saídas ----------------
ta_l, tz_l = L["t_ref_ascii"], L["t_var_zarr2_timeline"]
rows2 = []
for v in VARS + ["surface", "particles", "lithology"]:
    a, z = L["pervar"][v]
    rows2.append([NAMES[v], ta_l.get(v, ""), round(tz_l[v], 4), round(ta_l[v] / tz_l[v], 1) if v in ta_l else "",
                  round(mb(a), 3), round(mb(z), 3), round(a / z, 1)])
sa = sum(a for a, _ in L["pervar"].values()); sz_ = sum(z for _, z in L["pervar"].values())
fa = sum(ta_l[v] for v in VARS); fz = sum(tz_l[v] for v in VARS)
rows2.append(["TOTAL (tempo: campos)", round(fa, 3), round(fz, 3), round(fa / fz, 1),
              round(mb(sa), 2), round(mb(sz_), 2), round(sa / sz_, 1)])

# ---------------- formatos (local, saída do passo 1, layout por passo) ----------------
FMT = [("ascii", "ASCII (PETSc ASCII, padrão)"), ("bin", "PETSc binário (binary_output)"),
       ("h5", "HDF5 sem compressão"), ("zn0", "Zarr nativo, sem compressão"), ("zn1", "Zarr nativo, zlib 1 (adotado)"),
       ("zn5", "Zarr nativo, zlib 5"), ("zz1", "Zarr via z5 (lib externa), zlib 1")]
rows3 = []
base_t = sum(v for k, v in L["ascii"]["times"].items() if k in VARS)
base_s = L["ascii"]["size"]
for k, name in FMT:
    t = sum(v for kk, v in L[k]["times"].items() if kk in VARS)
    s = L[k]["size"]
    rows3.append([name, round(t, 4), round(base_t / t, 1), round(mb(s), 2), round(base_s / s, 1)])

H1 = ["variável", "tempo ASCII por saída (s)", "ASCII por saída (MB)", "Zarr zlib 1 estimado (MB)", "redução (×)"]
H2 = ["variável", "tempo ASCII (s)", "tempo Zarr (s)", "ganho (×)", "ASCII 2 saídas (MB)", "Zarr 2 saídas (MB)",
      "redução (×)"]
H3 = ["formato", "tempo dos campos (s)", "ganho vs ASCII (×)", "tamanho da saída (MB)", "redução vs ASCII (×)"]
for fn, h, r in [("1_craton_colision_lwr_hypatia.csv", H1, rows1), ("2_continental_rift_local.csv", H2, rows2),
                 ("3_formatos_local.csv", H3, rows3)]:
    with open(f"{OUT}/{fn}", "w", newline="") as f:
        w = csv.writer(f); w.writerow(h); w.writerows(r)


# ---------------- PDF ----------------
def fmt(x):
    if isinstance(x, float):
        return f"{x:.4f}".rstrip("0").rstrip(".") if x < 1 else f"{x:,.2f}".rstrip("0").rstrip(".")
    return str(x)


def wrap_h(h):
    import textwrap
    return "\n".join(textwrap.wrap(h, 16))


def table(ax, header, rows, col_w, bold_last=True, fs=8):
    ax.axis("off")
    cells = [[fmt(c) if c != "" else "–" for c in r] for r in rows]
    header = [wrap_h(h) for h in header]
    tb = ax.table(cellText=cells, colLabels=header, colWidths=col_w, loc="upper center", cellLoc="right",
                  colLoc="center")
    tb.auto_set_font_size(False); tb.set_fontsize(fs); tb.scale(1, 1.45)
    for (i, j), c in tb.get_celld().items():
        c.set_edgecolor(GRID); c.set_linewidth(0.6)
        if j == 0:
            c.set_text_props(ha="left"); c._loc = "left"
        if i == 0:
            c.set_facecolor("#f1f0ec"); c.set_text_props(weight="bold", color=INK, wrap=True); c.set_height(0.11)
        elif bold_last and i == len(rows):
            c.set_text_props(weight="bold")
    return tb


def page(pdf, title, subtitle=None):
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.text(0.07, 0.955, title, fontsize=15, weight="bold")
    if subtitle:
        fig.text(0.07, 0.932, subtitle, fontsize=9, color=INK2)
    return fig


def para(fig, y, text, fs=9, color=INK, x=0.07, w=0.86):
    import textwrap
    lines = []
    for p in text.split("\n"):
        lines += textwrap.wrap(p, int(w * 118 * 9 / fs)) or [""]
    fig.text(x, y, "\n".join(lines), fontsize=fs, color=color, va="top", linespacing=1.45)
    return y - len(lines) * fs * 1.45 / 72 / 11.69 - 0.012


def hbars(ax, labels, a, z, xlabel, log=False):
    import numpy as np
    y = np.arange(len(labels))[::-1]
    h = 0.36
    ax.barh(y + h / 2 + 0.01, a, h, color=ASCII_C, label="ASCII")
    ax.barh(y - h / 2 - 0.01, z, h, color=ZARR_C, label="Zarr zlib 1")
    ax.set_yticks(y, labels)
    if log:
        ax.set_xscale("log")
    ax.set_xlabel(xlabel)
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=8)


with PdfPages(f"{OUT}/comparacao_zarr_ascii.pdf") as pdf:
    # --- página 1: resumo + lwr
    fig = page(pdf, "Saídas do Mandyoc: Zarr x ASCII",
               "mandyoc_jpm, branch zarr_outputs (ab1561a) · zarr_output = True, zarr_compression = 1 · 26/11/2026")
    y = 0.905
    fig.text(0.07, y, "Resumo", fontsize=11, weight="bold"); y -= 0.02
    y = para(fig, y,
             f"• Tamanho: uma saída do craton_colision_lwr ocupa {mb(ta):.0f} MB em ASCII e ~{mb(tz):.0f} MB em Zarr "
             f"zlib 1 ({ta / tz:.0f}× menor). Até agora o lwr já gravou ~31 GB; o cenário completo "
             f"(~{N_FULL} saídas até {TMAX:.0f} Myr) chegaria a ~{mb(ta) * N_FULL / 1e3:.0f} GB em ASCII contra "
             f"~{mb(tz) * N_FULL / 1e3:.0f} GB em Zarr. 76% do ASCII é a litologia (litho_*).\n"
             f"• Tempo de escrita: no teste local, os campos levam {fa:.2f} s por saída em ASCII e {fz:.2f} s em Zarr "
             f"({fa / fz:.0f}× mais rápido); a saída Zarr inteira, com partículas, litologia e superfície, termina em "
             f"~0.16 s. No lwr, o ASCII gasta {t_ascii_out:.1f} s por saída, quase tudo nos campos, que passam pelo "
             f"rank 0 no viewer ASCII do PETSc.\n"
             f"• Impacto no tempo total: no lwr, a escrita ASCII soma {OUT_SO_FAR * t_ascii_out / 60:.0f} min em "
             f"{HOURS_SO_FAR:.1f} h de rodada ({100 * OUT_SO_FAR * t_ascii_out / 3600 / HOURS_SO_FAR:.1f}%). "
             f"O ganho em tempo de parede fica perto disso; o ganho grande é em disco e no pós-processamento.\n"
             f"• Pendente: o tempo do Zarr no hypatia. O job 336491 (craton_colision_zarr) ainda está na fila; os "
             f"tamanhos Zarr do lwr abaixo são estimados, não medidos.", fs=9)
    y -= 0.01
    fig.text(0.07, y, "1. craton_colision_lwr (hypatia, 96 ranks, malha 2501 × 351), por saída", fontsize=11,
             weight="bold"); y -= 0.008
    ax = fig.add_axes([0.07, y - 0.33, 0.86, 0.32])
    table(ax, H1, rows1, [0.30, 0.18, 0.17, 0.2, 0.15])
    y -= 0.345
    para(fig, y,
         "Tempo ASCII: campos = média das 68 saídas no FD.out; partículas e litologia = intervalo entre mtimes "
         "consecutivos dos arquivos (média de 67 saídas), porque o Mandyoc não as cronometra; a superfície não foi "
         "medida separadamente. ASCII: tamanhos reais da saída 3300. Zarr: os arquivos ASCII da saída 3300 "
         "comprimidos com zlib 1, nos tipos e chunks do zarr_output.cpp. Como os valores vêm com a precisão do texto "
         "(6 dígitos), os campos em float64 (viscosity, strain, strain_rate, heat, velocity) são um limite inferior; "
         "a litologia não tem o preenchimento de células vazias.", fs=7.5, color=INK2)
    pdf.savefig(fig); fig.savefig(f'{S}/pv/p{pdf.get_pagecount()}.png', dpi=110); plt.close(fig)

    # --- página 2: gráficos
    fig = page(pdf, "Tamanho e tempo por variável")
    ax = fig.add_axes([0.25, 0.53, 0.68, 0.36])
    hbars(ax, [NAMES[v] for v in ORDER], [mb(lwr_ascii[v]) for v in ORDER], [mb(lwr_zarr[v]) for v in ORDER],
          "MB por saída (escala log)", log=True)
    ax.set_title("craton_colision_lwr: tamanho de uma saída (Zarr estimado)", loc="left", fontsize=10,
                 weight="bold", pad=26)
    ax = fig.add_axes([0.25, 0.09, 0.68, 0.34])
    hbars(ax, [NAMES[v] for v in VARS], [ta_l[v] for v in VARS], [tz_l[v] for v in VARS],
          "segundos por saída (escala log)", log=True)
    ax.set_title("continental_rift (local): tempo de escrita de cada campo", loc="left", fontsize=10,
                 weight="bold", pad=26)
    pdf.savefig(fig); fig.savefig(f'{S}/pv/p{pdf.get_pagecount()}.png', dpi=110); plt.close(fig)

    # --- página 3: local
    fig = page(pdf, "2. continental_rift (local, 2 ranks, malha 801 × 151)",
               "mesmos inputs e parâmetros; ASCII com o mandyoc original, Zarr com o zarr_outputs; passos 0 e 1")
    ax = fig.add_axes([0.05, 0.60, 0.9, 0.31])
    table(ax, H2, rows2, [0.25, 0.13, 0.13, 0.1, 0.13, 0.13, 0.12], fs=7.5)
    y = para(fig, 0.585,
             "Saída do passo 1. ASCII: tempos do log (o viewer ASCII é coletivo, então o rank 0 mede tudo). "
             "Zarr: intervalo entre os mtimes de chunks consecutivos (temperature, o primeiro, vem do log), porque "
             "cada store é gravado por um rank diferente e o tempo impresso pelo rank 0 não inclui a compressão e a "
             "gravação feitas pelos outros. A saída Zarr completa (campos, partículas, litologia e superfície) leva "
             "~0.16 s. O ASCII não cronometra partículas, litologia e superfície. Tamanhos: soma das duas saídas (no Zarr, cada store guarda todas as "
             "saídas, incluindo eixos e metadados).", fs=7.5, color=INK2)
    fig.text(0.07, y - 0.01, "3. Formatos testados (local, saída do passo 1, campos + partículas)", fontsize=11,
             weight="bold")
    ax = fig.add_axes([0.07, y - 0.26, 0.86, 0.23])
    table(ax, H3, rows3, [0.36, 0.16, 0.15, 0.17, 0.16], bold_last=False)
    para(fig, y - 0.235,
         "Comparação feita antes da escolha do layout, com um mandyoc_<passo>.zarr por saída; a compressão é a "
         "mesma do layout atual por variável. O binário do PETSc é o mais rápido, mas não comprime e o tapioca não o "
         "lê. Tempos impressos pelo rank 0 (no Zarr e no HDF5, subestimam quando grava outro rank). zlib 5 custa ~35% a mais de tempo que zlib 1 para quase o mesmo tamanho, e o z5 empata com o "
         "nativo, por isso ficou o nativo com zlib 1. O ASCII inclui as partículas (step_*), mas não a litologia.",
         fs=7.5, color=INK2)
    pdf.savefig(fig); fig.savefig(f'{S}/pv/p{pdf.get_pagecount()}.png', dpi=110); plt.close(fig)

print(OUT)

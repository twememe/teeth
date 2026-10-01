#!/usr/bin/env python3
"""Generate dependency-free SVG figures for the Phase 2 release."""
import csv, json, html, math, argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument("--out", type=Path, default=ROOT/"results/figures");ap.add_argument("--mechanism", type=Path, default=ROOT/"results/16_mechanism");ap.add_argument("--project-root", type=Path, default=ROOT);args=ap.parse_args()
OUT = args.out
OUT.mkdir(parents=True, exist_ok=True)
W, H = 1000, 620

def esc(x): return html.escape(str(x))
def line(x1,y1,x2,y2,c="#cbd5e1",w=1): return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{c}" stroke-width="{w}"/>'
def txt(x,y,s,cls="label",anchor="start"): return f'<text x="{x}" y="{y}" class="{cls}" text-anchor="{anchor}">{esc(s)}</text>'
def doc(title, subtitle, body):
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial,sans-serif;fill:#1f2937}}.title{{font-size:25px;font-weight:700}}.sub{{font-size:14px;fill:#64748b}}.label{{font-size:13px}}.tick{{font-size:12px;fill:#64748b}}.legend{{font-size:13px}}</style>{txt(55,42,title,"title")}{txt(55,67,subtitle,"sub")}{body}</svg>'''
def save(name, body, title, subtitle): (OUT/name).write_text(doc(title,subtitle,body), encoding="utf-8")

# 1. Formal-MSA confidence ranking
rows = list(csv.DictReader((args.project_root/"results/14_quality/model_ranking_formal_msa.csv").open()))
rows.sort(key=lambda r: float(r["confidence_score"]), reverse=True); rows = rows[:15]
x0,x1,y0,y1=140,940,110,510
values=[float(r["confidence_score"]) for r in rows]
padding=max((max(values)-min(values))*.08,.002)
lo,hi=min(values)-padding,max(values)+padding
sx=lambda v: x0+(float(v)-lo)/(hi-lo)*(x1-x0)
body=line(x0,y0,x0,y1)+line(x0,y1,x1,y1)
for t in [lo+i*(hi-lo)/5 for i in range(6)]: body+=line(sx(t),y0,sx(t),y1)+txt(sx(t),y1+24,f"{t:.3f}","tick","middle")
for i,r in enumerate(rows):
    y=y0+i*(y1-y0)/max(1,len(rows)-1); col="#2563eb" if r["route"]=="native_blind" else "#f59e0b"; x=sx(r["confidence_score"])
    body+=line(x0,y,x1,y,"#f1f5f9")+f'<circle cx="{x}" cy="{y}" r="6" fill="{col}"/>'+txt(x0-10,y+4,f"{r['route'].replace('_blind','')} s{r['seed']}","tick","end")
body+=txt(530,585,"Confidence score","label","middle")+f'<circle cx="690" cy="565" r="6" fill="#2563eb"/>'+txt(705,570,"Native","legend")+f'<circle cx="800" cy="565" r="6" fill="#f59e0b"/>'+txt(815,570,"Immunogen","legend")
save("01_model_confidence_ranking.svg",body,"Formal-MSA model confidence ranking","Blue = Native; orange = Immunogen")

# 2. Route metric comparison from route summary
rows=list(csv.DictReader((args.project_root/"results/17_delivery/route_ensemble_summary.csv").open()))
body=""; metrics=[("iptm_mean","mean ipTM",0,1),("complex_plddt_mean","complex pLDDT",0,1),("hotspot_coverage_mean","hotspot coverage",0,1)]
for j,(key,label,a,b) in enumerate(metrics):
    y=150+j*115; body+=txt(180,y+6,label,"label","end")
    for i,r in enumerate(rows):
        val=float(r[key]); width=max(0,(val-a)/(b-a)*240); x=210+i*330; col="#2563eb" if "native" in r.get("route","") else "#f59e0b"
        body+=f'<rect x="{x}" y="{y-22}" width="{width:.1f}" height="28" fill="{col}"/>'+txt(x+width+8,y-2,f"{val:.3f}","tick")
body+=txt(500,445,"All bar lengths use the same 0–1 scale; values are dimensionless.","tick","middle")+f'<rect x="270" y="487" width="18" height="18" fill="#2563eb"/>'+txt(300,501,"Native","legend")+f'<rect x="465" y="487" width="18" height="18" fill="#f59e0b"/>'+txt(495,501,"Immunogen","legend")
save("02_native_immunogen_metrics.svg",body,"Native and Immunogen ensemble metrics","Empty-MSA formal ensembles; 15 seeds per construct; separate from the formal-MSA ranking")

# 3. ELISA curve
el=json.loads((args.mechanism/"ELISA_external_consistency.json").read_text())["rows"]; x0,x1,y0,y1=120,920,120,500
lx=lambda v: 0 if float(v)==0 else 1+(math.log2(float(v))+5)*.14
ex=lambda v:x0+v/2.1*(x1-x0); ey=lambda v:y1-float(v)/100*(y1-y0)
body=line(x0,y0,x0,y1)+line(x0,y1,x1,y1)
for v in [0,20,40,60,80,100]: body+=line(x0,ey(v),x1,ey(v))+txt(x0-10,ey(v)+4,v,"tick","end")
body+=txt(120,100,"Inhibition (%)","label")
for r in el:
    value=float(r['抗体浓度 (μg/mL)']);x=ex(lx(value))
    body+=line(x,y1,x,y1+5)+txt(x,y1+24,f"{value:g}","tick","middle")
body+=line(160,555,190,555,"#2563eb",3)+txt(198,559,"1A6-13-8","legend")+line(325,555,355,555,"#94a3b8",3)+txt(363,559,"Isotype IgG","legend")
for key,col in [("1A6-13-8 抑制率 (%)","#2563eb"),("同亚型 IgG 抑制率 (%)","#94a3b8")]:
    pts=" ".join(f"{ex(lx(r['抗体浓度 (μg/mL)'])):.1f},{ey(r[key]):.1f}" for r in el); body+=f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="3"/>'
    for r in el: body+=f'<circle cx="{ex(lx(r["抗体浓度 (μg/mL)"])):.1f}" cy="{ey(r[key]):.1f}" r="5" fill="{col}"/>'
body+=txt(520,597,"Antibody concentration (μg/mL; positive doses log-scaled; zero shown separately)","label","middle")
save("03_elisa_inhibition_curve.svg",body,"ELISA inhibition is concentration-dependent","Raw supplied data; assay metadata remain incomplete")

# 4. 6DF3 overlap
#    LOGIC-CLOSEOUT 20260930: read the corrected live result (results/16_mechanism/),
#    not the frozen bundle copy, so that this figure and the delivered table come from
#    one and the same computation. Superseded files remain only in the original project.
ov=list(csv.DictReader((args.mechanism/"6DF3_model_overlap.csv").open())); body=line(260,110,260,500)+line(260,500,930,500)
for i,r in enumerate(ov):
    y=150+i*65; val=float(r["jaccard"]); width=val/.40*650; label=r["model"].replace(".pdb","").replace("_representative","")
    body+=txt(250,y,label[:38],"tick","end")+f'<rect x="260" y="{y-20}" width="{width:.1f}" height="28" fill="#7c3aed"/>'+txt(270+width,y-2,f"{val:.3f}","tick")
body+=txt(600,555,"Jaccard: modelled IL-24 epitope vs 6DF3 IL-24 receptor-footprint residues (0-0.40)","label","middle")
save("04_6df3_receptor_overlap.svg",body,"IL-24 epitope overlap with the 6DF3 receptor footprint","Both sets in project numbering 1-181; human 6DF3 vs mouse project reference; hypothesis only")

(OUT/"README.md").write_text("""# Phase 2 visualizations\n\nThese SVG figures are generated from release CSV/JSON files (figure 4 reads the corrected live result in results/16_mechanism/) and are descriptive computational summaries, not experimental validation.\n\n- `01_model_confidence_ranking.svg`: formal-MSA confidence ranking by route.\n- `02_native_immunogen_metrics.svg`: route-level metric comparison.\n- `03_elisa_inhibition_curve.svg`: supplied ELISA inhibition data.\n- `04_6df3_receptor_overlap.svg`: Jaccard between the modelled IL-24 epitope and the IL-24-side residues of the 6DF3 receptor interface, both in project numbering. Corrected 2026-09-30 (logic closeout): protein contacts only, IL-24 side only, alignment-based numbering transfer.\n""",encoding="utf-8")
print("Wrote 4 SVG figures")

"""HTML pages rendered entirely by Python, with no templates or JavaScript."""
from html import escape
import json
from pathlib import Path
from urllib.parse import urlencode

from web.jobs import RESUMABLE, TERMINAL

STYLE = '''
:root{color-scheme:light;--ink:#192d38;--muted:#586b77;--line:#d9e1e5;--blue:#18576b}
*{box-sizing:border-box}body{margin:0;background:#f3f6f7;color:var(--ink);font:15px/1.5 system-ui,sans-serif}
header{background:#fff;border-bottom:1px solid var(--line);padding:22px max(24px,calc((100vw - 1220px)/2))}
header a{font-size:23px;font-weight:750;color:var(--ink);text-decoration:none}header p{margin:4px 0;color:var(--muted)}
main{max-width:1220px;margin:26px auto;padding:0 24px}h1{font-size:27px;margin:0 0 18px}h2{font-size:20px;margin:0 0 14px}h3{font-size:16px}
a{color:var(--blue)}.grid{display:grid;grid-template-columns:minmax(0,2fr) minmax(280px,1fr);gap:22px}
.card{background:#fff;border:1px solid var(--line);border-radius:12px;padding:23px;margin-bottom:22px}
.fields{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:15px}.wide{grid-column:1/-1}
label{display:block;font-weight:600;margin-bottom:5px}input,select{width:100%;padding:10px;border:1px solid #b9c8cf;border-radius:6px;font:inherit;background:#fff}
input[type=checkbox]{width:auto;margin-right:7px}button,.button{display:inline-block;background:var(--blue);color:#fff;border:0;border-radius:6px;padding:11px 16px;font:inherit;font-weight:650;text-decoration:none;cursor:pointer}
button.secondary,.secondary{background:#e6eef1;color:var(--ink)}button.danger{background:#963d37}form.inline{display:inline-block;margin:0 8px 8px 0}
small,.muted{color:var(--muted)}small{display:block;font-weight:400;margin-top:4px}.badge{display:inline-block;font-size:12px;padding:3px 9px;border-radius:20px;background:#e6eef1}
.succeeded{background:#ddf2e5;color:#246544}.failed,.partial{background:#fce7df;color:#963d37}.running{background:#dceff5;color:#18576b}
table{border-collapse:collapse;width:100%;font-size:14px}th,td{text-align:left;padding:10px 7px;border-bottom:1px solid var(--line);vertical-align:top}th{font-weight:650}td.path{overflow-wrap:anywhere;font-size:12px}
.scroll{overflow:auto}.notice{padding:13px 16px;background:#e6f0f4;border-radius:6px;margin-bottom:16px}.error{background:#fce7df;color:#87382e}
pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#132731;color:#dce9ed;padding:18px;border-radius:7px;max-height:560px;overflow:auto;font:12px/1.6 ui-monospace,monospace}
details{margin:18px 0}summary{cursor:pointer;font-weight:650}details .fields{margin-top:15px}.actions{margin-top:18px;display:flex;gap:10px;flex-wrap:wrap}img{max-width:100%;height:auto}.metric{font-size:28px;font-weight:750}.metrics{display:grid;grid-template-columns:repeat(4,1fr);gap:12px}.metrics small{margin:0}
@media(max-width:850px){.grid{grid-template-columns:1fr}.fields{grid-template-columns:1fr}.metrics{grid-template-columns:repeat(2,1fr)}}
'''


def e(value):
    return escape(str(value), quote=True)


def page(title, content, refresh=False):
    meta = '<meta http-equiv="refresh" content="5">' if refresh else ''
    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">{meta}<title>{e(title)} · PyMDA</title><style>{STYLE}</style></head><body><header><a href="/">PyMDA</a> <span class="badge">Local workstation</span><p>Diffraction processing and data assembly</p></header><main>{content}</main></body></html>'''


def field(name, title, value, kind='number', note='', attrs=''):
    return f'<div><label for="{e(name)}">{e(title)}</label><input id="{e(name)}" name="{e(name)}" type="{kind}" value="{e(value if value is not None else "")}" {attrs}>{"<small>" + e(note) + "</small>" if note else ""}</div>'


def select(name, title, choices, current):
    options = ''.join(f'<option value="{e(value)}" {"selected" if value == current else ""}>{e(label)}</option>' for value, label in choices)
    return f'<div><label for="{e(name)}">{e(title)}</label><select name="{e(name)}" id="{e(name)}">{options}</select></div>'


def parameters(mode):
    fields = []
    if mode == 'process':
        fields += [field('reso', 'Integration resolution (Å)', 1.8, attrs='step="0.01" min="0.1" required'),
                   field('reso_cchalf', 'Prefix assessment resolution (Å)', 2.5, attrs='step="0.01" min="0.1" required'),
                   field('thread', 'Image-processing workers', 2, attrs='min="1" max="2" required'),
                   field('wedges', 'Cumulative frame prefixes', 4, attrs='min="1" max="20" required')]
    fields += [field('resolutions', 'Merging resolution cutoffs (Å)', '2.0 2.2 2.26 2.3 2.5 3.0', 'text', 'Separate values with spaces or commas.'),
               field('merge_workers', 'Concurrent merge comparisons', 2, attrs='min="1" max="4" required')]
    advanced = [field('min_completeness', 'Minimum overall completeness (%)', 90, attrs='step="0.1" min="0" max="100"'),
                field('min_outer_completeness', 'Minimum outer-shell completeness (%)', 80, attrs='step="0.1" min="0" max="100"'),
                field('min_cc_half', 'Minimum outer-shell CC½', 0.3, attrs='step="0.01" min="0" max="1"'),
                field('min_isigi', 'Minimum outer-shell I/σ(I)', 1.0, attrs='step="0.1" min="0"'),
                select('prefix_policy', 'Frame-prefix comparison', [('all','I/σ(I), CC½ and full scans'),('selected','I/σ(I) prefixes'),('cchalf','CC½ prefixes'),('full','Full scans'),('both','I/σ(I) prefixes and full scans')], 'all'),
                select('crystal_identity', 'Crystal identity', [('unknown','Unknown'),('repeat','Repeat scans of one crystal'),('independent','Independent crystals')], 'unknown'),
                field('nres', 'Residues in the asymmetric unit (optional)', '', attrs='min="1"')]
    if mode == 'process':
        advanced += [field('frames', 'Maximum frames per scan (optional pilot)', '', attrs='min="1"'),
                     field('spg', 'Known space group (optional)', '', 'text'),
                     field('unit_cell', 'Known unit cell (optional)', '', 'text', 'a, b, c, α, β, γ'),
                     '<div><label><input type="checkbox" name="fix_detector" value="true">Fix detector geometry</label></div>']
    return '<div class="fields">' + ''.join(fields) + '</div><details><summary>Advanced parameters</summary><div class="fields">' + ''.join(advanced) + '</div></details>'


def job_list(jobs):
    if not jobs:
        return '<p class="muted">No jobs yet. Select data and submit a run.</p>'
    rows = ''.join(f'<tr><td><a href="/jobs/{e(j["id"])}">{e(j["id"][:10])}</a><small>{e(j["params"]["mode"])} · {e(j["created"][:16].replace("T"," "))} UTC</small></td><td><span class="badge {e(j["status"])}">{e(j["status"])}</span><small>{e(j["stage"])}</small></td></tr>' for j in jobs)
    return '<table><thead><tr><th>Job</th><th>Status</th></tr></thead><tbody>' + rows + '</tbody></table>'


def home(settings, jobs, input_dir, datasets, source_dir, sources, error=''):
    notice = f'<div class="notice error">{e(error)}</div>' if error else ''
    rows = ''
    for item in datasets:
        ready = item['ready']
        companions = ', '.join(x['name'] for x in item.get('companions', []) if not x['exists'])
        rows += f'<tr><td><input aria-label="Select {e(item["name"])}" type="checkbox" name="master_files" value="{e(item["path"])}" {"checked" if ready else "disabled"}></td><td class="path">{e(item["name"])}<small>{e(item["path"])}</small></td><td>{e(item.get("frames") or "—")}</td><td>{"Ready" if ready else e(item.get("error", "Unavailable"))}<small>{e(companions)}</small></td></tr>'
    choices = ''.join(f'<tr><td><input type="checkbox" name="dataset_ids" aria-label="Select {e(s["id"])}" value="{e(s["id"])}" checked></td><td class="path">{e(s["name"])}<small>{e(s["id"])}</small></td><td>{e(s.get("frame", "—"))}</td></tr>' for s in sources)
    raw = f'''<section class="card"><h2>Process diffraction data</h2><p class="muted">Select rotation master files. Companion HDF5 data stays on this workstation.</p><form method="get" action="/"><label for="input_dir">Input directory</label><input id="input_dir" name="input_dir" value="{e(input_dir)}"><input type="hidden" name="source_dir" value="{e(source_dir)}"><div class="actions"><button class="secondary">Discover scans</button></div></form><form method="post" action="/jobs"><input type="hidden" name="input_dir" value="{e(input_dir)}"><div class="scroll"><table><thead><tr><th></th><th>Master file</th><th>Frames</th><th>Companions</th></tr></thead><tbody>{rows or '<tr><td colspan="4">No usable master files found.</td></tr>'}</tbody></table></div><details open><summary>Processing parameters</summary>{parameters('process')}</details><div class="actions"><button name="mode" value="process">Process selected scans</button><button class="secondary" name="mode" value="both">Process and assemble</button></div></form></section>'''
    assembly = f'''<section class="card"><h2>Assemble processed datasets</h2><p class="muted">Use a completed processing job or an existing PyMDA output directory. A new job receives its own copy of the exported data.</p><form method="get" action="/"><label for="source_dir">Processed-data directory</label><input id="source_dir" name="source_dir" value="{e(source_dir)}"><input type="hidden" name="input_dir" value="{e(input_dir)}"><div class="actions"><button class="secondary">Find processed datasets</button></div></form><form method="post" action="/jobs"><input type="hidden" name="mode" value="assemble"><input type="hidden" name="source_dir" value="{e(source_dir)}"><div class="scroll"><table><thead><tr><th></th><th>Dataset</th><th>Selected last frame</th></tr></thead><tbody>{choices or '<tr><td colspan="3">No completed dataset manifests found.</td></tr>'}</tbody></table></div><details open><summary>Assembly parameters</summary>{parameters('assemble')}</details><div class="actions"><button>Assemble selected datasets</button></div></form></section>'''
    side = '<aside><section class="card"><h2>Jobs</h2><p class="muted">One processing job runs at a time. Additional jobs wait in the queue.</p>' + job_list(jobs) + '</section><section class="card"><h2>Output</h2><p class="path muted">' + e(settings.jobs_root) + '</p><p>Each job keeps its parameters, logs, plots, statistics and downloadable MTZ files.</p></section></aside>'
    return page('Workspace', '<h1>Run PyMDA</h1>' + notice + '<div class="grid"><div>' + raw + assembly + '</div>' + side + '</div>')


def load_json(path, default=None):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return default


def file_link(identifier, path, label=None, inline=False):
    query = urlencode({'path': str(path), 'inline': 'true' if inline else 'false'})
    return f'<a href="/jobs/{e(identifier)}/file?{e(query)}">{e(label or path)}</a>'


def tail(path, size=24000):
    try:
        with path.open('rb') as stream:
            stream.seek(max(0, path.stat().st_size - size))
            return stream.read().decode(errors='replace')
    except OSError:
        return 'The log will appear when the worker starts this job.'


def detail(job, refresh=True, log_file='job.stdout.log'):
    identifier = job['id']
    output = Path(job['output_dir'])
    active = job['status'] not in TERMINAL
    controls = ''
    if job['status'] in {'queued', 'running'}:
        controls += f'<form class="inline" method="post" action="/jobs/{identifier}/cancel"><button class="danger">Cancel job</button></form>'
    if job['status'] in RESUMABLE:
        controls += f'<form class="inline" method="post" action="/jobs/{identifier}/resume"><button>Resume completed stages</button></form>'
    if list(output.glob('crystals/*/manifest.json')):
        controls += f'<a class="button secondary" href="/?{e(urlencode({"source_dir": str(output)}))}">New assembly from these datasets</a>'
    top = f'<h1>Job {e(identifier[:10])}</h1><section class="card"><span class="badge {e(job["status"])}">{e(job["status"])}</span><h2>{e(job["stage"])}</h2><p>{e(job["message"])}</p><p class="muted">{e(job["params"]["mode"])} · Created {e(job["created"])}<br>Output: {e(output)}</p><div class="actions">{controls}</div></section>'
    manifests = list(output.glob('crystals/*/manifest.json'))
    trials = load_json(output / 'trials.json', [])
    top += f'<section class="card"><h2>Progress</h2><p>{len(manifests)} dataset manifests available · {len(trials)} merge configurations compared</p><p class="muted">Completed artifacts are retained if the job is cancelled or interrupted.</p></section>'
    best = load_json(output / 'best.json')
    results = ''
    if best:
        stats = best['statistics']
        results = '<section class="card"><h2>' + ('Results' if job['status'] == 'succeeded' else 'Last available results') + '</h2><div class="metrics">'
        for value, title in [(f"{stats['resolution'][0]:.2f} Å", 'Resolution'), (f"{stats['isigi'][0]:.1f}", 'I/σ(I), overall'), (f"{stats['cc_half'][2]:.3f}", 'CC½, outer shell'), (f"{stats['completeness'][0]:.1f}%", 'Completeness, overall')]:
            results += f'<div><div class="metric">{e(value)}</div><small>{e(title)}</small></div>'
        results += f'</div><p>{len(best["crystals"])} input datasets selected. Outer-shell I/σ(I): {stats["isigi"][2]:.1f}; outer-shell completeness: {stats["completeness"][2]:.1f}%.</p>'
        if best.get('crystal_identity', 'unknown') == 'unknown':
            results += '<p class="notice">Crystal identity is unverified. This result combines datasets; it does not establish that the scans came from independent crystals.</p>'
        alternatives = load_json(output/'alternatives.json', {})
        rows = ''
        labels = {'strongest_signal': 'Strongest signal', 'extended_resolution': 'Extended resolution', 'completeness_95': 'At least 95% complete'}
        for name, item in alternatives.items():
            s = item['statistics']
            rows += f'<tr><td>{e(labels.get(name,name))}</td><td>{s["resolution"][0]:.2f}</td><td>{s["isigi"][0]:.1f} / {s["isigi"][2]:.1f}</td><td>{s["completeness"][0]:.1f}% / {s["completeness"][2]:.1f}%</td><td>{file_link(identifier, "alternatives/"+name+"/scaled.mtz", "MTZ")}</td></tr>'
        if rows:
            results += '<div class="scroll"><table><thead><tr><th>Alternative</th><th>Å</th><th>I/σ(I), overall / outer</th><th>Completeness, overall / outer</th><th>Download</th></tr></thead><tbody>' + rows + '</tbody></table></div>'
        links = [('best/final.mtz','Final amplitudes and Free-R flags'),('best/scaled.mtz','Merged intensities'),('processing_report.md','Processing report'),('trials.csv','Comparison CSV'),('trials.json','Comparison JSON'),('pareto.json','Resolution/signal choices')]
        results += '<p>' + ' · '.join(file_link(identifier, p, label) for p, label in links if (output/p).is_file()) + '</p>'
        if (output/'resolution_tradeoff.png').is_file():
            url = f'/jobs/{identifier}/file?' + urlencode({'path':'resolution_tradeoff.png','inline':'true'})
            results += f'<img src="{e(url)}" alt="Resolution, signal and completeness comparison">'
        results += '</section>'
    if manifests:
        rows = ''
        for manifest in sorted(manifests):
            record = load_json(manifest, {})
            for prefix in record.get('prefixes', []):
                stats = prefix['statistics']
                mtz = Path(prefix['mtz'])
                if not mtz.is_relative_to(output):
                    continue
                rows += f'<tr><td>{e(manifest.parent.name)}<small>Through frame {e(prefix["last_frame"])}</small></td><td>{stats["resolution"][0]:.2f}</td><td>{stats["isigi"][0]:.1f} / {stats["isigi"][2]:.1f}</td><td>{stats["completeness"][0]:.1f}%</td><td>{file_link(identifier, mtz.relative_to(output), "Intensity MTZ")}</td></tr>'
        results += '<section class="card"><h2>Processed scan prefixes</h2><div class="scroll"><table><thead><tr><th>Dataset</th><th>Assessment Å</th><th>I/σ(I), overall / outer</th><th>Completeness</th><th>Download</th></tr></thead><tbody>' + rows + '</tbody></table></div></section>'
    paths = [Path('job.stdout.log')]
    if output.is_dir():
        paths += [p.relative_to(output) for p in sorted(output.rglob('*.console.log'), key=lambda p:p.stat().st_mtime, reverse=True)[:200]]
    if Path(log_file) not in paths:
        paths.insert(0, Path(log_file))
    options = ''.join(f'<option value="{e(p)}" {"selected" if str(p)==log_file else ""}>{e(p)}</option>' for p in paths)
    refresh_label = 'Pause refresh' if refresh else 'Resume refresh'
    query = urlencode({'refresh': '0' if refresh else '1', 'log_file': log_file})
    logs = f'<section class="card" id="log"><h2>Logs</h2><p>{"Page refreshes every five seconds while the job is active." if active and refresh else "Automatic refresh is paused or the job has finished."} <a href="?{e(query)}">{refresh_label}</a></p><form method="get"><label for="log_file">Log window</label><select id="log_file" name="log_file">{options}</select><input type="hidden" name="refresh" value="{int(refresh)}"><div class="actions"><button class="secondary">Open log</button></div></form><pre>{e(tail(output/log_file))}</pre><p>{file_link(identifier,log_file,"Download full log") if (output/log_file).is_file() else ""}</p></section>'
    params = '<section class="card"><details><summary>Submitted parameters</summary><pre>' + e(json.dumps(job['params'],indent=2)) + '</pre></details></section>'
    return page('Job ' + identifier[:10], top + results + logs + params, active and refresh)

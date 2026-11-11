"""HTML pages rendered entirely by Python, with no templates or JavaScript."""
from html import escape
import json
from math import ceil
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


def field(name, title, value, kind='number', note='', attrs='', identifier=None):
    identifier = identifier or name
    return f'<div><label for="{e(identifier)}">{e(title)}</label><input id="{e(identifier)}" name="{e(name)}" type="{kind}" value="{e(value if value is not None else "")}" {attrs}>{"<small>" + e(note) + "</small>" if note else ""}</div>'


def select(name, title, choices, current):
    options = ''.join(f'<option value="{e(value)}" {"selected" if value == current else ""}>{e(label)}</option>' for value, label in choices)
    return f'<div><label for="{e(name)}">{e(title)}</label><select name="{e(name)}" id="{e(name)}">{options}</select></div>'


def parameters(mode, integration_resolution=1.8):
    fields = []
    cutoffs = [2.0, 2.2, 2.26, 2.3, 2.5, 3.0]
    comparison_resolution = 2.5
    if mode == 'assemble':
        floor = ceil(integration_resolution*100)/100
        cutoffs = [value for value in cutoffs if value >= floor] or [floor]
        comparison_resolution = max(comparison_resolution, floor)
    if mode == 'process':
        fields += [field('reso', 'Integration resolution (Å)', 1.8, attrs='step="0.01" min="0.1" required'),
                   field('reso_cchalf', 'Prefix assessment resolution (Å)', 2.5, attrs='step="0.01" min="0.1" required'),
                   field('thread', 'Image-processing workers', 2, attrs='min="1" max="2" required'),
                   field('split_wedges', 'Consecutive wedges per scan', 1, note='1 keeps the whole scan. With 1,800 frames, 3 gives 1–600, 601–1200, 1201–1800; each wedge gets an MTZ and statistics.', attrs='min="1" max="100" step="1" required'),
                   field('wedges', 'Cumulative frame prefixes', 4, attrs='min="1" max="20" required')]
    fields += [field('resolutions', 'Merging resolution cutoffs (Å)', ' '.join(map(str,cutoffs)), 'text', 'Separate values with spaces or commas.'),
               field('merge_workers', 'Concurrent merge comparisons', 2, attrs='min="1" max="4" required')]
    advanced = [field('min_completeness', 'Minimum overall completeness (%)', 90, attrs='step="0.1" min="0" max="100"'),
                field('min_outer_completeness', 'Minimum outer-shell completeness (%)', 80, attrs='step="0.1" min="0" max="100"'),
                field('min_cc_half', 'Minimum outer-shell CC½', 0.3, attrs='step="0.01" min="0" max="1"'),
                field('min_isigi', 'Minimum outer-shell I/σ(I)', 1.0, attrs='step="0.1" min="0"'),
                select('subset_search', 'Dataset subset search', [('auto','Automatic: exhaustive through 10 datasets, bounded above 10'),('bounded','Bounded search'),('exhaustive','Exhaustive search (up to 10 datasets)')], 'auto'),
                field('max_subsets', 'Maximum bounded-search subsets', 64, note='Includes the full collection and a limited set of alternatives.', attrs='min="1" max="4096"'),
                select('prefix_policy', 'Frame-prefix comparison', [('all','I/σ(I), CC½ and full scans'),('selected','I/σ(I) prefixes'),('cchalf','CC½ prefixes'),('full','Full scans'),('both','I/σ(I) prefixes and full scans')], 'all'),
                select('crystal_identity', 'Crystal identity', [('unknown','Unknown'),('repeat','Repeat scans of one crystal'),('independent','Independent crystals')], 'unknown'),
                field('nres', 'Residues in the asymmetric unit (optional)', '', attrs='min="1"')]
    if mode == 'process':
        advanced += [field('frames', 'Maximum frames per scan (optional pilot)', '', attrs='min="1"'),
                     field('spg', 'Known space group (optional)', '', 'text'),
                     field('unit_cell', 'Known unit cell (optional)', '', 'text', 'a, b, c, α, β, γ'),
                     '<div><label><input type="checkbox" name="fix_detector" value="true">Fix detector geometry</label></div>']
    if mode == 'assemble':
        advanced += [field('cell_length_tolerance','Maximum cell length difference (%)',2,attrs='step="0.1" min="0.01"'), field('cell_angle_tolerance','Maximum cell angle difference (degrees)',2,attrs='step="0.1" min="0.01"'),field('min_pair_cc','Minimum pairwise intensity correlation',0.9,attrs='step="0.01" min="-1" max="1"'),field('min_common','Minimum common reflections',100,attrs='min="3"'),field('reso_cchalf','Comparison resolution (Å)',comparison_resolution,attrs='step="0.01" min="0.1"')]
    return '<div class="fields">' + ''.join(fields) + '</div><details><summary>Advanced parameters</summary><div class="fields">' + ''.join(advanced) + '</div></details>'


def job_list(jobs):
    if not jobs:
        return '<p class="muted">No jobs yet. Select data and submit a run.</p>'
    rows = ''.join(f'<tr><td><a href="/jobs/{e(j["id"])}">{e(j["id"][:10])}</a><small>{e(j["params"]["mode"])} · {e(j["created"][:16].replace("T"," "))} UTC</small></td><td><span class="badge {e(j["status"])}">{e(j["status"])}</span><small>{e(j["stage"])}</small></td></tr>' for j in jobs)
    return '<table><thead><tr><th>Job</th><th>Status</th></tr></thead><tbody>' + rows + '</tbody></table>'


def home(settings, jobs, input_dir, datasets, source_dir, sources, error='', output_dir='', comparison=None, selected=None):
    comparison = comparison or {}
    def browse(kind, label):
        query = urlencode(dict(kind=kind,input_dir=input_dir,source_dir=source_dir,output_dir=output_dir,path={'input_dir':input_dir,'source_dir':source_dir,'output_dir':output_dir}[kind]))
        return '<a class="button secondary" href="/browse?'+e(query)+'">'+label+'</a>'
    def find_outputs():
        query = urlencode(dict(input_dir=input_dir,source_dir=output_dir,output_dir=output_dir))
        return '<a class="button secondary" href="/?'+e(query)+'#processed-datasets">Find processed datasets in this output folder</a>'
    def output_field(identifier):
        return field('output_dir','Output directory',output_dir,'text',attrs='required',identifier=identifier)+'<div class="actions">'+browse('output_dir','Choose output directory')+find_outputs()+'</div><small>Each job gets its own folder inside the selected directory.</small>'
    notice = f'<div class="notice error">{e(error)}</div>' if error else ''
    rows = ''
    for item in datasets:
        ready = item['ready']
        companions = ', '.join(x['name'] for x in item.get('companions', []) if not x['exists'])
        rows += f'<tr><td><input aria-label="Select {e(item["name"])}" type="checkbox" name="master_files" value="{e(item["path"])}" {"checked" if ready else "disabled"}></td><td class="path">{e(item["name"])}<small>{e(item["path"])}</small></td><td>{e(item.get("frames") or "—")}</td><td>{"Ready" if ready else e(item.get("error", "Unavailable"))}<small>{e(companions)}</small></td></tr>'
    choices, selected_resolutions = '', []
    for item in sources:
        ready = item['ready']
        checked = item['default_selected'] if selected is None else item['id'] in selected
        attributes = 'disabled' if not ready else 'checked' if checked else ''
        if ready and checked:
            selected_resolutions.append(item['integration_resolution'])
        query = urlencode(dict(input_dir=input_dir,source_dir=item['source_dir'],output_dir=output_dir))
        location = '<small>'+e(item['source_dir'])+'<br><a href="/?'+e(query)+'#processed-datasets">Open this result directory</a></small>'
        status = ('Ready · selected by default' if item['default_selected'] else 'Ready') if ready else e(item['error'])
        choices += f'<tr><td><input type="checkbox" name="dataset_ids" aria-label="Select {e(item["id"])}" value="{e(item["id"])}" {attributes}></td><td class="path">{e(item["name"])}{location}</td><td>{e(item.get("frame") or "—")}</td><td>{e(item.get("integration_resolution") or "—")}</td><td>{status}</td></tr>'
    assembly_resolution = max(selected_resolutions, default=1.8)
    raw = f'''<section class="card"><h2>Process diffraction data</h2><p class="muted">Select rotation master files. Companion HDF5 data stays on this workstation. Scale one scan for single-crystal results; process several scans first to compare and select compatible datasets.</p><form method="get" action="/"><input type="hidden" name="output_dir" value="{e(output_dir)}"><label for="input_dir">Input directory</label><input id="input_dir" name="input_dir" value="{e(input_dir)}"><input type="hidden" name="source_dir" value="{e(source_dir)}"><div class="actions"><button class="secondary">Discover scans</button>{browse("input_dir", "Browse data directory")}</div></form><form method="post" action="/jobs"><input type="hidden" name="input_dir" value="{e(input_dir)}"><div class="scroll"><table><thead><tr><th></th><th>Master file</th><th>Frames</th><th>Companions</th></tr></thead><tbody>{rows or '<tr><td colspan="4">No usable master files found.</td></tr>'}</tbody></table></div>{output_field('process_output_dir')}<details open><summary>Processing parameters</summary>{parameters('process')}</details><div class="actions"><button name="mode" value="process">Process selected scans</button><button class="secondary" name="mode" value="both">Process and scale / assemble</button></div></form></section>'''
    assembly = f'''<section class="card" id="processed-datasets"><h2>Assemble processed datasets</h2><p class="muted">Choose your output folder, a job folder, a result directory, or a crystal folder. The search includes job subfolders and separate runs. One completed version per scan is selected by default; you can select other versions below. A new job receives its own copy of the exported data.</p><form method="get" action="/"><input type="hidden" name="output_dir" value="{e(output_dir)}"><label for="source_dir">Processed-data directory</label><input id="source_dir" name="source_dir" value="{e(source_dir)}"><input type="hidden" name="input_dir" value="{e(input_dir)}"><div class="actions"><button class="secondary">Find processed datasets</button>{browse("source_dir", "Browse processed results")}</div></form><form method="post" action="/jobs"><input type="hidden" name="source_dir" value="{e(source_dir)}"><div class="scroll"><table><thead><tr><th></th><th>Dataset / result directory</th><th>Selected last frame</th><th>Integration limit (Å)</th><th>Status</th></tr></thead><tbody>{choices or '<tr><td colspan="5">No processed datasets found in this folder or its job subfolders.</td></tr>'}</tbody></table></div>{output_field('assembly_output_dir')}<details open><summary>Assembly parameters</summary><p class="muted">Defaults use the selected results' integration limits. Check the cutoffs if you change the dataset selection.</p>{parameters('assemble',assembly_resolution)}</details><div class="actions"><button name="mode" value="assemble">Scale / assemble selected datasets</button><button class="secondary" name="mode" value="compare">Compare selected datasets</button></div></form></section>'''
    selection = '<section class="card"><h2>Output destination</h2><p>'+e(output_dir)+'</p><div class="actions">'+browse('output_dir','Choose output directory')+find_outputs()+'</div></section>'
    if comparison:
        assembly = comparison_html(comparison,source_dir,input_dir,output_dir) + assembly
    side = '<aside><section class="card"><h2>Jobs</h2><p class="muted">One processing job runs at a time. Additional jobs wait in the queue.</p>' + job_list(jobs) + '</section><section class="card"><h2>Output</h2><p class="path muted">' + e(settings.jobs_root) + '</p><p>Each job keeps its parameters, logs, plots, statistics and downloadable MTZ files.</p></section></aside>'
    return page('Workspace', '<h1>Run PyMDA</h1>' + notice + '<div class="grid"><div>' + selection + raw + assembly + '</div>' + side + '</div>')


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
    search = load_json(output/'search.json')
    if search:
        top += '<section class="card"><h2>Dataset subset search</h2><p>'+e(search['mode'])+' search: '+e(search['candidate_count'])+' candidate subsets for '+e(search['dataset_count'])+' datasets.</p><p>'+file_link(identifier,'search.json','Download evaluated subset plan')+'</p></section>'
    best = load_json(output / 'best.json')
    results = ''
    if best:
        stats = best['statistics']
        results = '<section class="card"><h2>' + ('Results' if job['status'] == 'succeeded' else 'Last available results') + '</h2><div class="metrics">'
        for value, title in [(f"{stats['resolution'][0]:.2f} Å", 'Resolution'), (f"{stats['isigi'][0]:.1f}", 'I/σ(I), overall'), (f"{stats['cc_half'][2]:.3f}", 'CC½, outer shell'), (f"{stats['completeness'][0]:.1f}%", 'Completeness, overall')]:
            results += f'<div><div class="metric">{e(value)}</div><small>{e(title)}</small></div>'
        results += f'</div><p>{len(best["crystals"])} input datasets selected. Outer-shell I/σ(I): {stats["isigi"][2]:.1f}; outer-shell completeness: {stats["completeness"][2]:.1f}%.</p>'
        if len(best['crystals']) > 1 and best.get('crystal_identity', 'unknown') == 'unknown':
            results += '<p class="notice">Crystal identity is unverified. This result combines datasets; it does not establish that the scans came from independent crystals.</p>'
        alternatives = load_json(output/'alternatives.json', {})
        rows = ''
        labels = {'strongest_signal': 'Strongest signal', 'extended_resolution': 'Extended resolution', 'completeness_95': 'At least 95% complete'}
        for name, item in alternatives.items():
            s = item['statistics']
            rows += f'<tr><td>{e(labels.get(name,name))}</td><td>{s["resolution"][0]:.2f}</td><td>{s["isigi"][0]:.1f} / {s["isigi"][2]:.1f}</td><td>{s["completeness"][0]:.1f}% / {s["completeness"][2]:.1f}%</td><td>{file_link(identifier, "alternatives/"+name+"/scaled.mtz", "MTZ")}</td></tr>'
        if rows:
            results += '<div class="scroll"><table><thead><tr><th>Alternative</th><th>Å</th><th>I/σ(I), overall / outer</th><th>Completeness, overall / outer</th><th>Download</th></tr></thead><tbody>' + rows + '</tbody></table></div>'
        links = [('best/final.mtz','Final amplitudes and Free-R flags'),('best/scaled.mtz','Merged intensities'),('processing_report.md','Processing report'),('trials.csv','Comparison CSV'),('trials.json','Comparison JSON'),('pareto.json','Resolution/signal choices'),('selection.json','Dataset selection'),('comparison.json','Compatibility comparison'),('search.json','Subset search details')]
        results += '<p>' + ' · '.join(file_link(identifier, p, label) for p, label in links if (output/p).is_file()) + '</p>'
        if (output/'resolution_tradeoff.png').is_file():
            url = f'/jobs/{identifier}/file?' + urlencode({'path':'resolution_tradeoff.png','inline':'true'})
            results += f'<img src="{e(url)}" alt="Resolution, signal and completeness comparison">'
        results += '</section>'
    comparison = load_json(output/'comparison.json', {})
    if trials and not best:
        rows = ''
        for trial in trials:
            s = trial['statistics']
            rows += f'<tr><td>{e(trial["prefix_policy"])}</td><td>{s["resolution"][0]:.2f}</td><td>{s["isigi"][0]:.1f} / {s["isigi"][2]:.1f}</td><td>{s["cc_half"][2]:.3f}</td><td>{s["completeness"][0]:.1f}% / {s["completeness"][2]:.1f}%</td><td>{"Passes" if trial["eligible"] else "Below thresholds"}</td></tr>'
        results += '<section class="card"><h2>Scaling trials</h2><p>'+('No trial met the configured quality thresholds. Inspect completeness and signal before changing the criteria.' if not any(t['eligible'] for t in trials) else 'Finalization is pending; these scaling statistics are available.')+'</p><div class="scroll"><table><tr><th>Prefixes</th><th>Å</th><th>I/σ(I), overall / outer</th><th>CC½, outer</th><th>Completeness, overall / outer</th><th>Quality</th></tr>'+rows+'</table></div><p>'+file_link(identifier,'trials.json','Trial statistics JSON')+'</p></section>'
    if comparison:
        results += '<p>'+file_link(identifier,'comparison.json','Download compatibility comparison')+'</p>'
        results += comparison_html(comparison,str(output),'','')
    if manifests:
        rows = ''
        wedge_rows = ''
        for manifest in sorted(manifests):
            record = load_json(manifest, {})
            for prefix in record.get('prefixes', []):
                stats = prefix['statistics']
                mtz = Path(prefix['mtz'])
                if not mtz.is_relative_to(output):
                    continue
                rows += f'<tr><td>{e(manifest.parent.name)}<small>Through frame {e(prefix["last_frame"])}</small></td><td>{stats["resolution"][0]:.2f}</td><td>{stats["isigi"][0]:.1f} / {stats["isigi"][2]:.1f}</td><td>{stats["completeness"][0]:.1f}%</td><td>{file_link(identifier, mtz.relative_to(output), "Intensity MTZ")}</td></tr>'
            for wedge in record.get('scan_wedges', []):
                first, last = wedge['first_frame'], wedge['last_frame']
                frames = f'{first}–{last}'
                used_first, used_last = wedge.get('used_first_frame'), wedge.get('used_last_frame')
                if used_first is not None and (used_first, used_last) != (first, last):
                    frames += f'<small>Integrated frames {used_first}–{used_last}; DIALS trimmed endpoints.</small>'
                stats = wedge.get('statistics')
                download = ''
                if wedge.get('mtz'):
                    mtz = Path(wedge['mtz'])
                    if mtz.is_relative_to(output):
                        download = file_link(identifier, mtz.relative_to(output), 'Wedge MTZ')
                values = (f'<td>{stats["resolution"][0]:.2f}</td><td>{stats["isigi"][0]:.1f} / {stats["isigi"][2]:.1f}</td><td>{stats["cc_half"][2]:.3f}</td><td>{stats["completeness"][0]:.1f}% / {stats["completeness"][2]:.1f}%</td>'
                          if stats else '<td colspan="4">'+e(wedge.get('error') or 'Statistics unavailable')+'</td>')
                wedge_rows += f'<tr><td>{e(manifest.parent.name)}<small>Wedge {e(wedge["number"])}</small></td><td>{frames}</td>{values}<td>{download}</td></tr>'
        if wedge_rows:
            results += '<section class="card"><h2>Consecutive scan wedges</h2><p>Each wedge covers its own frame range and has separate quality statistics.</p><div class="scroll"><table><thead><tr><th>Dataset / wedge</th><th>Frames</th><th>Assessment Å</th><th>I/σ(I), overall / outer</th><th>CC½, outer</th><th>Completeness, overall / outer</th><th>Download</th></tr></thead><tbody>'+wedge_rows+'</tbody></table></div></section>'
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


def comparison_html(data, source_dir, input_dir='', output_dir=''):
    rows = ''.join('<tr><td>'+e(d['name'])+'</td><td>'+e(d['space_group'])+'</td><td>'+e(', '.join(f'{v:.3f}' for v in d['cell']))+'</td></tr>' for d in data['datasets'])
    body = '<section class="card"><h2>Crystal compatibility</h2><p>'+e(data['method'])+'</p><p>Criteria: '+e(data['criteria'])+'</p><table><tr><th>Dataset</th><th>Symmetry</th><th>Unit cell (Å, degrees)</th></tr>'+rows+'</table><div class="scroll"><table><tr><th>Intensity correlation</th>'
    lookup = {frozenset((p['a'],p['b'])):p for p in data['pairs']}
    for d in data['datasets']:
        body += '<th>'+e(d['name'])+'</th>'
    body += '</tr>'
    for d in data['datasets']:
        body += '<tr><th>'+e(d['name'])+'</th>'
        for other in data['datasets']:
            pair = lookup.get(frozenset((d['id'],other['id'])))
            if pair is None:
                body += '<td>1.000</td>'; continue
            cc = pair['correlation']
            color = '#ddf2e5' if pair['compatible'] else '#fce7df'
            body += '<td style="background:'+color+'">'+(f'{cc:.3f}' if cc is not None else 'Unavailable')+'<small>'+str(pair['common_reflections'])+' common; cell '+('passes' if pair['cell_compatible'] else 'differs')+'; lengths '+f"{pair['length_difference_percent']:.2f}%"+'; angles '+f"{pair['angle_difference_degrees']:.2f}°"+'</small><small>'+e(pair['error'])+'</small></td>'
        body += '</tr>'
    body += '</table></div><p>Suggested groups (every pair passes the criteria):</p>'
    for index, group in enumerate(data['groups'],1):
        body += '<p><a class="button secondary" href="/?'+e(urlencode(dict(source_dir=source_dir,input_dir=input_dir,output_dir=output_dir,group=index)))+'">Select group '+str(index)+' ('+str(len(group))+' datasets)</a></p>'
    return body+'</section>'

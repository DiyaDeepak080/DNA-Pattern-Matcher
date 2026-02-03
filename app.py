from flask import Flask, request, jsonify, send_from_directory, Response
from motif.parsing import parse_motif_to_symbol_sets, read_fasta
from motif.automata import build_linear_nfa, nfa_to_dfa
from motif.matching import find_motif_positions_dfa
import graphviz


app = Flask(__name__)

@app.route('/')
def index():
    return send_from_directory('static', 'index.html')

@app.route('/match', methods=['POST'])
def match_api():
    motif = request.form.get('motif', '')
    fasta_file = request.files.get('fasta')
    if not motif or not fasta_file:
        return jsonify({'detail': 'Missing motif or FASTA file.'}), 400
    try:
        content = fasta_file.read().decode('utf-8', errors='ignore')
        entries = read_fasta(content)
        symbol_sets = parse_motif_to_symbol_sets(motif)
        motif_len = len(symbol_sets)
        effective_alphabet = set().union(*symbol_sets)
        nfa = build_linear_nfa(symbol_sets)
        dfa = nfa_to_dfa(nfa, effective_alphabet)

        results = []
        for header, seq in entries:
            positions = find_motif_positions_dfa(seq, dfa, motif_len, effective_alphabet)
            examples = []
            k = 8
            for p in positions[:10]:
                start = max(0, p-k)
                end = min(len(seq), p+motif_len+k)
                frag = seq[start:end]
                rel_start, rel_end = p - start, p - start + motif_len
                marked = frag[:rel_start] + "[" + frag[rel_start:rel_end] + "]" + frag[rel_end:]
                examples.append(marked)
            results.append({
                "header": header,
                "length": len(seq),
                "count": len(positions),
                "positions": positions,
                "examples": examples,
            })
        return jsonify({"results": results})
    except Exception as e:
        return jsonify({"detail": str(e)}), 400

# Helper function to create NFA SVG using graphviz
def nfa_to_dot(nfa):
    dot = graphviz.Digraph(format='svg')
    dot.attr(rankdir='LR')
    dot.node('start', shape='point')
    dot.edge('start', str(nfa.start))
    for state in range(nfa.num_states):
        shape = 'doublecircle' if state in nfa.accepts else 'circle'
        dot.node(str(state), shape=shape)
    edges = {}
    for (src, symbol), dst_set in nfa.delta.items():
        for dst in dst_set:
            edges.setdefault((str(src), str(dst)), []).append(symbol)
    for (src, dst), symbols in edges.items():
        label = ",".join(symbols)
        dot.edge(src, dst, label=label)
    return dot

# Helper function to create DFA SVG using graphviz    
def dfa_to_dot(dfa):
    dot = graphviz.Digraph(format='svg')
    dot.attr(rankdir='LR')
    dot.node('start', shape='point')
    dot.edge('start', str(dfa.start))
    state_names = {state: f'S{idx}' for idx, state in enumerate(sorted(dfa.accepts.union({dfa.start}), key=str))}
    for state in state_names.values():
        dot.node(state, shape='circle')
    for acc in dfa.accepts:
        dot.node(state_names[acc], shape='doublecircle')
    for (src, symbol), dst in dfa.delta.items():
        src_name = state_names.get(src, str(src))
        dst_name = state_names.get(dst, str(dst))
        dot.edge(src_name, dst_name, label=symbol)
    return dot

@app.route('/nfa.svg')
def get_nfa_svg():
    motif = request.args.get('motif', '')
    if not motif:
        return Response('Missing motif query param', status=400)
    try:
        symbol_sets = parse_motif_to_symbol_sets(motif)
        nfa = build_linear_nfa(symbol_sets)
        dot = nfa_to_dot(nfa)
        svg = dot.pipe().decode('utf-8')
        return Response(svg, mimetype='image/svg+xml')
    except Exception as e:
        return Response(f'Error: {e}', status=500)

@app.route('/dfa.svg')
def get_dfa_svg():
    motif = request.args.get('motif', '')
    if not motif:
        return Response('Missing motif query param', status=400)
    try:
        symbol_sets = parse_motif_to_symbol_sets(motif)
        nfa = build_linear_nfa(symbol_sets)
        dfa = nfa_to_dfa(nfa, set().union(*symbol_sets))
        dot = dfa_to_dot(dfa)
        svg = dot.pipe().decode('utf-8')
        return Response(svg, mimetype='image/svg+xml')
    except Exception as e:
        return Response(f'Error: {e}', status=500)


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=8000, debug=True)
